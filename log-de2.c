/* 
 *this function computes log, correctly rounded, 
 using experimental techniques based on  double-extended arithmetic

 THIS IS EXPERIMENTAL SOFTWARE

In particular it changes rounding modes all the time without warning
nor restoring.
 
 *
 * Author :  Florent de Dinechin
 * Florent.de.Dinechin at ens-lyon.fr
 *

 To have it replace the crlibm log, do:
on pentium, 
 gcc -DHAVE_CONFIG_H -I.  -fPIC  -O2 -c log-de2.c;   mv log-de2.o log_fast.o; make 
 
*/


#include <stdio.h>
#include <stdlib.h>
#include "crlibm.h"
#include "crlibm_private.h"
#include "double-extended.h"
#include "log-de2.h"


double log_rn(double x) {
  double wi;
  db_number y;
  long double r, logirh, logirl, z, t;
  long double  th, tl, eh,el;
  int E, i, index, roundtestmask;


   y.d=x;

   E=0;
   /* Filter cases */
   if (y.i[HI] < 0x00100000){        /* x < 2^(-1022)    */
     if (((y.i[HI] & 0x7fffffff)|y.i[LO])==0){
       return -1.0/0.0;     
     }                    		   /* log(+/-0) = -Inf */
     if (y.i[HI] < 0){ 
       return (x-x)/0;                      /* log(-x) = Nan    */
     }
     /* Subnormal number */
     E = -64; 		
     y.d *= two64; 	  /* make x a normal number    */ 
   }
    
   if (y.i[HI] >= 0x7ff00000){
     return  x+x;				 /* Inf or Nan       */
   }

   
   DOUBLE_EXTENDED_MODE;  /* This one should be overlapped with integer computation */

   if((y.i[HI]>MINYFAST) && (y.i[HI]<MAXYFAST)) {
     roundtestmask=0x7fc;
     z = x - 1 ; /* Sterbenz exact */
     t = c6;
     t = c5 + z*t;
     t = c4 + z*t;
     t = c3 + z*t;
     t = c2 + z*t;
     t = c1 + z*t;
     t = z*t;
     //printf("z= %1.20e,  t=%1.20e  \n  ", (double)z, (double)t);
     
   }    

   else {
   /* Extract exponent and mantissa */
     E += (y.i[HI]>>20)-1023;             /* extract the exponent */
     index = (y.i[HI] & 0x000fffff);
     y.i[HI] =  index | 0x3ff00000;	/* do exponent = 0 */
     index = index  >> (20-L); 
     /* now y.d holds 1+f, and E is the exponent */
     
     if(E>24 || E<-24)
       roundtestmask=0x7fe;
     else
       roundtestmask=0x7f0;

     logirh = argredtable[index].h;
     r = (long double) (argredtable[index].r); /* approx to 1/y.d */

     z = r*y.d - 1. ; /* even without an FMA, all exact */
     
     //   printf("  z=%1.20e\n  r=%1.20e\n  logirh=%1.20e\n  ", (double)z, (double)r, (double)logirh);
     /* Polynomial evaluation, unrolled to go through Gappa */
     t = c6;
     t = c5 + z*t;
     t = c4 + z*t;
     t = c3 + z*t;
     t = c2 + z*t;
     t = c1 + z*t;

     //printf("t=%1.20e  \n  ", (double)t);
     
     /* reconstruction */
     t = logirh + z*t;
     t = t + E*ln2h;
     
     //  printf("index=%d\ncr_libm    : ", index);
   }
#if 0 /* to time the first step only */
   BACK_TO_DOUBLE_MODE; return (double)t;
#endif

   TEST_AND_RETURN_RN(t, roundtestmask);


   /* Accurate phase */ 
#if EVAL_PERF
  crlibm_second_step_taken++;
#endif
   
   t = c13h;
   t = c12h + z*t;
   t = c11h + z*t;
   t = c10h + z*t;
   t = c9h  + z*t;
   t = c8h  + z*t;
   Mul12_ext(&th, &tl,   z,   t);
   Add22_ext(&th, &tl,   c7h,c7l,  th,tl);
   Mul22_ext(&th, &tl,   z,0,  th,tl);
   Add22_ext(&th, &tl,   c6h,c6l,  th,tl);
   Mul22_ext(&th, &tl,   z,0,  th,tl);
   Add22_ext(&th, &tl,   c5h,c5l,  th,tl);
   Mul22_ext(&th, &tl,   z,0,  th,tl);
   Add22_ext(&th, &tl,   c4h,c4l,  th,tl);
   Mul22_ext(&th, &tl,   z,0,  th,tl);
   Add22_ext(&th, &tl,   c3h,c3l,  th,tl);
   Mul22_ext(&th, &tl,   z,0,  th,tl);
   Add22_ext(&th, &tl,   c2h,c2l,  th,tl);
   Mul22_ext(&th, &tl,   z,0,  th,tl);
   Add22_ext(&th, &tl,   c1h,c1l,  th,tl);

   Mul22_ext(&th, &tl,   z,0,  th,tl);
   Add22_ext(&th, &tl,   logirh, argredtable[index].l,  th,tl);
   
   Mul22_ext(&eh, &el,   ln2h,ln2l, E, 0);
   Add22_ext(&th, &tl,   eh,el,  th,tl);
   
   BACK_TO_DOUBLE_MODE;
   return (double) (th+tl); /* The exact sum of these double-extended is rounded to the nearest */
}


double log_ru(double x) { return x;};
double log_rd(double x) { return x;};
double log_rz(double x) { return x;};
