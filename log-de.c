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
 gcc -DHAVE_CONFIG_H -I.  -fPIC  -O2 -c log-de.c;   mv log-de.o log_fast.o; make 
 
*/


#include <stdio.h>
#include <stdlib.h>
#include "crlibm.h"
#include "crlibm_private.h"
#include "double-extended.h"
#include "log-de.h"


double log_rn(double x) {
  double wi;
  db_number y;
  long double logwih, logwil, invwih, invwil, z, z0, t, p1,p2,p3,z2,z4;
  long double zh, zl, th, tl, eh,el;
  int E, i, index;


   E=0;
   y.d=x;

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

   /* Extract exponent and mantissa */
   E += (y.i[HI]>>20)-1023;             /* extract the exponent */
   index = (y.i[HI] & 0x000fffff);
   y.i[HI] =  index | 0x3ff00000;	/* do exponent = 0 */
   index = (index + (1<<(20-L-1))) >> (20-L); 

   /* reduce to  y.d such that sqrt(2)/2 < y.d < sqrt(2) */
   if (index >= MAXINDEX){ /* corresponds to y>sqrt(2)*/
     //y.d *= 0.5;
     y.i[HI] -= 0x00100000; 
     E++;
   }

   /* now y.d holds 1+f, and E is the exponent */

   /*read the tables */
   wi     = argredtable[index].wi;
   invwih = argredtable[index].invwih;
   logwih    = argredtable[index].logwih;
 

   z0 = y.d-wi; /* exact thanks to Sterbenz */
   z  = z0*invwih;

   /* Polynomial evaluation, unrolled to go through Gappa */
#if 0
   t = c8;
   t = c7 + z*t;
   t = c6 + z*t;
   t = c5 + z*t;
   t = c4 + z*t;
   t = c3 + z*t;
   t = c2 + z*t;
   t = c1 + z*t;
   t = logwih + z*t;
#else
   z2 = z*z;             p1 = c2+c3*z;    p2 = c4+c5*z;     p3 = c6+c7*z;
   t = logwih + c1*z;    z4 = z2*z2;      p1 = p1+p2*z2;    p3 = p3+c8*z2;
   p1 = p1 + z4*p3;
   t = t + z2 *p1;   
#endif
   /* reconstruction */
   t = t + E*ln2h;

#if 0 /* to time the first step only */
   BACK_TO_DOUBLE_MODE; return (double)t;
#endif

#if 0
   TEST_AND_RETURN_RN_ZIV(t,1.001);
#else
   TEST_AND_RETURN_RN(t, 0x7fe);
#endif


   /* Accurate phase */ 
#if EVAL_PERF
  crlibm_second_step_taken++;
#endif

   Mul22_ext(&zh,&zl, argredtable[index].invwih, argredtable[index].invwil,  z0, 0);
   
   t = c15h;
   t = c14h + zh*t;
   t = c13h + zh*t;
   t = c12h + zh*t;
   t = c11h + zh*t;
   t = c10h + zh*t;
   t = c9h  + zh*t;
   Mul22_ext(&th, &tl,   zh,zl,  t,0);
   Add22_ext(&th, &tl,   c8h,c8l,  th,tl);
   Mul22_ext(&th, &tl,   zh,zl,  th,tl);
   Add22_ext(&th, &tl,   c7h,c7l,  th,tl);
   Mul22_ext(&th, &tl,   zh,zl,  th,tl);
   Add22_ext(&th, &tl,   c6h,c6l,  th,tl);
   Mul22_ext(&th, &tl,   zh,zl,  th,tl);
   Add22_ext(&th, &tl,   c5h,c5l,  th,tl);
   Mul22_ext(&th, &tl,   zh,zl,  th,tl);
   Add22_ext(&th, &tl,   c4h,c4l,  th,tl);
   Mul22_ext(&th, &tl,   zh,zl,  th,tl);
   Add22_ext(&th, &tl,   c3h,c3l,  th,tl);
   Mul22_ext(&th, &tl,   zh,zl,  th,tl);
   Add22_ext(&th, &tl,   c2h,c2l,  th,tl);
   Mul22_ext(&th, &tl,   zh,zl,  th,tl);
   Add22_ext(&th, &tl,   c1h,c1l,  th,tl);

   Mul22_ext(&th, &tl,   zh,zl,  th,tl);
   Add22_ext(&th, &tl,   argredtable[index].logwih, argredtable[index].logwil,  th,tl);
   
   Mul22_ext(&eh, &el,   ln2h,ln2l, E, 0);
   Add22_ext(&th, &tl,   eh,el,  th,tl);
   
   BACK_TO_DOUBLE_MODE;
   return (double) (th+tl); /* The exact sum of these double-extended is rounded to the nearest */
}


double log_ru(double x) { return x;};
double log_rd(double x) { return x;};
double log_rz(double x) { return x;};
