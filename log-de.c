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

on Itanium:  
icc  -D__ICC__ -I/users/fdedinex/local/IA64/include -mcpu=itanium2\
 -Qoption,cpp,--extended_float_types -IPF_fp_speculationsafe -c log-de.c;\
 mv log-de.o log_fast.o; make

 
*/


#include <stdio.h>
#include <stdlib.h>
#include "crlibm.h"
#include "crlibm_private.h"
#include "double-extended.h"
#include "log-de.h"


double log_rn(double x) {
  double wi;
#if defined(CRLIBM_TYPECPU_X86) || defined(CRLIBM_TYPECPU_AMD64)
  db_number xdb;
  long double logwih, logwil, invwih, invwil, y, z, z0, t, p1,p2,p3,z2,z4, zh, zl, th, tl, eh,el;
  int E, i, index, roundtestmask;
#else
  long int E, i;
  unsigned long int index, roundtestmask;
  __fpreg logwih, logwil, invwih, invwil, y, z, z0, t, p1,p2,p3,z2,z4, zh, zl, th, tl, eh,el;
  __fpreg c1,c2,c3,c4,c5,c6,c7,c8;
  __fpreg c1h,c2h,c3h,c4h,c5h,c6h,c7h,c8h,c9h,c10h,c11h,c12h,c13h,c14h,c15h;
  __fpreg c1l,c2l,c3l,c4l,c5l,c6l,c7l,c8l;
#endif

   E=0;

#if defined(CRLIBM_TYPECPU_X86) || defined(CRLIBM_TYPECPU_AMD64)

   xdb.d=x;

   /* Filter cases */
   if (xdb.i[HI] < 0x00100000){        /* x < 2^(-1022)    */
     if (((xdb.i[HI] & 0x7fffffff)|xdb.i[LO])==0)    return -1.0/0.0;  /* log(+/-0) = -Inf */
     if (xdb.i[HI] < 0)                              return (x-x)/0;   /* log(-x) = Nan    */
     /* Else subnormal number */
     E = -64; 		
     xdb.d *= two64; 	  /* make x a normal number    */ 
   }
    
   if (xdb.i[HI] >= 0x7ff00000)                      return  x+x;      /* Inf or Nan       */
   
   DOUBLE_EXTENDED_MODE;  /* This one should be overlapped with following integer computation */

   /* Extract exponent and mantissa */
   E += (xdb.i[HI]>>20)-1023;             /* extract the exponent */
   index = (xdb.i[HI] & 0x000fffff);
   xdb.i[HI] =  index | 0x3ff00000;	/* do exponent = 0 */
   index = (index + (1<<(20-L-1))) >> (20-L); 
   /* reduce  such that sqrt(2)/2 < xdb.d < sqrt(2) */
   if (index >= MAXINDEX){ /* corresponds to y>sqrt(2)*/
     xdb.i[HI] -= 0x00100000; 
     E++;
   }
   y = xdb.d;
   /* now y holds 1+f, and E is the exponent */


#else  /*  Here come the code specific to Itanium processor */
   /* prefetch coefficients */
   PREFETCH_POLY_QUICK; /* defined in log-de.h */
 
   y=x;
   i =  _Asm_getf(2/*_FR_D*/, y);  /* Cast y to a 64-bit integer */

   /* Filter special cases */
   if (i<0x0010000000000000LL){   /* equivalent to : x < 2^(-1022)    */
     if ((i & 0x7fffffffffffffffULL)==0)  return -1.0/0.0;    /* log(+/-0) = -Inf */
     if (i<0)                             return (x-x)/0;     /* log(-x) = Nan    */
     /* Else subnormal number */
     y *= two64; 	  /* make x a normal number    */ 
     E = -64;
     i =  _Asm_getf(2/*_FR_D*/, y); /* and update i */ 
   }
    
   if (i >= 0x7ff0000000000000ULL)        return  x+x;	      /* Inf or Nan       */

   /* Extract exponent and mantissa */
     E += (i>>52)-1023;
     i = i & 0x000fffffffffffffULL;  /* keep only mantissa */
     index = (i + (1ULL<<(52-L-1))) >> (52-L);
     /* reduce  such that sqrt(2)/2 < xdb.d < sqrt(2) */
     if (index >= MAXINDEX){    /* corresponds to y>sqrt(2)*/
       y = _Asm_setf(2/*_FR_D*/, (i | 0x3ff0000000000000ULL) - 0x0010000000000000ULL ); /* exponent = -1 */
       E++;
     }
     else {
       y = _Asm_setf(2/*_FR_D*/, i | 0x3ff0000000000000ULL ); /* exponent = 0*/
     }

     /* now y holds 1+f, and E is the exponent */
   
#endif

   /*read the tables */
   wi     = argredtable[index].wi;
   invwih = argredtable[index].invwih;
   logwih    = argredtable[index].logwih;
 
   z0 = y-wi; /* exact thanks to Sterbenz */

   /* Polynomial evaluation, using Estrin scheme */


   roundtestmask=0x7fe;
   z  = z0*invwih;
   z2 = z*z;             p1 = c2+c3*z;    p2 = c4+c5*z;     p3 = c6+c7*z;
   t = logwih + c1*z;    z4 = z2*z2;      p1 = p1+p2*z2;    p3 = p3+c8*z2;
   p1 = p1 + z4*p3;      t = t + E*ln2h;
   t = t + z2 *p1;   


#if 0 /* to time the first step only */
   BACK_TO_DOUBLE_MODE; return (double)t;
#endif

#if 1
   TEST_AND_RETURN_RN(t, roundtestmask);
#else
   TEST_AND_RETURN_RN_ZIV(t,1.001);
#endif


   /* Accurate phase */ 
#if EVAL_PERF
  crlibm_second_step_taken++;
#endif

   Mul22_ext(&zh,&zl, argredtable[index].invwih, argredtable[index].invwil,  z0, 0);
   
   PREFETCH_POLY_ACCURATE;
   t = c15h;
   t = c14h + zh*t;
   t = c13h + zh*t;
   t = c12h + zh*t;
   t = c11h + zh*t;
   t = c10h + zh*t;
   t = c9h  + zh*t;
   FMA22_ext(&th, &tl,   zh,zl,  t,0,     c8h,c8l);
   FMA22_ext(&th, &tl,   zh,zl,  th,tl,   c7h,c7l);
   FMA22_ext(&th, &tl,   zh,zl,  th,tl,   c6h,c6l);
   FMA22_ext(&th, &tl,   zh,zl,  th,tl,   c5h,c5l);
   FMA22_ext(&th, &tl,   zh,zl,  th,tl,   c4h,c4l);
   FMA22_ext(&th, &tl,   zh,zl,  th,tl,   c3h,c3l);
   FMA22_ext(&th, &tl,   zh,zl,  th,tl,   c2h,c2l);
   FMA22_ext(&th, &tl,   zh,zl,  th,tl,   c1h,c1l);
   FMA22_ext(&th, &tl,   zh,zl,  th,tl,   argredtable[index].logwih, argredtable[index].logwil);
   
   Mul22_ext(&eh, &el,   ln2h,ln2l, E, 0);
   Add22_ext(&th, &tl,   eh,el,  th,tl);
   
   BACK_TO_DOUBLE_MODE;
   return (double) (th+tl); /* The exact sum of these double-extended is rounded to the nearest */
}


double log_ru(double x) { return x;};
double log_rd(double x) { return x;};
double log_rz(double x) { return x;};
