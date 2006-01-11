/* 
 *this function computes log, correctly rounded, 
 using  double-extended arithmetic

 THIS IS EXPERIMENTAL SOFTWARE

In particular it changes rounding modes all the time without warning
nor restoring.
 
 *
 * Author :  Florent de Dinechin
 * Florent.de.Dinechin at ens-lyon.fr
 *


This function compiles both on IA32 and IA64 architectures. On IA64,
it needs icc 8.1 or higher, with the following flags (which should be
set up by the autoconf).

icc -DHAVE_CONFIG_H  -Qoption,cpp,--extended_float_types \
                    -IPF_fp_speculationsafe -c log-de.c;\
 mv log-de.o log-td.o; make

 
*/


#include <stdio.h>
#include <stdlib.h>
#include "crlibm.h"
#include "crlibm_private.h"
#include "double-extended.h"
#include "log-de.h"



static void log_accurate(double_ext* prh, double_ext* prl, double_ext z, int E, int index) {

double_ext  th, tl, eh,el, t;

#if !(defined(CRLIBM_TYPECPU_X86) || defined(CRLIBM_TYPECPU_AMD64))
  double_ext c1h,c2h,c3h,c4h,c5h,c6h,c7h,c8h,c9h,c10h,c11h,c12h,c13h,c14h,c15h;
  double_ext c1l,c2l,c3l,c4l,c5l,c6l,c7l,c8l;
#endif


#if EVAL_PERF
  crlibm_second_step_taken++;
#endif
  
  /* TODO check the conditions for the double-double ops */

 
  PREFETCH_POLY_ACCURATE;
  t = c14h;
  t = c13h + z*t;
  t = c12h + z*t;
  t = c11h + z*t;
  t = c10h + z*t;
  t = c9h  + z*t;
  t = c8h  + z*t;
  Mul12_ext(&th, &tl,   z, t);
  Add22_ext(&th, &tl,   th,tl, c7h,c7l);
  FMA22_ext(&th, &tl,   z,0,   t,0,      c7h,c7l);
  FMA22_ext(&th, &tl,   z,0,   th,tl,    c6h,c6l);
  FMA22_ext(&th, &tl,   z,0,   th,tl,    c5h,c5l);
  FMA22_ext(&th, &tl,   z,0,   th,tl,    c4h,c4l);
  FMA22_ext(&th, &tl,   z,0,   th,tl,    c3h,c3l);
  FMA22_ext(&th, &tl,   z,0,   th,tl,    c2h,c2l);
  FMA22_ext(&th, &tl,   z,0,   th,tl,    c1h,c1l);
  FMA22_ext(&th, &tl,   z,0,   th,tl,    argredtable[index].logirh, argredtable[index].logirl);
  
  Mul22_ext(&eh, &el,   log2h,log2l, E, 0);
  Add22_ext(prh, prl,   eh,el,  th,tl);
}








double log_rn(double x) {
  double_ext logirh, r, y, z, z2,z4, th, tl, p01, p23, p45, p67, p03, p47, p07, logde;
#if defined(CRLIBM_TYPECPU_X86) || defined(CRLIBM_TYPECPU_AMD64)
  db_number xdb;
  int E, index, index0, roundtestmask;
#else
  int64_t  E, i;
  uint64_t index, roundtestmask;
  double_ext c1,c2,c3,c4,c5,c6,c7;
#endif


#if defined(CRLIBM_TYPECPU_X86) || defined(CRLIBM_TYPECPU_AMD64)
   xdb.d=x;

   index0 = (xdb.i[HI] & 0x000fffff);
   index = (index0 + (1<<(20-L-1))) >> (20-L); 
   E = (xdb.i[HI]>>20)-1023;             /* extract the exponent */

   /* Filter cases */
   if (xdb.i[HI] < 0x00100000){        /* x < 2^(-1022)    */
     if (((xdb.i[HI] & 0x7fffffff)|xdb.i[LO])==0)    return -1.0/0.0;  /* log(+/-0) = -Inf */
     if (xdb.i[HI] < 0)                              return (x-x)/0;   /* log(-x) = Nan    */
     /* Else subnormal number */
     xdb.d *= two64; 	  /* make x a normal number    */ 
     E = -64 + (xdb.i[HI]>>20)-1023;             /* extract the exponent */
     index0 = (xdb.i[HI] & 0x000fffff);
     index = (index0 + (1<<(20-L-1))) >> (20-L); 
   }
   if (xdb.i[HI] >= 0x7ff00000)                      return  x+x;      /* Inf or Nan       */
   
   DOUBLE_EXTENDED_MODE;  /* This one should be overlapped with following integer computation */

   /* Extract exponent and mantissa */
   xdb.i[HI] =  index0 | 0x3ff00000;	/* do exponent = 0 */
   /* reduce  such that sqrt(2)/2 < xdb.d < sqrt(2) */
   if (index >= MAXINDEX){ /* corresponds to y>sqrt(2)*/
     xdb.i[HI] -= 0x00100000; 
     index = index & INDEXMASK;
     E++;
}
   y = xdb.d;

#else /* defined(CRLIBM_TYPECPU_X86) || defined(CRLIBM_TYPECPU_AMD64) */
  /*  Here come the code specific to Itanium processor */
   E=0;
   PREFETCH_POLY_QUICK; /* defined in log-de.h */
   y=x;
   i =  _Asm_getf(2/*_FR_D*/, y);  /* Cast y to a 64-bit integer */

   /* Filter special cases */
   if (i<(int64_t)ULL(0010000000000000)){   /* equivalent to : x < 2^(-1022)    */
     if ((i & ULL(7fffffffffffffff))==0)  return -1.0/0.0;    /* log(+/-0) = -Inf */
     if (i<0)                             return (x-x)/0;     /* log(-x) = Nan    */
     /* Else subnormal number */
     y *= two64; 	  /* make x a normal number    */ 
     E = -64;
     i =  _Asm_getf(2/*_FR_D*/, y); /* and update i */ 
   }
   if (i >= ULL(7ff0000000000000))        return  x+x;	      /* Inf or Nan       */

   /* Extract exponent and mantissa */
   E += (i>>52)-1023;
   i = i & ULL(000fffffffffffff);  /* keep only mantissa */
   index = (i + (ULL(1)<<(52-L-1))) >> (52-L);
   /* reduce  such that sqrt(2)/2 < xdb.d < sqrt(2) */
   if (index >= MAXINDEX){    /* corresponds to y>sqrt(2)*/
     y = _Asm_setf(2/*_FR_D*/, (i | ULL(3ff0000000000000)) - ULL(0010000000000000) ); /* exponent = -1 */
     index = index & INDEXMASK;
     E++;
   }
   else 
     y = _Asm_setf(2/*_FR_D*/, i | ULL(3ff0000000000000) ); /* exponent = 0*/   
#endif  /* defined(CRLIBM_TYPECPU_X86) || defined(CRLIBM_TYPECPU_AMD64) */
 


   /* All the previous argument reduction was exact */
   /* now y holds 1+f, and E is the exponent */

     r = (double_ext) (argredtable[index].r); /* approx to 1/y.d */
     logirh = argredtable[index].logirh;
     z = y*r - 1. ; /* even without an FMA, all exact */

   if(E==0)
     roundtestmask=ACCURATE_TO_61_BITS;
    else
      roundtestmask=ACCURATE_TO_62_BITS;
      
   /* Estrin polynomial evaluation  */
   z2 = z*z;    p67 = c6 + z*c7;       p45 = c4 + z*c5;      p23 = c2 + z*c3;    p01 = logirh + z;
   z4 = z2*z2;  p47 = p45 + z2*p67;    p03 = p01 + z2*p23; 
   p07 = p03 + z4*p47;
   logde = p07 + E*log2h;
#if 0 /* to time the first step only */
   BACK_TO_DOUBLE_MODE; return (double)t;
#endif


   /* To test the second step only, comment out the following line */
   DE_TEST_AND_RETURN_RN(logde, roundtestmask);
   

   log_accurate(&th, &tl, z, E, index);
   
   BACK_TO_DOUBLE_MODE;
   
   return (double) (th+tl); /* The exact sum of these double-extended is rounded to the nearest */
}










double log_rd(double x) {
  double_ext logirh, r, y, z, z2,z4, th, tl, p01, p23, p45, p67, p03, p47, p07, logde;
#if defined(CRLIBM_TYPECPU_X86) || defined(CRLIBM_TYPECPU_AMD64)
  db_number xdb;
  int E, index, roundtestmask;
#else
  int64_t  E, i;
  uint64_t index, roundtestmask;
  double_ext c1,c2,c3,c4,c5,c6,c7;
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

#else /* defined(CRLIBM_TYPECPU_X86) || defined(CRLIBM_TYPECPU_AMD64) */
  /*  Here come the code specific to Itanium processor */
   PREFETCH_POLY_QUICK; /* defined in log-de.h */
   y=x;
   i =  _Asm_getf(2/*_FR_D*/, y);  /* Cast y to a 64-bit integer */

   /* Filter special cases */
   if (i<(int64_t)ULL(0010000000000000)){   /* equivalent to : x < 2^(-1022)    */
     if ((i & ULL(7fffffffffffffff))==0)  return -1.0/0.0;    /* log(+/-0) = -Inf */
     if (i<0)                             return (x-x)/0;     /* log(-x) = Nan    */
     /* Else subnormal number */
     y *= two64; 	  /* make x a normal number    */ 
     E = -64;
     i =  _Asm_getf(2/*_FR_D*/, y); /* and update i */ 
   }
   if (i >= ULL(7ff0000000000000))        return  x+x;	      /* Inf or Nan       */

   /* Extract exponent and mantissa */
   E += (i>>52)-1023;
   i = i & ULL(000fffffffffffff);  /* keep only mantissa */
   index = (i + (ULL(1)<<(52-L-1))) >> (52-L);
   /* reduce  such that sqrt(2)/2 < xdb.d < sqrt(2) */
   if (index >= MAXINDEX){    /* corresponds to y>sqrt(2)*/
     y = _Asm_setf(2/*_FR_D*/, (i | ULL(3ff0000000000000)) - ULL(0010000000000000) ); /* exponent = -1 */
     E++;
   }
   else 
     y = _Asm_setf(2/*_FR_D*/, i | ULL(3ff0000000000000) ); /* exponent = 0*/   
#endif  /* defined(CRLIBM_TYPECPU_X86) || defined(CRLIBM_TYPECPU_AMD64) */
 
   /* All the previous argument reduction was exact */
   /* now y holds 1+f, and E is the exponent */
   index = index & INDEXMASK;

   logirh = argredtable[index].logirh;
   r = (double_ext) (argredtable[index].r); /* approx to 1/y.d */
   z = y*r - 1. ; /* even without an FMA, all exact */
   
   if(E==0)
     roundtestmask=ACCURATE_TO_61_BITS;
   else
     roundtestmask=ACCURATE_TO_62_BITS;

   /* Estrin polynomial evaluation  */
   z2 = z*z;    p67 = c6 + z*c7;       p45 = c4 + z*c5;      p23 = c2 + z*c3;    p01 = logirh + z;
   z4 = z2*z2;  p47 = p45 + z2*p67;    p03 = p01 + z2*p23; 
   p07 = p03 + z4*p47;
   logde = p07 + E*log2h;
#if 0 /* to time the first step only */
   BACK_TO_DOUBLE_MODE; return (double)t;
#endif


   /* To test the second step only, comment out the following line */
   DE_TEST_AND_RETURN_RD(logde, roundtestmask);

   log_accurate(&th, &tl, z, E, index);

   RETURN_SUM_ROUNDED_DOWN(th, tl);

}










double log_ru(double x) {
  double_ext logirh, r, y, z, z2,z4, th, tl, p01, p23, p45, p67, p03, p47, p07, logde;
#if defined(CRLIBM_TYPECPU_X86) || defined(CRLIBM_TYPECPU_AMD64)
  db_number xdb;
  int E, index, roundtestmask;
#else
  int64_t  E, i;
  uint64_t index, roundtestmask;
  double_ext c1,c2,c3,c4,c5,c6,c7;
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
     index = index & INDEXMASK;
     xdb.i[HI] -= 0x00100000; 
     E++;
   }
   y = xdb.d;

#else /* defined(CRLIBM_TYPECPU_X86) || defined(CRLIBM_TYPECPU_AMD64) */
  /*  Here come the code specific to Itanium processor */
   PREFETCH_POLY_QUICK; /* defined in log-de.h */
   y=x;
   i =  _Asm_getf(2/*_FR_D*/, y);  /* Cast y to a 64-bit integer */

   /* Filter special cases */
   if (i<(int64_t)ULL(0010000000000000)){   /* equivalent to : x < 2^(-1022)    */
     if ((i & ULL(7fffffffffffffff))==0)  return -1.0/0.0;    /* log(+/-0) = -Inf */
     if (i<0)                             return (x-x)/0;     /* log(-x) = Nan    */
     /* Else subnormal number */
     y *= two64; 	  /* make x a normal number    */ 
     E = -64;
     i =  _Asm_getf(2/*_FR_D*/, y); /* and update i */ 
   }
   if (i >= ULL(7ff0000000000000))        return  x+x;	      /* Inf or Nan       */

   /* Extract exponent and mantissa */
   E += (i>>52)-1023;
   i = i & ULL(000fffffffffffff);  /* keep only mantissa */
   index = (i + (ULL(1)<<(52-L-1))) >> (52-L);
   /* reduce  such that sqrt(2)/2 < xdb.d < sqrt(2) */
   if (index >= MAXINDEX){    /* corresponds to y>sqrt(2)*/
     y = _Asm_setf(2/*_FR_D*/, (i | ULL(3ff0000000000000)) - ULL(0010000000000000) ); /* exponent = -1 */
     index = index & INDEXMASK;
     E++;
   }
   else 
     y = _Asm_setf(2/*_FR_D*/, i | ULL(3ff0000000000000) ); /* exponent = 0*/   
#endif  /* defined(CRLIBM_TYPECPU_X86) || defined(CRLIBM_TYPECPU_AMD64) */
 
   /* All the previous argument reduction was exact */
   /* now y holds 1+f, and E is the exponent */
   
   logirh = argredtable[index].logirh;
   r = (double_ext) (argredtable[index].r); /* approx to 1/y.d */
   z = y*r - 1. ; /* even without an FMA, all exact */
   
   if(E==0)
     roundtestmask=ACCURATE_TO_61_BITS;
   else
     roundtestmask=ACCURATE_TO_62_BITS;
   
   /* Estrin polynomial evaluation  */
   z2 = z*z;    p67 = c6 + z*c7;       p45 = c4 + z*c5;      p23 = c2 + z*c3;    p01 = logirh + z;
   z4 = z2*z2;  p47 = p45 + z2*p67;    p03 = p01 + z2*p23; 
   p07 = p03 + z4*p47;
   logde = p07 + E*log2h;
#if 0 /* to time the first step only */
   BACK_TO_DOUBLE_MODE; return (double)t;
#endif

   
   /* To test the second step only, comment out the following line */
   DE_TEST_AND_RETURN_RU(logde, roundtestmask);

   log_accurate(&th, &tl, z, E, index);

   RETURN_SUM_ROUNDED_UP(th, tl);

}


double log_rz(double x) {
  if (x>1.0)
    return log_rd(x);
  else
    return log_ru(x);
}


