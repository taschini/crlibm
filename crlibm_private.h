/*
 *  crlibm_private.h
 *  
 * This file contains useful tools and data for the crlibm functions.
 *
 */



#ifndef CRLIBM_PRIVATE_H
#define CRLIBM_PRIVATE_H 1

#include "scs_lib/scs.h"
#include "scs_lib/scs_private.h"


 /* undef all the variables that might have been defined in
    scs_lib/scs_private.h */
#undef VERSION 
#undef PACKAGE 
#undef HAVE_GMP_H
#undef HAVE_MPFR_H
#undef HAVE_MATHLIB_H
/* then include the proper definitions  */
#include "crlibm_config.h"


/* setting the following variable adds variables and code for
   monitoring the performance.
   Note that only round to nearest is instrumented */
#define EVAL_PERF  1


#if EVAL_PERF==1
/* counter of calls to the second step (accurate step) */
int crlibm_second_step_taken;
#endif


#ifdef WORDS_BIGENDIAN
 #define DB_ONE    {{0x3ff00000, 0x00000000}}
#else
 #define DB_ONE    {{0x00000000 ,0x3ff00000}}
#endif


#ifdef WORDS_BIGENDIAN
#define HI(x) (*((int*)(&x)))
#define LO(x) (*(1+(int*)(&x)))
#else
#define HI(x) (*(1+(int*)(&x)))
#define LO(x) (*((int*)(&x)))
#endif


#ifdef SCS_TYPECPU_SPARC
static const scs
/* 0   */
   scs_zer ={{0x00000000, 0x00000000, 0x00000000, 0x00000000},
             {{0, 0}},  0,   1 },
/* 1/2 */
   scs_half={{0x02000000, 0x00000000, 0x00000000, 0x00000000},
             DB_ONE, -1,   1 },
/*  1  */  
   scs_one ={{0x00000001, 0x00000000, 0x00000000, 0x00000000},
             DB_ONE,  0,   1 },
/*  2  */
   scs_two ={{0x00000002, 0x00000000, 0x00000000, 0x00000000},
             DB_ONE,  0,   1 },	

/* ~1.666667e-01 */ 
   scs_sixinv ={{0x0aaaaaaa, 0x2aaaaaaa, 0x2aaaaaaa, 0x2aaaaaaa},
	     DB_ONE,  -1,   1 };

#else
static const struct scs
/* 0   */
   scs_zer ={{0x00000000, 0x00000000, 0x00000000, 0x00000000,
             0x00000000, 0x00000000, 0x00000000, 0x00000000},
             {{0, 0}},  0,   1 },
/* 1/2 */
   scs_half={{0x20000000, 0x00000000, 0x00000000, 0x00000000,
             0x00000000, 0x00000000, 0x00000000, 0x00000000},
             DB_ONE, -1,   1 },
/*  1  */  
   scs_one ={{0x00000001, 0x00000000, 0x00000000, 0x00000000,
             0x00000000, 0x00000000, 0x00000000, 0x00000000},
             DB_ONE,  0,   1 },
/*  2  */
   scs_two ={{0x00000002, 0x00000000, 0x00000000, 0x00000000,
             0x00000000, 0x00000000, 0x00000000, 0x00000000},
             DB_ONE,  0,   1 },
/* 0.166666*/
   scs_sixinv ={{0x0aaaaaaa, 0x2aaaaaaa, 0x2aaaaaaa, 0x2aaaaaaa, 
	     0x2aaaaaaa, 0x2aaaaaaa, 0x2aaaaaaa, 0x2aaaaaaa},
	     DB_ONE,  -1,   1 };

#endif

#define SCS_ZERO    (scs_ptr)(&scs_zer)
#define SCS_HALF    (scs_ptr)(&scs_half)
#define SCS_ONE     (scs_ptr)(&scs_one)
#define SCS_TWO     (scs_ptr)(&scs_two)
#define SCS_SIXINV  (scs_ptr)(&scs_sixinv)

/*
 * Define rounding mode
 */
#define RNDN 0  /* to nearest  */
#define RNDD 1  /* toward -inf */
#define RNDU 2  /* toward +inf */







/* This sets round to the nearest and disables extended precision on
   the x86s. Should be done for the Itanii too */
#ifdef SCS_TYPECPU_X86
#include <fpu_control.h>
#ifndef __setfpucw
#define __setfpucw(cw) __asm__ ("fldcw %0" : : "m" (cw))
#endif 
#endif /*SCS_TYPECPU_X86*/










/*
 * In the following, when an operator is preceded by a '@' it means that we
 * are considering the IEEE-compliant machine operator, otherwise it
 * is the mathematical operator.
 *
 */


/*
 * computes s and r such that s + r = a + b,  with s = a @+ b exactly 
 */
#define Add12Cond(s, r, a, b)   \
        {double z;                 \
         db_number _a, _b;     \
          _a.d = a; _b.d = b;      \
         s = _a.d + _b.d;          \
         if ((_a.i[HI_ENDIAN]&0x7fffffff) > (_b.i[HI_ENDIAN]&0x7fffffff)){  \
           z = s - _a.d;           \
           r = _b.d - z;           \
         }else {                   \
           z = s - _b.d;           \
           r = _a.d - z;}}                          

/*
 *  computes s and r such that s + r = a + b,  with s = a @+ b exactly 
 * under the condition  a >= b
 */
#define Add12(s, r, a, b)      \
        {double z, _a=a, _b=b;    \
         s = _a + _b;             \
         z = s - _a;              \
         r = _b - z; }            


/*
 * computes r1, r2, r3 such that r1 + r2 + r3 = a + b + c exactly 
 */
#define Fast3Sum(r1, r2, r3, a, b, c) \
        {double u, v, w;              \
         Fast2Sum(u, v, b, c);        \
         Fast2Sum(r1, w, a, u);       \
         Fast2Sum(r2, r3, w, v); }





extern void Mul22(double * z, double * zz, double x, double xx, double y, double yy);
extern void Add22(double * z, double * zz, double x, double xx, double y, double yy);

/* Procedures like Dekker may be either defined as inline functions, or as #defines. 
Currently #define is better (smaller and faster) 
but it has to be used with more care  */

#define DEKKER_INLINED 1

#if DEKKER_INLINED
extern void Mul12(double *r1, double *r2, double u, double v);
extern void Mul12Cond(double *r1, double *r2, double a, double b);
#else


/*
 * computes r1 and r2 such that r1 + r2 = a * b with r1 = a @* b exactly
 * under the conditions : a < 2^970 et b < 2^970 
 */
#ifdef CRLIBM_TYPECPU_ITANIUM
/* One of the nice things with the fused multiply-and-add is that it
   greatly simplifies the Dekker : */
#define Dekker(r1,r2,u,v)                            \
{                                                     \
  *r1 = u*v;                                          \
  /* The following means: *r2 = FMS(u*v-*r1) */       \
  __asm__ __volatile__("fms %0 = %1, %2, %3\n ;;\n"   \
		       : "=f"(*r2)                     \
		       : "f"(u), "f"(v), "f"(*r1)     \
		       );                             \
}
#else
#define Dekker(r1,r2,u,v)                            \
{                                               \
  double c        = 134217729.;               /* 0x41A00000, 0x02000000 */ \
  double up, u1, u2, vp, v1, v2;                \
  double _u =u, _v=v;                           \
                                                \
  up = _u*c;        vp = _v*c;                  \
  u1 = (_u-up)+up;  v1 = (_v-vp)+vp;            \
  u2 = _u-u1;       v2 = _v-v1;                 \
                                                \
  *r1 = _u*_v;                                  \
  *r2 = (((u1*v1-*r1)+(u1*v2))+(u2*v1))+(u2*v2);\
}
#endif /*CRLIBM_TYPECPU_ITANIUM*/


/*
 * Computes r1 and r2 such that r1 + r2 = a * b and r1 = a @* b exactly
 */
#define DekkerCond(r1, r2, a,  b) \
{\
  double two_em53 = 1.1102230246251565404e-16; /* 0x3CA00000, 0x00000000 */\
  double two_e53  = 9007199254740992.;         /* 0x43400000, 0x00000000 */\
  double u, v;                                            \
  db_number _a=a, _b=b;                               \
                                                          \
  if (_a.i[HI_ENDIAN]>0x7C900000) u = _a*two_em53;        \
  else            u = _a;                                 \
  if (_b.i[HI_ENDIAN]>0x7C900000) v = _b*two_em53;        \
  else            v = _b;                                 \
                                                          \
  Dekker(r1, r2, u, v);                                   \
                                                          \
  if (_a.i[HI_ENDIAN]>0x7C900000) {*r1 *= two_e53; *r2 *= two_e53;} \
  if (_b.i[HI_ENDIAN]>0x7C900000) {*r1 *= two_e53; *r2 *= two_e53;} \
}


#endif /*DEKKERINLINED*/



/*
 * i = d in rounding to nearest
 */
#define DOUBLE2INT(ii, dd)       \
  {db_number t;              \
   t.d=(dd+6755399441055744.0);  \
   ii=t.i[LO_ENDIAN];}




/* A few prototypes that are not worth being in crlibm.h */

/*
 * Make a trigonometric range reduction with a relative
 * error less than 2^(-150) 
 */

extern int rem_pio2_scs(scs_ptr,  scs_ptr);



#endif /*CRLIBM_PRIVATE_H*/
