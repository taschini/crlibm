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




/* The Add22 and Add22 functions, as well as double-double
multiplications of the Dekker family may be either defined as
functions, or as #defines.  Which one is better depends on the
processor/compiler/OS.  As #define has to be used with more care (not
type-safe), the two following variables should  be set to 1 in the
development/debugging phase, until no type warning remains.  
Optimal values for Pentium 4:
#define ADD22_AS_FUNCTIONS 0
#define DEKKER_AS_FUNCTIONS 1

*/

#define ADD22_AS_FUNCTIONS 0
#define DEKKER_AS_FUNCTIONS 0



/* setting the following variable adds variables and code for
   monitoring the performance.
   Note that sometimes only round to nearest is instrumented */
#define EVAL_PERF  1


#if EVAL_PERF==1
/* counter of calls to the second step (accurate step) */
int crlibm_second_step_taken;
#endif





/*
 * i = d in rounding to nearest
 */
#define DOUBLE2INT(ii, dd)       \
  {db_number t;              \
   t.d=(dd+6755399441055744.0);  \
   ii=t.i[LO_ENDIAN];}








/* If the processor has a FMA, use it !   **/

#ifdef CRLIBM_TYPECPU_POWERPC
#define PROCESSOR_HAS_FMA 1
#define FMA(r, a,b,c)  /* r = a*b + c*/                \
do{                                                    \
  double _a, _b,_c,_r;                                 \
  _a=a; _b=b;_c=c;                                     \
  __asm__ __volatile__("fmadd %0, %1, %2, %3\n ;;\n"   \
		       : "=f"(_r)                      \
		       : "f"(_a), "f"(_b), "f"(_c)     \
		       );                              \
  r=_r;                                                \
} while(1+1==3)


#define FMS(r, a,b,c)   /* r = a*b - c*/               \
do{                                                    \
  double _a, _b,_c,_r;                                 \
  _a=a; _b=b;_c=c;                                     \
  __asm__ __volatile__("fmsub %0, %1, %2, %3\n ;;\n"   \
		       : "=f"(_r)                      \
		       : "f"(_a), "f"(_b), "f"(_c)     \
		       );                              \
  r=_r;                                                \
  } while(1+1==3)

#endif /*CRLIBM_TYPECPU_POWERPC*/




#ifdef CRLIBM_TYPECPU_ITANNIUM
#define PROCESSOR_HAS_FMA 1
#define FMA(r, a,b,c)\  /* r = a*b + c*/
do{                                                    \
  double _a, _b,_c,_r;                                 \
  _a=a; _b=b;_c=c;                                     \
  __asm__ __volatile__("fms %0 = %1, %2, %3\n ;;\n"    \
		       : "=f"(_r)                      \
		       : "f"(_a), "f"(_b), "f"(_c)     \
		       );                              \
  r=_r;                                                \
} while(1+1==3)


#define FMS(r, a,b,c)\  /* r = a*b - c*/
do{                                                    \
  double _a, _b, _c, _r;                               \
  _a=a; _b=b;_c=c;                                     \
  __asm__ __volatile__("fms %0 = %1, %2, %3\n ;;\n"    \
		       : "=f"(_r)                      \
		       : "f"(_a), "f"(_b), "f"(_c)     \
		       );                              \
  r=_r;                                                \
  } while(1+1==3)

#endif /*CRLIBM_TYPECPU_ITANIUM*/








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





#define ABS(x) (((x)>0) ? (x) : (-(x)))





/*
 * In the following, when an operator is preceded by a '@' it means that we
 * are considering the IEEE-compliant machine operator, otherwise it
 * is the mathematical operator.
 *
 */


/*
 * computes s and r such that s + r = a + b,  with s = a @+ b exactly 
 */
#define Add12Cond(s, r, a, b)     \
        {double z, _a=a, _b=b;    \
         s = _a + _b;             \
         if (ABS(a) > ABS(b)){    \
           z = s - _a;            \
           r = _b - z;            \
         }else {                  \
           z = s - _b;            \
           r = _a - z;}}                          

/*
 *  computes s and r such that s + r = a + b,  with s = a @+ b exactly 
 * under the condition  a >= b
 */
#define Add12(s, r, a, b)         \
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







/*
 * Functions to computes double-double addition: zh+zl = xh+xl + yh+yl
 * knowing that xh>yh
 * relative error is smaller than 2^-103 
 */


#if ADD22_AS_FUNCTIONS
extern void Add22(double *zh, double *zl, double xh, double xl, double yh, double yl);
extern void Add22Cond(double *zh, double *zl, double xh, double xl, double yh, double yl);

#else /* ADD22_AS_FUNCTIONS */

#define Add22Cond(zh,zl,xh,xl,yh,yl)                             \
do {                                                             \
  double r,s;                                                    \
  r = xh+yh;                                                     \
  s = ((ABS(xh)) > (ABS(yh)))? (xh-r+yh+yl+xl) : (yh-r+xh+xl+yl);\
  *zh = r+s;                                                     \
  *zl = r - (*zh) + s;                                           \
} while(2+2==5)

  
#define Add22(zh,zl,xh,xl,yh,yl)             \
do {                                         \
  double r,s;                                \
  r = xh+yh;                                 \
  s = xh-r+yh+yl+xl;                         \
  *zh = r+s;                                 \
  *zl = r - (*zh) + s;                       \
} while(2+2==5)
#endif /* ADD22_AS_FUNCTIONS */














#ifdef PROCESSOR_HAS_FMA
/* One of the nice things with the fused multiply-and-add is that it
   greatly simplifies the double-double multiplications : */
#define Mul12(rh,rl,u,v)                              \
{                                                     \
  *rh = u*v;                                          \
  FMS(*rl,   u,v, *rh);                               \
}


#define Mul22(pzh,pzl,xh,xl,yh,yl)                    \
{                                                     \
double mh, ml;                                        \
  mh = xh*yh;                                         \
  FMS(ml, xh, yh,  mh);                               \
  FMA(ml, xh,xl, ml);                                 \
  FMA(ml, xl,yh, ml);                                 \
  *pzh = mh+ml;					      \
  *pzl = mh - (*pzh) + ml;                            \
}

/* besides we don't care anymore about overflows in the mult  */
#define Mul12Cond Mul12    
#define Mul22cond Mul22


#else /* PROCESSOR_HAS_FMA */


#if DEKKER_AS_FUNCTIONS
extern void Mul12(double *rh, double *rl, double u, double v);
extern void Mul12Cond(double *rh, double *rl, double a, double b);
extern void Mul22(double *zh, double *zl, double xh, double xl, double yh, double yl);
#else /* if DEKKER_AS_FUNCTIONS  */
/*
 * computes rh and rl such that rh + rl = a * b with rh = a @* b exactly
 * under the conditions : a < 2^970 et b < 2^970 
 */
#define Mul12(rh,rl,u,v)                        \
{                                               \
  const double c  = 134217729.; /* 2^27 +1 */   \
  double up, u1, u2, vp, v1, v2;                \
  double _u =u, _v=v;                           \
                                                \
  up = _u*c;        vp = _v*c;                  \
  u1 = (_u-up)+up;  v1 = (_v-vp)+vp;            \
  u2 = _u-u1;       v2 = _v-v1;                 \
                                                \
  *rh = _u*_v;                                  \
  *rl = (((u1*v1-*rh)+(u1*v2))+(u2*v1))+(u2*v2);\
}


/*
 * Computes rh and rl such that rh + rl = a * b and rh = a @* b exactly
 */
#define Mul12Cond(rh, rl, a,  b) \
{\
  const double two_em53 = 1.1102230246251565404e-16; /* 0x3CA00000, 0x00000000 */\
  const double two_e53  = 9007199254740992.;         /* 0x43400000, 0x00000000 */\
  double u, v;                                            \
  db_number _a=a, _b=b;                               \
                                                          \
  if (_a.i[HI_ENDIAN]>0x7C900000) u = _a*two_em53;        \
  else            u = _a;                                 \
  if (_b.i[HI_ENDIAN]>0x7C900000) v = _b*two_em53;        \
  else            v = _b;                                 \
                                                          \
  Mul12(rh, rl, u, v);                                   \
                                                          \
  if (_a.i[HI_ENDIAN]>0x7C900000) {*rh *= two_e53; *rl *= two_e53;} \
  if (_b.i[HI_ENDIAN]>0x7C900000) {*rh *= two_e53; *rl *= two_e53;} \
}



/*
 * computes double-double multiplication: zh+zl = (xh+xl) *  (yh+yl)
 * relative error is smaller than 2^-102
 */
  

  
#define Mul22(zh,zl,xh,xl,yh,yl)                      \
{                                                     \
double mh, ml;                                        \
						      \
  const double c = 134217729.;			      \
  double up, u1, u2, vp, v1, v2;		      \
						      \
  up = (xh)*c;        vp = (yh)*c;			      \
  u1 = ((xh)-up)+up;  v1 = ((yh)-vp)+vp;		      \
  u2 = (xh)-u1;       v2 = (yh)-v1;                       \
  						      \
  mh = (xh)*(yh);					      \
  ml = (((u1*v1-mh)+(u1*v2))+(u2*v1))+(u2*v2);	      \
						      \
  ml += (xh)*(yl) + (xl)*(yh);				      \
  *zh = mh+ml;					      \
  *zl = mh - (*zh) + ml;                              \
}

#endif /* DEKKER_AS_FUNCTIONS */

#endif /* PROCESSOR_HAS_FMA */










/* A few prototypes that are not worth being in crlibm.h */

/*
 * Make a trigonometric range reduction with a relative
 * error less than 2^(-150) 
 */

extern int rem_pio2_scs(scs_ptr,  scs_ptr);



#endif /*CRLIBM_PRIVATE_H*/
