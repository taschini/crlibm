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


#if defined(CRLIBM_TYPECPU_X86) || defined(CRLIBM_TYPECPU_AMD64) 
#include <fpu_control.h>
#endif /* defined(CRLIBM_TYPECPU_X86) || defined(CRLIBM_TYPECPU_AMD64) */



/* The Add22 and Add22 functions, as well as double-double
multiplications of the Dekker family may be either defined as
functions, or as #defines.  Which one is better depends on the
processor/compiler/OS.  As #define has to be used with more care (not
type-safe), the two following variables should  be set to 1 in the
development/debugging phase, until no type warning remains.  

*/

#define ADD22_AS_FUNCTIONS 0
#define DEKKER_AS_FUNCTIONS 0



/* setting the following variable adds variables and code for
   monitoring the performance.
   Note that sometimes only round to nearest is instrumented */
#define EVAL_PERF  1


#if EVAL_PERF==1
/* counter of calls to the second step (accurate step) */
extern int crlibm_second_step_taken;
#endif



/* The prototypes of the second steps */
extern void exp_SC(scs_ptr res_scs, double x);
extern double scs_exp_rn(double);
extern double scs_exp_ru(double);
extern double scs_exp_rd(double);

extern void scs_log(scs_ptr,db_number, int);

extern double scs_atan_rn(double); 
extern double scs_atan_rd(double); 
extern double scs_atan_ru(double); 
 
extern double scs_sin_rn(double);
extern double scs_sin_ru(double);
extern double scs_sin_rd(double);
extern double scs_sin_rz(double);
extern double scs_cos_rn(double);
extern double scs_cos_ru(double);
extern double scs_cos_rd(double);
extern double scs_cos_rz(double);
extern double scs_tan_rn(double); 
extern double scs_tan_rd(double);  
extern double scs_tan_ru(double);  
extern double scs_tan_rz(double);  

extern int rem_pio2_scs(scs_ptr, scs_ptr);


/*
 * i = d in rounding to nearest
  The constant added is 2^52 + 2^51 
 */
#define DOUBLE2INT(_i, _d)       \
  {db_number _t;              \
   _t.d = (_d+6755399441055744.0);  \
   _i = _t.i[LO];}


/* Same idea but beware: works only for |_i| < 2^51 -1 */
#define DOUBLE2LONGINT(_i, _d)                                      \
  {                                                                 \
    db_number _t;                                                   \
    _t.d = (_d+6755399441055744.0);                                 \
    if (_d >= 0) /* sign extend */                                  \
      _i = _t.l & 0x0007FFFFFFFFFFFFLL;                             \
    else                                                            \
      _i = (_t.l & 0x0007FFFFFFFFFFFFLL) |  (0xFFF8000000000000LL); \
  }





/* Macros for the rounding tests in directed modes */
/* After Evgeny Gvozdev pointed out a bug in the rounding procedures I
   decided to centralize them here 

Note that these tests launch the accurate phase when yl=0, in
particular in the exceptional cases when the image of a double is a
double. See the chapter about the log for an example

*/


#define TEST_AND_RETURN_RU(__yh__, __yl__, __eps__)                    \
{                                                                      \
  db_number yh, yl, u53;  int yh_neg, yl_neg;                          \
  yh.d = __yh__;    yl.d = __yl__;                                     \
  yh_neg = (yh.i[HI] & 0x80000000);                                    \
  yl_neg = (yl.i[HI] & 0x80000000);                                    \
  yh.l = yh.l & 0x7fffffffffffffffLL;  /* compute the absolute value*/ \
  yl.l = yl.l & 0x7fffffffffffffffLL;  /* compute the absolute value*/ \
  u53.l     = (yh.l & 0x7ff0000000000000LL) +  0x0010000000000000LL;   \
  if(yl.d > __eps__ * u53.d){                                          \
    if(!yl_neg) {  /* The case yl==0 is filtered by the above test*/   \
      /* return next up */                                             \
      yh.d = __yh__;                                                   \
      if(yh_neg) yh.l--;  else yh.l++; /* Beware: fails for zero */    \
      return yh.d ;                                                    \
    }                                                                  \
    else  return __yh__;                                               \
  }                                                                    \
}


#define TEST_AND_RETURN_RD(__yh__, __yl__, __eps__)                    \
{                                                                      \
  db_number yh, yl, u53;  int yh_neg, yl_neg;                          \
  yh.d = __yh__;    yl.d = __yl__;                                     \
  yh_neg = (yh.i[HI] & 0x80000000);                                    \
  yl_neg = (yl.i[HI] & 0x80000000);                                    \
  yh.l = yh.l & 0x7fffffffffffffffLL;  /* compute the absolute value*/ \
  yl.l = yl.l & 0x7fffffffffffffffLL;  /* compute the absolute value*/ \
  u53.l     = (yh.l & 0x7ff0000000000000LL) +  0x0010000000000000LL;   \
  if(yl.d > __eps__ * u53.d){                                          \
    if(yl_neg) {   /* The case yl==0 is filtered by the above test*/   \
      /* return next down */                                           \
      yh.d = __yh__;                                                   \
      if(yh_neg) yh.l++;  else yh.l--; /* Beware: fails for zero */    \
      return yh.d ;                                                    \
    }                                                                  \
    else  return __yh__;                                               \
  }                                                                    \
}



#define TEST_AND_RETURN_RZ(__yh__, __yl__, __eps__)                    \
{                                                                      \
  db_number yh, yl, u53;  int yh_neg, yl_neg;                          \
  yh.d = __yh__;    yl.d = __yl__;                                     \
  yh_neg = (yh.i[HI] & 0x80000000);                                    \
  yl_neg = (yl.i[HI] & 0x80000000);                                    \
  yh.l = yh.l & 0x7fffffffffffffffLL;  /* compute the absolute value*/ \
  yl.l = yl.l & 0x7fffffffffffffffLL;  /* compute the absolute value*/ \
  u53.l     = (yh.l & 0x7ff0000000000000LL) +  0x0010000000000000LL;   \
  if(yl.d > __eps__ * u53.d){                                          \
    if(yl_neg!=yh_neg) {                                               \
      yh.d = __yh__;                                                   \
      yh.l--;                          /* Beware: fails for zero */    \
      return yh.d ;                                                    \
    }                                                                  \
    else  return __yh__;                                               \
  }                                                                    \
}


/* If the processor has a FMA, use it !   **/

/* All this probably works only with gcc. 
   See Markstein book for the case of HP's compiler */

#if defined(CRLIBM_TYPECPU_POWERPC) && defined(__GCC__)
#define PROCESSOR_HAS_FMA 1
#define FMA(a,b,c)  /* r = a*b + c*/                   \
({                                                     \
  double _a, _b,_c,_r;                                 \
  _a=a; _b=b;_c=c;                                     \
  __asm__ ("fmadd %0, %1, %2, %3\n ;;\n"               \
		       : "=f"(_r)                      \
		       : "f"(_a), "f"(_b), "f"(_c)     \
		       );                              \
 _r;                                                   \
})


#define FMS(a,b,c)   /* r = a*b - c*/                 \
({                                                    \
  double _a, _b,_c,_r;                                \
  _a=a; _b=b;_c=c;                                    \
  __asm__ ("fmsub %0, %1, %2, %3\n ;;\n"              \
		       : "=f"(_r)                     \
		       : "f"(_a), "f"(_b), "f"(_c)    \
		       );                             \
  _r;                                                 \
  })

#endif /* defined(CRLIBM_TYPECPU_POWERPC) && defined(__GCC__) */




/* On the Itanium 1 / gcc we lose 10 cycles when using the FMA !?! 
   It probably breaks the scheduling algorithms somehow... 
*/ 

/* __ICC__ is set by Makefile.am, there must be a cleaner way to detect icc compiler.
   Beware icc sets __GCC__*/
#if defined(CRLIBM_TYPECPU_ITANIUM) && defined(__GCC__)  && !defined(__ICC__) 
#define PROCESSOR_HAS_FMA 1
#define FMA(a,b,c)  /* r = a*b + c*/                 \
({                                                   \
  double _a, _b,_c,_r;                               \
  _a=a; _b=b;_c=c;                                   \
  __asm__ ("fma %0 = %1, %2, %3\n ;;\n"              \
		       : "=f"(_r)                    \
		       : "f"(_a), "f"(_b), "f"(_c)   \
		       );                            \
  _r;                                                \
})


#define FMS(a,b,c)  /* r = a*b - c*/                 \
({                                                   \
  double _a, _b, _c, _r;                             \
  _a=a; _b=b;_c=c;                                   \
  __asm__ ("fms %0 = %1, %2, %3\n ;;\n"              \
		       : "=f"(_r)                    \
		       : "f"(_a), "f"(_b), "f"(_c)   \
		       );                            \
  _r;                                                \
  })
#endif /* defined(CRLIBM_TYPECPU_ITANIUM) && defined(__GCC__) && !defined(__ICC) */


/* __ICC__ is set by Makefile.am, there must be a cleaner way to detect icc compiler.*/
#if defined(CRLIBM_TYPECPU_ITANIUM) && defined(__ICC__) 
#define PROCESSOR_HAS_FMA 1
/* Table 1-17: legal floating-point precision completers (.pc) */
typedef enum {
    _PC_S        = 1        /* single .s */
   ,_PC_D        = 2        /* double .d */
   ,_PC_NONE     = 3        /* dynamic   */
} _Asm_pc;

/* Table 1-22: legal getf/setf floating-point register access completers */
typedef enum {
    _FR_S        = 1        /* single form      .s   */
   ,_FR_D        = 2        /* double form      .d   */
   ,_FR_EXP      = 3        /* exponent form    .exp */
   ,_FR_SIG      = 4        /* significand form .sig */
} _Asm_fr_access;

/* Table 1-24: legal floating-point FPSR status field completers (.sf) */
typedef enum {
    _SF0         = 0        /* FPSR status field 0 .s0 */
   ,_SF1         = 1        /* FPSR status field 1 .s1 */
   ,_SF2         = 2        /* FPSR status field 2 .s2 */
   ,_SF3         = 3        /* FPSR status field 3 .s3 */
} _Asm_sf;

#define FMA(a,b,c)  /* r = a*b + c*/                 \
   _Asm_fma( _PC_D, a, b, c, _SF0 );              


#define FMS(a,b,c)  /* r = a*b - c*/                 \
   _Asm_fms( _PC_D, a, b, c, _SF0 );              

#endif /*defined(CRLIBM_TYPECPU_ITANIUM) && defined(__ICC)*/








#ifdef WORDS_BIGENDIAN
 #define DB_ONE    {{0x3ff00000, 0x00000000}}
#else
 #define DB_ONE    {{0x00000000 ,0x3ff00000}}
#endif






extern const scs scs_zer, scs_half, scs_one, scs_two, scs_sixinv;


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
#define Add12Cond(s, r, a, b)      \
        {double _z, _a=a, _b=b;    \
         s = _a + _b;              \
         if (ABS(a) > ABS(b)){     \
           _z = s - _a;            \
           r = _b - _z;            \
         }else {                   \
           _z = s - _b;            \
           r = _a - _z;}}                          

/*
 *  computes s and r such that s + r = a + b,  with s = a @+ b exactly 
 * under the condition  a >= b
 */
#define Add12(s, r, a, b)          \
        {double _z, _a=a, _b=b;    \
         s = _a + _b;              \
         _z = s - _a;              \
         r = _b - _z;   }            


/*
 * computes r1, r2, r3 such that r1 + r2 + r3 = a + b + c exactly 
 */
#define Fast3Sum(r1, r2, r3, a, b,  c) \
        {double u, v, w;               \
         Fast2Sum(u, v, b, c);         \
         Fast2Sum(r1, w, a, u);        \
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

#define Add22Cond(zh,zl,xh,xl,yh,yl)                                                   \
do {                                                                                   \
  double _r,_s;                                                                        \
  _r = (xh)+(yh);                                                                      \
  _s = ((ABS(xh)) > (ABS(yh)))? ((xh)-_r+(yh)+(yl)+(xl)) : ((yh)-_r+(xh)+(xl)+(yl));   \
  *zh = _r+_s;                                                                         \
  *zl = (_r - (*zh)) + _s;                                                             \
} while(2+2==5)

  

#define Add22(zh,zl,xh,xl,yh,yl)         \
do {                                     \
double _r,_s;                            \
_r = (xh)+(yh);                          \
_s = ((((xh)-_r) +(yh)) + (yl)) + (xl);  \
*zh = _r+_s;                             \
*zl = (_r - (*zh)) + _s;                 \
} while(0)

#endif /* ADD22_AS_FUNCTIONS */



#ifdef PROCESSOR_HAS_FMA
/* One of the nice things with the fused multiply-and-add is that it
   greatly simplifies the double-double multiplications : */
#define Mul12(rh,rl,u,v)                             \
{                                                    \
  *rh = u*v;                                         \
  *rl = FMS(u,v, *rh);                               \
}

#define Mul22(pzh,pzl, xh,xl, yh,yl)                  \
{                                                     \
double ph, pl;                                        \
  ph = xh*yh;                                         \
  pl = FMS(xh, yh,  ph);                              \
  pl = FMA(xh,yl, pl);                                \
  pl = FMA(xl,yh,pl);                                 \
  *pzh = ph+pl;					      \
  *pzl = ph - (*pzh);                                 \
  *pzl += pl;                                         \
}


/* besides we don't care anymore about overflows in the mult  */
#define Mul12Cond Mul12    
#define Mul22cond Mul22


#else /* ! PROCESSOR_HAS_FMA */


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
  double u, v;                                               \
  db_number _a=a, _b=b;                                      \
                                                             \
  if (_a.i[HI]>0x7C900000) u = _a*two_em53;                  \
  else            u = _a;                                    \
  if (_b.i[HI]>0x7C900000) v = _b*two_em53;                  \
  else            v = _b;                                    \
                                                             \
  Mul12(rh, rl, u, v);                                       \
                                                             \
  if (_a.i[HI]>0x7C900000) {*rh *= two_e53; *rl *= two_e53;} \
  if (_b.i[HI]>0x7C900000) {*rh *= two_e53; *rl *= two_e53;} \
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
  up = (xh)*c;        vp = (yh)*c;		      \
  u1 = ((xh)-up)+up;  v1 = ((yh)-vp)+vp;	      \
  u2 = (xh)-u1;       v2 = (yh)-v1;                   \
  						      \
  mh = (xh)*(yh);				      \
  ml = (((u1*v1-mh)+(u1*v2))+(u2*v1))+(u2*v2);	      \
						      \
  ml += (xh)*(yl) + (xl)*(yh);			      \
  *zh = mh+ml;					      \
  *zl = mh - (*zh) + ml;                              \
}



#endif /* DEKKER_AS_FUNCTIONS */

#endif /* PROCESSOR_HAS_FMA */



/* In the following the one-line computation of _cl was split so that
   icc(8.1) would compile it properly. It's a bug of icc */

#if DEKKER_AS_FUNCTIONS
extern void Div22(double *z, double *zz, double x, double xx, double y, double yy);
#else
#define  Div22(pzh,pzl,xh,xl,yh,yl)  {           \
  double _ch,_cl,_uh,_ul;                        \
  _ch=(xh)/(yh);   Mul12(&_uh,&_ul,_ch,(yh));    \
  _cl=((xh)-_uh);                                \
  _cl -= _ul;                                    \
  _cl += (xl);                                   \
  _cl -= _ch*(yl);                               \
  _cl /= (yh);                                   \
  *pzh=_ch+_cl;   *pzl=(_ch-(*pzh))+_cl;         \
}
#endif /* DEKKER_AS_FUNCTIONS */

/* A few prototypes that are not worth being in crlibm.h */

/*
 * Make a trigonometric range reduction with a relative
 * error less than 2^(-150) 
 */

extern int rem_pio2_scs(scs_ptr,  scs_ptr);



#endif /*CRLIBM_PRIVATE_H*/
