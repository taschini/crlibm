
#ifndef __DOUBLE_EXTENDED_H
#define __DOUBLE_EXTENDED_H

/* For debugging */
typedef union {
  int i[3];                 
  long double d;
} db_ext_number;


/**************************************************************************************/
/*********************************Rounding tests***************************************/
/**************************************************************************************/

#if (defined(CRLIBM_TYPECPU_X86) || defined(CRLIBM_TYPECPU_AMD64))

static const unsigned short RN_Double=(_FPU_DEFAULT & ~_FPU_EXTENDED)|_FPU_DOUBLE;
static const unsigned short RN_DoubleExt =_FPU_DEFAULT;

#define DOUBLE_EXTENDED_MODE  _FPU_SETCW(RN_DoubleExt)
#define BACK_TO_DOUBLE_MODE   _FPU_SETCW(RN_Double)

#define DE_EXP 2
#define DE_MANTISSA_HI 1
#define DE_MANTISSA_LO 0


/*  
Two rounding tests to the nearest.  On Pentium 3, gcc3.3, the second
is faster by 12 cycles (and also improves the worst-case time by 60
cycles since it doesn't switch processor rounding mode in this
case). However it uses a coarser error estimation.
*/

#define TEST_AND_RETURN_RN_ZIV(y,rncst)  \
{ double yh, yl;                         \
  yh = (double) y;                       \
  yl = y-yh;                             \
  BACK_TO_DOUBLE_MODE;                   \
  if(yh==yh + yl*rncst)   return yh;     \
  DOUBLE_EXTENDED_MODE;                  \
}


/* This test works by observing the bits of your double-extended after the 53rd.

   mask should be  7ff   if you trust your 64 bits (hum)
                   7fe   if you trust 63 (if you have proven that maxepsilon<2^(-63) )
                   7fc                62
                   7f8                61
                   7f0                60   etc
 */

#define TEST_AND_RETURN_RN(_y, _mask)                       \
{                                                           \
  db_ext_number _z;   double _yh;                           \
  int _lo;                                                  \
  _z.d = _y;                                                \
  _yh = (double) _y;                                        \
  _lo = _z.i[DE_MANTISSA_LO] &(_mask);                      \
  if((_lo!=(0x3ff&(_mask))) && (_lo!= (0x400&(_mask)))) {   \
    BACK_TO_DOUBLE_MODE;                                    \
    return _yh;                                             \
  }                                                         \
}


/* Use this one if you want a final computation step to overlap with
   the rounding test. Examples: multiplication by a sign or by a power of 2 */

#define TEST_AND_RETURN_RN2(_ytest, _yreturn, _mask)         \
{                                                            \
  db_ext_number _z;   double _y_return_d;                    \
  int _lo;                                                   \
  _z.d = _ytest;                                             \
  _y_return_d = (double) (_yreturn);                         \
  _lo = _z.i[DE_MANTISSA_LO] &(_mask);                       \
  if((_lo!=(0x3ff&(_mask))) && (_lo!= (0x400&(_mask)))) {    \
    BACK_TO_DOUBLE_MODE;                                     \
    return _y_return_d;                                      \
  }                                                          \
}




#else /* defined(CRLIBM_TYPECPU_X86) || defined(CRLIBM_TYPECPU_AMD64) */





#if !defined(CRLIBM_TYPECPU_ITANIUM)
#error "This file should be compiled only for IA32 or IA64 architecture "
#endif
#if !defined(__ICC__)
#error "Use icc, version 8.1 or higher to compile for IA64 architecture"
#endif


#define DOUBLE_EXTENDED_MODE {}
#define BACK_TO_DOUBLE_MODE {}


#define TEST_AND_RETURN_RN(_y, _mask)                                                \
{   unsigned long int _i1, _m;                                                              \
    _i1 = _Asm_getf(4/*_FR_SIG*/, _y);                                                    \
    _m =  _i1 & (0x7ff&(_mask));                                                      \
    if(__builtin_expect((_m!=(0x3ff&(_mask))) && (_m != (0x400&(_mask))), 1+1==2))    \
      return (double)_y;                                                             \
}


/* Use this one if you want a final computation step to overlap with
   the rounding test. Examples: multiplication by a sign or by a power of 2 */

#define TEST_AND_RETURN_RN2(_ytest, _yreturn, _mask)         \
{                                                            \
printf("Rounding test not yet implemented\n\n");             \
}


#endif /* defined(CRLIBM_TYPECPU_X86) || defined(CRLIBM_TYPECPU_AMD64) */




/**************************************************************************************/
/************************Double double-extended arithmetic*****************************/
/**************************************************************************************/



#define Add12_ext(prh, prl, a, b)       \
{                                       \
  long double _z, _a, _b;               \
  _a = a;   _b = b;                     \
  *prh = _a + _b;                       \
  _z = *prh - _a;                       \
  *prl = _b - _z;                       \
}




#if (defined(CRLIBM_TYPECPU_X86) || defined(CRLIBM_TYPECPU_AMD64))
#define Mul12_ext(prh,prl,u,v)                         \
{                                                      \
  const long double c  = 4294967297.L; /* 2^32 +1 */   \
  long double up, u1, u2, vp, v1, v2;                  \
  long double _u =u, _v=v;                             \
                                                       \
  up = _u*c;        vp = _v*c;                         \
  u1 = (_u-up)+up;  v1 = (_v-vp)+vp;                   \
  u2 = _u-u1;       v2 = _v-v1;                        \
                                                       \
  *prh = _u*_v;                                        \
  *prl = u1*v1 - *prh;                                 \
  *prl = *prl + u1*v2;                                 \
  *prl = *prl + u2*v1;                                 \
  *prl = *prl + u2*v2;                                 \
}

#define Mul22_ext(prh,prl, ah,al, bh,bl)               \
{                                                      \
  long double mh, ml;                                  \
  Mul12_ext(&mh,&ml,(ah),(bh));		               \
  ml += (ah)*(bl) + (al)*(bh);			       \
  Add12_ext(prh,prl, mh,ml);                           \
}

#define FMA22_ext(prh,prl, ah,al, bh,bl, ch,cl)        \
{                                                      \
  Mul22_ext(prh,prl, (ah),(al), (bh),(bl));            \
  Add22_ext(prh,prl, ch,cl, *prh, *prl);               \
}



#else  /* defined(CRLIBM_TYPECPU_X86) || defined(CRLIBM_TYPECPU_AMD64) */



#if 0 /* shouldn't be here */ 
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
#endif



#define ULL(bits) 0x##bits##uLL


/* The following is macros by Alexey Ershov in a more readable form */

#define __X_SQR_X( __resx__, __x_a__ ) \
    {                                                                                               \
        L_FLOAT_TYPE __xsqrx_r_hi__,__xsqrx_r_lo__,__xsqrx_t__,__xsqrx_two__;                       \
        __xsqrx_r_hi__ = (__x_a__).hi * (__x_a__).hi;                                               \
        __xsqrx_t__ = (__x_a__).hi * (__x_a__).lo;                                                  \
        __xsqrx_r_lo__ = _Asm_fms( _PC_NONE, (__x_a__).hi, (__x_a__).hi, __xsqrx_r_hi__, _SF1 );    \
        __xsqrx_two__ = _Asm_setf( _FR_EXP, 0xffff + 1 );                                           \
        __xsqrx_r_lo__ = __xsqrx_two__ * __xsqrx_t__ + __xsqrx_r_lo__;                              \
        (__resx__).hi = __xsqrx_r_hi__; (__resx__).lo = __xsqrx_r_lo__;                             \
    }

/* res = a*a
 *   res in X format
 *   a   in L format
 */
#define __X_SQR_L( __resx__, __l_a__ ) \
    {                                                                                               \
        L_FLOAT_TYPE __xsqrl_r_hi__,__xsqrl_r_lo__;                                                 \
        __xsqrl_r_hi__ = (__l_a__) * (__l_a__);                                                     \
        __xsqrl_r_lo__ = _Asm_fms( _PC_NONE, (__l_a__), (__l_a__), __xsqrl_r_hi__, _SF1 );          \
        (__resx__).hi = __xsqrl_r_hi__; (__resx__).lo = __xsqrl_r_lo__;                             \
    }

/* res = a*b
 *   res, a and b in X format
 */
#define __X_MUL_XX( __resx__, __x_a__, __x_b__ ) \
    {                                                                                               \
        L_FLOAT_TYPE __xmulxx_r_hi__,__xmulxx_r_lo__;                                               \
        L_FLOAT_TYPE __xmulxx_t1__,__xmulxx_t2__,__xmulxx_t3__;                                     \
        __xmulxx_r_hi__ = (__x_a__).hi * (__x_b__).hi;                                              \
        __xmulxx_t1__ = (__x_a__).hi * (__x_b__).lo;                                                \
        __xmulxx_t2__ = _Asm_fms( _PC_NONE, (__x_a__).hi, (__x_b__).hi, __xmulxx_r_hi__, _SF1 );    \
        __xmulxx_t3__ = (__x_a__).lo * (__x_b__).hi + __xmulxx_t1__;                                \
        __xmulxx_r_lo__ = (__xmulxx_t2__ + __xmulxx_t3__);                                          \
        (__resx__).hi = __xmulxx_r_hi__; (__resx__).lo = __xmulxx_r_lo__;                           \
    }



#define Mul12_ext( prh,prl, a, b )                              \
    {                                                           \
      *prh = (a) * (b);                                         \
      *prl = _Asm_fms( 3/*_PC_NONE*/, (a), (b), *prh, 1 );      \
    }


#if 0 /* transcription of Alexey's */
#define Mul22_ext( prh,prl, ah,al, bh,bl ) \
    {                                                            \
        long double _t1,_t2,_t3;                                 \
        *prh = (ah) * (bh);                                      \
        _t1 = (ah)*(bl);                                         \
        _t2 = _Asm_fms( 3/*_PC_NONE*/, (ah), (bh), *prh, 1 );    \
        _t3 = (al) * (bh) + _t1;                                 \
        *prl = (_t2 + _t3);                                      \
    }
#else
#define Mul22_ext( prh,prl, ah,al, bh,bl ) \
{                                                     \
  __fpreg ph, pl;                                   \
  ph = (ah)*(bh);                                         \
  pl = _Asm_fms( 3/*_PC_NONE*/, ah, bh, ph, 1/*_SF1*/ );;  \
  pl = (ah)*(bl) + pl;                                    \
  pl = (al)*(bh) + pl;                                    \
  Add12_ext(prh,prl, ph,pl); \
}
#endif

/* res = a*b + c, assume |a*b| <= |c|
 *   res, a, b and c in X format
 */
#define __X_FMA_GREATER_XXX( __resx__, __x_a__, __x_b__, __x_c__ ) \
    {                                                                                               \
        L_FLOAT_TYPE __xfmagxxx_r_hi__,__xfmagxxx_r_lo__;                                           \
        L_FLOAT_TYPE __xfmagxxx_t1__,__xfmagxxx_t2__;                                               \
        L_FLOAT_TYPE __xfmagxxx_t3__,__xfmagxxx_t4__;                                               \
        __xfmagxxx_r_hi__ = (__x_a__).hi * (__x_b__).hi + (__x_c__).hi;                             \
        __xfmagxxx_t1__ = (__x_a__).lo * (__x_b__).hi + (__x_c__).lo;                               \
        __xfmagxxx_t2__ = (__xfmagxxx_r_hi__ - (__x_c__).hi);                                       \
        __xfmagxxx_t3__ = (__x_a__).hi * (__x_b__).lo + __xfmagxxx_t1__;                            \
        __xfmagxxx_t4__ = _Asm_fms( _PC_NONE, (__x_a__).hi, (__x_b__).hi, __xfmagxxx_t2__, _SF1 );  \
        __xfmagxxx_r_lo__ = (__xfmagxxx_t3__ + __xfmagxxx_t4__);                                    \
        (__resx__).hi = __xfmagxxx_r_hi__; (__resx__).lo = __xfmagxxx_r_lo__;                       \
    }

/* res = a*b + c, assume |a*b| <= |c|
 *   res, b and c in X format
 *   a in L format
 */
#define __X_FMA_GREATER_LXX( __resx__, __l_a__, __x_b__, __x_c__ ) \
    {                                                                                               \
        L_FLOAT_TYPE __xfmaglxx_r_hi__,__xfmaglxx_r_lo__;                                           \
        L_FLOAT_TYPE __xfmaglxx_t2__;                                                               \
        L_FLOAT_TYPE __xfmaglxx_t3__,__xfmaglxx_t4__;                                               \
        __xfmaglxx_r_hi__ = (__l_a__) * (__x_b__).hi + (__x_c__).hi;                                \
        __xfmaglxx_t2__ = (__xfmaglxx_r_hi__ - (__x_c__).hi);                                       \
        __xfmaglxx_t3__ = (__l_a__) * (__x_b__).lo + (__x_c__).lo;                                  \
        __xfmaglxx_t4__ = _Asm_fms( _PC_NONE, (__l_a__), (__x_b__).hi, __xfmaglxx_t2__, _SF1 );     \
        __xfmaglxx_r_lo__ = (__xfmaglxx_t3__ + __xfmaglxx_t4__);                                    \
        (__resx__).hi = __xfmaglxx_r_hi__; (__resx__).lo = __xfmaglxx_r_lo__;                       \
    }

/* res = a*b + c, assume |a*b| <= |c|
 *   res, a and b in X format
 *   c in L format
 */
#define __X_FMA_GREATER_XXL( __resx__, __x_a__, __x_b__, __l_c__ ) \
    {                                                                                               \
        L_FLOAT_TYPE __xfmagxxl_r_hi__,__xfmagxxl_r_lo__;                                           \
        L_FLOAT_TYPE __xfmagxxl_t1__,__xfmagxxl_t2__;                                               \
        L_FLOAT_TYPE __xfmagxxl_t3__,__xfmagxxl_t4__;                                               \
        __xfmagxxl_r_hi__ = (__x_a__).hi * (__x_b__).hi + (__l_c__);                                \
        __xfmagxxl_t1__ = (__x_a__).lo * (__x_b__).hi;                                              \
        __xfmagxxl_t2__ = (__xfmagxxl_r_hi__ - (__l_c__));                                          \
        __xfmagxxl_t3__ = (__x_a__).hi * (__x_b__).lo + __xfmagxxl_t1__;                            \
        __xfmagxxl_t4__ = _Asm_fms( _PC_NONE, (__x_a__).hi, (__x_b__).hi, __xfmagxxl_t2__, _SF1 );  \
        __xfmagxxl_r_lo__ = (__xfmagxxl_t3__ + __xfmagxxl_t4__);                                    \
        (__resx__).hi = __xfmagxxl_r_hi__; (__resx__).lo = __xfmagxxl_r_lo__;                       \
    }

#if 0
#define FMA22_ext(prh,prl, ah,al, bh,bl, ch,cl)        \
{                                                      \
  Mul22_ext(prh,prl, (ah),(al), (bh),(bl));            \
  Add22_ext(prh,prl, ch,cl, *prh, *prl);               \
}
#else
#define FMA22_ext( prh,prl, ah,al,  bh,bl, ch,cl) \
    {                                                                                               \
        __fpreg __xfmagxxx_r_hi__,__xfmagxxx_r_lo__,                                            \
                __xfmagxxx_t1__,__xfmagxxx_t2__,                                               \
                __xfmagxxx_t3__,__xfmagxxx_t4__;                                               \
        __xfmagxxx_r_hi__ = ah * bh + ch;                             \
        __xfmagxxx_t1__ = al * bh + cl;                               \
        __xfmagxxx_t2__ = __xfmagxxx_r_hi__ - ch;                                       \
        __xfmagxxx_t3__ = ah * bl + __xfmagxxx_t1__;                            \
        __xfmagxxx_t4__ = _Asm_fms( 3/*_PC_NONE*/, ah, bh, __xfmagxxx_t2__, 1/*_SF1*/ );  \
        __xfmagxxx_r_lo__ = (__xfmagxxx_t3__ + __xfmagxxx_t4__);                                    \
        *prh = __xfmagxxx_r_hi__; *prl = __xfmagxxx_r_lo__;                       \
    }
#endif

#endif    /* defined(CRLIBM_TYPECPU_X86) || defined(CRLIBM_TYPECPU_AMD64) */





/* Computes a*b+c under the condition that a*b << c  
   and under the condition there will be no overflow
   which is easy to ensure as the inputs of crlibm functions are doubles */




#define  Div22_ext(prh,prl,xh,xl,yh,yl)             \
{                                                   \
  long double ch,cl,uh,ul;                          \
  ch = (xh)/(yh);                                   \
  Mul12_ext(&uh,&ul,ch,(yh));                       \
  cl = (xh)-uh;                                     \
  cl = cl - ul;                                     \
  cl = cl + (xl);                                   \
  cl = cl - ch*(yl);                                \
  cl = cl / (yh);                                   \
  Add12(prh,prl, ch, cl) ;                          \
}


#define Add22_ext(prh,prl,xh,xl,yh,yl)   \
do {                                     \
  long double _r,_s;                     \
  _r = (xh)+(yh);                        \
  _s = (xh)-_r;                          \
  _s = _s + (yh);                        \
  _s = _s + (yl);                        \
  _s = _s + (xl);                        \
  Add12_ext(prh,prl,_r,_s);              \
} while(0)

#endif /* ifndef __DOUBLE_EXTENDED_H*/
