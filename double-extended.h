
#ifndef __DOUBLE_EXT_H
#define __DOUBLE_EXT_H

#ifdef CRLIBM_TYPECPU_X86

static const unsigned short RN_Double=(_FPU_DEFAULT & ~_FPU_EXTENDED)|_FPU_DOUBLE;
static const unsigned short RN_DoubleExt =_FPU_DEFAULT;

#define DOUBLE_EXTENDED_MODE  _FPU_SETCW(RN_DoubleExt)
#define BACK_TO_DOUBLE_MODE   _FPU_SETCW(RN_Double)


typedef union {
  int i[3];                 
  long double d;
} db_ext_number;

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
#endif /* CRLIBM_TYPECPU_X86*/



#ifdef ITANIUMICC
typedef          __int64  INT64;
typedef   signed __int64 SINT64;
typedef unsigned __int64 UINT64;

/* FP register type */
typedef __fpreg L_FLOAT_TYPE;

/* Almost the same as the previous, except exponent field smaller, and morally in memory */
typedef long double LC_FLOAT_TYPE;

/* The double-double-ext type, using registers */
typedef struct __X_FLOAT_TYPE_TAG {
    L_FLOAT_TYPE hi,lo; /* order is critical! */
} X_FLOAT_TYPE;

/* The double-double-ext type, in memory */
typedef struct __XC_FLOAT_TYPE_TAG {
    LC_FLOAT_TYPE hi,lo; /* order is critical! */
} XC_FLOAT_TYPE;

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




#define ULL(bits) 0x##bits##uLL

#if (!defined(EM64T) && defined(__linux__) && defined(IA32))
# define LDOUBLE_ALIGN 12   /* IA32 Linux: 12-byte alignment */
#else
# define LDOUBLE_ALIGN 16   /* EM64T, IA32 Win or IPF Win/Linux: 16-byte alignm\
			       ent */
#endif

#if (LDOUBLE_ALIGN == 16)
#define _XPD_ ,0x0000,0x0000,0x0000
#else /*12*/
#define _XPD_ ,0x0000
#endif

#define LDOUBLE_HEX(w4,w3,w2,w1,w0) 0x##w0,0x##w1,0x##w2,0x##w3,0x##w4 _XPD_ /*LITTLE_ENDIAN*/



/* Load XC constant data and convert to X format */
#define __X_CONVERT_XC2X( __x__, __xc__ ) \
    (((__x__).hi = (__xc__).hi), ((__x__).lo = (__xc__).lo))

/* res = a*a
 *   res and a in X format
 */
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


#define Mul12_ext(_prh,_prl,_u,_v)                    \
{                                               \
  *_prh = _u*_v;                                  \
  *_prl = *_prh - _u*_v;                             \
}


#define Mul22_ext(pzh,pzl, xh,xl, yh,yl)              \
{                                                     \
long double ph, pl;                                   \
  ph = xh*yh;                                         \
  pl = xh*yh - ph;                                    \
  pl = xh*yl + pl;                                    \
  pl = xl*yh + pl;                                    \
  *pzh = ph+pl;					      \
  *pzl = ph - (*pzh);                                 \
  *pzl += pl;                                         \
}




#else /*ITANIUMICC*/
#define Mul12_ext(rh,rl,u,v)                        \
{                                               \
  const long double c  = 4294967297.L; /* 2^32 +1 */   \
  long double up, u1, u2, vp, v1, v2;                \
  long double _u =u, _v=v;                           \
                                                \
  up = _u*c;        vp = _v*c;                  \
  u1 = (_u-up)+up;  v1 = (_v-vp)+vp;            \
  u2 = _u-u1;       v2 = _v-v1;                 \
                                                \
  *rh = _u*_v;                                  \
  *rl = (((u1*v1-*rh)+(u1*v2))+(u2*v1))+(u2*v2);\
}

#define Mul22_ext(zh,zl,xh,xl,yh,yl)                  \
{                                                     \
long double mh, ml;                                   \
						      \
  const long double c = 4294967297.L /*2^32+1*/;      \
  long double up, u1, u2, vp, v1, v2;		      \
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


#endif /*ITANIUMICC*/

#define  Div22_ext(zh,zl,xh,xl,yh,yl)\
{long double ch,cl,uh,ul;  \
           ch=(xh)/(yh);   Mul12_ext(&uh,&ul,ch,(yh));  \
           cl=(((((xh)-uh)-ul)+(xl))-ch*(yl))/(yh);   zh=ch+cl;   zl=(ch-zh)+cl;\
}



#define Add12_ext(s, r, a, b)         \
        { long double _z, _a=a, _b=b;    \
         s = _a + _b;             \
         _z = s - _a;              \
         r = _b - _z; }            


#define Add22_ext(zh,zl,xh,xl,yh,yl) \
do {\
long double r,s;\
r = (xh)+(yh);\
s = (xh)-r+(yh)+(yl)+(xl);\
*zh = r+s;\
*zl = r - (*zh) + s;\
} while(0)


#endif // ifndef __DOUBLE_EXT_H
