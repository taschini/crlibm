/* file: exp_itanium.c


   Author Florent de Dinechin, ENS-Lyon
     Heavily inspired by code from Alexey V.Ershov, intel
     and Christoph Lauter, Technische Universitaet Muenchen

WARNING : This code is dirty and experimental, and remains here for
history. A cleaner, portable version of an exponential using
double-extended arithmetic will be available as exp-de.c

   To test within crlibm: (tested with Intel icc compiler version 8.1)

icc -mcpu=itanium  -Qoption,cpp,--extended_float_types \
    -IPF_fp_speculationsafe -c exp-itanium.c;\
    mv exp-itanium.o exp_fast.o; make


icc -mcpu=itanium2  -Qoption,cpp,--extended_float_types \
    -IPF_fp_speculationsafe -c exp-itanium.c;
    mv exp-itanium.o exp_fast.o; make

This file is completely self-contained so that we can change the crlibm infrastructure without bothering maintaining it.

*/


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



/* T1[i1] = 2^(i1/(K^1))  i1 = [-L/2..L/2]
 * T2[i2] = 2^(i2/(K^2))  i2 = [-L/2..L/2]
 */
#define TX ((const XC_FLOAT_TYPE *)_TXL)
__declspec(align(16)) static const unsigned short _TXL[] = {
/* -64*/  LDOUBLE_HEX(3ffe,b504,f333,f9de,6484),LDOUBLE_HEX(3fbd,b2fb,1366,ea95,7d3e)   /* ~0.70710678118654752438189403651591646848828531801700 T1[  0] hi,lo */
/* -64*/ ,LDOUBLE_HEX(3ffe,ff4e,cb59,511e,c8a5),LDOUBLE_HEX(3fbc,c06e,885f,bc63,75f1)   /* ~0.99729605608547012624747266085201147234329255297780 T2[  0] hi,lo */
/* -63*/ ,LDOUBLE_HEX(3ffe,b600,93a8,5ed5,f76c),LDOUBLE_HEX(bfbd,9566,7ff0,b0cc,0214)   /* ~0.71094630108458277992049267890450892082299105823040 T1[  1] hi,lo */
/* -63*/ ,LDOUBLE_HEX(3ffe,ff51,8f3a,34cb,1c70),LDOUBLE_HEX(3fbd,d414,fb77,db11,52bc)   /* ~0.99733824893045803280239303134635520109441131353378 T2[  1] hi,lo */
/* -62*/ ,LDOUBLE_HEX(3ffe,b6fd,91e3,28d1,7791),LDOUBLE_HEX(3fb9,e2cb,e1bb,aa83,4b3f)   /* ~0.71480666919598500561603207636984791406575823202729 T1[  2] hi,lo */
/* -62*/ ,LDOUBLE_HEX(3ffe,ff54,5322,c329,f86f),LDOUBLE_HEX(3fba,e9b1,4513,ba49,7f30)   /* ~0.99738044356050881746869260502386111966188764199614 T2[  2] hi,lo */
/* -61*/ ,LDOUBLE_HEX(3ffe,b7fb,efca,8ca4,1e7c),LDOUBLE_HEX(3fbc,fc36,9e7c,4277,ff37)   /* ~0.71868799872449116278826097836684994035749696195125 T1[  3] hi,lo */
/* -61*/ ,LDOUBLE_HEX(3ffe,ff57,1712,fc50,9e7f),LDOUBLE_HEX(3fbd,8700,fef5,d3b7,14da)   /* ~0.99742263997569800127026770630855878607690101489424 T2[  3] hi,lo */
/* -60*/ ,LDOUBLE_HEX(3ffe,b8fb,af47,62fb,9ee9),LDOUBLE_HEX(3fbb,dc3c,bbc2,b35b,2d0d)   /* ~0.72259040348852331001267348953298608194018015637993 T1[  4] hi,lo */
/* -60*/ ,LDOUBLE_HEX(3ffe,ff59,db0a,e054,50ba),LDOUBLE_HEX(3fbb,f679,bcc2,0658,8201)   /* ~0.99746483817610110853783128570526628209336195141077 T2[  4] hi,lo */
/* -59*/ ,LDOUBLE_HEX(3ffe,b9fc,d245,2c0b,9deb),LDOUBLE_HEX(bfbb,d96c,65d5,3b9f,5478)   /* ~0.72651399792452628242879139608412231154943583533167 T1[  5] hi,lo */
/* -59*/ ,LDOUBLE_HEX(3ffe,ff5c,9f0a,6f4a,5174),LDOUBLE_HEX(bfbd,b0bf,02f2,53a4,00e8)   /* ~0.99750703816179366674628259392676454808679409325122 T2[  5] hi,lo */
/* -58*/ ,LDOUBLE_HEX(3ffe,baff,5ab2,133e,45fb),LDOUBLE_HEX(3fbd,e9aa,33a4,8b27,0718)   /* ~0.73045889709032349430091121678110255288629559800028 T1[  6] hi,lo */
/* -58*/ ,LDOUBLE_HEX(3ffe,ff5f,6311,a947,e33b),LDOUBLE_HEX(bfbc,a528,cf59,14b2,d1a3)   /* ~0.99754923993285120651470718189379738305433420464396 T2[  6] hi,lo */
/* -57*/ ,LDOUBLE_HEX(3ffe,bc03,4a7e,f2e9,fb0d),LDOUBLE_HEX(bfbc,a3fa,fef4,e9a8,9a84)   /* ~0.73442521666849096343950356602370277414593147113919 T1[  7] hi,lo */
/* -57*/ ,LDOUBLE_HEX(3ffe,ff62,2720,8e62,48da),LDOUBLE_HEX(bfbd,d1c7,8694,72f5,3da0)   /* ~0.99759144348934926176900722660789710971585009247064 T2[  7] hi,lo */
/* -56*/ ,LDOUBLE_HEX(3ffe,bd08,a39f,580c,36bf),LDOUBLE_HEX(bfbd,aefd,c093,25e0,a10c)   /* ~0.73841307296974965571198165625865783567860489711165 T1[  8] hi,lo */
/* -56*/ ,LDOUBLE_HEX(3ffe,ff64,eb37,1eae,c555),LDOUBLE_HEX(bfbc,e4bf,fcf3,e4d6,86d6)   /* ~0.99763364883136336952506109665428368771245004609227 T2[  8] hi,lo */
/* -55*/ ,LDOUBLE_HEX(3ffe,be0f,6809,8609,93e2),LDOUBLE_HEX(3fbd,9334,4593,7562,b2dd)   /* ~0.74242258293637625025680054635657256767444778233766 T1[  9] hi,lo */
/* -55*/ ,LDOUBLE_HEX(3ffe,ff67,af55,5a42,9bec),LDOUBLE_HEX(3fbc,dd7f,3c8d,5596,1f78)   /* ~0.99767585595896907005135367807469037870760075747966 T2[  9] hi,lo */
/* -54*/ ,LDOUBLE_HEX(3ffe,bf17,99b6,7a73,1083),LDOUBLE_HEX(bfbb,bf51,7aa1,a07a,3d7b)   /* ~0.74645386414563242460538591838137278955400688573718 T1[ 10] hi,lo */
/* -54*/ ,LDOUBLE_HEX(3ffe,ff6a,737b,4133,101b),LDOUBLE_HEX(bfbd,a2df,4f03,c8fb,22c7)   /* ~0.99771806487224190686897637436736374638712732121348 T2[ 10] hi,lo */
/* -53*/ ,LDOUBLE_HEX(3ffe,c021,3aa1,f0d0,8db0),LDOUBLE_HEX(3fbd,de67,649a,354e,a707)   /* ~0.75050703481321276010901355491000686015468090772628 T1[ 11] hi,lo */
/* -53*/ ,LDOUBLE_HEX(3ffe,ff6d,37a8,d395,6596),LDOUBLE_HEX(bfbb,c802,1a75,cc6d,fa67)   /* ~0.99776027557125742653478667198996276965772267431020 T2[ 11] hi,lo */
/* -52*/ ,LDOUBLE_HEX(3ffe,c12c,4cca,6670,9456),LDOUBLE_HEX(3fbd,f88a,fab3,4a01,0f6b)   /* ~0.75458221379671136985669427366119066391547676175832 T1[ 12] hi,lo */
/* -52*/ ,LDOUBLE_HEX(3ffe,ff6f,fbde,117e,e04f),LDOUBLE_HEX(bfbd,de6d,fc11,1a38,5732)   /* ~0.99780248805609117896666879210521017284918343648314 T2[ 12] hi,lo */
/* -51*/ ,LDOUBLE_HEX(3ffe,c238,d231,1e3d,6673),LDOUBLE_HEX(bfbd,d094,6839,bf6a,c8f1)   /* ~0.75867952059910734942697538296663140044984174892306 T1[ 13] hi,lo */
/* -51*/ ,LDOUBLE_HEX(3ffe,ff72,c01a,fb04,c470),LDOUBLE_HEX(3fbd,a44f,1ba6,a6fe,09f6)   /* ~0.99784470232681871711827303883524109551217406988143 T2[ 13] hi,lo */
/* -50*/ ,LDOUBLE_HEX(3ffe,c346,ccda,2497,6407),LDOUBLE_HEX(3fbc,83b2,1584,a2e0,e90a)   /* ~0.76279907537226915341865493558337618651421507820487 T1[ 14] hi,lo */
/* -50*/ ,LDOUBLE_HEX(3ffe,ff75,845f,903c,5662),LDOUBLE_HEX(bfbc,df95,80fe,a0cf,dc67)   /* ~0.99788691838351559735848655963152964432083535939455 T2[ 14] hi,lo */
/* -49*/ ,LDOUBLE_HEX(3ffe,c456,3ecc,5334,cb33),LDOUBLE_HEX(bfbd,cf43,20d2,b162,9ed8)   /* ~0.76694099892047800092224549284303236618143273517489 T1[ 15] hi,lo */
/* -49*/ ,LDOUBLE_HEX(3ffe,ff78,48ab,d13a,dac5),LDOUBLE_HEX(bfbd,f1ee,64ae,cf61,1403)   /* ~0.99792913622625737903775247628068711946980329230427 T2[ 15] hi,lo */
/* -48*/ ,LDOUBLE_HEX(3ffe,c567,2a11,5506,dadd),LDOUBLE_HEX(3fbc,f8ab,4325,9376,7cde)   /* ~0.77110541270397041179298150415100110421917634084820 T1[ 16] hi,lo */
/* -48*/ ,LDOUBLE_HEX(3ffe,ff7b,0cff,be15,9675),LDOUBLE_HEX(3fbb,db1c,1100,661d,aecd)   /* ~0.99797135585511962475912042802583812317607225850224 T2[ 16] hi,lo */
/* -47*/ ,LDOUBLE_HEX(3ffe,c679,90b5,aa24,5f79),LDOUBLE_HEX(3fbd,aa1c,d161,c55d,84aa)   /* ~0.77529243884249997954350233642095702180085936561226 T1[ 17] hi,lo */
/* -47*/ ,LDOUBLE_HEX(3ffe,ff7d,d15b,56e1,ce8b),LDOUBLE_HEX(3fba,a4a7,c651,158b,a410)   /* ~0.99801357727017790043245668019089578137936769053339 T2[ 17] hi,lo */
/* -46*/ ,LDOUBLE_HEX(3ffe,c78d,74c8,abb9,b15d),LDOUBLE_HEX(bfbc,fb17,471a,24ff,6207)   /* ~0.77950220011891848353015668782184377505473094061017 T1[ 18] hi,lo */
/* -46*/ ,LDOUBLE_HEX(3ffe,ff80,95be,9bb4,c859),LDOUBLE_HEX(bfbc,b742,b767,0074,c32e)   /* ~0.99805580047150777505760368968346085694065550342202 T2[ 18] hi,lo */
/* -45*/ ,LDOUBLE_HEX(3ffe,c8a2,d85c,8ffe,2c45),LDOUBLE_HEX(3fbc,c368,d3ed,6e1c,0385)   /* ~0.78373481998277644652211096354399444408045383170247 T1[ 19] hi,lo */
/* -45*/ ,LDOUBLE_HEX(3ffe,ff83,5a29,8ca3,c96c),LDOUBLE_HEX(bfbb,a438,6792,f0f4,537c)   /* ~0.99809802545918482083280032224337219304288737475872 T2[ 19] hi,lo */
/* -44*/ ,LDOUBLE_HEX(3ffe,c9b9,bd86,6e2f,27a3),LDOUBLE_HEX(bfbd,fe3c,0dab,f5dd,2d04)   /* ~0.78799042255394324325455321367250860475905938073992 T1[ 20] hi,lo */
/* -44*/ ,LDOUBLE_HEX(3ffe,ff86,1e9c,29c4,178d),LDOUBLE_HEX(bfbc,db63,47ba,fee3,03cb)   /* ~0.99814025223328461320889196106698193489137338474392 T2[ 20] hi,lo */
/* -43*/ ,LDOUBLE_HEX(3ffe,cad2,265e,4290,774e),LDOUBLE_HEX(bfbd,b7c9,6a5f,0390,8383)   /* ~0.79226913262624686507939325563754096037882845848798 T1[ 21] hi,lo */
/* -43*/ ,LDOUBLE_HEX(3ffe,ff88,e316,732a,f8bf),LDOUBLE_HEX(bfbc,9a27,1ddb,a851,3860)   /* ~0.99818248079388273078091028955860508631303673610091 T2[ 21] hi,lo */
/* -42*/ ,LDOUBLE_HEX(3ffe,cbec,14fe,f272,7c5d),LDOUBLE_HEX(bfba,b6f8,370b,a140,9231)   /* ~0.79657107567113344897104590391911926872126059606671 T1[ 22] hi,lo */
/* -42*/ ,LDOUBLE_HEX(3ffe,ff8b,a798,68ed,b341),LDOUBLE_HEX(bfbd,f2fb,91eb,2cdf,6fb7)   /* ~0.99822471114105475539649350857906995315715903416275 T2[ 22] hi,lo */
/* -41*/ ,LDOUBLE_HEX(3ffe,cd07,8b86,503d,cdd2),LDOUBLE_HEX(bfbd,ef64,73b9,8c84,14e6)   /* ~0.80089637784134667692226977608882521053601521998643 T1[ 23] hi,lo */
/* -41*/ ,LDOUBLE_HEX(3ffe,ff8e,6c22,0b21,8d8b),LDOUBLE_HEX(3fbc,dbab,2022,5186,afc9)   /* ~0.99826694327487627199325601057289247819426236674189 T2[ 23] hi,lo */
/* -40*/ ,LDOUBLE_HEX(3ffe,ce24,8c15,1f84,80e4),LDOUBLE_HEX(bfbb,ee53,e383,5069,c895)   /* ~0.80524516597462715409606881511450637844973243772983 T1[ 24] hi,lo */
/* -40*/ ,LDOUBLE_HEX(3ffe,ff91,30b3,59db,ce53),LDOUBLE_HEX(3fbd,9ced,2076,73bc,c0c0)   /* ~0.99830917719542286886983892268965234961797250434756 T2[ 24] hi,lo */
/* -39*/ ,LDOUBLE_HEX(3ffe,cf43,18cf,1919,18c1),LDOUBLE_HEX(3fbc,994f,1cc9,8dc0,21f2)   /* ~0.80961756759743187464176475409693978235736722126603 T1[ 25] hi,lo */
/* -39*/ ,LDOUBLE_HEX(3ffe,ff93,f54c,5531,bc88),LDOUBLE_HEX(3fbc,badb,2327,d830,4c8e)   /* ~0.99835141290277013741485956366261689254315569996833 T2[ 25] hi,lo */
/* -38*/ ,LDOUBLE_HEX(3ffe,d063,33da,ef2b,2595),LDOUBLE_HEX(bfbc,a4ae,8e6a,996c,abf8)   /* ~0.81401371092867388343282744611606460694019915536046 T1[ 26] hi,lo */
/* -38*/ ,LDOUBLE_HEX(3ffe,ff96,b9ec,fd38,9f54),LDOUBLE_HEX(3fbc,deab,3cf6,1b93,876b)   /* ~0.99839365039699367221533166105729151240666396915912 T2[ 26] hi,lo */
/* -37*/ ,LDOUBLE_HEX(3ffe,d184,df62,5169,9ac6),LDOUBLE_HEX(3fba,b8fb,b86d,56aa,3fd1)   /* ~0.81843372488348224388140378371403471646772231906652 T1[ 27] hi,lo */
/* -37*/ ,LDOUBLE_HEX(3ffe,ff99,7e95,5205,be1d),LDOUBLE_HEX(bfbc,de71,2326,2868,54d2)   /* ~0.99843588967816907111087545989569491666770773008465 T2[ 27] hi,lo */
/* -36*/ ,LDOUBLE_HEX(3ffe,d2a8,1d91,f12a,e45a),LDOUBLE_HEX(3fbb,9124,72be,1ef2,0143)   /* ~0.82287773907698242225553647566371751054248306900262 T1[ 28] hi,lo */
/* -36*/ ,LDOUBLE_HEX(3ffe,ff9c,4345,53ae,6082),LDOUBLE_HEX(3fbc,fe9f,77a5,1aa2,3b01)   /* ~0.99847813074637193503108739678353344970673788338899 T2[ 28] hi,lo */
/* -35*/ ,LDOUBLE_HEX(3ffe,d3cc,f099,859a,c379),LDOUBLE_HEX(3fbd,dfb2,b158,f1a9,8797)   /* ~0.82734588382809719876243362279844006934581557288765 T1[ 29] hi,lo */
/* -35*/ ,LDOUBLE_HEX(3ffe,ff9f,07fd,0247,ce60),LDOUBLE_HEX(bfb8,da60,27f4,240a,2f75)   /* ~0.99852037360167786826659064303157720132730901241302 T2[ 29] hi,lo */
/* -34*/ ,LDOUBLE_HEX(3ffe,d4f3,5aab,cfed,fa1f),LDOUBLE_HEX(3fbd,b243,bdff,4c4c,58b5)   /* ~0.83183829016336821750429372790414106475509470328688 T1[ 30] hi,lo */
/* -34*/ ,LDOUBLE_HEX(3ffe,ffa1,ccbc,5de7,4fcc),LDOUBLE_HEX(3fbd,be96,9b04,c127,0421)   /* ~0.99856261824416247814377445291000867655384354293346 T2[ 30] hi,lo */
/* -33*/ ,LDOUBLE_HEX(3ffe,d61b,5dfe,9f9b,ce07),LDOUBLE_HEX(bfbc,8d32,b9db,3407,1837)   /* ~0.83635508982079828681687957980805947499902686104178 T1[ 31] hi,lo */
/* -33*/ ,LDOUBLE_HEX(3ffe,ffa4,9183,66a2,2d19),LDOUBLE_HEX(bfbc,d9f7,dc7c,8136,fec5)   /* ~0.99860486467390137535005481539407412583386758342385 T2[ 31] hi,lo */
/* -32*/ ,LDOUBLE_HEX(3ffe,d744,fcca,d69d,6af4),LDOUBLE_HEX(3fbc,e69a,2ee6,40b4,ff78)   /* ~0.84089641525371454301891749327779734812793321907520 T1[ 32] hi,lo */
/* -32*/ ,LDOUBLE_HEX(3ffe,ffa7,5652,1c8d,aed2),LDOUBLE_HEX(bfbd,c18b,c96e,08d6,74ee)   /* ~0.99864711289097017360861380241843221483577508479356 T2[ 32] hi,lo */
/* -31*/ ,LDOUBLE_HEX(3ffe,d870,394c,6db3,2c84),LDOUBLE_HEX(3fbc,8559,bf8d,ed94,1cbc)   /* ~0.84546239963465259097986914671452041147858835756778 T1[ 33] hi,lo */
/* -31*/ ,LDOUBLE_HEX(3ffe,ffaa,1b28,7fbf,1dbe),LDOUBLE_HEX(3fbd,8eb4,6838,bb56,a087)   /* ~0.99868936289544448989524000337425491125031840056180 T2[ 33] hi,lo */
/* -30*/ ,LDOUBLE_HEX(3ffe,d99d,15c2,78af,d7b6),LDOUBLE_HEX(bfb7,bc61,09ae,0f6a,2a20)   /* ~0.85005317685926173475099293375123465921205934137105 T1[ 34] hi,lo */
/* -30*/ ,LDOUBLE_HEX(3ffe,ffac,e006,904b,c2e1),LDOUBLE_HEX(3fbc,b14f,99ba,f83d,8dd4)   /* ~0.99873161468739994449253863373350270649098092690110 T2[ 34] hi,lo */
/* -29*/ ,LDOUBLE_HEX(3ffe,dacb,946f,2ac9,cc72),LDOUBLE_HEX(bfbc,efdd,dd36,f187,1d1c)   /* ~0.85466888155023141356459587258953547461715061217546 T1[ 35] hi,lo */
/* -29*/ ,LDOUBLE_HEX(3ffe,ffaf,a4ec,4e48,e778),LDOUBLE_HEX(bfbd,8671,6e2d,8219,1444)   /* ~0.99877386826691216077309110055182372889248654246330 T2[ 35] hi,lo */
/* -28*/ ,LDOUBLE_HEX(3ffe,dbfb,b797,daf2,3755),LDOUBLE_HEX(3fbc,f610,356a,78a6,a991)   /* ~0.85930964906123895780164570479264796176721574738621 T1[ 36] hi,lo */
/* -28*/ ,LDOUBLE_HEX(3ffe,ffb2,69d9,b9cb,d4fa),LDOUBLE_HEX(3fbd,d84d,2ee6,e8d4,d3b8)   /* ~0.99881612363405676525366511109282896541117224842309 T2[ 36] hi,lo */
/* -27*/ ,LDOUBLE_HEX(3ffe,dd2d,8185,0832,4c20),LDOUBLE_HEX(3fbd,cb3c,6af5,b47f,2972)   /* ~0.86397561548091878110000596535655859042890369892120 T1[ 37] hi,lo */
/* -27*/ ,LDOUBLE_HEX(3ffe,ffb5,2ece,d2e9,d51d),LDOUBLE_HEX(bfbd,994c,52bb,c51e,a42f)   /* ~0.99885838078890938786626521594946837012685136869549 T2[ 37] hi,lo */
/* -26*/ ,LDOUBLE_HEX(3ffe,de60,f482,5e0e,9124),LDOUBLE_HEX(bfbc,8be1,7498,5ee6,5e9c)   /* ~0.86866691763685312450450626275966214961954392492771 T1[ 38] hi,lo */
/* -26*/ ,LDOUBLE_HEX(3ffe,ffb7,f3cb,99b8,31cd),LDOUBLE_HEX(3fbb,8b13,c75c,9907,9d08)   /* ~0.99890063973154566147024183142555386893945978954434 T2[ 38] hi,lo */
/* -25*/ ,LDOUBLE_HEX(3ffe,df96,12de,b8f0,4420),LDOUBLE_HEX(3fbd,8d70,2518,e344,9fa0)   /* ~0.87338369309958447002373294054677899111993610858917 T1[ 39] hi,lo */
/* -25*/ ,LDOUBLE_HEX(3ffe,ffba,b8d0,0e4c,3535),LDOUBLE_HEX(3fba,bff0,5a1f,2fe7,8cc7)   /* ~0.99894290046204122234018221715423635487240971997380 T2[ 39] hi,lo */
/* -24*/ ,LDOUBLE_HEX(3ffe,e0cc,deec,2a94,e111),LDOUBLE_HEX(3fb9,cb12,a091,ba66,7944)   /* ~0.87812608018664974155473648687042498295340919867157 T1[ 40] hi,lo */
/* -24*/ ,LDOUBLE_HEX(3ffe,ffbd,7ddc,30bb,29b9),LDOUBLE_HEX(3fbc,c13b,06cc,24fb,adcf)   /* ~0.99898516298047170984064982435235435787035385146737 T2[ 40] hi,lo */
/* -23*/ ,LDOUBLE_HEX(3ffe,e205,5aff,fe83,d369),LDOUBLE_HEX(bfbd,b207,df0e,7d75,ac99)   /* ~0.88289421796663641054054086354518915413791546598076 T1[ 41] hi,lo */
/* -23*/ ,LDOUBLE_HEX(3ffe,ffc0,42f0,011a,59f9),LDOUBLE_HEX(3fbb,9baa,d16f,ffef,fe31)   /* ~0.99902742728691276658881462169325970990030327811837 T2[ 41] hi,lo */
/* -22*/ ,LDOUBLE_HEX(3ffe,e33f,8972,be8a,5a51),LDOUBLE_HEX(3fba,9bfe,9079,5980,eecf)   /* ~0.88768824626326062627321503151200943193543935194611 T1[ 42] hi,lo */
/* -22*/ ,LDOUBLE_HEX(3ffe,ffc3,080b,7f7f,10cf),LDOUBLE_HEX(3fbd,a67d,211d,b089,b842)   /* ~0.99906969338144003834603287805826710155088221654295 T2[ 42] hi,lo */
/* -21*/ ,LDOUBLE_HEX(3ffe,e47b,6ca0,373d,a88d),LDOUBLE_HEX(3fbd,cbc4,8805,c442,ddb5)   /* ~0.89250830565946749005053553749355899071815656498074 T1[ 43] hi,lo */
/* -21*/ ,LDOUBLE_HEX(3ffe,ffc5,cd2e,abfe,9952),LDOUBLE_HEX(bfbd,d27c,2e60,9c60,e59f)   /* ~0.99911196126412917418047748840947974713344592601060 T2[ 43] hi,lo */
/* -20*/ ,LDOUBLE_HEX(3ffe,e5b9,06e7,7c83,48a8),LDOUBLE_HEX(3fbb,f2f4,7a52,76dd,8765)   /* ~0.89735453750155359320742065598608405707636848092079 T1[ 44] hi,lo */
/* -20*/ ,LDOUBLE_HEX(3ffe,ffc8,9259,86ae,3ed1),LDOUBLE_HEX(bfbd,e1f3,4d80,330a,7c7d)   /* ~0.99915423093505582619608743066841327618021750822663 T2[ 44] hi,lo */
/* -19*/ ,LDOUBLE_HEX(3ffe,e6f8,5aaa,ee1f,ce22),LDOUBLE_HEX(3fbd,f895,8fac,51be,5160)   /* ~0.90222708390331194012751875321853844980068970471620 T1[ 45] hi,lo */
/* -19*/ ,LDOUBLE_HEX(3ffe,ffcb,578c,0fa3,4cd8),LDOUBLE_HEX(bfba,85a0,3df0,1bbe,480e)   /* ~0.99919650239429564980361830883737184194615110754966 T2[ 45] hi,lo */
/* -18*/ ,LDOUBLE_HEX(3ffe,e839,6a50,3c4b,dc68),LDOUBLE_HEX(3fbd,f22f,21a1,58e1,8fbc)   /* ~0.90712608775019937809927506044971323717618361115455 T1[ 46] hi,lo */
/* -18*/ ,LDOUBLE_HEX(3ffe,ffce,1cc6,46f3,0f2f),LDOUBLE_HEX(bfbd,d11c,79ec,fe71,2113)   /* ~0.99923877564192430366643224437517289970855927094817 T2[ 46] hi,lo */
/* -17*/ ,LDOUBLE_HEX(3ffe,e97c,3840,6c4f,8c57),LDOUBLE_HEX(bfba,f6e3,3b0a,efed,259d)   /* ~0.91205169270352665550133028360413334212353220209479 T1[ 47] hi,lo */
/* -17*/ ,LDOUBLE_HEX(3ffe,ffd0,e208,2cb2,d1d6),LDOUBLE_HEX(3fbd,822a,581b,5246,bdac)   /* ~0.99928105067801744948365744170004631996562238782644 T2[ 47] hi,lo */
/* -16*/ ,LDOUBLE_HEX(3ffe,eac0,c6e7,dd24,392f),LDOUBLE_HEX(bfbc,bf4a,2932,3e46,ac15)   /* ~0.91700404320467123175366838672140090693574165925383 T1[ 48] hi,lo */
/* -16*/ ,LDOUBLE_HEX(3ffe,ffd3,a751,c0f7,e10c),LDOUBLE_HEX(bfbc,b118,1d47,fb41,07e6)   /* ~0.99932332750265075236965894855956094033899717032909 T2[ 48] hi,lo */
/* -15*/ ,LDOUBLE_HEX(3ffe,ec07,18b6,4c1c,bddc),LDOUBLE_HEX(3fbc,9f3a,0910,0bf0,97d9)   /* ~0.92198328447931296252514099531794045105925761163234 T1[ 49] hi,lo */
/* -15*/ ,LDOUBLE_HEX(3ffe,ffd6,6ca3,03d7,8947),LDOUBLE_HEX(3fbc,eeaf,9a2f,a3c1,3789)   /* ~0.99936560611589988042035778703642279197083553299307 T2[ 49] hi,lo */
/* -14*/ ,LDOUBLE_HEX(3ffe,ed4f,301e,d994,2b84),LDOUBLE_HEX(3fbd,c01a,5b6d,4c97,f624)   /* ~0.92698956254169278417588684337857785067171789705753 T1[ 50] hi,lo */
/* -14*/ ,LDOUBLE_HEX(3ffe,ffd9,31fb,f567,173b),LDOUBLE_HEX(3fbd,e5f8,d3c8,1e5c,b0b2)   /* ~0.99940788651784050509270171391840165142639307305216 T2[ 50] hi,lo */
/* -13*/ ,LDOUBLE_HEX(3ffe,ee99,0f98,0da3,025b),LDOUBLE_HEX(3fbd,95de,3c06,30a3,9321)   /* ~0.93202302419889452238894664293233915941527811810374 T1[ 51] hi,lo */
/* -13*/ ,LDOUBLE_HEX(3ffe,ffdb,f75c,95bb,d7d7),LDOUBLE_HEX(bfbd,fbb3,2c62,3fb8,60df)   /* ~0.99945016870854830104203489482550537559291115030646 T2[ 51] hi,lo */
/* -12*/ ,LDOUBLE_HEX(3ffe,efe4,b99b,dcda,f5cb),LDOUBLE_HEX(3fbd,8cac,39ed,291b,7226)   /* ~0.93708381705514995065010516728243317174928961321711 T1[ 52] hi,lo */
/* -12*/ ,LDOUBLE_HEX(3ffe,ffde,bcc4,e4eb,1842),LDOUBLE_HEX(bfbd,fe9c,0372,21ba,c7cf)   /* ~0.99949245268809894595946757833715423657849896699190 T2[ 52] hi,lo */
/* -11*/ ,LDOUBLE_HEX(3ffe,f132,30a7,ad09,4509),LDOUBLE_HEX(3fbc,ec3f,42f5,b4c8,cfcf)   /* ~0.94217208951616722483130353382740906909020850434899 T1[ 53] hi,lo */
/* -11*/ ,LDOUBLE_HEX(3ffe,ffe1,8234,e30a,25e1),LDOUBLE_HEX(bfbd,d44f,a015,3fba,570a)   /* ~0.99953473845656812089713674773783225191436940804123 T2[ 53] hi,lo */
/* -10*/ ,LDOUBLE_HEX(3ffe,f281,773c,59ff,b13a),LDOUBLE_HEX(bfbb,bb3f,ab19,b85c,2da7)   /* ~0.94728799079348282067506575865323270591034088283777 T1[ 54] hi,lo */
/* -10*/ ,LDOUBLE_HEX(3ffe,ffe4,47ac,902e,4e53),LDOUBLE_HEX(bfb9,db28,c775,3bab,3e53)   /* ~0.99957702601403151005136568651998629775334848091006 T2[ 54] hi,lo */
/*  -9*/ ,LDOUBLE_HEX(3ffe,f3d2,8fde,3a64,1a5b),LDOUBLE_HEX(bfbd,b74d,7cdc,8798,a758)   /* ~0.95243167090883710184474539306442864017299143597483 T1[ 55] hi,lo */
/*  -9*/ ,LDOUBLE_HEX(3ffe,ffe7,0d2b,ec6c,df73),LDOUBLE_HEX(bfba,cdc3,3c89,4bd3,b790)   /* ~0.99961931536056480087108419563257655227062059566378 T2[ 55] hi,lo */
/*  -8*/ ,LDOUBLE_HEX(3ffe,f525,7d15,2486,cc2c),LDOUBLE_HEX(3fbd,f73a,18f5,db30,1f87)   /* ~0.95760328069857364691012946433801289458642713725566 T1[ 56] hi,lo */
/*  -8*/ ,LDOUBLE_HEX(3ffe,ffe9,d2b2,f7db,2756),LDOUBLE_HEX(bfbc,8838,b5d7,21d6,f1ce)   /* ~0.99966160649624368394940837623252605226298328489065 T2[ 56] hi,lo */
/*  -7*/ ,LDOUBLE_HEX(3ffe,f67a,416c,733f,846e),LDOUBLE_HEX(bfbd,fced,046b,6310,b9e0)   /* ~0.96280297181806246450529868097234498236502986401319 T1[ 57] hi,lo */
/*  -7*/ ,LDOUBLE_HEX(3ffe,ffec,9841,b28e,744c),LDOUBLE_HEX(3fbc,97cc,e653,acf3,2bde)   /* ~0.99970389942114385307785073830899591484921984374523 T2[ 57] hi,lo */
/*  -6*/ ,LDOUBLE_HEX(3ffe,f7d0,df73,0ad1,3bb9),LDOUBLE_HEX(bfb7,b795,b494,f824,8a8b)   /* ~0.96803089674614722529933166783600029248191276565194 T1[ 58] hi,lo */
/*  -6*/ ,LDOUBLE_HEX(3ffe,ffef,5dd8,1c9c,14e1),LDOUBLE_HEX(3fbd,94dc,70ee,035c,489a)   /* ~0.99974619413534100530053030930766055917047196999192 T2[ 58] hi,lo */
/*  -5*/ ,LDOUBLE_HEX(3ffe,f929,59bb,5dd4,ba74),LDOUBLE_HEX(3fbc,d2df,86c7,21a9,8d5b)   /* ~0.97328720878961664316093849169675422672298736870288 T1[ 59] hi,lo */
/*  -5*/ ,LDOUBLE_HEX(3ffe,fff2,2376,3619,57dc),LDOUBLE_HEX(bfbc,873d,93df,ba07,4a3e)   /* ~0.99978849063891084085996252550643248468986712396144 T2[ 59] hi,lo */
/*  -4*/ ,LDOUBLE_HEX(3ffe,fa83,b2db,722a,033a),LDOUBLE_HEX(3fbd,f84b,7628,62ba,ff99)   /* ~0.97857206208770013448287189428143051372899208217859 T1[ 60] hi,lo */
/*  -4*/ ,LDOUBLE_HEX(3ffe,fff4,e91b,ff1b,8c3e),LDOUBLE_HEX(bfbd,ef98,e3e2,81af,6b64)   /* ~0.99983078893192906314284912339118704949214588850736 T2[ 60] hi,lo */
/*  -3*/ ,LDOUBLE_HEX(3ffe,fbdf,ed6c,e5f0,9c49),LDOUBLE_HEX(bfbd,c4b4,018d,426a,3a32)   /* ~0.98388561161658788907719355720793430464254925027489 T1[ 61] hi,lo */
/*  -3*/ ,LDOUBLE_HEX(3ffe,fff7,aec9,77b8,0143),LDOUBLE_HEX(3fb9,87ed,a7ba,da2f,9976)   /* ~0.99987308901447137873428824828003769198403460904955 T2[ 61] hi,lo */
/*  -2*/ ,LDOUBLE_HEX(3ffe,fd3e,0c0c,f486,c175),LDOUBLE_HEX(bfbd,f581,8b4d,9c3e,23fa)   /* ~0.98922801319397548415511892128115789546427549794316 T1[ 62] hi,lo */
/*  -2*/ ,LDOUBLE_HEX(3ffe,fffa,747e,a004,0664),LDOUBLE_HEX(3fbc,8e3e,4bde,4901,6016)   /* ~0.99991539088661349752619467157188637429499067366123 T2[ 62] hi,lo */
/*  -1*/ ,LDOUBLE_HEX(3ffe,fe9e,115c,7b8f,884c),LDOUBLE_HEX(bfbd,a45b,4cd4,30c5,a1ed)   /* ~0.99459942348363317566987962248248322794097475707530 T1[ 63] hi,lo */
/*  -1*/ ,LDOUBLE_HEX(3ffe,fffd,3a3b,7814,eb54),LDOUBLE_HEX(bfbc,ca27,58a3,c057,ba57)   /* ~0.99995769454843113255466946487359791717608459293842 T2[ 63] hi,lo */
/*   0*/ ,LDOUBLE_HEX(3fff,8000,0000,0000,0000),LDOUBLE_HEX(0000,0000,0000,0000,0000)   /*  1.0                                                  T1[ 64] hi,lo */
/*   0*/ ,LDOUBLE_HEX(3fff,8000,0000,0000,0000),LDOUBLE_HEX(0000,0000,0000,0000,0000)   /*  1.0                                                  T2[ 64] hi,lo */
/*   1*/ ,LDOUBLE_HEX(3fff,80b1,ed4f,d999,ab6c),LDOUBLE_HEX(3fbd,94cd,5c66,db9b,f480)   /* ~1.00542990111280282133562891111466797156026586890220 T1[ 65] hi,lo */
/*   1*/ ,LDOUBLE_HEX(3fff,8001,62e6,1bed,4a49),LDOUBLE_HEX(bfbc,bd9e,8f2d,ce5c,6131)   /* ~1.00004230724139581934929027484670882586215157061815 T2[ 65] hi,lo */
/*   2*/ ,LDOUBLE_HEX(3fff,8164,d1f3,bc03,0773),LDOUBLE_HEX(3fbe,f7ca,ca4f,7a29,bde9)   /* ~1.01088928605170045996793776188482638644927646964788 T1[ 66] hi,lo */
/*   2*/ ,LDOUBLE_HEX(3fff,8002,c5d0,0fdc,fcb7),LDOUBLE_HEX(bfbe,9353,2b4e,7f6e,83c2)   /* ~1.00008461627269431323383058751730345647956710308790 T2[ 66] hi,lo */
/*   3*/ ,LDOUBLE_HEX(3fff,8218,af43,73fc,25ec),LDOUBLE_HEX(bfbe,c706,5df2,5b81,9130)   /* ~1.01637831491095303798263826955405875196447595953941 T1[ 67] hi,lo */
/*   3*/ ,LDOUBLE_HEX(3fff,8004,28bd,dbd9,bf99),LDOUBLE_HEX(3fbb,d289,3e85,affc,a646)   /* ~1.00012692709397120742909753632332581219088751822710 T2[ 67] hi,lo */
/*   4*/ ,LDOUBLE_HEX(3fff,82cd,8698,ac2b,a1d7),LDOUBLE_HEX(3fbd,f8a9,1d6d,1948,2ffd)   /* ~1.02189714865411667820815216911611855721275787800550 T1[ 68] hi,lo */
/*   4*/ ,LDOUBLE_HEX(3fff,8005,8baf,7fee,3b5d),LDOUBLE_HEX(3fbc,e38c,59c7,2a4e,5c9a)   /* ~1.00016923970530223107159445428138155875785741955041 T2[ 68] hi,lo */
/*   5*/ ,LDOUBLE_HEX(3fff,8383,594e,efb6,ee37),LDOUBLE_HEX(bfbc,eff1,589e,1360,4be1)   /* ~1.02744594911876369655156365778836402569140773266553 T1[ 69] hi,lo */
/*   5*/ ,LDOUBLE_HEX(3fff,8006,eea4,fc25,188e),LDOUBLE_HEX(bfbe,be51,bfd6,da9c,50ec)   /* ~1.00021155410676311644201097461603922056383453309535 T2[ 69] hi,lo */
/*   6*/ ,LDOUBLE_HEX(3fff,843a,28c3,acde,4046),LDOUBLE_HEX(3fbc,d7c9,7650,9fe8,ac10)   /* ~1.03302487902122842248868461734190304923686198890209 T1[ 70] hi,lo */
/*   6*/ ,LDOUBLE_HEX(3fff,8008,519e,5088,ffd3),LDOUBLE_HEX(bfbc,ddbf,fb2a,c38e,3196)   /* ~1.00025387029842959885680281351127973721304442733526 T2[ 70] hi,lo */
/*   7*/ ,LDOUBLE_HEX(3fff,84f1,f656,379c,1a29),LDOUBLE_HEX(3fbb,f030,62c2,6b5b,a5d1)   /* ~1.03863410196137879060607894787793270552356261759996 T1[ 71] hi,lo */
/*   7*/ ,LDOUBLE_HEX(3fff,8009,b49b,7d24,99f3),LDOUBLE_HEX(bfbe,8364,adc5,5d2d,55af)   /* ~1.00029618828037741710187263910469823713356163352727 T2[ 71] hi,lo */
/*   8*/ ,LDOUBLE_HEX(3fff,85aa,c367,cc48,7b15),LDOUBLE_HEX(bfbd,e8da,91cf,7aac,f938)   /* ~1.04427378242741384034662083246658426105568651109933 T1[ 72] hi,lo */
/*   8*/ ,LDOUBLE_HEX(3fff,800b,179c,8202,8fd1),LDOUBLE_HEX(bfbe,d743,563a,a3ce,1a20)   /* ~1.00033850805268231299888920249330226397432852536439 T2[ 72] hi,lo */
/*   9*/ ,LDOUBLE_HEX(3fff,8664,915b,923f,ba04),LDOUBLE_HEX(bfbd,91f4,8ed8,4742,eaa7)   /* ~1.04994408580068726609749174549790495802881196141242 T1[ 73] hi,lo */
/*   9*/ ,LDOUBLE_HEX(3fff,800c,7aa1,5f2d,8a6d),LDOUBLE_HEX(3fbe,b45d,df92,07f3,cad6)   /* ~1.00038082961542003151370755498206222000590059906244 T2[ 73] hi,lo */
/*  10*/ ,LDOUBLE_HEX(3fff,871f,6196,9e8d,1010),LDOUBLE_HEX(3fbd,e85c,9f15,ed4a,a559)   /* ~1.05564517836055715878373906235765389283187687397003 T1[ 74] hi,lo */
/*  10*/ ,LDOUBLE_HEX(3fff,800d,ddaa,14b0,32e8),LDOUBLE_HEX(bfbe,dbf5,ee63,a474,2c9a)   /* ~1.00042315296866632119004991707811313972342759370803 T2[ 74] hi,lo */
/*  11*/ ,LDOUBLE_HEX(3fff,87db,357f,f698,d792),LDOUBLE_HEX(bfbe,df6e,2275,ebd9,aeb2)   /* ~1.06137722728926208099788086602899284116574563086032 T1[ 75] hi,lo */
/*  11*/ ,LDOUBLE_HEX(3fff,800f,40b6,a295,327b),LDOUBLE_HEX(3fbe,9a11,8af4,78f3,787d)   /* ~1.00046547811249693317372372325380069923994597047567 T2[ 75] hi,lo */
/*  12*/ ,LDOUBLE_HEX(3fff,8898,0e80,92da,8527),LDOUBLE_HEX(3fbe,bbf1,aed9,318c,eac6)   /* ~1.06714040067682361812972241521535465835768263787031 T1[ 76] hi,lo */
/*  12*/ ,LDOUBLE_HEX(3fff,8010,a3c7,08e7,3282),LDOUBLE_HEX(3fbd,ae5b,58b5,4705,681e)   /* ~1.00050780504698762240524401168073609369457699358463 T2[ 76] hi,lo */
/*  13*/ ,LDOUBLE_HEX(3fff,8955,ee03,618e,5fdd),LDOUBLE_HEX(bfbe,d452,cdb2,971d,08d9)   /* ~1.07293486752597555142999669053338607227487955242395 T1[ 77] hi,lo */
/*  13*/ ,LDOUBLE_HEX(3fff,8012,06db,47b0,dc73),LDOUBLE_HEX(3fbe,bb24,d6a1,6ac0,de44)   /* ~1.00055013377221414664405146899284204664581920951604 T2[ 77] hi,lo */
/*  14*/ ,LDOUBLE_HEX(3fff,8a14,d575,496e,fd9a),LDOUBLE_HEX(3fbb,80ca,1d92,c368,0c22)   /* ~1.07876079775711979373727100739444040300440974533557 T1[ 78] hi,lo */
/*  14*/ ,LDOUBLE_HEX(3fff,8013,69f3,5efc,d9e4),LDOUBLE_HEX(bfbe,ba91,a929,afe6,0b3f)   /* ~1.00059246428825226711903373377765547047602012753486 T2[ 78] hi,lo */
/*  15*/ ,LDOUBLE_HEX(3fff,8ad4,c645,2c72,8924),LDOUBLE_HEX(3fba,d573,dd56,13bf,92a3)   /* ~1.08461836221330923781328015031988343253033235669136 T1[ 79] hi,lo */
/*  15*/ ,LDOUBLE_HEX(3fff,8014,cd0f,4ed5,d485),LDOUBLE_HEX(3fbe,95d6,5430,ec64,68fa)   /* ~1.00063479659517774787800409308502480598690453916788 T2[ 79] hi,lo */
/*  16*/ ,LDOUBLE_HEX(3fff,8b95,c1e3,ea8b,d6e7),LDOUBLE_HEX(bfba,8373,af14,eb58,6dfd)   /* ~1.09050773266525765920875040704274283598351757973432 T1[ 80] hi,lo */
/*  16*/ ,LDOUBLE_HEX(3fff,8016,302f,1746,7628),LDOUBLE_HEX(3fbd,da43,7f91,3447,4021)   /* ~1.00067713069306635665506322041551356960553675889968 T2[ 80] hi,lo */
/*  17*/ ,LDOUBLE_HEX(3fff,8c57,c9c4,646f,4dde),LDOUBLE_HEX(bfba,8f46,5c3d,afa3,6840)   /* ~1.09642908181637682338803452264386351089342497289180 T1[ 81] hi,lo */
/*  17*/ ,LDOUBLE_HEX(3fff,8017,9352,b859,68ba),LDOUBLE_HEX(bfbd,a49b,73de,ad19,46f2)   /* ~1.00071946658199386411165765498054724957910366356372 T2[ 81] hi,lo */
/*  18*/ ,LDOUBLE_HEX(3fff,8d1a,df5b,7e5b,a9e6),LDOUBLE_HEX(bfbe,9670,96d2,e37c,a594)   /* ~1.10238258330784094358827107651421783884870819747447 T1[ 82] hi,lo */
/*  18*/ ,LDOUBLE_HEX(3fff,8018,f67a,3219,5645),LDOUBLE_HEX(3fbd,b854,8436,4510,af0d)   /* ~1.00076180426203604405342023619951419277640525251626 T2[ 82] hi,lo */
/*  19*/ ,LDOUBLE_HEX(3fff,8ddf,0420,22e6,9cd6),LDOUBLE_HEX(bfbe,e18d,4bbd,81ca,0653)   /* ~1.10836841172367863805718612990602878198842518031597 T1[ 83] hi,lo */
/*  19*/ ,LDOUBLE_HEX(3fff,801a,59a5,8490,e8f3),LDOUBLE_HEX(bfbd,bfde,9abb,c89b,4566)   /* ~1.00080414373326867375543075544541693489009048789739 T2[ 83] hi,lo */
/*  20*/ ,LDOUBLE_HEX(3fff,8ea4,398b,45cd,53c0),LDOUBLE_HEX(3fbd,b700,5132,1e0f,5317)   /* ~1.11438674259589253628943694707231770735234022140502 T1[ 84] hi,lo */
/*  20*/ ,LDOUBLE_HEX(3fff,801b,bcd4,afca,cb09),LDOUBLE_HEX(bfbc,ee2b,3ca1,60ce,c87d)   /* ~1.00084648499576753342011486980211998343293089419603 T2[ 84] hi,lo */
/*  21*/ ,LDOUBLE_HEX(3fff,8f6a,8117,e6c8,e5c4),LDOUBLE_HEX(3fbb,cffb,0890,e8f2,826a)   /* ~1.12043775240960668442349867923724104912253096699714 T1[ 85] hi,lo */
/*  21*/ ,LDOUBLE_HEX(3fff,801d,2007,b3d1,a6eb),LDOUBLE_HEX(3fbc,e24f,a3f1,8819,c368)   /* ~1.00088882804960840661092497105855159134080167859792 T2[ 85] hi,lo */
/*  22*/ ,LDOUBLE_HEX(3fff,9031,dc43,1466,b1dc),LDOUBLE_HEX(3fbe,eeb0,2950,929d,0fc5)   /* ~1.12652161860824189974425446614247903198702260851860 T1[ 86] hi,lo */
/*  22*/ ,LDOUBLE_HEX(3fff,801e,833e,90b0,271b),LDOUBLE_HEX(bfbd,df0c,d686,c810,ce9e)   /* ~1.00093117289486708014391996846015331357193645089864 T2[ 86] hi,lo */
/*  23*/ ,LDOUBLE_HEX(3fff,90fa,4c8b,eee4,b12b),LDOUBLE_HEX(bfbe,d02d,6d6b,424b,49e1)   /* ~1.13263851959871922803115701361420519788225647062063 T1[ 87] hi,lo */
/*  23*/ ,LDOUBLE_HEX(3fff,801f,e679,4670,f637),LDOUBLE_HEX(bfbd,f4b1,1157,9a3c,c09c)   /* ~1.00097351953161934387092485421177912030543666332960 T2[ 87] hi,lo */
/*  24*/ ,LDOUBLE_HEX(3fff,91c3,d373,ab11,c336),LDOUBLE_HEX(3fbb,fd6d,8e0a,e5ac,9d82)   /* ~1.13878863475669165369712210189589995934511534869670 T1[ 88] hi,lo */
/*  24*/ ,LDOUBLE_HEX(3fff,8021,49b7,d51e,befb),LDOUBLE_HEX(3fbe,f7b7,5b79,1115,d652)   /* ~1.00101586795994099089637113797479628374276217073202 T2[ 88] hi,lo */
/*  25*/ ,LDOUBLE_HEX(3fff,928e,727d,9531,f9ac),LDOUBLE_HEX(3fbc,aadf,7a7a,5204,6a72)   /* ~1.14497214443180421938536794890239889355143532156944 T1[ 89] hi,lo */
/*  25*/ ,LDOUBLE_HEX(3fff,8022,acfa,3cc4,2c43),LDOUBLE_HEX(bfbb,9803,679c,fccf,bd52)   /* ~1.00105821817990781779413728136418626490922179073095 T2[ 89] hi,lo */
/*  26*/ ,LDOUBLE_HEX(3fff,935a,2b2f,13e6,e92c),LDOUBLE_HEX(bfbd,b319,afc5,89b6,c463)   /* ~1.15118922995298270583672262112884254747768864035606 T1[ 90] hi,lo */
/*  26*/ ,LDOUBLE_HEX(3fff,8024,1040,7d6b,e905),LDOUBLE_HEX(bfbd,d992,0edb,1ae2,6862)   /* ~1.00110057019159562395702739445724205324950162321329 T2[ 90] hi,lo */
/*  27*/ ,LDOUBLE_HEX(3fff,9426,ff0f,ab1c,04b6),LDOUBLE_HEX(3fbe,f15c,f03c,a096,7fdb)   /* ~1.15744007363375102956197515435832201546872965991497 T1[ 91] hi,lo */
/*  27*/ ,LDOUBLE_HEX(3fff,8025,738a,9720,a056),LDOUBLE_HEX(3fbd,8031,e4e2,c90f,b56a)   /* ~1.00114292399508021213887232203632038363139145076274 T2[ 91] hi,lo */
/*  28*/ ,LDOUBLE_HEX(3fff,94f4,efa8,fef7,0961),LDOUBLE_HEX(3fbd,ba2b,eb44,9547,7951)   /* ~1.16372485877757751379386191858955612588033545762300 T1[ 92] hi,lo */
/*  28*/ ,LDOUBLE_HEX(3fff,8026,d6d8,89ec,fd6a),LDOUBLE_HEX(bfbe,8df6,8809,7e58,ba93)   /* ~1.00118527959043738845452964358884173634578473865985 T2[ 92] hi,lo */
/*  29*/ ,LDOUBLE_HEX(3fff,95c3,fe86,d6cc,7fef),LDOUBLE_HEX(bfbb,adcd,6381,aa3b,dde8)   /* ~1.17004376968325018808485954435738563006452750414609 T1[ 93] hi,lo */
/*  29*/ ,LDOUBLE_HEX(3fff,8028,3a2a,55db,ab90),LDOUBLE_HEX(bfbc,bd1b,9f76,51ee,96be)   /* ~1.00122763697774296194620280431308856350369751453399 T2[ 93] hi,lo */
/*  30*/ ,LDOUBLE_HEX(3fff,9694,2d37,2018,5a00),LDOUBLE_HEX(3fbe,91d5,36d0,7538,458a)   /* ~1.17639699165028127625376441756088752299547195434570 T1[ 94] hi,lo */
/*  30*/ ,LDOUBLE_HEX(3fff,8029,9d7f,faf7,5637),LDOUBLE_HEX(bfbe,d236,fdc3,f478,d426)   /* ~1.00126999615707274512554220136095750603999476879835 T2[ 94] hi,lo */
/*  31*/ ,LDOUBLE_HEX(3fff,9765,7d49,f17a,b08e),LDOUBLE_HEX(3fbe,a0f4,5d52,3833,af61)   /* ~1.18278471098434102989037375319725242661661468446254 T1[ 95] hi,lo */
/*  31*/ ,LDOUBLE_HEX(3fff,802b,00d9,794a,a8e9),LDOUBLE_HEX(3fbe,856f,102a,5311,acd9)   /* ~1.00131235712850255343154409759520717670966405421495 T2[ 95] hi,lo */
/*  32*/ ,LDOUBLE_HEX(3fff,9837,f051,8db8,a96f),LDOUBLE_HEX(3fbe,8d5a,4630,5c85,eded)   /* ~1.18920711500272106668756738612202639160386752337217 T1[ 96] hi,lo */
/*  32*/ ,LDOUBLE_HEX(3fff,802c,6436,d0e0,4f51),LDOUBLE_HEX(bfb6,e62d,6b30,d098,6397)   /* ~1.00135471989210820588107192508076082049228716641664 T2[ 96] hi,lo */
/*  33*/ ,LDOUBLE_HEX(3fff,990b,87e2,66c1,89aa),LDOUBLE_HEX(bfbd,c61c,79fe,e0f2,443a)   /* ~1.19566439203982737460481289293312556765158660709857 T1[ 97] hi,lo */
/*  33*/ ,LDOUBLE_HEX(3fff,802d,c798,01c2,f534),LDOUBLE_HEX(3fbe,ee14,9ecb,14b8,e2d1)   /* ~1.00139708444796552430991476434485321078682318329811 T2[ 97] hi,lo */
/*  34*/ ,LDOUBLE_HEX(3fff,99e0,4593,20b7,fa65),LDOUBLE_HEX(bfbc,de7b,c9a6,5a50,1a8c)   /* ~1.20215673145270314210817513833617908858286682516336 T1[ 98] hi,lo */
/*  34*/ ,LDOUBLE_HEX(3fff,802f,2afd,0bfd,4678),LDOUBLE_HEX(bfbc,b751,8633,d8a0,aeca)   /* ~1.00143945079615033413172886511688375321682542562484 T2[ 98] hi,lo */
/*  35*/ ,LDOUBLE_HEX(3fff,9ab6,2afc,94ff,864a),LDOUBLE_HEX(3fbd,c468,ec6e,75e7,1adb)   /* ~1.20868432362658157733399655331396616020356304943561 T1[ 99] hi,lo */
/*  35*/ ,LDOUBLE_HEX(3fff,8030,8e65,ef99,ef1d),LDOUBLE_HEX(3fbd,b09f,c1fc,23b2,9fab)   /* ~1.00148181893673846357909612558856338182522449642419 T2[ 99] hi,lo */
/*  36*/ ,LDOUBLE_HEX(3fff,9b8d,39b9,d54e,5539),LDOUBLE_HEX(bfbe,baaf,d0ba,b867,81c2)   /* ~1.21524735998046887815605271443430979161348659545183 T1[100] hi,lo */
/*  36*/ ,LDOUBLE_HEX(3fff,8031,f1d2,aca3,9b44),LDOUBLE_HEX(bfbe,a4c4,911a,dfc4,d295)   /* ~1.00152418886980574446246561315376766287954524159431 T2[100] hi,lo */
/*  37*/ ,LDOUBLE_HEX(3fff,9c65,7368,2ec3,2c2d),LDOUBLE_HEX(3fbe,9cb0,d9be,d0c8,53bd)   /* ~1.22184603297275751687071126960759670510014984756708 T1[101] hi,lo */
/*  37*/ ,LDOUBLE_HEX(3fff,8033,5543,4324,f728),LDOUBLE_HEX(3fbe,e415,5eab,b7e1,f081)   /* ~1.00156656059542801141121204366868369106668978929519 T2[101] hi,lo */
/*  38*/ ,LDOUBLE_HEX(3fff,9d3e,d9a7,2cff,b751),LDOUBLE_HEX(bfbd,86da,cc3e,bc59,93d4)   /* ~1.22848053610687000570828725232175315795757342129945 T1[102] hi,lo */
/*  38*/ ,LDOUBLE_HEX(3fff,8034,b8b7,b328,af26),LDOUBLE_HEX(3fbc,ac57,aea7,6d92,1bde)   /* ~1.00160893411368110274099751944021363669889979064464 T2[102] hi,lo */
/*  39*/ ,LDOUBLE_HEX(3fff,9e19,6e18,9d47,2420),LDOUBLE_HEX(3fb7,f914,5ac7,9bba,f035)   /* ~1.23515106393693330569250043993179133394733071327209 T1[103] hi,lo */
/*  39*/ ,LDOUBLE_HEX(3fff,8036,1c2f,fcb9,6fb5),LDOUBLE_HEX(bfbd,cfaa,e2c6,1949,76e1)   /* ~1.00165130942464085958640979123757119850779417902231 T2[103] hi,lo */
/*  40*/ ,LDOUBLE_HEX(3fff,9ef5,3260,91a1,11ae),LDOUBLE_HEX(bfbe,bedd,c1ec,288c,045d)   /* ~1.24185781207348404863409496723392066996893845498561 T1[104] hi,lo */
/*  40*/ ,LDOUBLE_HEX(3fff,8037,7fac,1fe1,e56a),LDOUBLE_HEX(3fbe,c39a,17ff,af9f,8d06)   /* ~1.00169368652838312633464312728648337724735029041767 T2[104] hi,lo */
/*  41*/ ,LDOUBLE_HEX(3fff,9fd2,2825,6400,dd06),LDOUBLE_HEX(bfba,8fe5,5be7,cd04,73e3)   /* ~1.24860097718920473662367054412669631346943788230419 T1[105] hi,lo */
/*  41*/ ,LDOUBLE_HEX(3fff,8038,e32c,1cac,bcfa),LDOUBLE_HEX(3fbd,ed51,9851,d0d3,59b7)   /* ~1.00173606542498375084233874776629136249539442360401 T2[105] hi,lo */
/*  42*/ ,LDOUBLE_HEX(3fff,a0b0,510f,b971,4fc2),LDOUBLE_HEX(3fbc,c96e,3cf6,d87e,cd4c)   /* ~1.25538075702469108956872700932905217996449209749698 T1[106] hi,lo */
/*  42*/ ,LDOUBLE_HEX(3fff,803a,46af,f324,a335),LDOUBLE_HEX(3fbe,ac3c,8131,843a,6e57)   /* ~1.00177844611451858389348373856719831564987543970346 T2[106] hi,lo */
/*  43*/ ,LDOUBLE_HEX(3fff,a18f,aeca,8544,b6e4),LDOUBLE_HEX(bfbe,fbbc,6bef,3313,7e1e)   /* ~1.26219735039425070806731743466855277802096679806709 T1[107] hi,lo */
/*  43*/ ,LDOUBLE_HEX(3fff,803b,aa37,a354,450a),LDOUBLE_HEX(3fbe,9180,bea1,677e,0605)   /* ~1.00182082859706347963309192028447114353184588253498 T2[107] hi,lo */
/*  44*/ ,LDOUBLE_HEX(3fff,a270,4303,0c49,6819),LDOUBLE_HEX(bfbe,c90b,f620,fe60,42b1)   /* ~1.26905095719173322259699238090391304467630106955766 T1[108] hi,lo */
/*  44*/ ,LDOUBLE_HEX(3fff,803d,0dc3,2d46,4f85),LDOUBLE_HEX(3fbe,868a,dee3,72d5,ffa8)   /* ~1.00186321287269429535036341372133961158397141844034 T2[108] hi,lo */
/*  45*/ ,LDOUBLE_HEX(3fff,a352,0f68,e802,bb93),LDOUBLE_HEX(bfbe,ed0b,a6dd,6268,20c0)   /* ~1.27594177839639210043139877504003720787295605987310 T1[109] hi,lo */
/*  45*/ ,LDOUBLE_HEX(3fff,803e,7152,9105,6fd0),LDOUBLE_HEX(3fbb,f7f1,9404,a7bf,1165)   /* ~1.00190559894148689158710485713754678727127611637115 T2[109] hi,lo */
/*  46*/ ,LDOUBLE_HEX(3fff,a435,15ae,09e6,809e),LDOUBLE_HEX(3fbb,d1db,4831,781e,1eec)   /* ~1.28287001607877828072111492385687370187952183187007 T1[110] hi,lo */
/*  46*/ ,LDOUBLE_HEX(3fff,803f,d4e5,ce9c,5332),LDOUBLE_HEX(3fbc,8ac0,9487,d439,3b90)   /* ~1.00194798680351713202930918900079859668039716780185 T2[110] hi,lo */
/*  47*/ ,LDOUBLE_HEX(3fff,a519,5786,be9e,f339),LDOUBLE_HEX(3fbe,d8bc,f46f,9586,461e)   /* ~1.28983587340666581218685121656974956749763805419206 T1[111] hi,lo */
/*  47*/ ,LDOUBLE_HEX(3fff,8041,387c,e615,a710),LDOUBLE_HEX(3fbe,8e4c,645a,3a8a,6df3)   /* ~1.00199037645886088361557586523531426792033016681671 T2[111] hi,lo */
/*  48*/ ,LDOUBLE_HEX(3fff,a5fe,d6a9,b151,38ea),LDOUBLE_HEX(3fbc,e5eb,fb10,b883,80d9)   /* ~1.29683955465100966592158215906493978764046914875507 T1[112] hi,lo */
/*  48*/ ,LDOUBLE_HEX(3fff,8042,9c17,d77c,18ed),LDOUBLE_HEX(3fbe,93f9,0835,f753,878b)   /* ~1.00203276790759401653711085922182633112242911010980 T2[112] hi,lo */
/*  49*/ ,LDOUBLE_HEX(3fff,a6e5,94cf,eee8,6b1e),LDOUBLE_HEX(bfbe,c910,e561,f33b,7b4e)   /* ~1.30388126519193589861710103061653853728785179555416 T1[113] hi,lo */
/*  49*/ ,LDOUBLE_HEX(3fff,8043,ffb6,a2da,5669),LDOUBLE_HEX(3fbe,9bfa,bd93,3942,b72a)   /* ~1.00207516114979240412930644454903017503966111689805 T2[113] hi,lo */
/*  50*/ ,LDOUBLE_HEX(3fff,a7cd,93b4,e965,356a),LDOUBLE_HEX(bfbe,c274,9655,f8c1,1aa2)   /* ~1.31096121152476434196416932298490110042621381580829 T1[114] hi,lo */
/*  50*/ ,LDOUBLE_HEX(3fff,8045,6359,483b,0d42),LDOUBLE_HEX(3fbd,8ab1,e480,8fb1,67b1)   /* ~1.00211755618553192298016141226213449044735170900821 T2[114] hi,lo */
/*  51*/ ,LDOUBLE_HEX(3fff,a8b6,d516,7b32,0e09),LDOUBLE_HEX(bfbe,d0ad,37b2,7dde,6f19)   /* ~1.31807960126606399473437464253677831038658041507005 T1[115] hi,lo */
/*  51*/ ,LDOUBLE_HEX(3fff,8046,c6ff,c7a8,eb53),LDOUBLE_HEX(3fbd,cd04,05c3,2637,6be1)   /* ~1.00215995301488845282186085361431082674243953078985 T2[115] hi,lo */
/*  52*/ ,LDOUBLE_HEX(3fff,a9a1,5ab4,ea7c,0ef8),LDOUBLE_HEX(3fbe,a83c,49d8,6a63,f4e6)   /* ~1.32523664315974129459391184227001758699771016836166 T1[116] hi,lo */
/*  52*/ ,LDOUBLE_HEX(3fff,8048,2aaa,212e,9e96),LDOUBLE_HEX(bfbe,f210,9561,2774,6f44)   /* ~1.00220235163793787674761659456379447874496690928936 T2[116] hi,lo */
/*  53*/ ,LDOUBLE_HEX(3fff,aa8d,2652,ec90,7629),LDOUBLE_HEX(3fbe,ec62,0243,4ca6,7264)   /* ~1.33243254708316144930158736459091528558928985148668 T1[117] hi,lo */
/*  53*/ ,LDOUBLE_HEX(3fff,8049,8e58,54d6,d520),LDOUBLE_HEX(bfbd,fa6c,16d0,f9cb,b862)   /* ~1.00224475205475608077798632677968271309509873390197 T2[117] hi,lo */
/*  54*/ ,LDOUBLE_HEX(3fff,ab7a,39b5,a93e,d337),LDOUBLE_HEX(3fbe,cb00,4764,eb3c,00f3)   /* ~1.33966752405330300531704351696404842186893802136182 T1[118] hi,lo */
/*  54*/ ,LDOUBLE_HEX(3fff,804a,f20a,62ac,3d26),LDOUBLE_HEX(3fbd,8b1b,de3c,962c,ab75)   /* ~1.00228715426541895440297469388468698525684885680675 T2[118] hi,lo */
/*  55*/ ,LDOUBLE_HEX(3fff,ac68,96a4,be3f,e929),LDOUBLE_HEX(3fbe,bc2b,7343,bcf2,ec93)   /* ~1.34694178623294583574832722350222979912359733134508 T1[119] hi,lo */
/*  55*/ ,LDOUBLE_HEX(3fff,804c,55c0,4ab9,84fb),LDOUBLE_HEX(bfbe,ef50,9f43,e508,28a0)   /* ~1.00232955827000239036519285695803205271658953279256 T2[119] hi,lo */
/*  56*/ ,LDOUBLE_HEX(3fff,ad58,3eea,42a1,4ac6),LDOUBLE_HEX(3fbe,9301,5191,eb34,5d89)   /* ~1.35425554693689272826688518858162524338695220649242 T1[120] hi,lo */
/*  56*/ ,LDOUBLE_HEX(3fff,804d,b97a,0d09,5b0c),LDOUBLE_HEX(3fbe,d93e,3efa,3df9,fcd1)   /* ~1.00237196406858228422617762554125420138007029891014 T2[120] hi,lo */
/*  57*/ ,LDOUBLE_HEX(3fff,ae49,3452,ca35,b80e),LDOUBLE_HEX(3fbd,9637,02d3,0d44,07b1)   /* ~1.36160902063822475556963131904097963342792354524135 T1[121] hi,lo */
/*  57*/ ,LDOUBLE_HEX(3fff,804f,1d37,a9a6,6de9),LDOUBLE_HEX(bfbc,c1d5,6430,8a2f,e68c)   /* ~1.00241437166123453534217341287515523617912549525499 T2[121] hi,lo */
/*  58*/ ,LDOUBLE_HEX(3fff,af3b,78ad,690a,4375),LDOUBLE_HEX(bfbd,8367,bf8c,d132,bf35)   /* ~1.36900242297459061194351420676085240302199963480234 T1[122] hi,lo */
/*  58*/ ,LDOUBLE_HEX(3fff,8050,80f9,209b,6c3b),LDOUBLE_HEX(bfbe,a036,e851,2e22,bdf0)   /* ~1.00245678104803504577993006341429804706422146409749 T2[122] hi,lo */
/*  59*/ ,LDOUBLE_HEX(3fff,b02f,0dcb,b6e0,4584),LDOUBLE_HEX(bfbe,90a6,d5b7,91c4,cb15)   /* ~1.37643597075453010024695399415861629677237942814826 T1[123] hi,lo */
/*  59*/ ,LDOUBLE_HEX(3fff,8051,e4be,71f3,04ca),LDOUBLE_HEX(3fbc,cada,b8ef,31ad,ee75)   /* ~1.00249919222905972096722415631830926940892823040485 T2[123] hi,lo */
/*  60*/ ,LDOUBLE_HEX(3fff,b123,f581,d2ac,2590),LDOUBLE_HEX(bfbe,f05f,902d,25bd,44e3)   /* ~1.38390988196383195492356055211757848155684769153594 T1[124] hi,lo */
/*  60*/ ,LDOUBLE_HEX(3fff,8053,4887,9db7,e67d),LDOUBLE_HEX(3fbc,b8f5,8e77,78e8,f943)   /* ~1.00254160520438446969285900545187928400991950184106 T2[124] hi,lo */
/*  61*/ ,LDOUBLE_HEX(3fff,b21a,31a6,6618,fe3b),LDOUBLE_HEX(3fbe,f871,4c4e,d9a4,e410)   /* ~1.39142437577192618709722576886278488927928265184164 T1[125] hi,lo */
/*  61*/ ,LDOUBLE_HEX(3fff,8054,ac54,a3f4,c057),LDOUBLE_HEX(3fbd,ec30,98f2,c3d2,f4bb)   /* ~1.00258401997408520378140400763911088688473682850599 T2[125] hi,lo */
/*  62*/ ,LDOUBLE_HEX(3fff,b311,c412,a911,2489),LDOUBLE_HEX(3fbd,fb3c,5371,e629,4670)   /* ~1.39897967253831114018292752776417842142109293490648 T1[126] hi,lo */
/*  62*/ ,LDOUBLE_HEX(3fff,8056,1025,84b4,417a),LDOUBLE_HEX(bfbe,9388,2575,b35d,d384)   /* ~1.00262643653823783841845529440917061947402544319629 T2[126] hi,lo */
/*  63*/ ,LDOUBLE_HEX(3fff,b40a,aea2,654b,9841),LDOUBLE_HEX(bfbc,ea37,6118,3363,e500)   /* ~1.40657599381901544249601904157387366467446554452180 T1[127] hi,lo */
/*  63*/ ,LDOUBLE_HEX(3fff,8057,73fa,4001,1923),LDOUBLE_HEX(3fbe,b2a6,1fe1,eb6f,b480)   /* ~1.00266885489691829171695486300208699503855314105749 T2[127] hi,lo */
/*  64*/ ,LDOUBLE_HEX(3fff,b504,f333,f9de,6484),LDOUBLE_HEX(3fbe,b2fb,1366,ea95,7d3e)   /* ~1.41421356237309504876378807303183293697657063603401 T1[128] hi,lo */
/*  64*/ ,LDOUBLE_HEX(3fff,8058,d7d2,d5e5,f6b1),LDOUBLE_HEX(bfbe,d654,ec13,ee23,6abc)   /* ~1.00271127505020248547613209710860360246442724019289 T2[128] hi,lo */
};
#define L 7     /* # bits in table index (i1,i2) */
#define K 128   /* table size = 2^L */
#define L_EXPAND    (((K/2) << L) + (K/2))
#define L_MASK      ((1 << L)-1)

#define constants_80  ((const LC_FLOAT_TYPE *)_constants_80)
# define _DE2IntCst   (constants_80[ 0])
# define _K2byLog2    (constants_80[ 1])
# define _Log2byK2Hi  (constants_80[ 2])
# define _Log2byK2Med (constants_80[ 3])
# define _p3_HI       (constants_80[ 4])

# define _p5          (constants_80[ 5])
# define _p4          (constants_80[ 6])
# define _p3_LO       (constants_80[ 7])

# define _Log2byK2Lo  (constants_80[ 8])
# define _p6          (constants_80[ 9])

__declspec(align(16)) static const unsigned short _constants_80[] = {
     LDOUBLE_HEX(403e,c000,0000,0000,0000)  /*  0 DE2IntCst      = 1.5*2^63 */
    ,LDOUBLE_HEX(400d,b8aa,3b29,5c17,f0bc)  /*  1 K2byLog2   = K^2/log(2) ~ 23637.11554992 */
    ,LDOUBLE_HEX(bff0,b172,17f7,d1cf,79ac)  /*  2 Log2byK2Hi ~ -4.23063464697232244524347213247665588919943502332898e-05 */
    ,LDOUBLE_HEX(3fae,d871,319f,f000,0000)  /*  3 Log2byK2Med ~  6.99362349047620977408942987412840832949912800870764e-25 */
    ,LDOUBLE_HEX(3ffc,aaaa,aaaa,aaaa,aaab)  /*  4 p3_HI  ~  1.66666666666666666671184175718689601808364386670291e-01 */

    ,LDOUBLE_HEX(3ff8,8888,8888,91d5,7943)  /*  5 p5     ~  8.33333333346550207291141649720844775117711833445355e-03 */
    ,LDOUBLE_HEX(3ffa,aaaa,aaaa,aaaa,aaab)  /*  6 p4     ~  4.16666666666666666677960439296724004520910966675728e-02 */
    ,LDOUBLE_HEX(bfbb,ab70,78a8,ae65,09cf)  /*  7 p3_LO  ~ -4.53796157421198095955516441874944572386592208568553e-21 */
 
    ,LDOUBLE_HEX(3f84,d095,0bf0,cbcd,98d6)  /*  8 Log2byK2Lo ~  1.53242008494748302860112907643687641670519495734080e-37 */
    ,LDOUBLE_HEX(3ff5,b60b,60b6,1286,f983)  /*  9 p6     ~  1.38888888890158889343717346078962981970050805102800e-03 */
};

//============================================================================================

double exp_rn( double xd )
{
    UINT64 x_val,x_abs,sign_mask,range,x_exp;
    SINT64 m,n,i1,i2,xsign,yesno;
    X_FLOAT_TYPE xx, rx, r2x, tp0x, tp1x, tpxx, tt1x, tt2x, tptx, resx;
    L_FLOAT_TYPE DE2IntCst,K2byLog2,Log2byK2Hi,Log2byK2Med,Log2byK2Lo;
    X_FLOAT_TYPE p3;
    L_FLOAT_TYPE p2,p4,p5,p6;
    L_FLOAT_TYPE tmp,w,mr,tmp1,tmp5,tmp6,tmp7,res,ef,rx2,x4,inf_d,minnorm_d,maxnorm_d;
    L_FLOAT_TYPE sc; // Don't know if it wouldn't be better as a double
    double resd;
    double volatile exception;

    /* load constants */
    DE2IntCst      = _DE2IntCst;
    //DE2IntCst = _Asm_setf(_FR_D, 0x43e8000000000000); /* 1 cycle slower */
    K2byLog2   = _K2byLog2;
    Log2byK2Hi = _Log2byK2Hi;
    Log2byK2Med = _Log2byK2Med;
    p3.hi  = _p3_HI;
    p2 = _Asm_setf(_FR_EXP, 0xfffe); /* 0.5 */
    p5     = _p5;
    p4     = _p4;
    p3.lo  = _p3_LO;


    /* get and classify input value; obtain range value */

    /* 40862e42fefa39ef (+) overflow range
     * 40874910d52d3051 (-) underflow range
     * 000167522bd709be (^)
     */

    x_val = _Asm_getf( _FR_D, xd );
    sign_mask = ((SINT64)x_val >> 63);
    x_exp = _Asm_extr_u( x_val, 52, 11 );
    x_abs = (x_val & ULL(7fffffffffffffff));
    range = ULL(40862e42fefa39ef) ^ (sign_mask & ULL(000167522bd709be));

    /* filter out special cases */

    if (__builtin_expect( x_exp < 0x3eb, 0 )) { /* |x| < 2^(-20) */
        if (__builtin_expect( x_exp < 0x3c8, 0 )) { /* |x| < 2^(-55) */
            return _Asm_fma( _PC_D, xd, 1, 1, _SF0 );
        }
        /* 2^(-55) <= |x| < 2^(-20) */
	//	printf("\nToto %1.30e\n", xd);
         tmp = p5*xd + p4;
        __X_SQR_L( r2x, xd );
        tmp = p6*r2x.hi + tmp;
        tp1x.hi = p2;
        tp1x.lo = tmp*r2x.hi;
        __X_FMA_GREATER_LXX( tp0x, xd, p3, tp1x );   /* tp0 = xd * p3   +  tp1 */
        __X_FMA_GREATER_XXL( tpxx, r2x, tp0x, xd );  /* tpx = r2 * tp0  +  xd*/
        /* res = (tpxx.lo + tpxx.hi) + 1; */
        res = _Asm_fma( _PC_D, tpxx.hi, 1, 1, _SF1 );
        tmp = (1 - res);
        tmp = (tmp + tpxx.hi);
        tmp = (tmp + tpxx.lo);
        res = _Asm_fma( _PC_D, res, 1, tmp, _SF0 );
        return res;
    }

    if (__builtin_expect( x_abs > range, 0 )) { /* nan, inf, large finite */
        /* NaNs */
        if (x_abs > ULL(7ff0000000000000)) {
            if (x_abs <= ULL(7ff7ffffffffffff)) {
                inf_d = _Asm_setf( _FR_D, ULL(7ff0000000000000) );
                exception = _Asm_fma( _PC_D, inf_d, 0, 0, _SF0 );
            }
            return _Asm_fma( _PC_D, xd, 1, 0, _SF0 );
        }
        /* INFs */
        if (x_abs == ULL(7ff0000000000000)) {
            if (sign_mask)  /* -inf */
                return 0;
            else            /* +inf */
                return xd;
        }
        /* large finite */
        if (sign_mask) { /* underflow */
            minnorm_d = _Asm_setf( _FR_D, ULL(0010000000000000) );
            return _Asm_fma( _PC_D, minnorm_d, minnorm_d, 0, _SF0 );
        } else {        /* overflow */
            maxnorm_d = _Asm_setf( _FR_D, ULL(7fefffffffffffff) );
            return _Asm_fma( _PC_D, maxnorm_d, maxnorm_d, 0, _SF0 );
        }
    }



    /* main path: 2^(-55) <= |x| < range */

    /* reduction: x = m*log(2)/(K^2) + r; m = n*K^2 + i1*K + i2; |r|=[0..log(2)/(2*K^2)] */

    w = xd * K2byLog2 + DE2IntCst;   /* Double2Int on a doubledouble */
    mr = (w - DE2IntCst);
    m = _Asm_getf( _FR_SIG, w ) + L_EXPAND; /* add L/2 to i1,i2 (make it unsigned table offset) */
 
    rx.hi = mr * Log2byK2Hi + xd;


    i1 = (m >> L) & L_MASK;
    i2 = m & L_MASK;

    // was:   n = _Asm_extr( m, 2*L, 16 ); The compiler seems to get it
    n = m >> 2*L;
    sc = _Asm_setf( _FR_EXP, 0xffff + n ); /* build 2^n */
    //sc = _Asm_setf( _FR_D, ULL(3ff0000000000000) + (n<<52) ); /* build 2^n */


    /* First step */
    /* Compute e^xred as 1+p(xred) 
       Then compute T*e^xred as T+T*p(xred) */ 

    /* table lookup and reconstruction: T = T1[i1]*T2[i2]*2^n */
    /* We read both vlues but it isn't slower  */

    __X_CONVERT_XC2X( tt1x, TX[i1*2+0] );   /* load T1[i1] */
    __X_CONVERT_XC2X( tt2x, TX[i2*2+1] );   /* load T2[i2] */

    tptx.hi = tt2x.hi * tt1x.hi ;

      
    /* degree 3 polynomial evaluation:  */
    rx2=rx.hi*rx.hi;
    tmp = p2 + p3.hi*rx.hi;  

    // We're not adding a correcting term to rx.hi, we're multiplying the final exp by (1+f)
    rx.hi = rx.hi + mr * Log2byK2Med;

    //ef = 1.0 + mr * Log2byK2Med;
    tmp = rx.hi +rx2*tmp; /* tmp = x+x^2/2 + ...*/

    /* Final reconstruction and rounding test */
    
    res = tmp*tptx.hi + tptx.hi;
    resd = res*sc;    /* do the scaling now, saves a few cycles */

#if 0 /* to time the first step alone */
    return resd;
#endif

#if 0 // test a la Ziv
    /* This test is 3 cycles slower */
    L_FLOAT_TYPE reshi,reslo,test,rn_constant;
    reshi = _Asm_fma( _PC_D, tmp, tptx.hi, tptx.hi, _SF0 );      /*  tmp*tptx.hi + tptx.hi, round to double */
    rn_constant = 1.005;
    reslo = res - reshi;  
    test=_Asm_fma( _PC_D, reslo, rn_constant, reshi, _SF0 );
  
    if(__builtin_expect( reshi==test, 1+1==2) ) {
      return resd;   /* do the scaling now, saves a few cycles */
    }
#else // else test a la Ziv

    i1 = _Asm_getf( _FR_SIG, res);
    m =  i1 & (0xff<<3);
    if(__builtin_expect((m!=(0x7f<<3) && m!=(0x80<<3)), 1+1==2)) {
      return resd;
    }      
#endif // end test a la Ziv
    else {
      /***************************************************************************************************/
      /**************************************** second step **********************************************/
      /***************************************************************************************************/

#if EVAL_PERF
        crlibm_second_step_taken++;
#endif

    Log2byK2Lo = _Log2byK2Lo;
    p6     = _p6;

    // Cut from the first step
    tmp1 = mr * Log2byK2Med;
    tmp7 = mr * Log2byK2Hi + xd;

    tmp6 = _Asm_famax( tmp7, tmp1, _SF1 );
    tmp5 = _Asm_famin( tmp1, tmp7, _SF1 );
    rx.lo = (tmp6 - rx.hi);
    rx.lo = (tmp5 + rx.lo);
    rx.lo = (mr * Log2byK2Lo + rx.lo);

      
      /* table lookup and reconstruction: T = T1[i1]*T2[i2]*2^n */
      
      __X_MUL_XX( tptx, tt1x, tt2x ); /* T = T1[i1]*T2[i2] */
      
      tptx.hi *= sc;  /* scale T' = T*2^n */
      tptx.lo *= sc;
      
      /* polynomial evaluation: px = P(r)*r */
      
      tmp = p5*rx.hi + p4;
      
      __X_SQR_X( r2x, rx );
      
      tmp = p6*r2x.hi + tmp;
      tp1x.hi = p2;
      tp1x.lo = tmp*r2x.hi;
      
      __X_FMA_GREATER_XXX( tp0x, rx, p3, tp1x );
      __X_FMA_GREATER_XXX( tpxx, r2x, tp0x, rx );
      
      /* final reconstruction: res = (P(r)*r + 1)*T*2^n */
      
      __X_FMA_GREATER_XXX( resx, tptx, tpxx, tptx );
      
      /* convert and return result */

      return (double)( resx.hi + resx.lo );
    }
}





/* Just to compile OK within crlibm */

double exp_ru(double x) {
  return exp_rn(x);
}
double exp_rd(double x) {
  return exp_rn(x);
}
double exp_rz(double x) {
  return exp_rn(x);
}
