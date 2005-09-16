/* 
 *this function computes a correctly rounded atan using double-extended arithmetic, FMAs and other dirty tricks 
 *
 * Author : Nicolas Gast, Florent de Dinechin
 * nicolas.gast@ens.fr
 *

WARNING : This code is dirty and experimental, and remains here for
history. A cleaner, portable version using double-extended arithmetic will be available some day as atan-de.c
For this reason there is only atan_rn so it fails the "make check" for all the other rounding modes


   To test within crlibm: (tested with Intel icc compiler version 8.1)
icc -Qoption,cpp,--extended_float_types -IPF_fp_speculationsafe -c atan-itanium.c; mv atan-itanium.o atan_fast.o; make




This file is completely self-contained so that we can change the crlibm infrastructure without bothering maintaining it.


 */

/* WARNING Due to some quantum effect not understood so far, 
   turning debugging on may change the result */
#define DEBUG 0



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


/* For debugging */
typedef union {
  int i[3];                 
  long double d;
} db_ext_number;


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

#define print_debug(msg, _z) {\
  db_ext_number dbg;\
  dbg.d=_z;\
  printf(msg);\
  printf(" %08x %08x %08x \n", (dbg.i[2]<<16)>>16, dbg.i[1], dbg.i[0]);\
}


#define Add12_ext(s, r, a, b)            \
    { L_FLOAT_TYPE _z, _a, _b, _s;       \
      _a= (a); _b=(b);                   \
      s = (_a + _b);                     \
      _z= ( a - s );                     \
      r = (_b + _z); }            


#define Add22_ext(zh,zl,xh,xl,yh,yl) \
do {\
L_FLOAT_TYPE r,s;\
r = (xh)+(yh);\
s = (xh)-r;\
s+= (yh);\
s+= (yl);\
s+= (xl);\
zh = r+s;\
zl = r - (zh);\
zl+=  s;\
} while(0)



#define Mul12_ext(_rh,_rl,_u,_v)                    \
{                                                     \
  _rh = _u*_v;                                      \
  _rl = _Asm_fms( 3/*_PC_NONE*/, _u, _v, _rh, 1/*_SF1*/ );\
}
#define Mul22_ext(zh,zl, xh,xl, yh,yl)              \
{                                                     \
L_FLOAT_TYPE ph, pl;                                   \
  ph = (xh)*(yh);                                         \
  pl = _Asm_fms( 3/*_PC_NONE*/, xh, yh, ph, 1/*_SF1*/ );;  \
  pl = (xh)*(yl) + pl;                                    \
  pl = (xl)*(yh) + pl;                                    \
  zh = ph+pl;					      \
  zl = ph - zh;                                 \
  zl += pl;                                         \
}

#define  Div22_ext(zh,zl,xh,xl,yh,yl) \
  {           \
L_FLOAT_TYPE _ch,_cl,_uh,_ul;                        \
  _ch=(xh)/(yh);  \
  Mul12_ext(_uh,_ul,_ch,(yh));    \
  _cl=(xh)-_uh;                                   \
  _cl -= _ul;                                       \
  _cl += (xl);                                      \
  _cl -= _ch*(yl);                               \
  _cl /= (yh);                                   \
  zh = _ch + _cl;                                  \
  zl=(_ch-(zh)); zl += _cl;                  \
}





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



double dde_atan_rn(double x) {
  return 0;
}

double atan_rd(double x) {
  return 0;
}

double atan_ru(double x) {
  return 0;
}

double atan_rz(double x) {
  return 0;
}


static const double  HALFPI = 1.57079632679489655799898173427209258079528808593750e+00;
#define MIN_REDUCTION_NEEDED ULL(3F89FDF8BCCE533D)
#define A 0
#define B 1
#define ATAN_BHI 0
#define ATAN_BLO 1
#define epsilon 2.04221581890623872536809598138553304900554884091659e-19
#define epsilon_no_red 1.56771350764719825686165002299335165493769973908433e-19
#define TWO_M_64 5.42101086242752217003726400434970855712890625000000e-20
#define TWO_10 1.02400000000000000000000000000000000000000000000000e+03

__declspec(align(16)) 

static const struct{long long int a; double b;} ab_table[62] = {
 { /*a[0] ~= 1.26914436930661800408e-02   */   ULL(3F89FDF8BCCE533D),
   /*b[0] = */   2.53869765124364009378776785297304741106927394866943e-02},
 { /*a[1] ~= 3.80906929270782388369e-02   */   ULL(3FA3809F90CEBC31),
   /*b[1] = */   5.08066978456951506837313559117319528013467788696289e-02},
 { /*a[2] ~= 6.35391122156262234502e-02   */   ULL(3FB0441968FBA526),
   /*b[2] = */   7.62920780032335793530151590857713017612695693969727e-02},
 { /*a[3] ~= 8.90697640843219481662e-02   */   ULL(3FB6CD46ABCDFA25),
   /*b[3] = */   1.01876371166982934712841313285025535151362419128418e-01},
 { /*a[4] ~= 1.14716138034642060814e-01   */   ULL(3FBD5E096D2EA546),
   /*b[4] = */   1.27593346472767293908745500630175229161977767944336e-01},
 { /*a[5] ~= 1.40512327929006382604e-01   */   ULL(3FC1FC4ED691E891),
   /*b[5] = */   1.53477468508642272970732278736250009387731552124023e-01},
 { /*a[6] ~= 1.66493216120905490981e-01   */   ULL(3FC54FA6531F610B),
   /*b[6] = */   1.79564085612852891715718328669026959687471389770508e-01},
 { /*a[7] ~= 1.92694666476959805056e-01   */   ULL(3FC8AA380550EAF1),
   /*b[7] = */   2.05889628199359991933548030829115305095911026000977e-01},
 { /*a[8] ~= 2.19153728611415840590e-01   */   ULL(3FCC0D3AB8975BD9),
   /*b[8] = */   2.32491819536184141092860500066308304667472839355469e-01},
 { /*a[9] ~= 2.45908855876056406352e-01   */   ULL(3FCF79F0FEE46885),
   /*b[9] = */   2.59409901651160901270287695297156460583209991455078e-01},
 { /*a[10] ~= 2.73000139926648314534e-01   */   ULL(3FD178D5943274CA),
   /*b[10] = */   2.86684879348826082701151563014718703925609588623047e-01},
 { /*a[11] ~= 3.00469565029600954026e-01   */   ULL(3FD33AE4B2CFB5F7),
   /*b[11] = */   3.14359785700871030567071784389554522931575775146484e-01},
 { /*a[12] ~= 3.28361285690481766972e-01   */   ULL(3FD503DF0DD40A5B),
   /*b[12] = */   3.42479972833279300292730340515845455229282379150391e-01},
 { /*a[13] ~= 3.56721931693259067415e-01   */   ULL(3FD6D4883998DD14),
   /*b[13] = */   3.71093432391343347465095803272561170160770416259766e-01},
 { /*a[14] ~= 3.85600945252912822931e-01   */   ULL(3FD8ADAF964ABFA5),
   /*b[14] = */   4.00251150738601846335029676993144676089286804199219e-01},
 { /*a[15] ~= 4.15050955725992373816e-01   */   ULL(3FDA9031E241114E),
   /*b[15] = */   4.30007504761513281721363455289974808692932128906250e-01},
 { /*a[16] ~= 4.45128198220858643198e-01   */   ULL(3FDC7CFAFB78B41D),
   /*b[16] = */   4.60420705138676944478959285333985462784767150878906e-01},
 { /*a[17] ~= 4.75892983535655022698e-01   */   ULL(3FDE7507D82B9DC6),
   /*b[17] = */   4.91553295129659728601723145402502268552780151367188e-01},
 { /*a[18] ~= 5.07410228170177493351e-01   */   ULL(3FE03CB45FF4B2AB),
   /*b[18] = */   5.23472714391912563591802154405741021037101745605469e-01},
 { /*a[19] ~= 5.39750054761637805872e-01   */   ULL(3FE145A1E826E4EA),
   /*b[19] = */   5.56251939105489867642972967587411403656005859375000e-01},
 { /*a[20] ~= 5.72988475252136329570e-01   */   ULL(3FE255EBED462BAC),
   /*b[20] = */   5.89970211851368997457711884635500609874725341796875e-01},
 { /*a[21] ~= 6.07208171494496387417e-01   */   ULL(3FE36E3FD4CDD9AC),
   /*b[21] = */   6.24713877348479162954220100800739601254463195800781e-01},
 { /*a[22] ~= 6.42499390954343656748e-01   */   ULL(3FE48F5AE1FB2991),
   /*b[22] = */   6.60577343433393693317157158162444829940795898437500e-01},
 { /*a[23] ~= 6.78960978813340497734e-01   */   ULL(3FE5BA0C5FE86E27),
   /*b[23] = */   6.97664190728041089251121320558013394474983215332031e-01},
 { /*a[24] ~= 7.16701572306941533027e-01   */   ULL(3FE6EF3822C19A5D),
   /*b[24] = */   7.36088459496464064812926153535954654216766357421875e-01},
 { /*a[25] ~= 7.55840988781748695010e-01   */   ULL(3FE82FD970F967BD),
   /*b[25] = */   7.75976148518263131315109148999908939003944396972656e-01},
 { /*a[26] ~= 7.96511846049556065643e-01   */   ULL(3FE97D0669351A0D),
   /*b[26] = */   8.17466968767843527032823658373672515153884887695312e-01},
 { /*a[27] ~= 8.38861462565995493716e-01   */   ULL(3FEAD7F3FE730FCD),
   /*b[27] = */   8.60716404767067566616844942473107948899269104003906e-01},
 { /*a[28] ~= 8.83054096327761096527e-01   */   ULL(3FEC41FAAA0A733E),
   /*b[28] = */   9.05898149317818313086547732382314279675483703613281e-01},
 { /*a[29] ~= 9.29273595909162105525e-01   */   ULL(3FEDBC9BFAEEEADF),
   /*b[29] = */   9.53206993785724487899813084368361160159111022949219e-01},
 { /*a[30] ~= 9.77726555752981254442e-01   */   ULL(3FEF498933AC790A),
   /*b[30] = */   1.00286227737052557884567249857354909181594848632812e+00},
 { /*a[31] ~= 1.02864609206350806308e+00   */   ULL(3FF075559AC922B4),
   /*b[31] = */   1.05511202646791502068879253783961758017539978027344e+00},
 { /*a[32] ~= 1.08229638730567912228e+00   */   ULL(3FF151160440E8D3),
   /*b[32] = */   1.11023795151925819268967643438372761011123657226562e+00},
 { /*a[33] ~= 1.13897819300824741364e+00   */   ULL(3FF23941329D3DD8),
   /*b[33] = */   1.16856151675095110142876819736557081341743469238281e+00},
 { /*a[34] ~= 1.19903553596580987055e+00   */   ULL(3FF32F3FE2DB7094),
   /*b[34] = */   1.23045136228081597451478046423289924860000610351562e+00},
 { /*a[35] ~= 1.26286394722716532198e+00   */   ULL(3FF434B0D38A35D7),
   /*b[35] = */   1.29633244442242001603915468876948580145835876464844e+00},
 { /*a[36] ~= 1.33092063388866265448e+00   */   ULL(3FF54B736F41F96D),
   /*b[36] = */   1.36669737760087572908673791971523314714431762695312e+00},
 { /*a[37] ~= 1.40373715148086145849e+00   */   ULL(3FF675B5165CA5E1),
   /*b[37] = */   1.44212062317890032936418265308020636439323425292969e+00},
 { /*a[38] ~= 1.48193532552453321547e+00   */   ULL(3FF7B601D0DEA3C6),
   /*b[38] = */   1.52327639603630871079076314345002174377441406250000e+00},
 { /*a[39] ~= 1.56624743831976717041e+00   */   ULL(3FF90F5979506F51),
   /*b[39] = */   1.61096147803441858137318831722950562834739685058594e+00},
 { /*a[40] ~= 1.65754207708184630948e+00   */   ULL(3FFA854AD74CF791),
   /*b[40] = */   1.70612458293084490179580825497396290302276611328125e+00},
 { /*a[41] ~= 1.75685758736121174681e+00   */   ULL(3FFC1C16B3972246),
   /*b[41] = */   1.80990457885083300126893846027087420225143432617188e+00},
 { /*a[42] ~= 1.86544587781964938190e+00   */   ULL(3FFDD8DDC6DB1831),
   /*b[42] = */   1.92368085119253517945026032975874841213226318359375e+00},
 { /*a[43] ~= 1.98483051718814034750e+00   */   ULL(3FFFC1DDA4F6D032),
   /*b[43] = */   2.04914055707593512067887786542996764183044433593750e+00},
 { /*a[44] ~= 2.11688487740990979279e+00   */   ULL(4000EF6156AEFAF2),
   /*b[44] = */   2.18836977316091063627823132264893501996994018554688e+00},
 { /*a[45] ~= 2.26393888595347935033e+00   */   ULL(40021C8BFD9A80C1),
   /*b[45] = */   2.34397906437763481335423421114683151245117187500000e+00},
 { /*a[46] ~= 2.42892740222016626128e+00   */   ULL(40036E717D67269C),
   /*b[46] = */   2.51927965826279764982587039412464946508407592773438e+00},
 { /*a[47] ~= 2.61560046981161264128e+00   */   ULL(4004ECBFF069F1E4),
   /*b[47] = */   2.71853573297491069027387311507482081651687622070312e+00},
 { /*a[48] ~= 2.82882779840766906527e+00   */   ULL(4006A170780169B7),
   /*b[48] = */   2.94733416149008720097413061012048274278640747070312e+00},
 { /*a[49] ~= 3.07505072362971616974e+00   */   ULL(400899B4319C3F02),
   /*b[49] = */   3.21314087722892072207514502224512398242950439453125e+00},
 { /*a[50] ~= 3.36297230191158715455e+00   */   ULL(400AE75E05B0834A),
   /*b[50] = */   3.52616384863255349912947167467791587114334106445312e+00},
 { /*a[51] ~= 3.70464601821196143254e+00   */   ULL(400DA31D739BD0E3),
   /*b[51] = */   3.90073973345466518125590482668485492467880249023438e+00},
 { /*a[52] ~= 4.11726034471856573100e+00   */   ULL(401078131886BC57),
   /*b[52] = */   4.35765668014056828383218089584261178970336914062500e+00},
 { /*a[53] ~= 4.62619989820137847648e+00   */   ULL(4012813A8BCE2241),
   /*b[53] = */   4.92824409985376998832862227573059499263763427734375e+00},
 { /*a[54] ~= 5.27059285056349616385e+00   */   ULL(401515164ACECE78),
   /*b[54] = */   5.66202526987798027136022938066162168979644775390625e+00},
 { /*a[55] ~= 6.11406930017863578891e+00   */   ULL(401874CE9526FAB9),
   /*b[55] = */   6.64216890962962569489036468439735472202301025390625e+00},
 { /*a[56] ~= 7.26750136287798241547e+00   */   ULL(401D11EBE094C913),
   /*b[56] = */   8.01990986231011859786121931392699480056762695312500e+00},
 { /*a[57] ~= 8.94284159107796650204e+00   */   ULL(4021E2BC220DFA19),
   /*b[57] = */   1.01020964280653942068965989165008068084716796875000e+01},
 { /*a[58] ~= 1.16023240149353498339e+01   */   ULL(40273463D0337C49),
   /*b[58] = */   1.36206610885392880305744256475009024143218994140625e+01},
 { /*a[59] ~= 1.64826377753716631495e+01   */   ULL(40307B8E26350916),
   /*b[59] = */   2.08587363260064613257327437167987227439880371093750e+01},
 { /*a[60] ~= 2.83859754493341325216e+01   */   ULL(403C62CF497BF2F2),
   /*b[60] = */   4.43908820444562195461912779137492179870605468750000e+01},
 { /*a[61] ~= 1.01699461607316896213e+02   */   ULL(40596CC3FA9E0EF4),
   /*b[61] = */   8.27932424540746438879068591631948947906494140625000e+01}
};


#define atanb_table ((const XC_FLOAT_TYPE *)_atanb_table)
__declspec(align(16)) static const unsigned short _atanb_table[] = {
 /*atan_b[0] ~= 2.5381524664e-02*/
  LDOUBLE_HEX(3FF9, CFEC, EA4B, 4FCB, 5DFD), 
   LDOUBLE_HEX(BFB7, CBBA, 8342, F523, 8BE7), 
 /*atan_b[1] ~= 5.0763049304e-02*/
  LDOUBLE_HEX(3FFA, CFEC, EA49, B131, 647C), 
   LDOUBLE_HEX(3FB6, D38B, A5E1, 4DEF, A6BD), 
 /*atan_b[2] ~= 7.6144573921e-02*/
  LDOUBLE_HEX(3FFB, 9BF1, AFB6, 0F03, 5D53), 
   LDOUBLE_HEX(3FB8, EF7C, 871F, DC70, BCA9), 
 /*atan_b[3] ~= 1.0152609851e-01*/
  LDOUBLE_HEX(3FFB, CFEC, EA46, 78CC, AECA), 
   LDOUBLE_HEX(BFB7, DCB7, 3BED, 3BD7, 633C), 
 /*atan_b[4] ~= 1.2690762308e-01*/
  LDOUBLE_HEX(3FFC, 81F4, 126B, 0C0A, B24C), 
   LDOUBLE_HEX(3FB8, 9C93, 50C6, 8748, 202B), 
 /*atan_b[5] ~= 1.5228914763e-01*/
  LDOUBLE_HEX(3FFC, 9BF1, AFB2, 77C1, F1F3), 
   LDOUBLE_HEX(BFBB, 9D89, 6B54, 2B43, C3D3), 
 /*atan_b[6] ~= 1.7767067216e-01*/
  LDOUBLE_HEX(3FFC, B5EF, 4CF9, 8121, 27D9), 
   LDOUBLE_HEX(BFBB, D8AB, 134C, C337, 1424), 
 /*atan_b[7] ~= 2.0305219666e-01*/
  LDOUBLE_HEX(3FFC, CFEC, EA40, 29FE, 3D0C), 
   LDOUBLE_HEX(BFBA, 964C, 23A5, 78A9, 286C), 
 /*atan_b[8] ~= 2.2843372114e-01*/
  LDOUBLE_HEX(3FFC, E9EA, 8786, 746E, CBDE), 
   LDOUBLE_HEX(3FBB, 95CE, 8C74, D4B3, 3D3D), 
 /*atan_b[9] ~= 2.5381524560e-01*/
  LDOUBLE_HEX(3FFD, 81F4, 1266, 3163, 58ED), 
   LDOUBLE_HEX(3FBB, B292, B8DC, 903F, C86D), 
 /*atan_b[10] ~= 2.7919677004e-01*/
  LDOUBLE_HEX(3FFD, 8EF2, E108, FBCB, 4839), 
   LDOUBLE_HEX(BFBC, C5E3, D3F8, 42F0, A001), 
 /*atan_b[11] ~= 3.0457829447e-01*/
  LDOUBLE_HEX(3FFD, 9BF1, AFAB, 9AD5, 051A), 
   LDOUBLE_HEX(3FBC, BE9C, AF21, 45D0, CBC5), 
 /*atan_b[12] ~= 3.2995981887e-01*/
  LDOUBLE_HEX(3FFD, A8F0, 7E4E, 1002, FE3F), 
   LDOUBLE_HEX(3FB9, ACDF, 4585, 84D5, 7EE8), 
 /*atan_b[13] ~= 3.5534134325e-01*/
  LDOUBLE_HEX(3FFD, B5EF, 4CF0, 5CF3, 3B2F), 
   LDOUBLE_HEX(BFB9, DAF1, E542, E461, 5C3F), 
 /*atan_b[14] ~= 3.8072286762e-01*/
  LDOUBLE_HEX(3FFD, C2EE, 1B92, 835E, 5241), 
   LDOUBLE_HEX(3FBC, F450, E872, E8D5, 5B89), 
 /*atan_b[15] ~= 4.0610439197e-01*/
  LDOUBLE_HEX(3FFD, CFEC, EA34, 8516, 3E60), 
   LDOUBLE_HEX(BFBC, 91DD, F6E6, 0680, E8AD), 
 /*atan_b[16] ~= 4.3148591630e-01*/
  LDOUBLE_HEX(3FFD, DCEB, B8D6, 6405, 31AA), 
   LDOUBLE_HEX(BFBC, 8502, E09D, 5663, 1B39), 
 /*atan_b[17] ~= 4.5686744062e-01*/
  LDOUBLE_HEX(3FFD, E9EA, 8778, 222C, 48BB), 
   LDOUBLE_HEX(BFBB, F51E, C2F3, 5A3E, F53D), 
 /*atan_b[18] ~= 4.8224896492e-01*/
  LDOUBLE_HEX(3FFD, F6E9, 5619, C1A2, 5014), 
   LDOUBLE_HEX(BFBB, E1E1, FABB, 35B7, 64D8), 
 /*atan_b[19] ~= 5.0763048922e-01*/
  LDOUBLE_HEX(3FFE, 81F4, 125D, A249, 1B96), 
   LDOUBLE_HEX(BFBB, FEB6, 20F5, A80E, ABD8), 
 /*atan_b[20] ~= 5.3301201350e-01*/
  LDOUBLE_HEX(3FFE, 8873, 79AE, 569C, E82C), 
   LDOUBLE_HEX(BFBD, 9333, CB85, 3253, A31F), 
 /*atan_b[21] ~= 5.5839353776e-01*/
  LDOUBLE_HEX(3FFE, 8EF2, E0FE, FEF4, 22DF), 
   LDOUBLE_HEX(3FBD, FBF4, E487, 2960, 19F2), 
 /*atan_b[22] ~= 5.8377506202e-01*/
  LDOUBLE_HEX(3FFE, 9572, 484F, 9C7E, 4569), 
   LDOUBLE_HEX(BFBD, ED41, 6021, 317B, 1548), 
 /*atan_b[23] ~= 6.0915658627e-01*/
  LDOUBLE_HEX(3FFE, 9BF1, AFA0, 3071, E801), 
   LDOUBLE_HEX(3FBD, C46B, 95C4, B736, D8A5), 
 /*atan_b[24] ~= 6.3453811052e-01*/
  LDOUBLE_HEX(3FFE, A271, 16F0, BC0B, F541), 
   LDOUBLE_HEX(3FBD, E479, 64B6, 873E, E8BE), 
 /*atan_b[25] ~= 6.5991963475e-01*/
  LDOUBLE_HEX(3FFE, A8F0, 7E41, 408E, DDC6), 
   LDOUBLE_HEX(3FBD, C200, D1A3, 7D02, 9DAA), 
 /*atan_b[26] ~= 6.8530115898e-01*/
  LDOUBLE_HEX(3FFE, AF6F, E591, BF41, BD98), 
   LDOUBLE_HEX(3FBC, AB83, 86B7, DBD3, 49B9), 
 /*atan_b[27] ~= 7.1068268321e-01*/
  LDOUBLE_HEX(3FFE, B5EF, 4CE2, 396F, 887A), 
   LDOUBLE_HEX(3FB9, 93C0, 6F69, 2472, DD13), 
 /*atan_b[28] ~= 7.3606420743e-01*/
  LDOUBLE_HEX(3FFE, BC6E, B432, B066, 2617), 
   LDOUBLE_HEX(BFBD, C5F2, 72DA, A216, 8845), 
 /*atan_b[29] ~= 7.6144573166e-01*/
  LDOUBLE_HEX(3FFE, C2EE, 1B83, 2575, A17C), 
   LDOUBLE_HEX(3FBA, FC52, 25AC, D135, 67B0), 
 /*atan_b[30] ~= 7.8682725588e-01*/
  LDOUBLE_HEX(3FFE, C96D, 82D3, 99EF, 4753), 
   LDOUBLE_HEX(3FBC, E6CB, 9CE5, F7DC, 32EF), 
 /*atan_b[31] ~= 8.1220878010e-01*/
  LDOUBLE_HEX(3FFE, CFEC, EA24, 0F24, C5A3), 
   LDOUBLE_HEX(BFBB, 9F94, 64A4, 0D49, 77DA), 
 /*atan_b[32] ~= 8.3759030433e-01*/
  LDOUBLE_HEX(3FFE, D66C, 5174, 8667, 5086), 
   LDOUBLE_HEX(BFBC, E480, 36A7, 98A0, E416), 
 /*atan_b[33] ~= 8.6297182855e-01*/
  LDOUBLE_HEX(3FFE, DCEB, B8C5, 0106, C115), 
   LDOUBLE_HEX(BFBB, AE5E, 111C, 0925, 5FC1), 
 /*atan_b[34] ~= 8.8835335278e-01*/
  LDOUBLE_HEX(3FFE, E36B, 2015, 8050, B874), 
   LDOUBLE_HEX(BFBC, 8DD3, E1A9, 67EE, B236), 
 /*atan_b[35] ~= 9.1373487702e-01*/
  LDOUBLE_HEX(3FFE, E9EA, 8766, 058F, C400), 
   LDOUBLE_HEX(BFBD, 994E, 5D94, 7944, 5BF2), 
 /*atan_b[36] ~= 9.3911640126e-01*/
  LDOUBLE_HEX(3FFE, F069, EEB6, 920A, 8756), 
   LDOUBLE_HEX(BFBD, F0FC, 830B, 5639, 9FED), 
 /*atan_b[37] ~= 9.6449792552e-01*/
  LDOUBLE_HEX(3FFE, F6E9, 5607, 2702, D403), 
   LDOUBLE_HEX(BFBD, B0EF, D9DB, FF7A, BBF3), 
 /*atan_b[38] ~= 9.8987944978e-01*/
  LDOUBLE_HEX(3FFE, FD68, BD57, C5B4, F372), 
   LDOUBLE_HEX(BFBD, 9706, 5831, 4248, 656E), 
 /*atan_b[39] ~= 1.0152609740e+00*/
  LDOUBLE_HEX(3FFF, 81F4, 1254, 37AB, 59C4), 
   LDOUBLE_HEX(3FBE, C83B, C3BE, 8160, FE56), 
 /*atan_b[40] ~= 1.0406424983e+00*/
  LDOUBLE_HEX(3FFF, 8533, C5FC, 928B, 5DCD), 
   LDOUBLE_HEX(3FBE, C025, 7DA6, 5435, CDA0), 
 /*atan_b[41] ~= 1.0660240226e+00*/
  LDOUBLE_HEX(3FFF, 8873, 79A4, F40D, D390), 
   LDOUBLE_HEX(BFBE, BB70, CBE8, FB3B, AA03), 
 /*atan_b[42] ~= 1.0914055469e+00*/
  LDOUBLE_HEX(3FFF, 8BB3, 2D4D, 5CC1, ADB6), 
   LDOUBLE_HEX(3FBE, 8161, 18FB, A932, 136B), 
 /*atan_b[43] ~= 1.1167870712e+00*/
  LDOUBLE_HEX(3FFF, 8EF2, E0F5, CD31, 1F80), 
   LDOUBLE_HEX(BFBC, BD96, 57B0, 5730, 7576), 
 /*atan_b[44] ~= 1.1421685956e+00*/
  LDOUBLE_HEX(3FFF, 9232, 949E, 45E1, 3E02), 
   LDOUBLE_HEX(BFBD, CDB1, 87A1, 5D56, 06EC), 
 /*atan_b[45] ~= 1.1675501199e+00*/
  LDOUBLE_HEX(3FFF, 9572, 4846, C751, B4C7), 
   LDOUBLE_HEX(BFBD, A1AB, 140B, 2B49, DF68), 
 /*atan_b[46] ~= 1.1929316443e+00*/
  LDOUBLE_HEX(3FFF, 98B1, FBEF, 51FC, 635A), 
   LDOUBLE_HEX(3FBE, CA64, 3ADC, 86D5, FB02), 
 /*atan_b[47] ~= 1.2183131687e+00*/
  LDOUBLE_HEX(3FFF, 9BF1, AF97, E655, 1527), 
   LDOUBLE_HEX(3FBE, CA1D, 3262, C2F9, D84C), 
 /*atan_b[48] ~= 1.2436946931e+00*/
  LDOUBLE_HEX(3FFF, 9F31, 6340, 84C9, 33A7), 
   LDOUBLE_HEX(3FBD, AF23, 2B16, BE75, 8B87), 
 /*atan_b[49] ~= 1.2690762175e+00*/
  LDOUBLE_HEX(3FFF, A271, 16E9, 2DBF, 7CA7), 
   LDOUBLE_HEX(3FBE, FDDA, 7599, 4DA2, 0F86), 
 /*atan_b[50] ~= 1.2944577420e+00*/
  LDOUBLE_HEX(3FFF, A5B0, CA91, E197, C307), 
   LDOUBLE_HEX(BFBC, D265, 9307, D567, 08BE), 
 /*atan_b[51] ~= 1.3198392664e+00*/
  LDOUBLE_HEX(3FFF, A8F0, 7E3A, A0AA, A7E2), 
   LDOUBLE_HEX(3FBE, BE3C, 4D06, 7D11, 0641), 
 /*atan_b[52] ~= 1.3452207909e+00*/
  LDOUBLE_HEX(3FFF, AC30, 31E3, 6B49, 6713), 
   LDOUBLE_HEX(BFBE, B9DD, 9D13, C459, 6F6C), 
 /*atan_b[53] ~= 1.3706023154e+00*/
  LDOUBLE_HEX(3FFF, AF6F, E58C, 41BD, 9EA8), 
   LDOUBLE_HEX(BFBD, 802F, 2153, DC49, 3698), 
 /*atan_b[54] ~= 1.3959838399e+00*/
  LDOUBLE_HEX(3FFF, B2AF, 9935, 2449, 1D44), 
   LDOUBLE_HEX(3FBE, CAFC, 43E2, 3F23, 5075), 
 /*atan_b[55] ~= 1.4213653645e+00*/
  LDOUBLE_HEX(3FFF, B5EF, 4CDE, 1325, B93A), 
   LDOUBLE_HEX(BFBA, 9155, 4FBC, 9598, FA3D), 
 /*atan_b[56] ~= 1.4467468891e+00*/
  LDOUBLE_HEX(3FFF, B92F, 0087, 0E85, 296B), 
   LDOUBLE_HEX(3FBE, C76A, DB5B, 6055, 9EA6), 
 /*atan_b[57] ~= 1.4721284137e+00*/
  LDOUBLE_HEX(3FFF, BC6E, B430, 1690, E405), 
   LDOUBLE_HEX(3FBA, A6CB, 4564, 7FF8, 4121), 
 /*atan_b[58] ~= 1.4975099383e+00*/
  LDOUBLE_HEX(3FFF, BFAE, 67D9, 2B6A, 02AA), 
   LDOUBLE_HEX(BFBD, B0AE, B984, 420B, 761D), 
 /*atan_b[59] ~= 1.5228914629e+00*/
  LDOUBLE_HEX(3FFF, C2EE, 1B82, 4D29, 2EBE), 
   LDOUBLE_HEX(BFBE, 9CBD, 26E8, 9FF8, E917), 
 /*atan_b[60] ~= 1.5482729876e+00*/
  LDOUBLE_HEX(3FFF, C62D, CF2B, 7BDE, 8EE3), 
   LDOUBLE_HEX(BFBE, AF45, EFD8, 2A64, 49A5), 
 /*atan_b[61] ~= 1.5587186337e+00*/
  LDOUBLE_HEX(3FFF, C784, 1799, 9E5D, D2A5), 
   LDOUBLE_HEX(BFBE, A231, BD90, F170, 34A5), 
};
 static const long double coef_poly[9][2] = {
{ -3.33333333333333333342368351437379203616728773340583e-01L,  9.03501810404587028364033466367082415937499719525463e-21L},
{ 2.00000000000000000002710505431213761085018632002175e-01L,  -2.71050543121376108505536620063805076318847614178820e-21L},
{ -1.42857142857142857140921067549133027796415262855589e-01L,  -1.93607530800982934641564128836546985281459293443700e-21L},
{ 1.11111111111111111109605274760436799397211871109903e-01L,  1.50583635067431171387883211317314321885579450456211e-21L},
{ -9.09090909090909090933731867556488737136533018201590e-02L, 0},
{ 7.69230769230769230779655790120052927250071661546826e-02L, 0},
{ -6.66666666666666666698289230030827212658550706692040e-02L, 0},
{ 5.88235294117647058825522430464127765503690170589834e-02L, 0},
{ -5.26315789473684210515616425929419364138084347359836e-02L, 0},
}; 




extern double atan_rn(double xd) {

  unsigned int hx;
  double sign;
  double u;
  double comp;

  int i, i1, m;
  UINT64 x_val,x_abs,sign_mask;
  L_FLOAT_TYPE xe, tmp, bi, atanbhi, xred, xred2,q;
  L_FLOAT_TYPE res,reshi,reslo,rn_constant,test;
  L_FLOAT_TYPE xred4,tmp2;    
  L_FLOAT_TYPE a,b,e0,e1,e2,e3,q0,q1,q2,y0,y1,y2,xred2coarse;
  L_FLOAT_TYPE C3,C5,C7,C9 ;

  
  x_val = _Asm_getf( _FR_D, xd );
  x_abs = (x_val & ULL(7fffffffffffffff));
  sign_mask = ((SINT64)x_val >> 63); /* either 00..00 or 11...11 */



  /* cast x to a DE register */
  if(sign_mask) 
    xe=-xd;
  else
    xe=xd;


  /* Filter cases */
  if (__builtin_expect( x_abs >= ULL(4350000000000000), 0)) {           /* x >= 2^54 */
    if (xd!=xd )
      return xd+xd;                /* NaN */
    else {/* atan(x) = +/- Pi/2 */
      if(sign_mask) return -HALFPI; else return HALFPI;
    }
  }
  else if (__builtin_expect( x_abs < ULL(3E40000000000000), 0))
    /* TODO Add stuff to raise inexact flag */ 
    return xd;                   /* x<2^-27 then atan(x) =~ x */


  /* Now there is something to compute*/
  
  /* load polynomial coeffs */
  C3=coef_poly[0][0];
  C5=coef_poly[1][0];
  C7=coef_poly[2][0];
  C9=coef_poly[3][0];
  
  if (__builtin_expect(x_abs > MIN_REDUCTION_NEEDED, 0)) /* test if reduction is necessary : */
    {
      /* 1) Argument reduction :  */
      /* This constant was found by dichotomy. I am very ashamed */
      rn_constant = 1.002;
      
      /* compute i so that a[i] < x < a[i+1] */

     if (x_abs>ab_table[61].a)
        i=61;
     else {
       i=31;
        if (x_abs < ab_table[i].a) i-= 16;
        else i+=16;
        if (x_abs < ab_table[i].a) i-= 8;
        else i+= 8;
        if (x_abs < ab_table[i].a) i-= 4;
        else i+= 4;
        if (x_abs < ab_table[i].a) i-= 2;
        else i+= 2;
        if (x_abs < ab_table[i].a) i-= 1;
        else i+= 1;
        if (x_abs < ab_table[i].a) i-= 1;
     }
          
     bi= ab_table[i].b;
     atanbhi = atanb_table[i].hi;
     
     /* the dividend and the divisor for the argument reduction */
     a = xe-bi;       b = 1 +  xe * bi;
     

#if 1
     /* now we want to compute  (xe - bi )/b as a DE, but 
	we will need the accurate quotient  only later on, 
	we can start the computation of the polynomial with a much coarser approximation.
	Saves 12 cycles.
      */
     /* Algo 8.11 in Markstein book */
     _Asm_frcpa(&y0, a, b, _SF1);
     
      e0 = 1 - b*y0;       q0 = a*y0;
      e2 = e0 + e0*e0;     e1 = e0*e0;
      e3 = e0 + e1*e1;     q1 = q0+q0*e2;        
      xred = q0 + q1*e3;   xred2coarse = q1*q1;   /* 62 bits in xred, more than enough */
      xred2 = xred*xred;   xred4 = xred2coarse*xred2coarse; 



      /*polynom  evaluation */

      tmp2 = C7 + xred2coarse * C9 ;

      /* here we need xred2,  xred2coarse loses a lot of precision to win 3 cycles. */
      tmp = C3 + xred2 * C5;

      q = tmp + xred4 * tmp2;


#else
      xred=a/b; 
      xred2=xred*xred;
      xred4=xred2*xred2;
      tmp2 = C7 + xred2 * C9 ;
      tmp = C3 + xred2 * C5;
      q = tmp + xred4 * tmp2;
#endif

      tmp = 1+q*xred2;
      /* reconstruction : atan(x) = atan(b[i]) + atan(x) */
      res = atanbhi+xred*tmp;
      /*      atan = atanbhi + tmp;  with round to double */
      reshi = _Asm_fma( _PC_D, xred, tmp, atanbhi, _SF0 );

    }
  else
    /* no reduction needed */
    {


      /* Polynomial evaluation */
      
      xred2 = xe*xe;
      /*poly eval */
      xred4=xred2*xred2;
      tmp2 = C7 + xred2 * C9 ;
      tmp = C3 + xred2 * C5;
      q = tmp + xred4 * tmp2;
      q *= xred2;
      

      /* compute q*xe+xe with round to double */
      res = _Asm_fma( _PC_NONE, q, xe, xe, _SF1 );
      reshi = _Asm_fma( _PC_D, q, xe, xe, _SF0 );
    }
  
#if 0 /* To time the first step only */
  if(sign_mask) 
    return -reshi; 
  else 
    return reshi;
#endif
  
#if 1
    i1 = _Asm_getf( _FR_SIG, res);
    m =  i1 & (0xff<<3);
    if(__builtin_expect((m!=(0x7f<<3) && m!=(0x80<<3)), 1+1==2)) {
      if(sign_mask) 
	return -reshi; 
      else 
	return reshi;
    }      
#else
  /* ROUNDING TEST à la Ziv */
      /* This constant was found by dichotomy. I am very ashamed */
      rn_constant = 1.01;
  reslo = res - reshi;
  test=_Asm_fma( _PC_D, reslo, rn_constant, reshi, _SF0 );
  
  if (__builtin_expect(reshi == test, 1+1==2)) {
    if(sign_mask) 
      return -reshi; 
    else 
      return reshi;
  }
#endif


  else {

    /******************************************************************/
    /*  Double-double-extended  */
    L_FLOAT_TYPE tmphi, tmplo, x0hi, x0lo, xmBihi, xmBilo, Xredhi, Xredlo, Xred2, qhi,qlo, q, Xred2hi,Xred2lo, atanhi,atanlo;
    int j;




#if EVAL_PERF
    crlibm_second_step_taken++;
#endif

#if DEBUG
       printf("Toto\n");
#endif

  if (__builtin_expect(x_abs > MIN_REDUCTION_NEEDED, 0))  {/* test if reduction is necessary : */
    if(i==61){
      Add12_ext( xmBihi , xmBilo , xe , -ab_table[61].b);
    }
    else {
      xmBihi = xe-ab_table[i].b;
      xmBilo = 0.0;
    }
    
    Mul12_ext(tmphi,tmplo, xe, (ab_table[i].b));

    if (xe > 1) /* TODO remplacer par xabs */
      Add22_ext(x0hi,x0lo,tmphi,tmplo, 1.0,0.0);
    else {
      Add22_ext(x0hi , x0lo , 1.0,0.0,tmphi,tmplo);
    }
    
#if 1
    Div22_ext(Xredhi, Xredlo, xmBihi , xmBilo , x0hi,x0lo);
#else
    Xredhi=1; Xredlo=0; /* to time the Div22*/
#endif
    
#if DEBUG
    printf("i=%d, num=%1.15e   den=%1.15e\n",i, (double)xmBihi, (double)x0hi);
    printf("Xred=%1.15e\n", (double)Xredhi);
#endif

    Xred2 = Xredhi*Xredhi;
    Mul22_ext(Xred2hi,Xred2lo,Xredhi,Xredlo,Xredhi, Xredlo);
    
      /*poly eval */
      
    q = (coef_poly[4][0]+Xred2*
	 (coef_poly[5][0]+Xred2*
	  (coef_poly[6][0]+Xred2*
	   (coef_poly[7][0]+
	    (Xred2*coef_poly[8][0])))));
      
    Mul12_ext(qhi, qlo, q, Xred2);
    
    for(j=3;j>=0;j--)
      {
	Add22_ext(qhi,qlo, (coef_poly[j][0]), (coef_poly[j][1]), qhi,qlo);
	Mul22_ext(qhi,qlo, qhi,qlo, Xred2hi,Xred2lo);
      }
      
    Mul22_ext(qhi,qlo, Xredhi,Xredlo, qhi,qlo);
    Add22_ext(qhi,qlo, Xredhi,Xredlo, qhi,qlo);
    
    /* reconstruction : atan(x) = atan(b[i]) + atan(x) */
    Add22_ext(atanhi,atanlo, atanb_table[i].hi, atanb_table[i].lo, qhi,qlo);
  }
  else
    /* no reduction needed */
    {

#if DEBUG
       printf("Tata\n");
#endif
      /* Polynomial evaluation */
      Mul12_ext( Xred2hi,Xred2lo,xe,xe);

      /*poly eval - don't take risks, keep plain Horner */

      q = coef_poly[8][0];
      q = coef_poly[7][0]+Xred2hi*q;
      q = coef_poly[6][0]+Xred2hi*q;
      q = coef_poly[5][0]+Xred2hi*q;
      
      Add12_ext(qhi,qlo, coef_poly[4][0], Xred2hi*q);
#if DEBUG
      printf("    qhi+ql = %1.50Le + %1.50Le\n",(long double)qhi, (long double)qlo);
      print_debug("qhi", qhi);
      print_debug("qlo", qlo);
#endif
      Mul22_ext(qhi,qlo, qhi,qlo, Xred2hi,Xred2lo);
#if DEBUG
      printf("    Xred2  = %1.50Le + %1.50Le\n",(long double)Xred2hi, (long double)Xred2lo);
      printf("    qhi+ql = %1.50Le + %1.50Le\n",(long double)qhi, (long double)qlo);
      print_debug("qhi", qhi);
      print_debug("qlo", qlo);
#endif
      
      for(j=3;j>=0;j--)
        {
          Add22_ext(qhi,qlo, (coef_poly[j][0]), (coef_poly[j][1]), qhi,qlo);
          Mul22_ext(qhi,qlo, qhi,qlo, Xred2hi,Xred2lo);
        }
      
      Mul22_ext (qhi,qlo, xe,0, qhi,qlo);

#if DEBUG
      printf("    qhi+ql = %1.50Le + %1.50Le\n",(long double)qhi, (long double)qlo);
      print_debug("qhi", qhi);
      print_debug("qlo", qlo);
#endif
      /* Now comes the addition sequence proven in the TOMS paper */
      Add12_ext(atanhi,atanlo,xe,qhi);
#if DEBUG
      print_debug("atanhi", atanhi);
      printf(" atan hi+lo  %1.50Le + %1.50Le\n",(long double)atanhi, (long double)atanlo);
#endif
      atanlo += qlo;

    }

#if DEBUG
  printf(" atan hi+lo  %1.50Le + %1.50Le\n",(long double)atanhi, (long double)atanlo);
  printf("             %1.50e + %1.50e\n",(double)atanhi,(double) atanlo);
  printf("             %1.50Le\n",(long double)(atanhi + atanlo));
  printf("             ");
#endif
  
    if(sign_mask) 
      res= -(double) (atanhi+atanlo);
    else 
      res= (double) (atanhi+atanlo);

    return res;
    
  }
} 

