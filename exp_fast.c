#include <stdio.h>
#include <stdlib.h>
#include "crlibm.h"
#include "crlibm_private.h"
#include "exp_fast.h"


/***********************************************************/
/*        Fast first step for the exponential              */
/***********************************************************/




#define BIAS 89



/* useful constants */
static const double 
 two1000  = 0x1p1000, /* 0x7E700000, 0x00000000, 1.07150860718626732095e301 */
 twom1000 = 0x1p-1000,/* 0x01700000, 0x00000000, 9.33263618503218878990e-302 */
 ln2_hi   = 6.93147180559890330187e-01, /* 0x3FE62E42, 0xFEFA3800 */ 
 ln2_lo   = 5.49792301870837115524e-14, /* 0x3D2EF357, 0x93C76730 */
 inv_ln2  = 1.44269504088896338700e+00;

static const db_number 
#ifdef WORDS_BIGENDIAN
 scs_huge ={{0x7fefffff, 0xffffffff}},
 scs_small={{0x00000000, 0x00000001}};
#else
 scs_huge ={{0xffffffff, 0x7fefffff}},
 scs_small={{0x00000001, 0x00000000}};
#endif





/***********************************************************/
/*                      First step                         */
/***********************************************************/

/* 1) part */
static const double 
 o_bound = 7.09782712893383973096e+02,  /* 0x40862E42, 0xFEFA39EF */
 u_bound =-7.45133219101941108420e+02;  /* 0xC0874910, 0xD52D3051 */



/* 2) part */
static const double two_44_43 = 26388279066624.;      /* 0x42B80000 0x00000000 */   /* 2^44 + 2^43 */
static const double two_m52_56 = 2.35922392732845764840e-16; /* 0x3CB10000 0x00000000 */
static const double two_m53_56 = 1.24900090270330110798e-16; /* 0x3CA20000 0x00000000 */

/* 3) part */
static const double 
 c_0 = 1.66666666666665769236e-01,      /* 0x3FC55555, 0x55555535 */
 c_1 = 4.16666666666664631257e-02,      /* 0x3FA55555, 0x55555538 */
 c_2 = 8.33333427943885873823e-03,      /* 0x3F811111, 0x31931950 */
 c_3 = 1.38888903080471677251e-03;      /* 0x3F56C16C, 0x3DC3DC5E */


/* 4) part */
static const double
 two_m75  = 2.6469779601696885595885078e-23; /* 0x3B400000, 0x00000000 */

static const db_number 
#ifdef WORDS_BIGENDIAN
 _errn ={{0x3FF00020, 0x00000000}};
#else
 _errn ={{0x00000000, 0x3FF00020}};
#endif
#define errn _errn.d


static const int errd  = 73400320;                   /* 70 * 2^20 */  


#define EVAL_EXP_DECLARATION \
  db_number db, R1, R8, R11;\
  double rp_hi, rp_lo, ex_hi, ex_lo;\
  double r_hi, r_lo;\
  double tmp;\
  double P_r;\
  double R2, R3_hi, R3_lo, R4, R5_hi, R5_lo, R6, R7, R9, R10, crp_hi;\
  unsigned int hx;\
  int k;\
  int indx;


#define EVAL_EXP_STEP_2_AND_3 \
{                                                                 \
  /* Arrondi au plus près */                                      \
  DOUBLE2INT(k, (x * inv_ln2));                                   \
                                                                  \
  if (k != 0){                                                    \
    /* r_hi+r_lo =  x - (ln2_hi + ln2_me + ln2_lo)*k */           \
    rp_hi =  x-ln2_hi*k;                                          \
    rp_lo =   -ln2_lo*k;                                          \
    Add12Cond(r_hi, r_lo, rp_hi, rp_lo);                       \
  }else {                                                         \
    r_hi = x;  r_lo = 0.;                                         \
  }                                                               \
                                                                  \
                                                                  \
  /*                                                              \
   * 2) Deuxième reduction d'argument                             \
   */                                                             \
                                                                  \
                                                                  \
  /* Arrondi au plus près */                                      \
  db.d  = (r_hi + two_44_43);                                     \
  indx  = db.i[LO_ENDIAN];                                        \
  indx += BIAS;                                                   \
  r_hi  -= (db.d - two_44_43);                                    \
                                                                  \
  /* Normalisation du résultat */                                 \
  Add12(rp_hi, rp_lo, r_hi, r_lo);                             \
                                                                  \
  /* Lecture de table */                                          \
  ex_hi = (tab_exp[indx][0]).d;                                   \
  ex_lo = (tab_exp[indx][1]).d;                                   \
                                                                  \
                                                                  \
  /*                                                              \
   * 3) Evaluation Polynomiale                                    \
   */                                                             \
                                                                  \
  P_r = (c_0 + rp_hi * (c_1 + rp_hi * (c_2 + (rp_hi * c_3))));    \
                                                                  \
                                                                  \
  R1.d = rp_hi * rp_hi;                                           \
                                                                  \
  crp_hi = R1.d * rp_hi;                                          \
  /*R1.i[HI_ENDIAN] -= 0x00100000;*/                              \
  R1.d *= 0.5;                                                    \
                                                                  \
  R2 =  P_r * crp_hi;                                             \
                                                                  \
  Mul12(&R3_hi, &R3_lo, ex_hi, rp_hi);                           \
  R4 = ex_hi * rp_lo;                                             \
                                                                  \
  Mul12(&R5_hi, &R5_lo, ex_hi, R1.d);                            \
  R6 = R4 + (ex_lo * (R1.d + rp_hi));                             \
                                                                  \
  R7  = ex_hi * R2;                                               \
  R7 += (R6 + R5_lo) + (R3_lo + ex_lo);                           \
                                                                  \
  Add12(R9, R8.d, R7, R5_hi);                                  \
                                                                  \
  Add12(R10, tmp, R3_hi, R9);                                  \
  R8.d += tmp;                                                    \
                                                                  \
  Add12(R11.d, tmp, ex_hi, R10);                               \
  R8.d += tmp;                                                    \
                                                                  \
  Add12(R11.d, R8.d, R11.d, R8.d);                             \
                                                                  \
} 


/***************************
 ***************************
 ***  ROUDING TO NEAREST ***
 ***************************
 ***************************/
double exp_rn(double x){
  EVAL_EXP_DECLARATION

  db_number st_to_mem;
  int    exp_R11;


  /*
   * 1) Première réduction d'argument
   */

  db.d = x;
  hx   = db.i[HI_ENDIAN]; 
  hx  &= 0x7fffffff;  

  /* Filtre les cas spéciaux */
  if (hx >= 0x40862E42){
    if (hx >= 0x7ff00000){
      if (((hx&0x000fffff)|db.i[LO_ENDIAN])!=0)
	return x+x;                                        /* Nan */ 
      else return ((db.i[HI_ENDIAN]&0x80000000)==0)? x:0.0;/* exp(+/-inf) = inf,0 */
    }
    if (x > o_bound) return scs_huge.d * scs_huge.d;       /* overflow  */ 
    if (x < u_bound) return scs_small.d * scs_small.d;     /* underflow */ 
  }

  if (hx < 0x3C900000)                                     /* if (hx <= 2^(-54)) */
    return ((hx == 0) && (db.i[LO_ENDIAN] == 0))? 1.: 1.+scs_small.d;  
    
  EVAL_EXP_STEP_2_AND_3;

  /* Résultat = (R11 + R8) */
  exp_R11 = (HI(R11) & 0x7ff00000) - errd;

  if (R11.d == (R11.d + (R8.d * errn))){
    if (k > -1020){                                               
      if (k < 1020){                                              
	R11.i[HI_ENDIAN] += (k<<20);                              
	return R11.d;                                             
      }else {                                                     
	/* On est proche de + Inf */                              
	R11.i[HI_ENDIAN] += ((k-1000)<<20);                       
	return R11.d*two1000;                                     
      }                                                           
    }else {                                                       
      /* On est dans les dénormalisés */                          
      R11.i[HI_ENDIAN] += ((k+1000)<<20);    
      db.d = R11.d * twom1000;

      /*
       * We are working on adress to force the compiler and data
       * to transit throught memory and avoid extra precision.
       */
      st_to_mem.i[HI_ENDIAN] = db.i[HI_ENDIAN];
      st_to_mem.i[LO_ENDIAN] = db.i[LO_ENDIAN];
      
      R11.d -= st_to_mem.d * two1000;
      R1.i[HI_ENDIAN] = R11.i[HI_ENDIAN] & 0x7fffffff;
      R1.i[LO_ENDIAN] = R11.i[LO_ENDIAN];

      if (R1.d == two_m75){
	if ((HI(R8) & 0x7ff00000) < exp_R11){
	  /* Arrondi difficile ! */
	  return scs_exp_rn(x);
	}
	/* The error  is exactly 1/2 ulp of the result */
	if ((R11.i[HI_ENDIAN] > 0)&&(R8.i[HI_ENDIAN] > 0)){
	  /* st_to_mem                   */
	  /*          R11.d    R8        */ 
	  /*   |---| |10--0| |0----01--| */
	  /* We need to add 1 ulp        */
	  db.l += 1;
	}else
	if ((R11.i[HI_ENDIAN] < 0)&&(R8.i[HI_ENDIAN] < 0)){
	  /* st_to_mem    R11.d         R8       */ 
	  /*  |--(+1)| |(-1)0--0| |0----0(-1)--| */
	  /*                                     */
	  /* We need to remove 1 ulp             */
	  db.l -= 1;
	}
      }
      return db.d;                                      
    }                                                             
  }else {
    /* Cas difficile */
    return scs_exp_rn(x);
  }
}


/***************************
 ***************************
 ***  ROUDING TO + INF   ***
 ***************************
 ***************************/
double exp_ru(double x){
  EVAL_EXP_DECLARATION

  db_number st_to_mem;
  int           exp_R11;

  /*
   * 1) Première réduction d'argument
   */

  db.d = x;
  hx   = db.i[HI_ENDIAN]; 
  hx  &= 0x7fffffff;  

  /* Filtre les cas spéciaux */
  if (hx >= 0x40862E42){
    if (hx >= 0x7ff00000){
      if (((hx&0x000fffff)|db.i[LO_ENDIAN])!=0)
	return x+x;                                        /* Nan */ 
      else return ((db.i[HI_ENDIAN]&0x80000000)==0)? x:0.0;/* exp(+/-inf) = inf,0 */
    }
    if (x > o_bound) return scs_huge.d * scs_huge.d;       /* overflow  */ 
    if (x < u_bound) return scs_small.d*(1.+scs_small.d);  /* 2^(-1074) */ 
   }

  if (hx < 0x3CA00000){                                    /* if (hx <= 2^(-53)) */ 
    if ((hx == 0) && (db.i[LO_ENDIAN] == 0))
      return 1.;                                           /* exp(0)=1. */
    if (db.i[HI_ENDIAN] < 0)
      return 1. + scs_small.d;                             /* 1 and inexact */
    else
      return 1. + two_m52_56;                              /* 1 + 2^(-52) and inexact */
  }

  EVAL_EXP_STEP_2_AND_3;

  exp_R11   = (R11.i[HI_ENDIAN] & 0x7ff00000) - errd;

  if ((R8.i[HI_ENDIAN] & 0x7ff00000) > exp_R11){
    /* On est capable d'arrondir */
    if (k > -1020){                                         
      if (k < 1020){                                              
	R11.i[HI_ENDIAN] += (k<<20);                              
      }else {                                                     
	/* On est proche de + Inf */                              
	R11.i[HI_ENDIAN] += ((k-1000)<<20);                       
	R11.d *= two1000;                                     
      }
      if (R8.i[HI_ENDIAN] > 0)
	R11.l += 1; /* Be carefull this work only if R11>0 see David PhD */
                    /* Which is the case for exp(x) !! */
      return R11.d;
    }else {                                                       
      /* On est dans les dénormalisés */     
      
      R11.i[HI_ENDIAN] += ((k+1000)<<20);                         
      db.d = R11.d * twom1000;

      /*
       * We are working on adress to force the compiler and data
       * to transit throught memory and avoid extra precision.
       */
      st_to_mem.i[HI_ENDIAN] = db.i[HI_ENDIAN];
      st_to_mem.i[LO_ENDIAN] = db.i[LO_ENDIAN];

      R11.d -= st_to_mem.d * two1000;

      /* Be carefull this work only if R11>0 see David PhD */
      /* Which is the case for exp(x) !! */
      if ((R11.i[HI_ENDIAN]  > 0)||
	  ((R11.i[HI_ENDIAN] == 0)&&(R8.i[HI_ENDIAN] > 0))) db.l += 1;
      
      return db.d;                                      
    }                                                             
  }else {
    /* Cas difficile */
    return scs_exp_ru(x);
  }
}




/***************************
 ***************************
 ***  ROUDING TO - INF   ***
 ***************************
 ***************************/

double exp_rd(double x){
  EVAL_EXP_DECLARATION
  
  db_number st_to_mem;
  int    exp_R11;

  /*
   * 1) Première réduction d'argument
   */
  db.d = x;
  hx   = db.i[HI_ENDIAN];
  hx  &= 0x7fffffff;  

  /* Filtre les cas spéciaux */
  if (hx >= 0x40862E42){
    if (hx >= 0x7ff00000){
      if (((hx&0x000fffff)|db.i[LO_ENDIAN])!=0)
	return x+x;                                      /* Nan */ 
      else return ((hx&0x80000000)==0)? x:0.0;           /* exp(+/-inf) = inf,0 */
    }
    if (x > o_bound) return scs_huge.d*(1.+scs_small.d); /* (1-2^-53).2^1024  */ 
    if (x < u_bound) return scs_small.d * scs_small.d;   /* underflow */ 
  }

  if (hx < 0x3CA00000){                                  /* if (hx<=2^(-53)) */
    if ((hx == 0) && (db.i[LO_ENDIAN] == 0)) 
      return 1.;                                         /* exp(0)=1. */
    if (db.i[HI_ENDIAN] < 0)
      return 1. - two_m53_56;                            /* 1-2^(-53) and inexact */
    else
      return 1. + scs_small.d;                           /* 1 and inexact         */
  }

  EVAL_EXP_STEP_2_AND_3;

  exp_R11   = (R11.i[HI_ENDIAN] & 0x7ff00000) - errd;

  if ((R8.i[HI_ENDIAN] & 0x7ff00000) > exp_R11){
    /* On est capable d'arrondir */
    if (k > -1020){                                         
      if (k < 1020){                                              
	R11.i[HI_ENDIAN] += (k<<20);                              
      }else {                                                     
	/* On est proche de + Inf */                              
	R11.i[HI_ENDIAN] += ((k-1000)<<20);                       
	R11.d *= two1000;                                     
      }
      if (R8.i[HI_ENDIAN] < 0)
	R11.l -= 1; /* Be carefull this work only if R11>0 see David PhD */
                    /* Which is the case for exp(x) !! */
      return R11.d;
    }else {                                                       
      /* On est dans les dénormalisés */     
      
      R11.i[HI_ENDIAN] += ((k+1000)<<20);                         
      db.d = R11.d * twom1000;

      /*
       * We are working on adress to force the compiler and data
       * to transit throught memory and avoid extra precision.
       */
      st_to_mem.i[HI_ENDIAN] = db.i[HI_ENDIAN];
      st_to_mem.i[LO_ENDIAN] = db.i[LO_ENDIAN];

      R11.d -= st_to_mem.d * two1000;

      /* Be carefull this work only if R11>0 see David PhD */
      /* Which is the case for exp(x) !! */
      if ((R11.i[HI_ENDIAN]  < 0)||
	  ((R11.i[HI_ENDIAN] == 0)&&(R8.i[HI_ENDIAN] < 0))) db.l -= 1;
      
      return db.d;                                      
    }
  }else {
    /* Cas difficile */
    return scs_exp_rd(x);
  }
}




/***************************
 ***************************
 ***  ROUDING TO ZERO    ***
 ***************************
 ***************************/
/*
 * exp(x) is a function >0 so rounding to 0
 * is equivalent to rounding to - Inf.
 * In crlibm.h we have :
 *
 *#define exp_rz(x)  exp_rd(x);
 */

/*
 * probleme en arrondi au plus près 
 *-7.084061284577078367874491959810256958007812500000000000000000e+02
 *-7.100496760636764292939915321767330169677734375000000000000000e+02
 *
 */
