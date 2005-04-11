/*
 * Correctly rounded exponential
 *
 * Author: David Defour
 *
 * This file is part of the crlibm library developed by the Arenaire
 * project at Ecole Normale Superieure de Lyon
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
*/
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


/* For non-C99 compilers maybe which don't handle hex floats, maybe we
   should have some #ifs there */
static const double largest_double = 0x1.fffffffffffffp1023;
static const double tiniest_double = 0x1.0p-1074;






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

/*
 * Evaluate the exponential of x 
 * (res_hi + res_lo)*2^k = exp(x)
 * used for the exp, pow function
 */
void exp_quick(db_number * res_hi, db_number * res_lo, int *k, double x){
  double R1, R2, R3_hi, R3_lo, R4, R5_hi, R5_lo, R6, R7, R9, R10, crp_hi;
  double rp_hi, rp_lo, ex_hi, ex_lo;
  double r_hi, r_lo;
  double P_r;
  db_number db;
  int indx;

  /* Arrondi au plus près */
  DOUBLE2INT(*k, (x * inv_ln2));

  if (k != 0){
    /* r_hi+r_lo =  x - (ln2_hi + ln2_me + ln2_lo)*k */
    rp_hi =  x-ln2_hi*(*k);
    rp_lo =   -ln2_lo*(*k);
    Add12Cond(r_hi, r_lo, rp_hi, rp_lo);
  }else {
    r_hi = x;  r_lo = 0.;
  }


  /*
   * 2) Deuxième reduction d'argument
   */

  /* Arrondi au plus près */
  db.d  = (r_hi + two_44_43);
  indx  = db.i[LO];
  indx += BIAS;
  r_hi  -= (db.d - two_44_43);

  /* Normalisation du résultat */
  Add12(rp_hi, rp_lo, r_hi, r_lo);

  /* Lecture de table */
  ex_hi = (tab_exp[indx][0]).d;
  ex_lo = (tab_exp[indx][1]).d;


  /*
   * 3) Evaluation Polynomiale
   */

  P_r = (c_0 + rp_hi * (c_1 + rp_hi * (c_2 + (rp_hi * c_3))));


  R1 = rp_hi * rp_hi;

  crp_hi = R1 * rp_hi;
  R1 *= 0.5;

  R2 =  P_r * crp_hi;

  Mul12(&R3_hi, &R3_lo, ex_hi, rp_hi);
  R4 = ex_hi * rp_lo;

  Mul12(&R5_hi, &R5_lo, ex_hi, R1);
  R6 = (ex_hi * rp_lo) + (ex_lo * (R1 + rp_hi));

  R7  = ex_hi * R2;
  R7 += (R6 + R5_lo) + (R3_lo + ex_lo);

  Add12(R9, (*res_lo).d, R7, R5_hi);

  Add12(R10, R1, R3_hi, R9);
  (*res_lo).d += R1;

  Add12((*res_hi).d, R1, ex_hi, R10);                                  
  (*res_lo).d += R1;                                                    

  Add12((*res_hi).d, (*res_lo).d, (*res_hi).d, (*res_lo).d);                                
}



/***************************
 ***************************
 ***  ROUDING TO NEAREST ***
 ***************************
 ***************************/
double exp_rn(double x){
  db_number db, R1, R8, R11;
  unsigned int hx;
  int k;
  
  db_number st_to_mem;
  int    exp_R11;



  /*
   * 1) First argument reduction
   */

  db.d = x;
  hx   = db.i[HI]; 
  hx  &= 0x7fffffff;  

  /* Filter special cases */
  if (hx >= 0x40862E42){
    if (hx >= 0x7ff00000){
      if (((hx&0x000fffff)|db.i[LO])!=0)
	return x+x;                                                  /* Nan */ 
      else return ((db.i[HI]&0x80000000)==0)? x:0.0; /* exp(+/-inf) = inf,0 */
    }
    if (x > o_bound) return largest_double * largest_double;   /* overflow  */ 
    if (x < u_bound) return tiniest_double * tiniest_double;   /* underflow */ 
  }

  if (hx < 0x3C900000)                                /* if (hx <= 2^(-54)) */
    return ((hx == 0) && (db.i[LO] == 0))? 1.: 1.+tiniest_double;  
    

  /* The evaluation */
  exp_quick(&R11, &R8, &k, x);


  /* Résultat = (R11 + R8) */
  exp_R11 = (R11.i[HI] & 0x7ff00000) - errd;

  if (R11.d == (R11.d + (R8.d * errn))){
    if (k > -1020){                                               
      if (k < 1020){                                              
	R11.i[HI] += (k<<20);                              
	return R11.d;                                             
      }else {                                                     
	/* On est proche de + Inf */                              
	R11.i[HI] += ((k-1000)<<20);                       
	return R11.d*two1000;                                     
      }                                                           
    }else {                                                       
      /* result is a subnorml number  */                          
      R11.i[HI] += ((k+1000)<<20);    
      db.d = R11.d * twom1000;

      /*
       * We are working on adress to force the compiler and data
       * to transit throught memory and avoid extra precision.
       */
      st_to_mem.i[HI] = db.i[HI];
      st_to_mem.i[LO] = db.i[LO];
      
      R11.d -= st_to_mem.d * two1000;
      R1.i[HI] = R11.i[HI] & 0x7fffffff;
      R1.i[LO] = R11.i[LO];

      if (R1.d == two_m75){
	if ((R8.i[HI] & 0x7ff00000) < exp_R11){
	  /* Arrondi difficile ! */
	  return scs_exp_rn(x);
	}
	/* The error  is exactly 1/2 ulp of the result */
	if ((R11.i[HI] > 0)&&(R8.i[HI] > 0)){
	  /* st_to_mem                   */
	  /*          R11.d    R8        */ 
	  /*   |---| |10--0| |0----01--| */
	  /* We need to add 1 ulp        */
	  db.l += 1;
	}else
	if ((R11.i[HI] < 0)&&(R8.i[HI] < 0)){
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
  db_number db, R8, R11;
  unsigned int hx;
  int k;
  
  db_number st_to_mem;
  int    exp_R11;

  /*
   * 1) Première réduction d'argument
   */

  db.d = x;
  hx   = db.i[HI]; 
  hx  &= 0x7fffffff;  

  /* Filtre les cas spéciaux */
  if (hx >= 0x40862E42){
    if (hx >= 0x7ff00000){
      if (((hx&0x000fffff)|db.i[LO])!=0)
	return x+x;                                        /* Nan */ 
      else return ((db.i[HI]&0x80000000)==0)? x:0.0;/* exp(+/-inf) = inf,0 */
    }
    if (x > o_bound) return largest_double * largest_double;       /* overflow  */ 
    if (x < u_bound) return tiniest_double*(1.+tiniest_double);  /* 2^(-1074) */ 
   }

  if (hx < 0x3CA00000){                                    /* if (hx <= 2^(-53)) */ 
    if ((hx == 0) && (db.i[LO] == 0))
      return 1.;                                           /* exp(0)=1. */
    if (db.i[HI] < 0)
      return 1. + tiniest_double;                             /* 1 and inexact */
    else
      return 1. + two_m52_56;                              /* 1 + 2^(-52) and inexact */
  }

  /* The evaluation */
  exp_quick(&R11, &R8, &k, x);

  exp_R11   = (R11.i[HI] & 0x7ff00000) - errd;

  if ((R8.i[HI] & 0x7ff00000) > exp_R11){
    /* On est capable d'arrondir */
    if (k > -1020){                                         
      if (k < 1020){                                              
	R11.i[HI] += (k<<20);                              
      }else {                                                     
	/* On est proche de + Inf */                              
	R11.i[HI] += ((k-1000)<<20);                       
	R11.d *= two1000;                                     
      }
      if (R8.i[HI] > 0)
	R11.l += 1; /* Be carefull this work only if R11>0 see David PhD */
                    /* Which is the case for exp(x) !! */
      return R11.d;
    }else {                                                       
      /* On est dans les dénormalisés */     
      
      R11.i[HI] += ((k+1000)<<20);                         
      db.d = R11.d * twom1000;

      /*
       * We are working on adress to force the compiler and data
       * to transit throught memory and avoid extra precision.
       */
      st_to_mem.i[HI] = db.i[HI];
      st_to_mem.i[LO] = db.i[LO];

      R11.d -= st_to_mem.d * two1000;

      /* Be carefull this work only if R11>0 see David PhD */
      /* Which is the case for exp(x) !! */
      if ((R11.i[HI]  > 0)||
	  ((R11.i[HI] == 0)&&(R8.i[HI] > 0))) db.l += 1;
      
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
  db_number db, R8, R11;
  unsigned int hx;
  int k;
  
  db_number st_to_mem;
  int    exp_R11;

  /*
   * 1) Première réduction d'argument
   */
  db.d = x;
  hx   = db.i[HI];
  hx  &= 0x7fffffff;  

  /* Filtre les cas spéciaux */
  if (hx >= 0x40862E42){
    if (hx >= 0x7ff00000){
      if (((hx&0x000fffff)|db.i[LO])!=0)
	return x+x;                                      /* Nan */ 
      else return ((db.i[HI]&0x80000000)==0)? x:0.0;           /* exp(+/-inf) = inf,0 */
    }
    if (x > o_bound) return largest_double * (1.+tiniest_double); /* (1-2^-53).2^1024  */ 
    if (x < u_bound) return tiniest_double * tiniest_double;   /* underflow */ 
  }

  if (hx < 0x3CA00000){                                  /* if (hx<=2^(-53)) */
    if ((hx == 0) && (db.i[LO] == 0)) 
      return 1.;                                         /* exp(0)=1. */
    if (db.i[HI] < 0)
      return 1. - two_m53_56;                            /* 1-2^(-53) and inexact */
    else
      return 1. + tiniest_double;                           /* 1 and inexact         */
  }

  /* The evaluation */
  exp_quick(&R11, &R8, &k, x);

  exp_R11   = (R11.i[HI] & 0x7ff00000) - errd;

  if ((R8.i[HI] & 0x7ff00000) > exp_R11){
    /* On est capable d'arrondir */
    if (k > -1020){                                         
      if (k < 1020){                                              
	R11.i[HI] += (k<<20);                              
      }else {                                                     
	/* On est proche de + Inf */                              
	R11.i[HI] += ((k-1000)<<20);                       
	R11.d *= two1000;                                     
      }
      if (R8.i[HI] < 0)
	R11.l -= 1; /* Be carefull this work only if R11>0 see David PhD */
                    /* Which is the case for exp(x) !! */
      return R11.d;
    }else {                                                       
      /* Subnormal here */     
      
      R11.i[HI] += ((k+1000)<<20);                         
      db.d = R11.d * twom1000;

      /*
       * We are working on adress to force the compiler and data
       * to transit throught memory and avoid extra precision.
       */
      st_to_mem.i[HI] = db.i[HI];
      st_to_mem.i[LO] = db.i[LO];

      R11.d -= st_to_mem.d * two1000;

      /* Be carefull this work only if R11>0 see David PhD */
      /* Which is the case for exp(x) !! */
      if ((R11.i[HI]  < 0)||
	  ((R11.i[HI] == 0)&&(R8.i[HI] < 0))) db.l -= 1;
      
      return db.d;                                      
    }
  }else {
    /* Difficult to round */
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
 
