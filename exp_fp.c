#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include "cr_tools.h"
#include "exp_fp_tbl.h"

#include <fpu_control.h>
#ifndef __setfpucw
#define __setfpucw(cw) __asm__ ("fldcw %0" : : "m" (cw))
#endif 
	

/*
 * Définition des constantes 
 */

/* 1) part */
double 
 huge    = 1.0e+300,
 small   = 1.0e-300,
 o_bound = 7.09782712893383973096e+02,  /* 0x40862E42, 0xFEFA39EF */
 u_bound =-7.45133219101941108420e+02,  /* 0xC0874910, 0xD52D3051 */
 ln2_hi  = 6.93147180559890330187e-01,  /* 0x3FE62E42, 0xFEFA3800 */ 
 ln2_me  = 5.49792301870720995198e-14,  /* 0x3D2EF357, 0x93C76000 */
 ln2_lo  = 1.16122272293625324218e-26,  /* 0x3A8CC01F, 0x97B57A08 */ 
 half[2] = {0.5,
           -0.5},
 inv_ln2 = 1.44269504088896338700e+00;


/* 2) part */
double two_44_43 = 26388279066624.;     /* 0x42B80000 0x00000000 */   /* 2^44 + 2^43 */
int         bias = 89;               


/* 3) part */
double 
 c_0 = 1.66666666666665769236e-01,      /* 0x3FC55555, 0x55555535 */
 c_1 = 4.16666666666664631257e-02,      /* 0x3FA55555, 0x55555538 */
 c_2 = 8.33333427943885873823e-03,      /* 0x3F811111, 0x31931950 */
 c_3 = 1.38888903080471677251e-03;      /* 0x3F56C16C, 0x3DC3DC5E */


/* 4) part */
double
 two1000  = 1.07150860718626732095e301, /* 0x7E700000 0x00000000 */
 twom1000 = 9.33263618503218878990e-302,/* 0x01700000, 0x00000000 */
 twom53   = 1.11022302462515654042e-16, /* 0x3CA00000, 0x00000000 */
 errn     = 1.00012207031250000000e0;   /* 0x3FF00080, 0x00000000 */

int errd  = 71303168;                   /* 60 * 2^20 */  



/***************************
 ***************************
 ***  ROUDING TO NEAREST ***
 ***************************
 ***************************/

double exp_rn(double x){
  /* 1) Variables temporaires */
  double r_hi, r_lo;
  double u, tmp;
  unsigned int hx;
  int k;


  /* 2) Variables temporaires */
  double rp_hi, rp_lo, ex_hi, ex_lo;
  int index;

  /* 3) Variables temporaires */
  double P_r;

  /* 4) Variables temporaires */
  double R1, R2, R3_hi, R3_lo, R4, R5_hi, R5_lo, R6, R7, R8, R9, R10, R11, crp_hi;
  


  /*
   * 1) Première réduction d'argument
   */
  hx  = HI(x);       
  hx &= 0x7fffffff;  

  /* Filtre les cas spéciaux */
  if (hx >= 0x40862E42){
    if (hx >= 0x7ff00000){
      if (((hx&0x000fffff)|LO(x))!=0)
	return x+x;                     /* Nan */ 
      else return ((hx&0x80000000)==0)? x:0.0;      /* exp(+/-inf) = inf,0 */
    }
    if (x > o_bound) return huge*huge;  /* overflow  */ 
    if (x < u_bound) return small*small;/* underflow */ 
  }

  if (hx < 0x3C900000) return 1.; 

  /* Arrondi au plus près */
  DOUBLE2INT(k, x * inv_ln2);

  if (k != 0){ 
    /* r_hi+r_lo =  x - (ln2_hi + ln2_me + ln2_lo)*k */
    rp_hi =  x-ln2_hi*k;
    rp_lo = -ln2_me*k; 
    Fast2SumCond(r_hi, u, rp_hi, rp_lo);
    r_lo = u - ln2_lo*k;                         
  }else {
    r_hi = x;  r_lo = 0.;                        
  }

  /*
   * 2) Deuxième reduction d'argument
   */
  

  /* Arrondi au plus près */
  tmp    = (r_hi + two_44_43);         
  index  = LO(tmp); 
  index += bias;
  r_hi  -= (tmp - two_44_43);
 
  /* Normalisation du résultat */
  Fast2Sum(rp_hi, rp_lo, r_hi, r_lo); 

  /* Lecture de table */
  ex_hi = tab_exp[index][0];         
  ex_lo = tab_exp[index][1];         

  
  /*
   * 3) Evaluation Polynomiale
   */

  P_r = (c_0 + rp_hi * (c_1 + rp_hi * (c_2 + (rp_hi * c_3)))); 


  R1 = rp_hi * rp_hi;              

  crp_hi = R1 * rp_hi;             
  /* Correspond à R1 /= 2; */
  HI(R1) = HI(R1)-0x00100000;      

  R2 =  P_r * crp_hi;              

  /*  (ex_hi+ex_lo)*(rp_hi+rp_lo) */
  Dekker(&R3_hi, &R3_lo, ex_hi, rp_hi);      
  R4 = ex_hi * rp_lo;                      

  Dekker(&R5_hi, &R5_lo, ex_hi, R1);         
  R6 = R4 + (ex_lo * (R1 + rp_hi));      

  R7  = ex_hi * R2;                      
  R7 += (R6 + R5_lo) + (R3_lo + ex_lo);  

  Fast2Sum(R9, R8, R7, R5_hi);           
  
  Fast2Sum(R10, tmp, R3_hi, R9);           
  R8 += tmp;                               

  Fast2Sum(R11, tmp, ex_hi, R10);          
  R8 += tmp;                               

  Fast2Sum(R11, R8, R11, R8);              

  if (R11 = (R11 + (R11 * errn))){
    if (k > -1020){
      if (k < 1020){
	HI(R11) += (k<<20);                  
	return R11;
      }else {
	/* On est proche de + Inf */
	HI(R11) += ((k-1000)<<20);
	return R11*two1000;
      }
    }else {
      /* On est dans les dénormalisés */
      HI(R11) += ((k+1000)<<20);           
      return R11*twom1000;
    }
  }else {
    /* Cas difficile */
    sn_exp(x);
    return 0;
  }
}


/***************************
 ***************************
 ***  ROUDING TO + INF   ***
 ***************************
 ***************************/

double exp_ru(double x){
  /* 1) Variables temporaires */
  double r_hi, r_lo;
  double u, tmp;
  unsigned int hx;
  int k;


  /* 2) Variables temporaires */
  double rp_hi, rp_lo, ex_hi, ex_lo;
  int index;

  /* 3) Variables temporaires */
  double P_r;

  /* 4) Variables temporaires */
  double R1, R2, R3_hi, R3_lo, R4, R5_hi, R5_lo, R6, R7, R8, R9, R10, R11, crp_hi;
  
  double ulp_R11;
  int    exp_R11;

  /*
   * 1) Première réduction d'argument
   */
  hx  = HI(x);       
  hx &= 0x7fffffff;  

  /* Filtre les cas spéciaux */
  if (hx >= 0x40862E42){
    if (hx >= 0x7ff00000){
      if (((hx&0x000fffff)|LO(x))!=0)
	return x+x;                     /* Nan */ 
      else return ((hx&0x80000000)==0)? x:0.0;      /* exp(+/-inf) = inf,0 */
    }
    if (x > o_bound) return huge*huge;  /* overflow  */ 
    if (x < u_bound) return small*small;/* underflow */ 
  }

  if (hx < 0x3C900000) return 1.; 

  /* Arrondi au plus près */
  DOUBLE2INT(k, x * inv_ln2);

  if (k != 0){ 
    /* r_hi+r_lo =  x - (ln2_hi + ln2_me + ln2_lo)*k */
    rp_hi =  x-ln2_hi*k;
    rp_lo = -ln2_me*k; 
    Fast2SumCond(r_hi, u, rp_hi, rp_lo);
    r_lo = u - ln2_lo*k;                         
  }else {
    r_hi = x;  r_lo = 0.;                        
  }

  /*
   * 2) Deuxième reduction d'argument
   */
  

  /* Arrondi au plus près */
  tmp    = (r_hi + two_44_43);         
  index  = LO(tmp); 
  index += bias;
  r_hi  -= (tmp - two_44_43);
 
  /* Normalisation du résultat */
  Fast2Sum(rp_hi, rp_lo, r_hi, r_lo); 

  /* Lecture de table */
  ex_hi = tab_exp[index][0];         
  ex_lo = tab_exp[index][1];         

  
  /*
   * 3) Evaluation Polynomiale
   */

  P_r = (c_0 + rp_hi * (c_1 + rp_hi * (c_2 + (rp_hi * c_3)))); 


  R1 = rp_hi * rp_hi;              

  crp_hi = R1 * rp_hi;             
  /* Correspond à R1 /= 2; */
  HI(R1) = HI(R1)-0x00100000;      

  R2 =  P_r * crp_hi;              

  /*  (ex_hi+ex_lo)*(rp_hi+rp_lo) */
  Dekker(&R3_hi, &R3_lo, ex_hi, rp_hi);      
  R4 = ex_hi * rp_lo;                      

  Dekker(&R5_hi, &R5_lo, ex_hi, R1);         
  R6 = R4 + (ex_lo * (R1 + rp_hi));      

  R7  = ex_hi * R2;                      
  R7 += (R6 + R5_lo) + (R3_lo + ex_lo);  

  Fast2Sum(R9, R8, R7, R5_hi);           
  
  Fast2Sum(R10, tmp, R3_hi, R9);           
  R8 += tmp;                               

  Fast2Sum(R11, tmp, ex_hi, R10);          
  R8 += tmp;                               

  Fast2Sum(R11, R8, R11, R8);              

  exp_R11 = (HI(R11) & 0x7ff00000) - errd;
  ulp_R11 = R11 * twom53;

  if ((HI(R8) & 0x7ff00000) > exp_R11){
    /* On est capable d'arrondir */
    if (HI(R8) > 0){     
      /* arrondi vers + inf */
      HI(ulp_R11) &= 0x7ff00000;
      R11 += ulp_R11;
    }     
    if (k > -1020){                    
      if (k < 1020){                   
	HI(R11) += (k<<20);             
	return R11;
      }else {
	/* On est proche de + Inf */
	HI(R11) += ((k-1000)<<20);      
	return R11*two1000;
      }
    }else {
      /* On est dans les dénormalisés */
      HI(R11) += ((k+1000)<<20);        
      return R11*twom1000;
    }
  }else {
    /* Cas difficile */
    su_exp(x);
  }
}














/***************************
 ***************************
 ***  ROUDING TO - INF   ***
 ***************************
 ***************************/

double exp_rd(double x){
  /* 1) Variables temporaires */
  double r_hi, r_lo;
  double u, tmp;
  unsigned int hx;
  int k;


  /* 2) Variables temporaires */
  double rp_hi, rp_lo, ex_hi, ex_lo;
  int index;

  /* 3) Variables temporaires */
  double P_r;

  /* 4) Variables temporaires */
  double R1, R2, R3_hi, R3_lo, R4, R5_hi, R5_lo, R6, R7, R8, R9, R10, R11, crp_hi;
  
  double ulp_R11;
  int    exp_R11;

  /*
   * 1) Première réduction d'argument
   */
  hx  = HI(x);       
  hx &= 0x7fffffff;  

  /* Filtre les cas spéciaux */
  if (hx >= 0x40862E42){
    if (hx >= 0x7ff00000){
      if (((hx&0x000fffff)|LO(x))!=0)
	return x+x;                     /* Nan */ 
      else return ((hx&0x80000000)==0)? x:0.0;      /* exp(+/-inf) = inf,0 */
    }
    if (x > o_bound) return huge*huge;  /* overflow  */ 
    if (x < u_bound) return small*small;/* underflow */ 
  }

  if (hx < 0x3C900000) return 1.; 

  /* Arrondi au plus près */
  DOUBLE2INT(k, x * inv_ln2);

  if (k != 0){ 
    /* r_hi+r_lo =  x - (ln2_hi + ln2_me + ln2_lo)*k */
    rp_hi =  x-ln2_hi*k;
    rp_lo = -ln2_me*k; 
    Fast2SumCond(r_hi, u, rp_hi, rp_lo);
    r_lo = u - ln2_lo*k;                         
  }else {
    r_hi = x;  r_lo = 0.;                        
  }

  /*
   * 2) Deuxième reduction d'argument
   */
  

  /* Arrondi au plus près */
  tmp    = (r_hi + two_44_43);         
  index  = LO(tmp); 
  index += bias;
  r_hi  -= (tmp - two_44_43);
 
  /* Normalisation du résultat */
  Fast2Sum(rp_hi, rp_lo, r_hi, r_lo); 

  /* Lecture de table */
  ex_hi = tab_exp[index][0];         
  ex_lo = tab_exp[index][1];         

  
  /*
   * 3) Evaluation Polynomiale
   */

  P_r = (c_0 + rp_hi * (c_1 + rp_hi * (c_2 + (rp_hi * c_3)))); 


  R1 = rp_hi * rp_hi;              

  crp_hi = R1 * rp_hi;             
  /* Correspond à R1 /= 2; */
  HI(R1) = HI(R1)-0x00100000;      

  R2 =  P_r * crp_hi;              

  /*  (ex_hi+ex_lo)*(rp_hi+rp_lo) */
  Dekker(&R3_hi, &R3_lo, ex_hi, rp_hi);      
  R4 = ex_hi * rp_lo;                      

  Dekker(&R5_hi, &R5_lo, ex_hi, R1);         
  R6 = R4 + (ex_lo * (R1 + rp_hi));      

  R7  = ex_hi * R2;                      
  R7 += (R6 + R5_lo) + (R3_lo + ex_lo);  

  Fast2Sum(R9, R8, R7, R5_hi);           
  
  Fast2Sum(R10, tmp, R3_hi, R9);           
  R8 += tmp;                               

  Fast2Sum(R11, tmp, ex_hi, R10);          
  R8 += tmp;                               

  Fast2Sum(R11, R8, R11, R8);              

  exp_R11 = (HI(R11) & 0x7ff00000) - errd;
  ulp_R11 = R11 * twom53;

  if ((HI(R8) & 0x7ff00000) > exp_R11){
    /* On est capable d'arrondir */
    if (HI(R8) < 0){     
      /* arrondi vers - inf */
      HI(ulp_R11) &= 0x7ff00000;
      R11 -= ulp_R11;
    }     
    if (k > -1020){                    
      if (k < 1020){                   
	HI(R11) += (k<<20);             
	return R11;
      }else {
	/* On est proche de + Inf */
	HI(R11) += ((k-1000)<<20);      
	return R11*two1000;
      }
    }else {
      /* On est dans les dénormalisés */
      HI(R11) += ((k+1000)<<20);        
      return R11*twom1000;
    }
  }else {
    /* Cas difficile */
    su_exp(x);
  }
}












/***************************
 ***************************
 ***  ROUDING TO ZERO    ***
 ***************************
 ***************************/

double exp_rz(double x){
  /* 1) Variables temporaires */
  double r_hi, r_lo;
  double u, tmp;
  unsigned int hx;
  int k;


  /* 2) Variables temporaires */
  double rp_hi, rp_lo, ex_hi, ex_lo;
  int index;

  /* 3) Variables temporaires */
  double P_r;

  /* 4) Variables temporaires */
  double R1, R2, R3_hi, R3_lo, R4, R5_hi, R5_lo, R6, R7, R8, R9, R10, R11, crp_hi;
  
  double ulp_R11;
  int    exp_R11;

  /*
   * 1) Première réduction d'argument
   */
  hx  = HI(x);       
  hx &= 0x7fffffff;  

  /* Filtre les cas spéciaux */
  if (hx >= 0x40862E42){
    if (hx >= 0x7ff00000){
      if (((hx&0x000fffff)|LO(x))!=0)
	return x+x;                     /* Nan */ 
      else return ((hx&0x80000000)==0)? x:0.0;      /* exp(+/-inf) = inf,0 */
    }
    if (x > o_bound) return huge*huge;  /* overflow  */ 
    if (x < u_bound) return small*small;/* underflow */ 
  }

  if (hx < 0x3C900000) return 1.; 

  /* Arrondi au plus près */
  DOUBLE2INT(k, x * inv_ln2);

  if (k != 0){ 
    /* r_hi+r_lo =  x - (ln2_hi + ln2_me + ln2_lo)*k */
    rp_hi =  x-ln2_hi*k;
    rp_lo = -ln2_me*k; 
    Fast2SumCond(r_hi, u, rp_hi, rp_lo);
    r_lo = u - ln2_lo*k;                         
  }else {
    r_hi = x;  r_lo = 0.;                        
  }

  /*
   * 2) Deuxième reduction d'argument
   */
  

  /* Arrondi au plus près */
  tmp    = (r_hi + two_44_43);         
  index  = LO(tmp); 
  index += bias;
  r_hi  -= (tmp - two_44_43);
 
  /* Normalisation du résultat */
  Fast2Sum(rp_hi, rp_lo, r_hi, r_lo); 

  /* Lecture de table */
  ex_hi = tab_exp[index][0];         
  ex_lo = tab_exp[index][1];         

  
  /*
   * 3) Evaluation Polynomiale
   */

  P_r = (c_0 + rp_hi * (c_1 + rp_hi * (c_2 + (rp_hi * c_3)))); 


  R1 = rp_hi * rp_hi;              

  crp_hi = R1 * rp_hi;             
  /* Correspond à R1 /= 2; */
  HI(R1) = HI(R1)-0x00100000;      

  R2 =  P_r * crp_hi;              

  /*  (ex_hi+ex_lo)*(rp_hi+rp_lo) */
  Dekker(&R3_hi, &R3_lo, ex_hi, rp_hi);      
  R4 = ex_hi * rp_lo;                      

  Dekker(&R5_hi, &R5_lo, ex_hi, R1);         
  R6 = R4 + (ex_lo * (R1 + rp_hi));      

  R7  = ex_hi * R2;                      
  R7 += (R6 + R5_lo) + (R3_lo + ex_lo);  

  Fast2Sum(R9, R8, R7, R5_hi);           
  
  Fast2Sum(R10, tmp, R3_hi, R9);           
  R8 += tmp;                               

  Fast2Sum(R11, tmp, ex_hi, R10);          
  R8 += tmp;                               

  Fast2Sum(R11, R8, R11, R8);              

  exp_R11 = (HI(R11) & 0x7ff00000) - errd;
  ulp_R11 = R11 * twom53;

  if ((HI(R8) & 0x7ff00000) > exp_R11){
    /* On est capable d'arrondir */
    if (HI(R8) < 0){     
      /* arrondi vers - inf */
      HI(ulp_R11) &= 0x7ff00000;
      R11 -= ulp_R11;
    }     
    if (k > -1020){                    
      if (k < 1020){                   
	HI(R11) += (k<<20);             
	return R11;
      }else {
	/* On est proche de + Inf */
	HI(R11) += ((k-1000)<<20);      
	return R11*two1000;
      }
    }else {
      /* On est dans les dénormalisés */
      HI(R11) += ((k+1000)<<20);        
      return R11*twom1000;
    }
  }else {
    /* Cas difficile */
    su_exp(x);
  }
}











main(int argc, char *argv[]){
  mpfr_t xx, yy, zz;
  double x, y, z;
  int i, e, n, k;

  if (argc != 2) exit;

  __setfpucw((_FPU_DEFAULT & ~_FPU_EXTENDED)|_FPU_DOUBLE );

  n = atoi(argv[1]);

  mpfr_init2(xx, 300); mpfr_init2(yy, 300); mpfr_init2(zz, 300);
  mpfr_set_ui(xx, 0, GMP_RNDN);
  
  k = 0;
  for(i=0; i<n; i++){
    mpfr_random(xx);
    e = (rand()*(1600./RAND_MAX))-800;
    mpfr_mul_2si(xx, xx, e, GMP_RNDN);

    x = mpfr_get_d(xx, GMP_RNDN);
    z = exp(x);

    mpfr_set_d(xx, x, GMP_RNDN);
    mpfr_exp(yy, xx, GMP_RNDN);
    y = mpfr_get_d(yy, GMP_RNDN);

    if (z != y){
      e = (z-y)/y;
      e = HI(e);
      k ++;

      printf("%.80e \n",x);
      
      //      printf("x:%e nous:%.60e mpfr:%.60e \n", x, z, y);
      printf("%d \n", ((e & 0x7fffffff )>>20)-1023);
      printf("%d\n\n",k);
    }
  } 

  mpfr_clear(xx); mpfr_clear(yy); mpfr_clear(zz);
}
