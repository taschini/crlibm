/*
 * Function to compute the hyperbolic cosine and the hyperbolic sine  with fully exact rounding
 *
 *
 * Date of creation : 16/06/2004   
 */

#include <stdio.h>
#include <stdlib.h>
#include <crlibm.h>
#include <crlibm_private.h>
#include "csh_fast.h"
/* switches on various printfs. Default 0 */
#define DEBUG 0


/*************************************************************
**************************************************************
*         COMMON PROCEDURES TO COSH & SINH                   *
**************************************************************
**************************************************************/

void ch_reducted(double *res_hi, double * res_lo, double x_hi, double x_lo) {

  /* we'll use cosh(x) == 1 + x²/2 + x^4 * P(x) */

  double square_x_hi, square_x_lo;

  double temp1_hi, temp1_lo, temp2_hi, temp2_lo;
  double poly;
  db_number x;
  x.d = x_hi;
#if DEBUG
  printf("x = hexa2ieee([\"%X\",\"%X\"]); \n", x.i[HI_ENDIAN],x.i[LO_ENDIAN]);

#endif
  if (((x.i[HI_ENDIAN])&(0x7FFFFFFF)) <= (min_absolute_input_ch.i[HI_ENDIAN])) {
    *res_hi = 1;
    *res_lo = 0;
  }
  else {

    /*   first, x² = square_x_hi + square_x_lo  */
    Mul22(&square_x_hi, &square_x_lo, x_hi, x_lo, x_hi, x_lo);
  
    /*  poly = (double) 1;*/
    poly = square_x_hi * c7.d;
    poly = poly + c6.d;
    poly = poly * square_x_hi;
    poly = poly + c5.d;
    poly = poly * square_x_hi;
    poly = poly + c4.d;
    poly = poly * square_x_hi;
    poly = poly + c3.d;
    poly = poly * square_x_hi;
    poly = poly + c2.d;
#if 0
    poly = poly * square_x_hi;
    poly = poly + c1.d;
    Mul22(&temp1_hi, &temp1_lo, square_x_hi, square_x_lo, poly, (double) 0);
#else
    Mul22(&temp1_hi, &temp1_lo, square_x_hi, square_x_lo, poly, (double) 0);
    Add22(&temp2_hi, &temp2_lo, temp1_hi, temp1_lo, (double) c1.d, (double) 0);
    Mul22(&temp1_hi, &temp1_lo, square_x_hi, square_x_lo, temp2_hi, temp2_lo);
#endif
    Add22(res_hi, res_lo, temp1_hi, temp1_lo, c0.d, (double) 0);
#if DEBUG
    db_number debug;
    debug.d = *res_hi;
    printf("hexa2ieee([\"%X\",\"%X\"]); \n", debug.i[HI_ENDIAN],debug.i[LO_ENDIAN]);
    debug.d = *res_lo;
    printf("hexa2ieee([\"%X\",\"%X\"]); \n", debug.i[HI_ENDIAN],debug.i[LO_ENDIAN]);
#endif
  }
}

void sh_reducted(double *res_hi, double *res_lo, double x_hi, double x_lo) {
  //we'll use sinh(x) == x * (1 + x² * P(x))
  /* TODO : un test pour les dénormalisés ici |x| < min_absolute_input_sh */

  double square_x_hi, square_x_lo;
  double temp1_hi, temp1_lo, temp2_hi, temp2_lo;
  double poly;
  db_number x;
  x.d = x_hi;
#if DEBUG
  printf("x = hexa2ieee([\"%X\",\"%X\"]); \n", x.i[HI_ENDIAN],x.i[LO_ENDIAN]);

#endif

  if (((x.i[HI_ENDIAN])&(0x7FFFFFFF)) <= (min_absolute_input_sh.i[HI_ENDIAN])) {
    *res_hi = x_hi;
    *res_lo = x_lo;
  }
  else {
    /*   first, x² = square_x_hi + square_x_lo  */
    Mul22(&square_x_hi, &square_x_lo, x_hi, x_lo, x_hi, x_lo);
    /*  poly = (double) 1;*/
    
    poly = square_x_hi * s6.d;
    poly = poly * square_x_hi;
    poly = poly + s5.d;
    poly = poly * square_x_hi;
    poly = poly + s4.d;
    poly = poly * square_x_hi;
    poly = poly + s3.d;
    poly = poly * square_x_hi;
    poly = poly + s2.d;
    poly = poly * square_x_hi;
    Add22(&temp2_hi, &temp2_lo, poly, (double) 0,  s1_hi.d, s1_lo.d);
    Mul22(&temp1_hi, &temp1_lo, square_x_hi, square_x_lo, temp2_hi, temp2_lo);
    Add22(&temp2_hi, &temp2_lo, temp1_hi, temp1_lo, s0_hi.d, s0_lo.d);
    Mul22(res_hi, res_lo, temp2_hi, temp2_lo, x_hi, x_lo);
#if DEBUG
    db_number debug;
    debug.d = *res_hi;
    printf("hexa2ieee([\"%X\",\"%X\"]); \n", debug.i[HI_ENDIAN],debug.i[LO_ENDIAN]);
    debug.d = *res_lo;
    printf("hexa2ieee([\"%X\",\"%X\"]); \n", debug.i[HI_ENDIAN],debug.i[LO_ENDIAN]);
#endif
  }
}

/*************************************************************
 *************************************************************
 *               ROUNDED  TO NEAREST			     *
 *************************************************************
 *************************************************************/

double cosh_rn(double x){ 
  //some variable declarations
  int k;
  db_number y;
  double res_hi, res_lo;
  db_number ch_hi, ch_lo, sh_hi, sh_lo;
  db_number temp1_hi, temp1_lo, temp2_hi, temp2_lo, temp3_hi, temp3_lo, temp4_hi, temp4_lo;
  double r_hi, r_lo; /*, temp1_hi, temp1_lo, temp2_hi, temp2_lo;*/
  // r_hi + r_lo will be the reducted argument on which we'll do all the calculus

  //first, we consider special cases (out of range, etc ...)
  y.d = x;
  int hx = y.i[HI_ENDIAN] & 0x7FFFFFFF; /* to get the absolute value of the input */


  if ((hx & 0x7FF00000) >= (0x7FF00000)) {/*particular cases : QNaN, SNaN, +- oo*/
    return (x);
  }
  if ((hx > max_inf_input_ch.i[HI_ENDIAN])&&(y.i[LO_ENDIAN] > max_inf_input_ch.i[LO_ENDIAN])) {
    y.i[LO_ENDIAN] = 0;
    y.i[HI_ENDIAN] = 0x7FF00000;
    return(y.d);
  }
  /*the infinite is obtained earlier with the exponential than with cosh !*/
  if ((hx > max_input_ch.i[HI_ENDIAN])&&(hx < max_inf_input_exp.i[HI_ENDIAN])) {
     /* we use cosh(x) = exp(x)//  2 since x is great enough*/
    y.i[HI_ENDIAN] = hx;/* x <- |x| */
    y.d = exp_rn(y.d);
     k = (y.i[HI_ENDIAN] & 0x7FF00000) - (1 << 20);/* exp <- exp -1 */
    y.i[HI_ENDIAN] = ((y.i[HI_ENDIAN] & 0x000FFFFF)) | k;
    return(y.d);
  }

  //Now we can do the first range reduction
  DOUBLE2INT(k, x * inv_ln_2.d)
  if (k != 0){ 
    /* r_hi+r_lo =  x - (ln2_hi + ln2_lo) * k */
    temp1_hi.d =  x - ln2_hi.d * k;                                         
    temp1_lo.d =   -ln2_lo.d * k;                                          
    Add12Cond(r_hi, r_lo, temp1_hi.d, temp1_lo.d);                      
  }
  else {                                                         
    r_hi = x;  r_lo = 0.;                                        
}                                                               
  // at this stage, we've done the first range reduction : we have r_hi+r_lo  between -ln(2)/2 and ln(2)/2
 
  ch_reducted(&(ch_hi.d), &(ch_lo.d), r_hi,r_lo);
  if(k != 0) {
    sh_reducted(&(sh_hi.d), &(sh_lo.d), r_hi, r_lo);
#if DEBUG
    printf("ch_hi := hexa2ieee([\"%X\",\"%X\"]); ", ch_hi.i[HI_ENDIAN],ch_hi.i[LO_ENDIAN]);
    printf("ch_lo := hexa2ieee([\"%X\",\"%X\"]); ", ch_lo.i[HI_ENDIAN],ch_lo.i[LO_ENDIAN]);
    printf("sh_hi := hexa2ieee([\"%X\",\"%X\"]); ", sh_hi.i[HI_ENDIAN],sh_hi.i[LO_ENDIAN]);
    printf("sh_lo := hexa2ieee([\"%X\",\"%X\"]); ", sh_lo.i[HI_ENDIAN],sh_lo.i[LO_ENDIAN]);
    printf("k := %d;\n", k);
#endif
    if( (k < 35) && (k > -35) ) {
      temp1_hi.i[HI_ENDIAN] = ch_hi.i[HI_ENDIAN] + ((k-1) << 20);
      temp1_hi.i[LO_ENDIAN] = ch_hi.i[LO_ENDIAN];
      temp1_lo.i[HI_ENDIAN] = ch_lo.i[HI_ENDIAN] + ((k-1) << 20);
      temp1_lo.i[LO_ENDIAN] = ch_lo.i[LO_ENDIAN];

      temp2_hi.i[HI_ENDIAN] = ch_hi.i[HI_ENDIAN] + ((-k-1) << 20);
      temp2_hi.i[LO_ENDIAN] = ch_hi.i[LO_ENDIAN];
      temp2_lo.i[HI_ENDIAN] = ch_lo.i[HI_ENDIAN] + ((-k-1) << 20);
      temp2_lo.i[LO_ENDIAN] = ch_lo.i[LO_ENDIAN];
      
      temp3_hi.i[HI_ENDIAN] = sh_hi.i[HI_ENDIAN] + ((k-1) << 20);
      temp3_hi.i[LO_ENDIAN] = sh_hi.i[LO_ENDIAN];
      temp3_lo.i[HI_ENDIAN] = sh_lo.i[HI_ENDIAN] + ((k-1) << 20);
      temp3_lo.i[LO_ENDIAN] = sh_lo.i[LO_ENDIAN];

      temp4_hi.i[HI_ENDIAN] = sh_hi.i[HI_ENDIAN] + ((-k-1) << 20) + (1 << 31);
      temp4_hi.i[LO_ENDIAN] = sh_hi.i[LO_ENDIAN];
      temp4_lo.i[HI_ENDIAN] = sh_lo.i[HI_ENDIAN] + ((-k-1) << 20) + (1 << 31);
      temp4_lo.i[LO_ENDIAN] = sh_lo.i[LO_ENDIAN];

      Add22(&res_hi, &res_lo, temp1_hi.d, temp1_lo.d, temp3_hi.d, temp3_lo.d);
      Add22(&temp1_hi.d, &temp1_lo.d, res_hi, res_lo, temp2_hi.d, temp2_lo.d);
      Add22(&res_hi, &res_lo, temp1_hi.d, temp1_lo.d, temp4_hi.d, temp4_lo.d);
    }
    else if (k >= 35) 
      {
      temp1_hi.i[HI_ENDIAN] = ch_hi.i[HI_ENDIAN] + ((k-1) << 20);
      temp1_hi.i[LO_ENDIAN] = ch_hi.i[LO_ENDIAN];
      temp1_lo.i[HI_ENDIAN] = ch_lo.i[HI_ENDIAN] + ((k-1) << 20);
      temp1_lo.i[LO_ENDIAN] = ch_lo.i[LO_ENDIAN];

      temp3_hi.i[HI_ENDIAN] = sh_hi.i[HI_ENDIAN] + ((k-1) << 20);
      temp3_hi.i[LO_ENDIAN] = sh_hi.i[LO_ENDIAN];
      temp3_lo.i[HI_ENDIAN] = sh_lo.i[HI_ENDIAN] + ((k-1) << 20);
      temp3_lo.i[LO_ENDIAN] = sh_lo.i[LO_ENDIAN];
      
      
      Add22(&res_hi, &res_lo, temp1_hi.d, temp1_lo.d, temp3_hi.d, temp3_lo.d);
      }
    else /* if (k <= -35) */
      {/* sinh(x + k.ln(2)) = cosh(x)*(-2^(k-1)) + (2^(k-1))*sinh(x) */
	temp1_hi.i[HI_ENDIAN] = ch_hi.i[HI_ENDIAN] + ((-k-1) << 20);
	temp1_hi.i[LO_ENDIAN] = ch_hi.i[LO_ENDIAN];
	temp1_lo.i[HI_ENDIAN] = ch_lo.i[HI_ENDIAN] + ((-k-1) << 20);
	temp1_lo.i[LO_ENDIAN] = ch_lo.i[LO_ENDIAN];
	
	temp3_hi.i[HI_ENDIAN] = sh_hi.i[HI_ENDIAN] + ((-k-1) << 20) + (1 << 31);
	temp3_hi.i[LO_ENDIAN] = sh_hi.i[LO_ENDIAN];
	temp3_lo.i[HI_ENDIAN] = sh_lo.i[HI_ENDIAN] + ((-k-1) << 20) + (1 << 31);
	temp3_lo.i[LO_ENDIAN] = sh_lo.i[LO_ENDIAN];
	
	Add22(&res_hi, &res_lo, temp1_hi.d, temp1_lo.d, temp3_hi.d, temp3_lo.d);
      }

  }
  else {
    res_hi = ch_hi.d;
    res_lo = ch_lo.d;
    }
  return(res_hi);
   /*

   // Call the actual computation
   log_quick(&res_hi, &res_lo, &rndcstindex, &y, E);
   roundcst = rncst[rndcstindex];


  // Test for rounding to the nearest 
  if(res_hi == (res_hi + (res_lo * roundcst)))
    return res_hi;
  else {
    scs_t res;
#if DEBUG
    printf("Going for Accurate Phase for x=%1.50e\n",x);
#endif
    scs_log(res, y, E);
    scs_get_d(&res_hi, res);
    return res_hi;
  }*/
 }

/*************************************************************
 *************************************************************
 *               ROUNDED  TO NEAREST			     *
 *************************************************************
 *************************************************************/

double sinh_rn(double x){ 
  //some variable declarations
  int k;
  db_number y;
  double res_hi, res_lo;/*, ch_hi, ch_lo, sh_hi, sh_lo;*/
  db_number ch_hi, ch_lo, sh_hi, sh_lo;
  db_number temp1_hi, temp1_lo, temp2_hi, temp2_lo, temp3_hi, temp3_lo, temp4_hi, temp4_lo;
  double r_hi, r_lo;/*, temp1_hi, temp1_lo, temp2_hi, temp2_lo;*/
  /* r_hi + r_lo will be the reducted argument on which we'll do all the computation */

  /*first, we consider special cases (out of range, etc ...)*/
  y.d = x;

  int hx = y.i[HI_ENDIAN] & 0x7FFFFFFF; /* to get the absolute value of the input */
  int sx = y.i[HI_ENDIAN] - hx;/* to get the sign */

  if ((hx & 0x7FF00000) >= (0x7FF00000)) {/*particular cases : QNaN, SNaN, +- oo*/
    return (x);
  }
  if ((hx > max_inf_input_sh.i[HI_ENDIAN])&&(y.i[LO_ENDIAN] > max_inf_input_sh.i[LO_ENDIAN])) {
    y.i[LO_ENDIAN] = 0;
    y.i[HI_ENDIAN] = 0x7FF00000 + sx;
    return(y.d);
  }
  /*the infinite is obtained earlier with the exponential than with cosh !*/
  if ((hx > max_input_sh.i[HI_ENDIAN])&&(hx < max_inf_input_exp.i[HI_ENDIAN])) {
     /* we use cosh(x) = exp(x)//  2 since x is great enough*/
    y.i[HI_ENDIAN] = hx;/* x <- |x| */
    y.d = exp_rn(y.d);
     k = (y.i[HI_ENDIAN] & 0x7FF00000) - (1 << 20);/* exp <- exp -1 */
    y.i[HI_ENDIAN] = ((y.i[HI_ENDIAN] & 0x000FFFFF) + sx) | k;
    return(y.d);
  }
#if DEBUG
  printf("On passe par chez moi ^^ \n ");
#endif
  //Now we can do the first range reduction
  DOUBLE2INT(k, x * inv_ln_2.d)
     
     if (k != 0){ 
    // r_hi+r_lo =  x - (ln2_hi + ln2_lo) * k 
    temp1_hi.d =  x - ln2_hi.d * k;                                         
    temp1_lo.d =   -ln2_lo.d * k;                                          
    Add12Cond(r_hi, r_lo, temp1_hi.d, temp1_lo.d);
  }
  else {                                                         
    r_hi = x;  r_lo = 0.;                                        
}                                                               
  // at this stage, we've done the first range reduction : we have r_hi+r_lo  between -ln(2)/2 and ln(2)/2
  sh_reducted(&sh_hi.d, &sh_lo.d, r_hi,r_lo);
  if(k != 0) {
        ch_reducted(&ch_hi.d, &ch_lo.d, r_hi, r_lo);
#if DEBUG
    printf("ch_hi := hexa2ieee([\"%X\",\"%X\"]); ", ch_hi.i[HI_ENDIAN],ch_hi.i[LO_ENDIAN]);
    printf("ch_lo := hexa2ieee([\"%X\",\"%X\"]); ", ch_lo.i[HI_ENDIAN],ch_lo.i[LO_ENDIAN]);
    printf("sh_hi := hexa2ieee([\"%X\",\"%X\"]); ", sh_hi.i[HI_ENDIAN],sh_hi.i[LO_ENDIAN]);
    printf("sh_lo := hexa2ieee([\"%X\",\"%X\"]); ", sh_lo.i[HI_ENDIAN],sh_lo.i[LO_ENDIAN]);
    printf("k := %d;\n", k);
#endif
    if( (k < 35) && (k > -35) ) {
      temp1_hi.i[HI_ENDIAN] = ch_hi.i[HI_ENDIAN] + ((k-1) << 20);
      temp1_hi.i[LO_ENDIAN] = ch_hi.i[LO_ENDIAN];
      temp1_lo.i[HI_ENDIAN] = ch_lo.i[HI_ENDIAN] + ((k-1) << 20);
      temp1_lo.i[LO_ENDIAN] = ch_lo.i[LO_ENDIAN];

      temp2_hi.i[HI_ENDIAN] = ch_hi.i[HI_ENDIAN] + ((-k-1) << 20) + (1 << 31);
      temp2_hi.i[LO_ENDIAN] = ch_hi.i[LO_ENDIAN];
      temp2_lo.i[HI_ENDIAN] = ch_lo.i[HI_ENDIAN] + ((-k-1) << 20) + (1 << 31);
      temp2_lo.i[LO_ENDIAN] = ch_lo.i[LO_ENDIAN];
      
      temp3_hi.i[HI_ENDIAN] = sh_hi.i[HI_ENDIAN] + ((k-1) << 20);
      temp3_hi.i[LO_ENDIAN] = sh_hi.i[LO_ENDIAN];
      temp3_lo.i[HI_ENDIAN] = sh_lo.i[HI_ENDIAN] + ((k-1) << 20);
      temp3_lo.i[LO_ENDIAN] = sh_lo.i[LO_ENDIAN];

      temp4_hi.i[HI_ENDIAN] = sh_hi.i[HI_ENDIAN] + ((-k-1) << 20);
      temp4_hi.i[LO_ENDIAN] = sh_hi.i[LO_ENDIAN];
      temp4_lo.i[HI_ENDIAN] = sh_lo.i[HI_ENDIAN] + ((-k-1) << 20);
      temp4_lo.i[LO_ENDIAN] = sh_lo.i[LO_ENDIAN];

      Add22(&res_hi, &res_lo, temp1_hi.d, temp1_lo.d, temp3_hi.d, temp3_lo.d);
      Add22(&temp1_hi.d, &temp1_lo.d, res_hi, res_lo, temp2_hi.d, temp2_lo.d);
      Add22(&res_hi, &res_lo, temp1_hi.d, temp1_lo.d, temp4_hi.d, temp4_lo.d);
    }
    else if (k >= 35) 
      {
      temp1_hi.i[HI_ENDIAN] = ch_hi.i[HI_ENDIAN] + ((k-1) << 20);
      temp1_hi.i[LO_ENDIAN] = ch_hi.i[LO_ENDIAN];
      temp1_lo.i[HI_ENDIAN] = ch_lo.i[HI_ENDIAN] + ((k-1) << 20);
      temp1_lo.i[LO_ENDIAN] = ch_lo.i[LO_ENDIAN];
      
      temp3_hi.i[HI_ENDIAN] = sh_hi.i[HI_ENDIAN] + ((k-1) << 20);
      temp3_hi.i[LO_ENDIAN] = sh_hi.i[LO_ENDIAN];
      temp3_lo.i[HI_ENDIAN] = sh_lo.i[HI_ENDIAN] + ((k-1) << 20);
      temp3_lo.i[LO_ENDIAN] = sh_lo.i[LO_ENDIAN];
      
      Add22(&res_hi, &res_lo, temp1_hi.d, temp1_lo.d, temp3_hi.d, temp3_lo.d);
      }
    else 
      {
	temp1_hi.i[HI_ENDIAN] = ch_hi.i[HI_ENDIAN] + ((-k-1) << 20) + (1 << 31);
	temp1_hi.i[LO_ENDIAN] = ch_hi.i[LO_ENDIAN];
	temp1_lo.i[HI_ENDIAN] = ch_lo.i[HI_ENDIAN] + ((-k-1) << 20) + (1 << 31);
	temp1_lo.i[LO_ENDIAN] = ch_lo.i[LO_ENDIAN];
	
	temp3_hi.i[HI_ENDIAN] = sh_hi.i[HI_ENDIAN] + ((-k-1) << 20);
	temp3_hi.i[LO_ENDIAN] = sh_hi.i[LO_ENDIAN];
	temp3_lo.i[HI_ENDIAN] = sh_lo.i[HI_ENDIAN] + ((-k-1) << 20);
	temp3_lo.i[LO_ENDIAN] = sh_lo.i[LO_ENDIAN];
	
	Add22(&res_hi, &res_lo, temp1_hi.d, temp1_lo.d, temp3_hi.d, temp3_lo.d);
      }
  }
    else {
  res_hi = sh_hi.d;
  res_lo = sh_lo.d;
    }
    db_number debug;
    debug.d = res_hi;
#if DEBUG
   printf("hexa2ieee([\"%X\",\"%X\"]); \n", debug.i[HI_ENDIAN],debug.i[LO_ENDIAN]);
#endif
  return(res_hi);

 }


