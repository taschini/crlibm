/*
 * Function to compute the hyperbolic cosine and the hyperbolic sine  with fully exact rounding
 *
 * Author : Matthieu Gallet
 * E-Mail : mgallet@ens-lyon.fr
 * Date of creation : 16/06/2004   
 */

#include <stdio.h>
#include <stdlib.h>
#include <crlibm.h>
#include <crlibm_private.h>
#include "csh_fast.h"
#include "exp.h"
/* switches on various printfs. Default 0 */
#define DEBUG 0

void scs_div_2(scs_t num) {
  /* small function to divide by 2 any SCS number */
  unsigned int carry, new_value, mask, old_value;
  int i;
  carry = 0x00000000;
  mask = ((0x1) << SCS_NB_BITS)-1;/*we now have a mask for the used bits in a word*/
  /* if it's a normal number, i.e. not zero nor NaN */
  if( (num->exception).d == (double) 1) {
    /* first, a loop to rotate all numbers to the right*/
    for(i = 0; i < SCS_NB_WORDS; i++) {
      old_value = num->h_word[i];
      new_value = (old_value & mask);/* to keep only used bits */
      num->h_word[i] = (old_value & !mask) | carry | ((old_value >> 1) & mask);
      carry = old_value & 0x00000001;/* it can be interesting to keep the last bit of each number =) */
      carry = carry << (SCS_NB_BITS-1);
    }
    /* In the SCS format, the first number can't be zero, so we must handle this particular case */
    if (num->h_word[0] == 0){
      num->index = num->index - 1;
      for(i = 1; i < SCS_NB_WORDS; i++) {
	num->h_word[i-1] = num->h_word[i];
      }
      num->h_word[SCS_NB_WORDS-1] = 0;
    }
  }
  else {
    (num->exception).d = (num->exception).d / 2;/* zero, NaN, ... */
  }
}


/*************************************************************
 *************************************************************
 *               ROUNDED  TO NEAREST			     *
 *************************************************************
 *************************************************************/
 enum {RN, RZ, RD, RU};
double cosh_quick(double x, int rounding_mode){

  /*some variable declarations */
  int k;
  db_number y;
  double res_hi, res_lo;
  double ch_hi, ch_lo, sh_hi, sh_lo;/* cosh(x) = (ch_hi + ch_lo)*(cosh(k*ln(2)) + (sh_hi + sh_lo)*(sinh(k*ln(2))) */
  db_number  table_index_float;
  int table_index;
  double temp_hi, temp_lo, temp;/* some temporary variables */
  double b_hi, b_lo,b_ca_hi, b_ca_lo, b_sa_hi, b_sa_lo;
  double ca_hi, ca_lo, sa_hi, sa_lo; /*will be the tabulated values */
  double tcb_hi, tsb_hi; /*results of polynomial approximations*/
  double square_y_hi;
  double ch_2_pk_hi, ch_2_pk_lo, ch_2_mk_hi, ch_2_mk_lo;
  double sh_2_pk_hi, sh_2_pk_lo, sh_2_mk_hi, sh_2_mk_lo;
  db_number two_p_plus_k, two_p_minus_k; /* 2^(k-1) + 2^(-k-1) */
  db_number absyh, absyl, u53, u;


  /* Now we can do the first range reduction*/
  DOUBLE2INT(k, x * inv_ln_2.d)
    if (k != 0){ /* b_hi+b_lo =  x - (ln2_hi + ln2_lo) * k */
      temp_hi = x - ln2_hi.d * k;                                         
      temp_lo = -ln2_lo.d * k;                                          
      Add12Cond(b_hi, b_lo, temp_hi, temp_lo); 
    }
    else {                                                         
      b_hi = x;  b_lo = 0.;
      if ((x < two_minus_30.d)&&(x > -two_minus_30.d)) {
	return((double) 1);
	  }
    }                                                               
  /*we'll construct 2 constants for the last reconstruction */
  two_p_plus_k.i[LO_ENDIAN] = 0;
  two_p_plus_k.i[HI_ENDIAN] = (k-1+1023) << 20;
  two_p_minus_k.i[LO_ENDIAN] = 0;
  two_p_minus_k.i[HI_ENDIAN] = (-k-1+1023) << 20;

  /* at this stage, we've done the first range reduction : we have b_hi + b_lo  between -ln(2)/2 and ln(2)/2 */
  /* now we can do the second range reduction */
  /* we'll get the 8 leading bits of b_hi */
  table_index_float.d = b_hi + two_43_44.d;
  /*this add do the float equivalent of a rotation to the right, since -0.5 <= b_hi <= 0.5*/
  table_index = LO(table_index_float.d);/* -89 <= table_index <= 89 */
  table_index_float.d -= two_43_44.d;
  table_index += bias; /* to have only positive values */
  b_hi -= table_index_float.d;/* to remove the 8 leading bits*/
  /* since b_hi was between -2^-1 and 2^1, we now have b_hi between -2^-9 and 2^-9 */


  
  y.d = b_hi;
  /*   first, y�  */
  square_y_hi = b_hi * b_hi;
  /* effective computation of the polynomial approximation */
  
  if (((y.i[HI_ENDIAN])&(0x7FFFFFFF)) < (two_minus_30.i[HI_ENDIAN])) {
    tcb_hi = 0;
    tsb_hi = 0;
  }
  else {
    /*   second, cosh(y) = y� * (1/2 + y� * (1/24 + y� * 1/720)) */
    /*tcb_hi = (square_y_hi)* (c2.d + square_y_hi * (c4.d + square_y_hi * (c6.d + square_y_hi * c8.d)));*/
    tcb_hi = (square_y_hi)* (c2.d + square_y_hi * (c4.d + square_y_hi * c6.d));
    tsb_hi = square_y_hi * (s3.d + square_y_hi * (s5.d + square_y_hi * s7.d));
  }
 

  if( table_index != bias) {
    /* we get the tabulated the tabulated values */
    ca_hi = cosh_sinh_table[table_index][0].d;
    ca_lo = cosh_sinh_table[table_index][1].d;
    sa_hi = cosh_sinh_table[table_index][2].d;
    sa_lo = cosh_sinh_table[table_index][3].d;
    
    /* first reconstruction of the cosh (corresponding to the second range reduction) */
    Mul12(&b_sa_hi,&b_sa_lo, sa_hi, b_hi);
    temp =  ((((((ca_lo + (b_hi * sa_lo)) + b_lo * sa_hi) + b_sa_lo) + (b_sa_hi * tsb_hi)) + ca_hi * tcb_hi) + b_sa_hi);
    Add12Cond(ch_hi, ch_lo, ca_hi, temp);
      /* first reconstruction for the sinh (corresponding to the second range reduction) */
  }
  else {
    Add12Cond(ch_hi, ch_lo, (double) 1, tcb_hi);
  }
  
  
  if(k != 0) {
    if( table_index != bias) {
      /* first reconstruction for the sinh (corresponding to the second range reduction) */
      Mul12(&b_ca_hi , &b_ca_lo, ca_hi, b_hi);
      temp = (((((sa_lo + (b_lo * ca_hi)) + (b_hi * ca_lo)) + b_ca_lo) + (sa_hi*tcb_hi)) + (b_ca_hi * tsb_hi));
      Add12Cond(temp_hi, temp_lo, b_ca_hi, temp);
      Add22Cond(&sh_hi, &sh_lo, sa_hi, (double) 0, temp_hi, temp_lo);
    }
    else {
      Add12Cond(sh_hi, sh_lo, b_hi, tsb_hi * b_hi + b_lo);
    }
    if((k < 35) && (k > -35) )
      {
	ch_2_pk_hi = ch_hi * two_p_plus_k.d;
	ch_2_pk_lo = ch_lo * two_p_plus_k.d;
	ch_2_mk_hi = ch_hi * two_p_minus_k.d;
	ch_2_mk_lo = ch_lo * two_p_minus_k.d;
	sh_2_pk_hi = sh_hi * two_p_plus_k.d;
	sh_2_pk_lo = sh_lo * two_p_plus_k.d;
	sh_2_mk_hi = -1 * sh_hi * two_p_minus_k.d;
	sh_2_mk_lo = -1 * sh_lo * two_p_minus_k.d;
	
	Add22Cond(&res_hi, &res_lo, ch_2_mk_hi, ch_2_mk_lo, sh_2_mk_hi, sh_2_mk_lo);
	Add22Cond(&ch_2_mk_hi, &ch_2_mk_lo , sh_2_pk_hi, sh_2_pk_lo, res_hi, res_lo);
	Add22Cond(&res_hi, &res_lo, ch_2_pk_hi, ch_2_pk_lo, ch_2_mk_hi, ch_2_mk_lo);
      } 
    else if (k >= 35) 
      {
	ch_2_pk_hi = ch_hi * two_p_plus_k.d;
	ch_2_pk_lo = ch_lo * two_p_plus_k.d;
	sh_2_pk_hi = sh_hi * two_p_plus_k.d;
	sh_2_pk_lo = sh_lo * two_p_plus_k.d;
	Add22Cond(&res_hi, &res_lo, ch_2_pk_hi, ch_2_pk_lo, sh_2_pk_hi, sh_2_pk_lo);
      }
    else /* if (k <= -35) */
      {
	ch_2_mk_hi = ch_hi * two_p_minus_k.d;
	ch_2_mk_lo = ch_lo * two_p_minus_k.d;
	sh_2_mk_hi = -1 * sh_hi * two_p_minus_k.d;
	sh_2_mk_lo = -1 * sh_lo * two_p_minus_k.d;
	Add22Cond(&res_hi, &res_lo, ch_2_mk_hi, ch_2_mk_lo, sh_2_mk_hi, sh_2_mk_lo);
      }
  }
  else {
    res_hi = ch_hi;
    res_lo = ch_lo;
  }

  switch(rounding_mode) {
  case RN:
    {  /* Test for rounding to the nearest */
      if (res_hi == (res_hi + (res_lo * round_cst_cosh.d))) return res_hi;
      break;
    }
  case RU:
    {
      /* Rounding test to + infinity */
      absyh.d = res_hi;
      absyl.d = res_lo;
      absyh.i[HI_ENDIAN] = absyh.i[HI_ENDIAN] & 0x7fffffff;/* to get the absolute value */
      absyl.i[HI_ENDIAN] = absyl.i[HI_ENDIAN] & 0x7fffffff;/* to get the absolute value */
      /*      absyl.l = absyl.l & 0x7fffffffffffffffLL;*/
      u53.l = (absyh.l & 0x7ff0000000000000LL) +  0x0010000000000000LL;
      u.l = u53.l - 0x0350000000000000LL;
      /*      printf("coucou");*/
      double delta_cst_cosh = 1e-19;
      if(absyl.d > delta_cst_cosh * u53.d){ 
	if(res_lo > 0.)  res_hi += u.d;
	return res_hi;
      }
      break;
    }
  case RD:
    {
      /* Rounding test to - infinity (or to zero) */
      absyh.d = res_hi;
      absyl.d = res_lo;
      absyh.l = absyh.l & 0x7fffffffffffffffLL;
      absyl.l = absyl.l & 0x7fffffffffffffffLL;
      u53.l = (absyh.l & 0x7ff0000000000000LL) +  0x0010000000000000LL;
      u.l = u53.l - 0x0350000000000000LL;
      double delta_cst_cosh = 1e-19;
      if(absyl.d >  delta_cst_cosh * u53.d){ 
	if(res_lo < 0.)  res_hi -= u.d;
	return res_hi;
      }
       break;
    }
  }

  /* Now, the slow pass ! */
  scs_t res_scs, exp_scs, exp_minus_scs;
#if EVAL_PERF==1
  crlibm_second_step_taken++;
#endif
  /* we'll use the cosh(x) == (exp(x) + 1/exp(x))/2 */
  
  if ((k > -35) && (k < 35)) {
    exp_SC(exp_scs, x);
    scs_inv(exp_minus_scs, exp_scs);
    scs_add(res_scs, exp_scs, exp_minus_scs);
    scs_div_2(res_scs);
  }
  else if (k >= 35) {
    exp_SC(res_scs, x);
    scs_div_2(res_scs);
  }
  else {
    exp_SC(res_scs, -x);
    scs_div_2(res_scs);
  }


  switch(rounding_mode) {
  case RN:
    scs_get_d(&res_hi, res_scs); break;
  case RU:
    scs_get_d_pinf(&res_hi, res_scs); break;
  case RD:
    scs_get_d_minf(&res_hi, res_scs); break;
  }
  return(res_hi); 
}

double cosh_rn(double x){ 
  db_number y;
  y.d = x;
  y.i[HI_ENDIAN] = y.i[HI_ENDIAN] & 0x7FFFFFFF;     /* to get the absolute value of the input */
  if (y.d > max_input_ch.d) { /* out of range */
    y.i[LO_ENDIAN] = 0; y.i[HI_ENDIAN] = 0x7FF00000; return (y.d);
  }
  if ((y.i[HI_ENDIAN] & 0x7FF00000) >= (0x7FF00000)) {    /*particular cases : QNaN, SNaN, +- oo*/
   return (y.d);
  }
  return(cosh_quick(x, RN));
}
double cosh_rz(double x){ 
  db_number y;
  y.d = x;
  y.i[HI_ENDIAN] = y.i[HI_ENDIAN] & 0x7FFFFFFF;     /* to get the absolute value of the input */
  if ((y.i[HI_ENDIAN] & 0x7FF00000) >= (0x7FF00000)) {    /*particular cases : QNaN, SNaN, +- oo*/
    return (y.d);
  }
  if (y.d > max_input_ch.d) { /* out of range */
    y.i[LO_ENDIAN] = 0xFFFFFFFF; y.i[HI_ENDIAN] = 0x7FEFFFFF; return (y.d);
  }
  return(cosh_quick(x, RD));/* cosh is always positive, so rounding to -infinite is equal to rounding to zero */
}
double cosh_ru(double x){ 
  db_number y;
  y.d = x;
  y.i[HI_ENDIAN] = y.i[HI_ENDIAN] & 0x7FFFFFFF;     /* to get the absolute value of the input */
  if (y.d > max_input_ch.d) { /* out of range */
    y.i[LO_ENDIAN] = 0; y.i[HI_ENDIAN] = 0x7FF00000; return (y.d);
  }
  if ((y.i[HI_ENDIAN] & 0x7FF00000) >= (0x7FF00000)) {    /*particular cases : QNaN, SNaN, +- oo*/
   return (y.d);
  }
  return(cosh_quick(x, RU));
}
double cosh_rd(double x){ 
  db_number y;
  y.d = x;
  y.i[HI_ENDIAN] = y.i[HI_ENDIAN] & 0x7FFFFFFF;     /* to get the absolute value of the input */
  if ((y.i[HI_ENDIAN] & 0x7FF00000) >= (0x7FF00000)) {    /*particular cases : QNaN, SNaN, +- oo*/
    return (y.d);
  }
  if (y.d > max_input_ch.d) { /* out of range */
    y.i[LO_ENDIAN] = 0xFFFFFFFF; y.i[HI_ENDIAN] = 0x7FEFFFFF; return (y.d);
  }
  return(cosh_quick(x, RD));/* cosh is always positive, so rounding to -infinite is equal to rounding to zero */
}




















/*************************************************************
 *************************************************************
 *               ROUNDED  TO NEAREST			     *
 *************************************************************
 *************************************************************/

double sinh_quick(double x, int rounding_mode){ 

  /*some variable declarations */
  int k;
  db_number y;
  double res_hi, res_lo,temp1;
  double ch_hi, ch_lo, sh_hi, sh_lo;/* cosh(x) = (sh_hi + sh_lo)*(cosh(k*ln(2)) + (ch_hi + ch_lo)*(sinh(k*ln(2))) */
  double  table_index_float;
  int table_index;
  double ch_2_pk_hi, ch_2_pk_lo, ch_2_mk_hi, ch_2_mk_lo;
  double sh_2_pk_hi, sh_2_pk_lo, sh_2_mk_hi, sh_2_mk_lo;
  double b_hi, b_lo;
  double ca_b_hi, ca_b_lo, temp_hi, temp_lo, sa_b_hi, sa_b_lo;
  double ca_hi, ca_lo, sa_hi, sa_lo; /*tabulated values */
  double tcb_hi,  tsb_hi; /*results of polynomial approximations*/
  db_number two_p_plus_k, two_p_minus_k; /* 2^(k-1) + 2^(-k-1) */
  double square_y_hi;
  db_number absyh, absyl, u53, u;
  /* b_hi + b_lo will be the reducted argument on which we'll do all the calculus */
  
  /* first, we consider special cases (out of range, etc ...) */
  /*  y.d = x;
  int hx = y.i[HI_ENDIAN] & 0x7FFFFFFF;
  if ((hx & 0x7FF00000) >= (0x7FF00000)) {
    return (y.d);
  }
  if ((hx > max_input_sh.i[HI_ENDIAN])&&(y.i[LO_ENDIAN] > max_input_sh.i[LO_ENDIAN])) {
    y.i[LO_ENDIAN] = 0; y.i[HI_ENDIAN] = 0x7FF00000; return (y.d);
  }
  */
  /* Now we can do the first range reduction*/
  DOUBLE2INT(k, x * inv_ln_2.d)
    if (k != 0){ /* b_hi + b_lo =  x - (ln2_hi + ln2_lo) * k */
      temp_hi = x - ln2_hi.d * k;                                         
      temp_lo = -ln2_lo.d * k;                                          
      Add12Cond(b_hi, b_lo, temp_hi, temp_lo); 
    }
    else {                                                         
      b_hi = x;  b_lo = 0.;
      if ((x < two_minus_30.d)&&(x > -two_minus_30.d)) return((double) x);
    }                                                               

  /*we'll construct 2 constants for the last reconstruction */
  two_p_plus_k.i[LO_ENDIAN] = 0;
  two_p_plus_k.i[HI_ENDIAN] = (k-1+1023) << 20;
  two_p_minus_k.i[LO_ENDIAN] = 0;
  two_p_minus_k.i[HI_ENDIAN] = (-k-1+1023) << 20;

  /* at this stage, we've done the first range reduction : we have b_hi + b_lo  between -ln(2)/2 and ln(2)/2 */
  /* now we can do the second range reduction */
  /* we'll get the 8 leading bits of r_hi */
  
  table_index_float = b_hi + two_43_44.d;
  /*this add do the float equivalent of a rotation to the right, since -0.5 <= b_hi <= 0.5*/
  table_index = LO(table_index_float);/* -89 <= table_index <= 89 */
  table_index_float -= two_43_44.d;
  table_index += bias; /* to have only positive values */
  b_hi -= table_index_float;/* to remove the 8 leading bits*/
  /* since b_hi was between -2^-1 and 2^1, we now have b_hi between -2^-9 and 2^-9 */
  
  y.d = b_hi;
  /*   first, y� = square_y_hi + square_y_lo  */
  square_y_hi = b_hi * b_hi;
  /* effective computation of the polyomial approximation */
  if (((y.i[HI_ENDIAN])&(0x7FFFFFFF)) <= (two_minus_30.i[HI_ENDIAN])) {
    tsb_hi = 0;
    tcb_hi = 0;
  }
  else {
    tsb_hi = square_y_hi * (s3.d + square_y_hi * (s5.d + square_y_hi * s7.d));
    /*   second, cosh(y) = y� * (1/2 + y� * (1/24 + y� * 1/720)) */
    /*   tcb_hi = (square_y_hi)* (c2.d + square_y_hi * (c4.d + square_y_hi * (c6.d + (square_y_hi * c8.d))));*/
    tcb_hi = (square_y_hi)* (c2.d + square_y_hi * (c4.d + square_y_hi * c6.d));
  }
  
  if( table_index != bias) {
    /* we get the tabulated the tabulated values*/
    ca_hi = cosh_sinh_table[table_index][0].d;
    ca_lo = cosh_sinh_table[table_index][1].d;
    sa_hi = cosh_sinh_table[table_index][2].d;
    sa_lo = cosh_sinh_table[table_index][3].d;

    /* first reconstruction for the sinh (corresponding to the second range reduction) */
    temp1 = sa_lo;
    temp1 += b_lo * ca_hi;
    temp1 += b_hi * ca_lo;
    Mul12(&ca_b_hi, &ca_b_lo, ca_hi, b_hi);
    temp1 += ca_b_lo;
    temp1 += sa_hi * tcb_hi;
    temp1 += ca_b_hi * tsb_hi;
    Add12Cond(temp_hi, temp_lo, ca_b_hi, temp1);
    Add22Cond(&sh_hi, &sh_lo, sa_hi, (double) 0, temp_hi, temp_lo);
    /* first reconstruction of the cosh (corresponding to the second range reduction) */
    temp1 = ca_lo;
    Mul12(&sa_b_hi,&sa_b_lo, sa_hi, b_hi);
    temp1 += b_hi * sa_lo;
    temp1 += b_lo * sa_hi;
    temp1 += sa_b_lo;
    temp1 += sa_b_hi * tsb_hi;
    temp1 += ca_hi * tcb_hi;
    temp1 += sa_b_hi;
    Add12Cond(ch_hi, ch_lo, ca_hi, temp1);
  }
  else {
    Add12Cond(sh_hi, sh_lo, b_hi, tsb_hi * b_hi + b_lo);
    Add12Cond(ch_hi, ch_lo, (double) 1, tcb_hi);
  }
    
  if(k != 0) {
    if( (k < 35) && (k > -35) ) {
	ch_2_pk_hi = ch_hi * two_p_plus_k.d;
	ch_2_pk_lo = ch_lo * two_p_plus_k.d;
	ch_2_mk_hi = -1 * ch_hi * two_p_minus_k.d;
	ch_2_mk_lo = -1 * ch_lo * two_p_minus_k.d;
	sh_2_pk_hi = sh_hi * two_p_plus_k.d;
	sh_2_pk_lo = sh_lo * two_p_plus_k.d;
	sh_2_mk_hi = sh_hi * two_p_minus_k.d;
	sh_2_mk_lo = sh_lo * two_p_minus_k.d;

	Add22Cond(&res_hi, &res_lo, ch_2_mk_hi, ch_2_mk_lo, sh_2_mk_hi, sh_2_mk_lo);
	Add22Cond(&ch_2_mk_hi, &ch_2_mk_lo , sh_2_pk_hi, sh_2_pk_lo, res_hi, res_lo);
	Add22Cond(&res_hi, &res_lo, ch_2_pk_hi, ch_2_pk_lo, ch_2_mk_hi, ch_2_mk_lo);
    }
    else if (k >= 35) 
      {
	ch_2_pk_hi = ch_hi * two_p_plus_k.d;
	ch_2_pk_lo = ch_lo * two_p_plus_k.d;
	sh_2_pk_hi = sh_hi * two_p_plus_k.d;
	sh_2_pk_lo = sh_lo * two_p_plus_k.d;
	Add22Cond(&res_hi, &res_lo, ch_2_pk_hi, ch_2_pk_lo, sh_2_pk_hi, sh_2_pk_lo);
      }
    else 
      {
	ch_2_mk_hi = -1 * ch_hi * two_p_minus_k.d;
	ch_2_mk_lo = -1 * ch_lo * two_p_minus_k.d;
	sh_2_mk_hi = sh_hi * two_p_minus_k.d;
	sh_2_mk_lo = sh_lo * two_p_minus_k.d;
	Add22Cond(&res_hi, &res_lo, ch_2_mk_hi, ch_2_mk_lo, sh_2_mk_hi, sh_2_mk_lo);
      }
  }
  else {
    res_hi = sh_hi;
    res_lo = sh_lo;
  }
  /*  double roundcst = 1.0020; */
  /* Test for rounding to the nearest  */
  switch(rounding_mode) {
  case RN:
    {  /* Test for rounding to the nearest */
      if (res_hi == (res_hi + (res_lo * round_cst_cosh.d))) return res_hi;
      break;
    }
  case RU:
    {
      /* Rounding test to + infinity */
      absyh.d = res_hi;
      absyl.d = res_lo;
      absyh.i[HI_ENDIAN] = absyh.i[HI_ENDIAN] & 0x7fffffff;/* to get the absolute value */
      absyl.i[HI_ENDIAN] = absyl.i[HI_ENDIAN] & 0x7fffffff;/* to get the absolute value */
      /*      absyl.l = absyl.l & 0x7fffffffffffffffLL;*/
      u53.l = (absyh.l & 0x7ff0000000000000LL) +  0x0010000000000000LL;
      u.l = u53.l - 0x0350000000000000LL;
      /*      printf("coucou");*/
      double delta_cst_cosh = 1e-19;
      if(absyl.d > delta_cst_cosh * u53.d){ 
	if(res_lo > 0.)  res_hi += u.d;
	return res_hi;
      }
      break;
    }
  case RD:
    {
      /* Rounding test to - infinity (or to zero) */
      absyh.d = res_hi;
      absyl.d = res_lo;
      absyh.l = absyh.l & 0x7fffffffffffffffLL;
      absyl.l = absyl.l & 0x7fffffffffffffffLL;
      u53.l = (absyh.l & 0x7ff0000000000000LL) +  0x0010000000000000LL;
      u.l = u53.l - 0x0350000000000000LL;
      double delta_cst_cosh = 1e-19;
      if(absyl.d >  delta_cst_cosh * u53.d){ 
	if(res_lo < 0.)  res_hi -= u.d;
	return res_hi;
      }
       break;
    }
  }

  /*  if(res_hi == (res_hi + (res_lo * round_cst_sinh.d))) {
    return res_hi;
    }*/
  
  /* Now, the slow pass ! */
  scs_t res_scs, exp_scs, exp_minus_scs;
#if EVAL_PERF==1
  crlibm_second_step_taken++;
#endif
  /* we'll use the sinh(x) == (exp(x) - 1/exp(x))/2 */
  if ((k > -35) && (k < 35)) {
    exp_SC(exp_scs, x);
    scs_inv(exp_minus_scs, exp_scs);
    scs_sub(res_scs, exp_scs, exp_minus_scs);
    scs_div_2(res_scs);
  }
  else if (k >= 35) {
    exp_SC(res_scs, x);
    scs_div_2(res_scs);
  }
  else {
    exp_SC(res_scs, -x);
    res_scs->sign = -1;
    scs_div_2(res_scs);
  }
  switch(rounding_mode) {
  case RN:
    scs_get_d(&res_hi, res_scs); break;
  case RU:
    scs_get_d_pinf(&res_hi, res_scs); break;
  case RD:
    scs_get_d_minf(&res_hi, res_scs); break;
  }
    /*  scs_get_d(&res_hi, res_scs);*/
  return(res_hi);
}

double sinh_rn(double x){ 
  db_number y;
  y.d = x;
  y.i[HI_ENDIAN] = y.i[HI_ENDIAN] & 0x7FFFFFFF;     /* to get the absolute value of the input */
  if (y.d > max_input_ch.d) { /* out of range */
    y.d = x;
    y.i[LO_ENDIAN] = 0; y.i[HI_ENDIAN] = 0x7FF00000 | (y.i[HI_ENDIAN] & 0x80000000); return (y.d);
  }
  if ((y.i[HI_ENDIAN] & 0x7FF00000) >= (0x7FF00000)) {    /*particular cases : QNaN, SNaN, +- oo*/
   return (y.d);
  }
  return(sinh_quick(x, RN));
}
double sinh_rz(double x){ 
  db_number y;
  y.d = x;
  y.i[HI_ENDIAN] = y.i[HI_ENDIAN] & 0x7FFFFFFF;     /* to get the absolute value of the input */
  if ((y.i[HI_ENDIAN] & 0x7FF00000) >= (0x7FF00000)) {    /*particular cases : QNaN, SNaN, +- oo*/
    y.d = x;
    return (y.d);
  }
  if (y.d > max_input_ch.d) { /* out of range */
    y.d = x;
    y.i[LO_ENDIAN] = 0xFFFFFFFF; y.i[HI_ENDIAN] = 0x7FEFFFFF | (y.i[HI_ENDIAN] & 0x80000000); return (y.d);
  }
  if( x > 0) {
    return(sinh_quick(x, RD));
  }
  else {
    return(sinh_quick(x, RU));
  }
}
double sinh_ru(double x){ 
  db_number y;
  y.d = x;
  y.i[HI_ENDIAN] = y.i[HI_ENDIAN] & 0x7FFFFFFF;     /* to get the absolute value of the input */
  if ((y.i[HI_ENDIAN] & 0x7FF00000) >= (0x7FF00000)) {    /*particular cases : QNaN, SNaN, +- oo*/
    y.d = x;
   return (y.d);
  }
  if (y.d > max_input_ch.d) { /* out of range */
    if(x>0) {
      y.i[LO_ENDIAN] = 0; y.i[HI_ENDIAN] = 0x7FF00000; return (y.d);
    }
    else {
      y.i[LO_ENDIAN] = 0xFFFFFFFF; y.i[HI_ENDIAN] = 0xFFEFFFFF ; return (y.d);
    }
  }
  return(sinh_quick(x, RU));
}
double sinh_rd(double x){ 
  db_number y;
  y.d = x;
  y.i[HI_ENDIAN] = y.i[HI_ENDIAN] & 0x7FFFFFFF;     /* to get the absolute value of the input */
  if ((y.i[HI_ENDIAN] & 0x7FF00000) >= (0x7FF00000)) {    /*particular cases : QNaN, SNaN, +- oo*/
    y.d = x;
   return (y.d);
  }
  if (y.d > max_input_ch.d) { /* out of range */
    if(x>0) {
      y.i[LO_ENDIAN] = 0xFFFFFFFF; y.i[HI_ENDIAN] = 0x7FEFFFFF ; return (y.d);
    }
    else {
      y.i[LO_ENDIAN] = 0; y.i[HI_ENDIAN] = 0xFFF00000; return (y.d);
    }
  }
  return(sinh_quick(x, RD));
}



#if DEBUG
  printf("index := %d  %.8X %.8X\n", table_index-bias, table_index_float.i[HI_ENDIAN], table_index_float.i[LO_ENDIAN]);
  temp2_hi.d = tsb_hi;
  printf("tsinhb_hi := hexa2ieee([\"%.8X\",\"%.8X\"]); \n", temp2_hi.i[HI_ENDIAN],temp2_hi.i[LO_ENDIAN]);
  temp2_lo.d = 0;
  printf("tsinhb_lo := hexa2ieee([\"%.8X\",\"%.8X\"]); \n", temp2_lo.i[HI_ENDIAN],temp2_lo.i[LO_ENDIAN]);
  
  temp2_hi.d = tcb_hi;
  printf("tcoshb_hi := hexa2ieee([\"%.8X\",\"%.8X\"]); \n", temp2_hi.i[HI_ENDIAN],temp2_hi.i[LO_ENDIAN]);
  temp2_lo.d = tcb_lo;
  printf("tcoshb_lo := hexa2ieee([\"%.8X\",\"%.8X\"]); \n", temp2_lo.i[HI_ENDIAN],temp2_lo.i[LO_ENDIAN]);
  temp2_hi.d = b_hi;
  printf("\nb_hi := hexa2ieee([\"%.8X\",\"%.8X\"]); \n", temp2_hi.i[HI_ENDIAN],temp2_hi.i[LO_ENDIAN]);
  temp2_lo.d = b_lo;
  printf("b_lo := hexa2ieee([\"%.8X\",\"%.8X\"]); \n", temp2_lo.i[HI_ENDIAN],temp2_lo.i[LO_ENDIAN]);
  printf("cosh_table_hi := hexa2ieee([\"%.8X\",\"%.8X\"]);\n", cosh_table[table_index][0].i[HI_ENDIAN],cosh_table[table_index][0].i[LO_ENDIAN]);
  printf("cosh_table_lo := hexa2ieee([\"%.8X\",\"%.8X\"]);\n", cosh_table[table_index][1].i[HI_ENDIAN],cosh_table[table_index][1].i[LO_ENDIAN]);
  printf("sinh_table_hi := hexa2ieee([\"%.8X\",\"%.8X\"]);\n", sinh_table[table_index][0].i[HI_ENDIAN],sinh_table[table_index][0].i[LO_ENDIAN]);
  printf("sinh_table_lo := hexa2ieee([\"%.8X\",\"%.8X\"]);\n", sinh_table[table_index][1].i[HI_ENDIAN],sinh_table[table_index][1].i[LO_ENDIAN]);

  printf("k = %d\n", k);
  printf("ch_hi := hexa2ieee([\"%.8X\",\"%.8X\"]);\n", ch_hi.i[HI_ENDIAN],ch_hi.i[LO_ENDIAN]);
  printf("ch_lo := hexa2ieee([\"%.8X\",\"%.8X\"]);\n", ch_lo.i[HI_ENDIAN],ch_lo.i[LO_ENDIAN]);
  printf("sh_hi := hexa2ieee([\"%.8X\",\"%.8X\"]);\n", sh_hi.i[HI_ENDIAN],sh_hi.i[LO_ENDIAN]);
  printf("sh_lo := hexa2ieee([\"%.8X\",\"%.8X\"]);\n", sh_lo.i[HI_ENDIAN],sh_lo.i[LO_ENDIAN]);
  printf("cr libm     : ");
#endif 
