/**
 * Function to compute the tan and the Cotan functions on [-pi/2,pi/2]
 *
 * Author : Daramy Catherine  (Catherine.Daramy@ens-lyon.fr)
 *
 * Date of creation : 11/03/2002   
 * Last Modified    : 15/03/2003
 */
#include <stdio.h>
#include <stdlib.h>
#include "tan.h"
#include <coefpi2.h>

int rem_pio2_scs(scs_ptr, scs_ptr);

/**
 *  WHAT WE CAN DO :
 *
 * We are computing tan and cotan as the same time, assuming that tan (x + k*Pi/2) = Cotan (x) 
 *	 
 *
 * 1) Range reduction if needed ... find x in [-Pi/2, +Pi/2] tan is Pi-periodic
 *
 * 2) Third range reduction ... find x in [ 0, Pi/2] tan(x) = -tan(-x)
 *                                
 * 3) Polynomial evaluation of P(X) degree 34
 *                                  (-149)  
 *   Approximation error: |err| < 2^ 
 *
 *
 * 4) "Reconstruction":
 *
 *   tan(x) = P(X)
 *
 */

/*
 * Compute the tan (or cotan) of x in double floating point arithmetic with
 * correct rounding. 
 *
 * - Input x is assumed to be bounded by ~pi/2 in magnitude.
 * - we consider each scs digit to store 30 bits 
 * - It computes the polynom in one time ...
 */
void scs_tan(double x, scs_ptr res_scs){
  scs_t x_scs;
  scs_t x2;
  int i;
   scs_t y_scs;
  int N;

  scs_set_d(x_scs, x);
  
  /* x < 2^-18  => tan(x)~x+x^3/3+x^5/15+x^7/315 with accuracy 2^-143 */
  
  if(x_scs->exception.i[HI_ENDIAN] < 0x5f76b88){	/* Test if x<2^(-18) */
    scs_square(x2, x_scs);
    scs_mul(res_scs, cste_poly_ptr[0], x2);
    scs_add(res_scs, cste_poly_ptr[1], res_scs);
    scs_mul(res_scs, res_scs, x2);
    scs_add(res_scs, cste_poly_ptr[2], res_scs);
    scs_mul(res_scs, res_scs, x2);
    scs_mul(res_scs, x_scs, res_scs);
    scs_add(x_scs, x_scs, res_scs);	
    return;
  }
    
  /* Polynomial evaluation of tan(x) */
  
  else {
    N = rem_pio2_scs(y_scs, x_scs); 	/* x (=sc2) is in [-Pi/4,Pi/4] */ 
    N = N & 1;		/* extract the last bit of  N */
    scs_square(x2, y_scs);

    scs_mul(res_scs, constant_poly_ptr[0], x2);
    
    for(i=1; i<33; i++){					/* accuracy 2^(-151) */
      scs_add(res_scs, constant_poly_ptr[i], res_scs);
      scs_mul(res_scs, res_scs, x2);
    }
    
    scs_mul(res_scs, res_scs, y_scs);
    scs_add(res_scs, y_scs, res_scs);
    
    if(N==1) {
      scs_inv(res_scs, res_scs);
      res_scs->sign = -res_scs->sign;
    }
  }
}

/*************************************************************
 *************************************************************
 *               ROUNDED  TO NEAREST
 *************************************************************
 *************************************************************/
double scs_tan_rn(double x){  
  scs_t res_scs;
  double resd;

#if EVAL_PERF
  crlibm_second_step_taken++;
#endif

  scs_tan(x,res_scs);
  scs_get_d(&resd, res_scs);
  return resd;
}


/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  -INFINITY
 *************************************************************
 *************************************************************/
double scs_tan_rd(double x){  
  scs_t res_scs;
  double resd;

#if EVAL_PERF
  crlibm_second_step_taken++;
#endif

  scs_tan(x,res_scs);
  scs_get_d_minf(&resd, res_scs);
  return resd;
}


/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  +INFINITY
 *************************************************************
 *************************************************************/
double scs_tan_ru(double x){  
  scs_t res_scs;
  double resd;

#if EVAL_PERF
  crlibm_second_step_taken++;
#endif

  scs_tan(x,res_scs);
  scs_get_d_pinf(&resd, res_scs);
  return resd;
}
/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  ZERO
 *************************************************************
 *************************************************************/
double scs_tan_rz(double x){  
  scs_t res_scs;
  double resd;

#if EVAL_PERF
  crlibm_second_step_taken++;
#endif

  scs_tan(x,res_scs);
  scs_get_d_zero(&resd, res_scs);
  return resd;
}


