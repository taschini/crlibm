/*
 * Function to compute the sine function on [-pi/4,pi/4]
 *
 * Author : Defour David  (David.Defour@ens-lyon.fr), Daramy Catherine (Catherine.Daramy@ens-lyon.fr)
 *
 * Date of creation : 11/03/2002   
 * Last Modified    : 10/07/2003
 */
#include <stdio.h>
#include <stdlib.h>
#include "sine.h"

void cosine(scs_ptr);
void sine(scs_ptr);



/*
 *  WHAT WE CAN DO :
 *
 * 1) Range reduction if needed ... x in [-Pi/4, +Pi/4]
 *
 * 2) call cosine or sine function
 *
 * 3) Polynomial evaluation of sin(z) or cos(z) degree 16 
 *                                  (-160)  
 *   Approximation error: |err| < 2^ 
 *
 */   

void sine(scs_ptr x){
  scs_t res_scs;
  scs_t x2;
  int i;
 
  /* x < 2^-75  => cos(x)~1 (with accuracy 2^-150,999), when rounding to nearest, here we consider x< 2^-90 */
  if(X_IND < -3){
    return; 
    }  

   
 /* Polynomial evaluation of sin(x) */

  scs_square(x2, x);
  scs_mul(res_scs, constant_poly_ptr[0], x2);

  for(i=1; i<16; i++){
    scs_add(res_scs, constant_poly_ptr[i], res_scs);
    scs_mul(res_scs, res_scs, x2);
  }
  scs_mul(res_scs, res_scs, x);
  scs_add(x, x, res_scs);

  return;
}



/************************************************************
 *************************************************************
 *               ROUNDED  TO NEAREST
 *************************************************************
 *************************************************************/
double scs_sin_rn(double x){  
  scs_t sc1, sc2;
  double resd;
  int N;


#if EVAL_PERF
	crlibm_second_step_taken++;
#endif


  scs_set_d(sc1, x);
  N = rem_pio2_scs(sc2, sc1);
  N = N & 0x0000003;		/* extract the 2 last bits of  N */
  switch (N){
  case 0:
    sine(sc2);
    scs_get_d(&resd, sc2);
    return resd;
  case 1:
    cosine(sc2);
    scs_get_d(&resd, sc2);
    return resd;
  case 2:
    sine(sc2);
    scs_get_d(&resd, sc2);
    return -resd;		
  case 3:
    cosine(sc2);
    scs_get_d(&resd, sc2);
    return -resd;
    default:
    fprintf(stderr,"ERREUR: %d is not a valid value in s_sine \n", N);
    return 0.0;
  }
}




/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  -INFINITY
 *************************************************************
 *************************************************************/
double scs_sin_rd(double x){  
  scs_t sc1, sc2;
  double resd;
  int N;
    
  scs_set_d(sc1, x);
  N = rem_pio2_scs(sc2, sc1);
  N = N & 0x0000003;		/* extract the 2 last bits of  N */
  switch (N){
  case 0:
    sine(sc2);
    scs_get_d_minf(&resd, sc2);
    return resd;
  case 1:
    cosine(sc2);
    scs_get_d_minf(&resd, sc2);
    return resd;
  case 2:  
    sine(sc2);
    scs_get_d_pinf(&resd, sc2);
    return -resd;
  case 3:
    cosine(sc2);
    scs_get_d_pinf(&resd, sc2);
    return -resd;
  default:
    fprintf(stderr,"ERREUR: %d is not a valid value in s_sine \n", N);
    exit(1);
  }
  return resd;
}



/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  +INFINITY
 *************************************************************
 *************************************************************/
double scs_sin_ru(double x){  
  scs_t sc1, sc2;
  double resd;
  int N;

  scs_set_d(sc1, x);
  N = rem_pio2_scs(sc2, sc1);
  N = N & 0x0000003;		/* extract the 2 last bits of  N */
  switch (N){
  case 0:
    sine(sc2);
    scs_get_d_pinf(&resd, sc2);
    return resd;
  case 1:
    cosine(sc2);
    scs_get_d_pinf(&resd, sc2);
    return resd;
   case 2:
    sine(sc2);
    scs_get_d_minf(&resd, sc2);
    return -resd;
  case 3:
    cosine(sc2);
    scs_get_d_minf(&resd, sc2);
    return -resd;
  default:
    fprintf(stderr,"ERREUR: %d is not a valid value in s_sine \n", N);
    exit(1);
  }
  return resd;
}


/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  ZERO
 *************************************************************
 *************************************************************/
double scs_sin_rz(double x){  
  scs_t sc1, sc2;
  double resd;
  int N;

  scs_set_d(sc1, x);
  N = rem_pio2_scs(sc2, sc1);
  N = N & 0x0000003;		/* extract the 2 last bits of  N */
  switch (N){
  case 0:
    sine(sc2);
    scs_get_d_zero(&resd, sc2);
    return resd;
  case 1:
    cosine(sc2);
    scs_get_d_zero(&resd, sc2);
    return resd;
   case 2:
    sine(sc2);
    scs_get_d_zero(&resd, sc2);
    return -resd;
  case 3:
    cosine(sc2);
    scs_get_d_zero(&resd, sc2);
    return -resd;
  default:
    fprintf(stderr,"ERREUR: %d is not a valid value in s_sine \n", N);
    exit(1);
  }
  return resd;
}



