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

void tan(scs_ptr);
void cotan(scs_ptr);
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
 * 3) Polynomial evaluation of P(X) degree 34 (? terms to computes)
 *                  
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
 * CALCUL JUSTE ....
 * Compute the tan (or cotan) of x in double floating point arithmetic with
 * correct rounding. 
 *
 * - Input x is assumed to be bounded by ~pi/2 (~ ????) in magnitude.
 * - we consider each scs digit to store 30 bits 
 * - It computes the polynom in one time ...
 */
void scs_tan(scs_ptr x_scs){
  scs_t res_scs;
  scs_t x2;
  int i;
 
scs_square(x2, x_scs);
   
  /* x < 2^-18  => tan(x)~x+x^3/3+x^5/15+x^7/315 with accuracy 2^-143 */
  
      if(x_scs->exception.i[HI_ENDIAN] < 0x5f76b88){	/* Test if x<2^(-18) */
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
    scs_mul(res_scs, constant_poly_ptr[0], x2);

    for(i=1; i<33; i++){					/* accuracy 2^(-151) */
	scs_add(res_scs, constant_poly_ptr[i], res_scs);
	scs_mul(res_scs, res_scs, x2);
    }
  
  scs_mul(res_scs, res_scs, x_scs);
  scs_add(x_scs, x_scs, res_scs);
  
    
  return;
  }
}



/*************************************************************
 *************************************************************
 *               ROUNDED  TO NEAREST
 *************************************************************
 *************************************************************/
double scs_tan_rn(double x){  
  scs_t sc1;
  scs_t sc2;
  double resd;
  int N;
#if DEBUG
    double deb1, deb2;
#endif


  scs_set_d(sc1, x);
  N = rem_pio2_scs(sc2, sc1); 	/* x (=sc2) is in [-Pi/4,Pi/4] */ 
  N = N & 0x0000001;		/* extract the 2 last bits of  N */

    switch (N){
    case 0:
	scs_tan(sc2);
	scs_get_d(&resd, sc2);
	return resd;
   	break;
    case 1:
	scs_tan(sc2);
	scs_inv(sc2, sc2);
	scs_get_d(&resd, sc2);
	return -(resd);
	break;
    default:
	fprintf(stderr,"ERREUR: %d is not a valid value in sn_tan \n", N);
	return 0.0;
    }
}


/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  -INFINITY
 *************************************************************
 *************************************************************/
double scs_tan_rd(double x){  
  scs_t sc1, sc2;
  double resd;
  int N;
 
  scs_set_d(sc1, x);
  N = rem_pio2_scs(sc2, sc1); 	/* x is in [-Pi/4,Pi/4] */ 
  N = N & 0x0000001;		/* extract the 2 last bits of  N */

    switch (N){
      case 0:
	scs_tan(sc2);
	scs_get_d_minf(&resd, sc2);
	return resd;
   	break;
      case 1:
	scs_tan(sc2);
	scs_inv(sc2, sc2);
	scs_get_d_pinf(&resd, sc2);
	return -(resd);
	break;
    default:
	fprintf(stderr,"ERREUR: %d is not a valid value in tan_rd \n", N);
	exit(1);
    }
  return resd;
}


/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  +INFINITY
 *************************************************************
 *************************************************************/
double scs_tan_ru(double x){  
  scs_t sc1, sc2;
  double resd;
  int N;

  scs_set_d(sc1, x);
  N = rem_pio2_scs(sc2, sc1); 	/* x is in [-Pi/4,Pi/4] */ 
  N = N & 0x0000001;		/* extract the 2 last bits of  N */

    switch (N){
    case 0:
	scs_tan(sc2);
	scs_get_d_pinf(&resd, sc2);
	return resd;
   	break;
    case 1:
	scs_tan(sc2);
	scs_inv(sc2, sc2);
	scs_get_d_minf(&resd, sc2);
	return -(resd);
	break;
    default:
    fprintf(stderr,"ERREUR: %d is not a valid value in su_tan \n", N);
    exit(1);
    }
  return resd;
}


