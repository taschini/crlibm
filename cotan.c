/*
 * Function to compute the tan and the Cotan functions on [-pi/2,pi/2]
 *
 * Author : Daramy Catherine  (Catherine.Daramy@ens-lyon.fr)
 *
 * Date of creation : 11/03/2002   
 * Last Modified    : 15/03/2002
 */
#include <stdio.h>
#include <stdlib.h>
#include "coefpi2.h"
#include "tan.h"

void tan(scs_ptr x);
int rem_pio2_scs (scs_ptr, scs_ptr);

/*
 *  WHAT WE CAN DO :
 *
 * We are computing tan and cotan as the same time, assuming that cotan (x) = tan (PI/2 - x) 
 *	 
 *
 * 1) Range reduction if needed ... find x in [-Pi/2, +Pi/2] tan is Pi-periodic
 *
 * 2) Third range reduction ... find x in [ 0, Pi/2] tan(x) = -tan(-x)
 *                                
 * 3) Polynomial evaluation of P(X) degree ??? (? terms to computes)
 *                  
 *                                  (-150)  
 *   Approximation error: |err| < 2^ 
 *
 *
 * 4) "Reconstruction":
 *
 *   tan(x) = P(X)
 *
 */


/*************************************************************
 *************************************************************
 *               ROUNDED  TO NEAREST
 *************************************************************
 *************************************************************/
double cotan_rn(double x){
  scs_t sc1, sc2;
  double resd;
  int N;
  
  scs_set_d(sc1, x);
  N = rem_pio2_scs(sc2, sc1); 	/* x is in [0,Pi/4] */ 
  
    switch (N){
    case 0:
	tan(sc2);
	scs_inv(sc2, sc2);
	scs_get_d(&resd, sc2);
	return resd;
   	break;
    case 1:
	tan(sc2);
	scs_get_d(&resd, sc2);
	return -(resd);
	break;
    case 2:
	scs_sub(sc2, Pio2_ptr, sc2);
	tan(sc2);
	scs_get_d(&resd, sc2);
	return resd;
	break;
    case 3:
	scs_sub(sc2, Pio2_ptr, sc2);
	tan(sc2);
	scs_inv(sc2, sc2);
	scs_get_d(&resd, sc2);
	return -(resd);
	break;
    default:
    fprintf(stderr,"ERREUR: %d is not a valid value in cotan_rn \n", N);
    return x;
    }
 }


/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  -INFINITY
 *************************************************************
 *************************************************************/
double cotan_rd(double x){  
  scs_t sc1, sc2;
  double resd;
  int N;
 
  scs_set_d(sc1, x);
  N = rem_pio2_scs(sc2, sc1); 	/* x is in [0,Pi/4] */ 
 
    switch (N){
    case 0:
	tan(sc2);
	scs_inv(sc2, sc2);
	scs_get_d_minf(&resd, sc2);
	return resd;
   	break;
    case 1:
	tan(sc2);
	scs_get_d_pinf(&resd, sc2);
	return -(resd);
	break;
    case 2:
	scs_sub(sc2, Pio2_ptr, sc2);
	tan(sc2);
	scs_get_d_minf(&resd, sc2);
	return resd;
	break;
    case 3:
	scs_sub(sc2, Pio2_ptr, sc2);
	tan(sc2);
	scs_inv(sc2, sc2);
	scs_get_d_pinf(&resd, sc2);
	return -(resd);
	break;
    default:
    fprintf(stderr,"ERREUR: %d is not a valid value in sn_cotan \n", N);
    return x;
    }
 
}



/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  +INFINITY
 *************************************************************
 *************************************************************/
double cotan_ru(double x){  
  scs_t sc1, sc2;
  double resd;
  int N;

  scs_set_d(sc1, x);
  N = rem_pio2_scs(sc2, sc1); 	/* x is in [0,Pi/4] */ 
   
   switch (N){
    case 0:
	tan(sc2);
	scs_inv(sc2, sc2);
	scs_get_d_pinf(&resd, sc2);
	return resd;
   	break;
    case 1:
	tan(sc2);
	scs_get_d_minf(&resd, sc2);
	return -(resd);
	break;
    case 2:
	scs_sub(sc2, Pio2_ptr, sc2);
	tan(sc2);
	scs_get_d_pinf(&resd, sc2);
	return resd;
	break;
    case 3:
	scs_sub(sc2, Pio2_ptr, sc2);
	tan(sc2);
	scs_inv(sc2, sc2);
	scs_get_d_minf(&resd, sc2);
	return -(resd);
	break;
    default:
    fprintf(stderr,"ERREUR: %d is not a valid value in sn_cotan \n", N);
    return x;
    }     
}

