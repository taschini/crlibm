/*
 *  sin_fast.h
 *  pbcrlibm
 *
 *  Created by Catherine Daramy on Thu Mar 04 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <crlibm.h>
#include <crlibm_private.h>
#include "sine_fast.h"

extern double scs_sin_rn(double);
extern double scs_sin_ru(double);
extern double scs_sin_rd(double);

#define DEBUG 0


 

static void cos_fast(double *prhi, double *prlo, double yh, double yl){
  double yy, t;
  int n, i;


  n = 10;
 
  yy = yh*yh;
  
  /*Taylor series */
  t=yy*(PolyCos[n]).d;
  for(i=n-1; i >= 0; i--){
    t = yy * (t + (PolyCos[i]).d);

  }

  Add12(*prhi, *prlo, 1.0, t);
}





static void sin_fast(double *prhi, double *prlo, double yh, double yl){
  double yy, t, rh, rl;
  int n, i;
  n = 10;
  

  yy = yh*yh;
  
  /*Taylor series */
  t=yy*(PolySin[n]).d;
  for(i=n-1; i >= 0; i--){
    t = yy * (t + (PolySin[i]).d);
  }

  Mul22(&rh, &rl,  t,0.0,  yh, yl);
  Add22(prhi, prlo,  rh, rl, yh, yl);
}




/*************************************************************
 *************************************************************
 *               ROUNDED  TO NEAREST			     *
 *************************************************************
 *************************************************************/

double sin_rn(double x){ 
  double reshi, reslo;
  double kd, yh,yl;
  double absx;
  int k;
  
  if (x<0) 
    absx=-x; 
  else 
    absx=x;

  if (absx > MAX_FAST_SIN)
    return scs_sin_rn(x);
  if (absx < MIN_RETURN_X_FOR_SINE)
    return x;

  DOUBLE2INT(k, x * invpio2.d);
  kd = (double) k;

  /* Argument reduction  by Cody & Waite algorithm */
  Add12 (yh,yl,  ( x - kd*pio2hi.d),  (-kd*pio2lo.d) ) ;
 
#if DEBUG
   printf("\nk=%d k&3=%d \t y=%5.10f\n", k, k&3, yh);
#endif

  switch (k&3){
  case 0:
    sin_fast(&reshi, &reslo, yh, yl);
    break;
  case 1:
    cos_fast(&reshi, &reslo, yh,yl);
    break;
  case 2:  
    sin_fast(&reshi, &reslo, yh,yl);
    reshi *= -1;
    reslo *= -1;
    break;
  case 3:
    cos_fast(&reshi, &reslo, yh,yl);
    reshi *= -1;
    reslo *= -1;
    break;
  default:
    fprintf(stderr,"ERREUR: %d is not a valid value in sin_rn \n", k&3);
    exit(1);
  }
  
  /* ROUNDING TO NEAREST */

  //if(reshi == (reshi + (reslo * (1.0078125)))){	/* 2^-7 = 0.0078125 */
     return reshi;
  /* }else{ */
/*      return scs_sin_rn(y); */
/*   } */

}


