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

#define DEBUG 1


 

static void cos_fast(double *prhi, double *prlo, double yh, double yl){
  double yyh,yyl, t;
  double zh,zl,rh,rl;
  int n, i;

  n = (DEG_COS-6)/2;

  yyh = yh*yh;
  /* TODO check it is OK */
  yyl = 2*yh*yl;
  
  /*Taylor series */
  t=yyh*(PolyCos[n]).d;
  for(i=n-1; i >= 0; i--){
    t = yyh * (t + (PolyCos[i]).d);
  }

  /* degree 4 */
  /* check the Cond */
  Add22Cond(&zh,&zl,  Coef4Cos[0].d, Coef4Cos[1].d,  t,0.);
  Mul22(&rh,&rl,  zh,zl,  yyh,yyl);

  /* degree 2 */ 
  Add22Cond(&zh,&zl,  -0.5, 0.0,  rh,rl);
  Mul22(&rh,&rl,  zh,zl,  yyh,yyl);

  /* degree 0 */
  Add22(prhi, prlo, 1.0, 0.0,  rh,rl);
} 





static void sin_fast(double *prhi, double *prlo, double yh, double yl){
  double yyh,yyl, t;
  double zh,zl,rh,rl;
  int n, i;

  n = (DEG_SIN-5)/2;
  
  yyh = yh*yh;
  yyl = 2*yh*yl;
  
  /*Taylor series */
  t=yyh*(PolySin[n]).d;
  for(i=n-1; i >= 0; i--){ 
    t = yyh * (t + (PolySin[i]).d);
  }

  /* degree 3 */
  /* check the Cond */
  Add22Cond(&zh,&zl,  Coef3Sin[0].d, Coef3Sin[1].d,  t,0.);
  Mul22(&rh,&rl,  zh,zl,  yyh,yyl);

  Add22(&zh,&zl,  1.0, 0.0,  rh, rl);
  Mul22(prhi, prlo,  zh,zl,  yh, yl);
}




/*************************************************************
 *************************************************************
 *               ROUNDED  TO NEAREST			     *
 *************************************************************
 *************************************************************/

double sin_rn(double x){ 
  double reshi, reslo;
  double kd, yh,yl;
  double absx, absyh;
  int k;
  
  if (x<0) 
    absx=-x; 
  else 
    absx=x;

  if (absx > MAX_FAST_SIN)
    return scs_sin_rn(x);
  if (absx < MIN_RETURN_X_FOR_SINE)
    return x;

  /* Compute k */
  DOUBLE2INT(k, x * invpio2.d);
  kd = (double) k;

  /* Argument reduction  by Cody & Waite algorithm */
  Add12 (yh,yl,  (x - kd*pio2hi.d),  (-kd*pio2lo.d) ) ;

#if DEBUG
   printf("\nx=%1.8e      k=%d (%d) \t y=%5.10f\n", x, k, k&3, yh);
#endif
   if(yh<0) absyh=-yh; else absyh=yh;
   if((k!=0) && (absyh<MINY_ACCURATE)){
#if DEBUG
   printf("SCS!\n");    
#endif
    return scs_sin_rn(x);
  }


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
 
  if(reshi == (reshi + (reslo * RND_CST))){	
     return reshi;
  }else{ 
    return scs_sin_rn(x); 
  } 

}


