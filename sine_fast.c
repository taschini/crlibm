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
#include "sine_fast.h"
#include <crlibm.h>

extern double scs_sin_rn(double);
extern double scs_sin_ru(double);
extern double scs_sin_rd(double);


/*************************************************************
 *************************************************************
 *               ARGUMENT REDUCTION			     *
 *************************************************************
 *************************************************************/
 
/* Argument reduction for trigonometric functions */
/* by Cody & Waite algorithm */

int rempio2_fast(double y){
int k;

return k;
}

/*************************************************************
 *************************************************************
 *               ROUNDED  TO NEAREST			     *
 *************************************************************
 *************************************************************/
double sin_rn(double x){ 
  double y, reshi = 0, reslo = 0;
  double  varhi, varlo;
  int k;
  
  y = x;
  varhi = y * twoopihi.d;
  varlo = y * twoopilo.d;
  k = (int)(varhi + varlo);

  y= ( x - ((double)k)*pio2hi.d )    -  ((double)k)*pio2lo.d ;
 
  switch (k&3){
  case 0:
    sin_fast(y, reshi, reslo);
  case 1:
    cos_fast(y, reshi, reslo);
  case 2:  
    sin_fast(y, reshi, reslo);
    reshi *= -1;
    reslo *= -1;
  case 3:
    cos_fast(y,reshi, reslo);
    reshi *= -1;
    reslo *= -1;
  default:
    fprintf(stderr,"ERREUR: %d is not a valid value in sin_rn \n", N);
    exit(1);
  }
  
  /* ROUNDING TO NEAREST */

  if(reshi == (reshi + (reslo * (1.0078125)))){	/* 2^-7 = 0.0078125 */
     return reshi;
  }else{
     return scs_sin_rn(y);
  }

}


/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  -INFINITY		     *
 *************************************************************
 *************************************************************/
double sin_rd(double x){
double y;
db_number reshi, reslo;
int N;

y = x;
N = rempio2_fast(y);
switch (N){
  case 0:
    sin_fast(y, reshi, reslo);
  case 1:
    cos_fast(y, reshi, reslo);
  case 2:  
    sin_fast(y, reshi, reslo);
    reshi.d *= -1;
    reslo.d *= -1;
  case 3:
    cos_fast(y,reshi, reslo);
    reshi.d *= -1;
    reslo.d *= -1;
  default:
    fprintf(stderr,"ERREUR: %d is not a valid value in sin_rn \n", N);
    exit(1);
  }

  /* ROUNDING TO - INFINITY */

 {int logA, err;
  err = 59*2^(20);
  logA = (reshi.i[HI_ENDIAN] & 0x7FF00000) - err;
 
  if((reslo.i[HI_ENDIAN] & 0x7FF00000) > logA){
    if((reshi.i[HI_ENDIAN])^(reslo.i[HI_ENDIAN]) < 0){
      reshi.l -= 1-((reshi.i[HI_ENDIAN] >> 31) << 1);
    }
    return reshi.d;
  }else{
    return scs_sin_rd(y);
  }
 }
}

/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  +INFINITY		     *
 *************************************************************
 *************************************************************/
double sin_ru(double x){
double y;
db_number reshi, reslo;
int N;

y = x;
N = rempio2_fast(y);
switch (N){
  case 0:
    sin_fast(y, reshi, reslo);
  case 1:
    cos_fast(y, reshi, reslo);
  case 2:  
    sin_fast(y, reshi, reslo);
    reshi.d *= -1;
    reslo.d *= -1;
  case 3:
    cos_fast(y,reshi, reslo);
    reshi.d *= -1;
    reslo.d *= -1;
  default:
    fprintf(stderr,"ERREUR: %d is not a valid value in sin_rn \n", N);
    exit(1);
  }

   
  /* ROUNDING TO + INFINITY */

  {int logA, err;
   err = 59*2^(20);
   logA = (reshi.i[HI_ENDIAN] & 0x7FF00000) - err;
 
   if((reslo.i[HI_ENDIAN] & 0x7FF00000) > logA){
    if(reslo.i[HI_ENDIAN] > 0){
      reshi.l += 1;
    }
      return reshi.d;
    }else{
      return scs_sin_ru(y);
    }
  }
}

/*************************************************************
 *************************************************************
 *               	FUNCTION SINE       		     *
 *************************************************************
 *************************************************************/


void sin_fast(double x, double rhi, double rlo){
db_number z;
double yy, y, t = 0;
int n, i;

y = x;
N = rempio2_fast(y);
if((N==1)|(N==3)){
    cos_fast(y, rhi, rlo);
}else{

z.d = y;
if (z.i[HI_ENDIAN] < 0x3e500000){	/* z < 2^ -26 */
    rhi = y;
    rlo = 0;
    return;
}
else if(z.i[HI_ENDIAN] < 0x3fd00000 ){	/* z < 2^ -2 */
    n = 7;
}
else if (y < 0x3fe921fb54442d18){	/* z < Pi/4 */
    n = 10;
}
else {
    printf("ERROR, Z must be lower or equal to PI/4\n");
    return;
}

yy = y*y;

    /*Taylor series */
for(i=n; i > 0; i--){
    t = (t + (poly_sin_fast[i]).d) * yy;
}
t = t*y;
rhi = 1 + t;
rlo = (1-rhi) + t;
}
if(N>1){
    rhi*=-1;	
    rlo*=-1;
}
return;	
}
