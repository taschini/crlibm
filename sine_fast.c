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
  double yyh,yyl, t;
  double zh,zl,rh,rl;
  int n, i;

  if ((yh<1./8.) && (yh>-1./8.)) {
    /* Fast cos */
    n = (DEG_COS-14)/2;
    yyh=yh*yh;
    t=yyh*(PolyCos[n]).d;
    for(i=n-1; i >= 0; i--){
      t = yyh * (t + (PolyCos[i]).d);
    }
    t = yyh * (t + Coef4Cos[0].d);
    t = yyh * (t -0.5);
    Add12((*prhi), (*prlo), 1.0, t);
  }
  else {
    n = (DEG_COS-6)/2;
    
    Mul22(&yyh,&yyl, yh,yl,  yh, yl);
    
    /*Taylor series */
    t=yyh*(PolyCos[n]).d;
    for(i=n-1; i >= 0; i--){
      t = yyh * (t + (PolyCos[i]).d);
    }

    /* degree 4 */
    /* check the Cond */
    Add22(&zh,&zl,  Coef4Cos[0].d, Coef4Cos[1].d,  t,0.);
    Mul22(&rh,&rl,  zh,zl,  yyh,yyl);

    /* degree 2 */ 
    Add22(&zh,&zl,  -0.5, 0.0,  rh,rl);
    Mul22(&rh,&rl,  zh,zl,  yyh,yyl);
    
    /* degree 0 */
    Add22(prhi, prlo, 1.0, 0.0,  rh,rl);
  } 
}




static void sin_fast(double *prhi, double *prlo, double yh, double yl){
  double yyh,yyl, t;
  double zh,zl,rh,rl;
  int n, i;

  n = (DEG_SIN-5)/2;
  
  Mul22(&yyh,&yyl, yh,yl,  yh, yl);
  
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
  double kd, yh,yl,t;
  double absx, absyh;
  int k;
  
  if (x<0) 
    absx=-x; 
  else 
    absx=x;


  if (absx < (1./16.)){
    if (absx < MIN_RETURN_X_FOR_SINE)
      return x;
    /* else fast polynomial */
    yh=x*x;
  
    /*Taylor series */
    t = yh * (Coef3Sin[0].d + yh*((PolySin[0]).d + yh*((PolySin[1]).d + yh*((PolySin[2]).d))));
    Add12(reshi,reslo, x, t*x);
    if(reshi == (reshi + (reslo * 1.006))){	
      return reshi;
    }else{ 
#if DEBUG
      printf("SCS!\n");    
#endif
      return scs_sin_rn(x); 
    } 
  }


    
  if (absx > MAX_FAST_SIN)
    return scs_sin_rn(x);

  /* Compute k */
  DOUBLE2INT(k, x * invpio2.d);
  kd = (double) k;

  if(k==0) {
    /* To optimize heavily */
    yh=x; 
    yl=0;
    sin_fast(&reshi, &reslo, yh, yl);
  }
  else {
    /* Argument reduction  by Cody & Waite algorithm */
    /* all this is exact */
    Add12 (yh,yl,  (x - kd*pio2hi.d),  (-kd*pio2med.d) ) ;
    /* now an error but really small */
    Add22 (&yh,&yl,  yh,yl,   (-kd*pio2lo.d), 0.0 ) ;


#if DEBUG
    printf("\nx=%1.30e      k=%d (%d)\n yh=%1.30e  yl=%1.30e \n", x, k, k&3, yh,yl);
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
    }
  }  
  /* ROUNDING TO NEAREST */
 
  if(reshi == (reshi + (reslo * RND_CST))){	
     return reshi;
  }else{ 
#if DEBUG
   printf("SCS!\n");    
#endif
    return scs_sin_rn(x); 
  } 

}


