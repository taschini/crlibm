/**
 * Function to compute the tan and the Cotan functions on [-pi/2,pi/2]
 *
 * Author : Daramy Catherine  (Catherine.Daramy@ens-lyon.fr)
 *
 * Date of creation : 11/05/2003   
 * Last Modified    : 15/07/2004
 */
#include <stdio.h>
#include <stdlib.h>
#include <crlibm.h>
#include <crlibm_private.h>
/*#include "sine_fast.h"*/
#include "tan_fast.h"
#include "coefpi2.h"


extern double scs_tan_rn(double);
extern double scs_tan_ru(double);
extern double scs_tan_rd(double);

#define DEBUG 1
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
 *              tan(y/2)                   
 * 3) tan(y) = ----------
 *	       1-tan(y/2)
 *
 * 4) Polynomial evaluation of P(y/2), degree   , assuming that |y/2| < Pi/8 < 0.4  
 *                  
 *                                  (-59)  
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
 * - It computes the polynom in one time ...
 */


/*************************************************************
 *************************************************************
 *               ROUNDED  TO NEAREST
 *************************************************************
 *************************************************************/
double tan_rn(double x){  
  double reshi, reslo, sh, sl, ch, cl, cahyh_h, cahyh_l, sahyh_h, sahyh_l, kd, yh, yl, yh2, tc, ts, th, tl, sah, sal, cah, cal;
  db_number y;
  double rnconstant = 1.0502;
  int k, quadrant;

  /* x < 2^-26  => tan(x)~x with accuracy 2^-53.2 */
  y.d = x;
    if((y.i[HI_ENDIAN]&0x7FFFFFFF) < 0x3E4BEAD3){	/* Test if |x| < (1+e)2^(-26) */
    #if DEBUG
      printf("x est plus petit que 2^-26(1+e)\n");
    #endif
      return x;
    }

  /* Argument reduction */  
  /* Compute k */
  DOUBLE2INT(k, x * invpio256.d);
  quadrant = (k>>7)&3;	/* Pi is divided in 4 quarters */	
  kd = (double) k;
  k=(k&127)<<2;
  
  #if DEBUG
    printf("k = %d\n", k);
  #endif
 if(k==0)
      { /* Here we risk a large cancellation on yh+yl; on the other hand we have sa=0 and ca=1*/
	/* all this is exact */
  	Add12 (th,tl,  (x - kd*pio256hi.d),  (kd*mpio256med1.d) ) ;
	/* error on the last multiplication, and on the Add22 */
	Add22 (&yh, &yl, th, tl, kd*mpio256med2.d, kd*mpio256lo2.d) ;
	yh2 = yh*yh;
	
	/* Sine computation */
	ts = yh2 * (s3.d + yh2*(s5.d + yh2*s7.d));
	/* (1+ts)*(yh+yl) is an approx to sin(yh+yl) */
	/* Now we need to compute (1+ts)*(yh+yl) */
	Add12(sh, sl,   yh, yl+ ts*yh);
	
	/* Cosine computation */
	tc = yh2 * (c2.d + yh2*(c4.d + yh2*(c6.d)));
	/* 1+ tc is an approx to cos(yh+yl) */
	Add12(ch, cl, 1., tc);
      }
  else{
 switch (quadrant){
   case(0):
    #if DEBUG
      printf("Case 0\n");
    #endif
   if(k<=(64<<2)) {  /* sah <= cah */
    sah=sincosTable[k].d; /* sin(a), high part */
    sal=sincosTable[k+1].d; /* sin(a), low part */
    cah=sincosTable[k+2].d; /* cos(a), high part */
    cal=sincosTable[k+3].d; /* cos(a), low part */
  } else { /* cah <= sah */
    int k1=(128<<2) - k;
    cah=sincosTable[k1].d; 
    cal=sincosTable[k1+1].d;
    sah=sincosTable[k1+2].d;
    sal=sincosTable[k1+3].d;
  }     
     break;
   case(1):
    #if DEBUG
      printf("Case 1\n");
    #endif
     if(k<=(64<<2)) {  /* sah <= cah */
    cah=-sincosTable[k].d; /* sin(a), high part */
    cal=-sincosTable[k+1].d; /* sin(a), low part */
    sah=sincosTable[k+2].d; /* cos(a), high part */
    sal=sincosTable[k+3].d; /* cos(a), low part */
  } else { /* cah <= sah */
    int k1=(128<<2) - k;
    sah=sincosTable[k1].d; 
    sal=sincosTable[k1+1].d;
    cah=-sincosTable[k1+2].d;
    cal=-sincosTable[k1+3].d;
  }    
     break;
   case(2):
    #if DEBUG
      printf("Case 2\n");
    #endif
      if(k<=(64<<2)) {  /* sah <= cah */
    sah=-sincosTable[k].d; /* sin(a), high part */
    sal=-sincosTable[k+1].d; /* sin(a), low part */
    cah=-sincosTable[k+2].d; /* cos(a), high part */
    cal=-sincosTable[k+3].d; /* cos(a), low part */
  } else { /* cah <= sah */
    int k1=(128<<2) - k;
    cah=-sincosTable[k1].d; 
    cal=-sincosTable[k1+1].d;
    sah=-sincosTable[k1+2].d;
    sal=-sincosTable[k1+3].d;
  }    
   break;
      case(3):
    #if DEBUG
      printf("Case 3\n");
    #endif
     if(k<=(64<<2)) {  /* sah <= cah */
    cah=sincosTable[k].d ; /* sin(a), high part */
    cal=sincosTable[k+1].d; /* sin(a), low part */
    sah=-sincosTable[k+2].d; /* cos(a), high part */
    sal=-sincosTable[k+3].d; /* cos(a), low part */
  } else { /* cah <= sah */
    int k1=(128<<2) - k;
    sah=-sincosTable[k1].d ; 
    sal=-sincosTable[k1+1].d;
    cah=sincosTable[k1+2].d;
    cal=sincosTable[k1+3].d;
  }    
     break;   
   default:
     fprintf(stderr,"ERREUR: %d is not a valid value in sn_tan \n", quadrant);
     return 0.0;
    }
  

      /* Argument reduction  by Cody & Waite algorithm */
      /* all this is exact but the rightmost multiplication */
      Add12 (yh,yl,  (x - kd*pio256hi.d),  (kd*mpio256lo.d) ) ;
      /* Now y_h is in -Pi/512, Pi/512 */
 
      yh2 = yh*yh ;
 
      /* Sine computation, our numerator */      
      Mul12( &cahyh_h, &cahyh_l, cah, yh);
      Add12(th, tl, sah, cahyh_h);
      ts = yh2 * (s3.d + yh2*(s5.d + yh2*s7.d));
      /* (1+ts)*(yh+yl) is an approx to sin(yh+yl) */
       
      tc = yh2 * (c2.d + yh2*(c4.d + yh2*c6.d ));
      /* 1+ tc is an approx to cos(yh+yl) */
      /* now we compute an approximation to cos(a)sin(x) + sin(a)cos(x)   */
      /* read the sine and cos */ 
      
      Add12( sl, sh, th,    tc*sah + (ts*cahyh_h  +(sal + (tl + (cahyh_l  + (cal*yh + cah*yl)) ) ) )  );

/* Cosine computation, our denominator */
    /* now we compute an approximation to cos(a)cos(x) - sin(a)sin(x)   */
       
    Mul12(&sahyh_h,&sahyh_l, sah, yh);

    ts = yh2 * (s3.d + yh2*(s5.d + yh2*s7.d));
    /* (1+ts)*(yh+yl) is an approx to sin(yh+yl) */
    tc = yh2 * (c2.d + yh2*(c4.d + yh2*(c6.d)));
    /* 1+ tc is an approx to cos(yh+yl) */

    Add12(th, tl,  cah, -sahyh_h);
    Add12(ch, cl, th, tc*cah - (ts*sahyh_h -  (cal + (tl  - (sahyh_l + (sal*yh + sah*yl)) )))   );
}

   Div22(&reshi, &reslo, sh, sl, ch, cl);
   printf("reshi = %0.20f\nreslo = %0.20f\n\n", reshi, reslo);
 /*  
   DIV2(reshi, reslo, sh, sl, ch, cl);
   printf("reshi = %0.20f\nreslo = %0.20f\n", reshi, reslo);*/


  /* ROUNDING TO NEAREST */
 
  if(reshi == (reshi + (reslo * 1.3008))){
#if DEBUG
  printf("1ere etape\n");
#endif   
    return reshi;
  }else{ 
#if DEBUG
   printf("SCS! 2eme etape\n");    
#endif
    return scs_tan_rn(x); 
  } 

}


/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  -INFINITY
 *************************************************************
 *************************************************************/
double tan_rd(double x){  
return 2.114;
 }

/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  +INFINITY
 *************************************************************
 *************************************************************/
double tan_ru(double x){  
return 2.114;
 }
