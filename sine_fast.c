#include <stdio.h>
#include <stdlib.h>
#include <crlibm.h>
#include <crlibm_private.h>
#include "sine_fast.h"

extern double scs_sin_rn(double);
extern double scs_sin_ru(double);
extern double scs_sin_rd(double);

#define DEBUG 0

 



/*************************************************************
 *************************************************************
 *               ROUNDED  TO NEAREST			     *
 *************************************************************
 *************************************************************/

double sin_rn(double x){ 
  double reshi, reslo, absx;
  double kd, yh,yl;
  int k,quadrant;
  double sah,sal,cah,cal, yh2, ts, tc;
  
  absx=ABS(x);
  /* TODO remplacer les tests sur absx par des tests sur HI(x) */

#if 1
  if (absx < ((1./16.))){
    if (absx <MIN_RETURN_X_FOR_SINE)
      return x;
    /*  fast polynomial */
  
    /*Taylor series */
    yh=x*x;
    ts = yh * (s3.d + yh*(s5.d + yh*(s7.d + yh*(s9.d))));
    Add12(reshi,reslo, x, ts*x);
    if(reshi == (reshi + (reslo * RND_CST_SINFAST))){	
      return reshi;
    }else{ 
#if DEBUG
      printf("SCS!\n");    
#endif
      return scs_sin_rn(x); 
    } 
  }

#endif
    
  if (absx > MAX_FAST_SIN)
    return scs_sin_rn(x);
  
  /* Compute k */
  DOUBLE2INT(k, x * invpio256.d);
  kd = (double) k;
  quadrant = (k>>7)&3;
  
  /* Argument reduction  by Cody & Waite algorithm */
  
  /* all this is exact but the rightmost multiplication */
  Add12 (yh,yl,  (x - kd*pio256hi.d),  (-kd*pio256lo.d) ) ;

  /* Now y_h is in -Pi/512, Pi/512 */
    

#if DEBUG
  printf("\nx=%1.30e      k=%d (%d)\n yh=%1.30e  yl=%1.30e \n", x, k, k&3, yh,yl);
#endif
  

  k=(k&127)<<2;
  

/* read the sine and cos */ 

  if(k<=(64<<2)) {
    sah=sincosTable[k].d ; /* sin(a), high part */
    sal=sincosTable[k+1].d; /* sin(a), low part */
    cah=sincosTable[k+2].d; /* cos(a), high part */
    cal=sincosTable[k+3].d; /* cos(a), low part */
  } else {
    k=(128<<2) - k;
    cah=sincosTable[k].d ; 
    cal=sincosTable[k+1].d;
    sah=sincosTable[k+2].d;
    sal=sincosTable[k+3].d;
  }    

#if DEBUG
	printf("quadrant=%d  \n", quadrant);
	printf("sah=%1.30e sal=%1.30e  \n", sah,sal);
	printf("cah=%1.30e cal=%1.30e  \n", cah,cal);
	printf("ts=%1.30e tc=%1.30e  \n", ts,tc);
#endif
  if (quadrant&1) { 
    double thi, tlo, sahyh_h,sahyh_l; 
    /*compute the cos  */
    yh2 = yh*yh ;
    ts = yh2 * (s3.d + yh2*(s5.d + yh2*s7.d));
    /* (1+ts)*(yh+yl) is an approx to sin(yh+yl) */
    tc = yh2 * (c2.d + yh2*(c4.d + yh2*(c6.d)));
    /* 1+ tc is an approx to cos(yh+yl) */

    /* now we compute an approximation to cos(a)cos(x) - sin(a)sin(x)   */
    
    /* cos(a)cos(x) ~ (cah+cal)*(1+tc) ~ cah + tc*cah + cal (in magnitude order) */
    /* sin(a)sin(x) ~ (sah+sal)*(yh+yl)*(1+ts)  
       ~  sah*yh + sah*yh*ts + sal*yh + sah*yl 
       + (sal*yh*ts + sah*yl*ts) + sal*yl*ts
       the 3 last terms may be neglected  */
      
 
    Mul12(&sahyh_h,&sahyh_l, sah, yh);
    Add12(thi, tlo,   -sahyh_h,   tc*cah - (ts*sah*yh + (sahyh_l + (sal*yh + sah*yl)) )   );
    /* TODO : check Condition true except when k=0 */ 
    Add22(&reshi,&reslo,  cah,cal, thi, tlo);
  }


  else { /* compute the sine */

    if(k==0)
      { /* Here we risk a large cancellation on yh+yl; on the other hand we have sa=0 and ca=1*/
	double th,tl;
	/* all this is exact */
	Add12 (th,tl,  (x - kd*pio256hi.d),  (-kd*pio256med1.d) ) ;
	/* error on the last multiplication, and on the Add22 */
	Add22 (&yh, &yl, th, tl, -kd*pio256med2.d, -kd*pio256lo2.d) ;

	yh2 = yh*yh;
	ts = yh2 * (s3.d + yh2*(s5.d + yh2*s7.d));
	/* (1+ts)*(yh+yl) is an approx to sin(yh+yl) */
	/* Now we need to compute (1+ts)*(yh+yl) */
	Add12(reshi,reslo, yh, yl+ ts*yh);
      }
    else {
      double thi, tlo, cahyh_h, cahyh_l; 
      yh2 = yh*yh ;
      ts = yh2 * (s3.d + yh2*(s5.d + yh2*s7.d));
      /* (1+ts)*(yh+yl) is an approx to sin(yh+yl) */
      
      tc = yh2 * (c2.d + yh2*(c4.d + yh2*(c6.d)));
      /* 1+ tc is an approx to cos(yh+yl) */
      /* now we compute an approximation to cos(a)sin(x) + sin(a)cos(x)   */
      
      /* sin(a)cos(x) ~ (sah+sal)*(1+tc) ~ sah + tc*sah + sal (in magnitude order) */
      /* cos(a)sin(x) ~ (cah+cal)*(yh+yl)*(1+ts)  
	 ~  cah*yh + cah*yh*ts + cal*yh + cah*yl 
	 + (cal*yh*ts + cah*yl*ts) + cal*yl*ts
	 the 3 last terms may be neglected  */
      Mul12(&cahyh_h,&cahyh_l, cah, yh);
      Add12(thi, tlo,   cahyh_h,   tc*sah + (ts*cah*yh + (cahyh_l + (cal*yh + cah*yl)) )   );
      /* sacy>casy except when k=0, but then the Add22 works as well */ 
      Add22(&reshi,&reslo, sah,sal ,   thi, tlo);
    }
  }

    if(quadrant>=2) {
      reshi *= -1;
      reslo *= -1;
    }



  /* ROUNDING TO NEAREST */

  if(reshi == (reshi + (reslo * 1.0002))){	
     return reshi;
  }else{ 
#if DEBUG
   printf("SCS!\n");    
#endif
    return scs_sin_rn(x); 
  } 

}


