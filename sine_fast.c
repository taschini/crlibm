#include <stdio.h>
#include <stdlib.h>
#include <crlibm.h>
#include <crlibm_private.h>
#include "sine_fast.h"
#include "coefpi2.h"


extern double scs_sin_rn(double);
extern double scs_sin_ru(double);
extern double scs_sin_rd(double);

#define DEBUG 0

 


 

static void do_sinCodyWaite(double* reshi, double* reslo, double x, int k) {
  double kd, yh,yl;
  double thi, tlo, cahyh_h, cahyh_l; 
  double sah,sal,cah,cal, yh2, ts, tc;

  kd = (double) k;

  k=(k&127)<<2;

  if(k<=(64<<2)) {  /* sah <= cah */
    sah=sincosTable[k].d ; /* sin(a), high part */
    sal=sincosTable[k+1].d; /* sin(a), low part */
    cah=sincosTable[k+2].d; /* cos(a), high part */
    cal=sincosTable[k+3].d; /* cos(a), low part */
  } else { /* cah <= sah */
    int k1=(128<<2) - k;
    cah=sincosTable[k1].d ; 
    cal=sincosTable[k1+1].d;
    sah=sincosTable[k1+2].d;
    sal=sincosTable[k1+3].d;
  }    

    if(k==0)
      { /* Here we risk a large cancellation on yh+yl; on the other hand we have sa=0 and ca=1*/
	double th,tl;
	/* all this is exact */
  	Add12 (th,tl,  (x - kd*pio256hi.d),  (kd*mpio256med1.d) ) ;
	/* error on the last multiplication, and on the Add22 */
	Add22 (&yh, &yl, th, tl, kd*mpio256med2.d, kd*mpio256lo2.d) ;
	yh2 = yh*yh;
	ts = yh2 * (s3.d + yh2*(s5.d + yh2*s7.d));
	/* (1+ts)*(yh+yl) is an approx to sin(yh+yl) */
	/* Now we need to compute (1+ts)*(yh+yl) */
	Add12(*reshi,*reslo,   yh, yl+ ts*yh);
      }
    else {
      /* Argument reduction  by Cody & Waite algorithm */
      /* all this is exact but the rightmost multiplication */
      Add12 (yh,yl,  (x - kd*pio256hi.d),  (kd*mpio256lo.d) ) ;
      /* Now y_h is in -Pi/512, Pi/512 */
      
      Mul12(&cahyh_h,&cahyh_l, cah, yh);

      yh2 = yh*yh ;

      Add12(thi, tlo,     sah,cahyh_h);

      ts = yh2 * (s3.d + yh2*(s5.d + yh2*s7.d));
      /* (1+ts)*(yh+yl) is an approx to sin(yh+yl) */
       
      tc = yh2 * (c2.d + yh2*(c4.d + yh2*c6.d ));
      /* 1+ tc is an approx to cos(yh+yl) */
      /* now we compute an approximation to cos(a)sin(x) + sin(a)cos(x)   */
      /* read the sine and cos */ 
      
      Add12(*reshi,*reslo, thi,    tc*sah + (ts*cahyh_h  +(sal + (tlo + (cahyh_l  + (cal*yh + cah*yl)) ) ) )  );
    }
}





static void do_cosCodyWaite(double* reshi, double* reslo, double x, int k) {
  double kd, yh,yl;
  double sah,sal,cah,cal, yh2, ts, tc;
  double thi, tlo, sahyh_h,sahyh_l; 

  kd = (double) k;
  k=(k&127)<<2;

  if(k<=(64<<2)) {  /* sah <= cah */
    sah=sincosTable[k].d ; /* sin(a), high part */
    sal=sincosTable[k+1].d; /* sin(a), low part */
    cah=sincosTable[k+2].d; /* cos(a), high part */
    cal=sincosTable[k+3].d; /* cos(a), low part */
  } else { /* cah <= sah */
    int k1=(128<<2) - k;
    cah=sincosTable[k1].d ; 
    cal=sincosTable[k1+1].d;
    sah=sincosTable[k1+2].d;
    sal=sincosTable[k1+3].d;
  }    

    /* Argument reduction  by Cody & Waite algorithm */
    /* all this is exact but the rightmost multiplication */
    Add12 (yh,yl,  (x - kd*pio256hi.d),  (kd*mpio256lo.d) ) ;
    
    /* Now y_h is in -Pi/512, Pi/512 */
    yh2 = yh*yh ;

    /* now we compute an approximation to cos(a)cos(x) - sin(a)sin(x)   */
       
    Mul12(&sahyh_h,&sahyh_l, sah, yh);

    ts = yh2 * (s3.d + yh2*(s5.d + yh2*s7.d));
    /* (1+ts)*(yh+yl) is an approx to sin(yh+yl) */
    tc = yh2 * (c2.d + yh2*(c4.d + yh2*(c6.d)));
    /* 1+ tc is an approx to cos(yh+yl) */

    Add12(thi, tlo,  cah, -sahyh_h);
    Add12(*reshi, *reslo, thi,    tc*cah - (ts*sahyh_h -  (cal + (tlo  - (sahyh_l + (sal*yh + sah*yl)) )))   );
}



/*************************************************************
 *************************************************************
 *               ROUNDED  TO NEAREST			     *
 *************************************************************
 *************************************************************/ 

double sin_rn(double x){ 
  double reshi, reslo, yh, ts;
  int k,quadrant;
  int absxhi;
  db_number xx;

  double rnconstant=1.0002;


  xx.d=x;
  absxhi = xx.i[HI_ENDIAN] & 0x7fffffff;

  if (absxhi < XMAX_SIN_FAST){
    if (absxhi <XMAX_RETURN_X_FOR_SIN)
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

    
  if (absxhi > XMAX_CODY_WAITE)
    return scs_sin_rn(x);
  
  /* Compute k */
  DOUBLE2INT(k, x * invpio256.d);
  quadrant = (k>>7)&3;
  
  if (quadrant&1)   /*compute the cos  */
    do_cosCodyWaite(&reshi, &reslo,  x,  k);
  else /* compute the sine */
    do_sinCodyWaite(&reshi, &reslo,  x,  k);
  
  
  if(quadrant>=2) {
    reshi = -reshi;
    reslo = -reslo;
  }
  
#if DEBUG
  printf("\nx=%1.30e      k=%d (%d)\n yh=%1.30e  yl=%1.30e \n", x, k, k&127, yh,yl);
	printf("quadrant=%d  \n", quadrant);
	printf("sah=%1.30e sal=%1.30e  \n", sah,sal);
	printf("cah=%1.30e cal=%1.30e  \n", cah,cal);
	printf("ts=%1.30e tc=%1.30e  \n", ts,tc);
#endif


  /* ROUNDING TO NEAREST */

  if(reshi == (reshi + (reslo * rnconstant))){	
     return reshi;
  }else{ 
#if DEBUG
   printf("SCS!\n");    
#endif
    return scs_sin_rn(x); 
  } 

}


