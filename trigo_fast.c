#include <stdio.h>
#include <stdlib.h>
#include <crlibm.h>
#include <crlibm_private.h>
#include "trigo_fast.h"
#include "coefpi2.h"


extern double scs_sin_rn(double);
extern double scs_sin_ru(double);
extern double scs_sin_rd(double);
extern double scs_sin_rz(double);
extern double scs_cos_rn(double);
extern double scs_cos_ru(double);
extern double scs_cos_rd(double);
extern double scs_cos_rz(double);
extern double scs_tan_rn(double); 
extern double scs_tan_rd(double);  
extern double scs_tan_ru(double);  
extern double scs_tan_rz(double);  

#define DEBUG 0


static double sah,sal,cah,cal;


static void do_sin(double* sh, double* sl, double yh, double yl) {
  double thi, tlo, cahyh_h, cahyh_l, yh2, ts, tc;

  /* Add optimizations for small yh / k  here */

  yh2 = yh*yh;

  if(sah==0.0)
    { /*  sa=0 and ca=1, which simplifies computations */
      ts = yh2 * (s3.d + yh2*(s5.d + yh2*s7.d));
      /* (1+ts)*(yh+yl) is an approx to sin(yh+yl) */
      /* Now we need to compute (1+ts)*(yh+yl) */
      Add12(*sh,*sl,   yh, yl+ ts*yh);
    }
  else {
   
    Mul12(&cahyh_h,&cahyh_l, cah, yh);
    Add12(thi, tlo, sah,cahyh_h);
    
    ts = yh2 * (s3.d + yh2*(s5.d + yh2*s7.d));
    /* (1+ts)*(yh+yl) is an approx to sin(yh+yl) */
    
    tc = yh2 * (c2.d + yh2*(c4.d + yh2*c6.d ));
    /* 1+ tc is an approx to cos(yh+yl) */

    /* now we compute an approximation to cos(a)sin(x) + sin(a)cos(x)   */
    tlo = tc*sah + (ts*cahyh_h  +(sal + (tlo + (cahyh_l  + (cal*yh + cah*yl))))) ;
    Add12(*sh,*sl,  thi, tlo );
  }
}


static void do_cos(double* ch, double* cl, double yh, double yl) {
  double yh2, ts, tc, thi, tlo, sahyh_h,sahyh_l; 

  yh2 = yh*yh ;

  if(sah==0.0)
    { /*  sa=0 and ca=1, which simplifies computations */
    tc = yh2 * (c2.d + yh2*(c4.d + yh2*c6.d ));
    /* 1+ tc is an approx to cos(yh+yl) */

      /* Now we need to compute 1+tc */
      Add12(*ch,*cl, 1., tc);
    }
  else {
  
  /* now we compute an approximation to cos(a)cos(x) - sin(a)sin(x)   */
  
  Mul12(&sahyh_h,&sahyh_l, sah, yh);

  ts = yh2 * (s3.d + yh2*(s5.d + yh2*s7.d));
  /* (1+ts)*(yh+yl) is an approx to sin(yh+yl) */

  tc = yh2 * (c2.d + yh2*(c4.d + yh2*(c6.d)));
  /* 1+ tc is an approx to cos(yh+yl) */
  
  Add12(thi, tlo,  cah, -sahyh_h);
  tlo = tc*cah - (ts*sahyh_h -  (cal + (tlo  - (sahyh_l + (sal*yh + sah*yl)) ))) ;
  Add12(*ch, *cl,    thi, tlo );
  }
}







int static trig_range_reduction(double* pyh, double* pyl, 
				double x, int absxhi, 
				double (*scs_fun)(double)   ) {
  int k;
  double kd;
  if  (absxhi < XMAX_CODY_WAITE_3) {
    DOUBLE2INT(k, x * INV_PIO256);
    kd = (double) k;
    if(((k&127) == 0)) { 
      /* Here we risk a large cancellation on yh+yl; 
	 on the other hand we will have sa=0 and ca=1*/
      double kch_h,kch_l, kcm_h,kcm_l,  th, tl;
      /* TODO : improve this code by pre-splitting CH,  CM and k (as an int) 
	 Then you can improve the precision by taking kmax into account */
      /* all this is exact */
      Mul12(&kch_h, &kch_l,   kd, RR_DD_MCH);
      Mul12(&kcm_h, &kcm_l,   kd, RR_DD_MCM);
      Add12 (th,tl,  kch_l, kcm_h) ;
      /* only rounding error in the last multiplication and addition */ 
      Add22 (pyh, pyl,    (x + kch_h) , (kcm_l - kd*RR_DD_CL),   th, tl) ;
    } 
    else {      
      /* Argument reduction  by Cody & Waite algorithm */ 
      /* Here we do not care about cancellations on *pyh+yl */
      if (absxhi < XMAX_CODY_WAITE_2) { 
	/* all this is exact but the rightmost multiplication */
	Add12 (*pyh,*pyl,  (x - kd*RR_CW2_CH),  (kd*RR_CW2_MCL) ) ;
      }
     else 
       /* all this is exact but the rightmost multiplication */
       Add12Cond(*pyh,*pyl,  (x - kd*RR_CW3_CH) -  kd*RR_CW3_CM,   kd*RR_CW3_MCL);
    }
  }
  else  if ( absxhi < XMAX_DDRR ) {
    long long int kl;
    double kch_h,kch_l, kcm_h,kcm_l,  th, tl;
    DOUBLE2LONGINT(kl, x*INV_PIO256);
    kd=(double)kl;
    k = (int) kl;
#if DEBUG
    printf("kl=%lld  \n", kl);
#endif
    if((k&127) == 0) { 
      scs_t X, Y,Yh,Yl;
      scs_set_d(X, x*128.0); 
      k= rem_pio2_scs(Y, X);
      /* TODO an optimized procedure for the following */
      scs_get_d(pyh, Y);
      scs_set_d(Yh, *pyh);
      scs_sub(Yl, Y,Yh);
      scs_get_d(pyl, Yl);
      *pyh = *pyh * (1./128.) ;
      *pyl = *pyl * (1./128.) ;
    } 
    else {
      /* all this is exact */
      Mul12(&kch_h, &kch_l,   kd, RR_DD_MCH);
      Mul12(&kcm_h, &kcm_l,   kd, RR_DD_MCM);
      Add12 (th,tl,  kch_l, kcm_h) ;
      /* only rounding error in the last multiplication and addition */ 
      Add22 (pyh, pyl,    (x + kch_h) , (kcm_l - kd*RR_DD_CL),   th, tl) ;
    }
  }
  else {
    scs_t X, Y,Yh,Yl;
    if (absxhi > 0x7F700000) /*2^(1023-7)*/
      return (*scs_fun)(x);
    else {      
      scs_set_d(X, x*128.0); 
      k= rem_pio2_scs(Y, X);
      /* TODO an optimized procedure for the following */
      scs_get_d(pyh, Y);
      scs_set_d(Yh, *pyh);
      scs_sub(Yl, Y,Yh);
      scs_get_d(pyl, Yl);
      *pyh = *pyh * (1./128.) ;
      *pyl = *pyl * (1./128.) ;
    } 
  }
 return k;
}








/*************************************************************
 *************************************************************
 *              SIN ROUNDED  TO NEAREST			     *
 *************************************************************
 *************************************************************/ 

double sin_rn(double x){ 
  double sh, sl, yh, yl, xx;
  int quadrant;
  int k;
  int absxhi;
  db_number x_split;
  double rx, sx, cx, th, tl, ts, gh, gl; 

  x_split.d=x;
  absxhi = x_split.i[HI_ENDIAN] & 0x7fffffff;

  if (absxhi < XMAX_SIN_FAST2){
    if (absxhi < XMAX_SIN_FAST){
      /* CASE 1 : x small enough sin(x)=x */
      if (absxhi <XMAX_RETURN_X_FOR_SIN)
	return x;

      /* CASE 2 : x < 2^-7
                  Fast polynomial evaluation */
      xx = x*x;
      ts = x * xx * (s3.d + xx*(s5.d + xx*s7.d ));
      Add12(sh,sl, x, ts);
      if(sh == (sh + (sl * RN_CST_SINFAST1))){	
	return sh;
      }else{ 
	return scs_sin_rn(x); 
      } 
    }
    /* CASE3 : 2^-1 > (x) > 2^-7 > Pi/512 
               easy range reduction (no dramatic cancellation)
	       + table look-up 
               + fast polynomial evaluation */

    /* Cody and Wayte range reduction */
    DOUBLE2INT(k, x * INV_PIO256);
    rx = (x - k*RR_CW2_CH) + k*RR_CW2_MCL;

    quadrant = (k>>7)&3;
    k=(k&127)<<2;
    
    if(k<=(64<<2)) {
      sah=sincosTable[k+0].d; /* sin(a), high part */
      sal=sincosTable[k+1].d; /* sin(a), low part */
      cah=sincosTable[k+2].d; /* cos(a), high part */
      cal=sincosTable[k+3].d; /* cos(a), low part */
    } else { /* cah <= sah */
      int k1=(128<<2) - k;
      cah=sincosTable[k1+0].d; /* cos(a), high part */
      cal=sincosTable[k1+1].d; /* cos(a), low part  */ 
      sah=sincosTable[k1+2].d; /* sin(a), high part */
      sal=sincosTable[k1+3].d; /* sin(a), low part */
    }
    xx = rx*rx;

    sx = rx * xx * (s3.d + xx*(s5.d + xx*s7.d )); // rx is missing to get sin
    /* !!!!!!!! c2.d = 1/2 exactement, à modifier pour augmenter la
       précision de l'evaluation polynomiale !!!!!!!!! */
    cx = xx * (c2.d + xx*(c4.d + xx*c6.d));       //  1 is missing to have cos
    
    ts = (((sal + cah*sx) + sah*cx) + cah*rx);
    Add12(sh, sl, sah, ts);
    /* !!!!!!!!!  ATTENTION: prendre en compte le fait que l'on oublie
       des termes et que l'on fait des additions bizarres !!!!!!!!!!
       d'apres mes calculs 2^(-61) et meme un peu moins */
    if (sh == (sh + (sl * RN_CST_SINFAST2))){	
      return sh; 
    }else{  
      /* Can improve the accuracy of the result by splitting */
      Mul12(&gh, &gl, cah, rx);
      ts = (((gl + cal*rx) + sal) + cah*sx) + sah*cx; 
      Add12(th, tl, gh, ts);
      Add22(&sh, &sl, sah, 0, th, tl);
      
      if (sh == (sh + (sl * RN_CST_SINFAST3)))	
	return sh; 
      else
	return scs_sin_rn(x); 
    }  
  } 
  /* CASE 4: x>2^(-1) */

  /* Otherwise : Range reduction then standard evaluation */
  k=trig_range_reduction(&yh, &yl,  x, absxhi, &scs_sin_rn);
    
  /* Now y_h is in -Pi/512, Pi/512 and k holds the 32 lower bits of an
     int such that x = yh+yl + kPi/256 */
  
  quadrant = (k>>7)&3;
  k=(k&127)<<2;
  
  if(k<=(64<<2)) {
    sah=sincosTable[k+0].d; /* sin(a), high part */
    sal=sincosTable[k+1].d; /* sin(a), low part */
    cah=sincosTable[k+2].d; /* cos(a), high part */
    cal=sincosTable[k+3].d; /* cos(a), low part */
  } else { /* cah <= sah */
    int k1=(128<<2) - k;
    cah=sincosTable[k1+0].d; /* cos(a), high part */
    cal=sincosTable[k1+1].d; /* cos(a), low part  */ 
    sah=sincosTable[k1+2].d; /* sin(a), high part */
    sal=sincosTable[k1+3].d; /* sin(a), low part  */
  }

#if DEBUG
	printf("sah=%1.30e sal=%1.30e  \n", sah,sal);
	printf("cah=%1.30e cal=%1.30e  \n", cah,cal);
#endif

  if (quadrant&1)   /*compute the cos  */
    do_cos(&sh, &sl,  yh,yl);
  else /* compute the sine */
    do_sin(&sh, &sl,  yh,yl);

  
  if(quadrant>=2) { 
    sh = -sh;
    sl = -sl;
  }
  
  if(sh == (sh + (sl * 1.0004))){	
     return sh;
  }else{
    return scs_sin_rn(x); 
  } 

}

/* TODO */
double sin_rd(double x){
return scs_sin_rd(x);
}

/* TODO */
double sin_ru(double x){ 
return scs_sin_ru(x);
}

/* TODO */
double sin_rz(double x){ 
return scs_sin_rz(x);
}

/*************************************************************
 *************************************************************
 *              COS ROUNDED  TO NEAREST			     *
 *************************************************************
 *************************************************************/
double cos_rn(double x){ 
  double ch, cl, yh, yl,  tc;
  int quadrant;
  int k;
  int absxhi;
  db_number xx;


  xx.d=x;
  absxhi = xx.i[HI_ENDIAN] & 0x7fffffff;

  if (absxhi < XMAX_COS_FAST){
    if (absxhi <XMAX_RETURN_1_FOR_COS)
      return 1.;
    /* Fast Taylor series */
    yh=x*x;
    tc = yh * (c2.d + yh*(c4.d + yh*c6.d));
    Add12(ch,cl, 1, tc);
    if(ch == (ch + (cl * RN_CST_COSFAST))){	
      return ch;
    }else{ 
      return scs_cos_rn(x); 
    } 
  }
  
  /* Otherwise : Range reduction then standard evaluation */
  k=trig_range_reduction(&yh, &yl,  x, absxhi, &scs_cos_rn);
    
  /* Now y_h is in -Pi/512, Pi/512 and k holds the 32 lower bits of an
     int such that x = yh+yl + kPi/256 */
  
  quadrant = (k>>7)&3;
  k=(k&127)<<2;
  
  if(k<=(64<<2)) {
    sah=sincosTable[k+0].d; /* sin(a), high part */
    sal=sincosTable[k+1].d; /* sin(a), low part */
    cah=sincosTable[k+2].d; /* cos(a), high part */
    cal=sincosTable[k+3].d; /* cos(a), low part */
  } else { /* cah <= sah */
    int k1=(128<<2) - k;
    cah=sincosTable[k1+0].d; /* cos(a), high part */
    cal=sincosTable[k1+1].d; /* cos(a), low part  */ 
    sah=sincosTable[k1+2].d; /* sin(a), high part */
    sal=sincosTable[k1+3].d; /* sin(a), low part  */
  }


  if (quadrant&1)   /*compute the cos  */
    do_sin(&ch, &cl,  yh,yl);
  else /* compute the sine */
    do_cos(&ch, &cl,  yh,yl);
  
  if((quadrant == 1)||(quadrant == 2)) { 
    ch = -ch;
    cl = -cl;
  }
  
  if(ch == (ch + (cl * 1.0004))){	
     return ch;
  }else{
    return scs_cos_rn(x); 
  } 

}


/* TODO */
double cos_rd(double x){
return scs_cos_rd(x);
}

/* TODO */
double cos_ru(double x){ 
return scs_cos_ru(x);
}

/* TODO */
double cos_rz(double x){ 
return scs_cos_rz(x);
}

/*************************************************************
 *************************************************************
 *              TAN ROUNDED  TO NEAREST			     *
 *************************************************************
 *************************************************************/ 
double tan_rn(double x){  
  double reshi, reslo, sh, sl, ch, cl, kd, yh, yl;
  db_number y;
  int k, quadrant;


  int absxhi;
  db_number xx;

#if INLINE_SINCOS
  double sah,sal,cah,cal,ts,tc;
#endif

  xx.d=x;
  absxhi = xx.i[HI_ENDIAN] & 0x7fffffff;

  /* x < 2^-26  => tan(x)~x with accuracy 2^-53.2 */
  y.d = x;
    if((y.i[HI_ENDIAN]&0x7FFFFFFF) < 0x3E4BEAD3){	/* Test if |x| < (1+e)2^(-26) */
    #if DEBUG
      printf("x est plus petit que 2^-26(1+e)\n");
    #endif
      return x;
    }

    /*TODO Add polynomial for small values here */ 
  
  /* Otherwise : Range reduction then standard evaluation */
  k=trig_range_reduction(&yh, &yl,  x, absxhi, &scs_cos_rn);

  quadrant = (k>>7)&3;	/* Pi is divided in 4 quarters */	
  kd = (double) k;
  k=(k&127)<<2;

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

#if INLINE_SINCOS
DO_SIN(sh,sl);
DO_COS(ch,cl);
#else  
  do_sin(&sh, &sl, yh, yl);
  do_cos(&ch, &cl, yh, yl);
#endif

   Div22(&reshi, &reslo, sh, sl, ch, cl);

  /* ROUNDING TO NEAREST */
 
  if(reshi == (reshi + (reslo * 1.0004))){
    return reshi;
  }else{ 
    return scs_tan_rn(x); 
  } 

}


/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  -INFINITY
 *************************************************************
 *************************************************************/
/* TODO */
double tan_rd(double x){  
return scs_tan_rd(x);
 }

/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  +INFINITY
 *************************************************************
 *************************************************************/
/* TODO */
double tan_ru(double x){  
return scs_tan_ru(x);
 }

/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  ZERO
 *************************************************************
 *************************************************************/
/* TODO */
double tan_rz(double x){  
return scs_tan_rz(x);
 }
