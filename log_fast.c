/*
 * Function to compute the logarithm with fully exact rounding
 *
 * Author : Daramy Catherine, Florent de Dinechin
 * (Catherine.Daramy,Florent.de.Dinechin@ens-lyon.fr)
 *
 * Date of creation : 26/08/2003   
 * Last Modified    : 16/02/2002
 */
#include <stdio.h>
#include <stdlib.h>
#include "log_fast.h"
#include <crlibm.h>
/*
 *  1) First reduction: exponent extraction      
 *         E  
 *  x = 2^   .(y)    with  1 <= y < 2
 *
 *  log(x) = E.log(2) + log(y) where:
 *     - log(2)   is tabulated
 *     - log(y) need to be evaluated 
 *  
 *
 *  2) Avoiding accuracy problem when E=-1 by testing
 *   
 *    if (ny >= sqrt(2)) then 
 *        y = z/2;  E = E+1; 
 *    and,
 *        log(x) = (E+1).log(2) + log(y/2)
 *
 *    so now:    11/16 <= sqrt(2)/2 <= y < sqrt(2) <= 23/16
 *
 *
 *  3) Second reduction: tabular reduction
 *                    
 *     The interval 1/sqrt(2) .. sqrt(2) is divided in 8 intervals.                                 
 *     So, find the interval X_i where y is.
 *     And compute z = y - middle(X_i);
 *                                    
 *  4) Computation:
 *   
 *     Polynomial evaluation of:
 *        - P(z) ~ log(z+middle(X_i))
 *
 *                   -4      -5
 *   with  |z| < 2^   or 2^   depending the considered interval.
 *
 *
 *  5) Reconstruction:
 *   log(x) = E.log(2) + P(z)
 *
 */

/* The prototypes of the second step */
extern double scs_log_rn(db_number, int);
extern double scs_log_ru(db_number, int);
extern double scs_log_rd(db_number, int);


/*  Compile-time switches */
/* These #defines are intended to provide several performance tradeoffs. */

/* switches on various printfs. Default 0 */
#define DEBUG 0

/* The code for the first step may be shared, or not. Default to 1 (shared)
   There is a speed/size tradeoff there.
   Example : sharing costs 30 cycles out of 1300 on a PIV/Linux/gcc3.3 */
#define SHARE_CODE 1

/* The following constant selects which path to take.  Currently, 8
is optimal, but if we improve the speed of the second step in the
future, we may want to try 8 or 7, which will trade off speed for the
first step and % of taking second step
*/
#define CONST_FASTPATH 8

/* constant for directed rounding, do not edit */
#define DRCST_FASTPATH (1. / ((double)(2<<CONST_FASTPATH)) )

/* constant for rounding to nearest, do not edit */
#define RNCST_FASTPATH (1.0000000000000009 +(4. / (((double)(4<<CONST_FASTPATH)) -1.) ))



#if SHARE_CODE
static void log_quick(double *pres_hi, double *pres_lo, double * proundcst, db_number * py, int * pE) {
   db_number z,y;
   double ln2_times_E_HI, ln2_times_E_LO, res_hi, res_lo, roundcst;
   double res, P_hi, P_lo;
   int k, i, E;

   y = *py;
   E=*pE;
   
    /* E belongs to {-52-1023 .. 2046-1023} */
    
    /* find the interval including y.d */
    i = (((y.i[HI_ENDIAN] & 0x001F0000)>>16)-6) ;  /* to know which polynom to evaluate */
    if (i < 10)
      i = i>>1;
    else
      i = ((i-1)>>1);

    z.d = y.d - (middle[i]).d; 	/* evaluate the value of x in the
				   ii-th interval (exact thanks to
				   Sterbenz Lemma) */
    
    /* sc_ln2_times_E = E*log(2)  */
    Mul22(&ln2_times_E_HI, &ln2_times_E_LO, ln2hi.d, ln2lo.d, E*1., 0.);
    
    /*
     * Polynomial evaluation of log(1 + R) 
     */

    res = (poly_log_fast_h[i][13]).d;
    for(k=12; k>1; k--){
      res *= z.d;
      res += (poly_log_fast_h[i][k]).d;
    }
   
    if(ABS(ln2_times_E_HI) < ((double)(1<<CONST_FASTPATH))+0.5) {
      
      /* Slow path */
      roundcst = rncst[i];
      
      /* Multiply S2 by x = P2 */
      Mul12(&P_hi, &P_lo, res, z.d);
      
      /* add S1 = a1_hi + a1_lo to P2 */ 
      Add22(&res_hi, &res_lo, (poly_log_fast_h[i][1]).d,  (poly_log_fast_l[i][1]).d, P_hi, P_lo);
    
      /* multiply S1 by x = P1 */ 
      Mul22(&P_hi, &P_lo, res_hi, res_lo, z.d, 0.); 
      
      /* add S0 = a0_hi + a0_lo to P1=P1_hi+P1_lo */
      Add22(&res_hi, &res_lo, (poly_log_fast_h[i][0]).d, (poly_log_fast_l[i][0]).d, P_hi, P_lo);
      
      /* REBUILDING */
      if (!(E==0)) /* This test saves two tests in Add22Cond */
	Add22(&res_hi, &res_lo, ln2_times_E_HI, ln2_times_E_LO, res_hi, res_lo);
    }

    else { /* Fast path */
      roundcst=RNCST_FASTPATH;
      res =  (poly_log_fast_h[i][0]).d + z.d*((poly_log_fast_h[i][1]).d + z.d*res);
      
      /* REBUILDING */
      /* As |ln2_times_E_HI| > CONST_FASTPATH and |res| < 0.5 we may use Add22 */
      Add22(&res_hi, &res_lo, ln2_times_E_HI, ln2_times_E_LO, res, 0.0);
#if DEBUG
      printf("\ni=%d    roundcst=%1.20e\n", i , roundcst);
      printf("\ res=%1.20e\n E=%d\n Eln2HI=%1.20e\n Eln2LO=%1.20e\n \n",res,E,ln2_times_E_HI,ln2_times_E_LO);
      printf("\n   reshi=%1.20e\n   reslo=%1.20e\n",res_hi,res_lo);
      printf("\n    reslo*cst=%1.20e\n",res_lo*roundcst);
      printf("\n roundcst=%1.20e\n",roundcst);
#endif /* DEBUG */
    }

    /* exit */
    *proundcst = roundcst;
    *pres_hi=res_hi;
    *pres_lo=res_lo;
}

#endif /* SHARE_CODE */

/*************************************************************
 *************************************************************
 *               ROUNDED  TO NEAREST			     *
 *************************************************************
 *************************************************************/
 double log_rn(double x){ 
   db_number y;
   double res_hi,res_lo,roundcst;
   int E;

   E=0;
   y.d=x;

 /* Filter cases */
   if (y.i[HI_ENDIAN] < 0x00100000){        /* x < 2^(-1022)    */
     if (((y.i[HI_ENDIAN] & 0x7fffffff)|y.i[LO_ENDIAN])==0){
       return 1.0/0.0;     
     }                    		   /* log(+/-0) = -Inf */
     if (y.i[HI_ENDIAN] < 0){ 
      return (x-x)/0;                      /* log(-x) = Nan    */
     }
     /* Subnormal number */
     E = -52; 		
     y.d *= two52.d; 	  /* make x as normal number = x's mantissa    */ 
   }
    
   if (y.i[HI_ENDIAN] >= 0x7ff00000){
     return  x+x;				    /* Inf or Nan       */
   }
   
    /* find y.d such that sqrt(2)/2 < y.d < sqrt(2) */
    E += (y.i[HI_ENDIAN]>>20)-1023;				/* extract the exponent */
    y.i[HI_ENDIAN] =  (y.i[HI_ENDIAN] & 0x000fffff) | 0x3ff00000;	/* do exponent = 0 */
    if (y.d > SQRT_2){
      y.d *= 0.5;
      E++;
    }

#if SHARE_CODE
  log_quick(&res_hi, &res_lo, &roundcst, &y, &E);
#else
  {
    db_number z;
    double ln2_times_E_HI, ln2_times_E_LO;
    double res, P_hi, P_lo;
    int k, i;

    E=0;
    
    y.d = x;
    /* Filter cases */
  if (y.i[HI_ENDIAN] < 0x00100000){        /* x < 2^(-1022)    */
    if (((y.i[HI_ENDIAN] & 0x7fffffff)|y.i[LO_ENDIAN])==0){
      return 1.0/0.0;     
    }                    		   /* log(+/-0) = -Inf */
    if (y.i[HI_ENDIAN] < 0){ 
      return (x-x)/0;                      /* log(-x) = Nan    */
    }
    /* Subnormal number */
    E = -52; 		
    y.d *= two52.d; 	  /* make x as normal number = x's mantissa    */ 
    }
    
  if (y.i[HI_ENDIAN] >= 0x7ff00000){
    return  x+x;				    /* Inf or Nan       */
  }

    /* find y.d such that sqrt(2)/2 < y.d < sqrt(2) */
    E += (y.i[HI_ENDIAN]>>20)-1023;				/* extract the exponent */
    y.i[HI_ENDIAN] =  (y.i[HI_ENDIAN] & 0x000fffff) | 0x3ff00000;	/* do exponent = 0 */
    if (y.d > SQRT_2){
      y.d *= 0.5;
      E++;
    }

    /* E belongs to {-52-1023 .. 2046-1023} */
    
    /* find the interval including y.d */
    i = (((y.i[HI_ENDIAN] & 0x001F0000)>>16)-6) ;  /* to know which polynom to evaluate */
    if (i < 10)
      i = i>>1;
    else
      i = ((i-1)>>1);

    z.d = y.d - (middle[i]).d; 	/* evaluate the value of x in the
				   ii-th interval (exact thanks to
				   Sterbenz Lemma) */
    
    /* sc_ln2_times_E = E*log(2)  */
    Mul22(&ln2_times_E_HI, &ln2_times_E_LO, ln2hi.d, ln2lo.d, E*1., 0.);
    
    /*
     * Polynomial evaluation of log(1 + R) 
     */

    res = (poly_log_fast_h[i][13]).d;
    for(k=12; k>1; k--){
      res *= z.d;
      res += (poly_log_fast_h[i][k]).d;
    }
   
    if(ABS(ln2_times_E_HI) < ((double)(1<<CONST_FASTPATH))+0.5) {
      
      /* Slow path */
      roundcst = rncst[i];
      
      /* Multiply S2 by x = P2 */
      Mul12(&P_hi, &P_lo, res, z.d);
      
      /* add S1 = a1_hi + a1_lo to P2 */ 
      Add22(&res_hi, &res_lo, (poly_log_fast_h[i][1]).d,  (poly_log_fast_l[i][1]).d, P_hi, P_lo);
    
      /* multiply S1 by x = P1 */ 
      Mul22(&P_hi, &P_lo, res_hi, res_lo, z.d, 0.); 
      
      /* add S0 = a0_hi + a0_lo to P1=P1_hi+P1_lo */
      Add22(&res_hi, &res_lo, (poly_log_fast_h[i][0]).d, (poly_log_fast_l[i][0]).d, P_hi, P_lo);
      
      /* REBUILDING */
      if (!(E==0)) /* This test saves two tests in Add22Cond */
	Add22(&res_hi, &res_lo, ln2_times_E_HI, ln2_times_E_LO, res_hi, res_lo);
    }

    else { /* Fast path */
      roundcst=RNCST_FASTPATH;
      res =  (poly_log_fast_h[i][0]).d + z.d*((poly_log_fast_h[i][1]).d + z.d*res);
      
      /* REBUILDING */
      /* As |ln2_times_E_HI| > CONST_FASTPATH and |res| < 0.5 we may use Add22 */
      Add22(&res_hi, &res_lo, ln2_times_E_HI, ln2_times_E_LO, res, 0.0);
#if DEBUG
      printf("\ni=%d    roundcst=%1.20e\n", i , roundcst);
      printf("\ res=%1.20e\n E=%d\n Eln2HI=%1.20e\n Eln2LO=%1.20e\n \n",res,E,ln2_times_E_HI,ln2_times_E_LO);
      printf("\n   reshi=%1.20e\n   reslo=%1.20e\n",res_hi,res_lo);
      printf("\n    reslo*cst=%1.20e\n",res_lo*roundcst);
      printf("\n roundcst=%1.20e\n",roundcst);
#endif /* DEBUG */
    }
  }
  
#endif /* SHARE_CODE */
  /* ROUNDING TO NEAREST */
  if(res_hi == (res_hi + (res_lo * roundcst)))
    return res_hi;
  else { 
#if DEBUG
    printf("Going for Accurate Phase");
#endif
    return scs_log_rn(y, E);    
  }
 }





/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  -INFINITY		     *
 *************************************************************
 *************************************************************/
 double log_rd(double x){ 
   db_number y;
   double res_hi,res_lo,roundcst;
   int E;

   E=0;
   y.d=x;

 /* Filter cases */
   if (y.i[HI_ENDIAN] < 0x00100000){        /* x < 2^(-1022)    */
     if (((y.i[HI_ENDIAN] & 0x7fffffff)|y.i[LO_ENDIAN])==0){
       return 1.0/0.0;     
     }                    		   /* log(+/-0) = -Inf */
     if (y.i[HI_ENDIAN] < 0){ 
      return (x-x)/0;                      /* log(-x) = Nan    */
     }
     /* Subnormal number */
     E = -52; 		
     y.d *= two52.d; 	  /* make x as normal number = x's mantissa    */ 
   }
    
   if (y.i[HI_ENDIAN] >= 0x7ff00000){
     return  x+x;				    /* Inf or Nan       */
   }
   
    /* find y.d such that sqrt(2)/2 < y.d < sqrt(2) */
    E += (y.i[HI_ENDIAN]>>20)-1023;				/* extract the exponent */
    y.i[HI_ENDIAN] =  (y.i[HI_ENDIAN] & 0x000fffff) | 0x3ff00000;	/* do exponent = 0 */
    if (y.d > SQRT_2){
      y.d *= 0.5;
      E++;
    }
 

#if SHARE_CODE
  log_quick(&res_hi, &res_lo, &roundcst, &y, &E);
#else
  
#endif /* SHARE_CODE */ 
  /* ROUNDING TO - INFINITY */

 { 
   db_number absyh, absyl, u, u53;

   absyh.d=res_hi;
   absyl.d=res_lo;

   absyh.l = absyh.l & 0x7fffffffffffffffLL;
   absyl.l = absyl.l & 0x7fffffffffffffffLL;
   u.l     = absyh.l & 0x7fff000000000000LL;
   u53.l   = u.l     + 0x0035000000000000LL; /* exp + 53  */

   if(absyl.d > 0.5*u53.d){ /* FIXME 0.5*/
     if(res_lo<0.)
       res_hi -= u.d;
    return res_hi;
  }else{
    return scs_log_rd(y, E);
  }
 }
}

/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  +INFINITY		     *
 *************************************************************
 *************************************************************/
double log_ru(double x){ 
 double ln2_times_E_HI, ln2_times_E_LO;
 double res, P_hi, P_lo;
 db_number y, z, reshi, reslo;
 int k, i = 0, E = 0;
    
  y.d = x;
  /* Filter cases */
  if (y.i[HI_ENDIAN] < 0x00100000){        /* x < 2^(-1022)    */
    if (((y.i[HI_ENDIAN] & 0x7fffffff)|y.i[LO_ENDIAN])==0){
      return 1.0/0.0;     
    }                    		   /* log(+/-0) = -Inf */
    if (y.i[HI_ENDIAN] < 0){ 
      return (x-x)/0;                      /* log(-x) = Nan    */
    }
    /* Subnormal number */
    E = -52; 		
    y.d *= two52.d; 	  /* make x as normal number = x's mantissa    */ 
    }
    
  if (y.i[HI_ENDIAN] >= 0x7ff00000){
    return  x+x;				    /* Inf or Nan       */
  }

  /* find y.d such that sqrt(2)/2 < y.d < sqrt(2) */
  E += (y.i[HI_ENDIAN]>>20)-1023;				/* extract the exponent */
  y.i[HI_ENDIAN] =  (y.i[HI_ENDIAN] & 0x000fffff) | 0x3ff00000;	/* do exponent = 0 */
  if (y.d > SQRT_2){
    y.d *= 0.5;
    E++;
  }
  
  
  
  /* find the interval including y.d */
  i = (((y.i[HI_ENDIAN] & 0x001F0000)>>16)-6) ;  /* 11<= i <= 21, then we know which polynom to evaluate */
  if (i < 10) {
    i = i>>1;
  }
  else{
    i = ((i-1)>>1);
  }
 
  z.d = y.d - (middle[i]).d; 	/* evaluate the value of x in the ii-th interval */ 						/* Sterbenz Lemma */
 
  /*
   * Polynomial evaluation of log(1 + R) with an error less than 2^(-60)
   */

  res = (poly_log_fast_h[i][13]).d;
  for(k=12; k>2; k--){
    res *= z.d;
    res += (poly_log_fast_h[i][k]).d;
  }
   
   /* Multiply S3 by x = P3 */
 Mul12(&P_hi, &P_lo, res, z.d);
 
  /* add S2 = a2_hi to P3 */ 
 Add22(&reshi.d, &reslo.d, (poly_log_fast_h[i][2]).d, 0., P_hi, P_lo);
 
  /* Multiply S2 by x = P2 */
 Mul22(&P_hi, &P_lo, reshi.d, reslo.d, z.d, 0.);
 
  /* add S1 = a1_hi + a1_lo to P2 */ 
  Add22(&reshi.d, &reslo.d, (poly_log_fast_h[i][1]).d,  (poly_log_fast_l[i][1]).d, P_hi, P_lo);
 
  /* multiply S1 by x = P1 */ 
  Mul22(&P_hi, &P_lo, reshi.d, reslo.d, z.d, 0.);
       
  /* add S0 = a0_hi + a0_lo to P1=P1_hi+P1_lo */
  Add22(&reshi.d, &reslo.d, (poly_log_fast_h[i][0]).d, (poly_log_fast_l[i][0]).d, P_hi, P_lo);
  
  if (!(E==0)){
  
  /* sc_ln2_times_E = E*log(2)  */
  Mul22(&ln2_times_E_HI, &ln2_times_E_LO, ln2hi.d, ln2lo.d, E*1., 0.);
   
   /* RECONSTRUCTION */
   Add22(&reshi.d, &reslo.d, ln2_times_E_HI, ln2_times_E_LO, reshi.d, reslo.d);
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
      return scs_log_ru(y, E);
    }
  }
}

