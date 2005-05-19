/* 
 *this function computes log, correctly rounded, 
 using experimental techniques based on  triple double arithmetics

 THIS IS EXPERIMENTAL SOFTWARE
 
 *
 * Author :  Christoph Lauter
 * christoph.lauter at ens-lyon.fr
 *

 To have it replace the crlibm log, do:

 gcc -DHAVE_CONFIG_H -I.  -fPIC  -O2 -c log-td.c;   mv log-de.o log_fast.o; make 
 
*/


#include <stdio.h>
#include <stdlib.h>
#include "crlibm.h"
#include "crlibm_private.h"
#include "triple-double.h"
#include "log-td.h"

#define AVOID_FMA 0


void log_td_accurate(double *logh, double *logm, double *logl, int E, double ed, int index, double zh, double zl, double logih, double logim) {
  double highPoly, t1h, t1l, t2h, t2l, t3h, t3l, t4h, t4l, t5h, t5l, t6h, t6l, t7h, t7l, t8h, t8l, t9h, t9l, t10h, t10l, t11h, t11l;
  double t12h, t12l, t13h, t13l, t14h, t14l, zSquareh, zSquarem, zSquarel, zCubeh, zCubem, zCubel, higherPolyMultZh, higherPolyMultZm;
  double higherPolyMultZl, zSquareHalfh, zSquareHalfm, zSquareHalfl, polyWithSquareh, polyWithSquarem, polyWithSquarel;
  double polyh, polym, polyl, logil, logyh, logym, logyl;


#if EVAL_PERF
  crlibm_second_step_taken++;
#endif


  /* Accurate phase:

     Argument reduction is already done. 
     We must return logh, logm and logl representing the intermediate result in 118 bits precision.

     We use a 14 degree polynomial, computing the first 3 (the first is 0) coefficients in triple double,
     calculating the next 7 coefficients in double double arithmetics and the last in double.

     We must account for zl starting with the monome of degree 4 (7^3 + 53 - 7 >> 118); so 
     double double calculations won't account for it.

  */


  /* Start of the horner scheme */

#if defined(PROCESSOR_HAS_FMA) && !defined(AVOID_FMA)
  highPoly = FMA(FMA(FMA(FMA(accPolyC14,zh,accPolyC13),zh,accPolyC12),zh,accPolyC11),zh,accPolyC10);
#else
  highPoly = accPolyC10 + zh * (accPolyC11 + zh * (accPolyC12 + zh * (accPolyC13 + zh * accPolyC14)));
#endif
  
  /* We want to write 

     accPolyC3 + zh * (accPoly4 + zh * (accPoly5 + zh * (accPoly6 + zh * (accPoly7 + zh * (accPoly8 + zh * (accPoly9 + zh * highPoly))))));
     (        t14  t13         t12  t11         t10   t9          t8   t7          t6   t5          t4   t3          t2   t1  )

     with all additions and multiplications in double double arithmetics
     but we will produce intermediate results labelled t1h/t1l thru t14h/t14l
  */

  /* A VERIFIER: est-ce qu'on aurait besoin de Mul22 et Add22 conditionnelles? 
             et: Ne vaut il pas la peine d ecrire un Mul212 pour multiplier un double double par un double? 
	         Entre temps on prend on compte quand meme zl puisque ca ne coute rien (sauf au premier z**9)
  */


  Mul12(&t1h, &t1l, zh, highPoly);
  Add22(&t2h, &t2l, accPolyC9h, accPolyC9l, t1h, t1l);
  Mul22(&t3h, &t3l, zh, zl, t2h, t2l);
  Add22(&t4h, &t4l, accPolyC8h, accPolyC8l, t3h, t3l);
  Mul22(&t5h, &t5l, zh, zl, t4h, t4l);
  Add22(&t6h, &t6l, accPolyC7h, accPolyC7l, t5h, t5l);
  Mul22(&t7h, &t7l, zh, zl, t6h, t6l);
  Add22(&t8h, &t8l, accPolyC6h, accPolyC6l, t7h, t7l);
  Mul22(&t9h, &t9l, zh, zl, t8h, t8l);
  Add22(&t10h, &t10l, accPolyC5h, accPolyC5l, t9h, t9l);
  Mul22(&t11h, &t11l, zh, zl, t10h, t10l);
  Add22(&t12h, &t12l, accPolyC4h, accPolyC4l, t11h, t11l);
  Mul22(&t13h, &t13l, zh, zl, t12h, t12l);
  Add22(&t14h, &t14l, accPolyC3h, accPolyC3l, t13h, t13l);

  /* We must now prepare (zh + zl)^2 and (zh + zl)^3 as triple doubles */

  Mul23(&zSquareh, &zSquarem, &zSquarel, zh, zl, zh, zl); // A VERIFIER: c est sous optimal puisqu on calcule deux fois zh*zl 
  Mul233(&zCubeh, &zCubem, &zCubel, zh, zl, zSquareh, zSquarem, zSquarel); // A VERIFIER: Ici aussi on recalule des choses
  
  /* We can now multiplicate the middle and higher polynomial by z^3 */

  Mul233(&higherPolyMultZh, &higherPolyMultZm, &higherPolyMultZl, t14h, t14l, zCubeh, zCubem, zCubel);
  
  /* Multiply now z^2 by -1/2 (exact op) and add to middle and higher polynomial */
  
  zSquareHalfh = zSquareh * -0.5;
  zSquareHalfm = zSquarem * -0.5;
  zSquareHalfl = zSquarel * -0.5;

  Add33(&polyWithSquareh, &polyWithSquarem, &polyWithSquarel, 
	zSquareHalfh, zSquareHalfm, zSquareHalfl, 
	higherPolyMultZh, higherPolyMultZm, higherPolyMultZl);

  /* Add now zh and zl to obtain the polynomial evaluation result */

  Add233(&polyh, &polym, &polyl, zh, zl, polyWithSquareh, polyWithSquarem, polyWithSquarel);

  /* Reconstruct now log(y) = log(1 + z) - log(ri) by adding logih, logim, logil
     logil has not been read to the time, do this first 
  */

  logil =  argredtable[index].logil;
  Add33(&logyh, &logym, &logyl, polyh, polym, polyl, logih, logim, logil);

  /* Multiply log2 with E, i.e. log2h, log2m, log2l by ed 
     ed is always less than 2^(12) and log2h and log2m are stored with at least 12 trailing zeros 
     So multiplying naively is correct (up to 134 bits at least)

     The final result is thus obtained by adding log2 * E to log(y)
  */

  Add33(logh, logm, logl, log2h * ed, log2m * ed, log2l * ed, logyh, logym, logyl);
}



/*************************************************************
 *************************************************************
 *               ROUNDED  TO NEAREST			     *
 *************************************************************
 *************************************************************/
 double log_rn(double x){ 
   db_number xdb, loghdb;
   double res_hi,res_lo, y, ed, ri, logih, logim, yrih, yril, th, zh, zl;
   double polyHorner, zhSquareh, zhSquarel, polyUpper, zhSquareHalfh, zhSquareHalfl;
   double t1h, t1l, t2h, t2l, ph, pl, log2edh, log2edl, logTabPolyh, logTabPolyl, logh, logm, logl, roundcst;
   double zlPlusZhSquareHalflPlusZhZl, polyUpperPluszlh, miulp;
   int E, index;

   E=0;
   xdb.d=x;

   /* Filter cases */
   if (xdb.i[HI] < 0x00100000){        /* x < 2^(-1022)    */
     if (((xdb.i[HI] & 0x7fffffff)|xdb.i[LO])==0){
       return -1.0/0.0;     
     }                    		   /* log(+/-0) = -Inf */
     if (xdb.i[HI] < 0){ 
       return (x-x)/0;                      /* log(-x) = Nan    */
     }
     /* Subnormal number */
     E = -52; 		
     xdb.d *= ((db_number) ((double) two52)).d; 	  /* make x a normal number    */ 
   }
    
   if (xdb.i[HI] >= 0x7ff00000){
     return  x+x;				 /* Inf or Nan       */
   }
   
   
   /* Extract exponent and mantissa 
      Do range reduction,
      yielding to E holding the exponent and
      y the mantissa between sqrt(2)/2 and sqrt(2)
   */
   E += (xdb.i[HI]>>20)-1023;             /* extract the exponent */
   index = (xdb.i[HI] & 0x000fffff);
   xdb.i[HI] =  index | 0x3ff00000;	/* do exponent = 0 */
   index = (index + (1<<(20-L-1))) >> (20-L);
 
   /* reduce  such that sqrt(2)/2 < xdb.d < sqrt(2) */
   if (index >= MAXINDEX){ /* corresponds to xdb>sqrt(2)*/
     xdb.i[HI] -= 0x00100000; 
     E++;
   }
   y = xdb.d;
   index = index & INDEXMASK;
   // Cast integer E into double ed for multiplication later
   ed = (double) E;

   /* Read tables:
      Read one float for ri
      Read the first two doubles for -log(r_i) (out of three)

      Organization of the table:

      one struct entry per index, the struct entry containing 
      r, logih, logim and logil in this order
   */
   

   ri = argredtable[index].ri;
   // Actually we don't need the logarithm entries now
   // Move the following two lines to the eventual reconstruction
   // As long as we don't have any if in the following code, we can overlap 
   // memory access with calculations (static const?)
   logih = argredtable[index].logih;
   logim = argredtable[index].logim;

   /* Do range reduction:

      zh + zl = y * ri - 1.0 correctly

      Correctness is assured by use of Mul12 and Add12
      even if we don't force ri to have its' LSBs set to zero

      Discard zl for higher monome degrees
   */

   Mul12(&yrih, &yril, y, ri);
   th = yrih - 1.0; // Should be exact (Sterbenz) A VERIFIER !!!
   Add12Cond(zh, zl, th, yril); // Works perhaps without the condition A VERIFIER !!!

   /* Polynomial evaluation

      Use a 7 degree polynomial
      Evaluate the higher 5 terms in double precision (-7 * 3 = -21) using Horner's scheme
      Evaluate the lower 3 terms (the last is 0) in double double precision accounting also for zl
      using an ad hoc method

   */



#if defined(PROCESSOR_HAS_FMA) && !defined(AVOID_FMA)
   polyHorner = FMA(FMA(FMA(FMA(c7,zh,c6),zh,c5),zh,c4),zh,c3);
#else
   polyHorner = c3 + zh * (c4 + zh * (c5 + zh * (c6 + zh * c7)));
#endif

   Mul12(&zhSquareh, &zhSquarel, zh, zh);
   polyUpper = polyHorner * (zh * zhSquareh);
   zhSquareHalfh = zhSquareh * -0.5;
   zhSquareHalfl = zhSquarel * -0.5;
   Add12(t1h, t1l, polyUpper, -1 * (zh * zl));
   Add22(&t2h, &t2l, zh, zl, zhSquareHalfh, zhSquareHalfl);
   Add22(&ph, &pl, t2h, t2l, t1h, t1l);

   /* Reconstruction 

      Read logih and logim in the tables (already done)
      
      Compute log(x) = E * log(2) + log(1+z) - log(ri)
      i.e. log(x) = ed * (log2h + log2m) + (ph + pl) + (logih + logim) + delta

      Carry out everything in double double precision

   */
   
   /* We store log2 as log2h + log2m + log2l where log2h and log2m have 12 trailing zeros
      Multiplication of ed (double E) and log2h is thus correct
      The overall accuracy of log2h + log2m + log2l is 53 * 3 - 24 = 135 which
      is enough for the accurate phase
      The accuracy suffices also for the qucik phase: 53 * 2 - 24 = 82
      Condition is verified as log2m is smaller than log2h and both are scaled by ed
   */

   Add12(log2edh, log2edl, log2h * ed, log2m * ed);

   /* Add logih and logim to ph and pl 

      We must use conditioned Add22 as logih can move over ph

      A VERIFIER: on peut quand meme se passer de la condition
                  puisque le polynome reste par definition plus 
		  petit que la valeur tabluee sauf en 0 et la, 
		  le Add22 sans conditions passe peut etre quand meme

   */

   Add22Cond(&logTabPolyh, &logTabPolyl, logih, logim, ph, pl);

   /* Add log2edh + log2edl to logTabPolyh + logTabPolyl */

   Add22Cond(&logh, &logm, log2edh, log2edl, logTabPolyh, logTabPolyl);

   /* Rounding test and eventual return or call to the accurate function 

      Question: pourquoi les tests d arrondi des modes diriges sont factorises 
      mais pas celui du RN?

      Et: ne perd-on pas plein de cycles si on appelle une fonction pour la
      phase exacte au lieu de la recopier ou de la mettre en DEFINE?

   */

   if(E==0)
      roundcst = ROUNDCST1;
   else
      roundcst = ROUNDCST2;


   if(logh == (logh + (logm * roundcst)))
     return logh;
   else 
     {
       
#if DEBUG
       printf("Going for Accurate Phase for x=%1.50e\n",x);
#endif

       log_td_accurate(&logh, &logm, &logl, E, ed, index, zh, zl, logih, logim); 

       /*
	 Our renormalization procedure does not guarantee that logh = round(logh + logm + logl)
	 Further we have worst cases that are worser than 106 bits, so we can be in a 
	 situation where logm = +/- mi-ulp(logh) and we can not guarantee that logl and logm do
	 not overlap, so we have to take logl into account in these cases.
	 
	 We start by generating a mi-ulp of logh. We will not know its sign.
       */
       
       loghdb.d = logh;
       loghdb.l++;
       miulp = (loghdb.d - logh) / 2.0;
       
       /* We determine now if we are in a case worser than 106 bits, i.e. if logm = mi-ulp(logh) */
       
       if (((logm + miulp) != 0.0) && ((logm - miulp) != 0.0)) return logh;
       
       /* If we are here, we are in a case worser than 106 bits or in which logm and logl do overlap. */
       
       if (logh > 0.0) {
	 if (logm > 0.0) {
	   if (logl > 0.0) {
	     /* + + +
		We should have rounded up as we are greater than the middle and coming upwards.
		We have miulp positive as logh is greater than 0.
	     */
	     return loghdb.d;
	   } else {
	     /* + + -
		We have already rounded correctly as we are lesser than the middle and coming upwards.
	     */
	     return logh;
	   }
	 } else {
	   if (logl > 0.0) {
	     /* + - + 
		We have already rounded correctly as we are greater than the middle and coming downwards.
	     */
	     return logl;
	   } else {
	     /* + - -
		We should have rounded down as we are lesser than the middle and coming downwards.
		We had miulp positive as logh is greater than 0.
	     */
	     loghdb.l = loghdb.l - 2;
	     return loghdb.d;
	   }
	 }
       } else {
	 if (logm > 0.0) {
	   if (logl > 0.0) {
	     /* - + + 
		We should have rounded up as we are greater than the middle and coming upwards.
		We had milup negative as logh is lesser than 0.
	     */
	     loghdb.l = loghdb.l - 2;
	     return loghdb.d;
	   } else {
	     /* - + -
		We have already rounded correctly as we are lesser than the middle and coming upwards.
	     */
	     return logh;
	   }
	 } else {
	   if (logl > 0.0) {
	     /* - - +
		We have already rounded correctly as we are greater than the middle and coming downwards.
	     */
	     return logh;
	   } else {
	     /* - - -
		We should have rounded down as we are lesser than the middle and coming downwards.
		We have miulp negative as logh is lesser than 0.
	     */ 
	     return loghdb.d;
	   }
	 }
       }
     } // Accurate phase launched
 }


/*************************************************************
 *************************************************************
 *               ROUNDED  UPWARDS			     *
 *************************************************************
 *************************************************************/
 double log_ru(double x) { 
   db_number xdb, loghdb;
   double res_hi,res_lo, y, ed, ri, logih, logim, yrih, yril, th, zh, zl;
   double polyHorner, zhSquareh, zhSquarel, polyUpper, zhSquareHalfh, zhSquareHalfl;
   double t1h, t1l, t2h, t2l, ph, pl, log2edh, log2edl, logTabPolyh, logTabPolyl, logh, logm, logl;
   double zlPlusZhSquareHalflPlusZhZl, polyUpperPluszlh, roundcst;
   int E, index;

   E=0;
   xdb.d=x;

   if (x == 1.0) return 0.0; // This the only case in which the image under log of a double is a double. 

   /* Filter cases */
   if (xdb.i[HI] < 0x00100000){        /* x < 2^(-1022)    */
     if (((xdb.i[HI] & 0x7fffffff)|xdb.i[LO])==0){
       return -1.0/0.0;     
     }                    		   /* log(+/-0) = -Inf */
     if (xdb.i[HI] < 0){ 
       return (x-x)/0;                      /* log(-x) = Nan    */
     }
     /* Subnormal number */
     E = -52; 		
     xdb.d *= ((db_number) ((double) two52)).d; 	  /* make x a normal number    */ 
   }
    
   if (xdb.i[HI] >= 0x7ff00000){
     return  x+x;				 /* Inf or Nan       */
   }
   
   /* Extract exponent and mantissa 
      Do range reduction,
      yielding to E holding the exponent and
      y the mantissa between sqrt(2)/2 and sqrt(2)
   */
   E += (xdb.i[HI]>>20)-1023;             /* extract the exponent */
   index = (xdb.i[HI] & 0x000fffff);
   xdb.i[HI] =  index | 0x3ff00000;	/* do exponent = 0 */
   index = (index + (1<<(20-L-1))) >> (20-L);
 
   /* reduce  such that sqrt(2)/2 < xdb.d < sqrt(2) */
   if (index >= MAXINDEX){ /* corresponds to xdb>sqrt(2)*/
     xdb.i[HI] -= 0x00100000; 
     E++;
   }
   y = xdb.d;
   index = index & INDEXMASK;
   // Cast integer E into double ed for multiplication later
   ed = (double) E;

   /* Read tables:
      Read one float for ri
      Read the first two doubles for -log(r_i) (out of three)

      Organization of the table:

      one struct entry per index, the struct entry containing 
      r, logih, logim and logil in this order
   */
   

   ri = argredtable[index].ri;
   // Actually we don't need the logarithm entries now
   // Move the following two lines to the eventual reconstruction
   // As long as we don't have any if in the following code, we can overlap 
   // memory access with calculations (static const?)
   logih = argredtable[index].logih;
   logim = argredtable[index].logim;

   /* Do range reduction:

      zh + zl = y * ri - 1.0 correctly

      Correctness is assured by use of Mul12 and Add12
      even if we don't force ri to have its' LSBs set to zero

      Discard zl for higher monome degrees
   */

   Mul12(&yrih, &yril, y, ri);
   th = yrih - 1.0; // Should be exact (Sterbenz) A VERIFIER !!!
   Add12Cond(zh, zl, th, yril); // Works perhaps without the condition A VERIFIER !!!

   /* Polynomial evaluation

      Use a 7 degree polynomial
      Evaluate the higher 5 terms in double precision (-7 * 3 = -21) using Horner's scheme
      Evaluate the lower 3 terms (the last is 0) in double double precision accounting also for zl
      using an ad hoc method

   */


#if defined(PROCESSOR_HAS_FMA) && !defined(AVOID_FMA)
   polyHorner = FMA(FMA(FMA(FMA(c7,zh,c6),zh,c5),zh,c4),zh,c3);
#else
   polyHorner = c3 + zh * (c4 + zh * (c5 + zh * (c6 + zh * c7)));
#endif

   Mul12(&zhSquareh, &zhSquarel, zh, zh);
   polyUpper = polyHorner * (zh * zhSquareh);
   zhSquareHalfh = zhSquareh * -0.5;
   zhSquareHalfl = zhSquarel * -0.5;
   Add12(t1h, t1l, polyUpper, -1 * (zh * zl));
   Add22(&t2h, &t2l, zh, zl, zhSquareHalfh, zhSquareHalfl);
   Add22(&ph, &pl, t2h, t2l, t1h, t1l);

   /* Reconstruction 

      Read logih and logim in the tables (already done)
      
      Compute log(x) = E * log(2) + log(1+z) - log(ri)
      i.e. log(x) = ed * (log2h + log2m) + (ph + pl) + (logih + logim) + delta

      Carry out everything in double double precision

   */
   
   /* We store log2 as log2h + log2m + log2l where log2h and log2m have 12 trailing zeros
      Multiplication of ed (double E) and log2h is thus correct
      The overall accuracy of log2h + log2m + log2l is 53 * 3 - 24 = 135 which
      is enough for the accurate phase
      The accuracy suffices also for the qucik phase: 53 * 2 - 24 = 82
      Condition is verified as log2m is smaller than log2h and both are scaled by ed
   */

   Add12(log2edh, log2edl, log2h * ed, log2m * ed);

   /* Add logih and logim to ph and pl 

      We must use conditioned Add22 as logih can move over ph

      A VERIFIER: on peut quand meme se passer de la condition
                  puisque le polynome reste par definition plus 
		  petit que la valeur tabluee sauf en 0 et la, 
		  le Add22 sans conditions passe peut etre quand meme

   */

   Add22Cond(&logTabPolyh, &logTabPolyl, logih, logim, ph, pl);

   /* Add log2edh + log2edl to logTabPolyh + logTabPolyl */

   Add22Cond(&logh, &logm, log2edh, log2edl, logTabPolyh, logTabPolyl);

   /* Rounding test and eventual return or call to the accurate function 

      Question: pourquoi les tests d arrondi des modes diriges sont factorises 
      mais pas celui du RN?

      Et: ne perd-on pas plein de cycles si on appelle une fonction pour la
      phase exacte au lieu de la recopier ou de la mettre en DEFINE?

   */

   if(E==0)
      roundcst = RDROUNDCST1;
   else
      roundcst = RDROUNDCST2;

   TEST_AND_RETURN_RU(logh, logm, roundcst);

#if DEBUG
  printf("Going for Accurate Phase for x=%1.50e\n",x);
#endif

    log_td_accurate(&logh, &logm, &logl, E, ed, index, zh, zl, logih, logim); 

    /* 
       We will suppose the following:
       - logh and logm do not overlap in any case
       - logm and logl may overlap but logl will not "overtake" logm for overlapping with logh
       - logh may be different from round-to-nearest(logh + logm + logl) 
       - logm may be equal to 0 while logl is not; 
         either one of them is supposed to be different from 0 assured by Lindemann and Lefevre
       - in no case, logh is equal to zero; in particular, the argument of log is not 1.0

       In order to overcome the fact that overlapping of logm and logl is possible or that logm is equal
       to 0, we will first add logl to logm by a double precision add. This operation will only affect the
       sign of logm if it is equal to 0 prior to the addition; it will change the value of logm only
       if there is overlapping. It will not generate an overlapping of logm and logh.
       We will now know that logm = round-to-nearest(logm + logl)
    */

    logm = logm + logl;

    /* We can now simply test the sign of logh and logm */

    if (logm < 0.0) {
      return logh;
    } else {
      loghdb.d = logh;
      if (logh > 0.0) {
	loghdb.l++;
	return loghdb.d;
      } else {
	loghdb.l--;
	return loghdb.d;
      }
    }
 } 


/*************************************************************
 *************************************************************
 *               ROUNDED  DOWNWARDS			     *
 *************************************************************
 *************************************************************/
 double log_rd(double x) { 
   db_number xdb, loghdb;
   double res_hi,res_lo, y, ed, ri, logih, logim, yrih, yril, th, zh, zl;
   double polyHorner, zhSquareh, zhSquarel, polyUpper, zhSquareHalfh, zhSquareHalfl;
   double t1h, t1l, t2h, t2l, ph, pl, log2edh, log2edl, logTabPolyh, logTabPolyl, logh, logm, logl;
   double zlPlusZhSquareHalflPlusZhZl, polyUpperPluszlh, roundcst;
   int E, index;

   E=0;
   xdb.d=x;

   if (x == 1.0) return 0.0; // This the only case in which the image under log of a double is a double. 

   /* Filter cases */
   if (xdb.i[HI] < 0x00100000){        /* x < 2^(-1022)    */
     if (((xdb.i[HI] & 0x7fffffff)|xdb.i[LO])==0){
       return -1.0/0.0;     
     }                    		   /* log(+/-0) = -Inf */
     if (xdb.i[HI] < 0){ 
       return (x-x)/0;                      /* log(-x) = Nan    */
     }
     /* Subnormal number */
     E = -52; 		
     xdb.d *= ((db_number) ((double) two52)).d; 	  /* make x a normal number    */ 
   }
    
   if (xdb.i[HI] >= 0x7ff00000){
     return  x+x;				 /* Inf or Nan       */
   }
   
   /* Extract exponent and mantissa 
      Do range reduction,
      yielding to E holding the exponent and
      y the mantissa between sqrt(2)/2 and sqrt(2)
   */
   E += (xdb.i[HI]>>20)-1023;             /* extract the exponent */
   index = (xdb.i[HI] & 0x000fffff);
   xdb.i[HI] =  index | 0x3ff00000;	/* do exponent = 0 */
   index = (index + (1<<(20-L-1))) >> (20-L);
 
   /* reduce  such that sqrt(2)/2 < xdb.d < sqrt(2) */
   if (index >= MAXINDEX){ /* corresponds to xdb>sqrt(2)*/
     xdb.i[HI] -= 0x00100000; 
     E++;
   }
   y = xdb.d;
   index = index & INDEXMASK;
   // Cast integer E into double ed for multiplication later
   ed = (double) E;

   /* Read tables:
      Read one float for ri
      Read the first two doubles for -log(r_i) (out of three)

      Organization of the table:

      one struct entry per index, the struct entry containing 
      r, logih, logim and logil in this order
   */
   

   ri = argredtable[index].ri;
   // Actually we don't need the logarithm entries now
   // Move the following two lines to the eventual reconstruction
   // As long as we don't have any if in the following code, we can overlap 
   // memory access with calculations (static const?)
   logih = argredtable[index].logih;
   logim = argredtable[index].logim;

   /* Do range reduction:

      zh + zl = y * ri - 1.0 correctly

      Correctness is assured by use of Mul12 and Add12
      even if we don't force ri to have its' LSBs set to zero

      Discard zl for higher monome degrees
   */

   Mul12(&yrih, &yril, y, ri);
   th = yrih - 1.0; // Should be exact (Sterbenz) A VERIFIER !!!
   Add12Cond(zh, zl, th, yril); // Works perhaps without the condition A VERIFIER !!!

   /* Polynomial evaluation

      Use a 7 degree polynomial
      Evaluate the higher 5 terms in double precision (-7 * 3 = -21) using Horner's scheme
      Evaluate the lower 3 terms (the last is 0) in double double precision accounting also for zl
      using an ad hoc method

   */

#if defined(PROCESSOR_HAS_FMA) && !defined(AVOID_FMA)
   polyHorner = FMA(FMA(FMA(FMA(c7,zh,c6),zh,c5),zh,c4),zh,c3);
#else
   polyHorner = c3 + zh * (c4 + zh * (c5 + zh * (c6 + zh * c7)));
#endif

   Mul12(&zhSquareh, &zhSquarel, zh, zh);
   polyUpper = polyHorner * (zh * zhSquareh);
   zhSquareHalfh = zhSquareh * -0.5;
   zhSquareHalfl = zhSquarel * -0.5;
   Add12(t1h, t1l, polyUpper, -1 * (zh * zl));
   Add22(&t2h, &t2l, zh, zl, zhSquareHalfh, zhSquareHalfl);
   Add22(&ph, &pl, t2h, t2l, t1h, t1l);

   /* Reconstruction 

      Read logih and logim in the tables (already done)
      
      Compute log(x) = E * log(2) + log(1+z) - log(ri)
      i.e. log(x) = ed * (log2h + log2m) + (ph + pl) + (logih + logim) + delta

      Carry out everything in double double precision

   */
   
   /* We store log2 as log2h + log2m + log2l where log2h and log2m have 12 trailing zeros
      Multiplication of ed (double E) and log2h is thus correct
      The overall accuracy of log2h + log2m + log2l is 53 * 3 - 24 = 135 which
      is enough for the accurate phase
      The accuracy suffices also for the qucik phase: 53 * 2 - 24 = 82
      Condition is verified as log2m is smaller than log2h and both are scaled by ed
   */

   Add12(log2edh, log2edl, log2h * ed, log2m * ed);

   /* Add logih and logim to ph and pl 

      We must use conditioned Add22 as logih can move over ph

      A VERIFIER: on peut quand meme se passer de la condition
                  puisque le polynome reste par definition plus 
		  petit que la valeur tabluee sauf en 0 et la, 
		  le Add22 sans conditions passe peut etre quand meme

   */

   Add22Cond(&logTabPolyh, &logTabPolyl, logih, logim, ph, pl);

   /* Add log2edh + log2edl to logTabPolyh + logTabPolyl */

   Add22Cond(&logh, &logm, log2edh, log2edl, logTabPolyh, logTabPolyl);

   /* Rounding test and eventual return or call to the accurate function 

      Question: pourquoi les tests d arrondi des modes diriges sont factorises 
      mais pas celui du RN?

      Et: ne perd-on pas plein de cycles si on appelle une fonction pour la
      phase exacte au lieu de la recopier ou de la mettre en DEFINE?

   */

   if(E==0)
      roundcst = RDROUNDCST1;
   else
      roundcst = RDROUNDCST2;

   TEST_AND_RETURN_RD(logh, logm, roundcst);

#if DEBUG
  printf("Going for Accurate Phase for x=%1.50e\n",x);
#endif

    log_td_accurate(&logh, &logm, &logl, E, ed, index, zh, zl, logih, logim); 

    /* 
       We will suppose the following:
       - logh and logm do not overlap in any case
       - logm and logl may overlap but logl will not "overtake" logm for overlapping with logh
       - logh may be different from round-to-nearest(logh + logm + logl) 
       - logm may be equal to 0 while logl is not; 
         either one of them is supposed to be different from 0 assured by Lindemann and Lefevre
       - in no case, logh is equal to zero; in particular, the argument of log is not 1.0

       In order to overcome the fact that overlapping of logm and logl is possible or that logm is equal
       to 0, we will first add logl to logm by a double precision add. This operation will only affect the
       sign of logm if it is equal to 0 prior to the addition; it will change the value of logm only
       if there is overlapping. It will not generate an overlapping of logm and logh.
       We will now know that logm = round-to-nearest(logm + logl)
    */

    logm = logm + logl;

    /* We can now simply test the sign of logh and logm */

    if (logm > 0.0) {
      return logh;
    } else {
      loghdb.d = logh;
      if (logh > 0.0) {
	loghdb.l--;
	return loghdb.d;
      } else {
	loghdb.l++;
	return loghdb.d;
      }
    }
 } 

/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARDS ZERO			     *
 *************************************************************
 *************************************************************/
 double log_rz(double x) { 
   db_number xdb, loghdb;
   double res_hi,res_lo, y, ed, ri, logih, logim, yrih, yril, th, zh, zl;
   double polyHorner, zhSquareh, zhSquarel, polyUpper, zhSquareHalfh, zhSquareHalfl;
   double t1h, t1l, t2h, t2l, ph, pl, log2edh, log2edl, logTabPolyh, logTabPolyl, logh, logm, logl;
   double zlPlusZhSquareHalflPlusZhZl, polyUpperPluszlh, roundcst;
   int E, index;

   E=0;
   xdb.d=x;

   if (x == 1.0) return 0.0; // This the only case in which the image under log of a double is a double. 

   /* Filter cases */
   if (xdb.i[HI] < 0x00100000){        /* x < 2^(-1022)    */
     if (((xdb.i[HI] & 0x7fffffff)|xdb.i[LO])==0){
       return -1.0/0.0;     
     }                    		   /* log(+/-0) = -Inf */
     if (xdb.i[HI] < 0){ 
       return (x-x)/0;                      /* log(-x) = Nan    */
     }
     /* Subnormal number */
     E = -52; 		
     xdb.d *= ((db_number) ((double) two52)).d; 	  /* make x a normal number    */ 
   }
    
   if (xdb.i[HI] >= 0x7ff00000){
     return  x+x;				 /* Inf or Nan       */
   }
   
   /* Extract exponent and mantissa 
      Do range reduction,
      yielding to E holding the exponent and
      y the mantissa between sqrt(2)/2 and sqrt(2)
   */
   E += (xdb.i[HI]>>20)-1023;             /* extract the exponent */
   index = (xdb.i[HI] & 0x000fffff);
   xdb.i[HI] =  index | 0x3ff00000;	/* do exponent = 0 */
   index = (index + (1<<(20-L-1))) >> (20-L);
 
   /* reduce  such that sqrt(2)/2 < xdb.d < sqrt(2) */
   if (index >= MAXINDEX){ /* corresponds to xdb>sqrt(2)*/
     xdb.i[HI] -= 0x00100000; 
     E++;
   }
   y = xdb.d;
   index = index & INDEXMASK;
   // Cast integer E into double ed for multiplication later
   ed = (double) E;

   /* Read tables:
      Read one float for ri
      Read the first two doubles for -log(r_i) (out of three)

      Organization of the table:

      one struct entry per index, the struct entry containing 
      r, logih, logim and logil in this order
   */
   

   ri = argredtable[index].ri;
   // Actually we don't need the logarithm entries now
   // Move the following two lines to the eventual reconstruction
   // As long as we don't have any if in the following code, we can overlap 
   // memory access with calculations (static const?)
   logih = argredtable[index].logih;
   logim = argredtable[index].logim;

   /* Do range reduction:

      zh + zl = y * ri - 1.0 correctly

      Correctness is assured by use of Mul12 and Add12
      even if we don't force ri to have its' LSBs set to zero

      Discard zl for higher monome degrees
   */

   Mul12(&yrih, &yril, y, ri);
   th = yrih - 1.0; // Should be exact (Sterbenz) A VERIFIER !!!
   Add12Cond(zh, zl, th, yril); // Works perhaps without the condition A VERIFIER !!!

   /* Polynomial evaluation

      Use a 7 degree polynomial
      Evaluate the higher 5 terms in double precision (-7 * 3 = -21) using Horner's scheme
      Evaluate the lower 3 terms (the last is 0) in double double precision accounting also for zl
      using an ad hoc method

   */

#if defined(PROCESSOR_HAS_FMA) && !defined(AVOID_FMA)
   polyHorner = FMA(FMA(FMA(FMA(c7,zh,c6),zh,c5),zh,c4),zh,c3);
#else
   polyHorner = c3 + zh * (c4 + zh * (c5 + zh * (c6 + zh * c7)));
#endif

   Mul12(&zhSquareh, &zhSquarel, zh, zh);
   polyUpper = polyHorner * (zh * zhSquareh);
   zhSquareHalfh = zhSquareh * -0.5;
   zhSquareHalfl = zhSquarel * -0.5;
   Add12(t1h, t1l, polyUpper, -1 * (zh * zl));
   Add22(&t2h, &t2l, zh, zl, zhSquareHalfh, zhSquareHalfl);
   Add22(&ph, &pl, t2h, t2l, t1h, t1l);

   /* Reconstruction 

      Read logih and logim in the tables (already done)
      
      Compute log(x) = E * log(2) + log(1+z) - log(ri)
      i.e. log(x) = ed * (log2h + log2m) + (ph + pl) + (logih + logim) + delta

      Carry out everything in double double precision

   */
   
   /* We store log2 as log2h + log2m + log2l where log2h and log2m have 12 trailing zeros
      Multiplication of ed (double E) and log2h is thus correct
      The overall accuracy of log2h + log2m + log2l is 53 * 3 - 24 = 135 which
      is enough for the accurate phase
      The accuracy suffices also for the qucik phase: 53 * 2 - 24 = 82
      Condition is verified as log2m is smaller than log2h and both are scaled by ed
   */

   Add12(log2edh, log2edl, log2h * ed, log2m * ed);

   /* Add logih and logim to ph and pl 

      We must use conditioned Add22 as logih can move over ph

      A VERIFIER: on peut quand meme se passer de la condition
                  puisque le polynome reste par definition plus 
		  petit que la valeur tabluee sauf en 0 et la, 
		  le Add22 sans conditions passe peut etre quand meme

   */

   Add22Cond(&logTabPolyh, &logTabPolyl, logih, logim, ph, pl);

   /* Add log2edh + log2edl to logTabPolyh + logTabPolyl */

   Add22Cond(&logh, &logm, log2edh, log2edl, logTabPolyh, logTabPolyl);

   /* Rounding test and eventual return or call to the accurate function 

      Question: pourquoi les tests d arrondi des modes diriges sont factorises 
      mais pas celui du RN?

      Et: ne perd-on pas plein de cycles si on appelle une fonction pour la
      phase exacte au lieu de la recopier ou de la mettre en DEFINE?

   */

   if(E==0)
      roundcst = RDROUNDCST1;
   else
      roundcst = RDROUNDCST2;

   TEST_AND_RETURN_RZ(logh, logm, roundcst);

#if DEBUG
  printf("Going for Accurate Phase for x=%1.50e\n",x);
#endif

    log_td_accurate(&logh, &logm, &logl, E, ed, index, zh, zl, logih, logim); 

    /* 
       We will suppose the following:
       - logh and logm do not overlap in any case
       - logm and logl may overlap but logl will not "overtake" logm for overlapping with logh
       - logh may be different from round-to-nearest(logh + logm + logl) 
       - logm may be equal to 0 while logl is not; 
         either one of them is supposed to be different from 0 assured by Lindemann and Lefevre
       - in no case, logh is equal to 0; in particular, the argument of log is not 1.0

       In order to overcome the fact that overlapping of logm and logl is possible or that logm is equal
       to 0, we will first add logl to logm by a double precision add. This operation will only affect the
       sign of logm if it is equal to 0 prior to the addition; it will change the value of logm only
       if there is overlapping. It will not generate an overlapping of logm and logh.
       We will now know that logm = round-to-nearest(logm + logl)
    */

    logm = logm + logl;

    /* We can now simply test the sign of logh and logm dependingly on rounding up or down */

    if (logh > 0.0) {
      // round down
      if (logm > 0) {
	return logh;
      } else {
	loghdb.d = logh;
	loghdb.l--;
	return loghdb.d;
      }
    } else {
      // round up
      if (logm < 0) {
	return logh;
      } else {
	loghdb.d = logh;
	loghdb.l--;
	return loghdb.d;
      }
    }
 } 



