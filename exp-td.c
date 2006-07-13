/* 
 * This function computes exp, correctly rounded, 
 * using experimental techniques based on triple double arithmetics

 THIS IS EXPERIMENTAL SOFTWARE
 
 *
 * Author :  Christoph Lauter
 * christoph.lauter at ens-lyon.fr
 *

 To have it replace the crlibm exp, do:

 gcc -DHAVE_CONFIG_H -I.  -fPIC  -O2 -c exp-td.c;   mv exp-td.o exp_fast.o; make 
 
*/


#include <stdio.h>
#include <stdlib.h>
#include "crlibm.h"
#include "crlibm_private.h"
#include "triple-double.h"
#include "exp-td.h"
#ifdef BUILD_INTERVAL_FUNCTIONS
#include "interval.h"
#endif

#define AVOID_FMA 0
#define EVAL_PERF 1







void exp_td_accurate(double *polyTblh, double *polyTblm, double *polyTbll, 
		     double rh, double rm, double rl, 
		     double tbl1h, double tbl1m, double tbl1l,
		     double tbl2h, double tbl2m, double tbl2l) {
  double highPoly, highPolyMulth, highPolyMultm, highPolyMultl;
  double rhSquareh, rhSquarel, rhSquareHalfh, rhSquareHalfl;
  double rhCubeh, rhCubem, rhCubel;
  double t1h, t1l, t2h, t2l, t3h, t3l, t4h, t4l, t5, t6;
  double lowPolyh, lowPolym, lowPolyl;
  double ph, pm, pl, phnorm, pmnorm, rmlMultPh, rmlMultPl;
  double qh, ql, fullPolyh, fullPolym, fullPolyl;
  double polyWithTbl1h, polyWithTbl1m, polyWithTbl1l;
  double polyAddOneh,polyAddOnem,polyAddOnel;
  double polyWithTablesh, polyWithTablesm, polyWithTablesl;


#if EVAL_PERF
  crlibm_second_step_taken++;
#endif

#if defined(PROCESSOR_HAS_FMA) && !defined(AVOID_FMA)
  highPoly = FMA(FMA(accPolyC7,rh,accPolyC6),rh,accPolyC5);
#else
  highPoly = accPolyC5 + rh * (accPolyC6 + rh * accPolyC7);
#endif

  Mul12(&t1h,&t1l,rh,highPoly);
  Add22(&t2h,&t2l,accPolyC4h,accPolyC4l,t1h,t1l);
  Mul22(&t3h,&t3l,rh,0,t2h,t2l);
  Add22(&t4h,&t4l,accPolyC3h,accPolyC3l,t3h,t3l);

  Mul12(&rhSquareh,&rhSquarel,rh,rh);
  Mul23(&rhCubeh,&rhCubem,&rhCubel,rh,0,rhSquareh,rhSquarel);

  rhSquareHalfh = 0.5 * rhSquareh;
  rhSquareHalfl = 0.5 * rhSquarel;  

  Renormalize3(&lowPolyh,&lowPolym,&lowPolyl,rh,rhSquareHalfh,rhSquareHalfl);

  Mul233(&highPolyMulth,&highPolyMultm,&highPolyMultl,t4h,t4l,rhCubeh,rhCubem,rhCubel);

  Add33(&ph,&pm,&pl,lowPolyh,lowPolym,lowPolyl,highPolyMulth,highPolyMultm,highPolyMultl);

  Add12(phnorm,pmnorm,ph,pm);
  Mul22(&rmlMultPh,&rmlMultPl,rm,rl,phnorm,pmnorm);
  Add22(&qh,&ql,rm,rl,rmlMultPh,rmlMultPl);

  Add233Cond(&fullPolyh,&fullPolym,&fullPolyl,qh,ql,ph,pm,pl);
  Add12(polyAddOneh,t5,1,fullPolyh);
  Add12Cond(polyAddOnem,t6,t5,fullPolym);
  polyAddOnel = t6 + fullPolyl;
  Mul33(&polyWithTbl1h,&polyWithTbl1m,&polyWithTbl1l,tbl1h,tbl1m,tbl1l,polyAddOneh,polyAddOnem,polyAddOnel);
  Mul33(&polyWithTablesh,&polyWithTablesm,&polyWithTablesl,
	tbl2h,tbl2m,tbl2l,
	polyWithTbl1h,polyWithTbl1m,polyWithTbl1l);

  Renormalize3(polyTblh,polyTblm,polyTbll,polyWithTablesh,polyWithTablesm,polyWithTablesl);
}



/*************************************************************
 *************************************************************
 *               ROUNDED  TO NEAREST			     *
 *************************************************************
 *************************************************************/
double exp_rn(double x){ 
  double rh, rm, rl, tbl1h, tbl1m, tbl1l;
  double tbl2h, tbl2m, tbl2l;
  double xMultLog2InvMult2L, shiftedXMult, kd;
  double msLog2Div2LMultKh, msLog2Div2LMultKm, msLog2Div2LMultKl;
  double t1, t2, t3, t4, polyTblh, polyTblm, polyTbll;
  db_number shiftedXMultdb, twoPowerMdb, xdb, t4db, t4db2, polyTblhdb, resdb;
  int k, M, index1, index2, xIntHi, mightBeDenorm;
  double t5, t6, t7, t8, t9, t10, t11, t12, t13;
  double rhSquare, rhSquareHalf, rhC3, rhFour, monomialCube;
  double highPoly, highPolyWithSquare, monomialFour;
  double tablesh, tablesl;
  double s1, s2, s3, s4, s5;
  double res;
   
  /* Argument reduction and filtering for special cases */

  /* Compute k as a double and as an int */
  xdb.d = x;
  xMultLog2InvMult2L = x * log2InvMult2L;
  shiftedXMult = xMultLog2InvMult2L + shiftConst;
  kd = shiftedXMult - shiftConst;
  shiftedXMultdb.d = shiftedXMult;
  
  /* Special cases tests */
  xIntHi = xdb.i[HI];
  mightBeDenorm = 0;
  /* Test if argument is a denormal or zero */
  if ((xIntHi & 0x7ff00000) == 0) {
    /* We are in the RN case, return 1.0 in all cases */
    return 1.0;
  }
 
  /* Test if argument is greater than approx. 709 in magnitude */
  if ((xIntHi & 0x7fffffff) >= OVRUDRFLWSMPLBOUND) {
    /* If we are here, the result might be overflowed, underflowed, inf, or NaN */

    /* Test if +/- Inf or NaN */
    if ((xIntHi & 0x7fffffff) >= 0x7ff00000) {
      /* Either NaN or Inf in this case since exponent is maximal */

      /* Test if NaN: mantissa is not 0 */
      if (((xIntHi & 0x000fffff) | xdb.i[LO]) != 0) {
	/* x = NaN, return NaN */
	return x + x;
      } else {
	/* +/- Inf */

	/* Test sign */
	if ((xIntHi & 0x80000000)==0) 
	  /* x = +Inf, return +Inf */
	  return x;
	else
	  /* x = -Inf, return 0 */
	  return 0;
      } /* End which in NaN, Inf */
    } /* End NaN or Inf ? */
    
    /* If we are here, we might be overflowed, denormalized or underflowed in the result 
       but there is no special case (NaN, Inf) left */

    /* Test if actually overflowed */
    if (x > OVRFLWBOUND) {
      /* We are actually overflowed in the result */
      return LARGEST * LARGEST;
    }

    /* Test if surely underflowed */
    if (x <= UNDERFLWBOUND) {
      /* We are actually sure to be underflowed and not denormalized any more 
	 So we return 0 and raise the inexact flag */
      return SMALLEST * SMALLEST;
    }
       
    /* Test if possibly denormalized */
    if (x <= DENORMBOUND) {
      /* We know now that we are not sure to be normalized in the result
	 We just set an internal flag for a further test 
      */
      mightBeDenorm = 1;
    }
  } /* End might be a special case */

  /* If we are here, we are sure to be neither +/- Inf nor NaN nor overflowed nor denormalized in the argument
     but we might be denormalized in the result 

     We continue the argument reduction for the quick phase and table reads for both phases
  */

  Mul12(&s1,&s2,msLog2Div2Lh,kd);
  s3 = kd * msLog2Div2Lm;
  s4 = s2 + s3; 
  s5 = x + s1;
  Add12Cond(rh,rm,s5,s4);

  k = shiftedXMultdb.i[LO];
  M = k >> L;
  index1 = k & INDEXMASK1;
  index2 = (k & INDEXMASK2) >> LHALF;

  /* Table reads */
  tbl1h = twoPowerIndex1[index1].hi;
  tbl1m = twoPowerIndex1[index1].mi;
  tbl2h = twoPowerIndex2[index2].hi;
  tbl2m = twoPowerIndex2[index2].mi;

  /* Test now if it is sure to launch the quick phase because no denormalized result is possible */
  if (mightBeDenorm == 1) {
    /* The result might be denormalized, we launch the accurate phase in all cases */

    /* Rest of argument reduction for accurate phase */

    Mul133(&msLog2Div2LMultKh,&msLog2Div2LMultKm,&msLog2Div2LMultKl,kd,msLog2Div2Lh,msLog2Div2Lm,msLog2Div2Ll);
    t1 = x + msLog2Div2LMultKh;
    Add12Cond(rh,t2,t1,msLog2Div2LMultKm);
    Add12Cond(rm,rl,t2,msLog2Div2LMultKl);

    /* Table reads for accurate phase */
     tbl1l = twoPowerIndex1[index1].lo;
     tbl2l = twoPowerIndex2[index2].lo;

     /* Call accurate phase */
     exp_td_accurate(&polyTblh, &polyTblm, &polyTbll, rh, rm, rl, tbl1h, tbl1m, tbl1l, tbl2h, tbl2m, tbl2l); 

     /* Final rounding and multiplication with 2^M 

        We first multiply the highest significant byte by 2^M in two steps
	and adjust it then depending on the lower significant parts.

	We cannot multiply directly by 2^M since M is less than -1022.
	We first multiply by 2^(-1000) and then by 2^(M+1000).

     */
     
     t3 = polyTblh * twoPowerM1000;

     /* Form now twoPowerM with adjusted M */
     twoPowerMdb.i[LO] = 0;
     twoPowerMdb.i[HI] = (M + 2023) << 20;


     /* Multiply with the rest of M, the result will be denormalized */
     t4 = t3 * twoPowerMdb.d;

     /* For x86, force the compiler to pass through memory for having the right rounding */

     t4db.d = t4;   /* Do not #if-ify this line, we need the copy */
#if defined(CRLIBM_TYPECPU_AMD64) || defined(CRLIBM_TYPECPU_X86) 
     t4db2.i[HI] = t4db.i[HI];
     t4db2.i[LO] = t4db.i[LO];
     t4 = t4db2.d;
#endif

     /* Remultiply by 2^(-M) for manipulating the rounding error and the lower significant parts */
     M *= -1;
     twoPowerMdb.i[LO] = 0;
     twoPowerMdb.i[HI] = (M + 23) << 20;
     t5 = t4 * twoPowerMdb.d;
     t6 = t5 * twoPower1000;
     t7 = polyTblh - t6;
     
     /* The rounding decision is made at 1/2 ulp of a denormal, i.e. at 2^(-1075)
	We construct this number and by comparing with it we get to know 
	whether we are in a difficult rounding case or not. If not we just return 
	the known result. Otherwise we continue with further tests.
     */

     twoPowerMdb.i[LO] = 0;
     twoPowerMdb.i[HI] = (M - 52) << 20;

     if (ABS(t7) != twoPowerMdb.d) return t4;

     /* If we are here, we are in a difficult rounding case */
     
     /* We have to adjust the result iff the sign of the error on 
	rounding 2^M * polyTblh (which must be an ulp of a denormal) 
	and polyTblm +arith polyTbll is the same which means that 
	the error made was greater than an ulp of an denormal.
     */

     polyTblm = polyTblm + polyTbll;

     if (t7 > 0.0) {
       if (polyTblm > 0.0) {
	 t4db.l++;
	 return t4db.d;
       } else return t4;
     } else {
       if (polyTblm < 0.0) {
	 t4db.l--;
	 return t4db.d;
       } else return t4;
     }
  } /* End accurate phase launched as there might be a denormalized result */

  /* No more underflow nor denormal is possible. There may be the case where
     M is 1024 and the value 2^M is to be multiplied may be less than 1
     So the final result will be normalized and representable by the multiplication must be 
     made in 2 steps
  */

  /* Quick phase starts here */

  rhSquare = rh * rh;
  rhC3 = c3 * rh;

  rhSquareHalf = 0.5 * rhSquare;
  monomialCube = rhC3 * rhSquare;
  rhFour = rhSquare * rhSquare;

  monomialFour = c4 * rhFour;
  
  highPoly = monomialCube + monomialFour;

  highPolyWithSquare = rhSquareHalf + highPoly;

  Mul22(&tablesh,&tablesl,tbl1h,tbl1m,tbl2h,tbl2m);

  t8 = rm + highPolyWithSquare;
  t9 = rh + t8;

  t10 = tablesh * t9;
  
  Add12(t11,t12,tablesh,t10);
  t13 = t12 + tablesl;
  Add12(polyTblh,polyTblm,t11,t13);
  
  /* Rounding test 
     Since we know that the result of the final multiplication with 2^M 
     will always be representable, we can do the rounding test on the 
     factors and multiply only the final result.
     We implement the multiplication in integer computations to overcome
     the problem of the non-representability of 2^1024 if M = 1024
  */

  if(polyTblh == (polyTblh + (polyTblm * ROUNDCST))) {
    polyTblhdb.d = polyTblh;
    polyTblhdb.i[HI] += M << 20;
    return polyTblhdb.d;
  } else 
    {
      /* Rest of argument reduction for accurate phase */

      Mul133(&msLog2Div2LMultKh,&msLog2Div2LMultKm,&msLog2Div2LMultKl,kd,msLog2Div2Lh,msLog2Div2Lm,msLog2Div2Ll);
      t1 = x + msLog2Div2LMultKh;
      Add12Cond(rh,t2,t1,msLog2Div2LMultKm);
      Add12Cond(rm,rl,t2,msLog2Div2LMultKl);

      /* Table reads for accurate phase */
      tbl1l = twoPowerIndex1[index1].lo;
      tbl2l = twoPowerIndex2[index2].lo;
      
      /* Call accurate phase */
      exp_td_accurate(&polyTblh, &polyTblm, &polyTbll, rh, rm, rl, tbl1h, tbl1m, tbl1l, tbl2h, tbl2m, tbl2l); 

      /* Since the final multiplication is exact, we can do the final rounding before multiplying
	 We overcome this way also the cases where the final result is not underflowed whereas the
	 lower parts of the intermediate final result are.
      */
      
      RoundToNearest3(&res,polyTblh,polyTblm,polyTbll);

      /* Final multiplication with 2^M 
	 We implement the multiplication in integer computations to overcome
	 the problem of the non-representability of 2^1024 if M = 1024
      */

      resdb.d = res;
      resdb.i[HI] += M << 20;
      return resdb.d;
    } /* Accurate phase launched after rounding test*/
}


/*************************************************************
 *************************************************************
 *               ROUNDED  UPWARDS			     *
 *************************************************************
 *************************************************************/
double exp_ru(double x) { 
  double rh, rm, rl, tbl1h, tbl1m, tbl1l;
  double tbl2h, tbl2m, tbl2l;
  double xMultLog2InvMult2L, shiftedXMult, kd;
  double msLog2Div2LMultKh, msLog2Div2LMultKm, msLog2Div2LMultKl;
  double t1, t2, t3, t4, polyTblh, polyTblm, polyTbll;
  db_number shiftedXMultdb, twoPowerMdb, xdb, t4db, t4db2, resdb;
  int k, M, index1, index2, xIntHi, mightBeDenorm, roundable;
  double t5, t6, t7, t8, t9, t10, t11, t12, t13;
  double rhSquare, rhSquareHalf, rhC3, rhFour, monomialCube;
  double highPoly, highPolyWithSquare, monomialFour;
  double tablesh, tablesl;
  double s1, s2, s3, s4, s5;
  double res;
 
  /* Argument reduction and filtering for special cases */

  /* Compute k as a double and as an int */
  xdb.d = x;
  xMultLog2InvMult2L = x * log2InvMult2L;
  shiftedXMult = xMultLog2InvMult2L + shiftConst;
  kd = shiftedXMult - shiftConst;
  shiftedXMultdb.d = shiftedXMult;
  
  /* Special cases tests */
  xIntHi = xdb.i[HI];
  mightBeDenorm = 0;
  /* Test if argument is a denormal or zero */
  if ((xIntHi & 0x7ff00000) == 0) {
    /* If the argument is exactly zero, we just return 1.0
       which is the mathematical image of the function
    */
    if (x == 0.0) return 1.0;

    /* If the argument is a negative denormal, we 
       must return 1.0 and raise the inexact flag.
    */

    if (x < 0.0) return 1.0 + SMALLEST;

    /* Otherwise, we return 1.0 + 1ulp since 
       exp(greatest denorm) < 1.0 + 1ulp
       We must do the addition dynamically for
       raising the inexact flag.
    */

    return 1.0 + twoM52;
  }
 
  /* Test if argument is greater than approx. 709 in magnitude */
  if ((xIntHi & 0x7fffffff) >= OVRUDRFLWSMPLBOUND) {
    /* If we are here, the result might be overflowed, underflowed, inf, or NaN */

    /* Test if +/- Inf or NaN */
    if ((xIntHi & 0x7fffffff) >= 0x7ff00000) {
      /* Either NaN or Inf in this case since exponent is maximal */

      /* Test if NaN: mantissa is not 0 */
      if (((xIntHi & 0x000fffff) | xdb.i[LO]) != 0) {
	/* x = NaN, return NaN */
	return x + x;
      } else {
	/* +/- Inf */

	/* Test sign */
	if ((xIntHi & 0x80000000)==0) 
	  /* x = +Inf, return +Inf */
	  return x;
	else
	  /* x = -Inf, return 0 (even in RU!) */
	  return 0;
      } /* End which in NaN, Inf */
    } /* End NaN or Inf ? */
    
    /* If we are here, we might be overflowed, denormalized or underflowed in the result 
       but there is no special case (NaN, Inf) left */

    /* Test if actually overflowed */
    if (x > OVRFLWBOUND) {
      /* We are actually overflowed in the result */
      return LARGEST * LARGEST;
    }

    /* Test if surely underflowed */
    if (x <= UNDERFLWBOUND) {
      /* We are actually sure to be underflowed and not denormalized any more 
	 (at least where computing makes sense); since we are in the round 
	 upwards case, we return the smallest denormal possible.
      */
      return SMALLEST;
    }
       
    /* Test if possibly denormalized */
    if (x <= DENORMBOUND) {
      /* We know now that we are not sure to be normalized in the result
	 We just set an internal flag for a further test 
      */
      mightBeDenorm = 1;
    }
  } /* End might be a special case */

  /* If we are here, we are sure to be neither +/- Inf nor NaN nor overflowed nor denormalized in the argument
     but we might be denormalized in the result 

     We continue the argument reduction for the quick phase and table reads for both phases
  */

  Mul12(&s1,&s2,msLog2Div2Lh,kd);
  s3 = kd * msLog2Div2Lm;
  s4 = s2 + s3; 
  s5 = x + s1;
  Add12Cond(rh,rm,s5,s4);

  k = shiftedXMultdb.i[LO];
  M = k >> L;
  index1 = k & INDEXMASK1;
  index2 = (k & INDEXMASK2) >> LHALF;

  /* Table reads */
  tbl1h = twoPowerIndex1[index1].hi;
  tbl1m = twoPowerIndex1[index1].mi;
  tbl2h = twoPowerIndex2[index2].hi;
  tbl2m = twoPowerIndex2[index2].mi;

  /* Test now if it is sure to launch the quick phase because no denormalized result is possible */
  if (mightBeDenorm == 1) {
    /* The result might be denormalized, we launch the accurate phase in all cases */

    /* Rest of argument reduction for accurate phase */

    Mul133(&msLog2Div2LMultKh,&msLog2Div2LMultKm,&msLog2Div2LMultKl,kd,msLog2Div2Lh,msLog2Div2Lm,msLog2Div2Ll);
    t1 = x + msLog2Div2LMultKh;
    Add12Cond(rh,t2,t1,msLog2Div2LMultKm);
    Add12Cond(rm,rl,t2,msLog2Div2LMultKl);

    /* Table reads for accurate phase */
     tbl1l = twoPowerIndex1[index1].lo;
     tbl2l = twoPowerIndex2[index2].lo;

     /* Call accurate phase */
     exp_td_accurate(&polyTblh, &polyTblm, &polyTbll, rh, rm, rl, tbl1h, tbl1m, tbl1l, tbl2h, tbl2m, tbl2l); 

     /* Final rounding and multiplication with 2^M 

        We first multiply the highest significant byte by 2^M in two steps
	and adjust it then depending on the lower significant parts.

	We cannot multiply directly by 2^M since M is less than -1022.
	We first multiply by 2^(-1000) and then by 2^(M+1000).

     */
     
     t3 = polyTblh * twoPowerM1000;

     /* Form now twoPowerM with adjusted M */
     twoPowerMdb.i[LO] = 0;
     twoPowerMdb.i[HI] = (M + 2023) << 20;


     /* Multiply with the rest of M, the result will be denormalized */
     t4 = t3 * twoPowerMdb.d;

     /* For x86, force the compiler to pass through memory for having the right rounding */

     t4db.d = t4;   /* Do not #if-ify this line, we need the copy */
#if defined(CRLIBM_TYPECPU_AMD64) || defined(CRLIBM_TYPECPU_X86) 
     t4db2.i[HI] = t4db.i[HI];
     t4db2.i[LO] = t4db.i[LO];
     t4 = t4db2.d;
#endif


     /* Remultiply by 2^(-M) for manipulating the rounding error and the lower significant parts */
     M *= -1;
     twoPowerMdb.i[LO] = 0;
     twoPowerMdb.i[HI] = (M + 23) << 20;
     t5 = t4 * twoPowerMdb.d;
     t6 = t5 * twoPower1000;
     t7 = polyTblh - t6;

     /* The rounding can be decided using the sign of the arithmetical sum of the
	round-to-nearest-error (i.e. t7) and the lower part(s) of the final result.
	We add first the lower parts and add the result to the error in t7. We have to 
	keep in mind that everything is scaled by 2^(-M).
	t8 can never be exactly 0 since we filter out the cases where the image of the 
	function is algebraic and the implementation is exacter than the TMD worst case.
     */
 
     polyTblm = polyTblm + polyTbll;
     t8 = t7 + polyTblm;

     /* Since we are rounding upwards, the round-to-nearest-rounding result in t4 is 
	equal to the final result if the rounding error (i.e. the error plus the lower parts)
	is negative, i.e. if the rounding-to-nearest was upwards.
     */
     
     if (t8 < 0.0) return t4;

     /* If we are here, we must adjust the final result by +1ulp 
	Relying on the fact that the exponential is always positive, we can simplify this
	adjustment 
     */

     t4db.l++;
     return t4db.d;
  } /* End accurate phase launched as there might be a denormalized result */

  /* No more underflow nor denormal is possible. There may be the case where
     M is 1024 and the value 2^M is to be multiplied may be less than 1
     So the final result will be normalized and representable by the multiplication must be 
     made in 2 steps
  */

  /* Quick phase starts here */

  rhSquare = rh * rh;
  rhC3 = c3 * rh;

  rhSquareHalf = 0.5 * rhSquare;
  monomialCube = rhC3 * rhSquare;
  rhFour = rhSquare * rhSquare;

  monomialFour = c4 * rhFour;
  
  highPoly = monomialCube + monomialFour;

  highPolyWithSquare = rhSquareHalf + highPoly;

  Mul22(&tablesh,&tablesl,tbl1h,tbl1m,tbl2h,tbl2m);

  t8 = rm + highPolyWithSquare;
  t9 = rh + t8;

  t10 = tablesh * t9;
  
  Add12(t11,t12,tablesh,t10);
  t13 = t12 + tablesl;
  Add12(polyTblh,polyTblm,t11,t13);
  
  /* Rounding test 
     Since we know that the result of the final multiplication with 2^M 
     will always be representable, we can do the rounding test on the 
     factors and multiply only the final result.
     We implement the multiplication in integer computations to overcome
     the problem of the non-representability of 2^1024 if M = 1024
  */

  TEST_AND_COPY_RU(roundable,res,polyTblh,polyTblm,RDROUNDCST);

  if (roundable) {
    resdb.d = res;
    resdb.i[HI] += M << 20;
    return resdb.d;
  } else 
    {
      /* Rest of argument reduction for accurate phase */

      Mul133(&msLog2Div2LMultKh,&msLog2Div2LMultKm,&msLog2Div2LMultKl,kd,msLog2Div2Lh,msLog2Div2Lm,msLog2Div2Ll);
      t1 = x + msLog2Div2LMultKh;
      Add12Cond(rh,t2,t1,msLog2Div2LMultKm);
      Add12Cond(rm,rl,t2,msLog2Div2LMultKl);

      /* Table reads for accurate phase */
      tbl1l = twoPowerIndex1[index1].lo;
      tbl2l = twoPowerIndex2[index2].lo;
      
      /* Call accurate phase */
      exp_td_accurate(&polyTblh, &polyTblm, &polyTbll, rh, rm, rl, tbl1h, tbl1m, tbl1l, tbl2h, tbl2m, tbl2l); 

      /* Since the final multiplication is exact, we can do the final rounding before multiplying
	 We overcome this way also the cases where the final result is not underflowed whereas the
	 lower parts of the intermediate final result are.
      */
      
      RoundUpwards3(&res,polyTblh,polyTblm,polyTbll);

      /* Final multiplication with 2^M 
	 We implement the multiplication in integer computations to overcome
	 the problem of the non-representability of 2^1024 if M = 1024
      */

      resdb.d = res;
      resdb.i[HI] += M << 20;
      return resdb.d;
    } /* Accurate phase launched after rounding test*/
} 


/*************************************************************
 *************************************************************
 *               ROUNDED  DOWNWARDS			     *
 *************************************************************
 *************************************************************/
double exp_rd(double x) { 
  double rh, rm, rl, tbl1h, tbl1m, tbl1l;
  double tbl2h, tbl2m, tbl2l;
  double xMultLog2InvMult2L, shiftedXMult, kd;
  double msLog2Div2LMultKh, msLog2Div2LMultKm, msLog2Div2LMultKl;
  double t1, t2, t3, t4, polyTblh, polyTblm, polyTbll;
  db_number shiftedXMultdb, twoPowerMdb, xdb, t4db, t4db2, resdb;
  int k, M, index1, index2, xIntHi, mightBeDenorm, roundable;
  double t5, t6, t7, t8, t9, t10, t11, t12, t13;
  double rhSquare, rhSquareHalf, rhC3, rhFour, monomialCube;
  double highPoly, highPolyWithSquare, monomialFour;
  double tablesh, tablesl;
  double s1, s2, s3, s4, s5;
  double res;
 
  /* Argument reduction and filtering for special cases */

  /* Compute k as a double and as an int */
  xdb.d = x;
  xMultLog2InvMult2L = x * log2InvMult2L;
  shiftedXMult = xMultLog2InvMult2L + shiftConst;
  kd = shiftedXMult - shiftConst;
  shiftedXMultdb.d = shiftedXMult;
  
  /* Special cases tests */
  xIntHi = xdb.i[HI];
  mightBeDenorm = 0;
  /* Test if argument is a denormal or zero */
  if ((xIntHi & 0x7ff00000) == 0) {
    /* If the argument is exactly zero, we just return 1.0
       which is the mathematical image of the function
    */
    if (x == 0.0) return 1.0;

    /* If the argument is a positive denormal, we 
       must return 1.0 and raise the inexact flag.
    */
    
    if (x > 0.0) return 1.0 + SMALLEST;

    /* Otherwise, we return 1.0 - 1ulp since 
       exp(-greatest denorm) > 1.0 - 1ulp
       We must do the addition dynamically for
       raising the inexact flag.
    */
    
    return 1.0 + mTwoM53;

  }
 
  /* Test if argument is greater than approx. 709 in magnitude */
  if ((xIntHi & 0x7fffffff) >= OVRUDRFLWSMPLBOUND) {
    /* If we are here, the result might be overflowed, underflowed, inf, or NaN */

    /* Test if +/- Inf or NaN */
    if ((xIntHi & 0x7fffffff) >= 0x7ff00000) {
      /* Either NaN or Inf in this case since exponent is maximal */

      /* Test if NaN: mantissa is not 0 */
      if (((xIntHi & 0x000fffff) | xdb.i[LO]) != 0) {
	/* x = NaN, return NaN */
	return x + x;
      } else {
	/* +/- Inf */

	/* Test sign */
	if ((xIntHi & 0x80000000)==0) 
	  /* x = +Inf, return +Inf */
	  return x;
	else
	  /* x = -Inf, return 0 */
	  return 0;
      } /* End which in NaN, Inf */
    } /* End NaN or Inf ? */
    
    /* If we are here, we might be overflowed, denormalized or underflowed in the result 
       but there is no special case (NaN, Inf) left */

    /* Test if actually overflowed */
    if (x > OVRFLWBOUND) {
      /* We would be overflowed but as we are rounding downwards
	 the nearest number lesser than the exact result is the greatest 
	 normal. In any case, we must raise the inexact flag.
      */
      return LARGEST * (1.0 + SMALLEST);
    }

    /* Test if surely underflowed */
    if (x <= UNDERFLWBOUND) {
      /* We are actually sure to be underflowed and not denormalized any more 
	 (at least where computing makes sense); since we are in the round 
	 upwards case, we return the smallest denormal possible.
      */
      return SMALLEST * SMALLEST;
    }
       
    /* Test if possibly denormalized */
    if (x <= DENORMBOUND) {
      /* We know now that we are not sure to be normalized in the result
	 We just set an internal flag for a further test 
      */
      mightBeDenorm = 1;
    }
  } /* End might be a special case */

  /* If we are here, we are sure to be neither +/- Inf nor NaN nor overflowed nor denormalized in the argument
     but we might be denormalized in the result 

     We continue the argument reduction for the quick phase and table reads for both phases
  */

  Mul12(&s1,&s2,msLog2Div2Lh,kd);
  s3 = kd * msLog2Div2Lm;
  s4 = s2 + s3; 
  s5 = x + s1;
  Add12Cond(rh,rm,s5,s4);

  k = shiftedXMultdb.i[LO];
  M = k >> L;
  index1 = k & INDEXMASK1;
  index2 = (k & INDEXMASK2) >> LHALF;

  /* Table reads */
  tbl1h = twoPowerIndex1[index1].hi;
  tbl1m = twoPowerIndex1[index1].mi;
  tbl2h = twoPowerIndex2[index2].hi;
  tbl2m = twoPowerIndex2[index2].mi;

  /* Test now if it is sure to launch the quick phase because no denormalized result is possible */
  if (mightBeDenorm == 1) {
    /* The result might be denormalized, we launch the accurate phase in all cases */

    /* Rest of argument reduction for accurate phase */

    Mul133(&msLog2Div2LMultKh,&msLog2Div2LMultKm,&msLog2Div2LMultKl,kd,msLog2Div2Lh,msLog2Div2Lm,msLog2Div2Ll);
    t1 = x + msLog2Div2LMultKh;
    Add12Cond(rh,t2,t1,msLog2Div2LMultKm);
    Add12Cond(rm,rl,t2,msLog2Div2LMultKl);

    /* Table reads for accurate phase */
     tbl1l = twoPowerIndex1[index1].lo;
     tbl2l = twoPowerIndex2[index2].lo;

     /* Call accurate phase */
     exp_td_accurate(&polyTblh, &polyTblm, &polyTbll, rh, rm, rl, tbl1h, tbl1m, tbl1l, tbl2h, tbl2m, tbl2l); 

     /* Final rounding and multiplication with 2^M 

        We first multiply the highest significant byte by 2^M in two steps
	and adjust it then depending on the lower significant parts.

	We cannot multiply directly by 2^M since M is less than -1022.
	We first multiply by 2^(-1000) and then by 2^(M+1000).

     */
     
     t3 = polyTblh * twoPowerM1000;

     /* Form now twoPowerM with adjusted M */
     twoPowerMdb.i[LO] = 0;
     twoPowerMdb.i[HI] = (M + 2023) << 20;


     /* Multiply with the rest of M, the result will be denormalized */
     t4 = t3 * twoPowerMdb.d;

     /* For x86, force the compiler to pass through memory for having the right rounding */

     t4db.d = t4;   /* Do not #if-ify this line, we need the copy */
#if defined(CRLIBM_TYPECPU_AMD64) || defined(CRLIBM_TYPECPU_X86) 
     t4db2.i[HI] = t4db.i[HI];
     t4db2.i[LO] = t4db.i[LO];
     t4 = t4db2.d;
#endif

     /* Remultiply by 2^(-M) for manipulating the rounding error and the lower significant parts */
     M *= -1;
     twoPowerMdb.i[LO] = 0;
     twoPowerMdb.i[HI] = (M + 23) << 20;
     t5 = t4 * twoPowerMdb.d;
     t6 = t5 * twoPower1000;
     t7 = polyTblh - t6;

     /* The rounding can be decided using the sign of the arithmetical sum of the
	round-to-nearest-error (i.e. t7) and the lower part(s) of the final result.
	We add first the lower parts and add the result to the error in t7. We have to 
	keep in mind that everything is scaled by 2^(-M).
	t8 can never be exactly 0 since we filter out the cases where the image of the 
	function is algebraic and the implementation is exacter than the TMD worst case.
     */
 
     polyTblm = polyTblm + polyTbll;
     t8 = t7 + polyTblm;

     /* Since we are rounding downwards, the round-to-nearest-rounding result in t4 is 
	equal to the final result if the rounding error (i.e. the error plus the lower parts)
	is positive, i.e. if the rounding-to-nearest was downwards.
     */
     
     if (t8 > 0.0) return t4;

     /* If we are here, we must adjust the final result by +1ulp 
	Relying on the fact that the exponential is always positive, we can simplify this
	adjustment 
     */

     t4db.l--;
     return t4db.d;
  } /* End accurate phase launched as there might be a denormalized result */

  /* No more underflow nor denormal is possible. There may be the case where
     M is 1024 and the value 2^M is to be multiplied may be less than 1
     So the final result will be normalized and representable by the multiplication must be 
     made in 2 steps
  */

  /* Quick phase starts here */

  rhSquare = rh * rh;
  rhC3 = c3 * rh;

  rhSquareHalf = 0.5 * rhSquare;
  monomialCube = rhC3 * rhSquare;
  rhFour = rhSquare * rhSquare;

  monomialFour = c4 * rhFour;
  
  highPoly = monomialCube + monomialFour;

  highPolyWithSquare = rhSquareHalf + highPoly;

  Mul22(&tablesh,&tablesl,tbl1h,tbl1m,tbl2h,tbl2m);

  t8 = rm + highPolyWithSquare;
  t9 = rh + t8;

  t10 = tablesh * t9;
  
  Add12(t11,t12,tablesh,t10);
  t13 = t12 + tablesl;
  Add12(polyTblh,polyTblm,t11,t13);
  
  /* Rounding test 
     Since we know that the result of the final multiplication with 2^M 
     will always be representable, we can do the rounding test on the 
     factors and multiply only the final result.
     We implement the multiplication in integer computations to overcome
     the problem of the non-representability of 2^1024 if M = 1024
  */

  TEST_AND_COPY_RD(roundable,res,polyTblh,polyTblm,RDROUNDCST);

  if (roundable) {
    resdb.d = res;
    resdb.i[HI] += M << 20;
    return resdb.d;
  } else {      
      /* Rest of argument reduction for accurate phase */

      Mul133(&msLog2Div2LMultKh,&msLog2Div2LMultKm,&msLog2Div2LMultKl,kd,msLog2Div2Lh,msLog2Div2Lm,msLog2Div2Ll);
      t1 = x + msLog2Div2LMultKh;
      Add12Cond(rh,t2,t1,msLog2Div2LMultKm);
      Add12Cond(rm,rl,t2,msLog2Div2LMultKl);

      /* Table reads for accurate phase */
      tbl1l = twoPowerIndex1[index1].lo;
      tbl2l = twoPowerIndex2[index2].lo;
      
      /* Call accurate phase */
      exp_td_accurate(&polyTblh, &polyTblm, &polyTbll, rh, rm, rl, tbl1h, tbl1m, tbl1l, tbl2h, tbl2m, tbl2l); 

      /* Since the final multiplication is exact, we can do the final rounding before multiplying
	 We overcome this way also the cases where the final result is not underflowed whereas the
	 lower parts of the intermediate final result are.
      */
      
      RoundDownwards3(&res,polyTblh,polyTblm,polyTbll);

      /* Final multiplication with 2^M 
	 We implement the multiplication in integer computations to overcome
	 the problem of the non-representability of 2^1024 if M = 1024
      */

      resdb.d = res;
      resdb.i[HI] += M << 20;
      return resdb.d;
    } /* Accurate phase launched after rounding test*/
} 
 

#ifdef BUILD_INTERVAL_FUNCTIONS
interval j_exp(interval x)
{
  interval res;
  double x_rd, x_ru;
  double rh_ru, rm_ru, rl_ru, tbl1h_ru, tbl1m_ru, tbl1l_ru;
  double tbl2h_ru, tbl2m_ru, tbl2l_ru;
  double xMultLog2InvMult2L_ru, shiftedXMult_ru, kd_ru;
  double msLog2Div2LMultKh_ru, msLog2Div2LMultKm_ru, msLog2Div2LMultKl_ru;
  double t1_ru, t2_ru, polyTblh_ru, polyTblm_ru, polyTbll_ru;
  db_number shiftedXMultdb_ru, xdb_ru, resdb_ru;
  int k_ru, M_ru, index1_ru, index2_ru, xIntHi_ru, mightBeDenorm_ru, roundable;
  double t8_ru, t9_ru, t10_ru, t11_ru, t12_ru, t13_ru;
  double rhSquare_ru, rhSquareHalf_ru, rhC3_ru, rhFour_ru, monomialCube_ru;
  double highPoly_ru, highPolyWithSquare_ru, monomialFour_ru;
  double tablesh_ru, tablesl_ru;
  double s1_ru, s2_ru, s3_ru, s4_ru, s5_ru;
  double res_ru;

  double rh_rd, rm_rd, rl_rd, tbl1h_rd, tbl1m_rd, tbl1l_rd;
  double tbl2h_rd, tbl2m_rd, tbl2l_rd;
  double xMultLog2InvMult2L_rd, shiftedXMult_rd, kd_rd;
  double msLog2Div2LMultKh_rd, msLog2Div2LMultKm_rd, msLog2Div2LMultKl_rd;
  double t1_rd, t2_rd, polyTblh_rd, polyTblm_rd, polyTbll_rd;
  db_number shiftedXMultdb_rd, xdb_rd, resdb_rd;
  int k_rd, M_rd, index1_rd, index2_rd, xIntHi_rd, mightBeDenorm_rd;
  double t8_rd, t9_rd, t10_rd, t11_rd, t12_rd, t13_rd;
  double rhSquare_rd, rhSquareHalf_rd, rhC3_rd, rhFour_rd, monomialCube_rd;
  double highPoly_rd, highPolyWithSquare_rd, monomialFour_rd;
  double tablesh_rd, tablesl_rd;
  double s1_rd, s2_rd, s3_rd, s4_rd, s5_rd;
  double res_rd;

  double res_simple_rd, res_simple_ru;
  int infDone=0; int supDone=0;

  x_rd=LOW(x);
  x_ru=UP(x);

  /* Argument reduction and filtering for special cases */

  /* Compute k as a double and as an int */
  xdb_ru.d = x_ru;
  xdb_rd.d = x_rd;
  xMultLog2InvMult2L_ru = x_ru * log2InvMult2L;
  xMultLog2InvMult2L_rd = x_rd * log2InvMult2L;
  shiftedXMult_ru = xMultLog2InvMult2L_ru + shiftConst;
  shiftedXMult_rd = xMultLog2InvMult2L_rd + shiftConst;
  kd_ru = shiftedXMult_ru - shiftConst;
  kd_rd = shiftedXMult_rd - shiftConst;
  shiftedXMultdb_ru.d = shiftedXMult_ru;
  shiftedXMultdb_rd.d = shiftedXMult_rd;


  /* Special cases tests */
  xIntHi_ru = xdb_ru.i[HI];
  mightBeDenorm_ru = 0;

  /* Special cases tests */
  xIntHi_rd = xdb_rd.i[HI];
  mightBeDenorm_rd = 0;

  if ( __builtin_expect(
       ((xIntHi_ru & 0x7ff00000) == 0)
    || (((xIntHi_ru & 0x7ff00000) == 0)  && (x_ru == 0.0)) 
    || (((xIntHi_ru & 0x7ff00000) == 0)  && (x_ru < 0.0))
    || (((xIntHi_ru & 0x7fffffff) >= OVRUDRFLWSMPLBOUND) && ((xIntHi_ru & 0x7fffffff) >= 0x7ff00000))
    || (((xIntHi_ru & 0x7fffffff) >= OVRUDRFLWSMPLBOUND) && ((xIntHi_ru & 0x7fffffff) >= 0x7ff00000) && (((xIntHi_ru & 0x000fffff) | xdb_ru.i[LO]) != 0))
    || (((xIntHi_ru & 0x7fffffff) >= OVRUDRFLWSMPLBOUND) && ((xIntHi_ru & 0x7fffffff) >= 0x7ff00000) && ((xIntHi_ru & 0x80000000)==0))
    || (((xIntHi_ru & 0x7fffffff) >= OVRUDRFLWSMPLBOUND) && (x_ru > OVRFLWBOUND))
    || (((xIntHi_ru & 0x7fffffff) >= OVRUDRFLWSMPLBOUND) && (x_ru <= UNDERFLWBOUND))
    || (((xIntHi_ru & 0x7fffffff) >= OVRUDRFLWSMPLBOUND) && (x_ru <= DENORMBOUND))
    || ((xIntHi_rd & 0x7ff00000) == 0)
    || (((xIntHi_rd & 0x7ff00000) == 0) && (x_rd == 0.0))
    || (((xIntHi_rd & 0x7ff00000) == 0) && (x_rd > 0.0))
    || (((xIntHi_rd & 0x7fffffff) >= OVRUDRFLWSMPLBOUND) && ((xIntHi_rd & 0x7fffffff) >= 0x7ff00000))
    || (((xIntHi_rd & 0x7fffffff) >= OVRUDRFLWSMPLBOUND) && ((xIntHi_rd & 0x7fffffff) >= 0x7ff00000) && (((xIntHi_rd & 0x000fffff) | xdb_rd.i[LO]) != 0))
    || (((xIntHi_rd & 0x7fffffff) >= OVRUDRFLWSMPLBOUND) && ((xIntHi_rd & 0x7fffffff) >= 0x7ff00000) && ((xIntHi_rd & 0x80000000)==0))
    || (((xIntHi_rd & 0x7fffffff) >= OVRUDRFLWSMPLBOUND) && (x_rd > OVRFLWBOUND))
    || (((xIntHi_rd & 0x7fffffff) >= OVRUDRFLWSMPLBOUND) && (x_rd <= UNDERFLWBOUND))
    || (((xIntHi_rd & 0x7fffffff) >= OVRUDRFLWSMPLBOUND) && (x_rd <= DENORMBOUND))
     ,FALSE))
  {
    ASSIGN_LOW(res,exp_rd(LOW(x)));
    ASSIGN_UP(res,exp_ru(UP(x)));
    return res;
  }

  /* Test if argument is a denormal or zero */
  /* If we are here, we are sure to be neither +/- Inf nor NaN nor overflowed nor denormalized in the argument
     but we might be denormalized in the result 

     We continue the argument reduction for the quick phase and table reads for both phases
  */

  Mul12(&s1_ru,&s2_ru,msLog2Div2Lh,kd_ru);
  Mul12(&s1_rd,&s2_rd,msLog2Div2Lh,kd_rd);
  s3_ru = kd_ru * msLog2Div2Lm;
  s3_rd = kd_rd * msLog2Div2Lm;
  s4_ru = s2_ru + s3_ru; 
  s4_rd = s2_rd + s3_rd; 
  s5_ru = x_ru + s1_ru;
  s5_rd = x_rd + s1_rd;
  Add12Cond(rh_ru,rm_ru,s5_ru,s4_ru);
  Add12Cond(rh_rd,rm_rd,s5_rd,s4_rd);
  k_ru = shiftedXMultdb_ru.i[LO];
  k_rd = shiftedXMultdb_rd.i[LO];
  M_ru = k_ru >> L;
  M_rd = k_rd >> L;
  index1_ru = k_ru & INDEXMASK1;
  index1_rd = k_rd & INDEXMASK1;
  index2_ru = (k_ru & INDEXMASK2) >> LHALF;
  index2_rd = (k_rd & INDEXMASK2) >> LHALF;

  /* Table reads */
  tbl1h_ru = twoPowerIndex1[index1_ru].hi;
  tbl1h_rd = twoPowerIndex1[index1_rd].hi;
  tbl1m_ru = twoPowerIndex1[index1_ru].mi;
  tbl1m_rd = twoPowerIndex1[index1_rd].mi;
  tbl2h_ru = twoPowerIndex2[index2_ru].hi;
  tbl2h_rd = twoPowerIndex2[index2_rd].hi;
  tbl2m_ru = twoPowerIndex2[index2_ru].mi;
  tbl2m_rd = twoPowerIndex2[index2_rd].mi;




  /* No more underflow nor denormal is possible. There may be the case where
     M is 1024 and the value 2^M is to be multiplied may be less than 1
     So the final result will be normalized and representable by the multiplication must be 
     made in 2 steps
  */

  /* Quick phase starts here */

  rhSquare_ru = rh_ru * rh_ru;
  rhSquare_rd = rh_rd * rh_rd;
  rhC3_ru = c3 * rh_ru;
  rhC3_rd = c3 * rh_rd;
  rhSquareHalf_ru = 0.5 * rhSquare_ru;
  rhSquareHalf_rd = 0.5 * rhSquare_rd;
  monomialCube_ru = rhC3_ru * rhSquare_ru;
  monomialCube_rd = rhC3_rd * rhSquare_rd;
  rhFour_ru = rhSquare_ru * rhSquare_ru;
  rhFour_rd = rhSquare_rd * rhSquare_rd;
  monomialFour_ru = c4 * rhFour_ru;
  monomialFour_rd = c4 * rhFour_rd;
  highPoly_ru = monomialCube_ru + monomialFour_ru;
  highPoly_rd = monomialCube_rd + monomialFour_rd;
  highPolyWithSquare_ru = rhSquareHalf_ru + highPoly_ru;
  highPolyWithSquare_rd = rhSquareHalf_rd + highPoly_rd;
  Mul22(&tablesh_ru,&tablesl_ru,tbl1h_ru,tbl1m_ru,tbl2h_ru,tbl2m_ru);
  Mul22(&tablesh_rd,&tablesl_rd,tbl1h_rd,tbl1m_rd,tbl2h_rd,tbl2m_rd);
  t8_ru = rm_ru + highPolyWithSquare_ru;
  t8_rd = rm_rd + highPolyWithSquare_rd;
  t9_ru = rh_ru + t8_ru;
  t9_rd = rh_rd + t8_rd;
  t10_ru = tablesh_ru * t9_ru;
  t10_rd = tablesh_rd * t9_rd;
  Add12(t11_ru,t12_ru,tablesh_ru,t10_ru);
  Add12(t11_rd,t12_rd,tablesh_rd,t10_rd);
  t13_ru = t12_ru + tablesl_ru;
  t13_rd = t12_rd + tablesl_rd;
  Add12(polyTblh_ru,polyTblm_ru,t11_ru,t13_ru);
  Add12(polyTblh_rd,polyTblm_rd,t11_rd,t13_rd);
  
  /* Rounding test 
     Since we know that the result of the final multiplication with 2^M 
     will always be representable, we can do the rounding test on the 
     factors and multiply only the final result.
     We implement the multiplication in integer computations to overcome
     the problem of the non-representability of 2^1024 if M = 1024
  */

  if (infDone==1) res_rd=res_simple_rd;
  if (supDone==1) res_ru=res_simple_ru;

//  TEST_AND_COPY_RDRU_EXP(roundable,infDone,supDone,res_rd,polyTblh_rd,polyTblm_rd,res_ru,polyTblh_ru,polyTblm_ru,RDROUNDCST);
  db_number yh_rd, yl_rd, u53_rd, yh_ru, yl_ru, u53_ru;
  int yh_rd_neg, yl_rd_neg, yh_ru_neg, yl_ru_neg;
  int rd_ok, ru_ok;
  double save_res_rd=res_rd;
  double save_res_ru=res_ru;
  yh_rd.d = polyTblh_rd;    yl_rd.d = polyTblm_rd;
  yh_rd_neg = (yh_rd.i[HI] & 0x80000000);
  yl_rd_neg = (yl_rd.i[HI] & 0x80000000);
  yh_rd.l = yh_rd.l & 0x7fffffffffffffffLL;  /* compute the absolute value*/
  yl_rd.l = yl_rd.l & 0x7fffffffffffffffLL;  /* compute the absolute value*/
  u53_rd.l     = (yh_rd.l & ULL(7ff0000000000000)) +  ULL(0010000000000000);
  yh_ru.d = polyTblh_ru;    yl_ru.d = polyTblm_ru;
  yh_ru_neg = (yh_ru.i[HI] & 0x80000000);
  yl_ru_neg = (yl_ru.i[HI] & 0x80000000);
  yh_ru.l = yh_ru.l & 0x7fffffffffffffffLL;  /* compute the absolute value*/
  yl_ru.l = yl_ru.l & 0x7fffffffffffffffLL;  /* compute the absolute value*/
  u53_ru.l     = (yh_ru.l & ULL(7ff0000000000000)) +  ULL(0010000000000000);
  roundable = 0;
  rd_ok=(yl_rd.d > RDROUNDCST * u53_rd.d);
  ru_ok=(yl_ru.d > RDROUNDCST * u53_ru.d);
     if(yl_rd_neg) {  /* The case yl==0 is filtered by the above test*/
      /* return next down */
       yh_rd.d = polyTblh_rd;
      if(yh_rd_neg) yh_rd.l++;  else yh_rd.l--; /* Beware: fails for zero */
      res_rd = yh_rd.d;
    }
    else {
      res_rd = polyTblh_rd;
    }
    if(!yl_ru_neg) {  /* The case yl==0 is filtered by the above test*/
      /* return next up */
      yh_ru.d = polyTblh_ru;
      if(yh_ru_neg) yh_ru.l--;  else yh_ru.l++; /* Beware: fails for zero */
      res_ru = yh_ru.d;
    }
    else {
      res_ru = polyTblh_ru;
    }
  if(infDone) res_rd=save_res_rd;
  if(supDone) res_ru=save_res_ru;
  if(rd_ok && ru_ok){
    roundable=3;
  }
  else if (rd_ok){
    roundable=1;
  }
  else if (ru_ok){
     roundable=2;
  }
  resdb_rd.d = res_rd;
  resdb_ru.d = res_ru;

  if (roundable==3)
  {
    if (infDone==0){
      resdb_rd.i[HI] += M_rd << 20;
    }
    ASSIGN_LOW(res,resdb_rd.d);
    if (supDone==0){
      resdb_ru.i[HI] += M_ru << 20;
    }
    ASSIGN_UP(res,resdb_ru.d);
    return res;
  }
  if(roundable==1)
  {
    if(infDone==0){
      resdb_rd.i[HI] += M_rd << 20;
    }
    ASSIGN_LOW(res,resdb_rd.d);
    if(supDone==0){
    /* Rest of argument reduction for accurate phase */
    Mul133(&msLog2Div2LMultKh_ru,&msLog2Div2LMultKm_ru,&msLog2Div2LMultKl_ru,kd_ru,msLog2Div2Lh,msLog2Div2Lm,msLog2Div2Ll);
    t1_ru = x_ru + msLog2Div2LMultKh_ru;
    Add12Cond(rh_ru,t2_ru,t1_ru,msLog2Div2LMultKm_ru);
    Add12Cond(rm_ru,rl_ru,t2_ru,msLog2Div2LMultKl_ru);
    /* Table reads for accurate phase */
    tbl1l_ru = twoPowerIndex1[index1_ru].lo;
    tbl2l_ru = twoPowerIndex2[index2_ru].lo;
    /* Call accurate phase */
    exp_td_accurate(&polyTblh_ru, &polyTblm_ru, &polyTbll_ru, rh_ru, rm_ru, rl_ru, tbl1h_ru, tbl1m_ru, tbl1l_ru, tbl2h_ru, tbl2m_ru, tbl2l_ru); 
    /* Since the final multiplication is exact, we can do the final rounding before multiplying
       We overcome this way also the cases where the final result is not underflowed whereas the
       lower parts of the intermediate final result are.
    */
    RoundUpwards3(&res_ru,polyTblh_ru,polyTblm_ru,polyTbll_ru);
    /* Final multiplication with 2^M 
       We implement the multiplication in integer computations to overcome
       the problem of the non-representability of 2^1024 if M = 1024
    */
    resdb_ru.d = res_ru;
    resdb_ru.i[HI] += M_ru << 20;
    }
    ASSIGN_UP(res,resdb_ru.d);
    return res;
  } /* Accurate phase launched after rounding test*/
    
  if (roundable==2) {
    if (infDone==0){
    /* Rest of argument reduction for accurate phase */
    Mul133(&msLog2Div2LMultKh_rd,&msLog2Div2LMultKm_rd,&msLog2Div2LMultKl_rd,kd_rd,msLog2Div2Lh,msLog2Div2Lm,msLog2Div2Ll);
    t1_rd = x_rd + msLog2Div2LMultKh_rd;
    Add12Cond(rh_rd,t2_rd,t1_rd,msLog2Div2LMultKm_rd);
    Add12Cond(rm_rd,rl_rd,t2_rd,msLog2Div2LMultKl_rd);
    /* Table reads for accurate phase */
    tbl1l_rd = twoPowerIndex1[index1_rd].lo;
    tbl2l_rd = twoPowerIndex2[index2_rd].lo;
    /* Call accurate phase */
    exp_td_accurate(&polyTblh_rd, &polyTblm_rd, &polyTbll_rd, rh_rd, rm_rd, rl_rd, tbl1h_rd, tbl1m_rd, tbl1l_rd, tbl2h_rd, tbl2m_rd, tbl2l_rd); 
    /* Since the final multiplication is exact, we can do the final rounding before multiplying
       We overcome this way also the cases where the final result is not underflowed whereas the
       lower parts of the intermediate final result are.
    */

    RoundDownwards3(&res_rd,polyTblh_rd,polyTblm_rd,polyTbll_rd);
    /* Final multiplication with 2^M 
       We implement the multiplication in integer computations to overcome
       the problem of the non-representability of 2^1024 if M = 1024
    */

    resdb_rd.d = res_rd;
    resdb_rd.i[HI] += M_rd << 20;
    }
    ASSIGN_LOW(res,resdb_rd.d);
    if(supDone==0){
      resdb_ru.i[HI] += M_ru << 20;
    }
    ASSIGN_UP(res,resdb_ru.d);    
    return res;
  } /* Accurate phase launched after rounding test*/
  if(roundable==0)
  {
    if(supDone==0){
    /* Rest of argument reduction for accurate phase */
    Mul133(&msLog2Div2LMultKh_ru,&msLog2Div2LMultKm_ru,&msLog2Div2LMultKl_ru,kd_ru,msLog2Div2Lh,msLog2Div2Lm,msLog2Div2Ll);
    t1_ru = x_ru + msLog2Div2LMultKh_ru;
    Add12Cond(rh_ru,t2_ru,t1_ru,msLog2Div2LMultKm_ru);
    Add12Cond(rm_ru,rl_ru,t2_ru,msLog2Div2LMultKl_ru);
    /* Table reads for accurate phase */
    tbl1l_ru = twoPowerIndex1[index1_ru].lo;
    tbl2l_ru = twoPowerIndex2[index2_ru].lo;
    /* Call accurate phase */
    exp_td_accurate(&polyTblh_ru, &polyTblm_ru, &polyTbll_ru, rh_ru, rm_ru, rl_ru, tbl1h_ru, tbl1m_ru, tbl1l_ru, tbl2h_ru, tbl2m_ru, tbl2l_ru); 
    /* Since the final multiplication is exact, we can do the final rounding before multiplying
       We overcome this way also the cases where the final result is not underflowed whereas the
       lower parts of the intermediate final result are.
    */
    RoundUpwards3(&res_ru,polyTblh_ru,polyTblm_ru,polyTbll_ru);
    /* Final multiplication with 2^M 
       We implement the multiplication in integer computations to overcome
       the problem of the non-representability of 2^1024 if M = 1024
    */
    resdb_ru.d = res_ru;
    resdb_ru.i[HI] += M_ru << 20;
    }
    ASSIGN_UP(res,resdb_ru.d);
    if (infDone==0){
    /* Rest of argument reduction for accurate phase */
    Mul133(&msLog2Div2LMultKh_rd,&msLog2Div2LMultKm_rd,&msLog2Div2LMultKl_rd,kd_rd,msLog2Div2Lh,msLog2Div2Lm,msLog2Div2Ll);
    t1_rd = x_rd + msLog2Div2LMultKh_rd;
    Add12Cond(rh_rd,t2_rd,t1_rd,msLog2Div2LMultKm_rd);
    Add12Cond(rm_rd,rl_rd,t2_rd,msLog2Div2LMultKl_rd);
    /* Table reads for accurate phase */
    tbl1l_rd = twoPowerIndex1[index1_rd].lo;
    tbl2l_rd = twoPowerIndex2[index2_rd].lo;
    /* Call accurate phase */
    exp_td_accurate(&polyTblh_rd, &polyTblm_rd, &polyTbll_rd, rh_rd, rm_rd, rl_rd, tbl1h_rd, tbl1m_rd, tbl1l_rd, tbl2h_rd, tbl2m_rd, tbl2l_rd); 
    /* Since the final multiplication is exact, we can do the final rounding before multiplying
       We overcome this way also the cases where the final result is not underflowed whereas the
       lower parts of the intermediate final result are.
    */

    RoundDownwards3(&res_rd,polyTblh_rd,polyTblm_rd,polyTbll_rd);
    /* Final multiplication with 2^M 
       We implement the multiplication in integer computations to overcome
       the problem of the non-representability of 2^1024 if M = 1024
    */

    resdb_rd.d = res_rd;
    resdb_rd.i[HI] += M_rd << 20;
    }
    ASSIGN_LOW(res,resdb_rd.d);
    return res;
  } /* Accurate phase launched after rounding test*/

  return res;
}
#endif

