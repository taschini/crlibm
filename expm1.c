/*
 * Correctly rounded expm1 = e^x - 1
 *
 * Author : Christoph Lauter (ENS Lyon)
 *
 * This file is part of the crlibm library developed by the Arenaire
 * project at Ecole Normale Superieure de Lyon
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
*/

#include <stdio.h>
#include <stdlib.h>
#include "crlibm.h"
#include "crlibm_private.h"
#include "triple-double.h"
#include "expm1.h"


void expm1_direct_td(double *expm1h, double *expm1m, double *expm1l, 
		     double x, double xSqHalfh, double xSqHalfl, double xSqh, double xSql, int expoX) {
  double highPoly, tt1h, t1h, t1l, t2h, t2l, t3h, t3l, t4h, t4l, t5h, t5l, t6h, t6l;
  double tt6h, tt6m, tt6l, t7h, t7m, t7l, lowPolyh, lowPolym, lowPolyl;
  double fullHighPolyh, fullHighPolym, fullHighPolyl, polyh, polym, polyl;
  double xCubeh, xCubem, xCubel, tt7h, tt7m, tt7l, t8h, t8m, t8l;
  double expm1hover, expm1mover, expm1lover;
  double r1h, r1m, r1l, r2h, r2m, r2l, r3h, r3m, r3l;
  double rr1h, rr1m, rr1l, rr2h, rr2m, rr2l, rr3h, rr3m, rr3l;

  /* Double precision evaluation steps */
#if defined(PROCESSOR_HAS_FMA) && !defined(AVOID_FMA)
    highPoly = FMA(FMA(FMA(FMA(accuDirectpolyC15h ,x,accuDirectpolyC14h),x,accuDirectpolyC13h),x,
                               accuDirectpolyC12h),x,accuDirectpolyC11h);
#else
    highPoly = accuDirectpolyC11h + x * (accuDirectpolyC12h + x * (accuDirectpolyC13h + x * (
	       accuDirectpolyC14h + x *  accuDirectpolyC15h)));
#endif

    tt1h = x * highPoly;

    /* Triple-double steps for x + x^2/2 and x^3*/

    Add123(&lowPolyh,&lowPolym,&lowPolyl,x,xSqHalfh,xSqHalfl);
    Mul123(&xCubeh,&xCubem,&xCubel,x,xSqh,xSql);


    /* Double-double evaluation steps */

    Add12(t1h,t1l,accuDirectpolyC10h,tt1h);
    
    MulAdd212(&t2h,&t2l,accuDirectpolyC9h,accuDirectpolyC9m,x,t1h,t1l);
    MulAdd212(&t3h,&t3l,accuDirectpolyC8h,accuDirectpolyC8m,x,t2h,t2l);
    MulAdd212(&t4h,&t4l,accuDirectpolyC7h,accuDirectpolyC7m,x,t3h,t3l);
    MulAdd212(&t5h,&t5l,accuDirectpolyC6h,accuDirectpolyC6m,x,t4h,t4l);
    MulAdd212(&t6h,&t6l,accuDirectpolyC5h,accuDirectpolyC5m,x,t5h,t5l);

    /* Triple-double evaluation steps */

    Mul123(&tt6h,&tt6m,&tt6l,x,t6h,t6l);
    Add233(&t7h,&t7m,&t7l,accuDirectpolyC4h,accuDirectpolyC4m,tt6h,tt6m,tt6l);
   
    Mul133(&tt7h,&tt7m,&tt7l,x,t7h,t7m,t7l);
    Add33(&t8h,&t8m,&t8l,accuDirectpolyC3h,accuDirectpolyC3m,accuDirectpolyC3l,tt7h,tt7m,tt7l);

    Mul33(&fullHighPolyh,&fullHighPolym,&fullHighPolyl,xCubeh,xCubem,xCubel,t8h,t8m,t8l);

    Add33(&polyh,&polym,&polyl,lowPolyh,lowPolym,lowPolyl,fullHighPolyh,fullHighPolym,fullHighPolyl);

    /* Reconstruction steps */

    /* If we have not performed any range reduction, we have no reconstruction to do */
    if (expoX >= 0) {
      /* If we are here, we must perform reconstruction */

      /* First reconstruction step */
      Add133(&r1h,&r1m,&r1l,2,polyh,polym,polyl);
      Mul33(&rr1h,&rr1m,&rr1l,r1h,r1m,r1l,polyh,polym,polyl);
	
      if (expoX >= 1) {

	/* Second reconstruction step */
	Add133(&r2h,&r2m,&r2l,2,rr1h,rr1m,rr1l);
	Mul33(&rr2h,&rr2m,&rr2l,r2h,r2m,r2l,rr1h,rr1m,rr1l);

	if (expoX >= 2) {

	  /* Third reconstruction step */
	  Add133(&r3h,&r3m,&r3l,2,rr2h,rr2m,rr2l);
	  Mul33(&rr3h,&rr3m,&rr3l,r3h,r3m,r3l,rr2h,rr2m,rr2l);

	  /* expoX may be maximally 2 */

	  expm1hover = rr3h;
	  expm1mover = rr3m;
	  expm1lover = rr3l;

	} else {
	  expm1hover = rr2h;
	  expm1mover = rr2m;
	  expm1lover = rr1l;
	}

      } else {
	expm1hover = rr1h;
	expm1mover = rr1m;
	expm1lover = rr1l;
      }

    } else {
      expm1hover = polyh;
      expm1mover = polym;
      expm1lover = polyl;
    }

    /* Renormalize before returning */

    Renormalize3(expm1h,expm1m,expm1l,expm1hover,expm1mover,expm1lover);
}

void expm1_common_td(double *expm1h, double *expm1m, double *expm1l, 
		     double rh, double rm, double rl, 
		     double tbl1h, double tbl1m, double tbl1l, 
		     double tbl2h, double tbl2m, double tbl2l, 
		     int M) {
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
  double exph, expm, expl, expm1hover, expm1mover, expm1lover;
  db_number polyWithTableshdb, polyWithTablesmdb, polyWithTablesldb;

  /* Polynomial approximation - double precision steps */

#if defined(PROCESSOR_HAS_FMA) && !defined(AVOID_FMA)
  highPoly = FMA(FMA(accuCommonpolyC7h,rh,accuCommonpolyC6h),rh,accuCommonpolyC5h);
#else
  highPoly = accuCommonpolyC5h + rh * (accuCommonpolyC6h + rh * accuCommonpolyC7h);
#endif

  /* Polynomial approximation - double-double precision steps */

  Mul12(&t1h,&t1l,rh,highPoly);
  Add22(&t2h,&t2l,accuCommonpolyC4h,accuCommonpolyC4m,t1h,t1l);
  Mul122(&t3h,&t3l,rh,t2h,t2l);
  Add22(&t4h,&t4l,accuCommonpolyC3h,accuCommonpolyC3m,t3h,t3l);

  Mul12(&rhSquareh,&rhSquarel,rh,rh);
  Mul123(&rhCubeh,&rhCubem,&rhCubel,rh,rhSquareh,rhSquarel);

  rhSquareHalfh = 0.5 * rhSquareh;
  rhSquareHalfl = 0.5 * rhSquarel;  

  /* Polynomial approximation - triple-double precision steps */

  Renormalize3(&lowPolyh,&lowPolym,&lowPolyl,rh,rhSquareHalfh,rhSquareHalfl);

  Mul233(&highPolyMulth,&highPolyMultm,&highPolyMultl,t4h,t4l,rhCubeh,rhCubem,rhCubel);

  Add33(&ph,&pm,&pl,lowPolyh,lowPolym,lowPolyl,highPolyMulth,highPolyMultm,highPolyMultl);

  /* Reconstruction */

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

  /* Multiplication by 2^(M) 

     We perform it in integer to overcome the non-representability of 2^(1024) 
     This case is possible for M = 1024 and polyWithTablesh < 1

     The overlap in the triple-double polyWithTables[hml] stays unchanged.

  */

  polyWithTableshdb.d = polyWithTablesh;
  polyWithTablesmdb.d = polyWithTablesm;
  polyWithTablesldb.d = polyWithTablesl;

  polyWithTableshdb.i[HI] += M << 20;
  polyWithTablesmdb.i[HI] += M << 20;
  polyWithTablesldb.i[HI] += M << 20;

  exph = polyWithTableshdb.d;
  expm = polyWithTablesmdb.d;
  expl = polyWithTablesldb.d;

  /* Substraction of -1 
     
     We use a conditional Add133 
  */

  Add133(&expm1hover,&expm1mover,&expm1lover,-1,exph,expm,expl);

  /* Renormalization */

  Renormalize3(expm1h,expm1m,expm1l,expm1hover,expm1mover,expm1lover);
}


double expm1_rn(double x) {
  db_number xdb, scaledb, shiftedXMultdb, polyTblhdb, polyTblmdb;
  int xIntHi, expoX, k, M, index1, index2;
  double highPoly, tt1h, t1h, t1l, xSqh, xSql, xSqHalfh, xSqHalfl, xCubeh, xCubel, t2h, t2l, templ, tt3h, tt3l;
  double polyh, polyl, expm1h, expm1m, expm1l;
  double r1h, r1l, r1t, rr1h, rr1l;
  double r2h, r2l, r2t, rr2h, rr2l;
  double r3h, r3l, r3t, rr3h, rr3l;
  double xMultLog2InvMult2L, shiftedXMult, kd, s1, s2, s3, s4, s5, rh, rm, rl;
  double rhSquare, rhC3, rhSquareHalf, monomialCube, rhFour, monomialFour;
  double tbl1h, tbl1m, tbl1l, tbl2h, tbl2m, tbl2l;
  double highPolyWithSquare, tablesh, tablesl, t8, t9, t10, t11, t12, t13;
  double exph, expm, t1, t2, t3;
  double msLog2Div2LMultKh, msLog2Div2LMultKm, msLog2Div2LMultKl;


  xdb.d = x; 

  /* Strip off the sign of x for the following tests */

  xIntHi = xdb.i[HI] & 0x7fffffff;

  /* Test if we are so small that we can return (a corrected) x as correct rounding */
  if (xIntHi < RETURNXBOUND) {
    return x;
  }


  /* Filter out special cases like overflow, -1 in result, infinities and NaNs 
     The filters are not sharp, we have positive arguments that flow through
  */
  if (xIntHi >= SIMPLEOVERFLOWBOUND) {
    /* Test if we are +/-inf or NaN */
    if (xIntHi >= 0x7ff00000) {
      /* Test if NaN */
      if (((xIntHi & 0x000fffff) | xdb.i[LO]) != 0) {
	/* NaN */
	return x+x;  /* return NaN */
      }
      /* Test if +inf or -inf */
      if (xdb.i[HI] > 0) {
	/* +inf */
	return x+x;  /* return +inf */
      }
      
      /* If we are here, we are -inf */
      return -1.0;
    }

    /* If we are here, we are overflowed or a common case that flows through */

    /* Test if we are actually overflowed */
    if (x > OVERFLOWBOUND) {
      return LARGEST * LARGEST;  /* return +inf and set flag */
    }
  }
  
  /* Test if we know already that we are -1.0 (+ correction depending on rounding mode) in result */
  if (x < MINUSONEBOUND) {
    return -1.0;
  }

  /* Test if we have |x| <= 1/4-1/2ulp(1/4) for knowing if we use exp(x) or approximate directly */

  if (xIntHi < DIRECTINTERVALBOUND) {
    /* We approximate expm1 directly after a range reduction as follows

       expm1(x) = (expm1(x/2) + 2) * expm1(x/2)

       We perform the range reduction in such a way that finally |x| < 1/32 
    */

    /* Extract the exponent of |x| and add 5 (2^5 = 32) */
    expoX = ((xIntHi & 0x7ff00000) >> 20) - (1023 - 5);
    
    /* If this particularily biased exponent expoX is negative, we are already less than 1/32 */
    if (expoX >= 0) {
      /* If we are here, we must perform range reduction */


      /* We start by producing 2^(-expoX-1) as a factor to x */
      scaledb.i[HI] = ((-expoX-1)+1023) << 20;
      scaledb.i[LO] = 0;
      
      /* We multiply x with the scale */
      x = scaledb.d * x;
    }
    
    /* Here, we have always |x| < 1/32 */


    /* Double precision evaluation steps */
#if defined(PROCESSOR_HAS_FMA) && !defined(AVOID_FMA)
    highPoly = FMA(FMA(FMA(FMA(FMA(quickDirectpolyC9h ,x,quickDirectpolyC8h),x,quickDirectpolyC7h),x,
                                   quickDirectpolyC6h),x,quickDirectpolyC5h),x,quickDirectpolyC4h);
#else
    highPoly = quickDirectpolyC4h + x * (quickDirectpolyC5h + x * (quickDirectpolyC6h + x * (
	       quickDirectpolyC7h + x * (quickDirectpolyC8h + x *  quickDirectpolyC9h))));
#endif

    /* Double-double evaluation steps */
    tt1h = x * highPoly;

    Mul12(&xSqh,&xSql,x,x);
    xSqHalfh = 0.5 * xSqh;
    xSqHalfl = 0.5 * xSql;
    Add12(t2h,templ,x,xSqHalfh);
    t2l = templ + xSqHalfl;
    
    Add12(t1h,t1l,quickDirectpolyC3h,tt1h);
    Mul122(&xCubeh,&xCubel,x,xSqh,xSql);
    Mul22(&tt3h,&tt3l,xCubeh,xCubel,t1h,t1l);

    Add22(&polyh,&polyl,t2h,t2l,tt3h,tt3l);

    /* Reconstruction */

    /* If we have not performed any range reduction, we have no reconstruction to do */
    if (expoX >= 0) {
      /* If we are here, we must perform reconstruction */

      /* First reconstruction step */
      Add12(r1h,r1t,2,polyh);
      r1l = r1t + polyl;
      Mul22(&rr1h,&rr1l,r1h,r1l,polyh,polyl);

      if (expoX >= 1) {

	/* Second reconstruction step */
	Add12(r2h,r2t,2,rr1h);
	r2l = r2t + rr1l;
	Mul22(&rr2h,&rr2l,r2h,r2l,rr1h,rr1l);

	if (expoX >= 2) {

	  /* Third reconstruction step */
	  Add12(r3h,r3t,2,rr2h);
	  r3l = r3t + rr2l;
	  Mul22(&rr3h,&rr3l,r3h,r3l,rr2h,rr2l);
	  
	  /* expoX may be maximally 2 */

	  expm1h = rr3h;
	  expm1m = rr3l;

	} else {
	  expm1h = rr2h;
	  expm1m = rr2l;
	}

      } else {
	expm1h = rr1h;
	expm1m = rr1l;
      }

    } else {
      expm1h = polyh;
      expm1m = polyl;
    }

    /* Rounding test */
    if(expm1h == (expm1h + (expm1m * ROUNDCSTDIRECTRN)))
     return expm1h;
   else 
     {
       expm1_direct_td(&expm1h, &expm1m, &expm1l, x, xSqHalfh, xSqHalfl, xSqh, xSql, expoX);
      
       ReturnRoundToNearest3(expm1h, expm1m, expm1l);

     } /* Accurate phase launched */

    /* We cannot be here, since we return in all cases before */
  }

  /* If we are here, we can use expm1(x) = exp(x) - 1 */

  /* Range reduction - exact part: compute k as double and as int */

  xMultLog2InvMult2L = x * log2InvMult2L;
  shiftedXMult = xMultLog2InvMult2L + shiftConst;
  kd = shiftedXMult - shiftConst;
  shiftedXMultdb.d = shiftedXMult;
  k = shiftedXMultdb.i[LO];
  M = k >> 12;
  index1 = k & INDEXMASK1;
  index2 = (k & INDEXMASK2) >> 6;

  /* Range reduction - part affected by error - must be redone in accurate phase */
  Mul12(&s1,&s2,msLog2Div2Lh,kd);
  s3 = kd * msLog2Div2Lm;
  s4 = s2 + s3; 
  s5 = x + s1;
  Add12Cond(rh,rm,s5,s4);

  /* Table reads - read only two double-doubles by now */
  tbl1h = twoPowerIndex1[index1].hi;
  tbl1m = twoPowerIndex1[index1].mi;
  tbl2h = twoPowerIndex2[index2].hi;
  tbl2m = twoPowerIndex2[index2].mi;

  /* Quick phase starts here */

  rhSquare = rh * rh;
  rhC3 = quickCommonpolyC3h * rh;

  rhSquareHalf = 0.5 * rhSquare;
  monomialCube = rhC3 * rhSquare;
  rhFour = rhSquare * rhSquare;

  monomialFour = quickCommonpolyC4h * rhFour;
  
  highPoly = monomialCube + monomialFour;

  highPolyWithSquare = rhSquareHalf + highPoly;
  
  /* Reconstruction: integration of table values */
  
  Mul22(&tablesh,&tablesl,tbl1h,tbl1m,tbl2h,tbl2m);

  t8 = rm + highPolyWithSquare;
  t9 = rh + t8;

  t10 = tablesh * t9;
  
  Add12(t11,t12,tablesh,t10);
  t13 = t12 + tablesl;
  Add12(polyTblhdb.d,polyTblmdb.d,t11,t13);
  
  /* Reconstruction: multiplication by 2^M */

  /* Implement the multiplication by multiplication to overcome the
     problem of the non-representability of 2^1024 (M = 1024)
     This case is possible if polyTblhdb.d < 1
  */
  
  polyTblhdb.i[HI] += M << 20;
  polyTblmdb.i[HI] += M << 20;

  exph = polyTblhdb.d;
  expm = polyTblmdb.d;

  /* Substraction of 1 

     Testing if the operation is necessary is more expensive than 
     performing it in any case.

     We may cancellate at most 2 bits in the subtraction for 
     arguments 1/4 <= x <= ln(2) (0.25 <= x <= 0.69) 
     We must therefore use conditional Add12s

     Since we perform a substraction, we may not have addition overflow towards +inf

  */

  Add12Cond(t1,t2,-1,exph);
  t3 = t2 + expm;
  Add12Cond(expm1h,expm1m,t1,t3);


  /* Rounding test */
  if(expm1h == (expm1h + (expm1m * ROUNDCSTCOMMONRN))) {
    return expm1h;
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
    expm1_common_td(&expm1h, &expm1m, &expm1l, rh, rm, rl, tbl1h, tbl1m, tbl1l, tbl2h, tbl2m, tbl2l, M); 
    
    /* Final rounding */

    ReturnRoundToNearest3(expm1h, expm1m, expm1l);
  } /* Accurate phase launched */
  
  /* We cannot be here since we return before in any case */
}

double expm1_rd(double x) {
  db_number xdb, scaledb, shiftedXMultdb, polyTblhdb, polyTblmdb;
  int xIntHi, expoX, k, M, index1, index2;
  double highPoly, tt1h, t1h, t1l, xSqh, xSql, xSqHalfh, xSqHalfl, xCubeh, xCubel, t2h, t2l, templ, tt3h, tt3l;
  double polyh, polyl, expm1h, expm1m, expm1l;
  double r1h, r1l, r1t, rr1h, rr1l;
  double r2h, r2l, r2t, rr2h, rr2l;
  double r3h, r3l, r3t, rr3h, rr3l;
  double xMultLog2InvMult2L, shiftedXMult, kd, s1, s2, s3, s4, s5, rh, rm, rl;
  double rhSquare, rhC3, rhSquareHalf, monomialCube, rhFour, monomialFour;
  double tbl1h, tbl1m, tbl1l, tbl2h, tbl2m, tbl2l;
  double highPolyWithSquare, tablesh, tablesl, t8, t9, t10, t11, t12, t13;
  double exph, expm, t1, t2, t3;
  double msLog2Div2LMultKh, msLog2Div2LMultKm, msLog2Div2LMultKl;


  xdb.d = x; 

  /* Strip off the sign of x for the following tests */

  xIntHi = xdb.i[HI] & 0x7fffffff;

  /* Test if we are so small that we can return (a corrected) x as correct rounding */
  if (xIntHi < RETURNXBOUND) {
    /* The only algebraic result is 0 for x = 0; in this case, we can return x
       The truncation rest x^2/2 + x^3/6 + ... is always positive 
       but less than 1 ulp in this case, so we round down by returning x
    */
    return x;
  }


  /* Filter out special cases like overflow, -1 in result, infinities and NaNs 
     The filters are not sharp, we have positive arguments that flow through
  */
  if (xIntHi >= SIMPLEOVERFLOWBOUND) {
    /* Test if we are +/-inf or NaN */
    if (xIntHi >= 0x7ff00000) {
      /* Test if NaN */
      if (((xIntHi & 0x000fffff) | xdb.i[LO]) != 0) {
	/* NaN */
	return x+x;  /* return NaN */
      }
      /* Test if +inf or -inf */
      if (xdb.i[HI] > 0) {
	/* +inf */
	return x+x;  /* return +inf */
      }
      
      /* If we are here, we are -inf */
      return -1.0;
    }

    /* If we are here, we are overflowed or a common case that flows through */

    /* Test if we are actually overflowed */
    if (x > OVERFLOWBOUND) {
      /* We would be overflowed but as we are rounding downwards
	 the nearest number lesser than the exact result is the greatest 
	 normal. In any case, we must raise the inexact flag.
      */
      return LARGEST * (1.0 + SMALLEST);
    }
  }
  
  /* Test if we know already that we are -1.0 (+ correction depending on rounding mode) in result */
  if (x < MINUSONEBOUND) {
    /* We round down, so we are -1.0 */
    return -1.0;
  }

  /* Test if we have |x| <= 1/4-1/2ulp(1/4) for knowing if we use exp(x) or approximate directly */

  if (xIntHi < DIRECTINTERVALBOUND) {
    /* We approximate expm1 directly after a range reduction as follows

       expm1(x) = (expm1(x/2) + 2) * expm1(x/2)

       We perform the range reduction in such a way that finally |x| < 1/32 
    */

    /* Extract the exponent of |x| and add 5 (2^5 = 32) */
    expoX = ((xIntHi & 0x7ff00000) >> 20) - (1023 - 5);
    
    /* If this particularily biased exponent expoX is negative, we are already less than 1/32 */
    if (expoX >= 0) {
      /* If we are here, we must perform range reduction */


      /* We start by producing 2^(-expoX-1) as a factor to x */
      scaledb.i[HI] = ((-expoX-1)+1023) << 20;
      scaledb.i[LO] = 0;
      
      /* We multiply x with the scale */
      x = scaledb.d * x;
    }
    
    /* Here, we have always |x| < 1/32 */


    /* Double precision evaluation steps */
#if defined(PROCESSOR_HAS_FMA) && !defined(AVOID_FMA)
    highPoly = FMA(FMA(FMA(FMA(FMA(quickDirectpolyC9h ,x,quickDirectpolyC8h),x,quickDirectpolyC7h),x,
                                   quickDirectpolyC6h),x,quickDirectpolyC5h),x,quickDirectpolyC4h);
#else
    highPoly = quickDirectpolyC4h + x * (quickDirectpolyC5h + x * (quickDirectpolyC6h + x * (
	       quickDirectpolyC7h + x * (quickDirectpolyC8h + x *  quickDirectpolyC9h))));
#endif

    /* Double-double evaluation steps */
    tt1h = x * highPoly;

    Mul12(&xSqh,&xSql,x,x);
    xSqHalfh = 0.5 * xSqh;
    xSqHalfl = 0.5 * xSql;
    Add12(t2h,templ,x,xSqHalfh);
    t2l = templ + xSqHalfl;
    
    Add12(t1h,t1l,quickDirectpolyC3h,tt1h);
    Mul122(&xCubeh,&xCubel,x,xSqh,xSql);
    Mul22(&tt3h,&tt3l,xCubeh,xCubel,t1h,t1l);

    Add22(&polyh,&polyl,t2h,t2l,tt3h,tt3l);

    /* Reconstruction */

    /* If we have not performed any range reduction, we have no reconstruction to do */
    if (expoX >= 0) {
      /* If we are here, we must perform reconstruction */

      /* First reconstruction step */
      Add12(r1h,r1t,2,polyh);
      r1l = r1t + polyl;
      Mul22(&rr1h,&rr1l,r1h,r1l,polyh,polyl);

      if (expoX >= 1) {

	/* Second reconstruction step */
	Add12(r2h,r2t,2,rr1h);
	r2l = r2t + rr1l;
	Mul22(&rr2h,&rr2l,r2h,r2l,rr1h,rr1l);

	if (expoX >= 2) {

	  /* Third reconstruction step */
	  Add12(r3h,r3t,2,rr2h);
	  r3l = r3t + rr2l;
	  Mul22(&rr3h,&rr3l,r3h,r3l,rr2h,rr2l);
	  
	  /* expoX may be maximally 2 */

	  expm1h = rr3h;
	  expm1m = rr3l;

	} else {
	  expm1h = rr2h;
	  expm1m = rr2l;
	}

      } else {
	expm1h = rr1h;
	expm1m = rr1l;
      }

    } else {
      expm1h = polyh;
      expm1m = polyl;
    }

    /* Rounding test */
    TEST_AND_RETURN_RD(expm1h, expm1m, ROUNDCSTDIRECTRD);
    {
      expm1_direct_td(&expm1h, &expm1m, &expm1l, x, xSqHalfh, xSqHalfl, xSqh, xSql, expoX);
      
      ReturnRoundDownwards3(expm1h, expm1m, expm1l);
      
    } /* Accurate phase launched */

    /* We cannot be here, since we return in all cases before */
  }

  /* If we are here, we can use expm1(x) = exp(x) - 1 */

  /* Range reduction - exact part: compute k as double and as int */

  xMultLog2InvMult2L = x * log2InvMult2L;
  shiftedXMult = xMultLog2InvMult2L + shiftConst;
  kd = shiftedXMult - shiftConst;
  shiftedXMultdb.d = shiftedXMult;
  k = shiftedXMultdb.i[LO];
  M = k >> 12;
  index1 = k & INDEXMASK1;
  index2 = (k & INDEXMASK2) >> 6;

  /* Range reduction - part affected by error - must be redone in accurate phase */
  Mul12(&s1,&s2,msLog2Div2Lh,kd);
  s3 = kd * msLog2Div2Lm;
  s4 = s2 + s3; 
  s5 = x + s1;
  Add12Cond(rh,rm,s5,s4);

  /* Table reads - read only two double-doubles by now */
  tbl1h = twoPowerIndex1[index1].hi;
  tbl1m = twoPowerIndex1[index1].mi;
  tbl2h = twoPowerIndex2[index2].hi;
  tbl2m = twoPowerIndex2[index2].mi;

  /* Quick phase starts here */

  rhSquare = rh * rh;
  rhC3 = quickCommonpolyC3h * rh;

  rhSquareHalf = 0.5 * rhSquare;
  monomialCube = rhC3 * rhSquare;
  rhFour = rhSquare * rhSquare;

  monomialFour = quickCommonpolyC4h * rhFour;
  
  highPoly = monomialCube + monomialFour;

  highPolyWithSquare = rhSquareHalf + highPoly;
  
  /* Reconstruction: integration of table values */
  
  Mul22(&tablesh,&tablesl,tbl1h,tbl1m,tbl2h,tbl2m);

  t8 = rm + highPolyWithSquare;
  t9 = rh + t8;

  t10 = tablesh * t9;
  
  Add12(t11,t12,tablesh,t10);
  t13 = t12 + tablesl;
  Add12(polyTblhdb.d,polyTblmdb.d,t11,t13);
  
  /* Reconstruction: multiplication by 2^M */

  /* Implement the multiplication by multiplication to overcome the
     problem of the non-representability of 2^1024 (M = 1024)
     This case is possible if polyTblhdb.d < 1
  */
  
  polyTblhdb.i[HI] += M << 20;
  polyTblmdb.i[HI] += M << 20;

  exph = polyTblhdb.d;
  expm = polyTblmdb.d;

  /* Substraction of 1 

     Testing if the operation is necessary is more expensive than 
     performing it in any case.

     We may cancellate at most 2 bits in the subtraction for 
     arguments 1/4 <= x <= ln(2) (0.25 <= x <= 0.69) 
     We must therefore use conditional Add12s

     Since we perform a substraction, we may not have addition overflow towards +inf

  */

  Add12Cond(t1,t2,-1,exph);
  t3 = t2 + expm;
  Add12Cond(expm1h,expm1m,t1,t3);


  /* Rounding test */
  TEST_AND_RETURN_RD(expm1h, expm1m, ROUNDCSTCOMMONRD);
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
    expm1_common_td(&expm1h, &expm1m, &expm1l, rh, rm, rl, tbl1h, tbl1m, tbl1l, tbl2h, tbl2m, tbl2l, M); 
    
    /* Final rounding */

    ReturnRoundDownwards3(expm1h, expm1m, expm1l);
  } /* Accurate phase launched */
  
  /* We cannot be here since we return before in any case */
}

double expm1_ru(double x) {
  db_number xdb, scaledb, shiftedXMultdb, polyTblhdb, polyTblmdb;
  int xIntHi, expoX, k, M, index1, index2;
  double highPoly, tt1h, t1h, t1l, xSqh, xSql, xSqHalfh, xSqHalfl, xCubeh, xCubel, t2h, t2l, templ, tt3h, tt3l;
  double polyh, polyl, expm1h, expm1m, expm1l;
  double r1h, r1l, r1t, rr1h, rr1l;
  double r2h, r2l, r2t, rr2h, rr2l;
  double r3h, r3l, r3t, rr3h, rr3l;
  double xMultLog2InvMult2L, shiftedXMult, kd, s1, s2, s3, s4, s5, rh, rm, rl;
  double rhSquare, rhC3, rhSquareHalf, monomialCube, rhFour, monomialFour;
  double tbl1h, tbl1m, tbl1l, tbl2h, tbl2m, tbl2l;
  double highPolyWithSquare, tablesh, tablesl, t8, t9, t10, t11, t12, t13;
  double exph, expm, t1, t2, t3;
  double msLog2Div2LMultKh, msLog2Div2LMultKm, msLog2Div2LMultKl;


  xdb.d = x; 

  /* Strip off the sign of x for the following tests */

  xIntHi = xdb.i[HI] & 0x7fffffff;

  /* Test if we are so small that we can return (a corrected) x as correct rounding */
  if (xIntHi < RETURNXBOUND) {
    /* The only algebraic result is 0 for x = 0; in this case, we can return x
       The truncation rest x^2/2 + x^3/6 + ... is always positive 
       but less than 1 ulp in this case, so we round by adding 1 ulp 
    */
    if (xdb.i[HI] & 0x80000000) {
      /* x is negative 
	 We add 1 ulp by subtracting 1 in long
      */
      xdb.l--;
    } else {
      /* x is positive 
	 We add 1 ulp by adding 1 in long
      */
      xdb.l++;
    }
    return xdb.d;
  }


  /* Filter out special cases like overflow, -1 in result, infinities and NaNs 
     The filters are not sharp, we have positive arguments that flow through
  */
  if (xIntHi >= SIMPLEOVERFLOWBOUND) {
    /* Test if we are +/-inf or NaN */
    if (xIntHi >= 0x7ff00000) {
      /* Test if NaN */
      if (((xIntHi & 0x000fffff) | xdb.i[LO]) != 0) {
	/* NaN */
	return x+x;  /* return NaN */
      }
      /* Test if +inf or -inf */
      if (xdb.i[HI] > 0) {
	/* +inf */
	return x+x;  /* return +inf */
      }
      
      /* If we are here, we are -inf */
      return -1.0;
    }

    /* If we are here, we are overflowed or a common case that flows through */

    /* Test if we are actually overflowed */
    if (x > OVERFLOWBOUND) {
      return LARGEST * LARGEST;  /* return +inf and set flag */
    }
  }
  
  /* Test if we know already that we are -1.0 (+ correction depending on rounding mode) in result */
  if (x < MINUSONEBOUND) {
    /* Round up so we are -1.0 + 1ulp */
    return MINUSONEPLUSONEULP;
  }

  /* Test if we have |x| <= 1/4-1/2ulp(1/4) for knowing if we use exp(x) or approximate directly */

  if (xIntHi < DIRECTINTERVALBOUND) {
    /* We approximate expm1 directly after a range reduction as follows

       expm1(x) = (expm1(x/2) + 2) * expm1(x/2)

       We perform the range reduction in such a way that finally |x| < 1/32 
    */

    /* Extract the exponent of |x| and add 5 (2^5 = 32) */
    expoX = ((xIntHi & 0x7ff00000) >> 20) - (1023 - 5);
    
    /* If this particularily biased exponent expoX is negative, we are already less than 1/32 */
    if (expoX >= 0) {
      /* If we are here, we must perform range reduction */


      /* We start by producing 2^(-expoX-1) as a factor to x */
      scaledb.i[HI] = ((-expoX-1)+1023) << 20;
      scaledb.i[LO] = 0;
      
      /* We multiply x with the scale */
      x = scaledb.d * x;
    }
    
    /* Here, we have always |x| < 1/32 */


    /* Double precision evaluation steps */
#if defined(PROCESSOR_HAS_FMA) && !defined(AVOID_FMA)
    highPoly = FMA(FMA(FMA(FMA(FMA(quickDirectpolyC9h ,x,quickDirectpolyC8h),x,quickDirectpolyC7h),x,
                                   quickDirectpolyC6h),x,quickDirectpolyC5h),x,quickDirectpolyC4h);
#else
    highPoly = quickDirectpolyC4h + x * (quickDirectpolyC5h + x * (quickDirectpolyC6h + x * (
	       quickDirectpolyC7h + x * (quickDirectpolyC8h + x *  quickDirectpolyC9h))));
#endif

    /* Double-double evaluation steps */
    tt1h = x * highPoly;

    Mul12(&xSqh,&xSql,x,x);
    xSqHalfh = 0.5 * xSqh;
    xSqHalfl = 0.5 * xSql;
    Add12(t2h,templ,x,xSqHalfh);
    t2l = templ + xSqHalfl;
    
    Add12(t1h,t1l,quickDirectpolyC3h,tt1h);
    Mul122(&xCubeh,&xCubel,x,xSqh,xSql);
    Mul22(&tt3h,&tt3l,xCubeh,xCubel,t1h,t1l);

    Add22(&polyh,&polyl,t2h,t2l,tt3h,tt3l);

    /* Reconstruction */

    /* If we have not performed any range reduction, we have no reconstruction to do */
    if (expoX >= 0) {
      /* If we are here, we must perform reconstruction */

      /* First reconstruction step */
      Add12(r1h,r1t,2,polyh);
      r1l = r1t + polyl;
      Mul22(&rr1h,&rr1l,r1h,r1l,polyh,polyl);

      if (expoX >= 1) {

	/* Second reconstruction step */
	Add12(r2h,r2t,2,rr1h);
	r2l = r2t + rr1l;
	Mul22(&rr2h,&rr2l,r2h,r2l,rr1h,rr1l);

	if (expoX >= 2) {

	  /* Third reconstruction step */
	  Add12(r3h,r3t,2,rr2h);
	  r3l = r3t + rr2l;
	  Mul22(&rr3h,&rr3l,r3h,r3l,rr2h,rr2l);
	  
	  /* expoX may be maximally 2 */

	  expm1h = rr3h;
	  expm1m = rr3l;

	} else {
	  expm1h = rr2h;
	  expm1m = rr2l;
	}

      } else {
	expm1h = rr1h;
	expm1m = rr1l;
      }

    } else {
      expm1h = polyh;
      expm1m = polyl;
    }

    /* Rounding test */
    TEST_AND_RETURN_RU(expm1h, expm1m, ROUNDCSTDIRECTRD);
    {
      expm1_direct_td(&expm1h, &expm1m, &expm1l, x, xSqHalfh, xSqHalfl, xSqh, xSql, expoX);
      
      ReturnRoundUpwards3(expm1h, expm1m, expm1l);

    } /* Accurate phase launched */

    /* We cannot be here, since we return in all cases before */
  }

  /* If we are here, we can use expm1(x) = exp(x) - 1 */

  /* Range reduction - exact part: compute k as double and as int */

  xMultLog2InvMult2L = x * log2InvMult2L;
  shiftedXMult = xMultLog2InvMult2L + shiftConst;
  kd = shiftedXMult - shiftConst;
  shiftedXMultdb.d = shiftedXMult;
  k = shiftedXMultdb.i[LO];
  M = k >> 12;
  index1 = k & INDEXMASK1;
  index2 = (k & INDEXMASK2) >> 6;

  /* Range reduction - part affected by error - must be redone in accurate phase */
  Mul12(&s1,&s2,msLog2Div2Lh,kd);
  s3 = kd * msLog2Div2Lm;
  s4 = s2 + s3; 
  s5 = x + s1;
  Add12Cond(rh,rm,s5,s4);

  /* Table reads - read only two double-doubles by now */
  tbl1h = twoPowerIndex1[index1].hi;
  tbl1m = twoPowerIndex1[index1].mi;
  tbl2h = twoPowerIndex2[index2].hi;
  tbl2m = twoPowerIndex2[index2].mi;

  /* Quick phase starts here */

  rhSquare = rh * rh;
  rhC3 = quickCommonpolyC3h * rh;

  rhSquareHalf = 0.5 * rhSquare;
  monomialCube = rhC3 * rhSquare;
  rhFour = rhSquare * rhSquare;

  monomialFour = quickCommonpolyC4h * rhFour;
  
  highPoly = monomialCube + monomialFour;

  highPolyWithSquare = rhSquareHalf + highPoly;
  
  /* Reconstruction: integration of table values */
  
  Mul22(&tablesh,&tablesl,tbl1h,tbl1m,tbl2h,tbl2m);

  t8 = rm + highPolyWithSquare;
  t9 = rh + t8;

  t10 = tablesh * t9;
  
  Add12(t11,t12,tablesh,t10);
  t13 = t12 + tablesl;
  Add12(polyTblhdb.d,polyTblmdb.d,t11,t13);
  
  /* Reconstruction: multiplication by 2^M */

  /* Implement the multiplication by multiplication to overcome the
     problem of the non-representability of 2^1024 (M = 1024)
     This case is possible if polyTblhdb.d < 1
  */
  
  polyTblhdb.i[HI] += M << 20;
  polyTblmdb.i[HI] += M << 20;

  exph = polyTblhdb.d;
  expm = polyTblmdb.d;

  /* Substraction of 1 

     Testing if the operation is necessary is more expensive than 
     performing it in any case.

     We may cancellate at most 2 bits in the subtraction for 
     arguments 1/4 <= x <= ln(2) (0.25 <= x <= 0.69) 
     We must therefore use conditional Add12s

     Since we perform a substraction, we may not have addition overflow towards +inf

  */

  Add12Cond(t1,t2,-1,exph);
  t3 = t2 + expm;
  Add12Cond(expm1h,expm1m,t1,t3);


  /* Rounding test */
  TEST_AND_RETURN_RU(expm1h, expm1m, ROUNDCSTCOMMONRD);
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
    expm1_common_td(&expm1h, &expm1m, &expm1l, rh, rm, rl, tbl1h, tbl1m, tbl1l, tbl2h, tbl2m, tbl2l, M); 
    
    /* Final rounding */

    ReturnRoundUpwards3(expm1h, expm1m, expm1l);
  } /* Accurate phase launched */
  
  /* We cannot be here since we return before in any case */
}

double expm1_rz(double x) {
  db_number xdb, scaledb, shiftedXMultdb, polyTblhdb, polyTblmdb;
  int xIntHi, expoX, k, M, index1, index2;
  double highPoly, tt1h, t1h, t1l, xSqh, xSql, xSqHalfh, xSqHalfl, xCubeh, xCubel, t2h, t2l, templ, tt3h, tt3l;
  double polyh, polyl, expm1h, expm1m, expm1l;
  double r1h, r1l, r1t, rr1h, rr1l;
  double r2h, r2l, r2t, rr2h, rr2l;
  double r3h, r3l, r3t, rr3h, rr3l;
  double xMultLog2InvMult2L, shiftedXMult, kd, s1, s2, s3, s4, s5, rh, rm, rl;
  double rhSquare, rhC3, rhSquareHalf, monomialCube, rhFour, monomialFour;
  double tbl1h, tbl1m, tbl1l, tbl2h, tbl2m, tbl2l;
  double highPolyWithSquare, tablesh, tablesl, t8, t9, t10, t11, t12, t13;
  double exph, expm, t1, t2, t3;
  double msLog2Div2LMultKh, msLog2Div2LMultKm, msLog2Div2LMultKl;


  xdb.d = x; 

  /* Strip off the sign of x for the following tests */

  xIntHi = xdb.i[HI] & 0x7fffffff;

  /* Test if we are so small that we can return (a corrected) x as correct rounding */
  if (xIntHi < RETURNXBOUND) {
    /* The only algebraic result is 0 for x = 0; in this case, we can return x
       expm1 is positive for positive x, negative for negative x
       The truncation rest x^2/2 + x^3/6 + ... is always positive 
       but less than 1 ulp in this case, so we round as follows:
       
       - x is positive => expm1 is positive => round downwards => truncate by returning x
       - x is negative => expm1 is negative => round upwards => add 1 ulp
    */
    if (xdb.i[HI] & 0x80000000) {
      /* x is negative 
	 We add 1 ulp by subtracting 1 in long
      */
      xdb.l--;
      return xdb.d;
    } else {
      /* x is positive 
	 We do nothing (see above)
      */
      return x;
    }
  }


  /* Filter out special cases like overflow, -1 in result, infinities and NaNs 
     The filters are not sharp, we have positive arguments that flow through
  */
  if (xIntHi >= SIMPLEOVERFLOWBOUND) {
    /* Test if we are +/-inf or NaN */
    if (xIntHi >= 0x7ff00000) {
      /* Test if NaN */
      if (((xIntHi & 0x000fffff) | xdb.i[LO]) != 0) {
	/* NaN */
	return x+x;  /* return NaN */
      }
      /* Test if +inf or -inf */
      if (xdb.i[HI] > 0) {
	/* +inf */
	return x+x;  /* return +inf */
      }
      
      /* If we are here, we are -inf */
      return -1.0;
    }

    /* If we are here, we are overflowed or a common case that flows through */

    /* Test if we are actually overflowed */
    if (x > OVERFLOWBOUND) {
      /* We would be overflowed but as we are rounding towards zero, i.e. downwards,
	 the nearest number lesser than the exact result is the greatest 
	 normal. In any case, we must raise the inexact flag.
      */
      return LARGEST * (1.0 + SMALLEST);
    }
  }
  
  /* Test if we know already that we are -1.0 (+ correction depending on rounding mode) in result */
  if (x < MINUSONEBOUND) {
    /* We round towards zero, i.e. upwards, so we return -1.0+1ulp */
    return MINUSONEPLUSONEULP;
  }

  /* Test if we have |x| <= 1/4-1/2ulp(1/4) for knowing if we use exp(x) or approximate directly */

  if (xIntHi < DIRECTINTERVALBOUND) {
    /* We approximate expm1 directly after a range reduction as follows

       expm1(x) = (expm1(x/2) + 2) * expm1(x/2)

       We perform the range reduction in such a way that finally |x| < 1/32 
    */

    /* Extract the exponent of |x| and add 5 (2^5 = 32) */
    expoX = ((xIntHi & 0x7ff00000) >> 20) - (1023 - 5);
    
    /* If this particularily biased exponent expoX is negative, we are already less than 1/32 */
    if (expoX >= 0) {
      /* If we are here, we must perform range reduction */


      /* We start by producing 2^(-expoX-1) as a factor to x */
      scaledb.i[HI] = ((-expoX-1)+1023) << 20;
      scaledb.i[LO] = 0;
      
      /* We multiply x with the scale */
      x = scaledb.d * x;
    }
    
    /* Here, we have always |x| < 1/32 */


    /* Double precision evaluation steps */
#if defined(PROCESSOR_HAS_FMA) && !defined(AVOID_FMA)
    highPoly = FMA(FMA(FMA(FMA(FMA(quickDirectpolyC9h ,x,quickDirectpolyC8h),x,quickDirectpolyC7h),x,
                                   quickDirectpolyC6h),x,quickDirectpolyC5h),x,quickDirectpolyC4h);
#else
    highPoly = quickDirectpolyC4h + x * (quickDirectpolyC5h + x * (quickDirectpolyC6h + x * (
	       quickDirectpolyC7h + x * (quickDirectpolyC8h + x *  quickDirectpolyC9h))));
#endif

    /* Double-double evaluation steps */
    tt1h = x * highPoly;

    Mul12(&xSqh,&xSql,x,x);
    xSqHalfh = 0.5 * xSqh;
    xSqHalfl = 0.5 * xSql;
    Add12(t2h,templ,x,xSqHalfh);
    t2l = templ + xSqHalfl;
    
    Add12(t1h,t1l,quickDirectpolyC3h,tt1h);
    Mul122(&xCubeh,&xCubel,x,xSqh,xSql);
    Mul22(&tt3h,&tt3l,xCubeh,xCubel,t1h,t1l);

    Add22(&polyh,&polyl,t2h,t2l,tt3h,tt3l);

    /* Reconstruction */

    /* If we have not performed any range reduction, we have no reconstruction to do */
    if (expoX >= 0) {
      /* If we are here, we must perform reconstruction */

      /* First reconstruction step */
      Add12(r1h,r1t,2,polyh);
      r1l = r1t + polyl;
      Mul22(&rr1h,&rr1l,r1h,r1l,polyh,polyl);

      if (expoX >= 1) {

	/* Second reconstruction step */
	Add12(r2h,r2t,2,rr1h);
	r2l = r2t + rr1l;
	Mul22(&rr2h,&rr2l,r2h,r2l,rr1h,rr1l);

	if (expoX >= 2) {

	  /* Third reconstruction step */
	  Add12(r3h,r3t,2,rr2h);
	  r3l = r3t + rr2l;
	  Mul22(&rr3h,&rr3l,r3h,r3l,rr2h,rr2l);
	  
	  /* expoX may be maximally 2 */

	  expm1h = rr3h;
	  expm1m = rr3l;

	} else {
	  expm1h = rr2h;
	  expm1m = rr2l;
	}

      } else {
	expm1h = rr1h;
	expm1m = rr1l;
      }

    } else {
      expm1h = polyh;
      expm1m = polyl;
    }

    /* Rounding test */
    TEST_AND_RETURN_RZ(expm1h, expm1m, ROUNDCSTDIRECTRD);
    {
      expm1_direct_td(&expm1h, &expm1m, &expm1l, x, xSqHalfh, xSqHalfl, xSqh, xSql, expoX);
      
      ReturnRoundTowardsZero3(expm1h, expm1m, expm1l);

    } /* Accurate phase launched */

    /* We cannot be here, since we return in all cases before */
  }

  /* If we are here, we can use expm1(x) = exp(x) - 1 */

  /* Range reduction - exact part: compute k as double and as int */

  xMultLog2InvMult2L = x * log2InvMult2L;
  shiftedXMult = xMultLog2InvMult2L + shiftConst;
  kd = shiftedXMult - shiftConst;
  shiftedXMultdb.d = shiftedXMult;
  k = shiftedXMultdb.i[LO];
  M = k >> 12;
  index1 = k & INDEXMASK1;
  index2 = (k & INDEXMASK2) >> 6;

  /* Range reduction - part affected by error - must be redone in accurate phase */
  Mul12(&s1,&s2,msLog2Div2Lh,kd);
  s3 = kd * msLog2Div2Lm;
  s4 = s2 + s3; 
  s5 = x + s1;
  Add12Cond(rh,rm,s5,s4);

  /* Table reads - read only two double-doubles by now */
  tbl1h = twoPowerIndex1[index1].hi;
  tbl1m = twoPowerIndex1[index1].mi;
  tbl2h = twoPowerIndex2[index2].hi;
  tbl2m = twoPowerIndex2[index2].mi;

  /* Quick phase starts here */

  rhSquare = rh * rh;
  rhC3 = quickCommonpolyC3h * rh;

  rhSquareHalf = 0.5 * rhSquare;
  monomialCube = rhC3 * rhSquare;
  rhFour = rhSquare * rhSquare;

  monomialFour = quickCommonpolyC4h * rhFour;
  
  highPoly = monomialCube + monomialFour;

  highPolyWithSquare = rhSquareHalf + highPoly;
  
  /* Reconstruction: integration of table values */
  
  Mul22(&tablesh,&tablesl,tbl1h,tbl1m,tbl2h,tbl2m);

  t8 = rm + highPolyWithSquare;
  t9 = rh + t8;

  t10 = tablesh * t9;
  
  Add12(t11,t12,tablesh,t10);
  t13 = t12 + tablesl;
  Add12(polyTblhdb.d,polyTblmdb.d,t11,t13);
  
  /* Reconstruction: multiplication by 2^M */

  /* Implement the multiplication by multiplication to overcome the
     problem of the non-representability of 2^1024 (M = 1024)
     This case is possible if polyTblhdb.d < 1
  */
  
  polyTblhdb.i[HI] += M << 20;
  polyTblmdb.i[HI] += M << 20;

  exph = polyTblhdb.d;
  expm = polyTblmdb.d;

  /* Substraction of 1 

     Testing if the operation is necessary is more expensive than 
     performing it in any case.

     We may cancellate at most 2 bits in the subtraction for 
     arguments 1/4 <= x <= ln(2) (0.25 <= x <= 0.69) 
     We must therefore use conditional Add12s

     Since we perform a substraction, we may not have addition overflow towards +inf

  */

  Add12Cond(t1,t2,-1,exph);
  t3 = t2 + expm;
  Add12Cond(expm1h,expm1m,t1,t3);


  /* Rounding test */
  TEST_AND_RETURN_RZ(expm1h, expm1m, ROUNDCSTCOMMONRD);
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
    expm1_common_td(&expm1h, &expm1m, &expm1l, rh, rm, rl, tbl1h, tbl1m, tbl1l, tbl2h, tbl2m, tbl2l, M); 
    
    /* Final rounding */

    ReturnRoundTowardsZero3(expm1h, expm1m, expm1l);
  } /* Accurate phase launched */
  
  /* We cannot be here since we return before in any case */
}
