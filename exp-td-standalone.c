#include "crlibm.h"
#include "crlibm_private.h"
#include "triple-double.h"
#include "exp-td.h"

#define AVOID_FMA 0


/* Function exp13

   Computes exp(x) with an accuracy of 113 bits as

   2^exponent * (exph + expm + expl) \approx exp(x)

   Unless the subnormal case for x, no special cases are 
   handled.

   The triple-double exph + expm + expl is non-overlapping.
   The domain for exph + expm + expl is 1/2..2
   The integer exponent is in the range -1024..1024. The 
   value 2^(exponent) may therefore be non-representable
   whereas 2^exponent * (exph + expm + expl) is.

*/


void exp13(int *exponent, double *exph, double *expm, double *expl, double x) { 
  double rh, rm, rl, tbl1h, tbl1m, tbl1l;
  double tbl2h, tbl2m, tbl2l;
  double xMultLog2InvMult2L, shiftedXMult, kd;
  double msLog2Div2LMultKh, msLog2Div2LMultKm, msLog2Div2LMultKl;
  double t1, t2;
  db_number shiftedXMultdb, xdb;
  int k, M, index1, index2, xIntHi;
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
  double polyTblh, polyTblm, polyTbll;
   
  /* Argument reduction and filtering for special cases */

  /* Compute k as a double and as an int */
  xdb.d = x;
  xMultLog2InvMult2L = x * log2InvMult2L;
  shiftedXMult = xMultLog2InvMult2L + shiftConst;
  kd = shiftedXMult - shiftConst;
  shiftedXMultdb.d = shiftedXMult;
  
  /* Special cases tests */
  xIntHi = xdb.i[HI];
  /* Test if argument is a denormal or zero */
  if ((xIntHi & 0x7ff00000) == 0) {
    /* We are in the RN case, return 1.0 in all cases */
    *exph = 1.0;
    *expm = 0.0;
    *expl = 0.0;
    return;
  }
 
  k = shiftedXMultdb.i[LO];
  M = k >> L;
  index1 = k & INDEXMASK1;
  index2 = (k & INDEXMASK2) >> LHALF;

  /* Table reads */
  tbl1h = twoPowerIndex1[index1].hi;
  tbl1m = twoPowerIndex1[index1].mi;
  tbl2h = twoPowerIndex2[index2].hi;
  tbl2m = twoPowerIndex2[index2].mi;
  tbl1l = twoPowerIndex1[index1].lo;
  tbl2l = twoPowerIndex2[index2].lo;

  /* Argument reduction */

  Mul133(&msLog2Div2LMultKh,&msLog2Div2LMultKm,&msLog2Div2LMultKl,kd,msLog2Div2Lh,msLog2Div2Lm,msLog2Div2Ll);
  t1 = x + msLog2Div2LMultKh;
  Add12Cond(rh,t2,t1,msLog2Div2LMultKm);
  Add12Cond(rm,rl,t2,msLog2Div2LMultKl);

  /* Polynomial approximation */

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

  Renormalize3(&polyTblh,&polyTblm,&polyTbll,polyWithTablesh,polyWithTablesm,polyWithTablesl);

  *exponent = M;
  *exph = polyTblh;
  *expm = polyTblm;
  *expl = polyTbll;
}
