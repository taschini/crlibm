/**
 * Variables and common functions shared by many functions

 * Author : Defour David  (David.Defour@ens-lyon.fr)
 *	    Daramy Catherine (Catherine.Daramy@ens-lyon.fr)
 *          Florent de Dinechin (Florent.de.Dinechin@ens-lyon.fr)
 * Date of creation : 14/03/2002   
 * Last Modified    : 28/02/2003
 */
#include <stdio.h>
#include <stdlib.h>
#include "crlibm.h"
#include "crlibm_private.h"
#ifdef HAVE_FENV_H
#include <fenv.h>
#endif
#ifdef CRLIBM_TYPECPU_X86
#include <fpu_control.h>
#endif

/* An init function which sets FPU flags when needed */
void crlibm_init() {
#ifdef CRLIBM_TYPECPU_X86
  /* Set FPU flags to use double, not double extended, 
     with rounding to nearest */
  unsigned int cw = (_FPU_DEFAULT & ~_FPU_EXTENDED)|_FPU_DOUBLE;
  _FPU_SETCW(cw);
#else

#endif
}


double ABS(double x){
  double y;
  
  y=(x>0) ? x : -x;
  return y;
}

/*
 * In the following, when an operator is preceded by a '@' it means that we
 * are considering the IEEE-compliant machine operator, otherwise it
 * is the mathematical operator.
 *
 */


#if  DEKKER_INLINED 
/*
 * computes r1 and r2 such that r1 + r2 = a * b with r1 = a @* b exactly
 * under the conditions : a < 2^970 et b < 2^970 
 */
void  Mul12(double *r1, double *r2, double u, double v){
  double c        = 134217729.;                /* 0x41A00000, 0x02000000 */ 
  double up, u1, u2, vp, v1, v2;

  up = u*c;        vp = v*c;
  u1 = (u-up)+up;  v1 = (v-vp)+vp;
  u2 = u-u1;       v2 = v-v1;
  
  *r1 = u*v;
  *r2 = (((u1*v1-*r1)+(u1*v2))+(u2*v1))+(u2*v2);
}

/*
 * Computes r1 and r2 such that r1 + r2 = a * b and r1 = a @* b exactly
 */
void Mul12Cond(double *r1, double *r2, double a, double b){
  double two_em53 = 1.1102230246251565404e-16; /* 0x3CA00000, 0x00000000 */
  double two_e53  = 9007199254740992.;         /* 0x43400000, 0x00000000 */
  double u, v, a, b;

  if (HI(a)>0x7C900000) u = a*two_em53; 
  else            u = a;
  if (HI(b)>0x7C900000) v = b*two_em53; 
  else            v = b;

  Mul12(r1, r2, u, v);

  if (HI(a)>0x7C900000) {*r1 *= two_e53; *r2 *= two_e53;} 
  if (HI(b)>0x7C900000) {*r1 *= two_e53; *r2 *= two_e53;} 
}

#endif  /*DEKKER_INLINED */

/* /\* */
/*  *computes the exact product of x and y. Return a double-double result z, zz */
/*  *\/ */

/* void mul12(double *z, double *zz, double x, double y){ */
/* static const double c = 134217729.; */
/* double hx, tx, hy, ty, p, q; */

/* p = x*c; */
/* hx = (x-p)+p; */
/* tx = x-hx; */
/* p = y*c; */
/* hy = (y-p)+p; */
/* ty = y-hy; */
/* p = hx*hy; */
/* q = hx*ty + tx*hy; */
/* *z = p+q; */
/* *zz = p - (*z) + q + tx*ty; */
/* } */

/*
 * Computes double-double multiplication z+zz = (x+xx)*(y+yy) 
 * under the conditions : x < 2^970 et y < 2^970 
 */
  
void Mul22(double * z, double * zz, double x, double xx, double y, double yy){
double mc, mcc;

Mul12(&mc, &mcc, x, y);
mcc = x*yy + xx*y + mcc;
*z = mc+mcc;
*zz = mc - (*z) + mcc;
}  
  

/*
 * computes double-double addition: z+zz = x+xx + y+yy
 */
  
void Add22(double * z, double * zz, double x, double xx, double y, double yy){
double r,s;

r = x+y;
s = (ABS(x) > ABS(y))? (x-r+y+yy+xx) : (y-r+x+xx+yy);
*z = r+s;
*zz = r - (*z) + s;
}



