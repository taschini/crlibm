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


/*
 * In the following, when an operator is preceded by a '@' it means that we
 * are considering the IEEE-compliant machine operator, otherwise it
 * is the mathematical operation.
 *
 */



#if  DEKKER_AS_FUNCTIONS 
/* else it is defined in crlibm_private.h */

/*
 * computes rh and rl such that rh + rl = a * b with rh = a @* b exactly
 * under the conditions : a < 2^970 et b < 2^970 
 */
#ifndef CRLIBM_TYPECPU_ITANIUM
/* otherwise Mul12 and Mul12cond are #defined using FMS in crlibm_private.h  */
void  Mul12(double *rh, double *rl, double u, double v){
  const double c = 134217729.;   /*  1+2^27 */ 
  double up, u1, u2, vp, v1, v2;

  up = u*c;        vp = v*c;
  u1 = (u-up)+up;  v1 = (v-vp)+vp;
  u2 = u-u1;       v2 = v-v1;
  
  *rh = u*v;
  *rl = (((u1*v1-*rh)+(u1*v2))+(u2*v1))+(u2*v2);
}

/*
 * Computes rh and rl such that rh + rl = a * b and rh = a @* b exactly
 */
void Mul12Cond(double *rh, double *rl, double a, double b){
  const double two_970 = 0.997920154767359905828186356518419283e292;
  const double two_em53 = 0.11102230246251565404236316680908203125e-15;
  const double two_e53  = 9007199254740992.;
  double u, v;

  if (a>two_970)  u = a*two_em53; 
  else            u = a;
  if (b>two_970)  v = b*two_em53; 
  else            v = b;

  Mul12(rh, rl, u, v);

  if (a>two_970) {*rh *= two_e53; *rl *= two_e53;} 
  if (b>two_970) {*rh *= two_e53; *rl *= two_e53;} 
}
#endif /*CRLIBM_TYPECPU_ITANIUM*/



/*
 * computes double-double addition: zh+zl = xh+xl + yh+yl
 * relative error is smaller than 2^-103 
 */
  
void Add22(double *zh, double *zl, double xh, double xl, double yh, double yl)
{
double r,s;

r = xh+yh;
s = (ABS(xh) > ABS(yh))? (xh-r+yh+yl+xl) : (yh-r+xh+xl+yl);
*zh = r+s;
*zl = r - (*zh) + s;
}



/*
 * computes double-double multiplication: zh+zl = (xh+xl) *  (yh+yl)
 * relative error is smaller than 2^-102
 */
  

  
void Mul22(double *zh, double *zl, double xh, double xl, double yh, double yl)
{
double mh, ml;

  const double c = 134217729.;                /* 0x41A00000, 0x02000000 */ 
  double up, u1, u2, vp, v1, v2;

  up = xh*c;        vp = yh*c;
  u1 = (xh-up)+up;  v1 = (yh-vp)+vp;
  u2 = xh-u1;       v2 = yh-v1;
  
  mh = xh*yh;
  ml = (((u1*v1-mh)+(u1*v2))+(u2*v1))+(u2*v2);

  ml += xh*yl + xl*yh;
  *zh = mh+ml;
  *zl = mh - (*zh) + ml;
}

#endif  /*DEKKER_AS_FUNCTIONS */
