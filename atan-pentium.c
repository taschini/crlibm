/* 
 *this function computes atan correctly rounded to the nearest, 
 using experimental techniques based on  double-extended arithmetic

 THIS IS EXPERIMENTAL SOFTWARE

In particular it changes rounding modes all the time without warning
nor restoring.
 
 *
 * Author : Nicolas Gast, Florent de Dinechin
 * nicolas.gast@ens.fr
 *

 To have it replace the crlibm atan, do:
 gcc -DHAVE_CONFIG_H -I.  -fPIC  -O2 -c atan-pentium.c;   mv atan-pentium.o atan_fast.o; make 
 
*/


#include <stdio.h>
#include <stdlib.h>
#include <crlibm.h>
#include <crlibm_private.h>
#include "double-extended.h"

#define DEBUG 0
#define NICOLASTEST 0

#ifdef HAVE_FENV_H
#include <fenv.h>
#endif


/* The following seems perfectly harmless */

#ifdef FENV_H
#pragma STDC FENV_ACCESS ON
#endif

#include "atan-pentium.h"
#include <fpu_control.h>



/* Dummy functions to compile OK */
extern double atan_rd(double x) {return 0;}
extern double atan_ru(double x) {return 0;}
extern double atan_rz(double x) {return 0;}


extern double atan_rn(double x) {
  db_number x_db;
  unsigned int hx;
  double sign;
  double u;
  double comp;
  double atanhi, atanlo, atanlo_u;

  long double Xred;
  long double Xred2;
  long double q;
  long double atan;
  long double eps;
  int i;
  
  if(x>=0)
    sign = 1;
  else
    {sign = -1;
    x=-x;}
  
  x_db.d = x;
  hx = x_db.i[HI] & 0x7FFFFFFF; 
  
  /* Filter cases */
  if ( hx >= 0x43500000)           /* x >= 2^54 */
    {
      if ( (x_db.i[LO] == 0) && (hx & 0x000fffff) == 0x00080000)
        return x+x;                /* NaN */
      else
        return sign*HALFPI.d;           /* atan(x) = Pi/2 */
    }
  else
    if ( hx < 0x3E400000 )
      {return sign*x;}                   /* x<2^-27 then atan(x) =~ x */

   DOUBLE_EXTENDED_MODE;


  if (x > MIN_REDUCTION_NEEDED) /* test if reduction is necessary : */
    {
      /* 1) Argument reduction :  */
      
      /* compute i so that a[i] < x < a[i+1] */
      if (x>arctan_table[61][A])
        i=61;
      else {
        i=31;
        if (x < arctan_table[i][A]) i-= 16;
        else i+=16;
        if (x < arctan_table[i][A]) i-= 8;
        else i+= 8;
        if (x < arctan_table[i][A]) i-= 4;
        else i+= 4;
        if (x < arctan_table[i][A]) i-= 2;
        else i+= 2;
        if (x < arctan_table[i][A]) i-= 1;
        else i+= 1;
        if (x < arctan_table[i][A]) i-= 1;
      }
      Xred = (x - arctan_table[i][B] )/(1.0L +  x * arctan_table[i][B] );
      
      Xred2 = Xred*Xred;
      /* Polynomial evaluation */
      q = Xred2*(coef_poly[0][0]+Xred2*
           (coef_poly[1][0]+Xred2*
            (coef_poly[2][0]+Xred2*
             (coef_poly[3][0]))));
      
      /* reconstruction : atan(x) = atan(b[i]) + atan(x) */
      atan = arctan_table[i][ATAN_BHI] + (Xred + q*Xred);
      

      atanhi = (double) atan;
      atanlo = atan-atanhi;
      
      
    }
  else
    
    // no reduction needed
    {
      
      Xred2 = x*x;
      
      /* Polynomial evaluation */
      q = Xred2*(coef_poly[0][0]+Xred2*
           (coef_poly[1][0]+Xred2*
            (coef_poly[2][0]+Xred2*
             (coef_poly[3][0]))));
      
      atan = q*x + x;
    
    }
  
  TEST_AND_RETURN_RN2(atan, sign*atan, 0x7fe);
  /* or : 
  TEST_AND_RETURN_RN_ZIV(sign*atan, 1.003);
  */
   {

    /*Second step, double-double  */
    long double tmphi, tmplo;
    long double x0hi, x0lo;
    long double xmBihi, xmBilo;
    long double Xredhi, Xredlo;
    long double Xred2;
    long double qhi,qlo; /* q = polynomial */
    long double q;
    long double Xred2hi,Xred2lo;
    long double atanhi,atanlo;
    int j;
  
#if EVAL_PERF
  crlibm_second_step_taken++;
#endif
  
  if (x > MIN_REDUCTION_NEEDED) /* test if reduction is necessary : */
    {
      /* 1) Argument reduction :  */
      
      if (i==61) 
        {
          Add12_ext( xmBihi , xmBilo , x , -arctan_table[61][B]);
        }
      else 
        {
          xmBihi = x-arctan_table[i][B];
          xmBilo = 0.0;
        }
      
      Mul12_ext(&tmphi,&tmplo, x, arctan_table[i][B]);
      
      if (x > 1)
        Add22_ext(&x0hi,&x0lo,tmphi,tmplo, 1.0,0.0);
      else {Add22_ext( &x0hi , &x0lo , 1.0,0.0,tmphi,tmplo);}
      
      Div22_ext( Xredhi, Xredlo, xmBihi , xmBilo , x0hi,x0lo);
      
      Xred2 = Xredhi*Xredhi;
      Mul22_ext( &Xred2hi,&Xred2lo,Xredhi,Xredlo,Xredhi, Xredlo);
      
      /*poly eval */
      
      q = (coef_poly[4][0]+Xred2*
           (coef_poly[5][0]+Xred2*
            (coef_poly[6][0]+Xred2*
             (coef_poly[7][0]+
              (Xred2*coef_poly[8][0])))));
      
      Mul12_ext( &qhi, &qlo, q, Xred2);
      
      for(j=3;j>=0;j--)
        {
          Add22_ext(&qhi,&qlo, coef_poly[j][0], coef_poly[j][1], qhi,qlo);
          Mul22_ext(&qhi,&qlo, qhi,qlo, Xred2hi,Xred2lo);
        }
      
      Mul22_ext(&qhi,&qlo, Xredhi,Xredlo, qhi,qlo);
      Add22_ext(&qhi,&qlo, Xredhi,Xredlo, qhi,qlo);
      
      /* reconstruction : atan(x) = atan(b[i]) + atan(x) */
      Add22_ext(&atanhi,&atanlo, arctan_table[i][ATAN_BHI], arctan_table[i][ATAN_BLO], qhi,qlo);
    }
  else
    
    // no reduction needed
    {
      /* Polynomial evaluation */
      Mul12_ext( &Xred2hi,&Xred2lo,x,x);
      
      /*poly eval */
      q = Xred2hi*(coef_poly[5][0]+Xred2hi*
		   (coef_poly[6][0]+Xred2hi*
		    (coef_poly[7][0]+Xred2hi*
		     (coef_poly[8][0]))));
      
      Add12_ext(qhi,qlo, coef_poly[4][0],q);
      Mul22_ext(&qhi,&qlo, qhi,qlo, Xred2hi,Xred2lo);
      
      for(j=3;j>=0;j--)
        {
          Add22_ext(&qhi,&qlo, coef_poly[j][0], coef_poly[j][1], qhi,qlo);
          Mul22_ext(&qhi,&qlo, qhi,qlo, Xred2hi,Xred2lo);
        }
      
      Mul22_ext (&qhi,&qlo, x,0, qhi,qlo);
      Add12_ext (atanhi,atanlo,x,qhi);
      atanlo += qlo;
    }
#define debug 0
#if debug
  printf("             %1.50Le + %1.50Le\n",atanhi, atanlo);
  printf("             %1.50e + %1.50e\n",(double)atanhi,(double) atanlo);
  printf("             %1.50Le\n",atanhi + atanlo);
#endif
  
   BACK_TO_DOUBLE_MODE;
  return sign*((double) (atanhi+atanlo));
    }
}
