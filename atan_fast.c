/*
 * Function to compute atan with fully exact rounding
 *
 * Author : Gast Nicolas
 * nicolas.gast@ens.fr
 *
 * Date : 24/06/2004
 */

#include <stdio.h>
#include <stdlib.h>
#include <crlibm.h>
#include <crlibm_private.h>
#include <atan_fast.h>

#define  DIV2(x,xx,y,yy,z,zz)\
do \
{double __c,__cc,__u,__uu;  \
           __c=(x)/(y);   Mul12(&__u,&__uu,__c,y);  \
           __cc = (((((x)-__u)-__uu)+(xx))-__c*(yy))/(y);   z=__c+__cc;   zz=(__c-z)+__cc;\
} \
while(0)

/* the second step :   A VOIR PLUS TARD */ 
double scs_atan_rn(double); 

/*
 * 1) Argument reduction : 
 * 
 *  tan(x) = tan( b(i) ) + tan ( (x-b(i)) / (1+x*b(i)))
 *
 *                                                     6.3
 * we choose 63 b(i) so that (x-b(i)) / (1+x*b(i)) < 2^
 */ 

static double atan_rn2(double x);

extern double atan_rn (double x) {
  if (x>=0)
    return atan_rn2(x);
  else
    return -atan_rn2(-x);
}

static double atan_rn2 (double x) {  
  
  double atanhi,atanlo;
  double tmphi,tmplo, x0hi,x0lo;
  double q,Xred2,x2;
  double Xredhi,Xredlo;
  
  
  int i;
  db_number x_db;
  x_db.d = x;
  unsigned int hx = x_db.i[HI_ENDIAN] & 0x7FFFFFFF; 
  
  /* Filter cases */
  if ( hx >= 0x43500000)           /* x >= 2^54 */
    {
      if ( ( (hx & 0x000fffff) | x_db.i[LO_ENDIAN] ) == 0)
        return x+x;                /* NaN */
      else
        return HALFPI.d;           /* atan(x) = Pi/2 */
    }
  else
    if ( hx < 0x3E400000 )
      {return x;                   /* x<2^-27 then atan(x) =~ x */}
  
  if (x > MIN_REDUCTION_NEEDED) /* test if reduction is necessary : */
    {
      double xmBihi, xmBilo;
      
      if (x > value[61][B].d) {
        i=61;
        Add12( xmBihi , xmBilo , x , -value[61][B].d);
      }
      else {
        /* determine i so that a[i] < x < a[i+1] */
        i=31;
        if (x < value[i][A].d) i-= 16;
        else i+=16;
        if (x < value[i][A].d) i-= 8;
        else i+= 8;
        if (x < value[i][A].d) i-= 4;
        else i+= 4;
        if (x < value[i][A].d) i-= 2;
        else i+= 2;
        if (x < value[i][A].d) i-= 1;
        else i+= 1;
        if (x < value[i][A].d) i-= 1;     
        xmBihi = x-value[i][B].d;
        xmBilo = 0.0;
      }
        
      /* we now compute Xred = ( x-b[i] ) / ( 1 + x*b[i] )
       * 
       * def : x0 := 1+x*b[i]
       *
       * 1st we compute an approximation of y = 1/x0
       * then we compute a better approx x' = y*(2-x0*y)
       * we can proove that :
       * if y = 1/x0*(1+e) 
       *     then x' = 1/x0 * (1-e^2)
       *                   
       */
      
      Mul12(&tmphi,&tmplo, x, value[i][B].d);

      if (x > 1)
        Add22(&x0hi,&x0lo,tmphi,tmplo, 1.0,0.0);
      else {Add22( &x0hi , &x0lo , 1.0,0.0,tmphi,tmplo);}

      DIV2( xmBihi , xmBilo , x0hi,x0lo, Xredhi,Xredlo);

      /* Polynomial evaluation : 
       *  
       *  1rt compute Q(x^2) = (1 - x^2/3 + ...)
       *      then P(x) = x * Q(x^2)
       *
       */

      Xred2 = Xredhi*Xredhi;
      
      q = Xred2*(coef_poly[3]+Xred2*
                 (coef_poly[2]+Xred2*
                  (coef_poly[1]+Xred2*
                   coef_poly[0]))) ;

      /* reconstruction : atan(x) = atan(b[i]) + atan(x) */
      double testlo = Xredlo+ value[i][ATAN_BLO].d + Xredhi*q;
      double tmphi2, tmplo2;
      Add12( tmphi2, tmplo2, value[i][ATAN_BHI].d, Xredhi);
      Add12( atanhi, atanlo, tmphi2, (tmplo2+testlo));

    }
  else 
    // no reduction needed
    {
      /* Polynomial evaluation : 
       *  
       *  1rt compute Q(x^2) = (1 - x^2/3 + ...)
       *      then P(x) = x * Q(x^2)
       *
       */
      
      x2 = x*x;
      q = x2*(coef_poly[3]+x2*
                 (coef_poly[2]+x2*
                  (coef_poly[1]+x2*
                   coef_poly[0]))) ;
      Add12(atanhi,atanlo, x , x*q);
      
      if ( (atanhi == (atanhi + (atanlo*e_no_reduction)))
           || ( (hx <= 0x3F500000 ) && 
                atanhi == (atanhi + (atanlo*e_no_reduction_m10) ) ) )
        /* there is two cases : x > 2^-10 and x < 2^-10 */
        {
          return atanhi;}
      else
        {/* more accuracy is needed , lauch accurate phase */ 
          return scs_atan_rn(x);
        }
    }
  
  /* test if rounding is possible */
  if (atanhi == (atanhi + (atanlo*e)) 
      || ( i>=10 && (atanhi == (atanhi + (atanlo*e_i_10)) )))
    {
      return atanhi;}
  else
    {
      /* more accuracy is needed , lauch accurate phase */ 
      return scs_atan_rn(x);
    }
}
