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
  double q,x2;
  double xhi,xlo;
  
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
  
  if (x > my_e) /* test if reduction is necessary : */
    {
      double xmBIhi, xmBIlo;
      
      if (x > value[61][B].d) {
        i=61;
        Add12( xmBIhi , xmBIlo , x , -value[61][B].d);
      }
      else {
        /* determine i so that x E [a[i],a[i+1]] */
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
        else if (i<61) i+= 1;
        if (x < value[i][A].d) i-= 1;
          
        xmBIhi = x-value[i][B].d;
        xmBIlo = 0.0;
      }
        
        
      /* we now compute X = ( x-b[i] ) / ( 1 + x*b[i] )
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

      if (tmphi > 1)
        Add22(&x0hi,&x0lo,tmphi,tmplo, 1.0,0.0);
      else {Add22( &x0hi , &x0lo , 1.0,0.0,tmphi,tmplo);}
      
      DIV2( xmBIhi , xmBIlo , x0hi,x0lo, xhi,xlo);
      
      /* Polynomial evaluation : 
       *  
       *  1rt compute Q(x^2) = (1 - x^2/3 + ...)
       *      then P(x) = x * Q(x^2)
       *
       */

      x2 = xhi*xhi;
      
#if 0
      q = coef_poly[0];
      int j=1;
      /* a derouler ? */
      for (;j<DEGREE;j++) {
        q *= x2;
        q += coef_poly[j];
      }
      q*=x2;
#else
      q = x2*(coef_poly[3]+x2*(coef_poly[2]+x2*(coef_poly[1]+x2*coef_poly[0]))) ;
#endif 
      double atanXhi,atanXlo;
      
      Add12(atanXhi,atanXlo,xhi,( (xhi*q)+xlo));
      
      /* reconstruction : atan(x) = atan(b[i]) + atan(x) */
      Add22 (&atanhi,&atanlo, value[i][ATAN_BHI].d,value[i][ATAN_BLO].d,atanXhi,atanXlo);  
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
#if 0
      q = coef_poly[0];
      int j=1;
      for (;j<DEGREE;j++) 
        {q *= x2;
        q += coef_poly[j];}
      q*=x2;
#else
      q = x2*(coef_poly[3]+x2*(coef_poly[2]+x2*(coef_poly[1]+x2*coef_poly[0]))) ;
#endif

      Add12(atanhi,atanlo, x , x*q);
    }
  
  /* test if rounding is possible */
  if (atanhi == (atanhi + (atanlo*e)))
    {
      return atanhi;}
  else
    {/* more accuracy is needed , lauch accurate phase */ 
      return scs_atan_rn(x);
    }
}
