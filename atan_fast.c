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
#include "atan_fast.h"

#define  DIV2(x,xx,y,yy,z,zz)\
{double c,cc,u,uu;  \
           c=(x)/(y);   Mul12(&u,&uu,c,y);  \
           cc=(((((x)-u)-uu)+(xx))-c*(yy))/(y);   z=c+cc;   zz=(c-z)+cc;\
} 

//#define dprintf(format, args...)  printf(format , ## args)
#define dprintf(format, args...) 

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

  
  dprintf("\n\n");
  
  /* first : test special cases : */
  if ( hx >= 0x43500000)
    { /* x >= 2^54 */
      if ( ( (hx & 0x000fffff) | x_db.i[LO_ENDIAN] ) == 0)
        return x+x;                /* NaN */
      else
        return HALFPI.d;           /* atan(x) = Pi/2 */
    }
  else
    if ( hx < 0x3E400000 )
      {return x; /* x<2^-27 then atan(x) =~ x */}
  
  /* test if reduction is necessary : */
  if (x > my_e)
    {
      double xmBIhi, xmBIlo;
      
      // a amelio ?
      if (x > b[61].d) {
        i=61;
        Add12( xmBIhi , xmBIlo , x , -b[61].d);
      }
      else {
        /* determine i so that x E [a[i],a[i+1]] */
        i=31;
        if (x < a[i]) i-= 16;
        else i+=16;
        if (x < a[i]) i-= 8;
        else i+= 8;
        if (x < a[i]) i-= 4;
        else i+= 4;
        if (x < a[i]) i-= 2;
        else i+= 2;
        if (x < a[i]) i-= 1;
        else if (i<61) i+= 1;
        if (x < a[i]) i-= 1;
        
        xmBIhi = x-b[i].d;
        xmBIlo = 0.0;
      }
      
      
      dprintf("i = %d (b[i] = %f,%f<%f<%f = )\n",
              i,b[i].d,a[i], x ,a[i+1]);
      
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
      
            
      Mul12(&tmphi,&tmplo, x, b[i].d);
      dprintf("xfoisb[i] := nearest(%1.50e)+nearest(%1.50e)\n",tmphi,tmplo);
      if (tmphi > 1)Add22(&x0hi,&x0lo,tmphi,tmplo, 1.0,0.0);   /* ""          */
      else Add22(&x0hi,&x0lo, 1.0,0.0,tmphi,tmplo);
      dprintf("x0:          nearest(%1.50e)+nearest(%1.50e) \n", x0hi,x0lo);
      dprintf("x-b[i]:  nearest(%1.50e) \n",x-b[i].d );
      
      DIV2( xmBIhi , xmBIlo , x0hi,x0lo, xhi,xlo);
      
      dprintf("reduction : X := %1.50f+ ( %1.50f )\n",xhi,xlo);
      dprintf("reduction : b[i] = %1.50e\n",b[i].d);
    
      /* Polynomial evaluation : 
       *  
       *  1rt compute Q(x^2) = (1 - x^2/3 + ...)
       *      then P(x) = x * Q(x^2)
       *
       */

      x2 = xhi*xhi;
      q = coef_poly[0];
      int j=1;
      for (;j<DEGREE;j++) {
        q *= x2;
        q += coef_poly[j];
      }
      q*=x2;
      
      double atanXhi,atanXlo;
      
      Add12(atanXhi,atanXlo,xhi,( (xhi*q)+xlo));
      
      /* reconstruction : atan(x) = atan(b[i]) + atan(x) */
      Add22 (&atanhi,&atanlo, atan_b[i][0].d,atan_b[i][1].d,atanXhi,atanXlo);  
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
      q = coef_poly[0];
      int j=1;
      for (;j<DEGREE;j++) 
        {q *= x2;
        q += coef_poly[j];}
      q*=x2;
      
      Add12(atanhi,atanlo, x , x*q);
    }
  
  /* test if rounding is possible */
  double e = 1.00095;
  if (atanhi == (atanhi + (atanlo*e)))
    return atanhi;
  else
    {/* more accuracy is needed , lauch accurate phase */ 
      return scs_atan_rn(x);
    }
}
