 /*
 * Function to compute the logarithm with fully exact rounding
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

#define dprintf(format, args...)  printf(format , ## args)
//#define dprintf(format, args...) 

/* the second step :   A VOIR PLUS TARD */ 
double scs_atan_rn(double); 


/*
 * 1) Argument reduction : 
 * 
 *  tan(x) = tan( b(i) ) + tan ( (x-b(i)) / (1-x*b(i)))
 *
 *                                                     6.3
 * we choose 63 b(i) so that (x-b(i)) / (1-x*b(i)) < 2^
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
  db_number x_db;
  x_db.d = x;
  unsigned int hx = x_db.i[HI_ENDIAN] & 0x7FFFFFFF; 
  
  
  /* first : test special cases : */
  if ( hx >= 0x43500000)
    { /* x >= 2^54 */
      #ifdef DEBUG 
      dprintf("grand !");
      #endif
      if ( ( (hx & 0x000fffff) | x_db.i[LO_ENDIAN] ) == 0)
        return x+x;                /* NaN */
      else
        return HALFPI.d;           /* atan(x) = Pi/2 */
    }
  else
    if ( hx < 0x3CA00000 )
      {dprintf("coucou");
      return x; /* x<2^-53 then atan(x) =~ x */}
  
  /* TODO : test if x < e !! */
  
  /* determine i so that x E [a[i],a[i+1]] */

  int i=31; int j;
  dprintf("%f < %f ? %d\n",x,a[i],i);
  if (x < a[i]) i-= 16;
  else i+=16;
  dprintf("%f < %f ? %d\n",x,a[i],i);
  if (x < a[i]) i-= 8;
  else i+= 8;
  dprintf("%f < %f ? %d\n",x,a[i],i);
  if (x < a[i]) i-= 4;
  else i+= 4;
  dprintf("%f < %f ? %d\n",x,a[i],i);
  if (x < a[i]) i-= 2;
  else i+= 2;
  dprintf("%f < %f ? %d\n",x,a[i],i);
  if (x < a[i]) i-= 1;
  else if (i<61) i+= 1;
  dprintf("%f < %f ? %d\n",x,a[i],i);
  if (x < a[i]) i-= 1;
  
  
  /* for (;j<log_of_nb_of_ai && i < nb_of_ai;j++) */
  /*     {    */
  /*       i+=k; */
  /*       if (x < a[i]) */
  /*         {dprintf("%f < %f donc on elveve %d\n",x,a[i],k); */
  /*         i-=k;} */
  /*       /\* a ameliorer  */
  /*          normalement ca marche car au depart k = 2^.*\/ */
  /*     dprintf("om teste si %f < %f et si oui on elveve %d (i=%d)\n",x,a[i],k,i); */
  /*     k = (k+1) >> 1; */
  /*     } */
  /*   i--; */
  /*   if (i>=nb_of_bi) */
  /*     i = nb_of_bi-1; */
  
  dprintf("i = %d (b[i] = %f,%f<%f<%f = )\n",
         i,b[i][0].d,a[i], x ,a[i+1]);
  
  /* we now compute X = ( x-b[i] ) / ( 1 - x*b[i] )
   * 
   * def : x0 := 1-x*b[i]
   *
   * 1st we compute an approximation of y = 1/x0
   * then we compute a better approx x' = y*(2-x0*y)
   * we can proove that :
   * if y = 1/x0*(1+e) 
   *     then x' = 1/x0 * (1-e^2)
   *                   
   * TODO : test if x as to be accurate !!
  */
  
  double x0 = ( 1 + x*b[i][0].d );
  double y = 1/x0;
  double xhi,xlo, temphi,templo, multhi,multlo;
  
  /*if (i<=1)
    {*/
      if (i==-1)
        {xhi=x;
        xlo=0.0;
        }
      else {
        double tmphi, tmplo;
        Mul12(&multhi,&multlo, y,x0);
        // a amelio
        
        
        tmphi = (1-multhi)-multlo;
        Mul12 (&temphi,&templo, 1+tmphi,y);
        
        /* Add22Cond(&tmphi,&tmplo, 2, 0, -multhi, - multlo); */
/*         Mul22( &temphi, & templo, y,0.0,tmphi, tmplo); */
        
        
        Mul22( & xhi , & xlo , (x-b[i][0].d), 0.0 , temphi, templo);
        dprintf("coucou (encore !!)\n");
      }
      /*}
        else
        {xhi=( x-b[i][0].d ) / ( 1 + x*b[i][0].d );
        xlo=0;}
      */
  dprintf("reduction : X= %f\n",xhi);
  
  /* Ponynomial evaluation : 
   *  
   *  1rt compute Q(x^2) = (1 - x^2/3 + ...)
   *      then P(x) = x * Q(x^2)
   *
   */
  double q,x2;
  x2 = xhi*xhi;
  q = coef_poly[0];
  for (j=1;j<DEGREE;j++) {
    q *= x2;
    q += coef_poly[j];
    dprintf("value of q: %f ( a la %diem iteration)\n",q,j);
  }
  q*=x2;
  
  double qhi,qlo;
  double atanXhi,atanXlo;
  //if (i<=1) {
  
  //A AMELIO

  Add12 (qhi,qlo, 1.0, q);
  Mul22(&atanXhi, &atanXlo, xhi,xlo, qhi,qlo);
  
  dprintf("atanXhi=%f,   atanXlo=%f\n",atanXhi,atanXlo);
  // **  }
  //else
  {atanXhi = xhi*(1+q); atanXlo = 0;}

  dprintf("and now : atan(X=%f) = %f\n",xhi,atanXhi);
  
  /* reconstruction : atan(x) = atan(b[i]) + atan(x) */
  if (i>=0)
    {Add22 (&atanhi,&atanlo, atan_b[i][0].d,atan_b[i][1].d,atanXhi,atanXlo);}
  else {atanhi=atanXhi;atanlo=atanXlo;}
  
  /* test if rounding is possible */
  /* (TODO) */ 
  
  int e = 2;
  
  if (atanhi == (atanhi + (atanlo*e)))
    {  dprintf ("cr_libm     : ");return atanhi;}
  else
    {/* more accuracy is needed , lauch accurate phase */ 
      dprintf("2em etape avec %f\n",x);
      dprintf ("cr_libm     : ");
      return scs_atan_rn(x);
    }
}
