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

/* the second step :   A VOIR PLUS TARD */ 
void scs_atan(scs_ptr,db_number, int); 


/*
 * 1) Argument reduction : 
 * 
 *  tan(x) = tan( b(i) ) + tan ( (x-b(i)) / (1-x*b(i)))
 *
 *                                                     6.3
 * we choose 63 b(i) so that (x-b(i)) / (1-x*b(i)) < 2^
 */ 


extern double atan_rn (double x) {
  
  double atanhi,atanlo;
  db_number x_db;
  x_db.d = x;
  unsigned int hx = x_db.i[HI_ENDIAN] & 0x7FFFFFFF; 
  
  //return 0.0;
  
  /* first : test special cases : */
  if ( hx >= 0x43500000)
    { /* x >= 2^54 */
      printf("grand !");
      if ( ( (hx & 0x000fffff) | x_db.i[LO_ENDIAN] ) == 0)
        return x+x;                /* NaN */
      else
        return HALFPI.d;           /* atan(x) = Pi/2 */
    }
  else
    if ( hx < 0x3CA00000 )
      {printf("coucou");
      return x; /* x<2^-53 then atan(x) =~ x */}
  
  /* TODO : test if x < e !! */
  
  /* determine i so that x E [a[i],a[i+1]] */

  int i=31; int j;
  printf("%f < %f ? %d\n",x,a[i],i);
  if (x < a[i]) i-= 16;
  else i+=16;
  printf("%f < %f ? %d\n",x,a[i],i);
  if (x < a[i]) i-= 8;
  else i+= 8;
  printf("%f < %f ? %d\n",x,a[i],i);
  if (x < a[i]) i-= 4;
  else i+= 4;
  printf("%f < %f ? %d\n",x,a[i],i);
  if (x < a[i]) i-= 2;
  else i+= 2;
  printf("%f < %f ? %d\n",x,a[i],i);
  if (x < a[i]) i-= 1;
  else if (i<61) i+= 1;
  printf("%f < %f ? %d\n",x,a[i],i);
  if (x < a[i]) i-= 1;
  
  
  /* for (;j<log_of_nb_of_ai && i < nb_of_ai;j++) */
  /*     {    */
  /*       i+=k; */
  /*       if (x < a[i]) */
  /*         {printf("%f < %f donc on elveve %d\n",x,a[i],k); */
  /*         i-=k;} */
  /*       /\* a ameliorer  */
  /*          normalement ca marche car au depart k = 2^.*\/ */
  /*     printf("om teste si %f < %f et si oui on elveve %d (i=%d)\n",x,a[i],k,i); */
  /*     k = (k+1) >> 1; */
  /*     } */
  /*   i--; */
  /*   if (i>=nb_of_bi) */
  /*     i = nb_of_bi-1; */
  
  printf("i = %d (b[i] = %f,%f<%f<%f = )\n",
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
        double tmphi,tmplo;
        Mul12(&multhi,&multlo, y,x0);
        // a amelio
        Add22Cond(&temphi,&templo, 2, 0, -multhi, - multlo);
        Mul22( &tmphi, & tmplo, y,0.0,temphi, templo);
        //
        Mul22( & xhi , & xlo , (x-b[i][0].d), 0.0 , tmphi, tmplo);
        printf("coucou (encore !!)\n");
      }
      /*}
        else
        {xhi=( x-b[i][0].d ) / ( 1 + x*b[i][0].d );
        xlo=0;}
      */
  printf("reduction : X= %f\n",xhi);
  
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
    printf("value of q: %f ( a la %diem iteration)\n",q,j);
  }
  q*=x2;
  
  double qhi,qlo;
  double atanXhi,atanXlo;
  //if (i<=1) {

  //A AMELIO

  Add12 (qhi,qlo, 1.0, q);
  Mul22(&atanXhi, &atanXlo, xhi,xlo, qhi,qlo);
  printf("atanXhi=%f,   atanXlo=%f\n",atanXhi,atanXlo);
  // **  }
  //else
  {atanXhi = xhi*(1+q); atanXlo = 0;}

  printf("and now : atan(X=%f) = %f\n",xhi,atanXhi);
  
  /* reconstruction : atan(x) = atan(b[i]) + atan(x) */
  if (i>=0)
    {Add22 (&atanhi,&atanlo, atan_b[i][0].d,atan_b[i][1].d,atanXhi,atanXlo);}
  else {atanhi=atanXhi;atanlo=atanXlo;}
  
  /* test if rounding is possible */
  /* (TODO) */ 

  printf ("cr_libm     : ");

  return atanlo+atanhi; 
}
