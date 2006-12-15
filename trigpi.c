#include <stdio.h>
#include <stdlib.h>
#include "crlibm.h"
#include "crlibm_private.h"
#include "triple-double.h"
#include "trigpi.h"

#define sinpiacc_coeff_1h 3.14159265358979311599796346854418516159057617187500000000000000000000000000000000e+00
#define sinpiacc_coeff_1m 1.22464679914735320717376402945839660462569212467758006379625612680683843791484833e-16
#define sinpiacc_coeff_1l -2.87889731599645993207191707893463395148177292198731390393739579574603302514349608e-33
#define sinpiacc_coeff_3h -5.16771278004997025590228076907806098461151123046875000000000000000000000000000000e+00
#define sinpiacc_coeff_3m 2.26656228257550136196266687046492287115561324595258696490418515168130397796630859e-16
#define sinpiacc_coeff_5h 2.55016403987734552316624103696085512638092041015625000000000000000000000000000000e+00
#define sinpiacc_coeff_5m -7.93098961936403945684716222915171282926664203267314023904077657789457589387893677e-17
#define sinpiacc_coeff_7h -5.99264529320792105338000510528218001127243041992187500000000000000000000000000000e-01
#define sinpiacc_coeff_9h 8.21458866130424236740026344705256633460521697998046875000000000000000000000000000e-02
#define sinpiacc_coeff_11h -7.37046804820839888960914976223648409359157085418701171875000000000000000000000000e-03


#define cospiacc_coeff_0h 1.00000000000000000000000000000000000000000000000000000000000000000000000000000000e+00
#define cospiacc_coeff_2h -4.93480220054467899615247006295248866081237792968750000000000000000000000000000000e+00
#define cospiacc_coeff_2m -3.13264775437072047133490817894057799839785556899468543790021612949203699827194214e-16
#define cospiacc_coeff_4h 4.05871212641676848420502210501581430435180664062500000000000000000000000000000000e+00
#define cospiacc_coeff_4m -2.66019969731660223662555032718185048048635055542576743903282476821914315223693848e-16
#define cospiacc_coeff_6h -1.33526276885458949905682857206556946039199829101562500000000000000000000000000000e+00
#define cospiacc_coeff_8h 2.35330630358513925859398341344785876572132110595703125000000000000000000000000000e-01
#define cospiacc_coeff_10h -2.58068327360992909313974763563237502239644527435302734375000000000000000000000000e-02

static void sincospiacc(double *sinpiacc_resh, double *sinpiacc_resm, double *sinpiacc_resl,
		 double *cospiacc_resh, double *cospiacc_resm, double *cospiacc_resl,
		 double x) {
  double x2h, x2m;
  double sinpiacc_t_1_0h;
  double sinpiacc_t_2_0h;
  double sinpiacc_t_3_0h;
  double sinpiacc_t_4_0h;
  double sinpiacc_t_5_0h, sinpiacc_t_5_0m;
  double sinpiacc_t_6_0h, sinpiacc_t_6_0m;
  double sinpiacc_t_7_0h, sinpiacc_t_7_0m;
  double sinpiacc_t_8_0h, sinpiacc_t_8_0m;
  double sinpiacc_t_9_0h, sinpiacc_t_9_0m, sinpiacc_t_9_0l;
  double sinpiacc_t_10_0h, sinpiacc_t_10_0m, sinpiacc_t_10_0l;
 
  double cospiacc_t_1_0h;
  double cospiacc_t_2_0h;
  double cospiacc_t_3_0h;
  double cospiacc_t_4_0h;
  double cospiacc_t_5_0h, cospiacc_t_5_0m;
  double cospiacc_t_6_0h, cospiacc_t_6_0m;
  double cospiacc_t_7_0h, cospiacc_t_7_0m;
  double cospiacc_t_8_0h, cospiacc_t_8_0m;
  double cospiacc_t_9_0h, cospiacc_t_9_0m, cospiacc_t_9_0l;

  Mul12(&x2h,&x2m,x,x);  

  sinpiacc_t_1_0h = sinpiacc_coeff_11h;
  sinpiacc_t_2_0h = sinpiacc_t_1_0h * x2h;
  sinpiacc_t_3_0h = sinpiacc_coeff_9h + sinpiacc_t_2_0h;
  sinpiacc_t_4_0h = sinpiacc_t_3_0h * x2h;
  Add12(sinpiacc_t_5_0h,sinpiacc_t_5_0m,sinpiacc_coeff_7h,sinpiacc_t_4_0h);
  MulAdd22(&sinpiacc_t_6_0h,&sinpiacc_t_6_0m,sinpiacc_coeff_5h,sinpiacc_coeff_5m,x2h,x2m,sinpiacc_t_5_0h,sinpiacc_t_5_0m);
  MulAdd22(&sinpiacc_t_7_0h,&sinpiacc_t_7_0m,sinpiacc_coeff_3h,sinpiacc_coeff_3m,x2h,x2m,sinpiacc_t_6_0h,sinpiacc_t_6_0m);
  Mul22(&sinpiacc_t_8_0h,&sinpiacc_t_8_0m,sinpiacc_t_7_0h,sinpiacc_t_7_0m,x2h,x2m);
  Add233Cond(&sinpiacc_t_9_0h,&sinpiacc_t_9_0m,&sinpiacc_t_9_0l,sinpiacc_t_8_0h,sinpiacc_t_8_0m,sinpiacc_coeff_1h,sinpiacc_coeff_1m,sinpiacc_coeff_1l);
  Mul133(&sinpiacc_t_10_0h,&sinpiacc_t_10_0m,&sinpiacc_t_10_0l,x,sinpiacc_t_9_0h,sinpiacc_t_9_0m,sinpiacc_t_9_0l);
  Renormalize3(sinpiacc_resh,sinpiacc_resm,sinpiacc_resl,sinpiacc_t_10_0h,sinpiacc_t_10_0m,sinpiacc_t_10_0l);

  
  cospiacc_t_1_0h = cospiacc_coeff_10h;
  cospiacc_t_2_0h = cospiacc_t_1_0h * x2h;
  cospiacc_t_3_0h = cospiacc_coeff_8h + cospiacc_t_2_0h;
  cospiacc_t_4_0h = cospiacc_t_3_0h * x2h;
  Add12(cospiacc_t_5_0h,cospiacc_t_5_0m,cospiacc_coeff_6h,cospiacc_t_4_0h);
  MulAdd22(&cospiacc_t_6_0h,&cospiacc_t_6_0m,cospiacc_coeff_4h,cospiacc_coeff_4m,x2h,x2m,cospiacc_t_5_0h,cospiacc_t_5_0m);
  MulAdd22(&cospiacc_t_7_0h,&cospiacc_t_7_0m,cospiacc_coeff_2h,cospiacc_coeff_2m,x2h,x2m,cospiacc_t_6_0h,cospiacc_t_6_0m);
  Mul22(&cospiacc_t_8_0h,&cospiacc_t_8_0m,cospiacc_t_7_0h,cospiacc_t_7_0m,x2h,x2m);
  Add123(&cospiacc_t_9_0h,&cospiacc_t_9_0m,&cospiacc_t_9_0l,cospiacc_coeff_0h,cospiacc_t_8_0h,cospiacc_t_8_0m);
  *cospiacc_resh = cospiacc_t_9_0h; *cospiacc_resm = cospiacc_t_9_0m; *cospiacc_resl = cospiacc_t_9_0l;

}








/* Comment on comparing sa, ca, sy and cy 
   either index=0, then sa=0 and ca=1, therefore t2=0, and the Add33 will be exact 
   or index !=0, and       
	-eps1 < sy < eps1  (but sy may be negative)
	sa > eps1 (sa>0)
	1-eps2 < cy < 1
	ca < 1-eps2
	therefore 
	sacy = t2 >0
	casy = t1 may be negative
	abs(t1) <=  abs(t2)
	Unfortunately we need a stronger condition to use the Add33 
*/




static void sinpi_accurate(double *rh, double *rm, double *rl,
			   double y, int index, int quadrant)
{
   double syh, sym, syl, cyh, cym, cyl, sah, sam, sal, cah, cam, cal;
   double t1h, t1m, t1l, t2h, t2m, t2l;

   sincospiacc(&syh, &sym, &syl, &cyh, &cym, &cyl, y);
   
   sah=sincosTable[index].sh;
   cah=sincosTable[index].ch;
   sam=sincosTable[index].sm;
   cam=sincosTable[index].cm;
   sal=sincosTable[index].sl;
   cal=sincosTable[index].cl;

   if(quadrant==0 || quadrant==2) {
     /* compute sy*ca+sa*cy   :    t1 = sy*ca,     t2 =  sa*cy*/
     Mul33(&t1h,&t1m,&t1l, syh,sym,syl, cah,cam,cal);
     Mul33(&t2h,&t2m,&t2l, sah,sam,sal, cyh,cym,cyl);
     Add33Cond(rh, rm, rl, t2h,t2m,t2l, t1h,t1m,t1l);
   }
   else {
     /* compute cy*ca - sa*sy : t1 = cy*ca,    t2 =  sa*sy */
     Mul33(&t1h,&t1m,&t1l, cyh,cym,cyl, cah,cam,cal);
     Mul33(&t2h,&t2m,&t2l, sah,sam,sal, syh,sym,syl);     
     Add33Cond(rh, rm, rl, t1h,t1m,t1l, -t2h,-t2m,-t2l);
   }

   if (quadrant>=2) {
     *rh = -*rh;
     *rm = -*rm;
     *rl = -*rl;
   }

};


static void cospi_accurate(double *rh, double *rm, double *rl,
			   double y, int index, int quadrant)
{
   double syh, sym, syl, cyh, cym, cyl, sah, sam, sal, cah, cam, cal;
   double t1h, t1m, t1l, t2h, t2m, t2l;

   sincospiacc(&syh, &sym, &syl, &cyh, &cym, &cyl, y);
   
   sah=sincosTable[index].sh;
   cah=sincosTable[index].ch;
   sam=sincosTable[index].sm;
   cam=sincosTable[index].cm;
   sal=sincosTable[index].sl;
   cal=sincosTable[index].cl;

   if(quadrant==0 || quadrant==2) {
     /* compute cy*ca - sa*sy : t1 = cy*ca,    t2 =  sa*sy */
     Mul33(&t1h,&t1m,&t1l, cyh,cym,cyl, cah,cam,cal);
     Mul33(&t2h,&t2m,&t2l, sah,sam,sal, syh,sym,syl);     
     Add33Cond(rh, rm, rl, t1h,t1m,t1l, -t2h,-t2m,-t2l);
   }
   else {
     /* compute sy*ca+sa*cy   :    t1 = sy*ca,     t2 =  sa*cy*/
     Mul33(&t1h,&t1m,&t1l, syh,sym,syl, cah,cam,cal);
     Mul33(&t2h,&t2m,&t2l, sah,sam,sal, cyh,cym,cyl);
     Add33Cond(rh, rm, rl, t2h,t2m,t2l, t1h,t1m,t1l);
   }

   if (quadrant==1 || quadrant==2) {
     *rh = -*rh;
     *rm = -*rm;
     *rl = -*rl;
   }

};





 double sinpi_rn(double x){
   double xs, y,u, rh, rm, rl, sign,absx;
   db_number xdb, t;
   int32_t xh, absxh, index, quadrant;
   
   if (x<0) absx = -x;   else absx = x; 

   xdb.d = x;

   xs = x*128.0;

   /* argument reduction */
   if(absx>  TWOTO42 ) {  /* x is very large, let us first subtract a large integer from it */
     t.d = xs;
     t.i[LO] =0; /* remove the low part. The point is somewhere there since x > 2^42. 
		    So what remains in t is an FP integer almost as large as x */
     xs = xs-t.d; /* we are going to throw away the int part anyway */ 
   }

   t.d = TWOTO5251 + xs;
   u = t.d - TWOTO5251;
   y = xs - u;
   index = t.i[LO] & 0x3f;
   quadrant = (t.i[LO] & 0xff) >>6;

   /* Special case tests come late because the conversion FP to int is slow */
   xh = xdb.i[HI];
   absxh = xh & 0x7fffffff;
   if (xh<0)  sign=-1.;   else sign=1.; /* consider the sign bit */

   if(index==0 && y==0.0 && ((quadrant&1)==0)) return sign*0.0; /*signed, inspired by LIA-2 */

   y = y * INV128;

   /* SPECIAL CASES: x=(Nan, Inf) sin(pi*x)=Nan */
   if (absxh>=0x7ff00000) {
     xdb.l=0xfff8000000000000LL;
     return xdb.d - xdb.d; 
   }
      
   if(absxh>=0x43300000) /* 2^52, which entails that x is an integer */
     return sign*0.0; /*signed */

   sinpi_accurate(&rh, &rm, &rl, y, index, quadrant);
   ReturnRoundToNearest3(rh,rm,rl);   

 }







 double sinpi_rd(double x){
   double xs, y,u, rh, rm, rl, sign,absx;
   db_number xdb, t;
   int32_t xh, absxh, index, quadrant;
   
   if (x<0) absx = -x;   else absx = x; 

   xdb.d = x;

   xs = x*128.0;

   /* argument reduction */
   if(absx>  TWOTO42 ) {  /* x is very large, let us first subtract a large integer from it */
     t.d = xs;
     t.i[LO] =0; /* remove the low part. The point is somewhere there since x > 2^42. 
		    So what remains in t is an FP integer almost as large as x */
     xs = xs-t.d; /* we are going to throw away the int part anyway */ 
   }

   t.d = TWOTO5251 + xs;
   u = t.d - TWOTO5251;
   y = xs - u;
   index = t.i[LO] & 0x3f;
   quadrant = (t.i[LO] & 0xff) >>6;

   /* Special case tests come late because the conversion FP to int is slow */
   xh = xdb.i[HI];
   absxh = xh & 0x7fffffff;
   if (xh<0)  sign=-1.;   else sign=1.; /* consider the sign bit */

   if(index==0 && y==0.0 && ((quadrant&1)==0)) return sign*0.0; /*signed, inspired by LIA-2 */

   y = y * INV128;

   /* SPECIAL CASES: x=(Nan, Inf) sin(pi*x)=Nan */
   if (absxh>=0x7ff00000) {
     xdb.l=0xfff8000000000000LL;
     return xdb.d - xdb.d; 
   }
      
   if(absxh>=0x43300000) /* 2^52, which entails that x is an integer */
     return sign*0.0; /*signed */

   sinpi_accurate(&rh, &rm, &rl, y, index, quadrant);
   ReturnRoundDownwards3(rh,rm,rl);   
}; 



 double sinpi_ru(double x){
   double xs, y,u, rh, rm, rl, sign,absx;
   db_number xdb, t;
   int32_t xh, absxh, index, quadrant;
   
   if (x<0) absx = -x;   else absx = x; 

   xdb.d = x;

   xs = x*128.0;

   /* argument reduction */
   if(absx>  TWOTO42 ) {  /* x is very large, let us first subtract a large integer from it */
     t.d = xs;
     t.i[LO] =0; /* remove the low part. The point is somewhere there since x > 2^42. 
		    So what remains in t is an FP integer almost as large as x */
     xs = xs-t.d; /* we are going to throw away the int part anyway */ 
   }

   t.d = TWOTO5251 + xs;
   u = t.d - TWOTO5251;
   y = xs - u;
   index = t.i[LO] & 0x3f;
   quadrant = (t.i[LO] & 0xff) >>6;

   /* Special case tests come late because the conversion FP to int is slow */
   xh = xdb.i[HI];
   absxh = xh & 0x7fffffff;
   if (xh<0)  sign=-1.;   else sign=1.; /* consider the sign bit */

   if(index==0 && y==0.0 && ((quadrant&1)==0)) return sign*0.0; /*signed, inspired by LIA-2 */

   y = y * INV128;

   /* SPECIAL CASES: x=(Nan, Inf) sin(pi*x)=Nan */
   if (absxh>=0x7ff00000) {
     xdb.l=0xfff8000000000000LL;
     return xdb.d - xdb.d; 
   }
      
   if(absxh>=0x43300000) /* 2^52, which entails that x is an integer */
     return sign*0.0; /*signed */

   sinpi_accurate(&rh, &rm, &rl, y, index, quadrant);
   ReturnRoundUpwards3(rh,rm,rl);   
   
};  




 double sinpi_rz(double x){
   double xs, y,u, rh, rm, rl, sign,absx;
   db_number xdb, t;
   int32_t xh, absxh, index, quadrant;
   
   if (x<0) absx = -x;   else absx = x; 

   xdb.d = x;

   xs = x*128.0;

   /* argument reduction */
   if(absx>  TWOTO42 ) {  /* x is very large, let us first subtract a large integer from it */
     t.d = xs;
     t.i[LO] =0; /* remove the low part. The point is somewhere there since x > 2^42. 
		    So what remains in t is an FP integer almost as large as x */
     xs = xs-t.d; /* we are going to throw away the int part anyway */ 
   }

   t.d = TWOTO5251 + xs;
   u = t.d - TWOTO5251;
   y = xs - u;
   index = t.i[LO] & 0x3f;
   quadrant = (t.i[LO] & 0xff) >>6;

   /* Special case tests come late because the conversion FP to int is slow */
   xh = xdb.i[HI];
   absxh = xh & 0x7fffffff;
   if (xh<0)  sign=-1.;   else sign=1.; /* consider the sign bit */

   if(index==0 && y==0.0 && ((quadrant&1)==0)) return sign*0.0; /*signed, inspired by LIA-2 */

   y = y * INV128;

   /* SPECIAL CASES: x=(Nan, Inf) sin(pi*x)=Nan */
   if (absxh>=0x7ff00000) {
     xdb.l=0xfff8000000000000LL;
     return xdb.d - xdb.d; 
   }
      
   if(absxh>=0x43300000) /* 2^52, which entails that x is an integer */
     return sign*0.0; /*signed */

   sinpi_accurate(&rh, &rm, &rl, y, index, quadrant);
   ReturnRoundTowardsZero3(rh,rm,rl);   
};









 /* to nearest  */
 double cospi_rn(double x){
   double xs, y,u, rh, rm, rl, absx;
   db_number xdb, t;
   int32_t xh, absxh, index, quadrant;
   
   if (x<0) absx =-x; else absx=x; 

   xdb.d=x;
   xs = x*128.0;

   /* argument reduction. 
      We do it before the special case tests for performance, 
      it might compute garbage for inf, very large inputs, etc */
   if(absx>  TWOTO42 ) {  /* 2^42, x is very large, let us first subtract a large integer from it */
     t.d = xs;
     t.i[LO] =0; /* remove the low part, in which was the coma since x > 2^42. 
		    So what remains in t is an FP integer almost as large as x */
     xs = xs-t.d; /* we are going to throw away the int part anyway */ 
   }

   t.d = TWOTO5251 + xs;
   u = t.d - TWOTO5251;
   y = xs - u;
   y = y * INV128;
   index = t.i[LO] & 0x3f;
   quadrant = (t.i[LO] & 0xff) >>6;

   if(index==0 && y==0. && ((quadrant&1)==1)) return +0.; 
   /* Always +0, inpired by LIA2; We do not have cos(x+pi) == - cos(x)
      in this case */

   /* Special case tests come late because the conversion FP to int is slow */
   xh = xdb.i[HI];
   absxh = xh & 0x7fffffff;

   /* SPECIAL CASES: x=(Nan, Inf) cos(pi*x)=Nan */

   if (absxh>=0x7ff00000) {
     xdb.l=0xfff8000000000000LL;
     return xdb.d - xdb.d; 
   }
      
   if(absxh>=0x43400000) /* 2^53, which entails that x is an even integer */
     return 1.0; 

   if(index==0 && y==0. && quadrant==0) return 1.; 
   if(index==0 && y==0. && quadrant==2) return -1.; 

   if (xh<0x3E200000) /* 2^-29 */
     return 1;


   //   printf("\n\nint part = %f    frac part = %f     index=%d   quadrant=%d   \n", u, y, index, quadrant);

   cospi_accurate(&rh, &rm, &rl, y, index, quadrant);
   ReturnRoundToNearest3(rh,rm,rl);   
};





 double cospi_rd(double x){
   double xs, y,u, rh, rm, rl, absx;
   db_number xdb, t;
   int32_t xh, absxh, index, quadrant;
   
   if (x<0) absx =-x; else absx=x; 

   xdb.d=x;
   xs = x*128.0;


   /* argument reduction. 
      We do it before the special case tests for performance, 
      it might compute garbage for inf, very large inputs, etc */
   if(absx>  TWOTO42 ) {  /* 2^42, x is very large, let us first subtract a large integer from it */
     t.d = xs;
     t.i[LO] =0; /* remove the low part, in which was the coma since x > 2^42. 
		    So what remains in t is an FP integer almost as large as x */
     xs = xs-t.d; /* we are going to throw away the int part anyway */ 
   }

   t.d = TWOTO5251 + xs;
   u = t.d - TWOTO5251;
   y = xs - u;
   y = y * INV128;
   index = t.i[LO] & 0x3f;
   quadrant = (t.i[LO] & 0xff) >>6;


   /* Special case tests come late because the conversion FP to int is slow */
   xh = xdb.i[HI];
   absxh = xh & 0x7fffffff;

   /* SPECIAL CASES: x=(Nan, Inf) cos(pi*x)=Nan */

   if (absxh>=0x7ff00000) {
     xdb.l=0xfff8000000000000LL;
     return xdb.d - xdb.d; 
   }
      
   if(absxh>=0x43400000) /* 2^53, which entails that x is an even integer */
     return 1.0; /*signed */

   if(index==0 && y==0. && ((quadrant&1)==1)) return +0.; 

   if(index==0 && y==0. && quadrant==0) return 1.; 
   if(index==0 && y==0. && quadrant==2) return -1.; 

   if (xh<0x3E200000) /* 2^-29 */
     return 0.9999999999999998889776975374843459576368331909179687500; /* 1-2^-53 */
   /* Always +0, inpired by LIA2; We do not have cos(x+pi) == - cos(x)
      in this case */

   //   printf("\n\nint part = %f    frac part = %f     index=%d   quadrant=%d   \n", u, y, index, quadrant);

   cospi_accurate(&rh, &rm, &rl, y, index, quadrant);
   ReturnRoundDownwards3(rh,rm,rl);  
 }; 



 
 double cospi_ru(double x){
   double xs, y,u, rh, rm, rl, absx;
   db_number xdb, t;
   int32_t xh, absxh, index, quadrant;
   
   if (x<0) absx =-x; else absx=x; 

   xdb.d=x;
   xs = x*128.0;


   /* argument reduction. 
      We do it before the special case tests for performance, 
      it might compute garbage for inf, very large inputs, etc */
   if(absx>  TWOTO42 ) {  /* 2^42, x is very large, let us first subtract a large integer from it */
     t.d = xs;
     t.i[LO] =0; /* remove the low part, in which was the coma since x > 2^42. 
		    So what remains in t is an FP integer almost as large as x */
     xs = xs-t.d; /* we are going to throw away the int part anyway */ 
   }

   t.d = TWOTO5251 + xs;
   u = t.d - TWOTO5251;
   y = xs - u;
   y = y * INV128;
   index = t.i[LO] & 0x3f;
   quadrant = (t.i[LO] & 0xff) >>6;


   /* Special case tests come late because the conversion FP to int is slow */
   xh = xdb.i[HI];
   absxh = xh & 0x7fffffff;

   /* SPECIAL CASES: x=(Nan, Inf) cos(pi*x)=Nan */
   if (absxh>=0x7ff00000) {
     xdb.l=0xfff8000000000000LL;
     return xdb.d - xdb.d; 
   }
      
   if(absxh>=0x43400000) /* 2^53, which entails that x is an even integer */
     return 1.0; /*signed */

   if(index==0 && y==0. && quadrant==0) return 1.; 
   if(index==0 && y==0. && quadrant==2) return -1.; 

   if(index==0 && y==0. && ((quadrant&1)==1)) return +0.; 
   /* Always +0, inpired by LIA2; We do not have cos(x+pi) == - cos(x)
      in this case */

   if (xh<0x3E200000) /* 2^-29 */
     return 1;


   //   printf("\n\nint part = %f    frac part = %f     index=%d   quadrant=%d   \n", u, y, index, quadrant);

   cospi_accurate(&rh, &rm, &rl, y, index, quadrant);
   ReturnRoundUpwards3(rh,rm,rl);  
}; 




double cospi_rz(double x){
   double xs, y,u, rh, rm, rl, absx;
   db_number xdb, t;
   int32_t xh, absxh, index, quadrant;
   
   if (x<0) absx =-x; else absx=x; 

   xdb.d=x;
   xs = x*128.0;


   /* argument reduction. 
      We do it before the special case tests for performance, 
      it might compute garbage for inf, very large inputs, etc */
   if(absx>  TWOTO42 ) {  /* 2^42, x is very large, let us first subtract a large integer from it */
     t.d = xs;
     t.i[LO] =0; /* remove the low part, in which was the coma since x > 2^42. 
		    So what remains in t is an FP integer almost as large as x */
     xs = xs-t.d; /* we are going to throw away the int part anyway */ 
   }

   t.d = TWOTO5251 + xs;
   u = t.d - TWOTO5251;
   y = xs - u;
   y = y * INV128;
   index = t.i[LO] & 0x3f;
   quadrant = (t.i[LO] & 0xff) >>6;


   /* Special case tests come late because the conversion FP to int is slow */
   xh = xdb.i[HI];
   absxh = xh & 0x7fffffff;

   /* SPECIAL CASES: x=(Nan, Inf) cos(pi*x)=Nan */

   if (absxh>=0x7ff00000) {
     xdb.l=0xfff8000000000000LL;
     return xdb.d - xdb.d; 
   }
      
   if(absxh>=0x43400000) /* 2^53, which entails that x is an even integer */
     return 1.0; /*signed */

   if(index==0 && y==0. && ((quadrant&1)==1)) return +0.; 
   /* Always +0, inpired by LIA2; We do not have cos(x+pi) == - cos(x)
      in this case */

   if(index==0 && y==0. && quadrant==0) return 1.; 
   if(index==0 && y==0. && quadrant==2) return -1.; 

   if (xh<0x3E200000) /* 2^-29 */
     return 0.9999999999999998889776975374843459576368331909179687500; /* 1-2^-53 */

   //   printf("\n\nint part = %f    frac part = %f     index=%d   quadrant=%d   \n", u, y, index, quadrant);

   cospi_accurate(&rh, &rm, &rl, y, index, quadrant);
   ReturnRoundTowardsZero3(rh,rm,rl);
  }; 

/*  tangent of pi times x */
 double tanpi_rn(double x){
   printf("ERROR:  function not yet implemented \n");
   return 0.0/0.0;
}; /* to nearest  */
 double tanpi_rd(double x){
   printf("ERROR:  function not yet implemented \n");
   return 0.0/0.0;
}; /* toward -inf */ 
 double tanpi_ru(double x){
   printf("ERROR:  function not yet implemented \n");
   return 0.0/0.0;
}; /* toward +inf */
 double tanpi_rz(double x){
   printf("ERROR:  function not yet implemented \n");
   return 0.0/0.0;
}; /* toward zero */


