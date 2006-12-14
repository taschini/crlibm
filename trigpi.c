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

void sinpiacc(double *sinpiacc_resh, double *sinpiacc_resm, double *sinpiacc_resl, double x) {
  double sinpiacc_x_0_pow2h, sinpiacc_x_0_pow2m;

  Mul12(&sinpiacc_x_0_pow2h,&sinpiacc_x_0_pow2m,x,x);  

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
 
  sinpiacc_t_1_0h = sinpiacc_coeff_11h;
  sinpiacc_t_2_0h = sinpiacc_t_1_0h * sinpiacc_x_0_pow2h;
  sinpiacc_t_3_0h = sinpiacc_coeff_9h + sinpiacc_t_2_0h;
  sinpiacc_t_4_0h = sinpiacc_t_3_0h * sinpiacc_x_0_pow2h;
  Add12(sinpiacc_t_5_0h,sinpiacc_t_5_0m,sinpiacc_coeff_7h,sinpiacc_t_4_0h);
  MulAdd22(&sinpiacc_t_6_0h,&sinpiacc_t_6_0m,sinpiacc_coeff_5h,sinpiacc_coeff_5m,sinpiacc_x_0_pow2h,sinpiacc_x_0_pow2m,sinpiacc_t_5_0h,sinpiacc_t_5_0m);
  MulAdd22(&sinpiacc_t_7_0h,&sinpiacc_t_7_0m,sinpiacc_coeff_3h,sinpiacc_coeff_3m,sinpiacc_x_0_pow2h,sinpiacc_x_0_pow2m,sinpiacc_t_6_0h,sinpiacc_t_6_0m);
  Mul22(&sinpiacc_t_8_0h,&sinpiacc_t_8_0m,sinpiacc_t_7_0h,sinpiacc_t_7_0m,sinpiacc_x_0_pow2h,sinpiacc_x_0_pow2m);
  Add233Cond(&sinpiacc_t_9_0h,&sinpiacc_t_9_0m,&sinpiacc_t_9_0l,sinpiacc_t_8_0h,sinpiacc_t_8_0m,sinpiacc_coeff_1h,sinpiacc_coeff_1m,sinpiacc_coeff_1l);
  Mul133(&sinpiacc_t_10_0h,&sinpiacc_t_10_0m,&sinpiacc_t_10_0l,x,sinpiacc_t_9_0h,sinpiacc_t_9_0m,sinpiacc_t_9_0l);
  Renormalize3(sinpiacc_resh,sinpiacc_resm,sinpiacc_resl,sinpiacc_t_10_0h,sinpiacc_t_10_0m,sinpiacc_t_10_0l);
}


void cospiacc(double *cospiacc_resh, double *cospiacc_resm, double *cospiacc_resl, double x) {
  double cospiacc_x_0_pow2h, cospiacc_x_0_pow2m;

  Mul12(&cospiacc_x_0_pow2h,&cospiacc_x_0_pow2m,x,x);
  
  double cospiacc_t_1_0h;
  double cospiacc_t_2_0h;
  double cospiacc_t_3_0h;
  double cospiacc_t_4_0h;
  double cospiacc_t_5_0h, cospiacc_t_5_0m;
  double cospiacc_t_6_0h, cospiacc_t_6_0m;
  double cospiacc_t_7_0h, cospiacc_t_7_0m;
  double cospiacc_t_8_0h, cospiacc_t_8_0m;
  double cospiacc_t_9_0h, cospiacc_t_9_0m, cospiacc_t_9_0l;
  
  cospiacc_t_1_0h = cospiacc_coeff_10h;
  cospiacc_t_2_0h = cospiacc_t_1_0h * cospiacc_x_0_pow2h;
  cospiacc_t_3_0h = cospiacc_coeff_8h + cospiacc_t_2_0h;
  cospiacc_t_4_0h = cospiacc_t_3_0h * cospiacc_x_0_pow2h;
  Add12(cospiacc_t_5_0h,cospiacc_t_5_0m,cospiacc_coeff_6h,cospiacc_t_4_0h);
  MulAdd22(&cospiacc_t_6_0h,&cospiacc_t_6_0m,cospiacc_coeff_4h,cospiacc_coeff_4m,cospiacc_x_0_pow2h,cospiacc_x_0_pow2m,cospiacc_t_5_0h,cospiacc_t_5_0m);
  MulAdd22(&cospiacc_t_7_0h,&cospiacc_t_7_0m,cospiacc_coeff_2h,cospiacc_coeff_2m,cospiacc_x_0_pow2h,cospiacc_x_0_pow2m,cospiacc_t_6_0h,cospiacc_t_6_0m);
  Mul22(&cospiacc_t_8_0h,&cospiacc_t_8_0m,cospiacc_t_7_0h,cospiacc_t_7_0m,cospiacc_x_0_pow2h,cospiacc_x_0_pow2m);
  Add123(&cospiacc_t_9_0h,&cospiacc_t_9_0m,&cospiacc_t_9_0l,cospiacc_coeff_0h,cospiacc_t_8_0h,cospiacc_t_8_0m);
  *cospiacc_resh = cospiacc_t_9_0h; *cospiacc_resm = cospiacc_t_9_0m; *cospiacc_resl = cospiacc_t_9_0l;
}





 /* to nearest  */

 double sinpi_rn(double x){
   double xs,xsign, y,u,syh, sym, syl, cyh, cym, cyl, sah, sam, sal, cah, cam, cal;
   double t1h, t1m, t1l, t2h, t2m, t2l, rh, rm, rl;
   db_number t;
   int index, quadrant;


   xs = x*128.0;

   if(x>TWOTO52) /* x is an integer */
     return 0;

   /* Réduction d'argument */
   if(x<TWOTO42){
     t.d = TWOTO5251 + xs;
     u = t.d - TWOTO5251;
     y = xs - u;
     y = y * INV128;
   }

   if(u<0) 
     xsign = -1;
   else 
     xsign = 1;

   index = t.i[LO] & 0x3f;
   quadrant = (t.i[LO] & 0xff) >>6;

   //   printf("\n\nint part = %f    frac part = %f     index=%d   quadrant=%d   \n", u, y, index, quadrant);
   
   sinpiacc(&syh, &sym, &syl, y);
   cospiacc(&cyh, &cym, &cyl, y);
   
   sah=sincosTable[index].sh;
   cah=sincosTable[index].ch;
   sam=sincosTable[index].sm;
   cam=sincosTable[index].cm;
   sal=sincosTable[index].sl;
   cal=sincosTable[index].cl;


   if(quadrant==0 || quadrant==2) {
     /* compute sy*ca+sa*cy */
     
     /* t1 = sy*ca */
     Mul33(&t1h,&t1m,&t1l, syh,sym,syl, cah,cam,cal);
     
     /* t2 =  sa*cy */
     Mul33(&t2h,&t2m,&t2l, sah,sam,sal, cyh,cym,cyl);
     
     /* either index=0, then sa=0 and ca=1, therefore t2=0, and the Add33 will be exact 
	or index !=0, and       
	What we know is: 
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
     
     Add33Cond(&rh, &rm, &rl, t2h,t2m,t2l, t1h,t1m,t1l);
     if (quadrant==2) {
       rh = -rh;
       rm = -rm;
       rl = -rl;
     }
     ReturnRoundToNearest3(rh,rm,rl);
   }


   if(quadrant==1 || quadrant==3) {
     /* compute cy*ca - sa*sy */
     
     /* t1 = cy*ca */
     Mul33(&t1h,&t1m,&t1l, cyh,cym,cyl, cah,cam,cal);
     
     /* t2 =  sa*sy */
     Mul33(&t2h,&t2m,&t2l, sah,sam,sal, syh,sym,syl);     
     
     Add33Cond(&rh, &rm, &rl, t1h,t1m,t1l, -t2h,-t2m,-t2l);

     if (quadrant==3) {
       rh = -rh;
       rm = -rm;
       rl = -rl;
     }
     ReturnRoundToNearest3(rh,rm,rl);   
   }

};


 double sinpi_rd(double x){
   printf("ERROR:  function not yet implemented \n");
   return 0.0/0.0;
}; /* toward -inf */ 
 double sinpi_ru(double x){
   printf("ERROR:  function not yet implemented \n");
   return 0.0/0.0;
}; /* toward +inf */ 
 double sinpi_rz(double x){
   printf("ERROR:  function not yet implemented \n");
   return 0.0/0.0;
}; /* toward zero */ 





/*  cosine of pi times x  */
 double cospi_rn(double x){
   printf("ERROR:  function not yet implemented \n");
   return 0.0/0.0;
};
 /* to nearest  */
 double cospi_rd(double x){
   printf("ERROR:  function not yet implemented \n");
   return 0.0/0.0;
}; /* toward -inf */ 
 double cospi_ru(double x){
   printf("ERROR:  function not yet implemented \n");
   return 0.0/0.0;
}; /* toward +inf */ 
 double cospi_rz(double x){
   printf("ERROR:  function not yet implemented \n");
   return 0.0/0.0;
}; /* toward zero */ 


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


