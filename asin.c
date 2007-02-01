/*
 * Correctly rounded arcsine
 *
 * Author : Christoph Lauter (ENS Lyon)
 *
 * This file is part of the crlibm library developed by the Arenaire
 * project at Ecole Normale Superieure de Lyon
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
*/

#include <stdio.h>
#include <stdlib.h>
#include "crlibm.h"
#include "crlibm_private.h"
#include "triple-double.h"
#include "asin.h"

static inline void p0_quick(double *p_resh, double *p_resm, double x, int32_t xhi) {
double p_x_0_pow2h, p_x_0_pow2m;





double p_t_1_0h;
double p_t_2_0h;
double p_t_3_0h;
double p_t_4_0h;
double p_t_5_0h;
double p_t_6_0h;
double p_t_7_0h;
double p_t_8_0h;
double p_t_9_0h;
double p_t_10_0h;
double p_t_11_0h;
double p_t_12_0h;
double p_t_13_0h;
double p_t_14_0h;
double p_t_15_0h;
double p_t_16_0h;
double p_t_17_0h, p_t_17_0m;
double p_t_18_0h, p_t_18_0m;
double p_t_19_0h, p_t_19_0m;
double p_t_20_0h, p_t_20_0m;

 if (xhi < EXTRABOUND2) {
   
   double t1, t2, t3;

   t2 = p0_quick_coeff_3h * x;
   t1 = x * x;
   t3 = t1 * t2;
   
   Add12(*p_resh,*p_resm,x,t3);
   

   return;
 } 

Mul12(&p_x_0_pow2h,&p_x_0_pow2m,x,x);

p_t_15_0h = p0_quick_coeff_5h;
if (xhi > EXTRABOUND) {

p_t_1_0h = p0_quick_coeff_19h;
p_t_2_0h = p_t_1_0h * p_x_0_pow2h;
p_t_3_0h = p0_quick_coeff_17h + p_t_2_0h;
p_t_4_0h = p_t_3_0h * p_x_0_pow2h;
p_t_5_0h = p0_quick_coeff_15h + p_t_4_0h;
p_t_6_0h = p_t_5_0h * p_x_0_pow2h;
p_t_7_0h = p0_quick_coeff_13h + p_t_6_0h;
p_t_8_0h = p_t_7_0h * p_x_0_pow2h;
p_t_9_0h = p0_quick_coeff_11h + p_t_8_0h;
p_t_10_0h = p_t_9_0h * p_x_0_pow2h;
p_t_11_0h = p0_quick_coeff_9h + p_t_10_0h;
p_t_12_0h = p_t_11_0h * p_x_0_pow2h;
p_t_13_0h = p0_quick_coeff_7h + p_t_12_0h;
p_t_14_0h = p_t_13_0h * p_x_0_pow2h;
p_t_15_0h = p_t_15_0h + p_t_14_0h;
}

p_t_16_0h = p_t_15_0h * p_x_0_pow2h;
Add12(p_t_17_0h,p_t_17_0m,p0_quick_coeff_3h,p_t_16_0h);

 Mul122(&p_t_18_0h,&p_t_18_0m,x,p_x_0_pow2h,p_x_0_pow2m);
 Mul22(&p_t_19_0h,&p_t_19_0m,p_t_17_0h,p_t_17_0m,p_t_18_0h,p_t_18_0m);

 Add122(&p_t_20_0h,&p_t_20_0m,x,p_t_19_0h,p_t_19_0m);

*p_resh = p_t_20_0h; *p_resm = p_t_20_0m;


}

static inline void p_quick(double *p_resh, double *p_resm, double x, int index) {
  

double p_t_1_0h;
double p_t_2_0h;
double p_t_3_0h;
double p_t_4_0h;
double p_t_5_0h;
double p_t_6_0h;
double p_t_7_0h;
double p_t_8_0h;
double p_t_9_0h;
double p_t_10_0h;
double p_t_11_0h;
double p_t_12_0h;
double p_t_13_0h;
double p_t_14_0h;
double p_t_15_0h;
double p_t_16_0h;
double p_t_17_0h;
double p_t_18_0h;
double p_t_19_0h;
double p_t_20_0h;
double p_t_21_0h, p_t_21_0m;
double p_t_22_0h, p_t_22_0m;
double p_t_23_0h, p_t_23_0m;
 
p_t_1_0h = p_quick_coeff_12h;
p_t_2_0h = p_t_1_0h * x;
p_t_3_0h = p_quick_coeff_11h + p_t_2_0h;
p_t_4_0h = p_t_3_0h * x;
p_t_5_0h = p_quick_coeff_10h + p_t_4_0h;
p_t_6_0h = p_t_5_0h * x;
p_t_7_0h = p_quick_coeff_9h + p_t_6_0h;
p_t_8_0h = p_t_7_0h * x;
p_t_9_0h = p_quick_coeff_8h + p_t_8_0h;
p_t_10_0h = p_t_9_0h * x;
p_t_11_0h = p_quick_coeff_7h + p_t_10_0h;
p_t_12_0h = p_t_11_0h * x;
p_t_13_0h = p_quick_coeff_6h + p_t_12_0h;
p_t_14_0h = p_t_13_0h * x;
p_t_15_0h = p_quick_coeff_5h + p_t_14_0h;
p_t_16_0h = p_t_15_0h * x;
p_t_17_0h = p_quick_coeff_4h + p_t_16_0h;
p_t_18_0h = p_t_17_0h * x;
p_t_19_0h = p_quick_coeff_3h + p_t_18_0h;
p_t_20_0h = p_t_19_0h * x;
Add12(p_t_21_0h,p_t_21_0m,p_quick_coeff_2h,p_t_20_0h);
MulAdd212(&p_t_22_0h,&p_t_22_0m,p_quick_coeff_1h,p_quick_coeff_1m,x,p_t_21_0h,p_t_21_0m);
MulAdd212(&p_t_23_0h,&p_t_23_0m,p_quick_coeff_0h,p_quick_coeff_0m,x,p_t_22_0h,p_t_22_0m);
*p_resh = p_t_23_0h; *p_resm = p_t_23_0m;

}



static inline void p9_quick(double *p_resh, double *p_resm, double x) {




double p_t_1_0h;
double p_t_2_0h;
double p_t_3_0h;
double p_t_4_0h;
double p_t_5_0h;
double p_t_6_0h;
double p_t_7_0h;
double p_t_8_0h;
double p_t_9_0h;
double p_t_10_0h;
double p_t_11_0h;
double p_t_12_0h;
double p_t_13_0h;
double p_t_14_0h;
double p_t_15_0h;
double p_t_16_0h;
double p_t_17_0h;
double p_t_18_0h;
double p_t_19_0h;
double p_t_20_0h;
double p_t_21_0h, p_t_21_0m;
double p_t_22_0h, p_t_22_0m;
double p_t_23_0h, p_t_23_0m;
 


p_t_1_0h = p9_quick_coeff_11h;
p_t_2_0h = p_t_1_0h * x;
p_t_3_0h = p9_quick_coeff_10h + p_t_2_0h;
p_t_4_0h = p_t_3_0h * x;
p_t_5_0h = p9_quick_coeff_9h + p_t_4_0h;
p_t_6_0h = p_t_5_0h * x;
p_t_7_0h = p9_quick_coeff_8h + p_t_6_0h;
p_t_8_0h = p_t_7_0h * x;
p_t_9_0h = p9_quick_coeff_7h + p_t_8_0h;
p_t_10_0h = p_t_9_0h * x;
p_t_11_0h = p9_quick_coeff_6h + p_t_10_0h;
p_t_12_0h = p_t_11_0h * x;
p_t_13_0h = p9_quick_coeff_5h + p_t_12_0h;
p_t_14_0h = p_t_13_0h * x;
p_t_15_0h = p9_quick_coeff_4h + p_t_14_0h;
p_t_16_0h = p_t_15_0h * x;
p_t_17_0h = p9_quick_coeff_3h + p_t_16_0h;
p_t_18_0h = p_t_17_0h * x;
p_t_19_0h = p9_quick_coeff_2h + p_t_18_0h;
p_t_20_0h = p_t_19_0h * x;
Add12(p_t_21_0h,p_t_21_0m,p9_quick_coeff_1h,p_t_20_0h);
Mul122(&p_t_22_0h,&p_t_22_0m,x,p_t_21_0h,p_t_21_0m);
Add122(&p_t_23_0h,&p_t_23_0m,p9_quick_coeff_0h,p_t_22_0h,p_t_22_0m);
*p_resh = p_t_23_0h; *p_resm = p_t_23_0m;


}



static void p0_accu(double *p_resh, double *p_resm, double *p_resl, double x) {
double p_x_0_pow2h, p_x_0_pow2m;


Mul12(&p_x_0_pow2h,&p_x_0_pow2m,x,x);


double p_t_1_0h;
double p_t_2_0h;
double p_t_3_0h;
double p_t_4_0h;
double p_t_5_0h;
double p_t_6_0h;
double p_t_7_0h, p_t_7_0m;
double p_t_8_0h, p_t_8_0m;
double p_t_9_0h, p_t_9_0m;
double p_t_10_0h, p_t_10_0m;
double p_t_11_0h, p_t_11_0m;
double p_t_12_0h, p_t_12_0m;
double p_t_13_0h, p_t_13_0m;
double p_t_14_0h, p_t_14_0m;
double p_t_15_0h, p_t_15_0m;
double p_t_16_0h, p_t_16_0m;
double p_t_17_0h, p_t_17_0m;
double p_t_18_0h, p_t_18_0m;
double p_t_19_0h, p_t_19_0m;
double p_t_20_0h, p_t_20_0m;
double p_t_21_0h, p_t_21_0m;
double p_t_22_0h, p_t_22_0m;
double p_t_23_0h, p_t_23_0m;
double p_t_24_0h, p_t_24_0m;
double p_t_25_0h, p_t_25_0m;
double p_t_26_0h, p_t_26_0m;
double p_t_27_0h, p_t_27_0m;
double p_t_28_0h, p_t_28_0m, p_t_28_0l;
double p_t_29_0h, p_t_29_0m, p_t_29_0l;
double p_t_30_0h, p_t_30_0m, p_t_30_0l;
double p_t_31_0h, p_t_31_0m, p_t_31_0l;
double p_t_32_0h, p_t_32_0m, p_t_32_0l;
double p_t_33_0h, p_t_33_0m, p_t_33_0l;
double p_t_34_0h, p_t_34_0m, p_t_34_0l;
double p_t_35_0h, p_t_35_0m, p_t_35_0l;
 


p_t_1_0h = p0_accu_coeff_37h;
p_t_2_0h = p_t_1_0h * p_x_0_pow2h;
p_t_3_0h = p0_accu_coeff_35h + p_t_2_0h;
p_t_4_0h = p_t_3_0h * p_x_0_pow2h;
p_t_5_0h = p0_accu_coeff_33h + p_t_4_0h;
p_t_6_0h = p_t_5_0h * p_x_0_pow2h;
Add12(p_t_7_0h,p_t_7_0m,p0_accu_coeff_31h,p_t_6_0h);
Mul22(&p_t_8_0h,&p_t_8_0m,p_t_7_0h,p_t_7_0m,p_x_0_pow2h,p_x_0_pow2m);
Add122(&p_t_9_0h,&p_t_9_0m,p0_accu_coeff_29h,p_t_8_0h,p_t_8_0m);
Mul22(&p_t_10_0h,&p_t_10_0m,p_t_9_0h,p_t_9_0m,p_x_0_pow2h,p_x_0_pow2m);
Add122(&p_t_11_0h,&p_t_11_0m,p0_accu_coeff_27h,p_t_10_0h,p_t_10_0m);
Mul22(&p_t_12_0h,&p_t_12_0m,p_t_11_0h,p_t_11_0m,p_x_0_pow2h,p_x_0_pow2m);
Add122(&p_t_13_0h,&p_t_13_0m,p0_accu_coeff_25h,p_t_12_0h,p_t_12_0m);
Mul22(&p_t_14_0h,&p_t_14_0m,p_t_13_0h,p_t_13_0m,p_x_0_pow2h,p_x_0_pow2m);
Add122(&p_t_15_0h,&p_t_15_0m,p0_accu_coeff_23h,p_t_14_0h,p_t_14_0m);
Mul22(&p_t_16_0h,&p_t_16_0m,p_t_15_0h,p_t_15_0m,p_x_0_pow2h,p_x_0_pow2m);
Add122(&p_t_17_0h,&p_t_17_0m,p0_accu_coeff_21h,p_t_16_0h,p_t_16_0m);
Mul22(&p_t_18_0h,&p_t_18_0m,p_t_17_0h,p_t_17_0m,p_x_0_pow2h,p_x_0_pow2m);
Add122(&p_t_19_0h,&p_t_19_0m,p0_accu_coeff_19h,p_t_18_0h,p_t_18_0m);
Mul22(&p_t_20_0h,&p_t_20_0m,p_t_19_0h,p_t_19_0m,p_x_0_pow2h,p_x_0_pow2m);
Add122(&p_t_21_0h,&p_t_21_0m,p0_accu_coeff_17h,p_t_20_0h,p_t_20_0m);
Mul22(&p_t_22_0h,&p_t_22_0m,p_t_21_0h,p_t_21_0m,p_x_0_pow2h,p_x_0_pow2m);
Add122(&p_t_23_0h,&p_t_23_0m,p0_accu_coeff_15h,p_t_22_0h,p_t_22_0m);
MulAdd22(&p_t_24_0h,&p_t_24_0m,p0_accu_coeff_13h,p0_accu_coeff_13m,p_x_0_pow2h,p_x_0_pow2m,p_t_23_0h,p_t_23_0m);
MulAdd22(&p_t_25_0h,&p_t_25_0m,p0_accu_coeff_11h,p0_accu_coeff_11m,p_x_0_pow2h,p_x_0_pow2m,p_t_24_0h,p_t_24_0m);
MulAdd22(&p_t_26_0h,&p_t_26_0m,p0_accu_coeff_9h,p0_accu_coeff_9m,p_x_0_pow2h,p_x_0_pow2m,p_t_25_0h,p_t_25_0m);
Mul22(&p_t_27_0h,&p_t_27_0m,p_t_26_0h,p_t_26_0m,p_x_0_pow2h,p_x_0_pow2m);
Add23(&p_t_28_0h,&p_t_28_0m,&p_t_28_0l,p0_accu_coeff_7h,p0_accu_coeff_7m,p_t_27_0h,p_t_27_0m);
Mul233(&p_t_29_0h,&p_t_29_0m,&p_t_29_0l,p_x_0_pow2h,p_x_0_pow2m,p_t_28_0h,p_t_28_0m,p_t_28_0l);
Add233(&p_t_30_0h,&p_t_30_0m,&p_t_30_0l,p0_accu_coeff_5h,p0_accu_coeff_5m,p_t_29_0h,p_t_29_0m,p_t_29_0l);
Mul233(&p_t_31_0h,&p_t_31_0m,&p_t_31_0l,p_x_0_pow2h,p_x_0_pow2m,p_t_30_0h,p_t_30_0m,p_t_30_0l);
Add233(&p_t_32_0h,&p_t_32_0m,&p_t_32_0l,p0_accu_coeff_3h,p0_accu_coeff_3m,p_t_31_0h,p_t_31_0m,p_t_31_0l);
Mul233(&p_t_33_0h,&p_t_33_0m,&p_t_33_0l,p_x_0_pow2h,p_x_0_pow2m,p_t_32_0h,p_t_32_0m,p_t_32_0l);
Add133(&p_t_34_0h,&p_t_34_0m,&p_t_34_0l,p0_accu_coeff_1h,p_t_33_0h,p_t_33_0m,p_t_33_0l);
Mul133(&p_t_35_0h,&p_t_35_0m,&p_t_35_0l,x,p_t_34_0h,p_t_34_0m,p_t_34_0l);
Renormalize3(p_resh,p_resm,p_resl,p_t_35_0h,p_t_35_0m,p_t_35_0l);


}


static void p_accu(double *p_resh, double *p_resm, double *p_resl, double x, int index) {
  

double p_t_1_0h;
double p_t_2_0h;
double p_t_3_0h;
double p_t_4_0h;
double p_t_5_0h;
double p_t_6_0h;
double p_t_7_0h;
double p_t_8_0h;
double p_t_9_0h;
double p_t_10_0h;
double p_t_11_0h;
double p_t_12_0h;
double p_t_13_0h;
double p_t_14_0h;
double p_t_15_0h;
double p_t_16_0h;
double p_t_17_0h, p_t_17_0m;
double p_t_18_0h, p_t_18_0m;
double p_t_19_0h, p_t_19_0m;
double p_t_20_0h, p_t_20_0m;
double p_t_21_0h, p_t_21_0m;
double p_t_22_0h, p_t_22_0m;
double p_t_23_0h, p_t_23_0m;
double p_t_24_0h, p_t_24_0m;
double p_t_25_0h, p_t_25_0m;
double p_t_26_0h, p_t_26_0m;
double p_t_27_0h, p_t_27_0m;
double p_t_28_0h, p_t_28_0m;
double p_t_29_0h, p_t_29_0m;
double p_t_30_0h, p_t_30_0m;
double p_t_31_0h, p_t_31_0m;
double p_t_32_0h, p_t_32_0m;
double p_t_33_0h, p_t_33_0m;
double p_t_34_0h, p_t_34_0m;
double p_t_35_0h, p_t_35_0m, p_t_35_0l;
double p_t_36_0h, p_t_36_0m, p_t_36_0l;
double p_t_37_0h, p_t_37_0m, p_t_37_0l;
double p_t_38_0h, p_t_38_0m, p_t_38_0l;
double p_t_39_0h, p_t_39_0m, p_t_39_0l;
 


p_t_1_0h = p_accu_coeff_22h;
p_t_2_0h = p_t_1_0h * x;
p_t_3_0h = p_accu_coeff_21h + p_t_2_0h;
p_t_4_0h = p_t_3_0h * x;
p_t_5_0h = p_accu_coeff_20h + p_t_4_0h;
p_t_6_0h = p_t_5_0h * x;
p_t_7_0h = p_accu_coeff_19h + p_t_6_0h;
p_t_8_0h = p_t_7_0h * x;
p_t_9_0h = p_accu_coeff_18h + p_t_8_0h;
p_t_10_0h = p_t_9_0h * x;
p_t_11_0h = p_accu_coeff_17h + p_t_10_0h;
p_t_12_0h = p_t_11_0h * x;
p_t_13_0h = p_accu_coeff_16h + p_t_12_0h;
p_t_14_0h = p_t_13_0h * x;
p_t_15_0h = p_accu_coeff_15h + p_t_14_0h;
p_t_16_0h = p_t_15_0h * x;
Add12(p_t_17_0h,p_t_17_0m,p_accu_coeff_14h,p_t_16_0h);
Mul122(&p_t_18_0h,&p_t_18_0m,x,p_t_17_0h,p_t_17_0m);
Add122(&p_t_19_0h,&p_t_19_0m,p_accu_coeff_13h,p_t_18_0h,p_t_18_0m);
Mul122(&p_t_20_0h,&p_t_20_0m,x,p_t_19_0h,p_t_19_0m);
Add122(&p_t_21_0h,&p_t_21_0m,p_accu_coeff_12h,p_t_20_0h,p_t_20_0m);
Mul122(&p_t_22_0h,&p_t_22_0m,x,p_t_21_0h,p_t_21_0m);
Add122(&p_t_23_0h,&p_t_23_0m,p_accu_coeff_11h,p_t_22_0h,p_t_22_0m);
Mul122(&p_t_24_0h,&p_t_24_0m,x,p_t_23_0h,p_t_23_0m);
Add122(&p_t_25_0h,&p_t_25_0m,p_accu_coeff_10h,p_t_24_0h,p_t_24_0m);
Mul122(&p_t_26_0h,&p_t_26_0m,x,p_t_25_0h,p_t_25_0m);
Add122(&p_t_27_0h,&p_t_27_0m,p_accu_coeff_9h,p_t_26_0h,p_t_26_0m);
MulAdd212(&p_t_28_0h,&p_t_28_0m,p_accu_coeff_8h,p_accu_coeff_8m,x,p_t_27_0h,p_t_27_0m);
MulAdd212(&p_t_29_0h,&p_t_29_0m,p_accu_coeff_7h,p_accu_coeff_7m,x,p_t_28_0h,p_t_28_0m);
MulAdd212(&p_t_30_0h,&p_t_30_0m,p_accu_coeff_6h,p_accu_coeff_6m,x,p_t_29_0h,p_t_29_0m);
MulAdd212(&p_t_31_0h,&p_t_31_0m,p_accu_coeff_5h,p_accu_coeff_5m,x,p_t_30_0h,p_t_30_0m);
MulAdd212(&p_t_32_0h,&p_t_32_0m,p_accu_coeff_4h,p_accu_coeff_4m,x,p_t_31_0h,p_t_31_0m);
MulAdd212(&p_t_33_0h,&p_t_33_0m,p_accu_coeff_3h,p_accu_coeff_3m,x,p_t_32_0h,p_t_32_0m);
Mul122(&p_t_34_0h,&p_t_34_0m,x,p_t_33_0h,p_t_33_0m);
Add23(&p_t_35_0h,&p_t_35_0m,&p_t_35_0l,p_accu_coeff_2h,p_accu_coeff_2m,p_t_34_0h,p_t_34_0m);
Mul133(&p_t_36_0h,&p_t_36_0m,&p_t_36_0l,x,p_t_35_0h,p_t_35_0m,p_t_35_0l);
Add233(&p_t_37_0h,&p_t_37_0m,&p_t_37_0l,p_accu_coeff_1h,p_accu_coeff_1m,p_t_36_0h,p_t_36_0m,p_t_36_0l);
Mul133(&p_t_38_0h,&p_t_38_0m,&p_t_38_0l,x,p_t_37_0h,p_t_37_0m,p_t_37_0l);
Add233(&p_t_39_0h,&p_t_39_0m,&p_t_39_0l,p_accu_coeff_0h,p_accu_coeff_0m,p_t_38_0h,p_t_38_0m,p_t_38_0l);
Renormalize3(p_resh,p_resm,p_resl,p_t_39_0h,p_t_39_0m,p_t_39_0l);


}


static void p9_accu(double *p_resh, double *p_resm, double *p_resl, double x) {




double p_t_1_0h;
double p_t_2_0h;
double p_t_3_0h;
double p_t_4_0h;
double p_t_5_0h;
double p_t_6_0h;
double p_t_7_0h;
double p_t_8_0h;
double p_t_9_0h;
double p_t_10_0h;
double p_t_11_0h;
double p_t_12_0h;
double p_t_13_0h;
double p_t_14_0h;
double p_t_15_0h;
double p_t_16_0h;
double p_t_17_0h, p_t_17_0m;
double p_t_18_0h, p_t_18_0m;
double p_t_19_0h, p_t_19_0m;
double p_t_20_0h, p_t_20_0m;
double p_t_21_0h, p_t_21_0m;
double p_t_22_0h, p_t_22_0m;
double p_t_23_0h, p_t_23_0m;
double p_t_24_0h, p_t_24_0m;
double p_t_25_0h, p_t_25_0m;
double p_t_26_0h, p_t_26_0m;
double p_t_27_0h, p_t_27_0m;
double p_t_28_0h, p_t_28_0m;
double p_t_29_0h, p_t_29_0m;
double p_t_30_0h, p_t_30_0m;
double p_t_31_0h, p_t_31_0m;
double p_t_32_0h, p_t_32_0m;
double p_t_33_0h, p_t_33_0m, p_t_33_0l;
double p_t_34_0h, p_t_34_0m, p_t_34_0l;
double p_t_35_0h, p_t_35_0m, p_t_35_0l;
 


p_t_1_0h = p9_accu_coeff_20h;
p_t_2_0h = p_t_1_0h * x;
p_t_3_0h = p9_accu_coeff_19h + p_t_2_0h;
p_t_4_0h = p_t_3_0h * x;
p_t_5_0h = p9_accu_coeff_18h + p_t_4_0h;
p_t_6_0h = p_t_5_0h * x;
p_t_7_0h = p9_accu_coeff_17h + p_t_6_0h;
p_t_8_0h = p_t_7_0h * x;
p_t_9_0h = p9_accu_coeff_16h + p_t_8_0h;
p_t_10_0h = p_t_9_0h * x;
p_t_11_0h = p9_accu_coeff_15h + p_t_10_0h;
p_t_12_0h = p_t_11_0h * x;
p_t_13_0h = p9_accu_coeff_14h + p_t_12_0h;
p_t_14_0h = p_t_13_0h * x;
p_t_15_0h = p9_accu_coeff_13h + p_t_14_0h;
p_t_16_0h = p_t_15_0h * x;
Add12(p_t_17_0h,p_t_17_0m,p9_accu_coeff_12h,p_t_16_0h);
Mul122(&p_t_18_0h,&p_t_18_0m,x,p_t_17_0h,p_t_17_0m);
Add122(&p_t_19_0h,&p_t_19_0m,p9_accu_coeff_11h,p_t_18_0h,p_t_18_0m);
Mul122(&p_t_20_0h,&p_t_20_0m,x,p_t_19_0h,p_t_19_0m);
Add122(&p_t_21_0h,&p_t_21_0m,p9_accu_coeff_10h,p_t_20_0h,p_t_20_0m);
Mul122(&p_t_22_0h,&p_t_22_0m,x,p_t_21_0h,p_t_21_0m);
Add122(&p_t_23_0h,&p_t_23_0m,p9_accu_coeff_9h,p_t_22_0h,p_t_22_0m);
Mul122(&p_t_24_0h,&p_t_24_0m,x,p_t_23_0h,p_t_23_0m);
Add122(&p_t_25_0h,&p_t_25_0m,p9_accu_coeff_8h,p_t_24_0h,p_t_24_0m);
MulAdd212(&p_t_26_0h,&p_t_26_0m,p9_accu_coeff_7h,p9_accu_coeff_7m,x,p_t_25_0h,p_t_25_0m);
MulAdd212(&p_t_27_0h,&p_t_27_0m,p9_accu_coeff_6h,p9_accu_coeff_6m,x,p_t_26_0h,p_t_26_0m);
MulAdd212(&p_t_28_0h,&p_t_28_0m,p9_accu_coeff_5h,p9_accu_coeff_5m,x,p_t_27_0h,p_t_27_0m);
MulAdd212(&p_t_29_0h,&p_t_29_0m,p9_accu_coeff_4h,p9_accu_coeff_4m,x,p_t_28_0h,p_t_28_0m);
MulAdd212(&p_t_30_0h,&p_t_30_0m,p9_accu_coeff_3h,p9_accu_coeff_3m,x,p_t_29_0h,p_t_29_0m);
MulAdd212(&p_t_31_0h,&p_t_31_0m,p9_accu_coeff_2h,p9_accu_coeff_2m,x,p_t_30_0h,p_t_30_0m);
Mul122(&p_t_32_0h,&p_t_32_0m,x,p_t_31_0h,p_t_31_0m);
Add23(&p_t_33_0h,&p_t_33_0m,&p_t_33_0l,p9_accu_coeff_1h,p9_accu_coeff_1m,p_t_32_0h,p_t_32_0m);
Mul133(&p_t_34_0h,&p_t_34_0m,&p_t_34_0l,x,p_t_33_0h,p_t_33_0m,p_t_33_0l);
Add233(&p_t_35_0h,&p_t_35_0m,&p_t_35_0l,p9_accu_coeff_0h,p9_accu_coeff_0m,p_t_34_0h,p_t_34_0m,p_t_34_0l);
Renormalize3(p_resh,p_resm,p_resl,p_t_35_0h,p_t_35_0m,p_t_35_0l);


}




double asin_rn(double x) {
  db_number xdb, zdb;
  double sign, z, zp;
  int index;
  double asinh, asinm, asinl;
  double p9h, p9m, p9l, sqrh, sqrm, sqrl;
  double t1h, t1m, t1l;
  double asin;
  double xabs;

  /* Start already computations for argument reduction */

  zdb.d = 1.0 + x * x;

  xdb.d = x;

  /* Special case handling */
  
  /* Remove sign of x in floating-point */
  xabs = ABS(x);
  xdb.i[HI] &= 0x7fffffff;

  /* If |x| < 2^(-28) we have
     
     arcsin(x) = x * ( 1 + xi ) 

     with 0 <= xi < 2^(-55) 
          
     So we can decide the rounding without any computation 
  */
  if (xdb.i[HI] < SIMPLEBOUND) {
    return x;
  }

  /* asin is defined on -1 <= x <= 1, elsewhere it is NaN */
  if (xdb.i[HI] >= 0x3ff00000) {
    if (x == 1.0) {
      return PIHALFH;
    }
    if (x == -1.0) {
      return - PIHALFH;
    }
    return (x-x)/0.0;    /* return NaN */
  }

  /* Argument reduction:

     We have 10 intervals and 3 paths:

     - interval 0   => path 1 using p0
     - interval 1-8 => path 2 using p
     - interval 9   => path 3 using p9

  */

  index = (0x000f0000 & zdb.i[HI]) >> 16;

  /* 0 <= index <= 15 

     index approximates roughly x^2 

     Map indexes to intervals as follows:

     0  -> 0 
     1  -> 1
     ... 
     8  -> 8
     9  -> 9
     ... 
     15 -> 9

     For this mapping, filter first the case 0 -> 0
     In consequence, 1 <= index <= 15, i.e. 
     0 <= index - 1 <= 14 with the mapping index - 1 -> interval as

     0  -> 1
     ... 
     7  -> 8
     8  -> 9
     ...
     15 -> 9

     Thus it suffices to check the 3rd bit of index - 1 after the first filter.
     
  */

  if (index == 0) {
    /* Path 1 using p0 */

    p0_quick(&asinh, &asinm, x, xdb.i[HI]);

    /* Rounding test */

    if(asinh == (asinh + (asinm * RNROUNDCST))) 
      return asinh;

    /* Rounding test failed, launch accurate phase */

#if EVAL_PERF
  crlibm_second_step_taken++;
#endif
    
    p0_accu(&asinh, &asinm, &asinl, x);

    /* Final rounding */

    RoundToNearest3(&asin,asinh,asinm,asinl);

    return asin;
  } 

  /* Strip off the sign of argument x */
  sign = 1.0;
  if (x < 0.0) sign = -sign;
  
  index--;
  if ((index & 0x8) != 0) {
    /* Path 3 using p9 */

    /* Do argument reduction using a MI_9 as a midpoint value 
       for the polynomial and compute exactly zp = 2 * (1 - x) 
       for the asymptotical approximation using a square root.
    */

    z = xabs - MI_9;
    zp = 2.0 * (1.0 - xabs);

    /* Polynomial approximation and square root extraction */

    p9_quick(&p9h, &p9m, z);
    p9h = -p9h;
    p9m = -p9m;

    sqrt12_64_unfiltered(&sqrh,&sqrm,zp);

    /* Reconstruction */

    Mul22(&t1h,&t1m,sqrh,sqrm,p9h,p9m);
    Add22(&asinh,&asinm,PIHALFH,PIHALFM,t1h,t1m);

    /* Rounding test */

    if(asinh == (asinh + (asinm * RNROUNDCST))) 
      return sign * asinh;

    /* Rounding test failed, launch accurate phase */

#if EVAL_PERF
  crlibm_second_step_taken++;
#endif
    
    p9_accu(&p9h, &p9m, &p9l, z);
    p9h = -p9h;
    p9m = -p9m;
    p9l = -p9l;

    Sqrt13(&sqrh,&sqrm,&sqrl,zp);

    /* Reconstruction */

    Mul33(&t1h,&t1m,&t1l,sqrh,sqrm,sqrl,p9h,p9m,p9l);
    Add33(&asinh,&asinm,&asinl,PIHALFH,PIHALFM,PIHALFL,t1h,t1m,t1l);

    /* Final rounding */    

    RoundToNearest3(&asin,asinh,asinm,asinl);

    return sign * asin;

  }

  /* Path 2 using p */

  /* Do argument reduction using a table value for 
     the midpoint value 
  */

  z = xabs - mi_i;

  p_quick(&asinh, &asinm, z, index);


  /* Rounding test */
  
  if(asinh == (asinh + (asinm * RNROUNDCST))) 
    return sign * asinh;
  
  /* Rounding test failed, launch accurate phase */
  
#if EVAL_PERF
  crlibm_second_step_taken++;
#endif
  
  p_accu(&asinh, &asinm, &asinl, z, index);
  
  /* Final rounding */
  
  RoundToNearest3(&asin,asinh,asinm,asinl);
  
  return sign * asin;
  
}


double asin_ru(double x) {
  return asin_rn(x);
}

double asin_rd(double x) {
  return asin_rn(x);
}

double asin_rz(double x) {
  return asin_rn(x);
}
