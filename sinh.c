/*
 *  sinh.c
 *  
 *
 *  Created by Catherine Daramy on Fri Mar 14 2003.
 *  Copyright (c) 2003 __MyCompanyName__. All rights reserved.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include "sinh.h"
#include "exp.h"

/** 
 *  
 *  
 *  
 *  
 *  
 *  
 *  
 *  
 *  
 *  
 *  
 *  
 *  
 *  
 *  
 *  
 *  
 *  
 *  
 *  
 *  
 *  
 */
 
void exp_SC(scs_ptr, scs_ptr);

 
 /************************************************************
 *************************************************************
 *               ROUNDED  TO NEAREST		             *
 *************************************************************
 *************************************************************/
double  sinh_rn(double x){ 
db_number res;
scs_t sh1, sh2, x_scs, res_scs;
int sign = 0;

if(x < 0){
    sign = 1;
    x *=-1;
}
res.d = x;
scs_set_d(x_scs, x);   


    if(res.i[HI_ENDIAN] > 0xC9F34F7925F50D07D80800000U){	/* look if x < 2^-40   */
	printf("calcul 1\n");
	scs_square(res_scs, x_scs);
	scs_mul(res_scs, x_scs, res_scs);
	scs_mul(sh2, SCS_SIXINV, res_scs);
	scs_add(res_scs, sh2, x_scs);			/* return with accuracy 2^-167  */
    }
    else if(res.i[HI_ENDIAN] < mediumconst){			/* look if x < 55   */
	exp_SC(sh1, x_scs);					
	printf("calcul 2\ncalcul de exp(x):\n");
	scs_get_std(sh1);	
	scs_inv(sh2, sh1);
	scs_sub(res_scs, sh1, sh2);
	scs_mul(res_scs, res_scs, SCS_HALF);
	printf("res_scs = \n");
	scs_get_std(res_scs);
    }
    else if(res.i[0] > larg_int){			/* return overflow */
	scs_set_d(res_scs, radix_rng_double.d*radix_rng_double.d);  
    }
    else{						/* 55 < x < 710.47... = bigconst */
	printf("calcul 3\n");
	exp_SC(res_scs, x_scs);				/* return with accuracy 2^-158  */
	scs_mul(res_scs, res_scs, SCS_HALF);
    }
    
scs_get_d(&res.d, res_scs);

if(sign == 1){
    res.d *= -1;
}

return res.d;

}
   
 /************************************************************
 *************************************************************
 *               ROUNDED  TO - INFINITY		             *
 *************************************************************
 *************************************************************/
double sinh_rd(double x){ 
db_number res;
scs_t sh1, sh2, x_scs, res_scs;
int sign = 0;

if(x < 0){
    sign = 1;
    x *=-1;
}
res.d = x;
scs_set_d(x_scs, x);   


    if((res.l & 0x7fffffffffffffff) < 0x50400000U){	/* look if x < 2^-40   */
	printf("calcul 1\n");
	scs_square(res_scs, x_scs);
	scs_mul(res_scs, x_scs, res_scs);
	scs_mul(sh2, SCS_SIXINV, res_scs);
	scs_add(res_scs, sh2, x_scs);			/* return with accuracy 2^-167  */
    }
    else if(res.d < 55){				/* look if x < 55   */
	exp_SC(sh1, x_scs);					
	printf("calcul 2\ncalcul de exp(x):\n");
	scs_get_std(sh1);	
	scs_inv(sh2, sh1);
	scs_sub(res_scs, sh1, sh2);
	scs_mul(res_scs, res_scs, SCS_HALF);
	printf("res_scs = \n");
	scs_get_std(res_scs);
    }
    else if(res.d > 710.475860073943942037){			/* return overflow */
	printf("calcul impossible, x trop grand\n");
	scs_set_d(res_scs, radix_rng_double.d*radix_rng_double.d);  
    }
    else{							/* 55 < x < 710.47... = bigconst */
	printf("calcul 3\n");
	exp_SC(res_scs, x_scs);					/* return with accuracy 2^-158  */
	scs_mul(res_scs, res_scs, SCS_HALF);
    }
    
if (sign == 1){
    scs_get_d_pinf(&res.d, res_scs);
     res.d *= -1;
}
else{
    scs_get_d_minf(&res.d, res_scs);
}

return res.d;
} 


 /************************************************************
 *************************************************************
 *               ROUNDED  TO + INFINITY		             *
 *************************************************************
 *************************************************************/
double sinh_ru(double x){ 
db_number res;
scs_t sh1, sh2, x_scs, res_scs;
int sign = 0;

if(x < 0){
    sign = 1;
    x *=-1;
}
res.d = x;
scs_set_d(x_scs, x);   


    if((res.l & 0x7fffffffffffffff) < 0x50400000U){			/* look if x < 2^-40   */
	printf("calcul 1\n");
	scs_square(res_scs, x_scs);
	scs_mul(res_scs, x_scs, res_scs);
	scs_mul(sh2, SCS_SIXINV, res_scs);
	scs_add(res_scs, sh2, x_scs);			/* return with accuracy 2^-167  */
    }
    else if(res.d < 55){				/* look if x < 35   */
	exp_SC(sh1, x_scs);					
	printf("calcul 2\ncalcul de exp(x):\n");
	scs_get_std(sh1);	
	scs_inv(sh2, sh1);
	scs_sub(res_scs, sh1, sh2);
	scs_mul(res_scs, res_scs, SCS_HALF);
	printf("res_scs = \n");
	scs_get_std(res_scs);
    }
    else if(res.d > 710.475860073943942037){			/* return overflow */
	printf("calcul impossible, x trop grand\n");
	scs_set_d(res_scs, radix_rng_double.d*radix_rng_double.d);  
    }
    else{							/* 55 < x < 710.47... = bigconst */
	printf("calcul 3\n");
	exp_SC(res_scs, x_scs);					/* return with accuracy 2^-158  */
	scs_mul(res_scs, res_scs, SCS_HALF);
    }
    
if (sign == 1){
    scs_get_d_minf(&res.d, res_scs);
     res.d *= -1;
}
else{
    scs_get_d_pinf(&res.d, res_scs);
}

return res.d;
}  
