/*
 *  atan.c
 *  
 *
 *  Created by cathydar on Thu Nov 07 2002.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <atan.h>

void atan(scs_ptr, scs_ptr);
double scs_atan_rd(double);
double scs_atan_ru(double);
double scs_atan_rn(double);

/*
 *  WHAT WE CAN DO :
 *
 * 1) Range reduction 
 *
 *	x > 0  because atan(-x) = - atan(x)
 *	
 *	we have built 50 intervals I(i), associated to a b(i) so that :
 *	
 *	For every x :
 *	
 *	we find the interval I(i) , as atan(x) = atan(b(i)) + atan( (x - b(i)) / (1 + x * b(i)) ) 
 *	
 *		so that X = (x - b(i)) / (1 + x * b(i))  be in interval [ -2^(-6) , 2^(-6) ] 
 *		There is no cancellation because :
 *		for every x in [ -2^(-6) , 2^(-6) ],
 *		
 *					     atan(x) <= 0.01562372862     in binary 0.000001111111111
 *		AND for the smallest b(i)    atan(b(i)) = 0.04687118592   in binary 0.00001011111111 	   
 *
 *
 * 2) Polynomial evaluation of atan(X), atan(b(i)) is tabulated.
 *                  
 *                                  (-???)  
 *   Approximation error: |err| < 2^ 
 *
 *
 * 3) Reconstruction:
 *
 *    atan(x) = atan(b(i)) + atan(X) 
 *
 *
 * 4) Rounding:
 *
 *    when |x| is too big, the result is always sign(x) * Pi/2,
 *    because Pi/2 is appromated by the biggest value smallest than Pi/2, 
 *    in order not to have an atan > Pi/2.
 */


 
 
void scs_atan(scs_ptr res_scs, scs_ptr x){
    scs_t X_scs, denom1_scs, denom2_scs, num_scs, poly_scs, X2;
    db_number db;
    int k, ind = 100, j, intmid = 32;
  
  scs_get_d(&db.d, x);  
			     	     
    /* Test if x need to be reduced */
    
    if( db.i[0] >= borne_I[48]){
      	if ( db.i[0]  < borne_I[49]){
	    j = 48;}
	else if( db.i[0] > borne_I[49]){
	    j = 49;}
	else {
	    res_scs =  atan_bi_ptr[49];
	    return;
	}
	ind = j;
	
	/* evaluate X = (x - b(j)) / (1 + x*b(j)) */
	    	    	    
	scs_mul(denom1_scs,bsc_ptr[ind],x);
	scs_add(denom2_scs,denom1_scs,SCS_ONE);
	scs_sub(X_scs,x,bsc_ptr[ind]);
	scs_div(X_scs,X_scs,denom2_scs);
    }
    else if (((int) db.i[0]) == borne_I[48]) {
	res_scs = atan_bi_ptr[48];
	return ;
    }
    else if ( db.i[0] < borne_I[0]){
	scs_set(X_scs, x);
	ind = 60;
    }
    else{	/* First reduction : find the interval including x, "save" j to have b(j)  and then being able to calculate X */
	j = 32;
	for (k=1;k<=5;k++){
	    if (db.i[0] < borne_I[j]){
		j -= (intmid >> k);
            }
            else if(db.i[0] > borne_I[j]){
		j += (intmid >>  k);	    
	    }
	    else {res_scs = atan_bi_ptr[j];
		  return ;
	    }
	}
	if (db.i[0] < borne_I[j]){
	    ind = j-1;
	}
	else{
	    ind = j;
	}
	
    /* evaluate X = (x - b(j)) / (1 + x*b(j)) */
		    	    
	scs_mul(denom1_scs,bsc_ptr[ind],x);
	scs_add(denom2_scs,denom1_scs,SCS_ONE);   
	scs_sub(num_scs,x,bsc_ptr[ind]);
	scs_div(X_scs,num_scs,denom2_scs);
    } 
		
/* Polynomial evaluation of atan(X) , X = (x-b(i)) / (1+ x*b(i)) */ 
	
    scs_square(X2, X_scs);
    scs_set(res_scs, constant_poly_ptr[0]);
    for(k=1; k < 13; k++){
	scs_mul(res_scs, res_scs, X2);		/* we use Horner expression */
    	scs_add(res_scs, constant_poly_ptr[k], res_scs);
	
    }
    scs_mul(poly_scs, res_scs, X_scs);
    
    if(ind == 60){ 
	scs_set(res_scs, poly_scs);
	return;
    }else{
	/*scs_set(denom_scs, res_scs);*/
	scs_add(res_scs,atan_bi_ptr[ind], poly_scs); 
    }
    return;
}	
       
/*************************************************************
 *************************************************************
 *               ROUNDED  TO NEAREST
 *************************************************************
 *************************************************************/

double scs_atan_rn(double x){ 
  scs_t sc1;
  scs_t res_scs;
  db_number res;
  int sign =1;
  
   res.d = x;

    /* Filter cases */
    if((res.l & 0x7fffffffffffffffULL) > larg_int.l) 
    {
	if (x > 0){
	    return pio2.d;
	}
	else{
            return mpio2.d;
	}
    }
    else if ((res.l & 0x7fffffffffffffffULL) < tiny_int.l) {
	    return 0.0;
    }
    else{
	 if (x < 0){
	    sign = -1;
	    x *= -1;
	}
	scs_set_d(sc1, x);
	scs_atan(res_scs, sc1);
	scs_get_d(&res.d, res_scs);
	if (sign == -1){
	    res.d *= -1;
	}
	return res.d;
    }
}

/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  -INFINITY
 *************************************************************
 *************************************************************/

double scs_atan_rd(double x){ 
  scs_t sc1;
  scs_t res_scs;
  db_number res;
  int sign = 1;
   
  res.d = x;

    /* Filter cases */
    if((res.l & 0x7fffffffffffffffULL) > larg_int.l) 
    {
	if (x < 0){
	    return mpio2.d;
	}
	else{
		return pio2.d;
	}
    }
    else {  
	if (x < 0){
	    sign = -1;
	    x *= -1;
	}
	scs_set_d(sc1, x);
	scs_atan(res_scs, sc1);
	if (sign == -1){
	    scs_get_d_pinf(&res.d, res_scs);
	    res.d *= -1;
	    return res.d;
	}
	else{
	    scs_get_d_minf(&res.d, res_scs);		
	    return res.d;
	}
    }
}

/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  +INFINITY
 *************************************************************
 *************************************************************/

double scs_atan_ru(double x){ 
  scs_t sc1;
  scs_t res_scs;
  db_number res;
  int sign = 1;
  
  res.d = x;

    /* Filter cases */
    if((res.l & 0x7fffffffffffffffULL) > larg_int.l) 
    {
	if (x > 0)
        {
	    return pio2.d;
	}
	else
	{
            return mpio2.d; 		
	}
    }
    else 
	{
	   if (x < 0){
		sign = -1;
		x *= -1;
	    }
 
	    scs_set_d(sc1, x);
	    scs_atan(res_scs, sc1);
	    	    if (sign == -1){
		scs_get_d_minf(&res.d, res_scs);
		res.d *= -1;
		return res.d;
	    }
	    else{
		scs_get_d_pinf(&res.d, res_scs);		
		return res.d;
	    }
	}
}

