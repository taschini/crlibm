/** Functions for SCS multiplication operations 
@file multiplication_scs.c

@author Defour David David.Defour@ens-lyon.fr
@author Florent de Dinechin Florent.de.Dinechin@ens-lyon.fr 

 This file is part of the SCS library.
*/

/*
Copyright (C) 2002  David Defour and Florent de Dinechin

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/

#include "scs.h"
#include "scs_private.h"

#if 0 /* used to help debugging */
void pr(char* s,double d) {
  db_number x;
  x.d=d;
  printf(s);printf("   ");
  printf("%8x%8x . 2^%d   (%8f %8x %8x)   \n",
	 (x.i[HI_ENDIAN]&0x000FFFFF)+0x00100000,
	 x.i[LO_ENDIAN],
	 (x.i[HI_ENDIAN]>>20)-1023,
	 x.d,
	 x.i[HI_ENDIAN],
	 x.i[LO_ENDIAN]);
}
#endif


/***************************************************************/
/* First a few #defines selecting between architectures where integer
   multiplication is faster than FP multiplication (the normal case)
   and those where it is the opposite (some Sparcs). The latter
   doesn't work but might at some point, so the code remains */

#ifdef SCS_USE_FLT_MULT
/* Compute the carry of r1, remove it from r1, and add it to r0. This
 doesn't work: how to get a value truncated using only
 round-to-nearest FP ? Besides it doesn't work on x86 for at least two
 reasons : extended precision FP arithmetic, and a #define below. */
#define SCS_CARRY_PROPAGATE(r1,r0,tmp) \
     {tmp=((r1-SCS_FLT_TRUNC_CST)+SCS_FLT_SHIFT_CST)-SCS_FLT_SHIFT_CST; r1-= tmp; r0+=tmp*SCS_RADIX_MONE_DOUBLE;}      

typedef double SCS_CONVERSION_MUL;
 
#else /* SCS_USE_FLT_MULT*/

/*  Compute the carry of r1, remove it from r1, and add it to r0 */
#define SCS_CARRY_PROPAGATE(r1,r0,tmp) \
      {tmp = r1>>SCS_NB_BITS; r0 += tmp; r1 -= (tmp<<SCS_NB_BITS);}     
typedef unsigned long long int SCS_CONVERSION_MUL;
#endif  /* SCS_USE_FLT_MULT */







/************************************************************/
/* We have unrolled the loops for SCS_NB_WORDS==4 or 8 

   We just wish gcc would do it for us ! There are option switches,
   but they don't lead to any performance improvement. When they do,
   this part of the source code will be removed.

   In the meantime, feel free to unroll for other values. */


/***************************/
#if (SCS_NB_WORDS==4)
/***************************/
void scs_mul(scs_ptr result, scs_ptr x, scs_ptr y){
  SCS_CONVERSION_MUL     val, tmp;
  SCS_CONVERSION_MUL     r0,r1,r2,r3,r4;
  SCS_CONVERSION_MUL     x0,x1,x2,x3;
  int                    y0,y1,y2,y3;
    
  R_EXP = X_EXP * Y_EXP;
  R_SGN = X_SGN * Y_SGN;
  R_IND = X_IND + Y_IND;

  /* Partial products computation */   
  x3=X_HW[3];  y3=Y_HW[3];  x2=X_HW[2];  y2=Y_HW[2];
  x1=X_HW[1];  y1=Y_HW[1];  x0=X_HW[0];  y0=Y_HW[0];

  r4 = x3*y1 + x2*y2 + x1*y3;
  r3 = x3*y0 + x2*y1 + x1*y2 + x0*y3;
  r2 = x2*y0 + x1*y1 + x0*y2;
  r1 = x1*y0 + x0*y1 ;
  r0 = x0*y0;

  val= 0;
  /* Carry Propagate */
  SCS_CARRY_PROPAGATE(r4,r3,tmp)
  SCS_CARRY_PROPAGATE(r3,r2,tmp)
  SCS_CARRY_PROPAGATE(r2,r1,tmp)
  SCS_CARRY_PROPAGATE(r1,r0,tmp)      
  SCS_CARRY_PROPAGATE(r0,val,tmp)      
 
  if(val != 0){
    /* shift all the digits ! */
    R_HW[0] = val; R_HW[1] = r0; R_HW[2] = r1;  R_HW[3] = r2;
    R_IND += 1;
  }
  else {
    R_HW[0] = r0; R_HW[1] = r1; R_HW[2] = r2; R_HW[3] = r3;
  }

}


void scs_square(scs_ptr result, scs_ptr x){
  SCS_CONVERSION_MUL  r0,r1,r2,r3,r4;
  SCS_CONVERSION_MUL  x0,x1,x2,x3;
  SCS_CONVERSION_MUL  val, tmp;


  R_EXP = X_EXP * X_EXP;
  R_IND = X_IND + X_IND;
  R_SGN = 1;
    
  /*
   * Calcul des PP
   */   
  x3=X_HW[3];  x2=X_HW[2];  x1=X_HW[1];  x0=X_HW[0];

  r0 =  x0*x0;
  r1 = (x0*x1)* 2 ;
  r2 =  x1*x1 + (x0*x2*2);
  r3 = (x1*x2 +  x0*x3)* 2;
  r4 =  x2*x2 + (x1*x3)* 2;

  val= 0;
  /* Propagation des retenues */
  SCS_CARRY_PROPAGATE(r4,r3,tmp)
  SCS_CARRY_PROPAGATE(r3,r2,tmp)
  SCS_CARRY_PROPAGATE(r2,r1,tmp)
  SCS_CARRY_PROPAGATE(r1,r0,tmp)      
  SCS_CARRY_PROPAGATE(r0,val,tmp)      
 
  if(val != 0){
    /* shift all the digits ! */
    R_HW[0] = val; R_HW[1] = r0; R_HW[2] = r1;  R_HW[3] = r2;
    R_IND += 1;
  }
  else {
    R_HW[0] = r0; R_HW[1] = r1; R_HW[2] = r2; R_HW[3] = r3;
  }
  
}







/***************************/
#elif (SCS_NB_WORDS==8)
/***************************/
void scs_mul(scs_ptr result, scs_ptr x, scs_ptr y){
  SCS_CONVERSION_MUL     val, tmp;
  SCS_CONVERSION_MUL     r0,r1,r2,r3,r4,r5,r6,r7,r8;
  SCS_CONVERSION_MUL     x0,x1,x2,x3,x4,x5,x6,x7;
  int                    y0,y1,y2,y3,y4,y5,y6,y7;
    
  R_EXP = X_EXP * Y_EXP;
  R_SGN = X_SGN * Y_SGN;
  R_IND = X_IND + Y_IND;

  /* Partial products computation */   
  x7=X_HW[7];  y7=Y_HW[7];  x6=X_HW[6];  y6=Y_HW[6];
  x5=X_HW[5];  y5=Y_HW[5];  x4=X_HW[4];  y4=Y_HW[4];
  x3=X_HW[3];  y3=Y_HW[3];  x2=X_HW[2];  y2=Y_HW[2];
  x1=X_HW[1];  y1=Y_HW[1];  x0=X_HW[0];  y0=Y_HW[0];

  r8 = x7*y1 + x6*y2 + x5*y3 + x4*y4 + x3*y5 + x2*y6 + x1*y7;
  r7 = x7*y0 + x6*y1 + x5*y2 + x4*y3 + x3*y4 + x2*y5 + x1*y6 + x0*y7;
  r6 = x6*y0 + x5*y1 + x4*y2 + x3*y3 + x2*y4 + x1*y5 + x0*y6;
  r5 = x5*y0 + x4*y1 + x3*y2 + x2*y3 + x1*y4 + x0*y5;
  r4 = x4*y0 + x3*y1 + x2*y2 + x1*y3 + x0*y4 ;
  r3 = x3*y0 + x2*y1 + x1*y2 + x0*y3;
  r2 = x2*y0 + x1*y1 + x0*y2;
  r1 = x1*y0 + x0*y1 ;
  r0 = x0*y0 ;
 
  val= 0;
  /* Carry Propagate */
  SCS_CARRY_PROPAGATE(r8,r7,tmp)
  SCS_CARRY_PROPAGATE(r7,r6,tmp)
  SCS_CARRY_PROPAGATE(r6,r5,tmp)
  SCS_CARRY_PROPAGATE(r5,r4,tmp)
  SCS_CARRY_PROPAGATE(r4,r3,tmp)
  SCS_CARRY_PROPAGATE(r3,r2,tmp)
  SCS_CARRY_PROPAGATE(r2,r1,tmp)
  SCS_CARRY_PROPAGATE(r1,r0,tmp)      
  SCS_CARRY_PROPAGATE(r0,val,tmp)      
 
  if(val != 0){
    /* shift all the digits ! */
    R_HW[0] = val; R_HW[1] = r0; R_HW[2] = r1;  R_HW[3] = r2;
    R_HW[4] = r3;  R_HW[5] = r4; R_HW[6] = r5;  R_HW[7] = r6;
    R_IND += 1;
  }
  else {
    R_HW[0] = r0; R_HW[1] = r1; R_HW[2] = r2; R_HW[3] = r3;
    R_HW[4] = r4; R_HW[5] = r5; R_HW[6] = r6; R_HW[7] = r7;
  }

}


void scs_square(scs_ptr result, scs_ptr x){
  SCS_CONVERSION_MUL  r0,r1,r2,r3,r4,r5,r6,r7,r8;
  SCS_CONVERSION_MUL  x0,x1,x2,x3,x4,x5,x6,x7;
  SCS_CONVERSION_MUL  val, tmp;


  R_EXP = X_EXP * X_EXP;
  R_IND = X_IND + X_IND;
  R_SGN = 1;
    
  /*
   * Partial products
   */   
  x7=X_HW[7];  x6=X_HW[6];  x5=X_HW[5];  x4=X_HW[4];
  x3=X_HW[3];  x2=X_HW[2];  x1=X_HW[1];  x0=X_HW[0];

  r0 =  x0*x0;
  r1 = (x0*x1)* 2 ;
  r2 =  x1*x1 + (x0*x2*2);
  r3 = (x1*x2 +  x0*x3)* 2;
  r4 =  x2*x2 + (x1*x3 + x0*x4)* 2;
  r5 = (x2*x3 +  x1*x4 + x0*x5)* 2;
  r6 =  x3*x3 + (x2*x4 + x1*x5 + x0*x6)* 2;
  r7 = (x3*x4 +  x2*x5 + x1*x6 + x0*x7)* 2;
  r8 =  x4*x4 + (x3*x5 + x2*x6 + x1*x7)* 2;

  val= 0;
  /* Carry propagation */
  SCS_CARRY_PROPAGATE(r8,r7,tmp)
  SCS_CARRY_PROPAGATE(r7,r6,tmp)
  SCS_CARRY_PROPAGATE(r6,r5,tmp)
  SCS_CARRY_PROPAGATE(r5,r4,tmp)
  SCS_CARRY_PROPAGATE(r4,r3,tmp)
  SCS_CARRY_PROPAGATE(r3,r2,tmp)
  SCS_CARRY_PROPAGATE(r2,r1,tmp)
  SCS_CARRY_PROPAGATE(r1,r0,tmp)      
  SCS_CARRY_PROPAGATE(r0,val,tmp)      
 
  if(val != 0){
    /* shift all the digits ! */
    R_HW[0] = val; R_HW[1] = r0; R_HW[2] = r1;  R_HW[3] = r2;
    R_HW[4] = r3;  R_HW[5] = r4; R_HW[6] = r5;  R_HW[7] = r6;
    R_IND += 1;
  }
  else {
    R_HW[0] = r0; R_HW[1] = r1; R_HW[2] = r2; R_HW[3] = r3;
    R_HW[4] = r4; R_HW[5] = r5; R_HW[6] = r6; R_HW[7] = r7;
  }
  
}



/***************************/
#else
/***************************/
/* From there on, the normal, unrolled case */


void scs_mul(scs_ptr result, scs_ptr x, scs_ptr y){
  SCS_CONVERSION_MUL RES[SCS_NB_WORDS+1];
  SCS_CONVERSION_MUL val, tmp;
  int i, j;    

  R_EXP = X_EXP * Y_EXP;
  R_SGN = X_SGN * Y_SGN;
  R_IND = X_IND + Y_IND;

  for(i=0; i<=SCS_NB_WORDS; i++)
    RES[i]=0;

  /* Compute only the first half of the partial product. See the
     unrolled code for an example of what we compute */

#ifdef SCS_TYPECPU_X86
  /* This is the only place where there is assembly code to force 64-bit
     arithmetic. Someday gcc will catch up here, too.
  */
  {
    db_number t;
    /* i=0 */
    for(j=0; j<(SCS_NB_WORDS); j++) {
      __asm__ volatile("mull %3" 
		       : "=a" (t.i[LO_ENDIAN]), "=d" (t.i[HI_ENDIAN])
		       : "a" (X_HW[0]) , "g" (Y_HW[j]));
      RES[j] += t.l;
    }
    /* i = 1..SCS_NB_WORDS-1 */
    for(i=1 ; i<SCS_NB_WORDS; i++){
      for(j=0; j<(SCS_NB_WORDS-i); j++){
	__asm__ volatile("mull %3" 
			 : "=a" (t.i[LO_ENDIAN]), "=d" (t.i[HI_ENDIAN])
			 : "a" (X_HW[i]) , "g" (Y_HW[j]));
	RES[i+j] += t.l;
      }
      __asm__ volatile("mull %3" 
		       : "=a" (t.i[LO_ENDIAN]), "=d" (t.i[HI_ENDIAN])
		       : "a" (X_HW[i]) , "g" (Y_HW[j])); 
      /* here j==SCS_NB_WORDS-i */
      RES[SCS_NB_WORDS] += t.l;
    }
 }

#else /* other architectures */

 /* i=0 */
 tmp = X_HW[0];
 for(j=0; j<(SCS_NB_WORDS); j++)
   RES[j] += tmp * Y_HW[j];
 /* i = 1..SCS_NB_WORDS-1 */
 for(i=1 ; i<SCS_NB_WORDS; i++){
      tmp = X_HW[i];
      for(j=0; j<(SCS_NB_WORDS-i); j++)
	RES[i+j] += tmp * Y_HW[j];
      RES[SCS_NB_WORDS] += tmp * Y_HW[j]; /* here j==SCS_NB_WORDS-i */
  }
#endif/* SCS_TYPECPU_X86 */

  val = 0;

  /* Carry propagate */
  for(i=SCS_NB_WORDS; i>0; i--)
    SCS_CARRY_PROPAGATE(RES[i],RES[i-1],tmp)
  SCS_CARRY_PROPAGATE(RES[0],val,tmp)


  /* Store the result */
  if(val != 0){
    /* shift all the digits ! */     
    R_HW[0] = val;
    for(i=1; i<SCS_NB_WORDS; i++)
      R_HW[i] = RES[i-1];
  
    R_IND += 1;
  }else {
    for(i=0; i<SCS_NB_WORDS; i++)
      R_HW[i] = RES[i];
   }
}







void scs_square(scs_ptr result, scs_ptr x){
  SCS_CONVERSION_MUL RES[SCS_NB_WORDS+1];
  SCS_CONVERSION_MUL val, tmp;
  int i, j;
  

  R_EXP = X_EXP * X_EXP;
  R_SGN = 1;
  R_IND = X_IND + X_IND;

  /* Set to 0 intermediate register     */
  for(i=0; i<=SCS_NB_WORDS; i++)
    RES[i] = 0;

  /* Compute all the double partial products: 2 x_i * x_j, i!=j */
  tmp = (SCS_CONVERSION_MUL)X_HW[0];
  for(j=1; j<SCS_NB_WORDS; j++)
    RES[j] += tmp * X_HW[j];
  for(i=1 ; i<(SCS_NB_WORDS+1)/2; i++){
    tmp = (SCS_CONVERSION_MUL)X_HW[i];
    for(j=i+1; j<(SCS_NB_WORDS-i); j++)
      RES[i+j] += tmp * X_HW[j];
    RES[SCS_NB_WORDS] += tmp * X_HW[SCS_NB_WORDS-i];
  }

  /* All these partial products are double */
  for(i=0; i<=SCS_NB_WORDS; i++)
    RES[i] *=2;

  /* Add partial product of the form x_i^2 */
  for(i=0, j=0; i<=SCS_NB_WORDS; i+=2, j++){
    RES[i]  += (SCS_CONVERSION_MUL)X_HW[j] * X_HW[j];
  }  

  val = 0;
  /* Carry propagate */
  for(i=SCS_NB_WORDS; i>0; i--)
      SCS_CARRY_PROPAGATE(RES[i],RES[i-1],tmp)
 
  SCS_CARRY_PROPAGATE(RES[0],val,tmp)


  /* Store the result */
  if(val != 0){
    /* shift all the digits ! */     
    R_HW[0] = val;
    for(i=1; i<SCS_NB_WORDS; i++)
      R_HW[i] = RES[i-1];
  
    R_IND += 1;
  }else {
    for(i=0; i<SCS_NB_WORDS; i++)
      R_HW[i] = RES[i];
   }

}


/* 
 * #endif corresponding to the test #if (SCS_NB_WORDS==8)
 */
#endif


/*
 Multiply x by an integer val; result is returned in x.
 */
 void scs_mul_ui(scs_ptr x, unsigned int val_int){
  SCS_CONVERSION_MUL val, tmp, vald, rr;
  int i;

  if (val_int == 0)
    X_EXP = 0;
  
  vald = val_int;

  val = 0; 
  rr  = 0;
  for(i=(SCS_NB_WORDS-1); i>=0; i--){
    val    += vald * X_HW[i];
    SCS_CARRY_PROPAGATE(val, rr, tmp)
    X_HW[i] = val;
    val     = rr;
    rr      = 0; 
  }

  if(val != 0){
    /* shift all the digits ! */ 
    for(i=(SCS_NB_WORDS-1); i>0; i--)
      X_HW[i] = X_HW[i-1];

    X_HW[0] = val;
    X_IND  += 1;
  }
  
  return;
}

