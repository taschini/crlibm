#include "crlibm.h"
#include "crlibm_private.h"
#include "test_common.h"

/* A variable equal to zero, stored here so that the compiler doesn't
   know its value in the other functions, which allows to prevent some
   optimizations  */

double zero ;

/* Here come the various random number generators. They all use the
   rand() function. We probably don't mind setting the seed (srand())
   for reproducibility, since the program outputs the error cases */


/* Return 'sizeof(int)' random bits    */
int rand_int(){
  int val;
  int i;
  val = (random() & 0x000000ff);
  for(i=0; i<(sizeof(int)); i++){
    val = val << 8;
    val += (random() & 0x000000ff ); /* we keep only 8 bits */
  }
  return val;
}




/* Return a completely random double  */

double rand_generic(){
  db_number result;
  
  result.i[LO_ENDIAN]=rand_int();
  result.i[HI_ENDIAN]=rand_int();
 
  return result.d;
}


/* Return a random double between 0 and 1, with a normal law on the
   exponent  */

double rand_double(){
  db_number result;
  int e;
  /*first the low bits of the mantissa*/
  result.i[LO_ENDIAN]=rand_int();
  /* then the high bits of the mantissa, and the sign bit */
  result.i[HI_ENDIAN]=  rand_int() & 0x000fffff;
  /* Now set the exponent (negative value) */
  e = rand() & 0x000003ff; 
  if (e>0) e-=1;
  result.i[HI_ENDIAN] += e<<20;
  return (result.d);
}


/* Return a random double between 1 and 2, with a normal law on the
   mantissa and a constant exponent  */

double rand_double_normal(){
  db_number result;
  int e;
  /*first the low bits of the mantissa*/
  result.i[LO_ENDIAN]=rand_int();
  /* then the high bits of the mantissa, and the sign bit */
  result.i[HI_ENDIAN]=  rand_int() & 0x000fffff;
  /* Now set the exponent */
  e = 1023; 
  result.i[HI_ENDIAN] += e<<20;
  return (result.d);
}


/* For exp we will soaktest with a random sign, a random mantissa, and
   a random exponent between -16 and 15 */
double rand_for_exp(){
  db_number result;
  int e;

  /*first the low bits of the mantissa*/
  result.i[LO_ENDIAN]=rand_int();
  /* then the high bits of the mantissa, and the sign bit */
  result.i[HI_ENDIAN]=  rand_int() & 0x800fffff;
  /* Now set the exponent between -10 and 10, enough to cover the useful range  */
  e =  (int) ( (rand_double_normal()-1) * 20 );
  result.i[HI_ENDIAN] += (1023 + e -10)<<20;
  return result.d;
}

/* For log we only test the positive numbers*/
double rand_for_log(){
  db_number result;
  int e;

  /*first the low bits of the mantissa*/
  result.i[LO_ENDIAN]=rand_int();
  /* then the high bits of the mantissa, and the sign bit */
  result.i[HI_ENDIAN]=  rand_int() & 0x7fffffff;
  /* printf("x = %1.5e\n", result.d);*/
  return result.d;
}


/* a number in the range which never produces over/underflow for the
   exp function */
double rand_for_exp_normal(){
  return((750+710)*(rand_double_normal()-1)-750);
}




void test_rand()  {
  int i;
  double min=1e300, max=0.0;
  db_number input;
  for(i=0; i< 1000; i++){ 
    input.d = rand_for_exp();
    if (input.d<min) min=input.d;
    if (input.d>max) max=input.d;
    printf("%1.5ex \t%.8x %.8x\t%1.5e \t%1.5e\n", input.d, input.i[HI_ENDIAN], input.i[LO_ENDIAN],min,max );
  }
}


