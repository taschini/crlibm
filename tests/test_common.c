#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "crlibm.h"
#include "crlibm_private.h"
#include "test_common.h"

#ifdef HAVE_MPFR_H
#include <gmp.h>
#include <mpfr.h>
#endif

#ifdef HAVE_MATHLIB_H
#include <MathLib.h>
#endif


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

#define rand_for_cosh rand_for_exp
#define rand_for_sinh rand_for_exp


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


/* For exp we will soaktest with a random sign, a random mantissa, and
   a random exponent between -5 and 15 */
double rand_for_sin(){
  db_number result;
  int e;

  /*first the low bits of the mantissa*/
  result.i[LO_ENDIAN]=rand_int();
  /* then the high bits of the mantissa, and the sign bit */
  result.i[HI_ENDIAN]=  rand_int() & 0x800fffff;
  /* Now set the exponent between -10 and 15, enough to cover the useful range  */
#if 0  
  e =  (int) ( (rand_double_normal()-1) * 25 );
  result.i[HI_ENDIAN] += (1023 + e -10)<<20;
#else
  e =  (int) ( (rand_double_normal()-1) * 40 ); /* never Payne Hanek : 34 */
  result.i[HI_ENDIAN] += (1023 + e -10)<<20;
#endif
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





/* general init function */


void test_init(/* pointers to returned value */
	       double (**randfun)(), 
	       double (**testfun_crlibm)(), 
	       int    (**testfun_mpfr)  (),
	       double (**testfun_libm)  (),
	       double (**testfun_ibm)   (),
	       double* worst_case,
	       /* arguments */
	       char *func_name,
	       char *rnd_mode)  {

  int crlibm_rnd_mode;
  
  if      (strcmp(rnd_mode,"RU")==0) crlibm_rnd_mode = 2;
  else if (strcmp(rnd_mode,"RD")==0) crlibm_rnd_mode = 3;
  else if (strcmp(rnd_mode,"RZ")==0) crlibm_rnd_mode = 4;
  else crlibm_rnd_mode = 1;


  *randfun        = rand_generic; /* the default random function */
  *testfun_mpfr   = NULL;
  *testfun_libm   = NULL;
  *worst_case=0.;

  if (strcmp (func_name, "exp") == 0)
    {
      *randfun        = rand_for_exp;
      *worst_case= .75417527749959590085206221024712557043923055744016892276704311370849609375e-9;
      *testfun_libm   = exp;
      switch(crlibm_rnd_mode){
      case 2:
	*testfun_crlibm = exp_ru;	break;
      case 3:
	*testfun_crlibm = exp_rd;	break;
      case 4:
	*testfun_crlibm = exp_rz;	break;
      default:
	*testfun_crlibm = exp_rn;
      }
#ifdef HAVE_MATHLIB_H
      *testfun_ibm    = uexp;
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_exp;
#endif
    }
  else if (strcmp (func_name, "atan") == 0)
    {
      *randfun        = rand_generic;
      *worst_case= .75417527749959590085206221024712557043923055744016892276704311370849609375e-9;
      *testfun_libm   = atan;
      switch(crlibm_rnd_mode){
      case 2:
	*testfun_crlibm = atan_rn;	break;
      case 3:
	*testfun_crlibm = atan_rn;	break;
      case 4:
	*testfun_crlibm = atan_rn;	break;
      default:
        *testfun_crlibm = atan_rn ;
      }
#ifdef HAVE_MATHLIB_H
      *testfun_ibm    = uatan;
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_atan;
#endif
    }

  else if (strcmp (func_name, "cosh") == 0)
    {
      *randfun        = rand_for_cosh;
      *worst_case= .75417527749959590085206221024712557043923055744016892276704311370849609375e-9;
      *testfun_libm   = cosh;
      switch(crlibm_rnd_mode){
      case 2:
	*testfun_crlibm = cosh_ru;	break;
      case 3:
	*testfun_crlibm = cosh_rd;	break;
      case 4:
	*testfun_crlibm = cosh_rz;	break;
      default:
	*testfun_crlibm = cosh_rn;
      }
#ifdef HAVE_MATHLIB_H
      /* No hyperbolic function in Ziv library */ 
      *testfun_ibm    = NULL;
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_cosh;
#endif
    }

  else if (strcmp (func_name, "sinh") == 0)
    {
      *randfun        = rand_for_sinh;
      *worst_case= .75417527749959590085206221024712557043923055744016892276704311370849609375e-9;
      *testfun_libm   = sinh;
      switch(crlibm_rnd_mode){
      case 2:
	*testfun_crlibm = sinh_ru;	break;
      case 3:
	*testfun_crlibm = sinh_rd;	break;
      case 4:
	*testfun_crlibm = sinh_rz;	break;
      default:
	*testfun_crlibm = sinh_rn;
      }
#ifdef HAVE_MATHLIB_H
      /* No hyperbolic function in Ziv library */ 
      *testfun_ibm    = NULL;
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_sinh;
#endif
    }

  else  if (strcmp (func_name, "log") == 0)
    {
      *randfun        = rand_for_log;
      *worst_case=0.4009793462309855760053830468258630076242931610568335144339734234840014178511334897967240437927437320e-115;
      *testfun_libm   = log;
      switch(crlibm_rnd_mode){
      case 2:
	*testfun_crlibm = log_ru;	break;
      case 3:
	*testfun_crlibm = log_rd;	break;
      case 4:
	*testfun_crlibm = log_rz;	break;
      default:
	*testfun_crlibm = log_rn;
      }
#ifdef HAVE_MATHLIB_H
      *testfun_ibm    = ulog;
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_log;
#endif
    }

  else  if (strcmp (func_name, "sin") == 0)
    {
      *randfun        = rand_for_sin;
      /* This worst case works fast for Ziv. Damn. */
      *worst_case=0.498498785880875427967140467444551177322864532470703125000000 ;
      *testfun_libm   = sin;
      switch(crlibm_rnd_mode){
      case 2:
	*testfun_crlibm = sin_rn;	break;
      case 3:
	*testfun_crlibm = sin_rn;	break;
      case 4:
	*testfun_crlibm = sin_rn;	break;
      default:
	*testfun_crlibm = sin_rn;
      }
#ifdef HAVE_MATHLIB_H
      *testfun_ibm    = usin;
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_sin;
#endif
    }

  else  if (strcmp (func_name, "tan") == 0)
    {
      *randfun        = rand_for_sin; 
      /* *worst_case=0.4009793462309855760053830468258630076242931610568335144339734234840014178511334897967240437927437320e-115;
      */ 
      *testfun_libm   = tan; 
      switch(crlibm_rnd_mode){
      case 2:
	*testfun_crlibm = tan_rn;	break;
      case 3:
	*testfun_crlibm = tan_rn;	break;
      case 4:
	*testfun_crlibm = tan_rn;	break;
      default:
	*testfun_crlibm = tan_rn;
      }
#ifdef HAVE_MATHLIB_H
      *testfun_ibm    = utan;
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_tan;
#endif
    }
  else  if (strcmp (func_name, "scs_tan") == 0)
    {
      *randfun        = rand_for_sin; 
      /* *worst_case=0.4009793462309855760053830468258630076242931610568335144339734234840014178511334897967240437927437320e-115;
      */ 
      *testfun_libm   = tan; 
      switch(crlibm_rnd_mode){
      case 2:
	*testfun_crlibm = scs_tan_rn;	break;
      case 3:
	*testfun_crlibm = scs_tan_rn;	break;
      case 4:
	*testfun_crlibm = scs_tan_rn;	break;
      default:
	*testfun_crlibm = scs_tan_rn;
      }
#ifdef HAVE_MATHLIB_H
      *testfun_ibm    = utan;
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_tan;
#endif
    }

  else
    {
      fprintf (stderr, "Unknown function: %s\n", func_name);
      exit (1);
    }
}

