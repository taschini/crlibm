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

#ifdef HAVE_LIBMCR_H
#include <libmcr.h>
#endif


/* A variable equal to zero, stored here so that the compiler doesn't
   know its value in the other functions, which allows to prevent some
   optimizations  */


#define RN 1
#define RU 2
#define RD 3
#define RZ 4  

double zero ;

/* Here come the various random number generators. They all use the
   rand() function.  

   We may have two rand functions for each function under
   test. The first is for the soaktest, the second for testing the
   performance under what is supposed the main domain of use the function. 

   Typical examples:

   log has identical functions for soaktest and performance: random
   positive numbers. This means that negative numbers are not tested
   by soaktest, though.

   sin soaktests on all the floats, but tests for perf on a small
   interval around zero, shamelessely chosen as the one on which crlibm
   is the fastest.
*/


/**/


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
  
  result.i[LO]=rand_int();
  result.i[HI]=rand_int();
 
  return result.d;
}


/* Return a random double between 0 and 1, with a normal law on the
   exponent  */

double rand_double(){
  db_number result;
  int e;
  /*first the low bits of the mantissa*/
  result.i[LO]=rand_int();
  /* then the high bits of the mantissa, and the sign bit */
  result.i[HI]=  rand_int() & 0x000fffff;
  /* Now set the exponent (negative value) */
  e = rand() & 0x000003ff; 
  if (e>0) e-=1;
  result.i[HI] += e<<20;
  return (result.d);
}


/* Return a random double between 1 and 2, with a normal law on the
   mantissa and a constant exponent  */

double rand_double_normal(){
  db_number result;
  int e;
  /*first the low bits of the mantissa*/
  result.i[LO]=rand_int();
  /* then the high bits of the mantissa, and the sign bit */
  result.i[HI]=  rand_int() & 0x000fffff;
  /* Now set the exponent */
  e = 1023; 
  result.i[HI] += e<<20;
  return (result.d);
}



/* For exp we will test perf on numbers with a random sign, a random mantissa, and
   a random exponent between -9 and 9. And we soaktest on all the doubles */

#define rand_for_exp_soaktest rand_generic

double rand_for_exp_perf(){
  db_number result;
  int e;

  /*first the low bits of the mantissa*/
  result.i[LO]=rand_int();
  /* then the high bits of the mantissa, and the sign bit */
  result.i[HI]=  rand_int() & 0x800fffff;
  /* Now set the exponent between -9 and 9, enough to cover the useful range  */
  e =  (int) ( (rand_double_normal()-1) * 18 );
  result.i[HI] += (1023 + e -9)<<20;
  return result.d;
}


/* a number in the range which never produces over/underflow for the
   exp function. I don't trust the randomness of this function */
double rand_for_exp_normal(){
  return((750+710)*(rand_double_normal()-1)-750);
}





#define rand_for_csh_perf rand_for_exp_perf
/* I wish we could soaktest using rand_generic, but current MPFR is
   very slow for small and large arguments (up to a few minutes per
   call !). To check regularly, this is bound to improve. */

#define rand_for_csh_soaktest  rand_for_exp_perf

/* For log we only test the positive numbers*/
double rand_for_log(){
  db_number result;
  int e;

  /*first the low bits of the mantissa*/
  result.i[LO]=rand_int();
  /* then the high bits of the mantissa, and the sign bit */
  result.i[HI]=  rand_int() & 0x7fffffff;
  /* printf("x = %1.5e\n", result.d);*/
  return result.d;
}


/* For trigonometric functions it is difficult to tell what the test function should be */ 

double rand_for_trig_perf(){
  db_number result;
  int e;
  /*first the low bits of the mantissa*/
  result.i[LO]=rand_int();
  /* then the high bits of the mantissa, and the sign bit */
  result.i[HI]=  rand_int() & 0x800fffff;
  /* Now set the exponent  between -20 and 40 */
    e =  (int) ( (rand_double_normal()-1) * 60 ); 
   result.i[HI] += (1023 + e -20)<<20;
  return result.d;

}

#if 0
#define rand_for_trig_soaktest rand_generic
#else
#define rand_for_trig_soaktest rand_for_trig_perf
#endif

double rand_for_atan_perf(){
  db_number result;
  int e;
  /*first the low bits of the mantissa*/
  result.i[LO]=rand_int();
  /* then the high bits of the mantissa, and the sign bit */
  result.i[HI]=  rand_int() & 0x800fffff;
  /* Now set the exponent between -20 and 50, enough to cover the useful range  */
  e =  (int) ( (rand_double_normal()-1) * 70 );
  result.i[HI] += (1023 + e -20)<<20;
  return result.d;

}

double rand_for_atan_soaktest(){
  db_number result;
  int e;

  /*first the low bits of the mantissa*/
  result.i[LO]=rand_int();
  /* then the high bits of the mantissa, and the sign bit */
  result.i[HI]=  rand_int() & 0x800fffff;
  /* Now set the exponent between -20 and 50, enough to cover the useful range  */
  e =  (int) ( (rand_double_normal()-1) * 70 );
  result.i[HI] += (1023 + e -20)<<20;
  return result.d;

}

/* For pow we need to test the whole range of floating point numbers
 * However these definition should be slightly modified (keep x^y fp).
 */
#define rand_for_pow_soaktest rand_for_exp_perf

double rand_for_pow_perf(){
  db_number result;
  int e;

  /*first the low bits of the mantissa*/
  result.i[LO]=rand_int();
  /* then the high bits of the mantissa, and the sign bit */
  result.i[HI]=  rand_int() & 0x800fffff;
  /* Now set the exponent between -7 and 7, 
     enough to cover the useful range (does not overflow) */
  e =  (int) ( (rand_double_normal()-1) * 14 );
  result.i[HI] += (1023 + e -7)<<20;
  return result.d;
}




void test_rand()  {
  int i;
  double min=1e300, max=0.0;
  db_number input;
  for(i=0; i< 1000; i++){ 
    input.d = rand_for_exp_perf();
    if (input.d<min) min=input.d;
    if (input.d>max) max=input.d;
    printf("%1.5ex \t%.8x %.8x\t%1.5e \t%1.5e\n", input.d, input.i[HI], input.i[LO],min,max );
  }
}





/* general init function */


void test_init(/* pointers to returned value */
	       double (**randfun_perf)(), 
	       double (**randfun_soaktest)(), 
	       double (**testfun_crlibm)(), 
	       int    (**testfun_mpfr)  (),
	       double (**testfun_libultim)   (),
	       double (**testfun_libmcr)  (),
	       double (**testfun_libm)  (),
	       double* worst_case,
	       /* arguments */
	       char *func_name,
	       char *rnd_mode)  {

  int crlibm_rnd_mode;

  /* We have added the rounding mode designation used in libmcr's test files */
  if      ((strcmp(rnd_mode,"RU")==0) || (strcmp(rnd_mode,"P")==0)) crlibm_rnd_mode = RU;
  else if ((strcmp(rnd_mode,"RD")==0) || (strcmp(rnd_mode,"M")==0)) crlibm_rnd_mode = RD;
  else if ((strcmp(rnd_mode,"RZ")==0) || (strcmp(rnd_mode,"Z")==0)) crlibm_rnd_mode = RZ;
  else if ((strcmp(rnd_mode,"RN")==0) || (strcmp(rnd_mode,"N")==0)) crlibm_rnd_mode = RN;
  else {
    fprintf(stderr, "Unknown rounding mode: %s, exiting\n", rnd_mode);
    exit(EXIT_FAILURE);
  }


  *randfun_perf     = rand_generic; /* the default random function */
  *randfun_soaktest = rand_generic; /* the default random function */
  *testfun_mpfr     = NULL;
  *testfun_libm     = NULL;
  *worst_case=0.;

  if (strcmp (func_name, "exp") == 0)
    {
      *randfun_perf     = rand_for_exp_perf;
      *randfun_soaktest = rand_for_exp_soaktest;
      *worst_case= .75417527749959590085206221024712557043923055744016892276704311370849609375e-9;
      *testfun_libm   = exp;
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = exp_ru;	break;
      case RD:
	*testfun_crlibm = exp_rd;	break;
      case RZ:
	*testfun_crlibm = exp_rz;	break;
      default:
	*testfun_crlibm = exp_rn;
      }
#ifdef HAVE_MATHLIB_H
      *testfun_libultim    = uexp;
#endif
#ifdef HAVE_LIBMCR_H
      *testfun_libmcr    = __libmcr_exp;
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_exp;
#endif
    }


  else  if (strcmp (func_name, "log") == 0)
    {
      *randfun_perf     = rand_for_log;
      *randfun_soaktest = rand_for_log;
      *worst_case=0.4009793462309855760053830468258630076242931610568335144339734234840014178511334897967240437927437320e-115;
      *testfun_libm   = log;
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = log_ru;	break;
      case RD:
	*testfun_crlibm = log_rd;	break;
      case RZ:
	*testfun_crlibm = log_rz;	break;
      default:
	*testfun_crlibm = log_rn;
      }
#ifdef HAVE_MATHLIB_H
      *testfun_libultim    = ulog;
#endif
#ifdef HAVE_LIBMCR_H
      *testfun_libmcr    = __libmcr_log;
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_log;
#endif
    }


  else  if (strcmp (func_name, "sin") == 0)
    {
      *randfun_perf     = rand_for_trig_perf;
      *randfun_soaktest = rand_for_trig_soaktest;
      *worst_case= 9.24898516520941904595076721307123079895973205566406e-01;
      *testfun_libm   = sin;
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = sin_ru;	break;
      case RD:
	*testfun_crlibm = sin_rd;	break;
      case RZ:
	*testfun_crlibm = sin_rz;	break;
      default:
	*testfun_crlibm = sin_rn;
      }
#ifdef HAVE_MATHLIB_H
      *testfun_libultim    = usin;
#endif
#ifdef HAVE_LIBMCR_H
      *testfun_libmcr    = __libmcr_sin;
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_sin;
#endif
    }

  else  if (strcmp (func_name, "cos") == 0)
    {
      *randfun_perf     = rand_for_trig_perf;
      *randfun_soaktest = rand_for_trig_soaktest;
      *worst_case=  8.87406081479789610177988379291491582989692687988281e-01;
      *testfun_libm   = cos;
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = cos_ru;	break;
      case RD:
	*testfun_crlibm = cos_rd;	break;
      case RZ:
	*testfun_crlibm = cos_rz;	break;
      default:
	*testfun_crlibm = cos_rn;
      }
#ifdef HAVE_MATHLIB_H
      *testfun_libultim    = ucos;
#endif
#ifdef HAVE_LIBMCR_H
      *testfun_libmcr    = __libmcr_cos;
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_cos;
#endif
    }

  else  if (strcmp (func_name, "tan") == 0)
    {
      *randfun_perf     = rand_for_trig_perf;
      *randfun_soaktest = rand_for_trig_soaktest;
      *worst_case= 1.18008664944477814628953638020902872085571289062500e-01;
      *testfun_libm   = tan; 
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = tan_ru;	break;
      case RD:
	*testfun_crlibm = tan_rd;	break;
      case RZ:
	*testfun_crlibm = tan_rz;	break;
      default:
	*testfun_crlibm = tan_rn;
      }
#ifdef HAVE_MATHLIB_H
      *testfun_libultim    = utan;
#endif
#ifdef HAVE_LIBMCR_H
      *testfun_libmcr    = __libmcr_tan;
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_tan;
#endif
    }

#if 0 /* No cotan in the standard math.h ? */
  else  if (strcmp (func_name, "cotan") == 0)
    {
      *randfun_perf     = rand_for_trig_perf;
      *randfun_soaktest = rand_for_trig_soaktest;
      *worst_case= 1.18008664944477814628953638020902872085571289062500e-01;
      *testfun_libm   = cotan; 
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = cotan_ru;	break;
      case RD:
	*testfun_crlibm = cotan_rd;	break;
      case RZ:
	*testfun_crlibm = cotan_rz;	break;
      default:
	*testfun_crlibm = cotan_rn;
      }
#ifdef HAVE_MATHLIB_H
      *testfun_libultim    = ucotan;
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_cotan;
#endif
    }
#endif /* no cotan */

  else if (strcmp (func_name, "atan") == 0)
    {
      *randfun_perf     = rand_for_atan_perf;
      *randfun_soaktest = rand_for_atan_soaktest;
      *worst_case= 9.54714164331460501955461950274184346199035644531250e-02; 
      *testfun_libm   = atan;
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = atan_ru;	break;
      case RD:
	*testfun_crlibm = atan_rd;	break;
      case RZ:
	*testfun_crlibm = atan_rz;	break;
      default:
        *testfun_crlibm = atan_rn ;
      }
#ifdef HAVE_MATHLIB_H
      *testfun_libultim    = uatan;
#endif
#ifdef HAVE_LIBMCR_H
      *testfun_libmcr    = __libmcr_atan;
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_atan;
#endif
    }



  else if (strcmp (func_name, "cosh") == 0)
    {
      *randfun_perf     = rand_for_csh_perf;
      *randfun_soaktest = rand_for_csh_soaktest;
      *worst_case= 3.76323248339103422210882854415103793144226074218750e+00;
      *testfun_libm   = cosh;
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = cosh_ru;	break;
      case RD:
	*testfun_crlibm = cosh_rd;	break;
      case RZ:
	*testfun_crlibm = cosh_rz;	break;
      default:
	*testfun_crlibm = cosh_rn;
      }
#ifdef HAVE_MATHLIB_H
      /* No hyperbolic function in Ziv library */ 
      *testfun_libultim    = NULL;
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_cosh;
#endif
    }

  else if (strcmp (func_name, "sinh") == 0)
    {
      *randfun_perf     = rand_for_csh_perf;
      *randfun_soaktest = rand_for_csh_soaktest;
      *worst_case= 5.81191276791475441854117889306508004665374755859375;
      *testfun_libm   = sinh;
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = sinh_ru;	break;
      case RD:
	*testfun_crlibm = sinh_rd;	break;
      case RZ:
	*testfun_crlibm = sinh_rz; 	break;
      default:
	*testfun_crlibm = sinh_rn;
      }
#ifdef HAVE_MATHLIB_H
      /* No hyperbolic function in Ziv library */ 
      *testfun_libultim    = NULL;
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr   = mpfr_sinh;
#endif
    }

  else if (strcmp (func_name, "pow") == 0)
    {
      *randfun_perf     = rand_for_pow_perf;
      *randfun_soaktest = rand_for_pow_soaktest;
      *worst_case= 1;
      *testfun_libm   = pow;
      switch(crlibm_rnd_mode){
      case RU:
	*testfun_crlibm = NULL;	break;
      case RD:
	*testfun_crlibm = NULL;	break;
      case RZ:
	*testfun_crlibm = NULL;	break;
      default:
	/*	*testfun_crlibm = pow_rn;*/
	*testfun_crlibm = NULL;

      }
#ifdef HAVE_MATHLIB_H
      *testfun_libultim = upow;
#endif
#ifdef HAVE_LIBMCR_H
      *testfun_libmcr    = __libmcr_pow;
#endif
#ifdef HAVE_MPFR_H
      *testfun_mpfr     = mpfr_pow; 
#endif
    }



  else
    {
      fprintf (stderr, "Unknown function: %s\n", func_name);
      exit (1);
    }
}

