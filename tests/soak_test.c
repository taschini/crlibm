#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "crlibm.h"
#include "crlibm_private.h"
#ifdef HAVE_MPFR_H  /* stop here if MPFR not present */
#include "test_common.h"
#include <gmp.h>
#include <mpfr.h>

#ifdef HAVE_MATHLIB_H
#include <MathLib.h>
#endif

/* Stupidely soak-tests a function against mpfr */

/* if set to 1, print out detailed errors (vith values of x and exp(x))
   if set to 0, only count the errors and print out the count */
#define DETAILED_REPORT 1

/* if set to 1, print out errors due to NaNs
   if set to 0, don't print NaNs' errors */
#define PRINT_NAN 0
/* What we are looking for here is misrounding. Therefore these tests
   concern only intervals on which the function does something
   useful. For example for exp, it is enough to soaktest on -1024,1024.

   To achieve this we have several random generator function, tailored
   for each function or group of function.

   The other cases (including all the special cases) are supposedly
   tested exhaustively by the other programs of this directory. */ 

/* Basic-like programming with global variables: */
db_number input, res_crlibm, res_mpfr, res_ibm, res_libm;
mpfr_t mp_res, mp_inpt; 

/* The random number generator*/
double (*randfun)       () = NULL;
/* The function we test */
double (*testfun_crlibm)() = NULL;
/* The function we trust */
int    (*testfun_mpfr)  () = NULL;
/* The function to show off against for accuracy  */
double (*testfun_libm)  () = NULL;
/* The function to show off against for performance */
double (*testfun_ibm)   () = NULL;


/*
 * Define the rounding mode to test
 * 1 : to nearest              * 2 : to + inf
 * 3 : to - inf                * 4 : to 0
 */
int scs_rnd_mode=1;

/*
 * Rounding mode to test
 */
mp_rnd_t rnd;

void init(char *function){
  mpfr_init2(mp_res,  153);
  mpfr_init2(mp_inpt, 53);

  randfun        = rand_generic; /* the default random function */

  if (strcmp (function, "exp") == 0)
    {
      randfun        = rand_for_exp;
      testfun_libm   = exp;
      switch(scs_rnd_mode){
      case 2:
	testfun_crlibm = exp_ru;	break;
      case 3:
	testfun_crlibm = exp_rd;	break;
      case 4:
	testfun_crlibm = exp_rz;	break;
      default:
	testfun_crlibm = exp_rn;
      }

#ifdef HAVE_MATHLIB_H
      testfun_ibm    = uexp;
#endif
      testfun_mpfr   = mpfr_exp;
    }

  else if (strcmp (function, "log") == 0)
    {
      testfun_libm   = log;
      switch(scs_rnd_mode){
      case 2:
	testfun_crlibm = log_ru;	break;
      case 3:
	testfun_crlibm = log_rd;	break;
      default:
	testfun_crlibm = log_rn;
      }
#ifdef HAVE_MATHLIB_H
      testfun_ibm    = ulog;
#endif
      testfun_mpfr   = mpfr_log;
    }

  else if (strcmp (function, "log2") == 0)
    {
      testfun_libm   = NULL /* doesn't exist ? */;
      testfun_crlibm = log2_rn;
#ifdef HAVE_MATHLIB_H
      testfun_ibm    = ulog2;
#endif
      testfun_mpfr   = mpfr_log2;
    }

  else if (strcmp (function, "log10") == 0)
    {
      testfun_libm   = log10;
      testfun_crlibm = log10_rn;
#ifdef HAVE_MATHLIB_H
      testfun_ibm    = NULL;
#endif
      testfun_mpfr   = mpfr_log10;
    }

  else if (strcmp (function, "sin") == 0)
    {
      testfun_libm   = sin;
      testfun_crlibm = sin_rn;
#ifdef HAVE_MATHLIB_H
      testfun_ibm    = usin;
#endif
      testfun_mpfr   = mpfr_sin;
    }

  else if (strcmp (function, "tan") == 0)
    {
      testfun_libm   = tan   ;
      testfun_crlibm = tan_rn;
#ifdef HAVE_MATHLIB_H
      testfun_ibm    = utan;
#endif
      testfun_mpfr   = mpfr_tan;
    }

  else if (strcmp (function, "cotan") == 0)
    {
      testfun_libm   = NULL;
      testfun_crlibm = cotan_rn;
#ifdef HAVE_MATHLIB_H
/*       testfun_ibm    = ucotan; doesn't exist ?*/
#endif
      testfun_mpfr   = NULL; /* should check */
    }

    
  else if (strcmp (function, "log") == 0)
    {
      testfun_crlibm = log_rn;
      testfun_mpfr   = mpfr_log;
      randfun        = rand_generic;
    }
  else if (strcmp (function, "log2") == 0)
    {
      testfun_crlibm = log2_rn;
      testfun_mpfr   = mpfr_log2;
      randfun        = rand_generic;
    }
  else if (strcmp (function, "log10") == 0)
    {
      testfun_crlibm = log10_rn;
      testfun_mpfr   = mpfr_log10;
      randfun        = rand_generic;
    }
  else if (strcmp (function, "sin") == 0)
    {
      testfun_crlibm = sin_rn;
      testfun_mpfr   = mpfr_sin;
      randfun        = rand_generic;
    }
  else if (strcmp (function, "cos") == 0)
    {
      testfun_crlibm = cos_rn;
      testfun_mpfr   = mpfr_cos;
      randfun        = rand_generic;
    }
  else if (strcmp (function, "tan") == 0)
    {
      testfun_crlibm = tan_rn;
      testfun_mpfr   = mpfr_tan;
    }
#if 0 /* mpfr has no cotan? */ 
  else if (strcmp (function, "cotan") == 0)
    {
      testfun_crlibm = cotan_rn;
      testfun_mpfr   = mpfr_cotan;
    }
#endif
  else
    {
      fprintf (stderr, "Unknown function: %s\n", function);
      exit (1);
    }
}




void test() {
  int counter=0;
  int failures=0;

  while(1+1==2) {
    counter++;

    input.d = randfun();

    res_crlibm.d = testfun_crlibm(input.d);

    mpfr_set_d(mp_inpt, input.d, GMP_RNDN);
    testfun_mpfr(mp_res, mp_inpt, rnd);
    res_mpfr.d = mpfr_get_d(mp_res, rnd);
    
    if( (res_crlibm.i[LO_ENDIAN] != res_mpfr.i[LO_ENDIAN]) 
	|| (res_crlibm.i[HI_ENDIAN] != res_mpfr.i[HI_ENDIAN]) )  {
      printf("*** E  x=%.70e \n            (%8x %8x) \n", 
	     input.d, 
	     input.i[HI_ENDIAN], 
	     input.i[LO_ENDIAN]);
      
      printf("We find    %.50e \n         (%8x %8x) \n", 
	     res_crlibm.d, 
	     res_crlibm.i[HI_ENDIAN], 
	     res_crlibm.i[LO_ENDIAN]);
      printf("MPFR gives %.50e \n         (%8x %8x) \n\n", 
	     res_mpfr.d, 
	     res_mpfr.i[HI_ENDIAN], 
	     res_mpfr.i[LO_ENDIAN]);
      failures++;
      printf("That's %d failures out of %d (ratio %e)\n",
	     failures, counter, 
	     ((double)failures)/(double)counter);
    }
  }
}





/*
 * Give the number of missrounded results for 
 * crlibm, libmultim and libm.
 */
void test_all() {
  int counter=0;
  long long int failures_crlibm=0,
                failures_libm=0,
                failures_ibm=0;
  long long int i;

  i=0; 
  while(1+1==2){
    input.d = randfun();
    res_crlibm.d = testfun_crlibm(input.d);
    res_libm.d = testfun_libm(input.d);
#ifdef HAVE_MATHLIB_H
    if(rnd==0)
      res_ibm.d = testfun_ibm(input.d);
#endif
    mpfr_set_d(mp_inpt, input.d, GMP_RNDN);
    testfun_mpfr(mp_res, mp_inpt, rnd);
    res_mpfr.d = mpfr_get_d(mp_res, rnd);


#if PRINT_NAN
    if(1){
#else
      if((res_mpfr.i[HI_ENDIAN] & 0x7ff00000) != 0x7ff00000){
#endif
	if( (res_crlibm.i[LO_ENDIAN] != res_mpfr.i[LO_ENDIAN]) 
	    || (res_crlibm.i[HI_ENDIAN] != res_mpfr.i[HI_ENDIAN]) ) 
	  {
#if DETAILED_REPORT	  
	  printf("CRLIBM ERROR  x=%.70e \n            (%8x %8x) \n", 
		 input.d, 
		 input.i[HI_ENDIAN], 
		 input.i[LO_ENDIAN]);
	  printf("crlibm gives    %.50e \n         (%8x %8x) \n", 
		 res_crlibm.d, 
		 res_crlibm.i[HI_ENDIAN], 
		 res_crlibm.i[LO_ENDIAN]);
	  printf("MPFR gives %.50e \n         (%8x %8x) \n\n", 
		 res_mpfr.d, 
		 res_mpfr.i[HI_ENDIAN], 
		 res_mpfr.i[LO_ENDIAN]);
#endif
	  failures_crlibm++;
	  }
	
	if( (res_libm.i[LO_ENDIAN] != res_mpfr.i[LO_ENDIAN]) 
	    || (res_libm.i[HI_ENDIAN] != res_mpfr.i[HI_ENDIAN]) ) failures_libm++;
	  
#ifdef HAVE_MATHLIB_H
	if(rnd==0
	   && ((res_ibm.i[LO_ENDIAN] != res_mpfr.i[LO_ENDIAN]) 
	       || (res_ibm.i[HI_ENDIAN] != res_mpfr.i[HI_ENDIAN]) )) 
	  {
#if DETAILED_REPORT
	      printf("IBM ULTIM ERROR  x=%.50e \n            (%8x %8x) \n", 
		     input.d, 
		     input.i[HI_ENDIAN], 
		     input.i[LO_ENDIAN]);
	      printf("libultim gives    %.50e \n         (%8x %8x) \n", 
		     res_ibm.d, 
		     res_ibm.i[HI_ENDIAN], 
		     res_ibm.i[LO_ENDIAN]);
	      printf("MPFR gives %.50e \n         (%8x %8x) \n\n", 
		     res_mpfr.d, 
	       res_mpfr.i[HI_ENDIAN], 
		     res_mpfr.i[LO_ENDIAN]);
#endif
	      failures_ibm++;
	  }
#endif
      }
      i++;
      if((i % 1000000)==0) {
	printf(" CRLIBM : %lld failures out of %lld (ratio %e) \n",failures_crlibm, i,
	       ((double)failures_crlibm)/(double)i);
	printf(" LIBM   : %lld failures out of %lld (ratio %e) \n",failures_libm, i,
	       ((double)failures_libm)/(double)i);
#ifdef HAVE_MATHLIB_H
	printf(" IBM    : %lld failures out of %lld (ratio %e) \n \n",failures_ibm, i,
	       ((double)failures_ibm)/(double)i);
#endif
      }
  }
}


void test_random_gen() {

  double x, min, max;
  min=0;max=0;

  while(1+1==2) {
    x = randfun();
    if(x<min) min=x;
    if(x>max) max=x;
    printf("x=%f    (%f %f)\n", 
	   x, min, max
	   );
  }
}
  


void usage(char *fct_name){
  fprintf (stderr, "\n\n Soak-test for crlibm \n");
  fprintf (stderr, "Usage: %s name [RN/RU/RD/RZ] \n", fct_name);
  fprintf (stderr, " name          : name of function to test \n");
  fprintf (stderr, " [RN/RU/RD/RZ] : rounding mode, (default RN) \n");
  exit (1);
}



int main (int argc, char *argv[]) 
{ 
  crlibm_init();

  if ((argc < 2)||(argc > 3)) usage(argv[0]);
    
  if (argc == 3)
    if      (strcmp(argv[2],"RU")==0) scs_rnd_mode = 2;
    else if (strcmp(argv[2],"RD")==0) scs_rnd_mode = 3;
    else if (strcmp(argv[2],"RZ")==0) scs_rnd_mode = 4;
    else scs_rnd_mode = 1;
    
  switch(scs_rnd_mode){
  case(2) :    rnd = GMP_RNDU;    break;
  case(3) :    rnd = GMP_RNDD;    break;
  case(4) :    rnd = GMP_RNDZ;    break;
  default :    rnd = GMP_RNDN;  }

  printf("Testing %s function with rounding mode : %d \n", argv[1], rnd);

  srand(33);

  init(argv[1]);
#if 0
  test_random_gen();
#endif

  //   test();
   test_all();

  return 0;
}


#else
int main (int argc, char *argv[]) 
{ 
  printf("Sorry, I need to be compiled against MPFR\n");
}
#endif
