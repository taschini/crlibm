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
db_number input, res_crlibm, res_mpfr, res_libultim, res_libmcr, res_libm;
mpfr_t mp_res, mp_inpt; 

/* The random number generator*/
double (*randfun)       () = NULL;
/* Another unused random number generator*/
double (*randfun_perf)       () = NULL;

/* The function we test */
double (*testfun_crlibm)() = NULL;
/* The function we trust */
int    (*testfun_mpfr)  () = NULL;
/* The function to show off against for accuracy  */
double (*testfun_libm)  () = NULL;
/* The function to show off against for performance */
double (*testfun_libultim)   () = NULL;
/*  */
double (*testfun_libmcr)   () = NULL;


/*
 * Rounding mode to test
 */
mp_rnd_t mpfr_rnd_mode;







/*
 * Give the number of missrounded results 
 */

void test_all() {
  int counter=0;
  long long int 
    failures_crlibm=0,
    failures_libm=0,
    failures_libultim=0,
    failures_libmcr=0;
  long long int i;

  i=0; 
  while(1+1==2){
    input.d = randfun();
    res_crlibm.d = testfun_crlibm(input.d);
    res_libm.d = testfun_libm(input.d);

#ifdef HAVE_MATHLIB_H
    if(mpfr_rnd_mode==GMP_RNDN && testfun_libultim != NULL) /* not all the functions are in libultim */
      res_libultim.d = testfun_libultim(input.d);
#endif
#ifdef HAVE_LIBMCR_H
    if(mpfr_rnd_mode==GMP_RNDN && testfun_libmcr != NULL) /* not all the functions are in libultim */
      res_libmcr.d = testfun_libmcr(input.d);
#endif
    mpfr_set_d(mp_inpt, input.d, GMP_RNDN);
    testfun_mpfr(mp_res, mp_inpt, mpfr_rnd_mode);
    res_mpfr.d = mpfr_get_d(mp_res, mpfr_rnd_mode);


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
	if(mpfr_rnd_mode==0  && testfun_libultim != NULL 
	   && ((res_libultim.i[LO_ENDIAN] != res_mpfr.i[LO_ENDIAN]) 
	       || (res_libultim.i[HI_ENDIAN] != res_mpfr.i[HI_ENDIAN]) )) 
	  {
#if DETAILED_REPORT
	      printf("IBM ULTIM ERROR  x=%.50e \n            (%8x %8x) \n", 
		     input.d, 
		     input.i[HI_ENDIAN], 
		     input.i[LO_ENDIAN]);
	      printf("libultim gives    %.50e \n         (%8x %8x) \n", 
		     res_libultim.d, 
		     res_libultim.i[HI_ENDIAN], 
		     res_libultim.i[LO_ENDIAN]);
	      printf("MPFR gives %.50e \n         (%8x %8x) \n\n", 
		     res_mpfr.d, 
	       res_mpfr.i[HI_ENDIAN], 
		     res_mpfr.i[LO_ENDIAN]);
#endif
	      failures_libultim++;
	  }
#endif
#ifdef HAVE_LIBMCR_H
	if(mpfr_rnd_mode==0  && testfun_libmcr != NULL 
	   && ((res_libmcr.i[LO_ENDIAN] != res_mpfr.i[LO_ENDIAN]) 
	       || (res_libmcr.i[HI_ENDIAN] != res_mpfr.i[HI_ENDIAN]) )) 
	  {
#if DETAILED_REPORT
	      printf("LIBMCR ERROR  x=%.50e \n            (%8x %8x) \n", 
		     input.d, 
		     input.i[HI_ENDIAN], 
		     input.i[LO_ENDIAN]);
	      printf("libmcr gives    %.50e \n         (%8x %8x) \n", 
		     res_libmcr.d, 
		     res_libmcr.i[HI_ENDIAN], 
		     res_libmcr.i[LO_ENDIAN]);
	      printf("MPFR gives %.50e \n         (%8x %8x) \n\n", 
		     res_mpfr.d, 
	       res_mpfr.i[HI_ENDIAN], 
		     res_mpfr.i[LO_ENDIAN]);
#endif
	      failures_libmcr++;
	  }
#endif
      }
      i++;
      if((i % 1000000)==0) {
	printf(" CRLIBM       : %lld failures out of %lld (ratio %e) \n",failures_crlibm, i,
	       ((double)failures_crlibm)/(double)i);
	printf(" LIBM         : %lld failures out of %lld (ratio %e) \n",failures_libm, i,
	       ((double)failures_libm)/(double)i);
#ifdef HAVE_MATHLIB_H
	printf(" IBM LIBULTIM : %lld failures out of %lld (ratio %e) \n",failures_libultim, i,
	       ((double)failures_libultim)/(double)i);
#endif
#ifdef HAVE_LIBMCR_H
	printf(" SUN LIBMCR   : %lld failures out of %lld (ratio %e) \n",failures_libmcr, i,
	       ((double)failures_libmcr)/(double)i);
#endif
	printf("\n");
	
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
  /* fprintf (stderr, "\n%s: Soak-test for crlibm and other mathematical libraries \n", fct_name); */
  fprintf (stderr, "\nUsage: %s function (RN|RU|RD|RZ) seed \n", fct_name);
  fprintf (stderr, " function      : name of function to test \n");
  fprintf (stderr, " (RN|RU|RD|RZ) : rounding mode, \n");
  fprintf (stderr, " seed          : integer seed for the random number generator \n");
  exit (1);
}



int main (int argc, char *argv[]) 
{ 
  char* rounding_mode;
  char* function_name;
  int seed;
  double worstcase;

  if ((argc != 4)) usage(argv[0]);
  else{
    function_name = argv[1];
    rounding_mode = argv[2];
    sscanf(argv[3],"%d", &seed);
    
    crlibm_init();

    test_init(/* pointers to returned value */
	      &randfun_perf, /* unused here*/ 
	      &randfun, 
	      &testfun_crlibm, 
	      &testfun_mpfr,
	      &testfun_libultim,
	      &testfun_libmcr,
	      &testfun_libm,
	      &worstcase,
	      /* arguments */
	      function_name,
	      rounding_mode ) ;
    
    mpfr_init2(mp_res,  153);
    mpfr_init2(mp_inpt, 53);
    if      (strcmp(rounding_mode,"RU")==0) mpfr_rnd_mode = GMP_RNDU;
    else if (strcmp(rounding_mode,"RD")==0) mpfr_rnd_mode = GMP_RNDD;
    else if (strcmp(rounding_mode,"RZ")==0) mpfr_rnd_mode = GMP_RNDZ;
    else {
      mpfr_rnd_mode = GMP_RNDN; 
      rounding_mode="RN" ;
    }

  printf("Testing %s function with rounding mode : %s \n", function_name, rounding_mode);

  //srand(seed);

#if 0
  test_random_gen();
#endif

  //   test();
  test_all();

  return 0;
  }
}


#else
int main (int argc, char *argv[]) 
{ 
  printf("Sorry, I need to be compiled against MPFR\n");
}
#endif
