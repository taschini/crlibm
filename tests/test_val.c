#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "crlibm.h"
#include "crlibm_private.h"

#ifdef HAVE_MATHLIB_H
#include <MathLib.h>
#endif

#ifdef HAVE_MPFR_H
#include <gmp.h>
#include <mpfr.h>
#endif












void usage(char *fct_name){
  fprintf (stderr, "\nCompares results between different library (mpfr, Ziv, libm, cr_libm) \n");
  fprintf (stderr, "Usage: %s fun mode val \n", fct_name);
  fprintf (stderr, " fun : name of function to test \n");
  fprintf (stderr, " mode : rounding mode [RN, RD, RU, RZ] \n");
  fprintf (stderr, " val  : double precision input number \n");
  exit (1);
}




int main (int argc, char *argv[]) 
{ 
  char* rounding_mode;
  char* function_name;
  double input;

  
  db_number res_crlibm, res_mpfr, res_ibm, res_libm;
#ifdef HAVE_MPFR_H
  mp_rnd_t mpfr_rnd_mode;
  mpfr_t mp_res, mp_input; 
#endif


/* The random number generator (unused here) */
double (*randfun)       () = NULL;
/* The function we test */
double (*testfun_crlibm)() = NULL;
/* The function we trust */
int    (*testfun_mpfr)  () = NULL;
/* The function to show off against for accuracy  */
double (*testfun_libm)  () = NULL;
/* The function to show off against for performance */
double (*testfun_ibm)   () = NULL;


  if ((argc != 4)) usage(argv[0]);
  else{
    function_name = argv[1];
    rounding_mode = argv[2];
    sscanf(argv[3],"%le", &input);
    
    crlibm_init();
#ifdef HAVE_MATHLIB_H
  Init_Lib(); /* we don't save the state, no need here */ 
#endif

    test_init(/* pointers to returned value */
	       &randfun, 
	       &testfun_crlibm, 
	       &testfun_mpfr,
	       &testfun_libm,
	       &testfun_ibm,
	       /* arguments */
	       function_name,
	       rounding_mode ) ;
    


#ifdef HAVE_MPFR_H  /* stop here if MPFR not present */
    mpfr_init2(mp_res,  153);
    mpfr_init2(mp_input, 53);
    if      (strcmp(rounding_mode,"RU")==0) mpfr_rnd_mode = GMP_RNDU;
    else if (strcmp(rounding_mode,"RD")==0) mpfr_rnd_mode = GMP_RNDD;
    else if (strcmp(rounding_mode,"RZ")==0) mpfr_rnd_mode = GMP_RNDZ;
    else {
      mpfr_rnd_mode = GMP_RNDN; 
      rounding_mode="RN" ;
    }
#endif



  printf("cr_libm    : "); 
  fflush(stdout); /* To help debugging */
  if(testfun_crlibm != NULL)   {
    res_crlibm.d = testfun_crlibm(input);
    printf("%.50e  %8x %8x\n", 
	   res_crlibm.d, 
	   res_crlibm.i[HI_ENDIAN], 
	   res_crlibm.i[LO_ENDIAN] );
  }
  else 
    printf("Not available\n");
  fflush(stdout);
  
#ifdef HAVE_MPFR_H
  printf("mpfr_libm  : ");
  fflush(stdout);
  if(testfun_mpfr != NULL){
    mpfr_set_d(mp_input, input,  GMP_RNDN);
    testfun_mpfr(mp_res, mp_input, mpfr_rnd_mode);
    res_mpfr.d = mpfr_get_d(mp_res, mpfr_rnd_mode);
    printf("%.50e  %8x %8x \n", 
	   res_mpfr.d, 
	   res_mpfr.i[HI_ENDIAN], 
	   res_mpfr.i[LO_ENDIAN] );
  }else 
    printf("Not available\n");
  fflush(stdout);
#endif



#ifdef HAVE_MATHLIB_H
  printf("ibm_libm   : ");
  fflush(stdout);
  if(testfun_ibm != NULL)  {
    res_ibm.d = testfun_ibm(input);
    printf("%.50e  %8x %8x \n", 
	   res_ibm.d, 
	   res_ibm.i[HI_ENDIAN], 
	   res_ibm.i[LO_ENDIAN] );
  }
  else 
    printf("Not available\n");
  fflush(stdout);
#endif



  /* Last in the list because it segfaults more often than the
     others.  */
  printf("System libm : ");
  fflush(stdout);
  if(testfun_libm != NULL) 
    {
      res_libm.d = testfun_libm(input);
      printf("%.50e  %8x %8x \n", 
	     res_libm.d, 
	     res_libm.i[HI_ENDIAN], 
	     res_libm.i[LO_ENDIAN]) ;
    }
  else
    printf("Not available\n");
  fflush(stdout);

  /* release memory */
#ifdef HAVE_MPFR_H
  mpfr_clear(mp_input);
  mpfr_clear(mp_res);
#endif

  return 0;
  }
}


