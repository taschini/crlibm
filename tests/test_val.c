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
  /*fprintf (stderr, "\n%s: single-value test for crlibm and other mathematical libraries\n", fct_name);*/
  fprintf (stderr, "\nUsage: %s [-x] function (RN|RU|RD|RZ) value\n", fct_name);
  fprintf (stderr, " function      : name of function to test \n");
  fprintf (stderr, " (RN|RU|RD|RZ) : rounding mode, \n");
  fprintf (stderr, " value         : double precision input number (in hexadecimal if option -x is given) \n");
  exit (1);
}



int main (int argc, char *argv[]) 
{ 
  char* option;
  char* function_name;
  char* rounding_mode;
  double worstcase;

  
  db_number input, res_crlibm, res_mpfr, res_ibm, res_libmcr, res_libm;
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
  double (*testfun_libultim)   () = NULL;
  /* The last to join the family */
  double (*testfun_libmcr)   () = NULL;


  if ((argc != 4) && argc!=5) usage(argv[0]);
  if(argc == 4) {
    function_name = argv[1];
    rounding_mode = argv[2];
    sscanf(argv[3],"%le", &input.d);
  }
  if(argc == 5) {
    option = argv[1];
    if(!(option[0]=='-' && option[1]=='x')) usage(argv[0]);
    function_name = argv[2];
    rounding_mode = argv[3];
    sscanf(argv[4],"%llx", &input.l);
  }
  crlibm_init();
#ifdef HAVE_MATHLIB_H
  Init_Lib(); /* we don't save the state, no need here */ 
#endif

    test_init(/* pointers to returned value */
	      &randfun, /* unused */ 
	      &randfun, /* unused */ 
	      &testfun_crlibm, 
	      &testfun_mpfr,
	      &testfun_libultim,
	      &testfun_libmcr,
	      &testfun_libm,
	      &worstcase,
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



  printf("Input      : %.50e  %08x %08x\n", input.d, input.i[HI], input.i[LO] ); 
  printf("cr_libm    : "); 
  fflush(stdout); /* To help debugging */
  if(testfun_crlibm != NULL)   {
    res_crlibm.d = testfun_crlibm(input.d);
    printf("%.50e  %08x %08x\n", 
	   res_crlibm.d, 
	   res_crlibm.i[HI], 
	   res_crlibm.i[LO] );
  }
  else 
    printf("Not available\n");
  fflush(stdout);
  
#ifdef HAVE_MPFR_H
  printf("mpfr_libm  : ");
  fflush(stdout);
  if(testfun_mpfr != NULL){
    mpfr_set_d(mp_input, input.d,  GMP_RNDN);
    testfun_mpfr(mp_res, mp_input, mpfr_rnd_mode);
    res_mpfr.d = mpfr_get_d(mp_res, mpfr_rnd_mode);
    printf("%.50e  %08x %08x \n", 
	   res_mpfr.d, 
	   res_mpfr.i[HI], 
	   res_mpfr.i[LO] );
  }else 
    printf("Not available\n");
  fflush(stdout);
#endif



#ifdef HAVE_MATHLIB_H
  printf("libultim   : ");
  fflush(stdout);
  if(testfun_libultim != NULL)  {
    res_ibm.d = testfun_libultim(input.d);
    printf("%.50e  %08x %08x \n", 
	   res_ibm.d, 
	   res_ibm.i[HI], 
	   res_ibm.i[LO] );
  }
  else 
    printf("Not available\n");
  fflush(stdout);
#endif

#ifdef HAVE_LIBMCR_H
  printf("libmcr     : ");
  fflush(stdout);
  if(testfun_libmcr != NULL)  {
    res_libmcr.d = testfun_libmcr(input.d);
    printf("%.50e  %08x %08x \n", 
	   res_libmcr.d, 
	   res_libmcr.i[HI], 
	   res_libmcr.i[LO] );
  }
  else 
    printf("Not available\n");
  fflush(stdout);
#endif



  /* Last in the list because it segfaults more often than the
     others.  */
  printf("System libm: ");
  fflush(stdout);
  if(testfun_libm != NULL) 
    {
      res_libm.d = testfun_libm(input.d);
      printf("%.50e  %08x %08x \n", 
	     res_libm.d, 
	     res_libm.i[HI], 
	     res_libm.i[LO]) ;
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



