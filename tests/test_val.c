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

/* log2:
 * x=4.33808546963214027000663771883653290038704313524560e-309
 *
 * Problematic value for IBM lib : 
 *  # log2 :
 *     * 3.93340915000562568776548176218927720232626410258349e-309
 *
 */
double (*testfun_libm)   () = NULL;
double (*testfun_crlibm) () = NULL;
double (*testfun_ibm)    () = NULL;
int    (*testfun_mpfr)   () = NULL;


void testfun(char *function, double inpt){
  db_number res_cr, res_ibm, res_libm, res_mpfr;

#ifdef HAVE_MPFR_H
  mpfr_t mp_res, mp_inpt; 

  mpfr_init2(mp_res,  160);
  mpfr_init2(mp_inpt, 53);
#endif
  res_cr.d = inpt;
  printf("Input number %.70e  %8x %8x \n\n", 
	 res_cr.d, 
	 res_cr.i[HI_ENDIAN], 
	 res_cr.i[LO_ENDIAN]);
  
  if (strcmp (function, "exp_rn") == 0)
    {
      testfun_libm   = exp;
      testfun_crlibm = exp_rn;
#ifdef HAVE_MATHLIB_H
      testfun_ibm    = uexp;
#endif
#ifdef HAVE_MPFR_H
      testfun_mpfr   = mpfr_exp;
#endif
    }

  else if (strcmp (function, "log") == 0)
    {
      testfun_libm   = log;
      testfun_crlibm = log_rn;
#ifdef HAVE_MATHLIB_H
      testfun_ibm    = ulog;
#endif
#ifdef HAVE_MPFR_H
      testfun_mpfr   = mpfr_log;
#endif
    }

  else if (strcmp (function, "log2") == 0)
    {
      testfun_libm   = NULL /* doesn't exist ? */;
      testfun_crlibm = log2_rn;
#ifdef HAVE_MATHLIB_H
      testfun_ibm    = ulog2;
#endif
#ifdef HAVE_MPFR_H
      testfun_mpfr   = mpfr_log2;
#endif
    }

  else if (strcmp (function, "log10") == 0)
    {
      testfun_libm   = log10;
      testfun_crlibm = log10_rn;
#ifdef HAVE_MATHLIB_H
      testfun_ibm    = NULL;
#endif
#ifdef HAVE_MPFR_H
      testfun_mpfr   = mpfr_log10;
#endif
    }

  else if (strcmp (function, "tan") == 0)
    {
      testfun_libm   = NULL ; /*segfaults on Linux, I don't know why*/
      testfun_crlibm = tan_rn;
#ifdef HAVE_MATHLIB_H
      testfun_ibm    = utan;
#endif
#ifdef HAVE_MPFR_H
      testfun_mpfr   = mpfr_tan;
#endif
    }

  else if (strcmp (function, "sin") == 0)
    {
      testfun_libm   = sin;
      testfun_crlibm = sin_rn;
#ifdef HAVE_MATHLIB_H
      testfun_ibm    = usin;
#endif
#ifdef HAVE_MPFR_H
      testfun_mpfr   = mpfr_sin;
#endif
    }

  else if (strcmp (function, "cos") == 0)
    {
      testfun_libm   = cos;
      testfun_crlibm = cos_rn;
#ifdef HAVE_MATHLIB_H
      testfun_ibm    = ucos;
#endif
#ifdef HAVE_MPFR_H
      testfun_mpfr   = mpfr_cos;
#endif
    }

  else
    {
      fprintf (stderr, "Unknown function: %s\n", function);
      exit (1);
    }


  printf("cr_libm    : "); 
  fflush(stdout); /* To help debugging */
  if(testfun_crlibm != NULL)   {
    res_cr.d = testfun_crlibm(inpt);
    printf("%.50e  %8x %8x\n", 
	   res_cr.d, 
	   res_cr.i[HI_ENDIAN], 
	   res_cr.i[LO_ENDIAN] );
  }
  else 
    printf("Not available\n");
  fflush(stdout);
  
#ifdef HAVE_MATHLIB_H
  printf("ibm_libm   : ");
  fflush(stdout);
  if(testfun_ibm != NULL)  {
    res_ibm.d = testfun_ibm(inpt);
    printf("%.50e  %8x %8x \n", 
	   res_ibm.d, 
	   res_ibm.i[HI_ENDIAN], 
	   res_ibm.i[LO_ENDIAN] );
  }
  else 
    printf("Not available\n");
  fflush(stdout);
#endif

#ifdef HAVE_MPFR_H
  printf("mpfr_libm  : ");
  fflush(stdout);
  if(testfun_mpfr != NULL){
    mpfr_set_d(mp_inpt, inpt, GMP_RNDN);
    testfun_mpfr(mp_res, mp_inpt, GMP_RNDN);
    res_mpfr.d = mpfr_get_d1(mp_res);
    printf("%.50e  %8x %8x \n", 
	   res_mpfr.d, 
	   res_mpfr.i[HI_ENDIAN], 
	   res_mpfr.i[LO_ENDIAN] );
  }else 
    printf("Not available\n");
  fflush(stdout);
  
  
#endif

  /* Last in the list because it segfaults more often than the
     others.  */
  printf("System libm : ");
  fflush(stdout);
  if(testfun_libm != NULL) 
    {
      res_libm.d = testfun_libm(inpt);
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
  mpfr_clear(mp_inpt);
  mpfr_clear(mp_res);
#endif
}



void usage(char *fct_name){
  fprintf (stderr, "\nCompares results between different library (mpfr, Ziv, libm, cr_libm) \n");
  fprintf (stderr, "Usage: %s name vals \n", fct_name);
  fprintf (stderr, " name : name of function to test \n");
  fprintf (stderr, " val  : double precision input number \n");
  exit (1);
}



int main (int argc, char *argv[]) 
{ 
  double input=0;
  
  crlibm_init();
#ifdef HAVE_MATHLIB_H
  Init_Lib(); /* we don't save the state, no need here */ 
#endif

  if (argc != 3) usage(argv[0]);
    
  sscanf(argv[2],"%le", &input);
  testfun(argv[1], input);



  return 0;
}

