
/* 
Beware to compile without optimizations

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#include "crlibm.h"
#include "crlibm_private.h"
#include "test_common.h"

#include "scs_lib/tests/tbx_timing.h"

#ifdef HAVE_MATHLIB_H
#include <MathLib.h>
#endif

#ifdef HAVE_MPFR_H
#include <gmp.h>
#include <mpfr.h>
#endif



#define DETAILED_REPORT 0

/* If set, the behaviour of the function with respect to cache memory
   will be tested*/
#define TEST_CACHE 0




#if EVAL_PERF==1  
/* counter of calls to the second step */
extern int crlibm_second_step_taken; 
int crlibm_first_step_taken;
#endif


/* Basic-like programming with global variables: */
db_number input, res_crlibm, res_mpfr, res_ibm, res_libm;

double worst_case; /* worst case for round to nearest only */

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


/* TESTSIZE doubles should be enough to flush the cache */
#define TESTSIZE 200000 
double inputs[TESTSIZE];



void init(char *function){
  crlibm_init();

  randfun        = rand_generic; /* the default random function */
  worst_case = 0.0;
  if (strcmp (function, "exp") == 0)
    {
      randfun        = rand_for_exp;
      worst_case= .75417527749959590085206221024712557043923055744016892276704311370849609375e-9;
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
      randfun        = rand_for_log;
      worst_case=0.4009793462309855760053830468258630076242931610568335144339734234840014178511334897967240437927437320e-115;
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

  else if (strcmp (function, "tan") == 0)
    {
      testfun_libm   = tan   ;
      testfun_crlibm = tan_rn;
#ifdef HAVE_MATHLIB_H
      testfun_ibm    = utan;
#endif
#ifdef HAVE_MPFR_H
      testfun_mpfr   = mpfr_tan;
#endif
    }

  else if (strcmp (function, "cotan") == 0)
    {
      testfun_libm   = NULL;
      testfun_crlibm = cotan_rn;
#ifdef HAVE_MATHLIB_H
/*       testfun_ibm    = ucotan; doesn't exist ?*/
#endif
#ifdef HAVE_MPFR_H
      testfun_mpfr   = NULL; /* should check */
#endif
    }

    
  else if (strcmp (function, "log") == 0)
    {
      testfun_crlibm = log_rn;
#ifdef HAVE_MPFR_H
      testfun_mpfr   = mpfr_log;
#endif
      randfun        = rand_generic;
    }
  else if (strcmp (function, "log2") == 0)
    {
      testfun_crlibm = log2_rn;
#ifdef HAVE_MPFR_H
      testfun_mpfr   = mpfr_log2;
#endif
      randfun        = rand_generic;
    }
  else if (strcmp (function, "log10") == 0)
    {
      testfun_crlibm = log10_rn;
#ifdef HAVE_MPFR_H
      testfun_mpfr   = mpfr_log10;
#endif
      randfun        = rand_generic;
    }
  else if (strcmp (function, "sin") == 0)
    {
      testfun_crlibm = sin_rn;
#ifdef HAVE_MPFR_H
      testfun_mpfr   = mpfr_sin;
#endif
      randfun        = rand_generic;
    }
  else if (strcmp (function, "cos") == 0)
    {
      testfun_crlibm = cos_rn;
#ifdef HAVE_MPFR_H
      testfun_mpfr   = mpfr_cos;
#endif
      randfun        = rand_generic;
    }
  else if (strcmp (function, "tan") == 0)
    {
      testfun_crlibm = tan_rn;
#ifdef HAVE_MPFR_H
      testfun_mpfr   = mpfr_tan;
#endif
    }
#if 0 /* mpfr has no cotan? */ 
  else if (strcmp (function, "cotan") == 0)
    {
      testfun_crlibm = cotan_rn;
#ifdef HAVE_MPFR_H
      testfun_mpfr   = mpfr_cotan;
#endif
    }
#endif
  else
    {
      fprintf (stderr, "Unknown function: %s\n", function);
      exit (1);
    }
}




void fill_and_flush(int seed) {
  int i;
  srandom(seed);
  for (i=0; i<TESTSIZE; i++)
    inputs[i] = randfun();
}






int main (int argc, char *argv[]) 
{ 
  
  int i, j, k, n=100000;
  double input, result;
  tbx_tick_t   t1, t2; 
  unsigned long long 
    dt, dtmin, 
    libm_dtmin, libm_dtmax, libm_dtsum, libm_dtwc,
    crlibm_dtmin, crlibm_dtmax, crlibm_dtsum, crlibm_dtwc,
    mpfr_dtmin, mpfr_dtmax, mpfr_dtsum, mpfr_dtwc,
    ibm_dtmin, ibm_dtmax, ibm_dtsum, ibm_dtwc,
    dtsum, min_dtsum;
  unsigned long seed = 42;
  int output_latex=0;
#ifdef HAVE_MPFR_H  
  mpfr_t mp_res, mp_inpt; 
#endif
  short Original_Mode;

  /************  INITIALISATIONS  *********************/

  if ((argc != 3)) {
    printf("Usage: %s function iterations\n", argv[0]);
    exit(1);
  }
 
  init(argv[1]);

  sscanf(argv[2],"%d", &n);


#ifdef HAVE_MPFR_H
  mpfr_init2(mp_res,  53);
  mpfr_init2(mp_inpt, 53);
#endif


  crlibm_dtmin=1<<30; crlibm_dtmax=0; crlibm_dtsum=0;
  libm_dtmin=1<<30;   libm_dtmax=0;   libm_dtsum=0;
  mpfr_dtmin=1<<30;   mpfr_dtmax=0;   mpfr_dtsum=0;
  ibm_dtmin=1<<30;    ibm_dtmax=0;    ibm_dtsum=0;



#if TEST_CACHE
  /************  TESTS WITH CACHES  *********************/
  /* First tests in real conditions, where cache considerations
     matter */

  /* libm */
  printf("\nLIBM\n");
    for(i=1; i<=10000; i*=10) { /* i=1,10,100...*/
      min_dtsum=1<<30; 
      for(k=0;k<10;k++) { /* do the whole test 10 times and take the min */
	fill_and_flush(n);
	dtsum=0;
	for (j=0; j<i; j++) {
	  input=inputs[i];
	  TBX_GET_TICK(t1);
	  result = testfun_libm(input);
	  TBX_GET_TICK(t2);
	  dt = TBX_TICK_RAW_DIFF(t1, t2); 
	  dtsum += dt;
	}
	if (dtsum < min_dtsum) min_dtsum=dtsum; 
      }
      printf("  %d loops: \t avg time = %f ticks\n",i, ((double)min_dtsum)/i);
    }

  printf("\nCRLIBM\n");
    for(i=1; i<=10000; i*=10) { /* i=1,10,100...*/
      min_dtsum=1<<30; 
      for(k=0;k<10;k++) { /* do the whole test 10 times and take the min */
	fill_and_flush(n);
	dtsum=0;
	for (j=0; j<i; j++) {
	  input=inputs[i];
	  TBX_GET_TICK(t1);
	  result = testfun_crlibm(input);
	  TBX_GET_TICK(t2);
	  dt = TBX_TICK_RAW_DIFF(t1, t2); 
	  dtsum += dt;
	}
	if (dtsum < min_dtsum) min_dtsum=dtsum; 
      }
      printf("  %d loops: \t avg time = %f ticks\n",i, ((double)min_dtsum)/i);
    }

#ifdef HAVE_MATHLIB_H
  printf("\nIBM\n");
    for(i=1; i<=10000; i*=10) { /* i=1,10,100...*/
      min_dtsum=1<<30; 
      for(k=0;k<10;k++) { /* do the whole test 10 times and take the min */
	fill_and_flush(n);
	dtsum=0;
	for (j=0; j<i; j++) {
	  input=inputs[i];
	  TBX_GET_TICK(t1);
	  result = testfun_ibm(input);
	  TBX_GET_TICK(t2);
	  dt = TBX_TICK_RAW_DIFF(t1, t2); 
	  dtsum += dt;
	}
	if (dtsum < min_dtsum) min_dtsum=dtsum; 
      }
      printf("  %d loops: \t avg time = %f ticks\n",i, ((double)min_dtsum)/i);
    }
#endif

#endif /* TEST_CACHE*/
  /************  TESTS WITHOUT CACHES  *********************/
  srandom(n);

#if EVAL_PERF==1  
  crlibm_second_step_taken=0; 
#endif


  /* take the min of 10 identical calls to leverage interruptions */
  /* As a consequence, the cache impact of these calls disappear...*/

  for(i=0; i< n; i++){ 
    input = randfun();

    /* libm timing */
    dtmin=1<<30;
    for(j=0; j<5; j++) {
      TBX_GET_TICK(t1);
#ifdef CRLIBM_TYPECPU_POWERPC
      for(k=0; k<50;k++)
#endif
	result = testfun_libm(input);
      TBX_GET_TICK(t2);
      dt = TBX_TICK_RAW_DIFF(t1, t2); 
      if (dt<dtmin)  dtmin=dt;
    }
    libm_dtsum+=dtmin;
    if (dtmin<libm_dtmin)  libm_dtmin=dtmin;
    if (dtmin>libm_dtmax)  libm_dtmax=dtmin;
#if DETAILED_REPORT
    printf("\n input=%1.15e\tTlibm=%lld", input, dtmin);
#endif /*DETAILED_REPORT*/

    /* crlibm timing */
    dtmin=1<<30;
    for(j=0; j<5; j++) {
      TBX_GET_TICK(t1);
#ifdef CRLIBM_TYPECPU_POWERPC
      for(k=0; k<50;k++)
#endif
      result = testfun_crlibm(input);
      TBX_GET_TICK(t2);
      dt = TBX_TICK_RAW_DIFF(t1, t2); 
      if (dt<dtmin)  dtmin=dt;
    }
    crlibm_dtsum+=dtmin;
    if (dtmin<crlibm_dtmin)  crlibm_dtmin=dtmin;
    if (dtmin>crlibm_dtmax)  crlibm_dtmax=dtmin;

#if DETAILED_REPORT
    printf("\tTcrlibm=%lld", dtmin);
#endif /*DETAILED_REPORT*/


#ifdef HAVE_MPFR_H
    /* mpfr timing */
    dtmin=1<<30;
    for(j=0; j<5; j++) {
      TBX_GET_TICK(t1);
#ifdef CRLIBM_TYPECPU_POWERPC
      for(k=0; k<50;k++) {
#endif
      mpfr_set_d(mp_inpt, input, GMP_RNDN);
      testfun_mpfr(mp_res, mp_inpt, GMP_RNDN);
      result = mpfr_get_d1(mp_res);
#ifdef CRLIBM_TYPECPU_POWERPC
      }
#endif
      TBX_GET_TICK(t2);
      dt = TBX_TICK_RAW_DIFF(t1, t2); 
      if (dt<dtmin)  dtmin=dt;
    }
    mpfr_dtsum+=dtmin;
    if (dtmin<mpfr_dtmin)  mpfr_dtmin=dtmin;
    if (dtmin>mpfr_dtmax)  mpfr_dtmax=dtmin;
#if DETAILED_REPORT
    printf("\tTmpfr=%lld",dtmin);
#endif /*DETAILED_REPORT*/
#endif /* HAVE_MPFR */

#ifdef HAVE_MATHLIB_H
    Original_Mode = Init_Lib(); 
    dtmin=1<<30;
    for(j=0; j<5; j++) {
      TBX_GET_TICK(t1);
#ifdef CRLIBM_TYPECPU_POWERPC
      for(k=0; k<50;k++)
#endif
      result = testfun_ibm(input);
      TBX_GET_TICK(t2);
      dt = TBX_TICK_RAW_DIFF(t1, t2); 
      if (dt<dtmin)  dtmin=dt;
    }
    Exit_Lib(Original_Mode);
    ibm_dtsum+=dtmin;
    if (dtmin<ibm_dtmin)  ibm_dtmin=dtmin;
    if (dtmin>ibm_dtmax)  ibm_dtmax=dtmin;
#if DETAILED_REPORT
    printf("\tTibm=%lld",dtmin);
#endif /*DETAILED_REPORT*/
#endif /*HAVE_MATHLIB_H*/

  } 


#if EVAL_PERF==1  
  printf("\nCRLIBM : Second step taken %d times out of %d\n",
	 crlibm_second_step_taken/20, n );
#endif


  /************  WORST CASE TESTS   *********************/
  /* worst case test */
  input = worst_case;
  /* libm timing */
  dtmin=1<<30;
  for(j=0; j<10; j++) {
    TBX_GET_TICK(t1);
    result = testfun_libm(input);
    TBX_GET_TICK(t2);
    dt = TBX_TICK_RAW_DIFF(t1, t2); 
    if (dt<dtmin)  dtmin=dt;
  }
  libm_dtwc = dtmin;
  /* crlibm timing */
  dtmin=1<<30;
  for(j=0; j<10; j++) {
    TBX_GET_TICK(t1);
    result = testfun_crlibm(input);
    TBX_GET_TICK(t2);
    dt = TBX_TICK_RAW_DIFF(t1, t2); 
    if (dt<dtmin)  dtmin=dt;
  }
  crlibm_dtwc = dtmin;

#ifdef HAVE_MPFR_H
    /* mpfr timing */
    dtmin=1<<30;
    for(j=0; j<10; j++) {
      TBX_GET_TICK(t1);
      mpfr_set_d(mp_inpt, input, GMP_RNDN);
      testfun_mpfr(mp_res, mp_inpt, GMP_RNDN);
      result = mpfr_get_d1(mp_res);
      TBX_GET_TICK(t2);
      dt = TBX_TICK_RAW_DIFF(t1, t2); 
      if (dt<dtmin)  dtmin=dt;
    }
    mpfr_dtwc = dtmin;
#endif /*HAVE_MPFR_H*/

#ifdef HAVE_MATHLIB_H
    Original_Mode = Init_Lib(); 
    dtmin=1<<30;
    for(j=0; j<10; j++) {
      TBX_GET_TICK(t1);
      result = testfun_ibm(input);
      TBX_GET_TICK(t2);
      dt = TBX_TICK_RAW_DIFF(t1, t2); 
      if (dt<dtmin)  dtmin=dt;
    }
    Exit_Lib(Original_Mode);
    ibm_dtwc = dtmin;
#endif /*HAVE_MATHLIB_H*/


    /*************Normal output*************************/

    printf("\nLIBM\n");
    printf("Tmin = %lld ticks,\t Tmax = %lld ticks\t avg = %f\tworst case = %lld\n ",  
	   libm_dtmin, libm_dtmax,
	   (((double)libm_dtsum) / ((double) n)),
	   libm_dtwc
	   );
    
#ifdef HAVE_MPFR_H
    printf("\nMPFR\nTmin = %lld ticks,\t Tmax = %lld ticks\t avg = %f\tworst case = %lld\n",
	   mpfr_dtmin, mpfr_dtmax,
	   ((double)mpfr_dtsum) / ((double) n),
	   mpfr_dtwc
	 );
#endif /*HAVE_MPFR_H*/
    
#ifdef HAVE_MATHLIB_H
    printf("\nIBM\nTmin = %lld ticks,\t Tmax = %lld ticks\t avg = %f\tworst case = %lld\n",
	 ibm_dtmin, ibm_dtmax,
	   ((double)ibm_dtsum) / ((double) n),
	   ibm_dtwc
	   );
#endif /*HAVE_MATHLIB_H*/
    
    printf("\nCRLIBM\n");
    printf("Tmin = %lld ticks,\t Tmax = %lld ticks\t avg = %f\tworst case = %lld\n ",  
	   crlibm_dtmin, crlibm_dtmax,
	   ((double)crlibm_dtsum) / ((double) n),
	   crlibm_dtwc
	   );


    /******************* Latex output ****************/
    printf("\\multicolumn{4}{|c|}{Processor / system / compiler}   \\\\ \n \\hline");
    printf("\n                   \t & min time \t & max time \t & avg time \\\\ \n \\hline\n");
    printf(" \\texttt{libm}     \t & %lld    \t& %lld     \t& %10.0f \\\\ \n \\hline\n ",  
	   libm_dtmin,  libm_dtmax,
	   (((double)libm_dtsum) / ((double) n))
	   );
#ifdef HAVE_MPFR_H
    if (mpfr_dtwc > mpfr_dtmax) mpfr_dtmax=mpfr_dtwc;
    printf(" \\texttt{mpfr}     \t & %lld    \t& %lld     \t& %10.0f \\\\ \n \\hline\n ",  
	 mpfr_dtmin, mpfr_dtmax,
	 ((double)mpfr_dtsum) / ((double) n)
	 );
#endif /*HAVE_MPFR_H*/

#ifdef HAVE_MATHLIB_H
    if (ibm_dtwc > ibm_dtmax) ibm_dtmax=ibm_dtwc;
    printf(" \\texttt{libultim}  \t & %lld    \t& %lld     \t& %10.0f \\\\ \n \\hline\\hline\n ",  
	 ibm_dtmin, ibm_dtmax,
	 ((double)ibm_dtsum) / ((double) n)
	 );
#endif /*HAVE_MATHLIB_H*/

    if (crlibm_dtwc > crlibm_dtmax) crlibm_dtmax=crlibm_dtwc;
    printf("\\texttt{crlibm}     \t & %lld    \t& %lld     \t& %10.0f \\\\ \n \\hline\n ",  
	 crlibm_dtmin, crlibm_dtmax,
	 ((double)crlibm_dtsum) / ((double) n)
	   );


  /* release memory */
#ifdef HAVE_MPFR_H
  mpfr_clear(mp_inpt);
  mpfr_clear(mp_res);
#endif /*HAVE_MPFR_H*/

  return 0;
}

