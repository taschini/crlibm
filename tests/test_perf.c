
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


#define N1 4

#define DETAILED_REPORT 0

/* If set, the behaviour of the function with respect to cache memory
   will be tested*/
#define TEST_CACHE 0

/*
 * Rounding mode to test
 */




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







void fill_and_flush(int seed) {
  int i;
  srandom(seed);
  for (i=0; i<TESTSIZE; i++)
    inputs[i] = randfun();
}



void usage(char *fct_name){
  fprintf (stderr, "\n\n Performance test for crlibm \n");
  fprintf (stderr, "Usage: %s name (RN|RU|RD|RZ) iterations \n", fct_name);
  fprintf (stderr, " name          : name of function to test \n");
  fprintf (stderr, " (RN/RU/RD/RZ) : rounding mode, \n");
  fprintf (stderr, " iterations    : integer, seed for the random number generator \n");
  exit (1);
}





int main (int argc, char *argv[]) 
{ 
  
  int i, j, k, n;
  int counter;
  double input, result;
  char* rounding_mode;
  char* function_name;
  double worstcase;
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
mp_rnd_t mpfr_rnd_mode;
#endif
  short Original_Mode;





  if ((argc !=4)) usage(argv[0]);
  else{
    function_name = argv[1];
    rounding_mode = argv[2];
    sscanf(argv[3],"%d", &n);
    
    
    crlibm_init();

    test_init(/* pointers to returned value */
	       &randfun, 
	       &testfun_crlibm, 
	       &testfun_mpfr,
	       &testfun_libm,
	       &testfun_ibm,
	       &worstcase,
	       /* arguments */
	       function_name,
	       rounding_mode ) ;
    
#ifdef HAVE_MPFR_H  
    mpfr_init2(mp_res,  153);
    mpfr_init2(mp_inpt, 53);
    if      (strcmp(rounding_mode,"RU")==0) mpfr_rnd_mode = GMP_RNDU;
    else if (strcmp(rounding_mode,"RD")==0) mpfr_rnd_mode = GMP_RNDD;
    else if (strcmp(rounding_mode,"RZ")==0) mpfr_rnd_mode = GMP_RNDZ;
    else {
      mpfr_rnd_mode = GMP_RNDN; 
      rounding_mode="RN" ;
    }
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

  counter=0;

  /* take the min of N1 identical calls to leverage interruptions */
  /* As a consequence, the cache impact of these calls disappear...*/

  for(i=0; i< n; i++){ 
    input = randfun();

    /* libm timing */
    dtmin=1<<30;
    /* take the min of N1 consecutive calls */
    for(j=0; j<N1; j++) {
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
    /* take the min of N1 consecutive calls */
    for(j=0; j<N1; j++) {
      counter++;
      TBX_GET_TICK(t1);
#ifdef CRLIBM_TYPECPU_POWERPC
      counter+=49; for(k=0; k<50;k++)
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
    /* take the min of N1 consecutive calls */
    for(j=0; j<N1; j++) {
      TBX_GET_TICK(t1);
#ifdef CRLIBM_TYPECPU_POWERPC
      for(k=0; k<50;k++) {
#endif
      mpfr_set_d(mp_inpt, input, GMP_RNDN);
      testfun_mpfr(mp_res, mp_inpt, GMP_RNDN);
      result = mpfr_get_d(mp_res,mpfr_rnd_mode);
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
    /* take the min of N1 consecutive calls */
    for(j=0; j<N1; j++) {
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
	 crlibm_second_step_taken/N1, counter/N1 );
#endif


  /************  WORST CASE TESTS   *********************/
  /* worst case test */
  input = worst_case;
  /* libm timing */
  dtmin=1<<30;
  for(j=0; j<N1; j++) {
    TBX_GET_TICK(t1);
#ifdef CRLIBM_TYPECPU_POWERPC
      for(k=0; k<50;k++)
#endif
    result = testfun_libm(input);
    TBX_GET_TICK(t2);
    dt = TBX_TICK_RAW_DIFF(t1, t2); 
    if (dt<dtmin)  dtmin=dt;
  }
  libm_dtwc = dtmin;
  /* crlibm timing */
  dtmin=1<<30;
  for(j=0; j<N1; j++) {
    TBX_GET_TICK(t1);
#ifdef CRLIBM_TYPECPU_POWERPC
      for(k=0; k<50;k++)
#endif
    result = testfun_crlibm(input);
    TBX_GET_TICK(t2);
    dt = TBX_TICK_RAW_DIFF(t1, t2); 
    if (dt<dtmin)  dtmin=dt;
  }
  crlibm_dtwc = dtmin;

#ifdef HAVE_MPFR_H
    /* mpfr timing */
    dtmin=1<<30;
    for(j=0; j<N1; j++) {
      TBX_GET_TICK(t1);
#ifdef CRLIBM_TYPECPU_POWERPC
      for(k=0; k<50;k++){
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
    mpfr_dtwc = dtmin;
#endif /*HAVE_MPFR_H*/

#ifdef HAVE_MATHLIB_H
    Original_Mode = Init_Lib(); 
    dtmin=1<<30;
    for(j=0; j<N1; j++) {
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
}

