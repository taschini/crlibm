/* timings for David's thesis */

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

#if EVAL_PERF==1  
/* defined in exp.c */
extern int crlibm_exp_taken;
/* counter of calls to the second step */
extern int crlibm_exp_second_step_taken; 
#endif

//#undef HAVE_MPFR_H
//#undef HAVE_MATHLIB_H


#define DETAILED_REPORT 0

/* TESTSIZE doubles should be enough to flush the cache */
#define TESTSIZE 200000 
double inputs[TESTSIZE];

void fill_and_flush(int seed) {
  int i;
  srandom(seed);
  for (i=0; i<TESTSIZE; i++)
    inputs[i] = rand_for_exp_normal();
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
#ifdef HAVE_MPFR_H  
  mpfr_t mp_res, mp_inpt; 
#endif
  short Original_Mode;

  /************  INITIALISATIONS  *********************/

  if (argc != 2) {
    printf("Usage: %s iterations \n", argv[0]);
    exit(1);
  }
 
  sscanf(argv[1],"%d", &n);


  crlibm_init();

#ifdef HAVE_MPFR_H
  mpfr_init2(mp_res,  53);
  mpfr_init2(mp_inpt, 53);
#endif

  crlibm_dtmin=1<<30; crlibm_dtmax=0; crlibm_dtsum=0;
  libm_dtmin=1<<30;   libm_dtmax=0;   libm_dtsum=0;
  mpfr_dtmin=1<<30;   mpfr_dtmax=0;   mpfr_dtsum=0;
  ibm_dtmin=1<<30;    ibm_dtmax=0;    ibm_dtsum=0;


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
	  result = exp(input);
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
	  result = exp_rn(input);
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
	  result = uexp(input);
	  TBX_GET_TICK(t2);
	  dt = TBX_TICK_RAW_DIFF(t1, t2); 
	  dtsum += dt;
	}
	if (dtsum < min_dtsum) min_dtsum=dtsum; 
      }
      printf("  %d loops: \t avg time = %f ticks\n",i, ((double)min_dtsum)/i);
    }
#endif

  /************  TESTS WITHOUT CACHES  *********************/
  srandom(n);

#if EVAL_PERF==1  
  crlibm_exp_taken=0; 
  crlibm_exp_second_step_taken=0; 
#endif
  /* take the min of 10 identical calls to leverage interruptions */
  /* As a consequence, the cache impact of these calls disappear...*/

  for(i=0; i< n; i++){ 
    input = rand_for_exp_normal();

    /* libm timing */
    dtmin=1<<30;
    for(j=0; j<10; j++) {
      TBX_GET_TICK(t1);
      result = exp(input);
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
    for(j=0; j<10; j++) {
      TBX_GET_TICK(t1);
      result = exp_rn(input);
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
    for(j=0; j<10; j++) {
      TBX_GET_TICK(t1);
      mpfr_set_d(mp_inpt, input, GMP_RNDN);
      mpfr_exp(mp_res, mp_inpt, GMP_RNDN);
      result = mpfr_get_d1(mp_res);
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
    for(j=0; j<10; j++) {
      TBX_GET_TICK(t1);
      result = uexp(input);
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



  /************  WORST CASE TESTS   *********************/
  /* worst case test */
  input = .75417527749959590085206221024712557043923055744016892276704311370849609375e-9;

  /* libm timing */
  dtmin=1<<30;
  for(j=0; j<10; j++) {
    TBX_GET_TICK(t1);
    result = exp(input);
    TBX_GET_TICK(t2);
    dt = TBX_TICK_RAW_DIFF(t1, t2); 
    if (dt<dtmin)  dtmin=dt;
  }
  libm_dtwc = dtmin;
  /* crlibm timing */
  dtmin=1<<30;
  for(j=0; j<10; j++) {
    TBX_GET_TICK(t1);
    result = exp_rn(input);
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
      mpfr_exp(mp_res, mp_inpt, GMP_RNDN);
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
      result = uexp(input);
      TBX_GET_TICK(t2);
      dt = TBX_TICK_RAW_DIFF(t1, t2); 
      if (dt<dtmin)  dtmin=dt;
    }
    Exit_Lib(Original_Mode);
    ibm_dtwc = dtmin;
#endif /*HAVE_MATHLIB_H*/

  printf("\nLIBM\n");
  printf("Tmin = %lld ticks,\t Tmax = %lld ticks\t avg = %f\tworst case = %lld\n ",  
	 libm_dtmin, libm_dtmax,
	 (((double)libm_dtsum) / ((double) n)),
	 libm_dtwc
	 );

  printf("\nCRLIBM\n");
  printf("Tmin = %lld ticks,\t Tmax = %lld ticks\t avg = %f\tworst case = %lld\n",
	 crlibm_dtmin, crlibm_dtmax,
	 ((double)crlibm_dtsum) / ((double) n),
	 crlibm_dtwc
	 );

#if EVAL_PERF==1  
  printf("Second step taken %d times out of %d\n",
	 crlibm_exp_second_step_taken, crlibm_exp_taken );
#endif

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
  
  /* release memory */
#ifdef HAVE_MPFR_H
  mpfr_clear(mp_inpt);
  mpfr_clear(mp_res);
#endif /*HAVE_MPFR_H*/

  return 0;
}

