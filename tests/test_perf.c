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

#ifdef HAVE_LIBMCR_H
#include <libmcr.h>
#endif

#ifdef HAVE_MPFR_H
# include <gmp.h>
# include <mpfr.h>
# ifdef MPFR_VERSION
#  if MPFR_VERSION < MPFR_VERSION_NUM(2,2,0)
#   define mpfr_subnormalize(mp_res, inexact, mpfr_rnd_mode) 
#  endif
# else
#  define mpfr_get_version() "<2.1.0"
# endif
#endif


#define N1 20

#define TIMING_ITER 100

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

#ifdef HAVE_MPFR_H  
/* The rounding mode for mpfr function */
mp_rnd_t mpfr_rnd_mode;
#endif

/* Global variable for time stamping */
static unsigned long long tbx_time;

/* Unused random number generator*/
double (*randfun_soaktest) () = NULL;
/* The random number generator*/
double (*randfun)       () = NULL;
/* The function we test */
double (*testfun_crlibm)() = NULL;
/* The function we trust */
int    (*testfun_mpfr)  () = NULL;
/* The function to show off against for accuracy  */
double (*testfun_libm)  () = NULL;
/* The function to show off against for performance */
double (*testfun_libultim)   () = NULL;
/* The last competitor  */
double (*testfun_libmcr)   () = NULL;


/* TESTSIZE doubles should be enough to flush the cache */
#define TESTSIZE 200000
 
#if TEST_CACHE
static double inputs[TESTSIZE];
static double inputs2[TESTSIZE];
#endif /* TEST_CACHE */


/* indicate the number of argument taken by the function */
static int nbarg;          






static void usage(char *fct_name){
  /*fprintf (stderr, "\n%s: Performance test for crlibm and other mathematical libraries\n", fct_name);*/
  fprintf (stderr, "\nUsage: %s function (RN|RU|RD|RZ) iterations \n", fct_name);
  fprintf (stderr, " function      : name of function to test \n");
  fprintf (stderr, " (RN|RU|RD|RZ) : rounding mode, \n");
  fprintf (stderr, " iterations    : number of iterations, also seed for the random number generator \n");
  exit (EXIT_FAILURE);
}


#if TEST_CACHE

static void fill_and_flush(int seed, int args) {
  int i;
  srandom(seed);
  for (i=0; i<TESTSIZE; i++){
    if (args==1) {
      inputs[i] = randfun();
      inputs2[i] = randfun();
    } else {
      inputs[i] = (*((double (*)(double *))randfun))(&(inputs2[i]));
    }
  }
}

static void test_with_cache(const char *name, double (*testfun)(), int n, int args){
  int i, j, k;
  double i1, i2, rd;
  tbx_tick_t   t1, t2; 
  unsigned long long dt, min_dtsum, dtsum;

  if(testfun!=NULL) { /* test if some functions are missing */
    printf("\n%s\n",name);
    for(i=1; i<=10000; i*=10) { /* i=1,10,100...*/
      min_dtsum=1<<30; 
      for(k=0;k<10;k++) { /* do the whole test 10 times and take the min */
	fill_and_flush(n,args);
	dtsum=0;
	for (j=0; j<i; j++) {
	  i1 = inputs[i];
	  i2 = inputs2[i];

	  if (nbarg==1){
	    TBX_GET_TICK(t1);
	    rd = testfun(i1);
	    TBX_GET_TICK(t2);
	  }else{
	    TBX_GET_TICK(t1);
	    rd = testfun(i1, i2);
	    TBX_GET_TICK(t2);
	  }

	  dt = TBX_TICK_RAW_DIFF(t1, t2)-tbx_time; 
	  dtsum += dt;
	}
	if (dtsum < min_dtsum) min_dtsum=dtsum; 
      }
      printf("  %d loops: \t avg time = %f ticks\n",i, ((double)min_dtsum)/i);
    }
  }
}
#endif /* TEST_CACHE */

static void test_without_cache(const char *name, 
			double (*testfun)(), 
			double i1,
			double i2,
			unsigned long long *lib_dtmin,
			unsigned long long *lib_dtmax,
			unsigned long long *lib_dtsum,
			double *lib_dtmini1,
			double *lib_dtmini2,
			double *lib_dtmaxi1,
			double *lib_dtmaxi2,
			int func_type){
  double result;
  unsigned long long dt, dtmin;
  tbx_tick_t   t1, t2; 
  int j;
#ifdef TIMING_USES_GETTIMEOFDAY
  int k;
#endif 
#ifdef HAVE_MPFR_H  
  mpfr_t mp_res, mp_inpt;
  mpfr_t mp_inpt2; /* For the pow function */
  int inexact;
  double dtmini1, dtmini2;

  mpfr_init2(mp_res,  53);
  mpfr_init2(mp_inpt, 53);
  mpfr_init2(mp_inpt2, 53);
#endif

  if(testfun!=NULL) { /* test if some functions are missing */
    dtmin=1<<30;
    /* take the min of N1 consecutive calls */
    for(j=0; j<N1; j++) {

      if (func_type == 0){              /* func_type = normal function */
	if (nbarg==1){
	  TBX_GET_TICK(t1);
#ifdef TIMING_USES_GETTIMEOFDAY /* use inaccurate timer, do many loops */
	  for(k=0; k<TIMING_ITER;k++)
#endif
	    result = testfun(i1);
	  TBX_GET_TICK(t2);
	}else{
	  TBX_GET_TICK(t1);
#ifdef TIMING_USES_GETTIMEOFDAY /* use inaccurate timer, do many loops */
	  for(k=0; k<TIMING_ITER;k++)
#endif
	    result = testfun(i1,i2);
	  TBX_GET_TICK(t2);	  
	}
      }else{                             /* func_type = MPFR function  */
#ifdef   HAVE_MPFR_H 
	if (nbarg==1){
	  TBX_GET_TICK(t1);
#ifdef TIMING_USES_GETTIMEOFDAY /* use inaccurate timer, do many loops */
	  for(k=0; k<TIMING_ITER;k++){
#endif    
	    mpfr_set_d(mp_inpt, i1, GMP_RNDN);
	    inexact = testfun(mp_res, mp_inpt, mpfr_rnd_mode);
            mpfr_subnormalize (mp_res, inexact, mpfr_rnd_mode);
	    result = mpfr_get_d(mp_res, GMP_RNDN);
#ifdef TIMING_USES_GETTIMEOFDAY /* use inaccurate timer, do many loops */
	  }
#endif
	  TBX_GET_TICK(t2);
	}else{
	  TBX_GET_TICK(t1);
#ifdef TIMING_USES_GETTIMEOFDAY /* use inaccurate timer, do many loops */
	  for(k=0; k<TIMING_ITER;k++){
#endif    
	    mpfr_set_d(mp_inpt, i1, GMP_RNDN);
	    mpfr_set_d(mp_inpt2, i2, GMP_RNDN);
	    inexact = testfun(mp_res, mp_inpt, mp_inpt2, mpfr_rnd_mode);
            mpfr_subnormalize (mp_res, inexact, mpfr_rnd_mode);
	    result = mpfr_get_d(mp_res, GMP_RNDN);
#ifdef TIMING_USES_GETTIMEOFDAY /* use inaccurate timer, do many loops */
	  }
#endif 
	  TBX_GET_TICK(t2);
	}
#endif /*HAVE_MPFR_H*/
      }

      dt = TBX_TICK_RAW_DIFF(t1, t2)-tbx_time;
      if (dt<dtmin)  {
	dtmin=dt;
	dtmini1 = i1;
	dtmini2 = i2;
      }
    }
    *lib_dtsum+=dtmin;
    if (dtmin<*lib_dtmin)  {
      *lib_dtmin=dtmin;
      *lib_dtmini1 = dtmini1;
      *lib_dtmini2 = dtmini2;
    }
    if (dtmin>*lib_dtmax)  {
      *lib_dtmax=dtmin;
      *lib_dtmaxi1 = dtmini1;
      *lib_dtmaxi2 = dtmini2;
    }
#if      DETAILED_REPORT
    printf("\n input=%1.15e\tT%s=%lld", i1, name, dtmin);
#endif /*DETAILED_REPORT*/
  }

  /* release memory */
#ifdef   HAVE_MPFR_H
  mpfr_clear(mp_inpt2);
  mpfr_clear(mp_inpt);
  mpfr_clear(mp_res);
#endif /*HAVE_MPFR_H*/

}




static void test_worst_case(double (*testfun)(), 
		     double i1, 
		     double i2, 
		     unsigned long long *lib_dtwc, 
		     int func_type){


  double res;
  tbx_tick_t   t1, t2; 
  unsigned long long dtmin, dt;
  int j;
#ifdef TIMING_USES_GETTIMEOFDAY
  int k;
#endif 
#ifdef HAVE_MPFR_H  
  mpfr_t mp_res, mp_inpt;
  mpfr_t mp_inpt2; /* For the pow function */
  int inexact;

  mpfr_init2(mp_res,  53);
  mpfr_init2(mp_inpt, 53);
  mpfr_init2(mp_inpt2, 53);
#endif

  if(testfun!=NULL) { /* test if some functions are missing  */
    dtmin=1<<30;
    for(j=0; j<N1; j++) {
      if (func_type == 0){    /* func_type = normal function */
	if (nbarg==1){
	  TBX_GET_TICK(t1);
#ifdef TIMING_USES_GETTIMEOFDAY
	  for(k=0; k<TIMING_ITER;k++)
#endif
	    res = testfun(i1);
	  TBX_GET_TICK(t2);
	}else{
	  TBX_GET_TICK(t1);
#ifdef TIMING_USES_GETTIMEOFDAY
	  for(k=0; k<TIMING_ITER;k++)
#endif
	    res = testfun(i1,i2);
	  TBX_GET_TICK(t2);
	}
      }else{                   /* func_type = MPFR function  */
#ifdef HAVE_MPFR_H 
	if (nbarg==1){
	  TBX_GET_TICK(t1);
#ifdef TIMING_USES_GETTIMEOFDAY
	  for(k=0; k<TIMING_ITER;k++){
#endif
	    mpfr_set_d (mp_inpt, i1, GMP_RNDN);
	    inexact = testfun (mp_res, mp_inpt, mpfr_rnd_mode);
            mpfr_subnormalize (mp_res, inexact, mpfr_rnd_mode);
	    res = mpfr_get_d (mp_res, GMP_RNDN);
#ifdef TIMING_USES_GETTIMEOFDAY
	  }
#endif
	  TBX_GET_TICK(t2);
	}else{
	  TBX_GET_TICK(t1);
#ifdef TIMING_USES_GETTIMEOFDAY
	  for(k=0; k<TIMING_ITER;k++){
#endif
	    mpfr_set_d (mp_inpt, i1, GMP_RNDN);
	    mpfr_set_d (mp_inpt2, i2, GMP_RNDN);
	    inexact = testfun (mp_res, mp_inpt, mp_inpt2, mpfr_rnd_mode);
            mpfr_subnormalize (mp_res, inexact, mpfr_rnd_mode);
	    res = mpfr_get_d (mp_res, GMP_RNDN);
#ifdef TIMING_USES_GETTIMEOFDAY
	  }
#endif
	  TBX_GET_TICK(t2);
	}
#endif /*HAVE_MPFR_H*/
      }
      dt = TBX_TICK_RAW_DIFF(t1, t2)-tbx_time; 
      if (dt<dtmin)  dtmin=dt;
    }
    *lib_dtwc = dtmin;
  }

  /* release memory */
#ifdef   HAVE_MPFR_H
  mpfr_clear(mp_inpt2);
  mpfr_clear(mp_inpt);
  mpfr_clear(mp_res);
#endif /*HAVE_MPFR_H*/
}



static void normal_output(const char *name,
		   double (*testfun)(),
		   unsigned long long lib_dtmin,
		   unsigned long long lib_dtmax,
		   unsigned long long lib_dtsum,
		   unsigned long long lib_dtwc,
		   double lib_dtmini1, 
		   double lib_dtmini2, 
		   double lib_dtmaxi1,
		   double lib_dtmaxi2,  
		   int n, int args){
  db_number dtmini1, dtmini2, dtmaxi1, dtmaxi2;

  dtmini1.d = lib_dtmini1;
  dtmini2.d = lib_dtmini2;
  dtmaxi1.d = lib_dtmaxi1;
  dtmaxi2.d = lib_dtmaxi2;

  if(testfun!=NULL) { /* some functions are missing in libultim (cosh, ...  */
    if (args == 1) {
      printf("\n%s\nTmin = %lld ticks (0x%08x%08x),\t Tmax = %lld ticks (0x%08x%08x),\t avg = %f\tT worst case = %lld\n",
	     name, lib_dtmin, dtmini1.i[HI], dtmini1.i[LO], lib_dtmax, dtmaxi1.i[HI], dtmaxi1.i[LO], 
	     ((double)lib_dtsum) / ((double) n), lib_dtwc);
    } else {
      printf("\n%s\nTmin = %lld ticks (0x%08x%08x 0x%08x%08x),\t Tmax = %lld ticks (0x%08x%08x 0x%08x%08x),\t avg = %f\tT worst case = %lld\n",
	     name, lib_dtmin, dtmini1.i[HI], dtmini1.i[LO], dtmini2.i[HI], dtmini2.i[LO], 
	     lib_dtmax, dtmaxi1.i[HI], dtmaxi1.i[LO], dtmaxi2.i[HI], dtmaxi2.i[LO], 
	     ((double)lib_dtsum) / ((double) n), lib_dtwc);
    }
  }
}

static void latex_output(const char *name,
		  double (*testfun)(),
		  unsigned long long lib_dtmin,
		  unsigned long long lib_dtmax,
		  unsigned long long lib_dtsum,
		  unsigned long long lib_dtwc,
		  int n){
    if(testfun!=NULL) { /* some functions are missing in libultim (cosh, ...  */
      if (lib_dtwc > lib_dtmax) lib_dtmax=lib_dtwc;
      printf(" %s  \t& %lld    \t& %10.0f   \t& %lld      \\\\ \n \\hline\n",  
	     name, lib_dtmin, ((double)lib_dtsum) / ((double) n), lib_dtmax);
    }
}


void generate_pow_exact_case(double *x, double *y, int doSubnormals) {
  int ok,E,F,i,maskN,n,expoBound,maxM,maskM,valM,a,b,EP,EPP,maxA;
  db_number xdb, ydb, zdb;
  double m;

  /* Decide first if we generate a positive or negative y */

  if (rand_int() & 1) {

    /* Positive y */

    ok = 0;
    
    while (!ok) {

      /* First generate an appropriate exponent F for y 
	 The exponent must be between -5 and 5
      */
      do {
	F = (rand_int() & 0xf) - 8;
      } while ((F < -5) || (F > -5));
      
      /* Generate E 
	 If F is positive, E may be a random integer between 0 and 15 
	 If F is negative, E must be divisible by 2^(-F), and less than 
	 927
	 
	 So we generate an integer between 0 and 15 and multiply by 2^(-F) if 
	 F is negative
	 
      */
      
      E = rand_int() & 0x0f;
      
      EP = E;

      if (F < 0) {
	if (E != 0) { 
	  while ((E & 1) == 0) E >>= 1; 
	}
	E <<= -F;
      }
      
      EPP = E;

      
      /* Generate n 
	 If F is negative, n may be in the range 1..32
	 If F is positive, n may be in the range 1..(32 / 2^F)
	 
	 We will generate n - 1 in a range 0..31 or 0..((32 / 2^F) - 1)
	 We will than bring n to the next odd number lower than it.
	 
	 If F = 0, n may not be equal 1, we take 3 instead.
	 
	 We generate first an appropriate mask.
	 
      */
      
      if (F < 0) {
	maskN = 0x1f;
      } else {
	maskN = 0;
	for (i=0;i<F;i++) maskN = maskN * 2 + 1;
      }
      
      n = (rand_int() & maskN) + 1;
      
      if ((n & 1) != 1) n--;
      
      if ((F == 0) && (n == 1)) n = 3;

      
      /* Generate now ydb.d = 2^F * n */
      
      ydb.d = n;
      ydb.i[HI] += F << 20;
      
      /* Generate now m 
	 
         If F is negative generate j and take m = j^(2^-F)
	 If F is positive generate m 
	 
	 j must be such that j^(2^-F) is less than 2^53 - 1
	 m must be such that m^(2^F * n) is less than 2^53 - 1
	 
	 Call 2^-F and 2^F * n respectively expoBound. expoBound is less or
	 equal to 32.
	 Out of expoBound, compute the maximum value for m or j
	 If F negative take m = j^(2^(-F)) else take the value for m
	 
      */
      
      if (F < 0) {
	expoBound = 1 << -F;
      } else {
	expoBound = n << F;
      }
      
      switch (expoBound) {
      case 2: 
	maxM = 94906265; maskM = 134217727;
	break;
      case 3:
	maxM = 208063; maskM = 262143;
	break;
      case 4:
	maxM = 9741; maskM = 16383;
	break;
      case 5:
	maxM = 1552; maskM = 2047;
	break;
      case 6:
	maxM = 456; maskM = 511;
	break;
      case 7:
	maxM = 190; maskM = 255;
	break;
      case 8:
	maxM = 98; maskM = 127;
	break;
      case 9:
	maxM = 59; maskM = 63;
	break;
      case 10:
	maxM = 39; maskM = 63;
	break;
      case 11:
	maxM = 28; maskM = 31;
	break;
      case 12:
	maxM = 21; maskM = 31;
	break;
      case 13:
	maxM = 16; maskM = 15;
	break;
      case 14:
	maxM = 13; maskM = 15;
	break;
      case 15:
	maxM = 11; maskM = 15;
	break;
      case 16:
	maxM = 9; maskM = 15;
	break;
      case 17:
	maxM = 8; maskM = 7;
	break;
      case 18:
	maxM = 7; maskM = 7;
	break;
      case 19:
      case 20:
	maxM = 6; maskM = 7;
	break;
      case 21:
      case 22:
	maxM = 5; maskM = 7;
	break;
      case 23:
      case 24:
      case 25:
      case 26:
	maxM = 4; maskM = 3;
	break;
      case 27:
      case 28:
      case 29:
      case 30:
      case 31:
      case 32:
	maxM = 3; maskM = 3;
      default:
	maxM = 1; maskM = 0;
      }
      
      /* Compute now an odd value valM such that 1 <= valM <= maxM */
      
      if (maxM == 1) {
	valM = 1;
      } else {
	while ((valM = ((rand_int() & maskM) + 1)) > maxM);
      }
      
      if ((valM & 1) == 0) valM--;
      
      if (F >= 0) {
	m = valM;
      } else {
	m = valM;
	for (i=-F;i>0;i--) m *= m;
      }
      
      /* Generate now xdb.d = 2^E * m */
      
      xdb.d = m;
      xdb.i[HI] += E << 20;
      
      if (doSubnormals) {
	if ((xdb.d != 0.0)  && (xdb.d != 1.0) && (ydb.d != 1.0) && (ydb.d != 0.0)) ok = 1;
      } else {
    
	/* Test if we produce a subnormal */
     
	zdb.d = pow_rn(xdb.d,ydb.d);

	/* If we produce a subnormal, we restart the procedure */

	if ((zdb.i[HI] & 0xfff00000) != 0) {
	  /* We are not subnormal and thus okay */
	  if ((xdb.d != 0.0)  && (xdb.d != 1.0) && (ydb.d != 1.0) && (ydb.d != 0.0)) ok = 1;
	}
      }
    }

  } else {
    
    /* Negative y */

    ok = 0;
    
    while (!ok) {

      /* Generate an exponent b such that 

	 -32 <= b <= 31, b != 0

	 Take x = 2^b

	 abs(y) must be bounded by 32
      */
    

      b = (rand_int() & 0x3f) - 32;

      if (b == 0) b = -1;
             
      
      /* Compute a bound for a such that a * b is bounded by 32 in absolute value */
      
      if (b < 0) b = -b;
      
      maxA = 1;
      while (b < 16) {
	maxA *= 2;
	b *= 2;
      } 
      
      do {
	a = (rand_int() & 0x1f);
      } while (a > maxA);
      
      a *= b;
      
      if (a > 0) a = -a;

      if (a * b >= 1024) {
	b >>= 1;
      }
     
      xdb.i[HI] = (b + 1023) << 20;
      xdb.i[LO] = 0;

      ydb.d = a;

      if (doSubnormals) {
	if ((xdb.d != 0.0)  && (xdb.d != 1.0) && (ydb.d != 1.0) && (ydb.d != 0.0)) ok = 1;
      } else {
    
	/* Test if we produce a subnormal */
     
	zdb.d = pow_rn(xdb.d,ydb.d);

	/* If we produce a subnormal, we restart the procedure */

	if ((zdb.i[HI] & 0xfff00000) != 0) {
	  /* We are not subnormal and thus okay */
	  if ((xdb.d != 0.0)  && (xdb.d != 1.0) && (ydb.d != 1.0) && (ydb.d != 0.0)) ok = 1;
	}
      }
    }
  }
    
  /* Assign now x and y */

  *x = xdb.d;
  *y = ydb.d;
    
}



int main (int argc, char *argv[]){ 
  int i, j, n;
  double i1, i2;
  char* rounding_mode;
  char* function_name;
  double worstcase;
  tbx_tick_t   t1, t2; 
  unsigned long long 
    dt,
    libm_dtmin, libm_dtmax, libm_dtsum, libm_dtwc,
    crlibm_dtmin, crlibm_dtmax, crlibm_dtsum, crlibm_dtwc,
    mpfr_dtmin, mpfr_dtmax, mpfr_dtsum, mpfr_dtwc,
    libultim_dtmin, libultim_dtmax, libultim_dtsum, libultim_dtwc,
    libmcr_dtmin, libmcr_dtmax, libmcr_dtsum, libmcr_dtwc,
    libm_dtecSample, crlibm_dtecSample, mpfr_dtecSample, libmcr_dtecSample,
    libm_dtecMin, crlibm_dtecMin, mpfr_dtecMin, libmcr_dtecMin,
    libm_dtecMax, crlibm_dtecMax, mpfr_dtecMax, libmcr_dtecMax,
    libm_dtecSum, crlibm_dtecSum, mpfr_dtecSum, libmcr_dtecSum;
  double libm_dtmini1,libm_dtmini2,crlibm_dtmini1,crlibm_dtmini2,libmcr_dtmini1,libmcr_dtmini2,
    mpfr_dtmini1, mpfr_dtmini2, libm_dtmaxi1, libm_dtmaxi2,crlibm_dtmaxi1,crlibm_dtmaxi2,libmcr_dtmaxi1,libmcr_dtmaxi2,
    mpfr_dtmaxi1, mpfr_dtmaxi2, libultim_dtmini1,libultim_dtmini2, libultim_dtmaxi1,libultim_dtmaxi2;
  double libm_dtecMinX, crlibm_dtecMinX, mpfr_dtecMinX, libmcr_dtecMinX;
  double libm_dtecMinY, crlibm_dtecMinY, mpfr_dtecMinY, libmcr_dtecMinY;
  double libm_dtecMaxX, crlibm_dtecMaxX, mpfr_dtecMaxX, libmcr_dtecMaxX;
  double libm_dtecMaxY, crlibm_dtecMaxY, mpfr_dtecMaxY, libmcr_dtecMaxY;


#ifdef   HAVE_MATHLIB_H
  short Original_Mode;
#endif


  if ((argc !=4)) usage(argv[0]);

  function_name = argv[1];
  rounding_mode = argv[2];
  sscanf(argv[3],"%d", &n);
    

#ifdef HAVE_MPFR_H  
  if      (strcmp(rounding_mode,"RU")==0) mpfr_rnd_mode = GMP_RNDU;
  else if (strcmp(rounding_mode,"RD")==0) mpfr_rnd_mode = GMP_RNDD;
  else if (strcmp(rounding_mode,"RZ")==0) mpfr_rnd_mode = GMP_RNDZ;
  else {
    mpfr_rnd_mode = GMP_RNDN; 
    rounding_mode="RN" ;
  }
#endif
    
  if (strcmp(function_name,"pow")==0) nbarg=2;
  else nbarg=1;

  crlibm_init();

  test_init(/* pointers to returned value */
	    &randfun, 
	    &randfun_soaktest, /* unused here */ 
	    &testfun_crlibm, 
	    &testfun_mpfr,
	    &testfun_libultim,
	    &testfun_libmcr,
	    &testfun_libm,
	    &worstcase,
	    /* arguments */
	    function_name,
	    rounding_mode ) ;
  


  crlibm_dtmin=1<<30; crlibm_dtmax=0; crlibm_dtsum=0; crlibm_dtwc=0;
  libm_dtmin=1<<30;   libm_dtmax=0;   libm_dtsum=0; libm_dtwc=0;
  libultim_dtmin=1<<30;    libultim_dtmax=0;    libultim_dtsum=0; libultim_dtwc=0;
  libmcr_dtmin=1<<30;    libmcr_dtmax=0;    libmcr_dtsum=0; libmcr_dtwc=0;
  mpfr_dtmin=1<<30;   mpfr_dtmax=0;   mpfr_dtsum=0; mpfr_dtwc=0;

  
  /* take the min of N1 consecutive calls */
  tbx_time=1<<30;
  for(j=0; j<1000*N1; j++) {
    TBX_GET_TICK(t1);
    TBX_GET_TICK(t2);
    dt = TBX_TICK_RAW_DIFF(t1, t2);
    if(dt<tbx_time) tbx_time = dt;
  }
#if HAVE_MPFR_H
  printf ("GMP version %s MPFR version %s ",
          gmp_version, mpfr_get_version ());
#endif
  printf("tbx_time=%llu\n", tbx_time);

#if TEST_CACHE
  /************  TESTS WITH CACHES  *********************/
  /* First tests in real conditions, where cache considerations
     matter */

  /* libm */
  printf("TEST WITH CACHE CONSIDERATION \n");
  test_with_cache("LIBM", testfun_libm, n, nbarg);
  test_with_cache("CRLIBM", testfun_crlibm, n, nbarg);
#ifdef HAVE_MATHLIB_H
  Original_Mode = Init_Lib();
  test_with_cache("IBM", testfun_libultim, n, nbarg);
  Exit_Lib(Original_Mode);
#endif
#ifdef HAVE_LIBMCR_H
  test_with_cache("SUN", testfun_libmcr, n, nbarg);
#endif

#endif /* TEST_CACHE*/

  /************  TESTS WITHOUT CACHES  *******************/
  srandom(n);
#if EVAL_PERF==1  
  crlibm_second_step_taken=0; 
#endif

  /* take the min of N1 identical calls to leverage interruptions */
  /* As a consequence, the cache impact of these calls disappear...*/
  for(i=0; i< n; i++){ 
    if (nbarg==1) {
      i1 = randfun();
      i2 = randfun();
    } else {
      i1 = (*((double (*)(double *))randfun))(&i2);
    }
    
    test_without_cache("libm", testfun_libm, i1, i2, &libm_dtmin, &libm_dtmax, &libm_dtsum, &libm_dtmini1, &libm_dtmini2, &libm_dtmaxi1, &libm_dtmaxi2, 0);
    test_without_cache("crlibm", testfun_crlibm, i1, i2, &crlibm_dtmin, &crlibm_dtmax, &crlibm_dtsum, &crlibm_dtmini1, &crlibm_dtmini2, &crlibm_dtmaxi1, &crlibm_dtmaxi2, 0);
#ifdef   HAVE_MATHLIB_H
    Original_Mode = Init_Lib(); 
    test_without_cache("ultim", testfun_libultim, i1, i2, &libultim_dtmin, &libultim_dtmax, &libultim_dtsum, &libultim_dtmini1, &libultim_dtmini2, &libultim_dtmaxi1, &libultim_dtmaxi2, 0);
    Exit_Lib(Original_Mode);
#endif /*HAVE_MATHLIB_H*/
#ifdef   HAVE_LIBMCR_H
    test_without_cache("libmcr", testfun_libmcr, i1, i2, &libmcr_dtmin, &libmcr_dtmax, &libmcr_dtsum, &libmcr_dtmini1, &libmcr_dtmini2, &libmcr_dtmaxi1, &libmcr_dtmaxi2, 0);
#endif /*HAVE_LIBMCR_H*/
#ifdef   HAVE_MPFR_H
    test_without_cache("mpfr", (double(*)()) testfun_mpfr, i1, i2, &mpfr_dtmin, &mpfr_dtmax, &mpfr_dtsum, &mpfr_dtmini1, &mpfr_dtmini2, &mpfr_dtmaxi1, &mpfr_dtmaxi2, 1);
#endif /*HAVE_MPFR_H*/
  } 

#if EVAL_PERF==1  
#ifdef TIMING_USES_GETTIMEOFDAY /* use inaccurate timer, do many loops */
	 printf("\nCRLIBM : Second step taken %d times out of %d\n",
		crlibm_second_step_taken/(N1 * TIMING_ITER), n );
#else
	 printf("\nCRLIBM : Second step taken %d times out of %d\n",
		crlibm_second_step_taken/N1, n );
#endif

#endif


  /************  WORST CASE TESTS   *********************/
  /* worst case test */
  i1 = worstcase;
  i2 = 0.156250000000000000000000000000000000000000000000000000e1; /* TODO when we have worst cases for power...*/

  test_worst_case(testfun_libm, i1, i2, &libm_dtwc, 0);
  test_worst_case(testfun_crlibm, i1, i2, &crlibm_dtwc, 0);
#ifdef   HAVE_MPFR_H
  test_worst_case((double(*)())testfun_mpfr, i1, i2, &mpfr_dtwc, 1);
#endif /*HAVE_MPFR_H*/
#ifdef   HAVE_MATHLIB_H
  Original_Mode = Init_Lib(); 
  test_worst_case(testfun_libultim, i1, i2, &libultim_dtwc, 0);
  Exit_Lib(Original_Mode);
#endif /*HAVE_MATHLIB_H*/
#ifdef   HAVE_LIBMCR_H
  test_worst_case(testfun_libmcr, i1, i2, &libmcr_dtwc, 0);
#endif /*HAVE_LIBMCR_H*/


  /************* Exact cases of power ****************/

  if (strcmp(function_name,"pow")==0) {
    /* Generate n exact cases of power, take the timings */

    printf("Power exact cases:\n");

    libm_dtecMin = 4294967295;
    crlibm_dtecMin = 4294967295;
    mpfr_dtecMin = 4294967295;
    libmcr_dtecMin = 4294967295;

    libm_dtecMax = 0;
    crlibm_dtecMax = 0;
    mpfr_dtecMax = 0;
    libmcr_dtecMax = 0;

    libm_dtecSum = 0;
    crlibm_dtecSum = 0;
    mpfr_dtecSum = 0;
    libmcr_dtecSum = 0;

    libm_dtecMaxX = 0.0;
    crlibm_dtecMaxX = 0.0;
    mpfr_dtecMaxX = 0.0;
    libmcr_dtecMaxX = 0.0;

    libm_dtecMinX = 0.0;
    crlibm_dtecMinX = 0.0;
    mpfr_dtecMinX = 0.0;
    libmcr_dtecMinX = 0.0;


    for (i=0;i<2* n;i++) {
      generate_pow_exact_case(&i1,&i2,1);

      test_worst_case(testfun_libm, i1, i2, &libm_dtecSample, 0);
      if (libm_dtecSample > libm_dtecMax) {
	libm_dtecMax = libm_dtecSample;
	libm_dtecMaxX = i1; libm_dtecMaxY = i2;
      }
      if (libm_dtecSample < libm_dtecMin) {
	libm_dtecMin = libm_dtecSample;
	libm_dtecMinX = i1; libm_dtecMinY = i2;
      }
      libm_dtecSum += libm_dtecSample;
      test_worst_case(testfun_crlibm, i1, i2, &crlibm_dtecSample, 0);
      if (crlibm_dtecSample > crlibm_dtecMax) {
	crlibm_dtecMax = crlibm_dtecSample;
	crlibm_dtecMaxX = i1; crlibm_dtecMaxY = i2;
      }
      if (crlibm_dtecSample < crlibm_dtecMin) {
	crlibm_dtecMin = crlibm_dtecSample;
	crlibm_dtecMinX = i1; crlibm_dtecMinY = i2;
      }
      crlibm_dtecSum += crlibm_dtecSample;
#ifdef   HAVE_MPFR_H
      test_worst_case((double(*)())testfun_mpfr, i1, i2, &mpfr_dtecSample, 1);
      if (mpfr_dtecSample > mpfr_dtecMax) {
	mpfr_dtecMax = mpfr_dtecSample;
	mpfr_dtecMaxX = i1; mpfr_dtecMaxY = i2;
      }
      if (mpfr_dtecSample < mpfr_dtecMin) {
	mpfr_dtecMin = mpfr_dtecSample;
	mpfr_dtecMinX = i1; mpfr_dtecMinY = i2;
      }
      mpfr_dtecSum += mpfr_dtecSample;
#endif /*HAVE_MPFR_H*/
#ifdef   HAVE_LIBMCR_H
      test_worst_case(testfun_libmcr, i1, i2, &libmcr_dtecSample, 0);
      if (libmcr_dtecSample > libmcr_dtecMax) {
	libmcr_dtecMax = libmcr_dtecSample;
	libmcr_dtecMaxX = i1; libmcr_dtecMaxY = i2;
      }
      if (libmcr_dtecSample < libmcr_dtecMin) {
	libmcr_dtecMin = libmcr_dtecSample;
	libmcr_dtecMinX = i1; libmcr_dtecMinY = i2;
      }

      libmcr_dtecSum += libmcr_dtecSample;
#endif /*HAVE_LIBMCR_H*/
 
    }

    printf("LIBM: min = %lld, avg = %f, max = %lld\n",libm_dtecMin,((double) libm_dtecSum)/((double) (2 * n)),libm_dtecMax);
    printf("Minimum value on:\n");
    printHexa("x",libm_dtecMinX);
    printHexa("y",libm_dtecMinY);
    printf("Maximum value on:\n");
    printHexa("x",libm_dtecMaxX);
    printHexa("y",libm_dtecMaxY);
    printf("CRLIBM: min = %lld, avg = %f, max = %lld\n",crlibm_dtecMin,((double) crlibm_dtecSum)/((double) (2 * n)),crlibm_dtecMax);
    printf("Minimum value on:\n");
    printHexa("x",crlibm_dtecMinX);
    printHexa("y",crlibm_dtecMinY);
    printf("Maximum value on:\n");
    printHexa("x",crlibm_dtecMaxX);
    printHexa("y",crlibm_dtecMaxY);
#ifdef   HAVE_MPFR_H
    printf("MPFR: min = %lld, avg = %f, max = %lld\n",mpfr_dtecMin,((double) mpfr_dtecSum)/((double) (2 * n)),mpfr_dtecMax);
    printf("Minimum value on:\n");
    printHexa("x",mpfr_dtecMinX);
    printHexa("y",mpfr_dtecMinY);
    printf("Maximum value on:\n");
    printHexa("x",mpfr_dtecMaxX);
    printHexa("y",mpfr_dtecMaxY);
#endif /*HAVE_MPFR_H*/
#ifdef   HAVE_LIBMCR_H
    printf("LIBMCR: min = %lld, avg = %f, max = %lld\n",libmcr_dtecMin,((double) libmcr_dtecSum)/((double) (2 * n)),libmcr_dtecMax);
    printf("Minimum value on:\n");
    printHexa("x",libmcr_dtecMinX);
    printHexa("y",libmcr_dtecMinY);
    printf("Maximum value on:\n");
    printHexa("x",libmcr_dtecMaxX);
    printHexa("y",libmcr_dtecMaxY);
#endif /*HAVE_LIBMCR_H*/

    printf("Power exact cases excluding subnormal results:\n");

    libm_dtecMin = 4294967295;
    crlibm_dtecMin = 4294967295;
    mpfr_dtecMin = 4294967295;
    libmcr_dtecMin = 4294967295;

    libm_dtecMax = 0;
    crlibm_dtecMax = 0;
    mpfr_dtecMax = 0;
    libmcr_dtecMax = 0;

    libm_dtecSum = 0;
    crlibm_dtecSum = 0;
    mpfr_dtecSum = 0;
    libmcr_dtecSum = 0;

    libm_dtecMaxX = 0.0;
    crlibm_dtecMaxX = 0.0;
    mpfr_dtecMaxX = 0.0;
    libmcr_dtecMaxX = 0.0;

    libm_dtecMinX = 0.0;
    crlibm_dtecMinX = 0.0;
    mpfr_dtecMinX = 0.0;
    libmcr_dtecMinX = 0.0;


    for (i=0;i<2 * n;i++) {
      generate_pow_exact_case(&i1,&i2,0);

      test_worst_case(testfun_libm, i1, i2, &libm_dtecSample, 0);
      if (libm_dtecSample > libm_dtecMax) {
	libm_dtecMax = libm_dtecSample;
	libm_dtecMaxX = i1; libm_dtecMaxY = i2;
      }
      if (libm_dtecSample < libm_dtecMin) {
	libm_dtecMin = libm_dtecSample;
	libm_dtecMinX = i1; libm_dtecMinY = i2;
      }
      libm_dtecSum += libm_dtecSample;
      test_worst_case(testfun_crlibm, i1, i2, &crlibm_dtecSample, 0);
      if (crlibm_dtecSample > crlibm_dtecMax) {
	crlibm_dtecMax = crlibm_dtecSample;
	crlibm_dtecMaxX = i1; crlibm_dtecMaxY = i2;
      }
      if (crlibm_dtecSample < crlibm_dtecMin) {
	crlibm_dtecMin = crlibm_dtecSample;
	crlibm_dtecMinX = i1; crlibm_dtecMinY = i2;
      }
      crlibm_dtecSum += crlibm_dtecSample;
#ifdef   HAVE_MPFR_H
      test_worst_case((double(*)())testfun_mpfr, i1, i2, &mpfr_dtecSample, 1);
      if (mpfr_dtecSample > mpfr_dtecMax) {
	mpfr_dtecMax = mpfr_dtecSample;
	mpfr_dtecMaxX = i1; mpfr_dtecMaxY = i2;
      }
      if (mpfr_dtecSample < mpfr_dtecMin) {
	mpfr_dtecMin = mpfr_dtecSample;
	mpfr_dtecMinX = i1; mpfr_dtecMinY = i2;
      }
      mpfr_dtecSum += mpfr_dtecSample;
#endif /*HAVE_MPFR_H*/
#ifdef   HAVE_LIBMCR_H
      test_worst_case(testfun_libmcr, i1, i2, &libmcr_dtecSample, 0);
      if (libmcr_dtecSample > libmcr_dtecMax) {
	libmcr_dtecMax = libmcr_dtecSample;
	libmcr_dtecMaxX = i1; libmcr_dtecMaxY = i2;
      }
      if (libmcr_dtecSample < libmcr_dtecMin) {
	libmcr_dtecMin = libmcr_dtecSample;
	libmcr_dtecMinX = i1; libmcr_dtecMinY = i2;
      }

      libmcr_dtecSum += libmcr_dtecSample;
#endif /*HAVE_LIBMCR_H*/
 
    }

    printf("LIBM: min = %lld, avg = %f, max = %lld\n",libm_dtecMin,((double) libm_dtecSum)/((double) (2 * n)),libm_dtecMax);
    printf("Minimum value on:\n");
    printHexa("x",libm_dtecMinX);
    printHexa("y",libm_dtecMinY);
    printf("Maximum value on:\n");
    printHexa("x",libm_dtecMaxX);
    printHexa("y",libm_dtecMaxY);
    printf("CRLIBM: min = %lld, avg = %f, max = %lld\n",crlibm_dtecMin,((double) crlibm_dtecSum)/((double) (2 * n)),crlibm_dtecMax);
    printf("Minimum value on:\n");
    printHexa("x",crlibm_dtecMinX);
    printHexa("y",crlibm_dtecMinY);
    printf("Maximum value on:\n");
    printHexa("x",crlibm_dtecMaxX);
    printHexa("y",crlibm_dtecMaxY);
#ifdef   HAVE_MPFR_H
    printf("MPFR: min = %lld, avg = %f, max = %lld\n",mpfr_dtecMin,((double) mpfr_dtecSum)/((double) (2 * n)),mpfr_dtecMax);
    printf("Minimum value on:\n");
    printHexa("x",mpfr_dtecMinX);
    printHexa("y",mpfr_dtecMinY);
    printf("Maximum value on:\n");
    printHexa("x",mpfr_dtecMaxX);
    printHexa("y",mpfr_dtecMaxY);
#endif /*HAVE_MPFR_H*/
#ifdef   HAVE_LIBMCR_H
    printf("LIBMCR: min = %lld, avg = %f, max = %lld\n",libmcr_dtecMin,((double) libmcr_dtecSum)/((double) (2 * n)),libmcr_dtecMax);
    printf("Minimum value on:\n");
    printHexa("x",libmcr_dtecMinX);
    printHexa("y",libmcr_dtecMinY);
    printf("Maximum value on:\n");
    printHexa("x",libmcr_dtecMaxX);
    printHexa("y",libmcr_dtecMaxY);
#endif /*HAVE_LIBMCR_H*/

  }

    /*************Normal output*************************/
  normal_output("LIBM", testfun_libm, libm_dtmin, libm_dtmax, libm_dtsum, libm_dtwc, libm_dtmini1, libm_dtmini2, libm_dtmaxi1, libm_dtmaxi2, n, nbarg);
#ifdef   HAVE_MPFR_H
  normal_output("MPFR", (double(*)())testfun_mpfr, mpfr_dtmin, mpfr_dtmax, mpfr_dtsum, mpfr_dtwc, mpfr_dtmini1, mpfr_dtmini2, mpfr_dtmaxi1, mpfr_dtmaxi2, n, nbarg);
#endif /*HAVE_MPFR_H*/
#ifdef   HAVE_MATHLIB_H
  normal_output("IBM", testfun_libultim, libultim_dtmin, libultim_dtmax, libultim_dtsum, libultim_dtwc, libultim_dtmini1, libultim_dtmini2, libultim_dtmaxi1, libultim_dtmaxi2, n, nbarg);
#endif /*HAVE_MATHLIB_H*/
#ifdef   HAVE_LIBMCR_H
  normal_output("SUN", testfun_libmcr, libmcr_dtmin, libmcr_dtmax, libmcr_dtsum, libmcr_dtwc, libmcr_dtmini1, libmcr_dtmini2, libmcr_dtmaxi1, libmcr_dtmaxi2, n, nbarg);
#endif /*HAVE_LIBMCR_H*/
  normal_output("CRLIBM", testfun_crlibm, crlibm_dtmin, crlibm_dtmax, crlibm_dtsum, crlibm_dtwc, crlibm_dtmini1, crlibm_dtmini2, crlibm_dtmaxi1, crlibm_dtmaxi2, n, nbarg);


  /******************* Latex output ****************/
  printf("\\multicolumn{4}{|c|}{Processor / system / compiler}   \\\\ \n \\hline");
  printf("\n                             & min time \t & avg time \t& max time \t  \\\\ \n \\hline\n");
  latex_output("default \\texttt{libm}  ", testfun_libm, libm_dtmin, libm_dtmax, libm_dtsum, libm_dtwc, n);
#ifdef   HAVE_MPFR_H
  latex_output("MPFR                   ", (double(*)())testfun_mpfr, mpfr_dtmin, mpfr_dtmax, mpfr_dtsum, mpfr_dtwc, n);
#endif /*HAVE_MPFR_H*/
#ifdef   HAVE_MATHLIB_H
  latex_output("IBM's \\texttt{libultim}", testfun_libultim, libultim_dtmin, libultim_dtmax, libultim_dtsum, libultim_dtwc, n);
#endif /*HAVE_MATHLIB_H*/
#ifdef   HAVE_LIBMCR_H
  latex_output("Sun's \\texttt{libmcr}  ", testfun_libmcr, libmcr_dtmin, libmcr_dtmax, libmcr_dtsum, libmcr_dtwc, n);
#endif /*HAVE_LIBMCR_H*/
  latex_output("\\texttt{crlibm}        ", testfun_crlibm, crlibm_dtmin, crlibm_dtmax, crlibm_dtsum, crlibm_dtwc, n);

  return 0;
}


