#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* Ziv library */
#include <MathLib.h>
#include <gmp.h>
#include <mpfr.h>
/* Our library */
#include "crlibm.h"

#define TIRET "\n========================================================\n"

/* To adapt to the tested precision */
#define FPPREC 53       /* 24       53         64*/

/* stolen from mpfr-impl.h */
#define MPFR_EXP(x) ((x)->_mpfr_exp)

#ifdef __i386
#include <fpu_control.h>
#ifndef __setfpucw
#define __setfpucw(cw) __asm__ ("fldcw %0" : : "m" (cw))
#endif /* ifndef __setfpucw */
#endif

#ifdef __alpha
#ifndef FP_RND_RN
#define FP_RND_RZ       0
#define FP_RND_RN       1
#define FP_RND_RP       2
#define FP_RND_RM       3
#endif
#define TONEAREST write_rnd(FP_RND_RN)
#define TOZERO    write_rnd(FP_RND_RZ)
#define TOINFP    write_rnd(FP_RND_RP)
#define TOINFM    write_rnd(FP_RND_RM)
#else /* ifdef __alpha */
#include <fenv.h>
#define TONEAREST fesetround(FE_TONEAREST)
#define TOZERO    fesetround(FE_TOWARDZERO)
#define TOINFP    fesetround(FE_UPWARD)
#define TOINFM    fesetround(FE_DOWNWARD)
#endif

#ifndef MAX_RND
#define MAX_RND 4
#endif

void test _PROTO ((char *, mp_exp_t, unsigned long, unsigned long));
void test2 _PROTO ((char *, mp_exp_t, mp_exp_t, unsigned long, unsigned long));
void testall _PROTO ((unsigned long, unsigned long));
void check_fp _PROTO ((void));

/* sets the machine rounding mode to the value rnd_mode */
double (*testfun_libm) () = NULL;
double (*testfun_ibm)  () = NULL;
double (*testfun_sn_cr)() = NULL;
double (*testfun_sd_cr)() = NULL;
double (*testfun_su_cr)() = NULL;
int    (*testfun_mpfr) () = NULL;

double MAX_ERR_NEAR = 0.0;
double MAX_ERR_DIR = 0.0;

void 
mpfr_set_machine_rnd_mode (mp_rnd_t rnd_mode)
{
  switch (rnd_mode) {
  case GMP_RNDN: TONEAREST; break;
  case GMP_RNDZ: TOZERO; break;
  case GMP_RNDU: TOINFP; break;
  case GMP_RNDD: TOINFM; break;
  default: fprintf(stderr, "invalid rounding mode\n"); exit(1);
  }
}

void
check_fp ()
{
  double x, y, c, d, dj;
  int j;

#ifdef __i386
  /* sets the precision to double */
  __setfpucw((_FPU_DEFAULT & (~_FPU_EXTENDED)) | _FPU_DOUBLE);
#endif
  x = DBL_MIN;
  if (2.0 * (x / 2.0) != x)
    fprintf (stderr, "Warning: no denormalized numbers\n");

  c = 1.46484375e-3;
  dj = 1.0;
  for (j=0; j<54; j++) dj *= 0.5;
  d = 0.75 + dj;
  d /= 1 << 9;
  if (c != d)
    {
      fprintf (stderr, "default seems to use extended precision\n");
      exit (1);
    }

  mpfr_set_machine_rnd_mode (GMP_RNDD);
  x = 2.0; /* don't write x = sqrt (2.0) and y = sqrt (2.0) otherwise the
              compiler may optimize too much */
  x = sqrt (x);
  mpfr_set_machine_rnd_mode (GMP_RNDU);
  y = 2.0;
  y = sqrt (y);
  if (x == y)
    {
      fprintf (stderr, "setting rounding modes does not work, "
               "you may have to add an option to the C compiler\n");
      fprintf (stderr, "   on Alpha, try:\n");
      fprintf (stderr, "   '-fprm d -ieee_with_inexact' with cc\n");
      fprintf (stderr, "   '-mfp-rounding-mode=d -mieee-with-inexact' with gcc\n");
      exit (1);
    }
}



void test_fct(double (*testfct)(), mp_exp_t e, mp_rnd_t rnd_mode, unsigned long N, unsigned long seed){
  double xd, yd, r, u;
  double umax, usum, xmax;
  mpfr_t x, y, z;
  mpfr_t xx, yy, xin;
  mp_exp_t expy;
  gmp_randstate_t state;
  unsigned long i, tot;


  mpfr_init2 (x, FPPREC);
  mpfr_init2 (y, FPPREC);
  mpfr_init2 (z, FPPREC);
  mpfr_init2 (xx, 4*FPPREC);
  mpfr_init2 (xin, 4*FPPREC);
  mpfr_init2 (yy, 4*FPPREC);

  tot = 0; usum = 0.0; umax = 0.0; xmax = 0.0;
  mpfr_set_machine_rnd_mode (GMP_RNDN); 
  /* reset the seed to test the same sequence of numbers */
  //mpfr_set_machine_rnd_mode (rnd_mode); 
  gmp_randinit (state, GMP_RAND_ALG_LC, 128);
  gmp_randseed_ui (state, seed);

  for (i=0; i<N; i++){
    mpfr_urandomb (x, state);
    MPFR_EXP(x) = e;
    
    /* Conversion from a mpfr to a fp number : x */
    xd = (double) mpfr_get_d1 (x);
    mpfr_set_d(x, xd, GMP_RNDN);
    
    testfun_mpfr (y, x, rnd_mode);
    /* Conversion from a mpfr to a fp number : y */
    yd = (double) mpfr_get_d1 (y);
    
    /* Set the input number with many extra bits */
    mpfr_set_d (xin, xd, GMP_RNDN); 
  
   
    if (testfct != NULL) {
      /* Tested library */
      r = testfct (xd);
      /* check for correct directed rounding */
      
      if (yd != r) {
	/* Computation with extra bits */
	testfun_mpfr (yy, xin, rnd_mode);
	expy = MPFR_EXP(yy);
	mpfr_set_d (xx, r, GMP_RNDN); /* exact */
	mpfr_sub (xx, xx, yy, GMP_RNDN);
	MPFR_EXP(xx) += FPPREC - expy;
	u = mpfr_get_d1 (xx);
	
	tot ++;
	usum += fabs(u);
	if (fabs(u) > fabs(umax))
	  {
	    umax = fabs(u);
	    xmax = xd;
	  }
      }
    }else {
      umax = 0; xmax = 0;
    }
  }

  printf ("(nb errors) / (max ulp diff) / (mean ulp diff):\n %lu/%f/%f/\n", 
	  tot, fabs(umax), usum/(double)N);
  printf ("%e ulp(s) =(%d bits) for x=%1.50e\n", 
	  umax, (umax!=0)? (int)ulog2(umax)+1:0, xmax);
  if (umax != 0){
    db_number db;

    db.d = testfct (xmax);
    printf("Lib  Result : %8x %8x \n", db.i[HI_ENDIAN], db.i[LO_ENDIAN]);
    
    mpfr_set_d(x, xmax, GMP_RNDN);
    testfun_mpfr (y, x, rnd_mode);
    db.d = (double) mpfr_get_d1 (y);
    printf("True Result : %8x %8x \n", db.i[HI_ENDIAN], db.i[LO_ENDIAN]);
  }

  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (z);
  mpfr_clear (xx);
  mpfr_clear (yy);

}


void test (char *foo, mp_exp_t e, unsigned long N, unsigned long seed)
{
 
   printf ("Testing function %s for exponent %ld.\n", foo, e);
   fflush (stdout);
  
  if (strcmp (foo, "exp") == 0)
    {
      testfun_libm = exp;
      testfun_ibm  = uexp;
      testfun_sn_cr= sn_exp;
      testfun_sd_cr= sd_exp;
      testfun_su_cr= su_exp;
      testfun_mpfr = mpfr_exp;
    }
  else if (strcmp (foo, "log") == 0)
    {
      testfun_libm = log;
      testfun_ibm  = ulog;
      testfun_sn_cr= sn_log;
      testfun_sd_cr= sd_log;
      testfun_su_cr= su_log;
      testfun_mpfr = mpfr_log;
    }
  else if (strcmp (foo, "log2") == 0)
    {
      testfun_libm = NULL;   /* Doesn't exist */
      testfun_ibm  = ulog2;
      testfun_sn_cr= sn_log2;
      testfun_sd_cr= sd_log2;
      testfun_su_cr= su_log2;
      testfun_mpfr = mpfr_log2;
    }
  else if (strcmp (foo, "log10") == 0)
    {
      testfun_libm = log10;
      testfun_ibm  = NULL;   /* Doesn't exist */
      testfun_sn_cr= sn_log10;
      testfun_sd_cr= sd_log10;
      testfun_su_cr= su_log10;
      testfun_mpfr = mpfr_log10;
    }
  else if (strcmp (foo, "sin") == 0)
    {
      testfun_libm = sin;
      testfun_ibm  = usin;
      testfun_sn_cr= sn_sin;
      testfun_sd_cr= sd_sin;
      testfun_su_cr= su_sin;
      testfun_mpfr = mpfr_sin;
    }
  else if (strcmp (foo, "cos") == 0)
    {
      testfun_libm = cos;
      testfun_ibm  = ucos;
      testfun_sn_cr= sn_cos;
      testfun_sd_cr= sd_cos;
      testfun_su_cr= su_cos;
      testfun_mpfr = mpfr_cos;
    }
  /*
  else if (strcmp (foo, "atan") == 0)
    {
      testfun_libm = atan;
      testfun_ibm  = uatan;
      testfun_sn_cr= sn_atan;
      testfun_sd_cr= sd_atan;
      testfun_su_cr= su_atan;
      testfun_mpfr = mpfr_atan;
    }
  */
  else
    {
      fprintf (stderr, "Unknown function: %s\n", foo);
      exit (1);
    }
  
  printf(TIRET);
  printf("LIBM To nearest: \n");
  test_fct(testfun_libm, e, GMP_RNDN, N, seed); 
  
  printf(TIRET);
  printf("IBM To nearest: \n");
  test_fct(testfun_ibm, e, GMP_RNDN, N, seed);
  
  printf(TIRET);
  printf("ENS To nearest: \n");
  test_fct(testfun_sn_cr, e, GMP_RNDN, N, seed);
  
  printf(TIRET);
  printf("ENS Toward -inf: \n");
  test_fct(testfun_sd_cr, e, GMP_RNDD, N, seed);

  printf(TIRET);
  printf("ENS Toward +inf: \n");
  test_fct(testfun_su_cr, e, GMP_RNDU, N, seed);
  

  fflush (stdout);
    
}


 
void
test2 (char *foo, mp_exp_t e, mp_exp_t f, unsigned long N, unsigned long seed)
{
  unsigned long i, wrong, tot;
  mpfr_t x, y, z, t;
  mpfr_t xx, yy, tt;
  mp_exp_t expy;
  double xd, td, yd, r;
  double u, usum, umax, xmax, tmax, max_err_near, max_err_dir;
  mp_rnd_t rnd;
  gmp_randstate_t state;

  printf ("Testing function %s for exponents %ld and %ld.\n", foo, e, f);
  fflush (stdout);

  if (strcmp (foo, "pow") == 0)
    {
      testfun_libm = pow;
      testfun_mpfr = mpfr_pow;
    }
  else
    {
      fprintf (stderr, "Unknown function: %s\n", foo);
      exit (1);
    }

  mpfr_init2 (x, FPPREC);
  mpfr_init2 (y, FPPREC);
  mpfr_init2 (z, FPPREC);
  mpfr_init2 (t, FPPREC);
  mpfr_init2 (xx, 106);
  mpfr_init2 (tt, 106);
  mpfr_init2 (yy, 106);

  max_err_near = 0.0;
  max_err_dir = 0.0;

  for (rnd=0; rnd<MAX_RND; rnd++)
    {
      printf ("   rounding mode %s:\n", mpfr_print_rnd_mode (rnd));
      fflush (stdout);
      mpfr_set_machine_rnd_mode (rnd);

      gmp_randinit (state, GMP_RAND_ALG_LC, 128);
      gmp_randseed_ui (state, seed);

      umax = 0.0;
      xmax = 0.0;
      tmax = 0.0;
      tot = 0;
      wrong = 0;
      for (i=0; i<N; i++)
	{
	  mpfr_urandomb (x, state);
	  MPFR_EXP(x) = e;
	  mpfr_urandomb (t, state);
	  MPFR_EXP(t) = f;
	  testfun_mpfr (y, x, t, rnd);
	  /* Conversion from a mpfr to a fp number : x */
	  xd = (double) mpfr_get_d1 (x);
	  mpfr_set_d(yy, (double)xd, rnd);
	  mpfr_sub(xx, x, yy, rnd);
	  xd += (double) mpfr_get_d1(xx);
	  /* Conversion from a mpfr to a fp number : t */
	  td = (double) mpfr_get_d1 (t);
	  mpfr_set_d(yy, (double)td, rnd);
	  mpfr_sub(xx, t, yy, rnd);
	  td += (double) mpfr_get_d1(xx);
	  /* Conversion from a mpfr to a fp number : y */
	  yd = (double) mpfr_get_d1 (y);
	  mpfr_set_d(yy, (double)yd, rnd);
	  mpfr_sub(xx, y, yy, rnd);
	  yd += (double) mpfr_get_d1(xx);

	  r = testfun_libm (xd, td);
	  
	  /* check for correct directed rounding */
	  if (yd != r) {
	    
	    mpfr_set_d (xx, xd, GMP_RNDN); /* exact */
	    mpfr_set_d (tt, td, GMP_RNDN); /* exact */
	    testfun_mpfr (yy, xx, tt, rnd);
	    expy = MPFR_EXP (yy);
	    mpfr_set_d (xx, r, GMP_RNDN); /* exact */
	    mpfr_sub (xx, xx, yy, GMP_RNDN);
	    MPFR_EXP(xx) += FPPREC - expy;
	    u = mpfr_get_d1 (xx);
	    
	    if (rnd != GMP_RNDN)
	      {
		mp_rnd_t rnd2;
		
		if (rnd == GMP_RNDZ)
		  rnd2 = (mpfr_cmp_ui (y, 0) > 0) ? GMP_RNDD : GMP_RNDU;
		else
		  rnd2 = rnd;
		if ((rnd2 == GMP_RNDU && u < 0.0) || (rnd2 == GMP_RNDD && u > 0.0))
		  {
		    wrong++;
		    if (wrong == 1)
		      {
			printf ("      wrong directed rounding for x=%1.20e "
				" t=%1.20e [%f]\n", xd, td, u);
			fflush (stdout);
		      }
		  }
	      }
	    tot ++;
	    usum += fabs(u);
	    if (fabs(u) > fabs(umax))
	      {
		umax = u;
		xmax = xd;
		tmax = td;
	      }
	  }
	}

      if (umax != 0.0)
	printf ("      %f ulp(s) for x=%1.20e t=%1.20e\n", umax, xmax, tmax);

      umax = fabs (umax);

      printf ("   nb errors/max ulp diff/mean ulp diff/wrong directed: %lu/%f/%f/%lu\n", 
	      tot, umax, usum/(double)N, wrong);
      fflush (stdout);

      if (rnd == GMP_RNDN)
	{
	  if (umax > max_err_near)
	    max_err_near = umax;
	}
      else
	{
	  if (umax > max_err_dir)
	    max_err_dir = umax;
	}
    }

  printf ("Maximal errors for %s: %f (nearest), %f (directed)\n", foo,
          max_err_near, max_err_dir);
  fflush (stdout);

  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (z);
  mpfr_clear (t);
  mpfr_clear (xx);
  mpfr_clear (tt);
  mpfr_clear (yy);

  if (max_err_near > MAX_ERR_NEAR)
    MAX_ERR_NEAR = max_err_near;

  if (max_err_dir > MAX_ERR_DIR)
    MAX_ERR_DIR = max_err_dir;
}

void
testall (unsigned long N, unsigned long seed)
{
  test ("exp",     0, N, seed);
  test ("exp",     4, N, seed);
  test ("exp",     8, N, seed);
  test ("exp",     9, N, seed);
  test ("log",     0, N, seed);
  test ("log",  1024, N, seed);
  test ("log2",    0, N, seed);
  test ("log2", 1024, N, seed);
  test ("log10",   0, N, seed);
  test ("log10",1024, N, seed);
  test ("sin",     0, N, seed);
  test ("sin",    10, N, seed); /* mpfr-2.0.1 is too slow for 1024 */
  test ("cos",     0, N, seed);
  test ("cos",    10, N, seed);
  test ("atan",    0, N, seed);
  test ("atan", 1024, N, seed);
  test2 ("pow", 0, 0, N, seed);
  test2 ("pow", 8, 7, N, seed);

  printf ("Maximal errors for all functions: %f (nearest), %f (directed)\n",
          MAX_ERR_NEAR, MAX_ERR_DIR);
}

int
main (int argc, char *argv[])
{
  mp_exp_t exponent, exp2 = 0;
  unsigned long N, nargs, seed = 0;

  check_fp ();

  if (argc == 1 || argc == 3)
    {
      fprintf (stderr, "Usage: mpfr-test [-seed s] N\n");
      fprintf (stderr, "Usage: mprf-test [-seed s] <function> <exponent> [N]\n");
      fprintf (stderr, "Usage: mpfr-test [-seed s] <function> <exp1> <exp2> [N]\n");
      exit (1);
    }

  if (strcmp (argv[1], "-seed") == 0)
    {
      seed = atoi (argv[2]);
      argc -= 2;
      argv += 2;
    }


  if (argc == 2)
    {
      N = atoi(argv[1]);
      testall (N, seed);
    }
  else
    {
      nargs = 1;
      if (strcmp (argv[1], "pow") == 0)
        {
          nargs = 2;
          exp2  = atoi (argv[3]);
        }
  
      exponent = atoi (argv[2]);
      N = (argc > 2 + nargs) ? atoi(argv[2 + nargs]) : 1;

      if (nargs == 1)
        test (argv[1], exponent, N, seed);
      else
        test2 (argv[1], exponent, exp2, N, seed);
    }
  
  return 0;
}
