#include <stdio.h>

#include <gmp.h>
#include <mpfr.h>


/*
 * This program builds the tables used in the first step of the exponential
 * Auteur : David Defour in a hurry             Date : 15/06/2003  
 *
 * build command:
  gcc -I/usr/local/include/ -L/usr/local/lib/  exp_firststep_table.c -o exp_firststep_table -lgmp  -lmpfr 
 */


/* don't forget to check that */
#define LITTLE_ENDIAN

#ifdef LITTLE_ENDIAN
#define HI(x) *(1+(int*)&x)
#define LO(x) *(int*)&x
#else
#define HI(x) *(int*)&x
#define LO(x) *(1+(int*)&x)
#endif


/* This is the number of bits considered to address the table. It can
   be changed by the command-line second argument, but should stay to
   8 for the exp as distributed */

int PREC=8;



#define CC_MPFR(d1, d2, a, res, rest, f, A, B){ \
    mpfr_exp(res, a, GMP_RNDN); \
    d1 = mpfr_get_d(res, GMP_RNDN);\
    mpfr_set_d(rest, d1, GMP_RNDN);\
    mpfr_sub(res, res, rest, GMP_RNDN);\
    d2 = mpfr_get_d(res, GMP_RNDN); \
    fprintf(f,"\t {{{0x%.8x, 0x%.8x}}, {{0x%.8x, 0x%.8x}}}", A(d1), B(d1), A(d2), B(d2));}


void compute(FILE *f){
  mpfr_t a, b, step, res, rest;
  double d, d1, d2;
  int i, end;

  /* Memory allocation */
  mpfr_init2(a,   300);
  mpfr_init2(b,   300);
  mpfr_init2(res, 300);
  mpfr_init2(rest,300);
  mpfr_init2(step,300);


  /* Set step to 2^-PREC */
  mpfr_set_ui(step, 1, GMP_RNDN);
  mpfr_div_2ui(step, step, PREC, GMP_RNDN);

  /* Compute bound of the interval [a, b] */
  mpfr_set_ui(res, 2, GMP_RNDN);
  mpfr_log(a, res, GMP_RNDN); 
  mpfr_div_2ui(a, a, 1, GMP_RNDN);

  mpfr_mul_2ui(a, a, PREC, GMP_RNDN);
  mpfr_round(a, a);
  d = mpfr_get_d(a, GMP_RNDN);
  mpfr_div_2ui(a, a, PREC, GMP_RNDN);

  mpfr_neg(a, a, GMP_RNDN);
  mpfr_set(b, a, GMP_RNDN);
  
  end = (2 * d) + 1; /* +1 because this range is centered in 0 and '0' need to be taken into account */

  fprintf(f,"static const scs_db_number tab_exp[%d][2] = { \n",end);
  fprintf(f,"#ifdef WORDS_BIGENDIAN\n");

  for(i=0; i<(end-1) ; i++){
    CC_MPFR(d1, d2, a, res, rest, f, HI, LO);
    fprintf(f,", \n");
    mpfr_add(a, a, step, GMP_RNDN);
  }
  CC_MPFR(d1, d2, a, res, rest, f, HI, LO);

  fprintf(f,"\n#else \n");
  mpfr_set(a, b, GMP_RNDN);
  for(i=0; i<(end-1) ; i++){
    CC_MPFR(d1, d2, a, res, rest, f, LO, HI);
    fprintf(f,", \n");
    mpfr_add(a, a, step, GMP_RNDN);
  }
  CC_MPFR(d1, d2, a, res, rest, f, LO, HI);

  fprintf(f,"\n#endif \n}; \n\n");


  mpfr_clear(a);  
  mpfr_clear(b);
  mpfr_clear(res); 
  mpfr_clear(rest); 
  mpfr_clear(step);
} 



int main (int argc, char *argv[]) 
{ 
  FILE *f;

  if (argc == 3) 
    PREC = atoi(argv[2]);
  
  f = fopen (argv[1], "w");
  if (f== NULL){
    fprintf(stderr," Cannot open file %s  \n", argv[1]);    exit(1);
  }

  compute(f);

  fclose(f);
  return 0;
}

