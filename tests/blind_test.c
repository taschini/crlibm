#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "crlibm.h"
#include "crlibm_private.h"
#include "test_common.h"

void usage(char *fct_name){
  fprintf (stderr, "\nUsage: %s [-v] file\n", fct_name);
  exit (EXIT_SUCCESS);
}



char* skip_comments(FILE* f, char* line) {
  char* r; 
  int i;
  do {
    r=fgets(line, 200, f);
  }
    /* Look if the first char was a #, a newline or a blank character */
  while((line[0]=='#' || line[0]=='\n' || line[0]==' ' || line[0]=='\t' || line[0]==0) && r!=NULL);

  /* look for a comment in the line, and remove it */
  i=0;
  while(line[i]!=0 && line[i]!='#') i++;
  line[i]=0;
  return r;
}


int main (int argc, char *argv[]) 
{ 
  int verbose=0;
  int failures=0;
  char* option;
  char* filename;
  char rounding_mode[10];
  char function_name[10];
  char line[400];
  char* r;
  int count=0;
  double worstcase;
  db_number input, output, expected, res_mpfr;
  double (*testfun_crlibm)() = NULL;
  double (*unused)() = NULL;
  int (*mpfr_fun)() = NULL;
  mpfr_t mp_res, mp_inpt; 
  mp_rnd_t mpfr_rnd_mode;

  FILE* f;

  if ((argc != 2) && argc!=3) usage(argv[0]);
  if (argc == 3) {
    option = argv[1];
    filename=argv[2];
    if(!(option[0]=='-' && option[1]=='v')) usage(argv[0]);
    else verbose=1;
  }
  else
    filename=argv[1];


  f=fopen(filename, "r");
  if(f==NULL) {
    fprintf(stderr, "%s: problem opening %s, exiting\n", argv[0], filename);
    exit(EXIT_FAILURE);
  }


  crlibm_init();
   
  fflush(stdout); /* To help debugging */

  skip_comments(f, line);
  sscanf(line, "%s", function_name);

  if(verbose)  printf("Testing function: %s\n", function_name);

  r=skip_comments(f, line);
  while(r!=0) { 
    sscanf(line, "%s %x %x %x %x\n", 
	   rounding_mode, 
	   &input.i[HI], &input.i[LO],
	   &expected.i[HI], &expected.i[LO] );
      /* Centralized test initialization function */
    test_init(
	      &unused, &unused, 
	      &testfun_crlibm, 
	      &mpfr_fun, &unused, &unused, &unused, &worstcase,
	      function_name,
	      rounding_mode ) ;
    
    output.d = testfun_crlibm(input.d);
    count++;

    if(verbose){
      printf("Input:      %08x %08x  (%0.50e)\n", input.i[HI], input.i[LO], input.d ); 
      printf("    Output: %08x %08x  (%0.50e)", output.i[HI], output.i[LO], output.d ); 
      if( output.l==expected.l)
	printf("   ...  OK \n");
      else
	printf(" \n");

    }
    /* The following tests is a little bit contrived because of NaN */
    if(    ((expected.d != expected.d) && (output.d == output.d)) /* expected NaN, got non-NaN */
        || ((expected.d == expected.d) && (output.l != expected.l))    ) { /* or expected non-NaN, got something different */
      failures ++;
      printf("ERROR for %s with rounding %s\n", function_name, rounding_mode);
      printf("       Input: %08x %08x  (%0.50e)\n", input.i[HI], input.i[LO], input.d ); 
      printf("      Output: %08x %08x  (%0.50e)\n", output.i[HI], output.i[LO], output.d ); 
      printf("    Expected: %08x %08x  (%0.50e)\n", expected.i[HI], expected.i[LO], expected.d ); 

      if      ((strcmp(rounding_mode,"RU")==0) || (strcmp(rounding_mode,"P")==0)) mpfr_rnd_mode = GMP_RNDU;
      else if ((strcmp(rounding_mode,"RD")==0) || (strcmp(rounding_mode,"M")==0)) mpfr_rnd_mode = GMP_RNDD;
      else if ((strcmp(rounding_mode,"RZ")==0) || (strcmp(rounding_mode,"Z")==0)) mpfr_rnd_mode = GMP_RNDZ;
      else if ((strcmp(rounding_mode,"RN")==0) || (strcmp(rounding_mode,"N")==0)) mpfr_rnd_mode = GMP_RNDN;
      else {
	fprintf(stderr, "Unknown rounding mode: %s, exiting\n", rounding_mode);
	exit(EXIT_FAILURE);
      }

      mpfr_init2(mp_res,  160);
      mpfr_init2(mp_inpt, 53);
      mpfr_set_d(mp_inpt, input.d, GMP_RNDN);
      mpfr_fun(mp_res, mp_inpt, mpfr_rnd_mode);
      res_mpfr.d = mpfr_get_d(mp_res, mpfr_rnd_mode);

      printf(" MPFR result: %08x %08x  (%0.50e)\n", res_mpfr.i[HI], res_mpfr.i[LO], res_mpfr.d ); 

    }
    fflush(stdout); /* To help debugging */
    
    r=skip_comments(f, line);
  } 
  printf("Test completed for %s, %d failures in %d tests\n", function_name, failures, count);
  return failures;
  
}



