/*
 * Function to compute the logarithm base 10 with fully exact rounding
 *
 * Author : Defour David  (David.Defour@ens-lyon.fr)
 *
 * Date of creation : 27/06/2002   
 * Last Modified    : 27/06/2002
 */
#include "log10.h"

/*
 *  1) First reduction: exponent extraction      
 *         E  
 *  x = 2^   .(1+f)    with  0 <= f < 1
 *
 *  log(x) = E.log10(2) + log10(1+f) where:
 *     - log10(2)   is tabulated
 *     - log10(1+f) need to be evalute 
 *  
 *
 *  2) Avoiding accuracy problem when E=-1 by testing
 *   
 *    if (1+f >= sqrt(2)) then 
 *        1+f = (1+f)/2;  E = E+1; 
 *    and,
 *        log10(x) = (E+1).log10(2) + log10((1+f)/2)
 *
 *    so now:      sqrt(2)/2 <= (1+f) < sqrt(2)
 *
 *
 *  3) Second reduction: tabular reduction
 *                   -4  
 *   wi = 1 + i. 2^
 *                                         1  
 *   log10(1+f) = log10(wi) + log10 ( 1 + --- . (1 + f - wi) ) 
 *                                        wi 
 *
 *   then |(1+f-wi)/wi| <= 2^-5 if we use rounded to nearest.
 *
 *  4) Computation:
 *   a) Table lookup of: 
 *        - ti    = log10(wi)
 *        - inv_wi = 1/(wi)
 *   b) Polynomial evaluation of:
 *        - P(R) ~ log10(1 + R), where R = (1+f-wi) * inv_wi 
 *
 *                 -5 
 *   with  |R| < 2^
 *
 *
 *  5) Reconstruction:
 *   log10(x) = E.log10(2) + t_i + P(R)
 *
 *
 *   Note 1:
 *	To guarantee log10(10^n)=n, where 10^n is normal, the rounding 
 *	mode must set to Round-to-Nearest.
 *
 *   Special cases:
 *	log10(x) is NaN with signal if x < 0; 
 *	log10(+INF) is +INF with no signal; log10(0) is -INF with signal;
 *	log10(NaN) is that NaN with no signal;
 *	log10(10^N) = N  
 *
 */
#define SQRT_2 1.4142135623730950489e0 


/*************************************************************
 *************************************************************
 *               ROUNDED  TO NEAREST
 *************************************************************
 *************************************************************/
double log10_rn(double x)
{
  scs_t R, res1;
  scs_t sc_ln2_r10_times_E;
  scs_ptr inv_wi, ti;
 
  db_number nb, nb2, wi, resd;
  int i, E=0;

  nb.d = x;
  /* Filter cases */
  if (nb.i[HI_ENDIAN] < 0x00100000){        /* x < 2^(-1022)    */
    if (((nb.i[HI_ENDIAN] & 0x7fffffff)|nb.i[LO_ENDIAN])==0)
      return 1.0/0.0;                       /* log(+/-0) = -Inf */
    if (nb.i[HI_ENDIAN] < 0) 
      return (x-x)/0;                       /* log(-x) = Nan    */

    /* Subnormal number */
    E    -= (SCS_NB_BITS*2); /* keep in mind that x is a subnormal number */ 
    nb.d *=SCS_RADIX_TWO_DOUBLE;  /* make x as normal number     */         
    /* We may just want add 2 to the scs number.index */
    /* may be .... we will see */
  }
  if (nb.i[HI_ENDIAN] >= 0x7ff00000)
    return x+x;                             /* Inf or Nan       */

  /* find n, nb.d such that sqrt(2)/2 < nb.d < sqrt(2) */
  E += (nb.i[HI_ENDIAN]>>20)-1023;
  nb.i[HI_ENDIAN] =  (nb.i[HI_ENDIAN] & 0x000fffff) | 0x3ff00000;
  if (nb.d > SQRT_2){
    nb.d *= 0.5;
    E++;
  }


  /* to normalize nb.d and round to nearest      */
  /* +((2^4 - trunc(sqrt(2)/2) *2^4 )*2 + 1)/2^5 */ 
  nb2.d = nb.d + norm_number.d; 
  i = (nb2.i[HI_ENDIAN] & 0x000fffff);
  i = i >> 16; /* 0<= i <=11 */
  
  wi.d = (11+i)*(double)0.6250e-1;

  /* (1+f-w_i) */
  nb.d -= wi.d; 
  
  /* Table reduction */
  ti     = table_ti_ptr[i]; 
  inv_wi = table_inv_wi_ptr[i];

 
  /* R = (1+f-w_i)/w_i */
  scs_set_d(R, nb.d);
  scs_mul(R, R, inv_wi);
 

  /* sc_ln2_r10_times_E = E*log10(2)  */
  scs_set(sc_ln2_r10_times_E, sc_ln2_r10_ptr);
  
  if (E >= 0){
    scs_mul_ui(sc_ln2_r10_times_E, (unsigned int) E);
  }else{
    scs_mul_ui(sc_ln2_r10_times_E, (unsigned int) -E);
    sc_ln2_r10_times_E->sign = -1;
  }


  /*
   * Polynomial evaluation of log10(1 + R) with an error less than 2^(-130)
   */
  scs_mul(res1, constant_poly_ptr[0], R);
  for(i=1; i<20; i++){
    scs_add(res1, constant_poly_ptr[i], res1);
    scs_mul(res1, res1, R);
  }
  scs_add(res1, res1, ti);
  scs_add(res1, res1, sc_ln2_r10_times_E);  


  scs_get_d(&resd.d, res1);  
  
  return resd.d;
}









/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  -INFINITY
 *************************************************************
 *************************************************************/
double log10_rd(double x)
{
  scs_t R, res1;
  scs_t sc_ln2_r10_times_E;
  scs_ptr inv_wi, ti;

  db_number nb, nb2, wi, resd;
  int i, E=0;

  nb.d = x;
  /* Filter cases */
  if (nb.i[HI_ENDIAN] < 0x00100000){        /* x < 2^(-1022)    */
    if (((nb.i[HI_ENDIAN] & 0x7fffffff)|nb.i[LO_ENDIAN])==0)
      return 1.0/0.0;                       /* log(+/-0) = -Inf */
    if (nb.i[HI_ENDIAN] < 0) 
      return (x-x)/0;                       /* log(-x) = Nan    */

    /* Subnormal number */
    E    -= (SCS_NB_BITS*2); /* keep in mind that x is a subnormal number */ 
    nb.d *=SCS_RADIX_TWO_DOUBLE;  /* make x as normal number     */         
    /* We may just want add 2 to the scs number.index */
    /* may be .... we will see */
  }
  if (nb.i[HI_ENDIAN] >= 0x7ff00000)
    return x+x;                             /* Inf or Nan       */

  /* find n, nb.d such that sqrt(2)/2 < nb.d < sqrt(2) */
  E += (nb.i[HI_ENDIAN]>>20)-1023;
  nb.i[HI_ENDIAN] =  (nb.i[HI_ENDIAN] & 0x000fffff) | 0x3ff00000;
  if (nb.d > SQRT_2){
    nb.d *= 0.5;
    E++;
  }


  /* to normalize nb.d and round to nearest      */
  /* +((2^4 - trunc(sqrt(2)/2) *2^4 )*2 + 1)/2^5 */ 
  nb2.d = nb.d + norm_number.d; 
  i = (nb2.i[HI_ENDIAN] & 0x000fffff);
  i = i >> 16; /* 0<= i <=11 */
  
  wi.d = (11+i)*(double)0.6250e-1;

  /* (1+f-w_i) */
  nb.d -= wi.d; 
  
  /* Table reduction */
  ti     = table_ti_ptr[i]; 
  inv_wi = table_inv_wi_ptr[i];

 
  /* R = (1+f-w_i)/w_i */
  scs_set_d(R, nb.d);
  scs_mul(R, R, inv_wi);
 

  /* sc_ln2_r10_times_E = E*log10(2)  */
  scs_set(sc_ln2_r10_times_E, sc_ln2_r10_ptr);
  
  if (E >= 0){
    scs_mul_ui(sc_ln2_r10_times_E, (unsigned int) E);
  }else{
    scs_mul_ui(sc_ln2_r10_times_E, (unsigned int) -E);
    sc_ln2_r10_times_E->sign = -1;
  }


  /*
   * Polynomial evaluation of log10(1 + R) with an error less than 2^(-130)
   */
  scs_mul(res1, constant_poly_ptr[0], R);
  for(i=1; i<20; i++){
    scs_add(res1, constant_poly_ptr[i], res1);
    scs_mul(res1, res1, R);
  }
  scs_add(res1, res1, ti);
  scs_add(res1, res1, sc_ln2_r10_times_E);  


  scs_get_d_minf(&resd.d, res1);  
  
  return resd.d;
}



/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  +INFINITY
 *************************************************************
 *************************************************************/
double log10_ru(double x)
{
  scs_t R, res1;
  scs_t sc_ln2_r10_times_E;
  scs_ptr inv_wi, ti;

  db_number nb, nb2, wi, resd;
  int i, E=0;

  nb.d = x;
  /* Filter cases */
  if (nb.i[HI_ENDIAN] < 0x00100000){        /* x < 2^(-1022)    */
    if (((nb.i[HI_ENDIAN] & 0x7fffffff)|nb.i[LO_ENDIAN])==0)
      return 1.0/0.0;                       /* log(+/-0) = -Inf */
    if (nb.i[HI_ENDIAN] < 0) 
      return (x-x)/0;                       /* log(-x) = Nan    */

    /* Subnormal number */
    E    -= (SCS_NB_BITS*2); /* keep in mind that x is a subnormal number */ 
    nb.d *=SCS_RADIX_TWO_DOUBLE;  /* make x as normal number     */         
    /* We may just want add 2 to the scs number.index */
    /* may be .... we will see */
  }
  if (nb.i[HI_ENDIAN] >= 0x7ff00000)
    return x+x;                             /* Inf or Nan       */

  /* find n, nb.d such that sqrt(2)/2 < nb.d < sqrt(2) */
  E += (nb.i[HI_ENDIAN]>>20)-1023;
  nb.i[HI_ENDIAN] =  (nb.i[HI_ENDIAN] & 0x000fffff) | 0x3ff00000;
  if (nb.d > SQRT_2){
    nb.d *= 0.5;
    E++;
  }


  /* to normalize nb.d and round to nearest      */
  /* +((2^4 - trunc(sqrt(2)/2) *2^4 )*2 + 1)/2^5 */ 
  nb2.d = nb.d + norm_number.d; 
  i = (nb2.i[HI_ENDIAN] & 0x000fffff);
  i = i >> 16; /* 0<= i <=11 */
  
  wi.d = (11+i)*(double)0.6250e-1;

  /* (1+f-w_i) */
  nb.d -= wi.d; 
  
  /* Table reduction */
  ti     = table_ti_ptr[i]; 
  inv_wi = table_inv_wi_ptr[i];

 
  /* R = (1+f-w_i)/w_i */
  scs_set_d(R, nb.d);
  scs_mul(R, R, inv_wi);
 

  /* sc_ln2_r10_times_E = E*log10(2)  */
  scs_set(sc_ln2_r10_times_E, sc_ln2_r10_ptr);
  
  if (E >= 0){
    scs_mul_ui(sc_ln2_r10_times_E, (unsigned int) E);
  }else{
    scs_mul_ui(sc_ln2_r10_times_E, (unsigned int) -E);
    sc_ln2_r10_times_E->sign = -1;
  }


  /*
   * Polynomial evaluation of log10(1 + R) with an error less than 2^(-130)
   */
  scs_mul(res1, constant_poly_ptr[0], R);
  for(i=1; i<20; i++){
    scs_add(res1, constant_poly_ptr[i], res1);
    scs_mul(res1, res1, R);
  }
  scs_add(res1, res1, ti);
  scs_add(res1, res1, sc_ln2_r10_times_E);  


  scs_get_d_pinf(&resd.d, res1);  
  
  return resd.d;
}



