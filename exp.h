


/***********************************************************/
/*                      Second step                        */
/***********************************************************/

#ifdef WORDS_BIGENDIAN
static const db_number 
  iln2_o512 = {{0x40871547, 0x652B82FE}}, /* 512/ln(2)        */
  tiny_int  = {{0x3c8fffff, 0xffffffff}}, /* 2^(-54)-2^(-107) */
  smll_int  = {{0xc0874910, 0xd52d3051}}, /* -745,13...       */
  larg_int  = {{0x40862e42, 0xfefa39ef}}; /* 709,78...        */
#else
static const db_number 
  iln2_o512 = {{0x652B82FE, 0x40871547}}, /* 512/ln(2)        */
  tiny_int  = {{0xffffffff, 0x3c8fffff}}, /* 2^(-54)-2^(-107) */
  smll_int  = {{0xd52d3051, 0xc0874910}}, /* -745,13...       */
  larg_int  = {{0xfefa39ef, 0x40862e42}}; /* 709,78...        */ 
#endif

static const struct scs 
sc_ln2_o512_1 = {{0x00162e42, 0x3fbe8e7b, 0x335793c7, 0x19cc01f9,
		  0x1ed5e81e, 0x1a193394, 0x316c5b14, 0x00000000},
		 DB_ONE,  -1,   1 },
sc_ln2_o512_2 = {{0x068badc5, 0x355f457c, 0x3dc3b103, 0x1bd75930, 
		  0x2acaa97d, 0x295f4362, 0x07697571, 0x2b827042},
		 DB_ONE,  -8,   1 };

#define sc_ln2_o512_1_ptr ((scs_ptr)(&sc_ln2_o512_1))
#define sc_ln2_o512_2_ptr ((scs_ptr)(&sc_ln2_o512_2))


/* We should add 1+x to this polynom and multiply it by x^2 */
static const scs constant_poly[10]=  		/* 2^-164 */
/* ~2.505211e-08 */ 
{{{0x0000001a, 0x39915ae6, 0x271b3b26, 0x3bc3c8d6, 
0x327f2894, 0x1a80bba6, 0x2d7be149, 0x149ab262},
DB_ONE,  -1,   1 } 
,
/* ~2.755732e-07 */ 
{{0x00000127, 0x393ee8a7, 0x2647a705, 0x037af958, 
0x0c753fb3, 0x0941e129, 0x13162578, 0x1773558f},
DB_ONE,  -1,   1 } 
,
/* ~2.755732e-06 */ 
{{0x00000b8e, 0x3c74aad8, 0x399abcfd, 0x1225e04b, 
0x2e8dea13, 0x038b8517, 0x3212cece, 0x2cb577c5},
DB_ONE,  -1,   1 } 
,
/* ~2.480159e-05 */ 
{{0x00006806, 0x201a01a0, 0x066e4bb6, 0x156666d6, 
0x333d1ab0, 0x06772281, 0x044e4578, 0x04d8d60d},
DB_ONE,  -1,   1 } 
,
/* ~1.984127e-04 */ 
{{0x00034034, 0x00d00d00, 0x34034034, 0x0249047b, 
0x38e5a9d7, 0x30805847, 0x335f9062, 0x206ba2e1},
DB_ONE,  -1,   1 } 
,
/* ~1.388889e-03 */ 
{{0x0016c16c, 0x05b05b05, 0x2c16c16c, 0x121e8acb, 
0x25da960e, 0x2999616d, 0x2b3ae80c, 0x1a17471d},
DB_ONE,  -1,   1 } 
,
/* ~8.333333e-03 */ 
{{0x00888888, 0x22222222, 0x08888888, 0x22222220, 
0x072a661e, 0x15f0a5cf, 0x3fc53e81, 0x3da68844},
DB_ONE,  -1,   1 } 
,
/* ~4.166667e-02 */ 
{{0x02aaaaaa, 0x2aaaaaaa, 0x2aaaaaaa, 0x2aaaaa9e, 
0x2e389977, 0x07fcaf81, 0x07fee303, 0x13eab39c},
DB_ONE,  -1,   1 } 
,
/* ~1.666667e-01 */ 
{{0x0aaaaaaa, 0x2aaaaaaa, 0x2aaaaaaa, 0x2aaaaaaa, 
0x2aaaaadc, 0x1902174f, 0x32f46829, 0x13f2855b},
DB_ONE,  -1,   1 } 
,
/* ~5.000000e-01 */ 
{{0x20000000, 0x00000000, 0x00000000, 0x00000000, 
0x00000075, 0x21e5e9c2, 0x1e95e039, 0x17dad7b5},
DB_ONE,  -1,   1 } 
};


#define constant_poly_ptr (scs_ptr)&constant_poly



