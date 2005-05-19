


/* Mul23
   Arguments: 2 double double numbers
   Returns: 1 triple double number
   Relative error: less than 2^(-150) // A VERIFIER

   We do not use the conditional versions of the Mul12 routines.
   ah and bh must therefore be bounded by abs(ah), abs(bh) < 2^(970)

   The arguments are supposed to be non-overlapping but may contain gaps.
   The result will be non-overlapping but may contain gaps.

   The use of the conditional Add?2 seems to be necessary as we cannot simply 
   bound and compare the intermediate lower parts.
*/

/* On n utilise pas la methode de Alexey Ershov qui consiste a separer les arguments en deux doubles
   et de les multiplier "a la main".
   La methode semble devenir suboptimale quand on a un fma; c est bizarre qu Alexey n en tienne pas compte

   On se pose donc sur les Mul12 et Mul22 qui a leur tour vont utiliser le fma
*/
/*

This is the version using a function. 

void Mul23(double *resh, double *resm, double *resl, double ah, double al, double bh, double bl) {
  double ahbhh, ahbhl, ahblh, ahbll, albhh, albhl, albl, mth, mtl, mtth, mttl, mttth, mtttl, mt4h, mt4l, mt5;

  Mul12(&ahbhh, &ahbhl, ah, bh);   
  Mul12(&ahblh, &ahbll, ah, bl);   
  Mul12(&albhh, &albhl, al, bh);   
  albl = al * bl;                  
  Add22Cond(&mth, &mtl, ahblh, ahbll, albhh, albhl);  
  Add12Cond(mtth, mttl, ahbhl, albl);  // A VERIFIER: on devrait pouvoir prouver une condition ici  
  Add22Cond(&mttth, &mtttl, mth, mtl, mtth, mttl);  

  Add12Cond(mt4h, mt4l, ahbhh, mttth); // A VERIFIER: la renormalisation est certainement sous optimale
  Add12Cond(mt5, *resl, mt4l, mtttl);
  Add12Cond(*resh, *resm, mt4h, mt5); 
}  
*/


#define Mul23(resh, resm, resl, ah, al, bh, bl)                                                              \
{                                                                                                            \
  double _ahbhh, _ahbhl, _ahblh, _ahbll, _albhh, _albhl;                                                     \
  double _albl, _mth, _mtl, _mtth, _mttl, _mttth, _mtttl;                                                    \
  double _mt4h, _mt4l, _mt5;                                                                                 \
  Mul12(&_ahbhh, &_ahbhl, (ah), (bh));                                                                       \
  Mul12(&_ahblh, &_ahbll, (ah), (bl));                                                                       \
  Mul12(&_albhh, &_albhl, (al), (bh));                                                                       \
  _albl = (al) * (bl);                                                                                       \
  Add22Cond(&_mth, &_mtl, _ahblh, _ahbll, _albhh, _albhl);                                                   \
  Add12Cond(_mtth, _mttl, _ahbhl, _albl);                                                                    \
  Add22Cond(&_mttth, &_mtttl, _mth, _mtl, _mtth, _mttl);                                                     \
  Add12Cond(_mt4h, _mt4l, _ahbhh, _mttth);                                                                   \
  Add12Cond(_mt5, (*(resl)), _mt4l, _mtttl);                                                                 \
  Add12Cond((*(resh)), (*(resm)), _mt4h, _mt5);                                                              \
} 

/* Mul233
   Arguments: 1 double double number, 1 triple double number
   Returns: 1 triple double number
   Relative error: less than 2^(-150) // A VERIFIER

   We do not use the conditional versions of the Mul12 routines.
   ah and bh must therefore be bounded by abs(ah), abs(bh) < 2^(970)
   
   The arguments are supposed to be non-overlapping but may contain gaps.
   The result will be non-overlapping but may also contain gaps.

   The use of the conditional Add?2 seems to be necessary as we cannot simply 
   bound and compare the intermediate lower parts.
*/
/*

This is the version using a function.

void Mul233(double *resh, double *resm, double *resl, double ah, double al, double bh, double bm, double bl) {
  double ahbhh, ahbhl, ahbmh, ahbml, ahblh, ahbll, albhh, albhl, albmh, albml, albl;
  double t1h, t1l, t2h, t2l, t3h, t3l, t4h, t4l;

  Mul12(&ahbhh, &ahbhl, ah, bh); 
  Mul12(&ahbmh, &ahbml, ah, bm); 
  Mul12(&ahblh, &ahbll, ah, bl); 
  Mul12(&albhh, &albhl, al, bh);  
  Mul12(&albmh, &albml, al, bm);  
  albl = al * bl;  // A VERIFIER: on peut peut etre s en passer  
  Add22Cond(&t1h, &t1l, albmh, albml, albhh, albhl);  
  Add22Cond(&t2h, &t2l, ahbmh, ahbml, ahblh, ahbll);  
  Add22Cond(&t3h, &t3l, t1h, t1l, t2h, t2l);  


  Add22Cond(&t4h, resl, ahbhl, albl, t3h, t3l);  // A VERIFIER: la renormalisation est-elle correcte ?
  Add12Cond(*resh, *resm, ahbhh, t4h);  
}  
*/

#define Mul233(resh, resm, resl, ah, al, bh, bm, bl)                                               \
{                                                                                                  \
  double _ahbhh, _ahbhl, _ahbmh, _ahbml, _ahblh, _ahbll, _albhh, _albhl, _albmh, _albml, _albl;    \
  double _t1h, _t1l, _t2h, _t2l, _t3h, _t3l, _t4h, _t4l;                                           \
  Mul12(&_ahbhh, &_ahbhl, (ah), (bh));                                                             \
  Mul12(&_ahbmh, &_ahbml, (ah), (bm));                                                             \
  Mul12(&_ahblh, &_ahbll, (ah), (bl));                                                             \
  Mul12(&_albhh, &_albhl, (al), (bh));                                                             \
  Mul12(&_albmh, &_albml, (al), (bm));                                                             \
  _albl = (al) * (bl);                                                                             \
  Add22Cond(&_t1h, &_t1l, _albmh, _albml, _albhh, _albhl);                                         \
  Add22Cond(&_t2h, &_t2l, _ahbmh, _ahbml, _ahblh, _ahbll);                                         \
  Add22Cond(&_t3h, &_t3l, _t1h, _t1l, _t2h, _t2l);                                                 \
  Add22Cond(&_t4h, (resl), _ahbhl, _albl, _t3h, _t3l);                                             \
  Add12Cond((*(resh)), (*(resm)), _ahbhh, _t4h);                                                   \
} 



/* Add33
   Arguments: 2 triple double numbers
   Returns: 1 triple double number
   Relative error: less than 2^(-150) // A VERIFIER

   The arguments are supposed to be non overlapping but may contain gaps.
   The result will be non-overlapping but may contain gaps.

   We currently use the conditional versions of Add?2 since we have no methodology of
   bounding the lower part values.
*/

/*

This is the version using a function.

void Add33(double *resh, double *resm, double *resl, double ah, double am, double al, double bh, double bm, double bl) {
  double ahPbhh, ahPbhl, amPbmh, amPbml, alPbl, t1h, t1l, t2, t3h, t3l, t4;
  

  // A VERIFIER: il suffirait en fait de faire un seul test de grandeur par composante, donc les versions
  //   conditionnelles sont peut etre de l overkill

  // Alexey Ershov fait pareil (inclus les tests de conditions dans les fonctions pas specialisees) mais il ne renormalise pas 

  Add12Cond(ahPbhh, ahPbhl, ah, bh);
  Add12Cond(amPbmh, amPbml, am, bm);
  Add12Cond(t1h, t1l, ahPbhl, amPbmh);
  t2 = (amPbml + t1l) + (al + bl);

  // A VERIFIER: peut etre pas de conditions dans ce qui suit 
  Add12Cond(t3h, t3l, ahPbhh, t1h);
  Add12Cond(t4, *resl, t3l, t2);
  Add12Cond(*resh, *resm, t3h, t4);
}  
*/

#define Add33(resh, resm, resl, ah, am, al, bh, bm, bl)                                  \
{                                                                                        \
  double _ahPbhh, _ahPbhl, _amPbmh, _amPbml, _alPbl, _t1h, _t1l, _t2, _t3h, _t3l, _t4;   \
  Add12Cond(_ahPbhh, _ahPbhl, (ah), (bh));                                               \
  Add12Cond(_amPbmh, _amPbml, (am), (bm));                                               \
  Add12Cond(_t1h, _t1l, _ahPbhl, _amPbmh);                                               \
  _t2 = (_amPbml + _t1l) + ((al) + (bl));                                                \
  Add12Cond(_t3h, _t3l, _ahPbhh, _t1h);                                                  \
  Add12Cond(_t4, (*(resl)), _t3l, _t2);                                                  \
  Add12Cond((*(resh)), (*(resm)), _t3h, _t4);                                            \
} 

/* Add233
   Arguments: 1 double double number, 1 triple double number
   Returns: 1 triple double number
   Relative error: less than 2^(-150) // A VERIFIER
 
   The arguments are supposed to be non-overlapping but may contain gaps.
   The result will be non-overlapping but may also contain gaps.

   The use of the conditional Add?2 seems to be necessary as we cannot simply 
   bound and compare the intermediate lower parts.
*/

/* 

 This is the version using a function

void Add233(double *resh, double *resm, double *resl, double ah, double al, double bh, double bm, double bl) {
  double ahPbhh, ahPbhl, alPbmh, alPbml, t1h, t1l, t2h, t2l, t3, t4, t5, t6h, t6l, t7;

  Add12Cond(ahPbhh, ahPbhl, ah, bh);
  Add12Cond(alPbmh, alPbml, al, bm);
  Add12Cond(t1h, t1l, ahPbhl, alPbmh);
  t3 = alPbml + bl;
  t4 = t3 + t1l;

  // A VERIFIER: peut etre pas de condition dans ce qui suit
  Add12Cond(t6h, t6l, ahPbhh, t1h);
  Add12Cond(t7, *resl, t6l, t4);
  Add12Cond(*resh, *resm, t6h, t7);
}

*/

#define Add233(resh, resm, resl, ah, al, bh, bm, bl)                 \
{                                                                    \
  double _ahPbhh, _ahPbhl, _alPbmh, _alPbml, _t1h, _t1l;             \
  double _t2h, _t2l, _t3, _t4, _t5, _t6h, _t6l, _t7;                 \
  Add12Cond(_ahPbhh, _ahPbhl, (ah), (bh));                           \
  Add12Cond(_alPbmh, _alPbml, (al), (bm));                           \
  Add12Cond(_t1h, _t1l, _ahPbhl, _alPbmh);                           \
  _t3 = _alPbml + (bl);                                              \
  _t4 = _t3 + _t1l;                                                  \
  Add12Cond(_t6h, _t6l, _ahPbhh, _t1h);                              \
  Add12Cond(_t7, (*(resl)), _t6l, _t4);                              \
  Add12Cond((*(resh)), (*(resm)), _t6h, _t7);                        \
}
