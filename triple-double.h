/*
 *  triple_double.h
 *  
 * This file contains useful tools and data for triple double data representation.
 *
 */

#ifndef TRIPLE_DOUBLE_H
#define TRIPLE_DOUBLE_H 1

#include "scs_lib/scs.h"
#include "scs_lib/scs_private.h"
 /* undef all the variables that might have been defined in
    scs_lib/scs_private.h */
#undef VERSION 
#undef PACKAGE 
#undef HAVE_GMP_H
#undef HAVE_MPFR_H
#undef HAVE_MATHLIB_H
/* then include the proper definitions  */
#include "crlibm_config.h"

#ifdef HAVE_INTTYPES_H
#include <inttypes.h>
#endif


/* Renormalize3

   Procedure for renormalizing a triple double number, i.e.
   computing exactly an equivalent sum of three non-overlapping
   double numbers


   Arguments:       a triple double number ah, am, al
   
   Results:         a triple double number resh, resm, resl

   Preconditions:   abs(ah) > abs(am) > abs(al)
                    ah and am are overlapping not more than 51 bits
                    am and al are overlapping not more than 51 bits

   Guarantees:      abs(resh) > abs(resm) > abs(resl)
                    resh and resm are non-overlapping
		    resm and resl are non-overlapping
		    resm = round-to-nearest(resm + resl)

   Details:         resh, resm and resl are considered to be pointers

*/
#define Renormalize3(resh, resm, resl, ah, am, al)     \
{                                                      \
    double _t1h, _t1l, _t2l;                           \
                                                       \
    Add12(_t1h, _t1l, (am), (al));                     \
    Add12((*(resh)), _t2l, (ah), (_t1h));              \
    Add12((*(resm)), (*(resl)), _t2l, _t1l);           \
}


/* Mul23

   Procedure for multiplying two double double numbers resulting
   in a triple double number


   Arguments:       two double double numbers:
                    ah, al and
		    bh, bl
   
   Results:         a triple double number resh, resm, resl

   Preconditions:   abs(ah) > abs(al) 
                    ah and al do not overlap
		    ah = round-to-nearest(ah + al)
		    abs(bh) > abs(bl) 
                    bh and bl do not overlap
		    bh = round-to-nearest(bh + bl)
		    
   Guarantees:      resm and resl are non-overlapping
                    resm = round-to-nearest(resm + resl)
		    abs(resm) <= 2^(-49) * abs(resh)
		    resh+resm+resl = (ah+al) * (bh+bl) * (1 + eps)
		    where
		    abs(eps) <= 2^(-149)

   Details:         resh, resm and resl are considered to be pointers
*/
#define Mul23(resh, resm, resl, ah, al, bh, bl)                \
{                                                              \
    double _t1, _t2, _t3, _t4, _t5, _t6, _t7, _t8, _t9, _t10;  \
                                                               \
    Mul12((resh),&_t1,(ah),(bh));                              \
    Mul12(&_t2,&_t3,(ah),(bl));                                \
    Mul12(&_t4,&_t5,(al),(bh));                                \
    _t6 = (al) * (bl);                                         \
    Add22Cond(&_t7,&_t8,_t2,_t3,_t4,_t5);                      \
    Add12(_t9,_t10,_t1,_t6);                                   \
    Add22Cond((resm),(resl),_t7,_t8,_t9,_t10);                 \
}

/* Mul233

   Procedure for multiplying a double double number by 
   a triple double number resulting in a triple double number


   Arguments:       a double double number ah, al
                    a triple double number bh, bm, bl
   
   Results:         a triple double number resh, resm, resl

   Preconditions:   abs(ah) > abs(al)
                    ah and al do not overlap
		    ah = round-to-nearest(ah + al)
		    abs(bm) <= 2^(-b_o) * abs(bh)
		    abs(bl) <= 2^(-b_u) * abs(bm)
		    where
		    b_o >= 2
		    b_u >= 1
		    
   Guarantees:      resm and resl are non-overlapping
                    resm = round-to-nearest(resm + resl)
		    abs(resm) <= 2^(\gamma) * abs(resh)
		    where
		    \gamma >= min(48,b_o-4,b_o+b_u-4)
		    resh+resm+resl=(ah+al) * (bh+bm+bl) * (1+eps)
		    where
		    abs(eps) <= 
                       (2^(-99-b_o) + 2^(-99-b_o-b_u) + 2^(-152)) / 
		         (1 - 2^(-53) - 2^(-b_o+1) - 2^(-b_o-b_u+1))

   Details:         resh, resm and resl are considered to be pointers
*/
#define Mul233(resh, resm, resl, ah, al, bh, bm, bl)            \
{                                                               \
    double _t1, _t2, _t3, _t4, _t5, _t6, _t7, _t8, _t9, _t10;   \
    double _t11, _t12, _t13, _t14, _t15, _t16, _t17, _t18;      \
                                                                \
    Mul12((resh),&_t1,(ah),(bh));                               \
    Mul12(&_t2,&_t3,(ah),(bm));                                 \
    Mul12(&_t4,&_t5,(ah),(bl));                                 \
    Mul12(&_t6,&_t7,(al),(bh));                                 \
    Mul12(&_t8,&_t9,(al),(bm));                                 \
    _t10 = (al) * (bl);                                         \
    Add22Cond(&_t11,&_t12,_t2,_t3,_t4,_t5);                     \
    Add22Cond(&_t13,&_t14,_t6,_t7,_t8,_t9);                     \
    Add22Cond(&_t15,&_t16,_t11,_t12,_t13,_t14);                 \
    Add12Cond(_t17,_t18,_t1,_t10);                              \
    Add22Cond((resm),(resl),_t17,_t18,_t15,_t16);               \
}




/* Add33

   Procedure for adding two triple double numbers resulting
   in a triple double number


   Arguments:       two triple double numbers:
                    ah, am, al and
		    bh, bm, bl
   
   Results:         a triple double number resh, resm, resl

   Preconditions:   abs(bh) <= 0.75 * abs(ah) (i)
                    abs(am) <= 2^(-a_o) * abs(ah)
		    abs(al) <= 2^(-a_u) * abs(am)
		    abs(bm) <= 2^(-b_o) * abs(bh)
		    abs(bl) <= 2^(-b_u) * abs(bm)
		    where
		    b_o >= a_o >= 4
		    b_u >= a_u >= 4

		    Condition (i) may not be respected if 
		    one can assume in this case that ah=am=al
		    
   Guarantees:      resm and resl are non-overlapping
                    resm = round-to-nearest(resm + resl)
		    abs(resm) <= 2^(-a_o + 5) * abs(resh)
		    resh+resm+resl = (ah+am+al + bh+bm+bl) * (1+eps)
                    where 
		    abs(eps) <= 2^(-a_o-a_u-47) + 2^(-a_o-98)

   Details:         resh, resm and resl are considered to be pointers
*/
#define Add33(resh, resm, resl, ah, am, al, bh, bm, bl)      \
{                                                            \
    double _t1, _t2, _t3, _t4, _t5, _t6, _t7, _t8;           \
                                                             \
    Add12((*(resh)),_t1,(ah),(bh));                          \
    Add12Cond(_t2,_t3,(am),(bm));                            \
    _t6 = (al) + (bl);                                       \
    Add12Cond(_t7,_t4,_t1,_t2);                              \
    _t5 = _t3 + _t4;                                         \
    _t8 = _t5 + _t6;                                         \
    Add12Cond((*(resm)),(*(resl)),_t7,_t8);                  \
}



/* Add233

   Procedure for adding a double double number to a triple 
   double number resulting in a triple double number


   Arguments:       a double double number ah, al
                    a triple double number bh, bm, bl
   
   Results:         a triple double number resh, resm, resl

   Preconditions:   abs(ah) > abs(al)
                    ah and al do not overlap
		    ah = round-to-nearest(ah + al)
		    abs(bh) <= 2^(-2) * abs(ah)
		    abs(bm) <= 2^(-b_o) * abs(bh)
		    abs(bl) <= 2^(-b_u) * abs(bm)
		    where
		    b_o >= 2
		    b_u >= 1
		    
   Guarantees:      resm and resl are non-overlapping
                    resm = round-to-nearest(resm + resl)
		    abs(resm) <= 2^(\gamma) * abs(resh)
		    where
		    \gamma >= min(45,b_o-4,b_o+b_u-2)
		    resh+resm+resl=((ah+al) + (bh+bm+bl)) * (1+eps)
		    where
		    abs(eps) <= 
                       <= 2^(-b_o-b_u-52) + 2^(-b_o-104) + 2^(-153)

   Details:         resh, resm and resl are considered to be pointers
*/
#define Add233(resh, resm, resl, ah, al, bh, bm, bl)            \
{                                                               \
    double _t1, _t2, _t3, _t4, _t5, _t6, _t7;                   \
                                                                \
    Add12((*(resh)),_t1,(ah),(bh));                             \
    Add12Cond(_t2,_t3,(al),(bm));                               \
    Add12Cond(_t4,_t5,_t1,_t2);                                 \
    _t6 = _t3 + (bl);                                           \
    _t7 = _t6 + _t5;                                            \
    Add12Cond((*(resm)),(*(resl)),_t4,_t7);                     \
}

/* ReturnRoundToNearest3

   Procedure for rounding a triple to a double number
   in round-to-nearest-ties-to-even mode.


   Arguments:       a triple double number xh, xm, xl
   
   Results:         a double number xprime 
                    returned by a return-statement

   Preconditions:   xh, xm and xl are non-overlapping
                    xm = RN(xm +math xl)
		    xh != 0, xm != 0
		    xl = 0 iff xm != +/- 0.5 * ulp(xh) (0.25 if xh = 2^e)
		    		    
   Guarantees:      xprime = RN(xh + xm + xl)

   Sideeffects:     returns, i.e. leaves the function

*/
#define ReturnRoundToNearest3(xh,xm,xl)                       \
{                                                             \
    double _t1, _t2, _t3, _t4, _t5, _t6;                      \
    db_number _xp, _xn;                                       \
                                                              \
    _xp.d = (xh);                                             \
    _xn.i[HI] = _xp.i[HI];                                    \
    _xn.i[LO] = _xp.i[LO];                                    \
    _xn.l--;                                                  \
    _t1 = _xn.d;                                              \
    _xp.l++;                                                  \
    _t4 = _xp.d;                                              \
    _t2 = (xh) - _t1;                                         \
    _t3 = _t2 * -0.5;                                         \
    _t5 = _t4 - (xh);                                         \
    _t6 = _t5 * 0.5;                                          \
    if (((xm) != _t3) && ((xm) != _t6)) return ((xh) + (xm)); \
    if ((xm) * (xl) > 0.0) {                                  \
      if ((xh) * (xl) > 0.0)                                  \
        return _t4;                                           \
      else                                                    \
        return _t1;                                           \
    } else return (xh);                                       \
}

/* ReturnRoundUpwards3

   Procedure for rounding a triple to a double number
   in round-upwards mode.


   Arguments:       a triple double number xh, xm, xl
   
   Results:         a double number xprime 
                    returned by a return-statement

   Preconditions:   xh, xm and xl are non-overlapping
                    xm = RN(xm +math xl)
		    xh != 0, xm != 0
		    		    
   Guarantees:      xprime = RU(xh + xm + xl)

   Sideeffects:     returns, i.e. leaves the function

*/
#define ReturnRoundUpwards3(xh,xm,xl)                         \
{                                                             \
    double _t1, _t2, _t3;                                     \
    db_number _tdb;                                           \
                                                              \
    Add12(_t1,_t2,(xh),(xm));                                 \
    _t3 = _t2 + (xl);                                         \
    if (_t3 > 0.0) {                                          \
      if (_t1 > 0.0) {                                        \
         _tdb.d = _t1;                                        \
         _tdb.l++;                                            \
         return _tdb.d;                                       \
      } else {                                                \
         _tdb.d = _t1;                                        \
         _tdb.l--;                                            \
         return _tdb.d;                                       \
      }                                                       \
    } else return _t1;                                        \
}


/* ReturnRoundDownwards3

   Procedure for rounding a triple to a double number
   in round-downwards mode.


   Arguments:       a triple double number xh, xm, xl
   
   Results:         a double number xprime 
                    returned by a return-statement

   Preconditions:   xh, xm and xl are non-overlapping
                    xm = RN(xm +math xl)
		    xh != 0, xm != 0
		    		    
   Guarantees:      xprime = RD(xh + xm + xl)

   Sideeffects:     returns, i.e. leaves the function

*/
#define ReturnRoundDownwards3(xh,xm,xl)                       \
{                                                             \
    double _t1, _t2, _t3;                                     \
    db_number _tdb;                                           \
                                                              \
    Add12(_t1,_t2,(xh),(xm));                                 \
    _t3 = _t2 + (xl);                                         \
    if (_t3 < 0.0) {                                          \
      if (_t1 > 0.0) {                                        \
         _tdb.d = _t1;                                        \
         _tdb.l--;                                            \
         return _tdb.d;                                       \
      } else {                                                \
         _tdb.d = _t1;                                        \
         _tdb.l++;                                            \
         return _tdb.d;                                       \
      }                                                       \
    } else return _t1;                                        \
}


/* ReturnRoundTowardsZero3

   Procedure for rounding a triple to a double number
   in round-towards-zero mode.


   Arguments:       a triple double number xh, xm, xl
   
   Results:         a double number xprime 
                    returned by a return-statement

   Preconditions:   xh, xm and xl are non-overlapping
                    xm = RN(xm +math xl)
		    xh != 0, xm != 0
		    		    
   Guarantees:      xprime = RZ(xh + xm + xl)

   Sideeffects:     returns, i.e. leaves the function

*/
#define ReturnRoundTowardsZero3(xh,xm,xl)                     \
{                                                             \
    double _t1, _t2, _t3;                                     \
    db_number _tdb;                                           \
                                                              \
    Add12(_t1,_t2,(xh),(xm));                                 \
    _t3 = _t2 + (xl);                                         \
    if (_t1 > 0.0) {                                          \
       if (_t3 < 0.0) {                                       \
         _tdb.d = _t1;                                        \
         _tdb.l--;                                            \
         return _tdb.d;                                       \
       } else return _t1;                                     \
    } else {                                                  \
       if (_t3 > 0.0) {                                       \
         _tdb.d = _t1;                                        \
         _tdb.l--;                                            \
         return _tdb.d;                                       \
       } else return _t1;                                     \
    }                                                         \
}


#endif /*TRIPLE_DOUBLE_H*/
