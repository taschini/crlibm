
#####################################################################
# Useful procedures for IEEE doubles


#---------------------------------------------------------------------

log2:=proc(x) evalf( log[2](x)) end proc:


#---------------------------------------------------------------------
# ieeedouble converts a number to IEEE double format.
# returns sign (-1 or 1), exponent between -1022 and 1023, mantissa as a fraction between 0.5 and 1.
ieeedouble:=proc(xx)
local x, sgn, logabsx, exponent, mantissa, infmantissa,powermin,powermax,expmin,expmax,expmiddle,powermiddle;
Digits := 100;
x := evalf(xx);
if (x=0) then sgn, exponent, mantissa := 1, -1022, 0
else
  if (x < 0) then sgn := -1
  else sgn := 1
  fi:
  x := abs(x);
  if x >=  2^(1023)*(2-2^(-53)) then mantissa := infinity; exponent := 1023
  else if x <= 2^(-1075) then mantissa := 0; exponent := -1022
      else
         if x <= 2^(-1022) then exponent := -1022
         else
# x is between 2^(-1022) and 2^(1024)
         powermin := 2^(-1022); expmin := -1022;
         powermax := 2^1024; expmax := 1024;
         while (expmax-expmin > 1) do
            expmiddle := round((expmax+expmin)/2);
            powermiddle := 2^expmiddle;
            if x >= powermiddle then
                powermin := powermiddle;
                expmin := expmiddle
            else
                powermax := powermiddle;
                expmax := expmiddle
            fi
          od;
# now, expmax - expmin = 1 and powermin <= x < powermax,
# powermin = 2^expmin and powermax = 2^expmax, so expmin is the exponent of x
         exponent := expmin;
         fi;
         infmantissa := x*2^(52-exponent);
	 if frac(infmantissa) <> 0.5 then mantissa := round(infmantissa)
            else
              mantissa := floor(infmantissa);
               if type(mantissa,odd) then mantissa := mantissa+1 fi
            fi;
         mantissa := mantissa*2^(-52);
      fi;
  fi;
fi;
sgn,exponent,mantissa;
end:

#---------------------------------------------------------------------
# ieeedoubleRU converts a number to IEEE double format rounding upwards.
# returns sign (-1 or 1), exponent between -1022 and 1023, mantissa as a fraction between 0.5 and 1.
ieeedoubleRU:=proc(xx)
local x, sgn, logabsx, exponent, mantissa, infmantissa,powermin,powermax,expmin,expmax,expmiddle,powermiddle;
Digits := 100;
x := evalf(xx);
if (x=0) then sgn, exponent, mantissa := 1, -1022, 0
else
  if (x < 0) then sgn := -1
  else sgn := 1
  fi:
  x := abs(x);
  if x >=  2^(1023)*(2-2^(-53)) then mantissa := infinity; exponent := 1023
  else if x <= 2^(-1075) then mantissa := 0; exponent := -1022
      else
         if x <= 2^(-1022) then exponent := -1022
         else
# x is between 2^(-1022) and 2^(1024)
         powermin := 2^(-1022); expmin := -1022;
         powermax := 2^1024; expmax := 1024;
         while (expmax-expmin > 1) do
            expmiddle := round((expmax+expmin)/2);
            powermiddle := 2^expmiddle;
            if x >= powermiddle then
                powermin := powermiddle;
                expmin := expmiddle
            else
                powermax := powermiddle;
                expmax := expmiddle
            fi
          od;
# now, expmax - expmin = 1 and powermin <= x < powermax,
# powermin = 2^expmin and powermax = 2^expmax, so expmin is the exponent of x
         exponent := expmin;
         fi;
         infmantissa := x*2^(52-exponent);
	 if frac(infmantissa) <> 0 then
         if (sgn > 0) then
             mantissa := ceil(infmantissa);
         else
             mantissa := floor(infmantissa);
         fi;
            else
              mantissa := infmantissa;
            fi;
         mantissa := mantissa*2^(-52);
      fi;
  fi;
fi;
sgn,exponent,mantissa;
end:

#---------------------------------------------------------------------
# ieeedoubleRD converts a number to IEEE double format rounding downwards.
# returns sign (-1 or 1), exponent between -1022 and 1023, mantissa as a fraction between 0.5 and 1.
ieeedoubleRD:=proc(xx)
local x, sgn, logabsx, exponent, mantissa, infmantissa,powermin,powermax,expmin,expmax,expmiddle,powermiddle;
Digits := 100;
x := evalf(xx);
if (x=0) then sgn, exponent, mantissa := 1, -1022, 0
else
  if (x < 0) then sgn := -1
  else sgn := 1
  fi:
  x := abs(x);
  if x >=  2^(1023)*(2-2^(-53)) then mantissa := infinity; exponent := 1023
  else if x <= 2^(-1075) then mantissa := 0; exponent := -1022
      else
         if x <= 2^(-1022) then exponent := -1022
         else
# x is between 2^(-1022) and 2^(1024)
         powermin := 2^(-1022); expmin := -1022;
         powermax := 2^1024; expmax := 1024;
         while (expmax-expmin > 1) do
            expmiddle := round((expmax+expmin)/2);
            powermiddle := 2^expmiddle;
            if x >= powermiddle then
                powermin := powermiddle;
                expmin := expmiddle
            else
                powermax := powermiddle;
                expmax := expmiddle
            fi
          od;
# now, expmax - expmin = 1 and powermin <= x < powermax,
# powermin = 2^expmin and powermax = 2^expmax, so expmin is the exponent of x
         exponent := expmin;
         fi;
         infmantissa := x*2^(52-exponent);
	 if frac(infmantissa) <> 0 then
         if (sgn < 0) then
             mantissa := ceil(infmantissa);
         else
             mantissa := floor(infmantissa);
         fi;
            else
              mantissa := infmantissa;
            fi;
         mantissa := mantissa*2^(-52);
      fi;
  fi;
fi;
sgn,exponent,mantissa;
end:




#---------------------------------------------------------------------
# pulp returns the precision of the ulp of x

pulp:=proc(x)
local flt, ulpy:
flt:=ieeedouble(x):
ulpy:=-52+flt[2]:
end proc:



#---------------------------------------------------------------------
# ulp returns the absolute value of the ulp of x
ulp:=proc(x)
2**(pulp(x)):
end proc:




#---------------------------------------------------------------------
# Returns nearest IEEE double:
nearest := proc(x)
  local sgn, exponent, mantissa:

  sgn,exponent,mantissa := ieeedouble(x):
  sgn*mantissa*2^(exponent):

end:

#---------------------------------------------------------------------
# Returns RU IEEE double:
roundUp := proc(x)
  local sgn, exponent, mantissa:

  sgn,exponent,mantissa := ieeedoubleRU(x):
  sgn*mantissa*2^(exponent):

end:

#---------------------------------------------------------------------
# Returns RD IEEE double:
roundDown := proc(x)
  local sgn, exponent, mantissa:

  sgn,exponent,mantissa := ieeedoubleRD(x):
  sgn*mantissa*2^(exponent):

end:

#---------------------------------------------------------------------
# Returns RZ IEEE double:
roundToZero := proc(x)
    if evalf(x) > 0 then roundDown(x) else roundUp(x) fi:
end:




#---------------------------------------------------------------------
# ieehexa returns a string containing the hexadecimal representation of the double nearest to input x.

ieeehexa:= proc(x)
local  hex2, xx, longint, expo, sgn, frac, resultat:
    if(x=0) then resultat:=["00000000","00000000"]:
    elif(x=-0) then resultat:=["80000000","00000000"]:   # nice try
    else
        xx:=ieeedouble(x):
        sgn:=xx[1]:
        expo:=xx[2]:
        frac:=xx[3]:
        if (expo = -1023) then
            longint := (frac)*2^51 :   # subnormal
        else
            longint := (frac-1)*2^52 +   (expo+1023)*2^52:
        fi:
        if (sgn=-1) then
            longint := longint + 2^63:
        fi:
        longint := longint + 2^64:  # to get all the hexadecimal digits when we'll convert to string
        hex2:=convert(longint, hex):
        hex2:=convert(hex2, string):

        resultat:=[substring(hex2,2..9), substring(hex2,10..18)]:
    fi:
    resultat:
end proc:

ieeehexaString := proc(x)
	local hex, result:
	hex := ieeehexa(x):
	result := cat(hex[1],hex[2]):
	return result:
end proc:



#---------------------------------------------------------------------
# reciprocal of the previous
hexa2ieee:= proc(hexa)
local dec, bin, expo, mantis, sgn, hex1, hex2, hexcat, res:

    hex1:= op(1, hexa):
    hex2:= op(2, hexa):
    hexcat:= cat(hex1, hex2):
    dec:= convert(hexcat, decimal, hex):

    if(dec >= 2^63) then
        dec := dec - 2^63:
        sgn:= -1:
    else
        sgn:= 1:
    fi:
    expo:= trunc(dec/(2^52)) - 1023:
    if(expo=-1023) then
        mantis:= frac(dec/(2^51)): # denormal
    else
        mantis:= 1+frac(dec/(2^52)):
    fi:
    res:= evalf(sgn*2^(expo)*mantis):
    res:
end proc:

#---------------------------------------------------------------------



# Print a number x in Low or Big Endian representation in opened file "fd":
printendian:=proc(fd,x,isbig)
local xhl:
xhl:=ieeehexa(x):

if(isbig=0 or isbig=1) then
  fprintf(fd,"{{0x%+0.8s,0x%+0.8s}} /* %+0.10e */", xhl[2-isbig], xhl[isbig+1], x):
else
  print("ERROR, isbig must be equal to 0 or 1"):
end if:
end proc:



#---------------------------------------------------------------------
# hi_lo takes an arbitrary precision number x and returns two doubles such that:
# x ~ x_hi + x_lo
hi_lo:= proc(x)
local x_hi, x_lo, res:
x_hi:= nearest(evalf(x)):
res:=x-x_hi:
if (res = 0) then
  x_lo:=0:
else
  x_lo:=nearest(evalf(res)):
end if:
x_hi,x_lo:
end:



#---------------------------------------------------------------------
# same as hi_lo, but returns hexadecimal strings
ieeehexa2:=proc(x)
local reshi, reslo, hexhi, hexlo:
reshi:=nearest(x):
hexhi:=ieee2hexa(reshi):
reslo:= nearest(x-reshi):
hexlo:=ieee2hexa(reslo):
reshi, reslo:
end proc:




#---------------------------------------------------------------------
# Computes the constant for the round-to-nearest test.
# epsilon is the overall relative error of the approximation scheme
# See the documentation for the explanation of these formulae
compute_rn_constant := proc(epsilon)
  local k, constGenCase, constPowerOf2:
  k := trunc(-log[2](epsilon)) - 53:
  (1 +     2**54*epsilon / (1 - epsilon - 2**(-k+1) ) )  / (1-2**(-53))
end proc:



#---------------------------------------------------------------------
# Takes a real number, and prints the bits after the 53th of its nearest IEEE floating-point number

showHowDifficultToRound:=proc(x)
local xb,xs,s,e,m:
    Digits:=200:
    s,e,m := ieeedouble(x):
    xb:=convert(evalf(x*2^(-e)),binary):
    xs:=convert(xb, string):
    substring(xs,55..153)
end proc:

#---------------------------------------------------------------------
# Computes the floating point successor of x
succDouble:=proc(x)
local s,he,hexcat,hehi,helo,castx,shex,neg;
he := ieeehexa(x);
hehi:= op(1, he);
helo:= op(2, he);
hexcat := cat(hehi, helo);
neg := 1;
castx := convert(hexcat, decimal, hex);
if (castx >= 2^(63)) then
	castx := castx - 2^(63);
	neg := -1;
end if;
castx := castx + neg;
shex := convert(convert(castx, hex),string);
s := neg * hexa2ieee([substring(shex,1..8), substring(shex,9..16)]);
s;
end proc:







#####################################################################

# Stuff about truncated polynomials

# Truncate a polynomial

#---------------------------------------------------------------------
#poly_exact takes a polynomial in x with arbitrary precision
# coefficients, and returns a truncated polynomial where coefficients
# are IEEE doubles.

poly_exact:=proc(P)
local deg,i, coef, coef_t, Q:
Q:= 0:
convert(Q, polynom):
deg:=degree(P,x):
  for i from 0 to deg do
    coef:=coeff(P,x,i):
    coef_t:=nearest(coef):
    Q:= Q + coef_t*x^i:
  od:
return(Q):
end:


#---------------------------------------------------------------------
#Like poly_exact, but the n first coefficients are exactly
# representable as the sum of two doubles.  (to actually get the two
# doubles, use procedure hi_lo)

poly_exact2:=proc(P,n)
local deg,i, coef, coef_hi, coef_lo, Q:
Q:= 0:
convert(Q, polynom):
deg:=degree(P,x):
  for i from 0 to deg do
    coef :=coeff(P,x,i):
    coef_hi, coef_lo:=hi_lo(coef):
    Q:= Q + coef_hi*x^i:
    if(i<n) then
        Q := Q + coef_lo*x^i:
    fi:
  od:
  return(Q):
end:


#---------------------------------------------------------------------
#  OBSOLETE use compute_horner_rounding_error below
# Compute a bound on the accumulated rounding error caused by the Horner evaluation of a truncated polynomial
# P is the polynomial.
# xmax is the max value of |x|.
# n is the degree when P is computed in double double. The first double-double operation is an addition.

# returns max absolute error, min of the function, max of the function.

# This procedure also checks on the fly that the fast (test-free) versions of the double-double addition can be used, i.e. that for all x, at each Horner step i computing ci+x*Si, we have |ci|>|x*Si|. It prints warnings if it not the case.

compute_abs_rounding_error:=proc(poly,xmax, nn)
local n, deg, delta, deltap, i, S, P, Snorm, Smin, Smax, prec:
deltap:=0:
delta:=0:
deg:=degree(poly):

prec:=53: # precision of the first iterations

S:=coeff(poly, x, deg):
Smax:=abs(S):
Smin:=Smax:

if nn<0 then n:=0: else n:=nn: fi:# sometimes called by compute_rel_rounding_error with n=-1

for i from (deg-1) to 0 by -1 do
  P:= convert(S*x, polynom):
  Smin := abs(coeff(poly,x,i)) - xmax*Smax :
  if(Smin<=0) then
    printf("Warning! in compute_abs_rounding_error, Smin<=0 at iteration %d, consider decreasing xmax\n",i):
  fi:
  delta:= evalf(xmax*deltap + 2**(-prec)*xmax*Smax):
  if i<n then
    # fast Add22 ?
    if abs(coeff(poly,x,i)) < xmax*Smax  # may be improved to xmax*Smax/2
    then printf("WARNING Add22 cannot be used at step %d, use Add22Cond\n" , i  ):
         printf("    coeff=%1.20e,  xmax*Smax=%1.20e"  ,  abs(coeff(poly,x,i)),  xmax*Smax ):
    fi:
  fi:
  S:=convert(P+coeff(poly,x,i), polynom):
  Snorm:=evalf(numapprox[infnorm](S, x=-xmax..xmax)):
  if i=n-1 then prec:=100: fi:  # from the addition of the n-1-th iteration
  deltap:= evalf(delta + 2**(-prec)*(delta + Snorm)):
  Smax := Snorm + deltap:
od:
deltap, Smin, Smax:
end proc:



#---------------------------------------------------------------------
#  OBSOLETE use compute_horner_rounding_error below
# Computes the total relative rounding error
compute_rel_rounding_error:=proc(poly,xmax, n)
local deg, p, rho, deltap, Smin, Smax:

deg:=degree(poly):
if(n>0) then p:=100: else p:=53: fi:

if coeff(poly,x, 0) = 0 then
   deltap, Smin, Smax := compute_abs_rounding_error(poly/x,xmax, n-1):
   rho :=  (2^(-p)*(Smax+deltap) +deltap ) / Smin :
else
   deltap, Smin, Smax := compute_abs_rounding_error(poly,xmax, n):
   rho := deltap /  Smin:
fi:
rho:
end proc:




#---------------------------------------------------------------------
#  OBSOLETE use compute_horner_rounding_error below
# Computes the accumulated rounding error during the polynomial evaluation.
# P is the polynomial.
# xmax is the max value of |x|.
# n is the degree when P is computed in double double. The first double-double operation is a multiplication (probably less useful).

# returns max absolute error, min of the function, max of the function.

# This procedure also checks on the fly that the fast (test-free) versions of the double-double addition can be used, i.e. that for all x, at each Horner step i computing ci+x*Si, we have |ci|>|x*Si|. It prints warnings if it not the case.

compute_abs_rounding_error_firstmult:=proc(poly,xmax, nn)
local n, deg, delta, deltap, i, S, P, Snorm, Smin, Smax, prec:
deltap:=0:
delta:=0:
deg:=degree(poly):

prec:=53: # precision of the first iterations

S:=coeff(poly, x, deg):
Smax:=abs(S):
Smin:=Smax:

if nn<0 then n:=0: else n:=nn: fi:# sometimes called by compute_rel_rounding_error with n=-1

for i from (deg-1) to 0 by -1 do
  if i=n-1 then prec:=100: fi:  # from the mult of the n-1-th iteration
  P:= convert(S*x, polynom):
  Smin := abs(coeff(poly,x,i)) - xmax*Smax :
  if(Smin<=0) then
    printf("Warning! in compute_abs_rounding_error, Smin<=0 at iteration %d, consider decreasing xmax\n",i):
  fi:
  delta:= evalf(xmax*deltap + 2**(-prec)*xmax*Smax):
  if i<n then
    # fast Add22 ?
    if abs(coeff(poly,x,i)) < xmax*Smax  # may be improved to xmax*Smax/2
    then printf("WARNING Add22 cannot be used at step %d, use Add22Cond\n" , i  ):
         printf("    coeff=%1.20e,  xmax*Smax=%1.20e"  ,  abs(coeff(poly,x,i)),  xmax*Smax ):
    fi:
  fi:
  S:=convert(P+coeff(poly,x,i), polynom):
  Snorm:=evalf(numapprox[infnorm](S, x=-xmax..xmax)):
  deltap:= evalf(delta + 2**(-prec)*(delta + Snorm)):
  Smax := Snorm + deltap:
od:
deltap, Smin, Smax:
end proc:




#---------------------------------------------------------------------
#  OBSOLETE use compute_horner_rounding_error below
# Computes the total relative rounding error
compute_rel_rounding_error_firstmult:=proc(poly,xmax, n)
local deg, p, rho, deltap, Smin, Smax:

deg:=degree(poly):
if(n>0) then p:=100: else p:=53: fi:

if coeff(poly,x, 0) = 0 then
   deltap, Smin, Smax := compute_abs_rounding_error_firstmult(poly/x,xmax, n-1):
   rho :=  (2^(-p)*(Smax+deltap) +deltap ) / Smin :
else
   deltap, Smin, Smax := compute_abs_rounding_error_firstmult(poly,xmax, n):
   rho := deltap /  Smin:
fi:
rho:
end proc:




# Compute a good truncated polynomial approximation for a function
# Computes an approximation to a function of x f, as a truncated polynomial of deegree deg with the n first coefficients exactly representable as double-double.
# The function f(x) must have as input interval xmin..xmax
# returns [ truncated polynomial, relative approx error of trunc. poly.  ,  infinite precision polynomial,   rel. error of inf. prec. poly ]
poly_trunc_classic:=proc(f,deg,xmin,xmax,n)
  local pe, repe, pt, ppe, rept, maxpt:
  pe:=numapprox[minimax]( f,  x=xmin..xmax,  [deg,0],  1,  'err'):
  pt := poly_exact2(pe,n):
  rept := numapprox[infnorm]( 1-pt/f, x=xmin..xmax) :
  maxpt := numapprox[infnorm]( pt, x=xmin..xmax) :
  pt,rept, maxpt:
end proc:



#---------------------------------------------------------------------
# Computes a truncated polynomial of degree deg with the two first coefficients stored as double-doubles.
# The function f(x) must have as input interval xmin..xmax
# returns [ truncated polynomial, relative error of trunc. poly.  ,  infinite precision polynomial,   rel. error of inf. prec. poly ]
poly_trunc_f2d_2:=proc(f,deg,xmin,xmax)
  local pe, repe, pt, c0, c1, c2, ppe, abserr, relerr, maxpt, err:
  pe:=numapprox[minimax](  f,  x=xmin..xmax,  [deg,0],  1,  'err'):
  pt := poly_exact2(pe,2):
  c0:=coeff(pt,x,0):
  c1:=coeff(pt,x,1):
  c2:=coeff(pt,x,2):
  ppe:=numapprox[minimax](  (f - c0 - c1 * x - c2*x*x) ,  x=xmin..xmax,  [deg,0],  1,  'err'):
  ppe:=expand(ppe) - coeff(ppe,x,0) - coeff(ppe,x,1)*x- coeff(ppe,x,2)*x*x:
  pt := poly_exact2(ppe,2) + c0 + c1*x + c2*x*x:
  abserr := numapprox[infnorm]( pt-f, x=xmin..xmax) :
  relerr := numapprox[infnorm]( 1-pt/f, x=xmin..xmax) :
  maxpt := numapprox[infnorm]( pt, x=xmin..xmax) :
  pt,abserr,relerr,maxpt:
end proc:




#---------------------------------------------------------------------
# Computes a truncated polynomial of degree deg with the first coefficient stored as double-double.
#  The function f(x) must have as input interval xmin..xmax
# returns [ truncated polynomial, relative error of trunc. poly.  ,  infinite precision polynomial,   rel. error of inf. prec. poly ]
poly_trunc_f2d_1:=proc(f,deg,xmin,xmax)
  local pe, repe, pt, c0, c1, ppe, relerr, abserr, maxpt, err:
  pe:=numapprox[minimax](  f ,  x=xmin..xmax,  [deg,0],  1,  'err'):
  pt := poly_exact2(pe,1):
  c0:=coeff(pt,x,0):
  c1:=coeff(pt,x,1):
  ppe:=numapprox[minimax](  (f - c0 - c1 * x) ,  x=xmin..xmax,  [deg,0],  1,  'err'):
  ppe:=ppe - coeff(ppe,x,0) - coeff(ppe,x,1)*x:
  pt := poly_exact2(ppe,1) + c0 + c1*x:
  abserr := numapprox[infnorm]( pt-f, x=xmin..xmax) :
  relerr := numapprox[infnorm]( 1-pt/f, x=xmin..xmax) :
  maxpt := numapprox[infnorm]( pt, x=xmin..xmax) :
  pt,abserr,relerr,maxpt:
end proc:



#---------------------------------------------------------------------
#  compute_horner_rounding_error

# Computes a bound on the accumulated rounding error caused by the Horner evaluation of a truncated polynomial
# It is designed to allow evaluating the error for various schemes:
#   - with or without an error on x
#   - using SCS operators
#   - using double and double-double operators.
# Arguments:
#   P is the polynomial.
#   xmax is  the max value of |x|.
#   errors is a list of size n where n is the degree of the polynomial.
#      Each element of this list is a triple (epsx,epsmul,epsadd) where
#        epsx is the relative error on x at each step,
#        epsadd is the max relative error on the addition
#        epsmul is the max relative error on the multiplication.
#      This allows to handle SCS Horner, as well as evaluation starting in double and ending in double-double.
#   check_dd is a flag, if set to 1 the procedure also checks on the fly that the fast (test-free) versions
#       of the double-double addition can be used, i.e. that for all x, at each Horner step i computing ci+x*Si,
#       we have |ci|>|x*Si|. It prints warnings if it not the case.
# returns (epsprimek, deltak, minP, maxP) where (see the doc)
#   epsprimek is the max rel error of the last multiplication (useful if the coeff of degree 0 is 0)
#   deltak is the max absolute error of the last addition (useful if the reconstruction adds something to result of the polynomial)
#   minP  min of the evaluated polynomial (useful to compute a relative error out of deltak)
#   maxP  max of the evaluated polynomial.

compute_horner_rounding_error:=proc(poly, x, xmax, errors, check_dd)
local deg, Sk, maxSk, minSk, epsaddk, epsmulk, deltaaddk, k, ck, epsx, epsmul, epsadd, deltaadd, Pk, maxPk:

  if assigned(x) then
    printf("Error in compute_horner_rounding_error, polynomial variable is assigned\n"):
    return 'procname(args)':
  fi:

  deg:=degree(poly,x):
  if(deg<0) then  printf("ERROR: negative degree in compute_abs_rounding_error"): return 'procname(args)': fi:

  Sk:=coeff(poly, x, deg):
  maxSk:=abs(Sk):
  minSk:=maxSk:
  epsmulk:=0:
  deltaaddk:=0:
  epsaddk:=0:

  for k from (deg) to 1 by -1 do

    # the errors to consider for this step
    epsx := errors[k][1]:
    epsadd := errors[k][2]:
    epsmul := errors[k][3]:

    # multiplication operation
    Pk:= convert(Sk*x, polynom):
    maxPk:=numapprox[infnorm](Pk, x=-xmax..xmax):

    ck:=coeff(poly,x,k-1):
    epsmulk:=evalf( (1+epsx)*(1+epsaddk)*(1+epsmul)-1  + 10^(-Digits+2)   ):

    #addition
    if(ck=0)  then
      Sk:=Pk:
      maxSk := maxPk:
      minSk:=0:
      deltaaddk:= evalf(epsmulk*maxPk):
      epsaddk:=epsmulk:
    else
      Sk:=convert(Pk+ck , polynom):
      maxSk:=numapprox[infnorm](Sk, x=-xmax..xmax):
      minSk:=minimize(abs(Sk), x=-xmax..xmax):
      if(epsadd=2^(-53)) then   # compute deltadd exactly as the max half ulp of the result
	deltaadd := 0.5*ulp(maxSk+epsmulk*maxSk):
      else    # compute deltaadd out of the relative error
        deltaadd := epsadd * (maxSk+epsmulk*maxSk):
      fi:
      deltaaddk := evalf(  epsmulk*maxPk + deltaadd  + 10^(-Digits+2)  ):
      epsaddk := deltaaddk/minSk + 10^(-Digits+2) :
      # warnings
      if (minSk=0) then
        printf("Warning! in compute_abs_rounding_error, minSk=0 at iteration %d, consider decreasing xmax\n",k):
      fi:
    fi:
    printf("step %d   epsmulk=%1.4e  deltaaddk=%1.4e   minSk=%1.4e   maxSk=%1.4e\n", k, epsmulk, deltaaddk, minSk, maxSk):


#  if (epsadd=2**(-103)) then
#    # fast Add22 ?
#    if abs(coeff(poly,x,k)) < xmax*maxSk  # may be improved to xmax*Smax/2
#    then printf("WARNING Add22 cannot be used at step %d, use Add22Cond\n" , k  ):
#         printf("    coeff=%1.20e,  xmax*Smax=%1.20e"  ,  abs(coeff(poly,x,i)),  xmax*maxSk ):
#    fi:
#  fi:

od:

return (epsmulk, deltaaddk, minSk, maxSk)
end proc:


#---------------------------------------------------------------------
# An helper function to build an errlist for the previous procedure.
# Arguments:
#   n is the degree of the polynomial,
#   ddadd, ddmul number of (final) double-double operations
#   depsx  error on x in the double steps
#   ddepsx error on x in the double-double steps (when x is represented by a double-double)

errlist_quickphase_horner := proc(n,ddadd, ddmul,depsx, ddepsx)
 local nddadd, epsadd, nddmul, epsmul, epsx:
 if n=0
  then []
  else
    if ddadd>0 then
      nddadd:=ddadd-1:
      epsadd:=2**(-103):
    else
      nddadd:=ddadd:
      epsadd:=2**(-53):
    fi:
    if ddmul>0 then
      nddmul:=ddmul-1:
      epsmul:=2**(-102):
      epsx:=ddepsx:
    else
      nddmul:=ddmul:
      epsmul:=2**(-53):
      epsx:=depsx:
    fi:
    [  [epsx,epsadd,epsmul] ,   op(errlist_quickphase_horner(n-1, nddadd, nddmul, depsx,ddepsx))]
  fi:
end proc:



#---------------------------------------------------------------------
# Finding the worst cases for additive range reduction
# cut and paste from the web page of Muller's book
# Much faster if called with Digits very large and an evalf() on C
# but then check

WorstCaseForAdditiveRangeReduction:=proc(B,n,emin,emax,C)
  local epsilonmin,powerofBoverC,e,a,Plast,r,Qlast, Q,P,NewQ,NewP,epsilon, numbermin,expmin,l:
  epsilonmin := 12345.0 :
  powerofBoverC := B^(emin-n)/C:
  for e from emin-n+1 to emax-n+1 do
    powerofBoverC := B*powerofBoverC:
    a := floor(powerofBoverC):
    Plast := a:
    r := 1/(powerofBoverC-a):
    a := floor(r):
    Qlast := 1:
    Q := a:
    P := Plast*a+1:
    while Q < B^n-1 do
      r := 1/(r-a):
      a := floor(r):
      NewQ := Q*a+Qlast:
      NewP := P*a+Plast:
      Qlast := Q:
      Plast := P:
      Q := NewQ:
      P := NewP
      od:
   epsilon := evalf(C*abs(Plast-Qlast*powerofBoverC)):
   if epsilon < epsilonmin then
     epsilonmin := epsilon: numbermin := Qlast:
     expmin := e
     fi
   od:
  print('mantissa',numbermin):
  print('exponent',expmin):
  print('epsilon',epsilonmin):
  l := evalf(log(epsilonmin)/log(B),10):
  print(numberofdigits,l):
  (numbermin, expmin, epsilonmin)
end proc:







#####################################################################

# Stuff for SCS
#####################################################################

# Global parameters
# Don´t forget to set all the parameters here
SCS_NB_WORDS := 8:
SCS_NB_BITS  := 30:





#---------------------------------------------------------------------
# This procedure convert a decimal number into it SCS representation.
#        x : input number to convert into it SCS representation
real_to_SCS := proc(x)
        local exception, index, sgn, mantissa, nb, i:

            if x <> 0 then
                exception := 1:
                if x > 0 then
                    sgn  := 1:
                    nb    := x:
                elif x < 0 then
                    sgn := -1:
                    nb   := -x:
                end if:

                index := 0:

                if nb >= 1 then
                    for i from 0 while nb > (2^(SCS_NB_BITS+1)-1) do
                        index := index+1:
                        nb    := nb * 2^(-SCS_NB_BITS):
                    end do:
                else
                    for i from 0 while nb < 1 do
                        index := index-1:
                        nb    := nb * 2^(SCS_NB_BITS):
                    end do:
                end if:

                for i from 0 by 1 to (SCS_NB_WORDS-1) do
                    mantissa[i] := trunc(nb):
                    nb          := (nb - mantissa[i]) * 2^(SCS_NB_BITS):
                end do:
            else
                for i from 0 by 1 to (SCS_NB_WORDS-1) do
                    mantissa[i] := 0:
                end do:

                index     := 1:
                exception := x:
                sgn      := 1:
            end if:
            mantissa[SCS_NB_WORDS]   := exception:
            mantissa[SCS_NB_WORDS+1] := index:
            mantissa[SCS_NB_WORDS+2] := sgn:

            return mantissa:
        end proc:




#---------------------------------------------------------------------
# Convert an SCS number into a rational number

SCS_to_real := proc(tab)
       local res, i:

           if (tab[SCS_NB_WORDS] <> 1) then
               return tab[SCS_NB_WORDS]:
           end if:

           res := 0:
           for i from (SCS_NB_WORDS-1) by -1 while i>=0 do
               res := 2^(-SCS_NB_BITS)*res + tab[i]
           end do:

           res := tab[SCS_NB_WORDS+2]*(res * 2.^(SCS_NB_BITS * tab[SCS_NB_WORDS+1])):

           return res:

       end proc:



#---------------------------------------------------------------------
# This procedure truncates the coefficients of a polynomial to SCS
# numbers, so that we can then evaluate its approximation error
# (equivalent to poly_exact for the doubles)
poly_exact_SCS:=proc(P)
local deg,i, coef, coef_t, Q:
Q:= 0:
convert(Q, polynom):
deg:=degree(P,x):
  for i from 0 to deg do
    coef:=coeff(P,x,i):
    coef_t:=SCS_to_real(real_to_SCS(evalf(coef))):
    Q:= Q + coef_t*x^i:
  od:
return(Q):
end:



#---------------------------------------------------------------------
# Write Into file fd the SCSS number stored into the table tab where
# tab[0..(SCS_NB_WORDS-1)] store the mantissa
# tab[SCS_NB_WORDS] store the exception
# tab[SCS_NB_WORDS+1] store the index
# tab[SCS_NB_WORDS+2] store the sign
# You probably want to use WriteSCS below !

WriteSCS_from_table := proc(fd, tab)
              local i:

                  fprintf(fd,"{{"):

                  fprintf(fd,"0x%+0.8x, ", tab[0]):
                  for i from 1 by 1 to (SCS_NB_WORDS-2) do
                      fprintf(fd,"0x%+0.8x, ", tab[i]):
                      if (i mod 4 = 3) then
                          fprintf(fd,"\n"):
                      fi:
                  end do:
                  fprintf(fd,"0x%+0.8x},\n", tab[SCS_NB_WORDS-1]):
                  if (tab[SCS_NB_WORDS]=1) then
                      fprintf(fd,"DB_ONE, %3d, %3d ", tab[SCS_NB_WORDS+1], tab[SCS_NB_WORDS+2]):
                  else
                      # the only other possible value is 0 so ...
                      fprintf(fd,"{0x00000000, 0x00000000}, %3d, %3d ", tab[SCS_NB_WORDS+1], tab[SCS_NB_WORDS+2]):
                  end if:

                  fprintf(fd, "} \n"):
              end proc:


#---------------------------------------------------------------------
# Write a real number as an SCS array to a file

WriteSCS := proc (fd,x)
      WriteSCS_from_table(fd , real_to_SCS (x)):
      end:



#---------------------------------------------------------------------
# A procedure to count the non-zero coefficients of a polynomial (to store it)

get_nb_terms := proc(poly)
       local i, deg_poly:

           deg_poly := degree(poly):
           for i from deg_poly by -1 while i>=0 do
               if coeff(poly, x, i)=0 then
                   deg_poly := deg_poly-1:
               end if:
           end do:

           return deg_poly:
       end proc:








#---------------------------------------------------------------------
#   Write a polynomial to an array of SCS coefficients

# fd : file where to put the result
# poly : input polynom
# name : name of the array

Write_SCS_poly := proc(fd,  name, poly)
local i, deg:
           #fclose(fd):
           try
           finally
               fprintf(fd,"static const scs %s [%d]=\n", name, get_nb_terms(poly)+1):
               deg := degree(poly):

               fprintf(fd,"/* ~%1.50e */ \n{", coeff(poly, x, deg)):
               WriteSCS(fd, coeff(poly, x, deg)):
               for i from (deg-1) by (-1) while i>=0 do
                   if (coeff(poly, x, i)<>0) then
                       fprintf(fd,",\n/* ~%1.50e */ \n", coeff(poly, x, i)):
                       WriteSCS(fd, coeff(poly, x, i), 0):
                   end if:
               end do:
               fprintf(fd,"};\n"):
            end try:
          end proc:

















