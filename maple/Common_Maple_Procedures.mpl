# Evaluate this worksheet to produce Common_maple_procedures.m, which is used by the other worksheets in this directory
restart:
Digits:=150:
with(numapprox):
with(orthopoly):

#####################################################################
# Useful procedures for IEEE doubles


#---------------------------------------------------------------------

log2:=proc(x) evalf( log(x)/log(2) ) end proc:


#---------------------------------------------------------------------
# ieeedouble converts a number to IEEE double format. 
# returns sign (-1 or 1), exponent between -1022 and 1023, mantissa as a fraction between 0.5 and 1.
# Handles denorms well, infinities more or less

ieeedouble:=proc(xx) 
local x, sign, logabsx, exponent, mantissa, infmantissa; 
x:=evalf(xx):
if (x=0) then 
  sign,exponent,mantissa := 0,0,0; 
else 
 if (x<0) then sign:=-1:
 else sign:=1:
 fi:
 exponent := floor(log2(sign*x));
 if (exponent>1023) then mantissa:=infinity: exponent=1023:
 elif (exponent<-1022) then 
    # denorm
    exponent := -1023
 fi:
 infmantissa := sign*x*2^(52-exponent);
 if frac(infmantissa) <> 0.5 then mantissa := round(infmantissa)
 else
    mantissa := floor(infmantissa);
    if type(mantissa,odd) then mantissa := mantissa+1 fi;
 fi;
 Digits := 53;
 mantissa := mantissa*2^(-52);
fi;
sign,exponent,mantissa;
end:


#---------------------------------------------------------------------
# pulp returns the precision of the ulp of x

pulp:=proc(x)
local flt, ulpy;
flt:=ieeedouble(x);
ulpy:=-53+flt[2];
end proc:



#---------------------------------------------------------------------
# ulp returns the absolute value of the ulp of x
ulp:=proc(x)
2**(pulp(x));
end proc:
# Returns nearest IEEE double:
nearest := proc(x)
  local sign, exponent, mantissa;

  sign,exponent,mantissa := ieeedouble(x);
  sign*mantissa*2^(exponent);
end:


#---------------------------------------------------------------------
# ieehexa returns a string containing the hexadecimal representation of the double nearest to input x.

ieeehexa:= proc(x)
  local signe, hex1, hex2, ma, manti, expo, expos, bina, bin1, bin2, dec1, dec2, resultat;
  if(x=0) then resultat:=["00000000","00000000"];
  elif(x=-0) then resultat:=["80000000","00000000"];
  elif(x=2) then resultat:=["40000000","00000000"];
  elif(x=-2) then resultat:=["C0000000","00000000"];
  else
   ma:=ieeedouble(x);
   expo:=ma[2]:
   if (ma[1]=-1) then 
    manti:=2**64 + 2**63+(ma[3]-1)*2**52+(expo+1023)*2**52;
   else 
     manti:=2**64 + (ma[3]-1)*2**52+(expo+1023)*2**52;
   fi:
   hex2:=convert(manti, hex); 
   hex2:=convert(hex2, string):  
  
   resultat:=[substring(hex2,2..9), substring(hex2,10..18)];
  end if;
  resultat;
end proc:



#---------------------------------------------------------------------
# reciprocal of the previous
hexa2ieee:= proc(hexa)
local dec, bin, expobin, expo, mantis, sign, hex1, hex2, hexcat;
global res;

hex1:= op(1, hexa):
hex2:= op(2, hexa):
hexcat:= cat(hex1, hex2);
dec:= convert(hexcat, decimal, hex):

if(dec >= 2**63) then
  dec = dec - 2**63:
  sign:= -1:else
  sign:= 1:
fi;  
expo:= trunc(dec/(2**52))-1023:
mantis:= 1+frac(dec/(2**52));
res:= evalf(sign*2**(expo)*mantis);
end proc:
# Print a number x in Low or Big Endian representation in opened file "fd":
printendian:=proc(fd,x,isbig)
local xhl:
xhl:=ieeehexa(x):

if(isbig=0 or isbig=1) then
  fprintf(fd,"{{0x%+0.8s,0x%+0.8s}} /* %+0.10e */", xhl[2-isbig], xhl[isbig+1], x):
else
  print("ERROR, isbig must be equal to 0 or 1");
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
end if;
x_hi,x_lo;
end:



#---------------------------------------------------------------------
# same as hi_lo, but returns hexadecimal strings
ieeehexa2:=proc(x)
local reshi, reslo, hexhi, hexlo;
reshi:=nearest(x);
hexhi:=ieee2hexa(reshi);
reslo:= nearest(x-reshi);
hexlo:=ieee2hexa(reslo);
reshi, reslo;
end proc:




#---------------------------------------------------------------------
# Computes the constant for the round-to-nearest test. 
# delta is the overall relative error of the approximation scheme
compute_rn_constant := proc(delta)
  local k;
  k := trunc(-log2(delta)) - 53: 
  nearest(  1+ 2**(-52) + (2**(54+k)*delta)  /  ( (2**k-1) * (1-2**(-53)) )  ):
end proc:



#---------------------------------------------------------------------
# Takes a real number, and prints the bits after the 53th of its nearest IEEE floating-point number 

testWorstCaseRN:=proc(x)
  local xh,xl,s,e,m:
  xh:=nearest(x):
  xl:=x-xh;
  s,e,m := ieeedouble(xl):
  convert(op(1,m),binary);
end proc:



#####################################################################

# Stuff about truncated polynomials

# Truncate a polynomial 

#---------------------------------------------------------------------
# poly_exact takes a polynomial in x with arbitrary precision coefficients, and returns a truncated polynomial where coefficients are IEEE doubles.
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
return(Q);
end:


#---------------------------------------------------------------------
# Like poly_exact, but the n first coefficients are exactly representable as the sum of two doubles.
# (to actually get the two doubles, use procedure hi_lo)

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
  return(Q);
end:


#---------------------------------------------------------------------
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

prec:=53; # precision of the first iterations

S:=coeff(poly, x, deg):
Smax:=abs(S):
Smin:=Smax:

if nn<0 then n:=0: else n:=nn: fi:# sometimes called by compute_rel_rounding_error with n=-1

for i from (deg-1) to 0 by -1 do
  P:= convert(S*x, polynom):
  Smin := abs(coeff(poly,x,i)) - xmax*Smax : 
  if(Smin<=0) then 
    printf("Warning! in compute_abs_rounding_error, Smin<=0 at iteration %d, consider decreasing xmax\n",i);
  fi:
  delta:= evalf(xmax*deltap + 2**(-prec)*xmax*Smax):
  if i<n then 
    # fast Add22 ?    
    if abs(coeff(poly,x,i)) < xmax*Smax  # may be improved to xmax*Smax/2
    then printf("WARNING Add22 cannot be used at step %d, use Add22Cond\n" , i  );   
         printf("    coeff=%1.20e,  xmax*Smax=%1.20e"  ,  abs(coeff(poly,x,i)),  xmax*Smax );
    fi:
  fi:
  S:=convert(P+coeff(poly,x,i), polynom):
  Snorm:=evalf(infnorm(S, x=-xmax..xmax)):
  if i=n-1 then prec:=100: fi:  # from the addition of the n-1-th iteration
  deltap:= evalf(delta + 2**(-prec)*(delta + Snorm)): 
  Smax := Snorm + deltap:  
od:
deltap, Smin, Smax;
end proc:



#---------------------------------------------------------------------
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
rho;
end proc:




#---------------------------------------------------------------------
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

prec:=53; # precision of the first iterations

S:=coeff(poly, x, deg):
Smax:=abs(S):
Smin:=Smax:

if nn<0 then n:=0: else n:=nn: fi:# sometimes called by compute_rel_rounding_error with n=-1

for i from (deg-1) to 0 by -1 do
  if i=n-1 then prec:=100: fi:  # from the mult of the n-1-th iteration
  P:= convert(S*x, polynom):
  Smin := abs(coeff(poly,x,i)) - xmax*Smax : 
  if(Smin<=0) then 
    printf("Warning! in compute_abs_rounding_error, Smin<=0 at iteration %d, consider decreasing xmax\n",i);
  fi:
  delta:= evalf(xmax*deltap + 2**(-prec)*xmax*Smax):
  if i<n then 
    # fast Add22 ?    
    if abs(coeff(poly,x,i)) < xmax*Smax  # may be improved to xmax*Smax/2
    then printf("WARNING Add22 cannot be used at step %d, use Add22Cond\n" , i  );   
         printf("    coeff=%1.20e,  xmax*Smax=%1.20e"  ,  abs(coeff(poly,x,i)),  xmax*Smax );
    fi:
  fi:
  S:=convert(P+coeff(poly,x,i), polynom):
  Snorm:=evalf(infnorm(S, x=-xmax..xmax)):
  deltap:= evalf(delta + 2**(-prec)*(delta + Snorm)): 
  Smax := Snorm + deltap:  
od:
deltap, Smin, Smax;
end proc:




#---------------------------------------------------------------------
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
rho;
end proc:
# Compute a good truncated polynomial approximation for a function
# Computes an approximation to a function of x f, as a truncated polynomial of deegree deg with the n first coefficients exactly representable as double-double.
# The function f(x) must have as input interval xmin..xmax
# returns [ truncated polynomial, relative approx error of trunc. poly.  ,  infinite precision polynomial,   rel. error of inf. prec. poly ]
poly_trunc_classic:=proc(f,deg,xmin,xmax,n) 
  local pe, repe, pt, ppe, rept, maxpt;
  pe:=minimax( f,  x=xmin..xmax,  [deg,0],  1,  'err'):
  pt := poly_exact2(pe,n):
  rept := infnorm( 1-pt/f, x=xmin..xmax) :
  maxpt := infnorm( pt, x=xmin..xmax) :
  pt,rept, maxpt:
end proc:
# Computes a truncated polynomial of degree deg with the two first coefficients stored as double-doubles. 
# The function f(x) must have as input interval xmin..xmax
# returns [ truncated polynomial, relative error of trunc. poly.  ,  infinite precision polynomial,   rel. error of inf. prec. poly ]
poly_trunc_f2d_2:=proc(f,deg,xmin,xmax) 
  local pe, repe, pt, c0, c1, c2, ppe, abserr, relerr, maxpt, err;
  pe:=minimax(  f,  x=xmin..xmax,  [deg,0],  1,  'err'):
  pt := poly_exact2(pe,2):
  c0:=coeff(pt,x,0):
  c1:=coeff(pt,x,1):
  c2:=coeff(pt,x,2):
  ppe:=minimax(  (f - c0 - c1 * x - c2*x*x) ,  x=xmin..xmax,  [deg,0],  1,  'err'):
  ppe:=ppe - coeff(ppe,x,0) - coeff(ppe,x,1)*x- coeff(ppe,x,1)*x*x:
  pt := poly_exact2(ppe,2) + c0 + c1*x + c2*x*x:
  abserr := infnorm( pt-f, x=xmin..xmax) :
  relerr := infnorm( 1-pt/f, x=xmin..xmax) :
  maxpt := infnorm( pt, x=xmin..xmax) :
  pt,abserr,relerr,maxpt:
end proc:




#---------------------------------------------------------------------
# Computes a truncated polynomial of degree deg with the first coefficient stored as double-double.
#  The function f(x) must have as input interval xmin..xmax
# returns [ truncated polynomial, relative error of trunc. poly.  ,  infinite precision polynomial,   rel. error of inf. prec. poly ]
poly_trunc_f2d_1:=proc(f,deg,xmin,xmax) 
  local pe, repe, pt, c0, c1, ppe, relerr, abserr, maxpt, err;
  pe:=minimax(  f ,  x=xmin..xmax,  [deg,0],  1,  'err'):
  pt := poly_exact2(pe,1):
  c0:=coeff(pt,x,0):
  c1:=coeff(pt,x,1):
  ppe:=minimax(  (f - c0 - c1 * x) ,  x=xmin..xmax,  [deg,0],  1,  'err'):
  ppe:=ppe - coeff(ppe,x,0) - coeff(ppe,x,1)*x:
  pt := poly_exact2(ppe,1) + c0 + c1*x:
  abserr := infnorm( pt-f, x=xmin..xmax) :
  relerr := infnorm( 1-pt/f, x=xmin..xmax) :
  maxpt := infnorm( pt, x=xmin..xmax) :
  pt,abserr,relerr,maxpt:
end proc:



#---------------------------------------------------------------------
# Compute a bound on the accumulated rounding error caused by the Horner evaluation of a truncated polynomial 

# Compute a bound on the accumulated rounding error caused by the Horner evaluation of a truncated polynomial 
# Computes the accumulated rounding error during the polynomial evaluation.
# It is designed to allow evaluating the error for various schemes:
#   - with or without an error on x
#   - using SCS operators
#   - using double and double-double operators.
# Arguments:
# P is the polynomial.
# xmax is  the max value of |x|.
# errors is a list of size n where n is the degree of the polynomial. Each element of this list is a triple (xerr,mulerr,adderr) where xerr is the relative error on x at each step, adderr is the relative error on the addition and mulerr is a relative error on the multiplication. This allows to handle SCS Horner as well as evaluation starting in double and ending in double-double.
# check_dd is a flag, if set to 1 the procedure also checks on the fly that the fast (test-free) versions of the double-double addition can be used, i.e. that for all x, at each Horner step i computing ci+x*Si, we have |ci|>|x*Si|. It prints warnings if it not the case. 

# returns eps, epsp, max absolute error, min of the function, max of the function. 

compute_horner_rounding_error:=proc(poly, x, xmax, errors, check_dd)
local deg, Sk, maxSk, minSk, epsk, epsprimek, deltak, k, ck, xerr, mulerr, adderr, Pk, maxPk;
 
  if assigned(x) then 
    printf("Error in compute_horner_rounding_error, polynomial variable is assigned\n");
    return 'procname(args)';
  fi;

deg:=degree(poly,x):
if(deg<0) then  printf("ERROR: negative degree in compute_abs_rounding_error"); return 'procname(args)'; fi:

Sk:=coeff(poly, x, deg):
maxSk:=abs(Sk);
minSk:=maxSk;
epsprimek:=0;
deltak:=0;
epsk:=0;

for k from (deg) to 1 by -1 do

  # the errors to consider for this step
  xerr := errors[k][1];
  mulerr:=errors[k][2];
  adderr:=errors[k][3];

  # multiplication operation
  Pk:= convert(Sk*x, polynom):
  epsprimek:=evalf( (1+xerr)*(1+epsk)*(1+mulerr)-1  );

  #addition
  ck:=coeff(poly,x,k-1);
  Sk:=convert(Pk+ck , polynom);
  maxPk:=evalf(maximize(abs(Pk), x=-xmax..xmax)):
  maxSk:=evalf(maximize(abs(Sk), x=-xmax..xmax)):
  minSk:=evalf(minimize(abs(Sk), x=-xmax..xmax));
  if(ck=0) 
  then
    deltak:=0; 
    epsk:=epsprimek;
  else
    if(adderr=2^(-53)) then
      deltak := evalf(  epsprimek*maxPk + 0.5*ulp(maxSk+epsprimek*maxSk)   ):
    else
      deltak := evalf(  epsprimek*maxPk + 2^(-adderr) * (maxSk+epsprimek*maxSk)):
    fi:
    epsk := deltak/minSk;
  fi:
  printf("step %d   epsprimek=%1.4e  deltak=%1.4e   minSk=%1.4e   maxSk=%1.4e\n", 
          k, epsprimek, deltak, minSk, maxSk);

  # warnings
  if(minSk=0) then 
    printf("Warning! in compute_abs_rounding_error, minSk=0 at iteration %d, consider decreasing xmax\n",k);
  fi:

#  if (adderr=2**(-103)) then 
#    # fast Add22 ?    
#    if abs(coeff(poly,x,k)) < xmax*maxSk  # may be improved to xmax*Smax/2
#    then printf("WARNING Add22 cannot be used at step %d, use Add22Cond\n" , k  );   
#         printf("    coeff=%1.20e,  xmax*Smax=%1.20e"  ,  abs(coeff(poly,x,i)),  xmax*maxSk );
#    fi:
#  fi:

od:

return (epsprimek, deltak, minSk, maxSk)
end proc:


#---------------------------------------------------------------------
# An helper function to build an errlist for the previous procedure.
# Arguments: 
#   n is the degree of the polynomial,
#   ddadd, ddmul number of (final) double-double operations
#   dxerr  error on x in the double steps
#   ddxerr error on x in the double-double steps (when x is represented by a double-double)

errlist_quickphase_horner := proc(n,ddadd, ddmul,dxerr, ddxerr)
 local nddadd, adderr, nddmul, mulerr, xerr; 
 if n=0
  then []
  else
    if ddadd>0 then 
      nddadd:=ddadd-1:
      adderr:=2**(-103):
    else
      nddadd:=ddadd:
      adderr:=2**(-53):
    fi:
    if ddmul>0 then 
      nddmul:=ddmul-1:
      mulerr:=2**(-102):
      xerr:=ddxerr;
    else
      nddmul:=ddmul:
      mulerr:=2**(-53):
      xerr:=dxerr:
    fi:
    [  [xerr,mulerr,adderr] ,   op(errlist_quickphase_horner(n-1, nddadd, nddmul, dxerr,ddxerr))]
  fi:
end proc:



#---------------------------------------------------------------------
# Finding the worst cases for additive range reduction
# cut and paste from the web page of Muller's book
# Much faster if called with Digits very large and an evalf() on C

WorstCaseForAdditiveRangeReduction:=proc(B,n,emin,emax,C)
  local epsilonmin,powerofBoverC,e,a,Plast,r,Qlast, Q,P,NewQ,NewP,epsilon, numbermin,expmin,l;

  epsilonmin := 12345.0 ;
  powerofBoverC := B^(emin-n)/C;
  for e from emin-n+1 to emax-n+1 do
    powerofBoverC := B*powerofBoverC;
    a := floor(powerofBoverC);
    Plast := a;
    r := 1/(powerofBoverC-a);
    a := floor(r);
    Qlast := 1;
    Q := a;
    P := Plast*a+1;
    while Q < B^n-1 do
      r := 1/(r-a);
      a := floor(r);
      NewQ := Q*a+Qlast;
      NewP := P*a+Plast;
      Qlast := Q;
      Plast := P;
      Q := NewQ;
      P := NewP
      od;
   epsilon :=
   evalf(C*abs(Plast-Qlast*powerofBoverC));
   print(e+n-1);
   if epsilon < epsilonmin then
     epsilonmin := epsilon; numbermin := Qlast;
     expmin := e
     fi
   od;
  print('mantissa',numbermin);
  print('exponent',expmin);
  print('epsilon',epsilonmin);
  l := evalf(log(epsilonmin)/log(B),10);
  print(numberofdigits,l)
end proc:







#####################################################################

# Stuff for SCS
# All the fllowing may be outdated

# Global parameters
# Don´t forget to set all the parameters here
NB_WORDS := 2:
NB_BITS  := 30:

# GetSCS_real
# This procedure convert a decimal number into it SCS representation.
#        x : input number to convert into it SCS representation
GetSCS_real := proc(x)
local exception, index, sign, mantissa, nb, i;

if x <> 0 then
  exception := 1;
  if x > 0 then
    sign  := 1;
    nb    := x;
  elif x < 0 then
    sign := -1;
    nb   := -x;
  end if;
  
  index := 0;

  if nb >= 1 then
    for i from 0 while nb > (2^(NB_BITS+1)-1) do
      index := index+1;
      nb    := nb * 2^(-NB_BITS);
    end do; 
  else
    for i from 0 while nb < 1 do
      index := index-1;
      nb    := nb * 2^(NB_BITS);
    end do; 
  end if;

  for i from 0 by 1 to (NB_WORDS-1) do
    mantissa[i] := trunc(nb);
    nb          := (nb - mantissa[i]) * 2^(NB_BITS);
  end do;
else
  for i from 0 by 1 to (NB_WORDS-1) do
    mantissa[i] := 0;
  end do;

  index     := 1;
  exception := x;
  sign      := 1;
end if;
mantissa[NB_WORDS]   := exception;
mantissa[NB_WORDS+1] := index;
mantissa[NB_WORDS+2] := sign;

return mantissa;
end proc:


# SetSCS_real
# Convert an SCS number into a rational number

SetSCS_real := proc(tab)
  local res, i;

  if (tab[NB_WORDS] <> 1) then
    return tab[NB_WORDS];
  end if;

  res := 0;
  for i from (NB_WORDS-1) by -1 while i>=0 do
    res := 2^(-NB_BITS)*res + tab[i] 
  end do;
  
res := tab[NB_WORDS+2]*(res * 2.^(NB_BITS * tab[NB_WORDS+1]));

return res;

end proc:


# WriteSCS
# Write Into file fd the SCSS number stored into the table tab where
# tab[0..(NB_WORDS-1)] store the mantissa
# tab[NB_WORDS] store the exception
# tab[NB_WORDS+1] store the index
# tab[NB_WORDS+2] store the sign

WriteSCS := proc(fd, tab)
 local i;

fprintf(fd,"{{");

fprintf(fd,"0x%+0.16x, ", tab[0]);
for i from 1 by 1 to (NB_WORDS-2) do
  fprintf(fd,"0x%+0.16x, ", tab[i]);
  if (i mod 4 = 3) then
    fprintf(fd,"\n"); 
  fi;
end do;
fprintf(fd,"0x%+0.16x},\n", tab[NB_WORDS-1]);
if (tab[NB_WORDS]=1) then
  fprintf(fd,"DB_ONE, %3d, %3d ", tab[NB_WORDS+1], tab[NB_WORDS+2]);
else
  # the only other possible value is 0 so ...
  fprintf(fd,"{0x00000000, 0x00000000}, %3d, %3d ", tab[NB_WORDS+1], tab[NB_WORDS+2]);
end if;

fprintf(fd, "} \n");
end proc:
# GetSCS_poly

get_nb_terms := proc(poly)
 local i, deg_poly;

 deg_poly := degree(poly);
 for i from deg_poly by -1 while i>=0 do
  if coeff(poly, x, i)=0 then
   deg_poly := deg_poly-1;
  end if;
 end do;

 return deg_poly;
end proc:

# Convert each coefficient of a polynom into it SCSS representation


# poly : input polynom 
#   file : name of the file where to put the result
# GetSCS_poly := proc(poly, file)
#   local i, fd, mantissa, deg;
#   #fclose(fd);
#   try
#     fd := fopen(file, WRITE,TEXT);
#   finally
#     fprintf(fd,"static const SCS constant_poly[%d]=\n",get_nb_terms(poly)+1);
#     deg := degree(poly); 

#     fprintf(fd,"/* ~%e */ \n{", coeff(poly, x, deg));   
#     mantissa := GetSCS_real(coeff(poly, x, deg));
#     WriteSCS(fd, mantissa); 
#     for i from (deg-1) by (-1) while i>=0 do
#       if (coeff(poly, x, i)<>0) then 
#         fprintf(fd,",\n/* ~%e */ \n", coeff(poly, x, i));
#         mantissa := GetSCS_real(coeff(poly, x, i));
#         WriteSCS(fd, mantissa, 0);
#       end if; 
#     end do;
#     fprintf(fd,"};\n");
#    
#   fclose(fd);
#   end try;
# end proc:


# hexa_scs

# Une dernière procédure pour écrire les nombres SCS (floatant exacts) en hexadécimal dans un fichier ".h"

# hexa_scs := proc (x) 
#   local i, elem, mantisse, exposant, resultat; 
#   global ulp, ulp_inv;
#   if (x :: list) then
#     [seq (hexa(elem), elem=x)];
#   else
#     Digits := 60;
#     if (x = 0) then 
#       resultat := 0;
#     else
#       exposant := floor(log[2](abs(evalf(x))));
#       if (exposant < -126) then exposant := -126; fi;
#       resultat := round (ulp_inv * (abs(x) * 2^(-exposant) -1)) + ulp_inv / 2 * (exposant + 2^(lgex-1) - 1);
#       if (x < 0) then resultat := resultat + ulp_inv * 2^(lgex-1); fi; 
#     fi;
#   fi;
#   convert(resultat, hex);
# end:

















