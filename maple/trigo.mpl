
Digits := 100:
interface(quiet=true):
read "common-procedures.mpl":


mkdir("TEMPTRIG"):


#################################
#Comments and todos :

#don't forget the error on k in Cody and Waite
# recheck the max_return_* : we shave off the lower bits, etc

# - Evaluation scheme :
# case 1 : return x
# case 2 (or Fast) : compute a simple polynomial
# case 3 : do an argument reduction...



########################################################
#  Case 1 : Small arguments
# return x for sine and tan, return 1 for cos
########################################################

xmax_return_x_for_sin := 2^(-26):
xmax_return_1_for_cos_RN := sqrt(2^(-53)):
xmax_return_1_for_cos_RDIR:=2^(-26):
one_rounded_down := evalf(1-ulp(1/2)):

xmax_return_x_for_tan := 2^(-26):




########################################################
# Case 2 : simple polynomial approximation
########################################################

# We want to use the same polynomial in case 2 and 3.
# So see after arg red

#################################################
#   Case 3 : Argument reduction
#################################################



# TODO  This value seems to work but needs proving. Besides maybe it
# should be the same for sine and cos

maxepstotalSinCase3:=2**(-65);
rnconstantSinCase3 := evalf(compute_rn_constant(maxepstotalSinCase3)):

maxepstotalCosCase3:=2**(-65);
rnconstantCosCase3 := evalf(compute_rn_constant(maxepstotalCosCase3)):


#The following is not finished

#################################################
# CODY and WAITE  Argument reduction

# TODO add the error due to rounding k


C := Pi/256;
invC:= nearest(1/C);
reminvC := evalf(1/C - invC);
expC:=ieeedouble(C)[2]:

# There are three sets of constants :
#  - split redC into two constants, for small values when we are concerned with absolute error
#  - split redC into three constants, for larger values
#  - split redC into three doubles, for the cases when we need
#     good relative precision on the result and fear cancellation




# Case2est reduction using two-part Cody and Waite (up to |k|=2^22)

bitsCh_0:=32:  # ensures at least 53+11 bits

# 1/2 <= C/2^(expC+1) <1
Ch:= round(evalf(  C * 2^(bitsCh_0-expC-1))) / (2^(bitsCh_0-expC-1)):
# recompute bitsCh in case we are lucky (and we are for bitsCh_0=32)
bitsCh:=1+log2(op(2,ieeedouble(Ch)[3])) :  # this means the log of the denominator

Cl:=nearest(C - Ch):
# Cody and Waite argument reduction will work for |k|<kmax_cw2
kmax_cw2:=2^(53-bitsCh):

# The constants to move to the .h file
RR_CW2_CH := Ch:
RR_CW2_MCL := -Cl:
XMAX_CODY_WAITE_2 := nearest(kmax_cw2*C):

# The error in this case (we need absolute error)
delta_repr_C_cw2   := abs(C-Ch-Cl);
delta_round_cw2    := kmax_cw2* 1/2 * ulp(Cl) ;
delta_cody_waite_2 := kmax_cw2 * delta_repr_C_cw2 + delta_round_cw2;
# This is the delta on y, the reduced argument

log2(%);




# Slower reduction using three-part Cody and Waite, up to |k|=2^31

bitsCh_0:=21:
Ch:= round(evalf(  C * 2^(bitsCh_0-expC-1))) / (2^(bitsCh_0-expC-1)):
# recompute bitsCh in case we are lucky (and we are for bitsCh_0=32)
bitsCh:=1+log2(op(2,ieeedouble(Ch)[3])) :  # this means the log of the denominator

r := C-Ch:
Cmed := round(evalf(  r * 2^(2*bitsCh-expC-1))) / (2^(2*bitsCh-expC-1)):
bitsCmed:=1+log2(op(2,ieeedouble(Cmed)[3])) :

Cl:=nearest(C - Ch - Cmed):

kmax_cw3 := 2^31:# Otherwise we have integer overflow



# The constants to move to the .h file
RR_CW3_CH  := Ch;
RR_CW3_CM  := Cmed:
RR_CW3_MCL := -Cl:
XMAX_CODY_WAITE_3 := nearest(kmax_cw3*C):

# The error in this case (we need absolute error)
delta_repr_C_cw3   := abs(C - Ch - Cmed - Cl):
delta_round_cw3    := kmax_cw3 * 1/2 * ulp(Cl) :
delta_cody_waite_3 := kmax_cw3 * delta_repr_C_cw3 + delta_round_cw3:
# This is the delta on y, the reduced argument

log2(%);





# Third range reduction, using double-double arithmetic, for |k| up to 2^51-1

# max int value that we can be produced by DOUBLE2LONGINT
kmax:=2^51-1:
XMAX_DDRR:=nearest(kmax*C);

#in this case we have C stored as 3 doubles
Ch   := nearest(C):
Cmed := nearest(C-Ch):
Cl   := nearest(C-Ch-Cmed):

RR_DD_MCH := -Ch:
RR_DD_MCM := -Cmed:
RR_DD_CL  := Cl:

delta_repr_C := abs(C - Ch - Cmed - Cl):

# and we have only exact Add12 and Mul12  operations. The only place
# with possible rounding errors is:
#       Add22 (pyh, pyl,    (x + kch_h) , (kcm_l - kd*RR_DD_CL),   th, tl) ;
# where (x + kch_h) is exact (Sterbenz) with up to kmax bits of cancellation
# and the error is simply the error in  (kcm_l - kd*RR_DD_CL)
# At the very worst :
delta_round :=
              kmax * 1/2 * ulp(Cl) # for   kd*RR_DD_CL
              + kmax*ulp(Cl)         # for the subtraction
              + 2^(-100) * Pi/512 :    # for the Add22
delta_RR_DD :=  kmax * delta_repr_C + delta_round:

#  the last case, Payne and Hanek reduction, gives a very small delta:
#  red arg is on 9*30 bits, then rounded to a double-double (106 bits)
# This should, of course, be optimized some day
delta_PayneHanek := 2^(-100):

# Finally the max delta on the reduced argument is
delta_ArgRed := max(delta_cody_waite_2, delta_cody_waite_3,
                    delta_RR_DD, delta_PayneHanek):
log2(delta_ArgRed);





# Now we use the above range reductions when k mod 256 <> 0
# otherwise we need to worry about relative accuracy of the result.

# First, what is the worst case for cancellation ?

emax := ieeedouble(XMAX_DDRR)[2] +1 :
# above emax, we will use Payne and Hanek so we do not worry

(wcn, wce, wceps) := WorstCaseForAdditiveRangeReduction(2,53,-8, emax, C):
wcx := wcn * 2^wce:
wck := round(wcx/C):
wcy := wcx - wck*C:

log2(wcy);   # y > 2^(-67);

# In these cases we use the double-double range reduction, for |k|<kmax_cw3
# and the relative precision in the worst case is for wcy

delta_round := kmax_cw3 * 1/2 * ulp(Cl)      # for   kd*RR_DD_CL
              + kmax_cw3 * ulp(Cl) :         # for the subtraction

delta_RR_DD :=  kmax_cw3 * delta_repr_C + delta_round:

eps_ArgRed := (1+delta_RR_DD/wcy)*(1+2^(-100)) -1:

log2(eps_ArgRed);

# In all the other cases we use Payne and Hanek, and eps_ArgRed is
# much smaller, so this is the max.











######################


# An attempt to use Chebychev polynomials, less accurate in our case after poly_exact2
# If we need extra precision for the polynomial evaluation use Chebpade as follows :
#with(orthopoly):
#Poly_P := series(sin(sqrt(x))/(x^(3/2))-1/x, x=0, degreeSinCase2*4):
#Poly_Q := convert(Poly_P,polynom):
#Poly_cheb := numapprox[chebpade](Poly_Q, x=0..xmaxSinCase2^2, [degreeSinCase2/2-2,0]):
#polySinCase2 := poly_exact2 ( expand(x + x^3 * subs(x=x^2, Poly_cheb)), 2);
#


###########
# Polynomials for do_sine and do_cos, and for the case 2
degreeSin := 8:
degreeCos := 7:
ymaxCase3  := Pi/512:
y2maxCase3 := ymaxCase3^2:
# These are the parameters to vary  (they should always be larger than Pi/512
xmaxCosCase2   := Pi/512:
xmaxSinCase2   := Pi/256:
x2maxSinCase2:= xmaxSinCase2^2:
x2maxCosCase2:= xmaxCosCase2^2:

# The difficulty here is to find polynomials which are good for case 2
# as well as for case 3. A simple solution is to set xmaxCosCase2 =ymaxCase3...


# For the sine: simple approach using Taylor is better than minimax
if(1+1=2) then
polySin:=  poly_exact(convert( series(sin(x), x=0, degreeSin+1), polynom)):
polyTs := expand(polySin/x-1):
polyTs2 := subs(x=sqrt(y), polyTs):
else
# More accurate. Compute the minimax up to the case 2 limit
# THIS MINIMAX MINIMIZES ABS ERROR AND NOT REL ERROR
polyTs2x := poly_exact(numapprox[minimax]((sin(sqrt(x))/sqrt(x)), x=2^(-2048)..x2maxSinCase2, [3,0]) ) -1:
polyTs2 := subs(x=y, polyTs2x):
polySin := expand(x*(1+ subs(y=x^2, polyTs2))):
end if:

# For the cos: compute a minimax
if(1+1=3) then
# simple approach using Taylor
polyCos  := poly_exact (convert( series(cos(x), x=0, degreeCos+1), polynom)):
polyTc2 := subs(x=sqrt(y), polyCos - 1):
else
# More accurate. Compute the minimax up to the case 2 limit
polyTc2x := poly_exact(numapprox[minimax](cos(sqrt(x)), x=2^(-2048)..x2maxCosCase2, [3,0]) ) -1:
polyTc2 := subs(x=y,polyTc2x):
polyCos := expand(1+ subs(y=x^2, polyTc2)):
end if:


eps_approx_Sin_Case2 := numapprox[infnorm]((x*polyTs+x -sin(x))/sin(x), x=0..xmaxSinCase2):
log2(%);
eps_approx_Sin_Case3 := numapprox[infnorm]((x*polyTs +x -sin(x))/sin(x), x=0..ymaxCase3):
log2(%);

delta_approx_Tc_Case2:= numapprox[infnorm](polyCos -  cos(x), x=0..xmaxCosCase2):
log2(%);
delta_approx_Tc_Case3:= numapprox[infnorm](polyCos -  cos(x), x=0..ymaxCase3):
log2(%);



########################## Case 2 for sine  ###########################



# evaluate this polynomial in double. The error on x*x is at most half an ulp
errlist:=errlist_quickphase_horner(degree(polyTs2),0,0,2^(-53), 0):
(eps_rounding_Ts, delta_rounding_Ts, minTs, maxTs):=
	compute_horner_rounding_error(polyTs2,y,x2maxSinCase2, errlist, true):

eps_poly_Ts_Case2 := numapprox[infnorm]((x*polyTs)/(sin(x)-x) -1, x=0..xmaxSinCase2):
maxeps2 := (1+eps_poly_Ts_Case2)*(1+eps_rounding_Ts)*(1+2^(-53))-1:

maxepstotalSinCase2 := maxeps2 * numapprox[infnorm](1-x/sin(x), x=0..xmaxSinCase2);
rnconstantSinCase2 := evalf(compute_rn_constant(maxepstotalSinCase2));





##################################### Case2 cos ###########################

# evaluate this polynomial in double. The error on x*x is at most half an ulp
errlist        := errlist_quickphase_horner(degree(polyTc2),0,0,2**(-53), 0):
(eps_rounding_Tc, delta_rounding_Tc, minTc, maxTc):=
              compute_horner_rounding_error(polyTc2,y,x2maxCosCase2, errlist, true):

# Then we have an Add12 which is exact. The result is greater then cos(xmaxCosCase2):
miny := cos(xmaxCosCase2);
maxepstotalCosCase2 :=  (delta_rounding_Tc + delta_approx_Tc_Case2) / miny ;
log2(%);
rnconstantCosCase2 := evalf(compute_rn_constant(maxepstotalCosCase2));




######################## Case2 Tangent #########################
#
xmaxTanCase2   := 2**(-3);
degreeTanCase2 := 16;
rnconstantTanCase21 := evalf(compute_rn_constant(2**(-62)));
rnconstantTanCase22 := evalf(compute_rn_constant(2**(-59)));

# Compute the Taylor series
with(orthopoly):
Poly_P := convert(series(tan(sqrt(x))/(x^(3/2))-1/x, x=0, degreeTanCase2*4),polynom):
Poly_cheb := numapprox[chebpade](Poly_P, x=0..xmaxTanCase2^2, [degreeTanCase2/2-2,0]):
polyTanCase2 :=  poly_exact2(expand(x + x^3 * subs(x=x^2, Poly_cheb)), 4);



approx_errorTanCase2:=numapprox[infnorm](1 - polyTanCase2 / tan(x), x=0..xmaxTanCase2):
log2(approx_errorTanCase2);




###############################################################################
#   Computing errors for Case3 : now we have an error due to arg red


# The error on yh*yh
# we had $y_h+y_l = y + \abserr{CodyWaite}$.
# Now we take only $y_h$ :  $y_h = (y_h + y_l)(1+2^{-53}) = (y+\abserr{CodyWaite})(1+2^{-53})$
# When squared we get $y_h\otimes y_h = (y+\abserr{CodyWaite})^2(1+2^{-53})^3 $

epsy2 := evalf(((1+eps_ArgRed)**2) *  (1+2**(-53))**2 - 1):

############### Errors in computing Ts and Tc


errlist:=errlist_quickphase_horner(degree(polyTs2),0,0, epsy2 , 0):
(eps_rounding_Ts, delta_rounding_Ts, minTs, maxTs):=
	compute_horner_rounding_error(polyTs2,y,y2maxCase3, errlist, true):

errlist:=errlist_quickphase_horner(degree(polyTc2),0,0, epsy2 , 0):
(eps_rounding_Tc, delta_rounding_Tc, minTc, maxTc):=
              compute_horner_rounding_error(polyTc2,y,y2maxCase3, errlist, true):



##############   Case  k&127 = 0
# See above. As the error will always be smaller than in the case
# k&127 != 0, no need to recompute it. It could be useful to have a
# better rounding constant, but it will be statistically
# unsignificant.





###### The extreme Cases for tabulated values for the Case k&127 != 0
minsca := sin(Pi/256);
minscah:= nearest(minsca);
maxsca := cos(Pi/256);
maxscah:= nearest(maxsca);

# Worst Case of approximation error

################ Summing everything up
# The only error is in
#  tlo =  tc*cah - (ts*sahyh_h -  (cal + (tlo  - (sahyh_l + (sal*yh + sah*yl)) )));

delta_round_tlo := 2*ulp(maxTc*maxsca) ; # TODO c'est à la louche


delta_do_cos := delta_round_tlo;
min_sin := sin(Pi/512); # TODO Add error on k here

epsilon_sincos := delta_round_tlo / min_sin; # TODO hum

rnconstant_sincos := compute_rn_constant(epsilon_sincos);





##############################################
## Compute constants for SCS arg red
oldDigits:=Digits;
Digits:=1000;
# for 2/Pi:
n:=round(2^(30*48)*evalf(2/Pi));
digitlist:=[]:
for i from 1 to 48 do
  r:=n mod (2^30):
  n:=floor(n/(2^30)):
  hexstring:= convert(convert(r,hex),string):
  digitlist:=[hexstring, op(digitlist)]:
end:
digitlist;

# for 256/Pi:
n:=round(2^(30*47)*evalf(256/Pi));
digitlist:=[]:
for i from 1 to 48 do
  r:=n mod (2^30):
  n:=floor(n/(2^30)):
  hexstring:= convert(convert(r,hex),string):
  digitlist:=[hexstring, op(digitlist)]:
end:
digitlist;
Digits:=oldDigits;



# an auxiliary output function:
# Outputs the high part of a double, and the double in comment.
# As all these high parts are used in code as
# if(absxhi < XMAX_COS_CASE2)
# we have to remove one LSB to the high part, or, divide var by
# (1+2^(-20))
# Now we have absxhi<highpart(var/(1+2^(-20))
# => absxhi*(1+2^(-20))<var
# => x<var
outputHighPart:=proc(cvarname, var)
  local varMinusLSB:
  Digits:=8:
  varMinusLSB:=var/(1+2^(-20)):
  ("#define " || cvarname || " 0x" || (ieeehexa(varMinusLSB)[1])
    ||  "        /* " || (convert(evalf(varMinusLSB),string)) ||  " */" )
end proc:





# Output:

filename:="TEMPTRIG/trigo_fast.h":
fd:=fopen(filename, WRITE, TEXT):

fprintf(fd, "#include \"crlibm.h\"\n#include \"crlibm_private.h\"\n"):

fprintf(fd, "\n/*File generated by maple/coef_sine.mw*/\n"):
fprintf(fd, "\n"):


fprintf(fd,  "%s\n",  outputHighPart("XMAX_RETURN_X_FOR_SIN", xmax_return_x_for_sin) ):
fprintf(fd,  "%s\n",  outputHighPart("XMAX_SIN_CASE2        ", xmaxSinCase2) ):

fprintf(fd,  "%s\n",  outputHighPart("XMAX_RETURN_1_FOR_COS_RN", xmax_return_1_for_cos_RN) ):
fprintf(fd,  "%s\n",  outputHighPart("XMAX_RETURN_1_FOR_COS_RDIR", xmax_return_1_for_cos_RDIR) ):
fprintf(fd,  "%s\n",  outputHighPart("XMAX_COS_CASE2        ", xmaxCosCase2) ):

fprintf(fd,  "%s\n",  outputHighPart("XMAX_RETURN_X_FOR_TAN", xmax_return_x_for_tan) ):
fprintf(fd,  "%s\n",  outputHighPart("XMAX_TAN_CASE2        ", xmaxTanCase2) ):

fprintf(fd, "\n"):

fprintf(fd, "#define EPS_SIN_CASE2     %1.25e \n", maxepstotalSinCase2):
fprintf(fd, "#define RN_CST_SIN_CASE2  %1.25f \n", rnconstantSinCase2):
fprintf(fd, "#define EPS_SIN_CASE3     %1.25e \n", maxepstotalSinCase3):
fprintf(fd, "#define RN_CST_SIN_CASE3  %1.25f \n", rnconstantSinCase3):

fprintf(fd, "#define EPS_COS_CASE2     %1.25e \n", maxepstotalCosCase2):
fprintf(fd, "#define RN_CST_COS_CASE2  %1.25f \n", rnconstantCosCase2):
fprintf(fd, "#define EPS_COS_CASE3     %1.25e \n", maxepstotalCosCase3):
fprintf(fd, "#define RN_CST_COS_CASE3  %1.25f \n", rnconstantCosCase3):
fprintf(fd, "#define ONE_ROUNDED_DOWN  %1.25e \n", one_rounded_down):

fprintf(fd, "#define RN_CST_TAN_CASE21 %1.25e \n", rnconstantTanCase21):
fprintf(fd, "#define RN_CST_TAN_CASE22 %1.25e \n", rnconstantTanCase22):

fprintf(fd, "\n"):

fprintf(fd, "#define INV_PIO256        %1.25f \n", 1/C):
fprintf(fd, "\n"):

fprintf(fd,  "%s\n",  outputHighPart("XMAX_CODY_WAITE_2", XMAX_CODY_WAITE_2) ):
fprintf(fd,  "%s\n",  outputHighPart("XMAX_CODY_WAITE_3", XMAX_CODY_WAITE_3) ):
fprintf(fd,  "%s\n",  outputHighPart("XMAX_DDRR        ", XMAX_DDRR) ):
#fprintf(fd,  "%s\n",  outputHighPart("", ) ):
fprintf(fd, "\n"):

fprintf(fd, "#define RR_CW2_CH  %1.25e\n", RR_CW2_CH):
fprintf(fd, "#define RR_CW2_MCL %1.25e\n", RR_CW2_MCL):
fprintf(fd, "\n"):

fprintf(fd, "#define RR_CW3_CH  %1.25e\n", RR_CW3_CH):
fprintf(fd, "#define RR_CW3_CM  %1.25e\n", RR_CW3_CM):
fprintf(fd, "#define RR_CW3_MCL %1.25e\n", RR_CW3_MCL):
fprintf(fd, "\n"):

fprintf(fd, "#define RR_DD_MCH  %1.25e\n", RR_DD_MCH):
fprintf(fd, "#define RR_DD_MCM  %1.25e\n", RR_DD_MCM):
fprintf(fd, "#define RR_DD_CL   %1.25e\n", RR_DD_CL):
fprintf(fd, "\n"):

fprintf(fd, "\n"):

  fprintf(fd, "\n\n"):

  # The 256/Pi SCS array
  fprintf(fd, "static const int digits_256_over_pi[] = \n{"):
  for i from 0 to 11 do
    for j from 1 to 4 do
      fprintf(fd, " 0x%s,  \t",digitlist[4*i+j]):
    end:
    fprintf(fd, "\n "):
  end:
  fprintf(fd, "};\n\n"):

  # The Pi/256 SCS constant
  fprintf(fd, "static const scs Pio256=\n"):
  WriteSCS(fd, evalf(C)):
  fprintf(fd, ";\n#define Pio256_ptr  (scs_ptr)(& Pio256)\n\n"):

fprintf(fd,"#ifdef WORDS_BIGENDIAN\n"):
for isbig from 1 to 0 by -1 do

  if(isbig=0) then
    fprintf(fd,"#else\n"):
  fi:

  # The sine polynomial

  fprintf(fd, "static db_number const s3 = "):
  printendian(fd, coeff(polySin,x,3), isbig):
  fprintf(fd, ";\n"):
  fprintf(fd, "static db_number const s5 = "):
  printendian(fd, coeff(polySin,x,5), isbig):
  fprintf(fd, ";\n"):
  fprintf(fd, "static db_number const s7 = "):
  printendian(fd, coeff(polySin,x,7), isbig):
  fprintf(fd, ";\n\n"):


  # the cos polynomial

  fprintf(fd, "static db_number const c2 = "):
  printendian(fd, coeff(polyCos,x,2), isbig):
  fprintf(fd, ";\n"):
  fprintf(fd, "static db_number const c4 = "):
  printendian(fd, coeff(polyCos,x,4), isbig):
  fprintf(fd, ";\n"):
  fprintf(fd, "static db_number const c6 = "):
  printendian(fd, coeff(polyCos,x,6), isbig):
  fprintf(fd, ";\n\n"):


  # the tan polynomial

  t3h, t3l := hi_lo(coeff(polyTanCase2,x,3)):
  fprintf(fd, "static db_number const t3h = "):
  printendian(fd, t3h, isbig):
  fprintf(fd, ";\n"):
  fprintf(fd, "static db_number const t3l = "):
  printendian(fd, t3l, isbig):
  fprintf(fd, ";\n"):
  fprintf(fd, "static db_number const t5 = "):
  printendian(fd, coeff(polyTanCase2,x,5), isbig):
  fprintf(fd, ";\n"):
  fprintf(fd, "static db_number const t7 = "):
  printendian(fd, coeff(polyTanCase2,x,7), isbig):
  fprintf(fd, ";\n"):
  fprintf(fd, "static db_number const t9 = "):
  printendian(fd, coeff(polyTanCase2,x,9), isbig):
  fprintf(fd, ";\n"):
  fprintf(fd, "static db_number const t11 = "):
  printendian(fd, coeff(polyTanCase2,x,11), isbig):
  fprintf(fd, ";\n"):
  fprintf(fd, "static db_number const t13 = "):
  printendian(fd, coeff(polyTanCase2,x,13), isbig):
  fprintf(fd, ";\n"):
  fprintf(fd, "static db_number const t15 = "):
  printendian(fd, coeff(polyTanCase2,x,15), isbig):
  fprintf(fd, ";\n\n"):


  # The sincos table

  fprintf(fd, "\n/*  sine and cos of kPi/256 in double-double */\n"):
  SinCosSize:= 128:
  fprintf(fd, "static db_number const sincosTable[%d] =\n{\n",  4*(SinCosSize/2+1)):
  for i from 0 to SinCosSize/2 do
      s:=hi_lo(sin(i*Pi/(2*SinCosSize))):
      c:=hi_lo(cos(i*Pi/(2*SinCosSize))):
      printendian(fd,s[1],isbig):
      fprintf(fd," ,\n"):
      printendian(fd,s[2],isbig):
      fprintf(fd," ,\n"):
      printendian(fd,c[1],isbig):
      fprintf(fd," ,\n"):
      printendian(fd,c[2],isbig):

    if i<SinCosSize-1 then fprintf(fd," ,\n"): fi:
  od:
  fprintf(fd, "\n};\n\n"):


od:
fprintf(fd,"#endif /* WORDS_BIGENDIAN */\n\n\n"):

fclose(fd):

print("************DONE************");



#################################################

if(1+1=3) then
# Modifications to hi_lo
# x ~ x_hi + x_lo
hi_lo_sincos:= proc(x)
  local x_hi, x_lo, res, s,m,e, den, num:
  if x=nearest(x) then
    x_hi := x; x_lo:= 0;
  else
    s,e,m := ieeedouble(x):

    num:=numer(m);
    den:=denom(m);
    num:=round(num/2^25) * 2^25;
    x_hi:=nearest(s*(num/den)*2^e);
    res:=x-x_hi:
    if (res = 0) then
      x_lo:=0:
    else
      x_lo:=nearest(evalf(res)):
    end if;
  end if;
x_hi,x_lo;
end:


fi:



#####################################################################
#Old stuff for SCS, cut from various old maple worksheets, not checked

if(1+2=4) then
Poly_P    := series(sin(sqrt(x))/(x^(3/2))-1/x, x=0, 40):
Poly_Q    := convert(Poly_P,polynom):
Poly_cheb := chebpade(Poly_Q, x=0..evalf(Pi/4), [17,0]);
Poly_Res  := x + x^3 * subs(x=x^2, Poly_cheb);
log(infnorm( 1 - (Poly_Res)/sin(x),x=0..evalf(Pi/4), err))/log(2.);

fi:
