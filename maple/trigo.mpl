
Digits := 100:
interface(quiet=true):
read "common-procedures.mpl":


mkdir("TEMPTRIG"):


#################################
#Comments and todos :

#don't forget the error on k in Cody and Waite


# - Evaluation scheme :
# case 1 : return x
# case 2 (or Fast) : compute a simple polynomial
# case 3 : do an argument reduction...



########################################################
#  Case 1 : Small arguments
# return x for sine and tan, return 1 for cos
########################################################

xmax_return_x_for_sin := 2^(-26):
xmax_return_1_for_cos := 2^(-27):
xmax_return_x_for_tan := 2^(-26):





########################################################
# Case 2 : simple polynomial approximation
########################################################



########################## Fast sine ###########################

# These are the parameters to vary for case 2 
xmaxSinCase2   := Pi/256; #  Should always be larger than Pi/512
degreeSinCase2 := 8;




# Compute the Taylor series
polySinCase2:=  poly_exact2 (convert( series(sin(x), x=0, degreeSinCase2), polynom),2);

# An attempt to use Chebychev polynomials, less accurate in our case after poly_exact2
# If we need extra precision for the polynomial evaluation use Chebpade as follows :
#with(orthopoly):
#Poly_P := series(sin(sqrt(x))/(x^(3/2))-1/x, x=0, degreeSinCase2*4):
#Poly_Q := convert(Poly_P,polynom):
#Poly_cheb := numapprox[chebpade](Poly_Q, x=0..xmaxSinCase2^2, [degreeSinCase2/2-2,0]):
#polySinCase2 := poly_exact2 ( expand(x + x^3 * subs(x=x^2, Poly_cheb)), 2);
#

approx_errorSinCase2:=numapprox[infnorm](1 - polySinCase2 / sin(x), x=-xmaxSinCase2..xmaxSinCase2):
log2(approx_errorSinCase2);

# remove the first x and compute the polynomial of x**2
polySinCase22 :=  subs(x=sqrt(y), expand(polySinCase2/x-1));
x2maxSinCase2:= xmaxSinCase2**2;

# evaluate this polynomial in double. The error on x*x is at most half an ulp
errlist:=errlist_quickphase_horner(degree(polySinCase22),0,0,2**(-53), 0):
rounding_error1:=compute_horner_rounding_error(polySinCase22,y,x2maxSinCase2, errlist, true):
maxeps1 := rounding_error1[1]:
log2(maxeps1);
# get the max value of this polynomial
maxP1:= rounding_error1[4];

# now we multiply this result by x
# Here x is exact, so $x\otimes P(y) = x P(y)(1+\epsilon_{-53})$
# therefore the relative error respective to the stored polynomial is

maxeps2 := (1+maxeps1)*(1+2**(-53))-1:
log2(maxeps2);
# the maximum value of maxP is scaled down, of course
maxP2 := xmaxSinCase2*maxP1;

# and use a fast2Sum to add the last x. Again x is exact

maxeps3 := (maxP1*maxeps2 + 2**(-100)*(1+maxP1)*(1+maxeps2) )   / (1-maxP1) ;
log2(maxeps3);


maxepstotalSinCase2 := (1+maxeps3) * (1+approx_errorSinCase2) - 1 ;
rnconstantSinCase2 := evalf(compute_rn_constant(maxepstotalSinCase2));





##################################### Case2 cos ###########################
# These are the parameters to vary
xmaxCosCase2   := Pi/256;
degreeCosCase2 := 7;


# Compute the Taylor series
polyCosCase2  := poly_exact2 (convert( series(cos(x), x=0, degreeCosCase2), polynom),2);
delta_approx := numapprox[infnorm](polyCosCase2 - cos(x), x=-xmaxCosCase2..xmaxCosCase2):
log2(%);

# remove the first 1 and compute the polynomial of x**2
polyCosCase22 := subs(x=sqrt(y), expand(polyCosCase2-1));
x2maxCosCase2 := xmaxCosCase2**2;

# evaluate this polynomial in double. The error on x*x is at most half an ulp
errlist        := errlist_quickphase_horner(degree(polyCosCase22),0,0,2**(-53), 0):
rounding_error1:= compute_horner_rounding_error(polyCosCase22,y,x2maxCosCase2, errlist, true):
delta_round    := rounding_error1[2]:
log2(%);

# Then we have an Add12 which is exact. The result is greater then cos(xmaxCosCase2):
miny := cos(xmaxCosCase2);
maxepstotalCosCase2 :=  (delta_round + delta_approx) / miny ;
log2(%);
rnconstantCosCase22 := evalf(compute_rn_constant(maxepstotalCosCase2));
rnconstantCosCase23 := evalf(compute_rn_constant(2**(-67))):




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








#################################################
#   Case 3 : Argument reduction
#################################################



# TODO  This value seems to work but needs proving

maxepstotalSinCase3:=2**(-65); 
rnconstantSinCase3 := evalf(compute_rn_constant(maxepstotalSinCase3)):


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







################################## Polynomials for do_sine and do_cos

ymax  := Pi/512;
y2max := ymax**2;

# The error on yh*yh
# we had $y_h+y_l = y + \abserr{CodyWaite}$.
# Now we take only $y_h$ :  $y_h = (y_h + y_l)(1+2^{-53}) = (y+\abserr{CodyWaite})(1+2^{-53})$
# When squared we get $y_h\otimes y_h = (y+\abserr{CodyWaite})^2(1+2^{-53})^3 $

epsy2 := evalf(((1+eps_ArgRed)**2) *  (1+2**(-53))**2 - 1):

############### Computing Ts

polySin:=  poly_exact (convert( series(sin(x), x=0, 9), polynom),x):
polyTs :=  expand(polySin/x-1):
delta_approx_Ts := numapprox[infnorm](polyTs - (sin(x)/x-1), x=-ymax..ymax):
maxTs   := numapprox[infnorm](polyTs, x=-ymax..ymax):
polyTs2 := subs(x=sqrt(y),polyTs):

errlist:=errlist_quickphase_horner(degree(polyTs2),0,0, epsy2 , 0):
(eps_rounding_Ts, delta_rounding_Ts, minTs, maxTs):=
	compute_horner_rounding_error(polyTs2,y,y2max, errlist, true):


############### Computing Tc

polyCos:= polyCosCase2 -x^8*coeff(polyCosCase2,x,8):
# More accurate
#poly_exact(subs(y=x^2, numapprox[minimax]((cos(sqrt(y))), y=2^(-2048)..y2max, [3,0])), x);

polyTc:=polyCos - 1:
delta_approx_Tc:= numapprox[infnorm](polyCos -  cos(x), x=-ymax..ymax):
maxTc   := numapprox[infnorm](polyTc, x=-ymax..ymax):
polyTc2 := subs(x=sqrt(y),polyTc):

errlist:=errlist_quickphase_horner(degree(polyTc2),0,0, epsy2 , 0):
(eps_rounding_Tc, delta_rounding_Tc, minTc, maxTc):=compute_horner_rounding_error(polyTc2,y,y2max, errlist, true):




###### The extreme cases for tabulated values for the case k&127 != 0
minsca := sin(Pi/256);
minscah:= nearest(minsca);
maxsca := cos(Pi/256);
maxscah:= nearest(maxsca);

# Worst case of approximation error

################ Summing everything up
# The only error is in
#  tlo =  tc*cah - (ts*sahyh_h -  (cal + (tlo  - (sahyh_l + (sal*yh + sah*yl)) )));

delta_round_tlo := 2*ulp(maxTc*maxsca) ; # TODO c'est à la louche


delta_do_cos := delta_round_tlo;
min_sin := sin(Pi/512); # TODO Add error on k here

delta_sincos := delta_round_tlo / min_sin; # TODO hum

rnconstant_sincos := compute_rn_constant(delta_sincos);



# an auxiliary output function:
# Outputs the high part of a double, and the double in comment
outputHighPart:=proc(cvarname, var)
  Digits:=8:
  ("#define " || cvarname || " 0x" || (ieeehexa(var)[1]) 
    ||  "    /* " || (convert(evalf(var),string)) ||  " */" )
end proc;





# Output:

filename:="TEMPTRIG/trigo_fast.h":
fd:=fopen(filename, WRITE, TEXT):

fprintf(fd, "#include \"crlibm.h\"\n#include \"crlibm_private.h\"\n"):

fprintf(fd, "\n/*File generated by maple/coef_sine.mw*/\n"):


fprintf(fd,  "%s\n",  outputHighPart("XMAX_RETURN_X_FOR_SIN", xmax_return_x_for_sin) ):
fprintf(fd,  "%s\n",  outputHighPart("XMAX_SIN_FAST", xmaxSinCase2) ):

fprintf(fd,  "%s\n",  outputHighPart("XMAX_RETURN_1_FOR_COS", xmax_return_1_for_cos) ):
fprintf(fd,  "%s\n",  outputHighPart("XMAX_COS_FAST", xmaxCosCase2) ):

fprintf(fd,  "%s\n",  outputHighPart("XMAX_RETURN_X_FOR_TAN", xmax_return_x_for_tan) ):
fprintf(fd,  "%s\n",  outputHighPart("XMAX_TAN_FAST", xmaxTanCase2) ):


# TODO remove the following
#fprintf(fd, "#define XMAX_RETURN_1_FOR_COS 0x%s\n", ieeehexa(xmax_return_1_for_cos)[1]):
#fprintf(fd, "#define XMAX_COS_FAST  0x%s\n", ieeehexa(xmaxCosCase2)[1]):
#fprintf(fd, "#define XMAX_COS_FAST2 0x%s\n", ieeehexa(xmaxCosCase22)[1]):

#fprintf(fd, "#define XMAX_RETURN_X_FOR_TAN 0x%s\n", ieeehexa(xmax_return_x_for_tan)[1]):
#fprintf(fd, "#define XMAX_TAN_FAST  0x%s\n", ieeehexa(xmaxTanCase2)[1]):


fprintf(fd, "\n"):

fprintf(fd, "#define EPS_SIN_CASE2 %e \n",    maxepstotalSinCase2);
fprintf(fd, "#define RN_CST_SIN_CASE2 %f \n", rnconstantSinCase2);
fprintf(fd, "#define EPS_SIN_CASE3 %e \n",    maxepstotalSinCase3);
fprintf(fd, "#define RN_CST_SIN_CASE3 %f \n", rnconstantSinCase3);

fprintf(fd, "#define EPS_COS_CASE2 %e \n",    maxepstotalCosCase2);
fprintf(fd, "#define RN_CST_COS_CASE2 %f \n", rnconstantCosCase2);
fprintf(fd, "#define EPS_COS_CASE3 %e \n",    maxepstotalCosCase3);
fprintf(fd, "#define RN_CST_COS_CASE3 %f \n", rnconstantCosCase3);

fprintf(fd, "#define RN_CST_TANFAST1 %e \n", rnconstantTanCase21);
fprintf(fd, "#define RN_CST_TANFAST2 %e \n", rnconstantTanCase22);
fprintf(fd, "\n"):

fprintf(fd, "#define INV_PIO256 %1.50f \n", 1/C);
fprintf(fd, "\n"):

fprintf(fd, "#define XMAX_CODY_WAITE_2 0x%s\n", ieeehexa(XMAX_CODY_WAITE_2)[1]):
fprintf(fd, "#define XMAX_CODY_WAITE_3 0x%s\n", ieeehexa(XMAX_CODY_WAITE_3)[1]):
fprintf(fd, "#define XMAX_DDRR 0x%s\n", ieeehexa(XMAX_DDRR)[1]):
fprintf(fd, "\n"):

fprintf(fd, "#define RR_CW2_CH %1.50e\n", RR_CW2_CH):
fprintf(fd, "#define RR_CW2_MCL %1.50e\n", RR_CW2_MCL):
fprintf(fd, "\n"):

fprintf(fd, "#define RR_CW3_CH %1.50e\n", RR_CW3_CH):
fprintf(fd, "#define RR_CW3_CM %1.50e\n", RR_CW3_CM):
fprintf(fd, "#define RR_CW3_MCL %1.50e\n", RR_CW3_MCL):
fprintf(fd, "\n"):

fprintf(fd, "#define RR_DD_MCH %1.50e\n", RR_DD_MCH):
fprintf(fd, "#define RR_DD_MCM %1.50e\n", RR_DD_MCM):
fprintf(fd, "#define RR_DD_CL %1.50e\n", RR_DD_CL):
fprintf(fd, "\n"):

fprintf(fd, "#define RN_CST_SINCOS %f \n", rnconstant_sincos);
fprintf(fd, "\n"):



fprintf(fd,"#ifdef WORDS_BIGENDIAN\n"):
for isbig from 1 to 0 by -1 do

  if(isbig=0) then
    fprintf(fd,"#else\n"):
  fi;

  # The sine polynomial

  fprintf(fd, "static db_number const s3 = "):
  printendian(fd, coeff(polySinCase2,x,3), isbig):
  fprintf(fd, ";\n"):
  fprintf(fd, "static db_number const s5 = "):
  printendian(fd, coeff(polySinCase2,x,5), isbig):
  fprintf(fd, ";\n"):
  fprintf(fd, "static db_number const s7 = "):
  printendian(fd, coeff(polySinCase2,x,7), isbig):
  fprintf(fd, ";\n\n"):


  # the cos polynomial

  fprintf(fd, "static db_number const c2 = "):
  printendian(fd, coeff(polyCosCase2,x,2), isbig):
  fprintf(fd, ";\n"):
  fprintf(fd, "static db_number const c4 = "):
  printendian(fd, coeff(polyCosCase2,x,4), isbig):
  fprintf(fd, ";\n"):
  fprintf(fd, "static db_number const c6 = "):
  printendian(fd, coeff(polyCosCase2,x,6), isbig):
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


  fprintf(fd, "\n\n"):

  #The sincos table

  fprintf(fd, "\n/*  sine and cos of kPi/256 in double-double */\n"):
  SinCosSize:= 128;
  fprintf(fd, "static db_number const sincosTable[%d] =\n{\n",  4*(SinCosSize/2+1)):
  for i from 0 to SinCosSize/2 do
      s:=hi_lo(sin(i*Pi/(2*SinCosSize)));
      c:=hi_lo(cos(i*Pi/(2*SinCosSize)));
      printendian(fd,s[1],isbig);
      fprintf(fd," ,\n"):
      printendian(fd,s[2],isbig);
      fprintf(fd," ,\n"):
      printendian(fd,c[1],isbig);
      fprintf(fd," ,\n"):
      printendian(fd,c[2],isbig);

    if i<SinCosSize-1 then fprintf(fd," ,\n"): fi:
  od:
  fprintf(fd, "\n};\n\n"):


od:
fprintf(fd,"#endif /* WORDS_BIGENDIAN */\n\n\n"):

fclose(fd):



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
