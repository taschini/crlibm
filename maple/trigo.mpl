
Digits := 100:
 interface(quiet=true):
read "Common_Maple_Procedures.mpl":
# read "temp_common.mpl";


mkdir("TEMPTRIG"):


#################################
#Comments and todos : 

#don't forget the error on k in Cody and Waite




#########################################################
# Small arguments

# When do we return x ?
xmax_return_x_for_sin := 2^(-26):
xmax_return_1_for_cos := 2^(-27):


#####################################Fast sine###########################
# These are the parameters to vary
xmaxSinFast:=2**(-4);
degreeSinFast:=10;

# Compute the Taylor series
polySinFast:=  poly_exact2 (convert( series(sin(x), x=0, degreeSinFast), polynom),2);
approx_errorSinFast:=numapprox[infnorm](1 - polySinFast / sin(x), x=-xmaxSinFast..xmaxSinFast):
log2(approx_errorSinFast);

# remove the first x and compute the polynomial of x**2
polySinFast2 :=  subs(x=sqrt(y), expand(polySinFast/x-1));
x2maxSinFast:= xmaxSinFast**2;

# evaluate this polynomial in double. The error on x*x is at most half an ulp
errlist:=errlist_quickphase_horner(degree(polySinFast2),0,0,2**(-53), 0):
rounding_error1:=compute_horner_rounding_error(polySinFast2,y,x2maxSinFast, errlist, true):
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
maxP2 := xmaxSinFast*maxP1;

# and use a fast2Sum to add the last x. Again x is exact

maxeps3 := (maxP1*maxeps2 + 2**(-100)*(1+maxP1)*(1+maxeps2) )   / (1-maxP1) ;
log2(maxeps3);


maxepstotalSinFast := (1+maxeps3) * (1+approx_errorSinFast) - 1 ;
rnconstantSinFast := evalf(compute_rn_constant(maxepstotalSinFast));





##################################### Fast cos ###########################
# These are the parameters to vary
xmaxCosFast:=2**(-5);
degreeCosFast:=10;

# Compute the Taylor series
polyCosFast:=  poly_exact2 (convert( series(cos(x), x=0, degreeCosFast), polynom),2);
delta_approx := numapprox[infnorm](polyCosFast - cos(x), x=-xmaxCosFast..xmaxCosFast):
log2(%); 

# remove the first 1 and compute the polynomial of x**2
polyCosFast2 :=  subs(x=sqrt(y), expand(polyCosFast-1));
x2maxCosFast:= xmaxCosFast**2;

# evaluate this polynomial in double. The error on x*x is at most half an ulp
errlist:=errlist_quickphase_horner(degree(polyCosFast2),0,0,2**(-53), 0):
rounding_error1:=compute_horner_rounding_error(polyCosFast2,y,x2maxCosFast, errlist, true):
delta_round := rounding_error1[2]:
log2(%); 

# Then we have an Add12 which is exact. The result is greater then cos(xmaxCosFast):
miny := cos(xmaxCosFast);
maxepstotalCosFast :=  (delta_round + delta_approx) / miny ;
log2(%);
rnconstantCosFast := evalf(compute_rn_constant(maxepstotalCosFast));















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




# Fastest reduction using two-part Cody and Waite

bitsCh_0:=32:  # ensures at least 53+11 bits

# 1/2 <= C/2^(expC+1) <1 
Ch:= round(evalf(  C * 2^(bitsCh_0-expC-1))) / (2^(bitsCh_0-expC-1)):
# recompute bitsCh in case we are lucky (and we are for bitsCh_0=32)
bitsCh:=1+log2(op(2,ieeedouble(Ch)[3])) :  # this means the log of the denominator

Cl:=nearest(C - Ch):
# Cody and Waite argument reduction will work for k<kmax
kmax:=2^(53-bitsCh):

# The constants to move to the .h file
CW2_CH := Ch;
CW2_MCL := -Cl:
XMAX_CODY_WAITE_2 := nearest(kmax*C):

# The error in this case (we need absolute error) 
delta_repr_C := abs(C-Ch-Cl);
delta_round := 1/2 * ulp(Cl) ;
delta_cody_waite_2 := kmax * (delta_repr_C + delta_round);

log2(%);




# Slower reduction using three-part Cody and Waite

bitsCh_0:=21:
Ch:= round(evalf(  C * 2^(bitsCh_0-expC-1))) / (2^(bitsCh_0-expC-1)):
# recompute bitsCh in case we are lucky (and we are for bitsCh_0=32)
bitsCh:=1+log2(op(2,ieeedouble(Ch)[3])) :  # this means the log of the denominator

r := C-Ch:
Cmed := round(evalf(  r * 2^(2*bitsCh-expC-1))) / (2^(2*bitsCh-expC-1)):
bitsCmed:=1+log2(op(2,ieeedouble(Cmed)[3])) : 

  (53-max(bitsCh, bitsCmed + 1));

Cl:=nearest(C - Ch - Cmed):
# Cody and Waite argument reduction will work for k<kmax
kmax:=2^(53-max(bitsCh, bitsCmed + 1)):

kmax:=2^31: # Otherwise we have integer overflow

cw2_kmax := kmax;


# The constants to move to the .h file
CW3_CH := Ch;
CW3_CM := Cmed:
CW3_MCL := -Cl:
XMAX_CODY_WAITE_3 := nearest(kmax*C):

# The error in this case (we need absolute error) 
delta_repr_C := abs(C - Ch - Cmed - Cl):
delta_round := 1/2 * ulp(Cl) :
delta_cody_waite_3 := kmax * (delta_repr_C + delta_round):

log2(%);


XMAX_DDRR:=nearest((2^51-1)*C);



delta_ArgRed := max(delta_cody_waite_2, delta_cody_waite_3); # TODO ajouter la troisième




# Now we use the above range reduction when k mod 256 <> 0
# otherwise we need to worry about relative accuracy of the result.
# what is the worst case for cancellation ?

emax := ieeedouble(XMAX_CODY_WAITE_3)[2] +1 ;
# above emax, we will use Payne and Hanek
(wcn, wce, wceps) := WorstCaseForAdditiveRangeReduction(2,53,-8, emax, C):
wcx := wcn * 2^wce:
wck := round(wcx/C): 
wcy := wcx - wck*C:

log2(wcy);   # y < 2^(-67);

# What is the precision we can expect if we do a stupid range reduction storing C as 3 doubles ? 
Ch := nearest(C);
Cmed := nearest(C-Ch);
Cl := nearest(C-Ch-Cmed);

delta_repr_C := abs(C - Ch - Cmed - Cl):
delta_round := 1/2 * ulp(Cl) :
delta_max :=  (delta_repr_C + delta_round) * cw2_kmax;
eps_ArgRed_k0 := delta_max/wcy :
eps_ArgRed_k0 := (1+eps_ArgRed_k0)*(1+2^(-102))-1:
log2(eps_ArgRed_k0);

RRK0_MCH := -Ch;
RRK0_MCM := -Cmed;
RRK0_CL := Cl;



eps_ArgRed := eps_ArgRed_k0; # TODO Correct this



###########################################################################################
#            Payne and Hanek argument reduction
#

# We use the SCS version...



################################## Polynomials for do_sine and do_cos 

ymax:=Pi/512;
y2max:= ymax**2; 

# The error on yh*yh 
# we had $y_h+y_l = y + \abserr{CodyWaite}$. 
# Now we take only $y_h$ :  $y_h = (y_h + y_l)(1+2^{-53}) = (y+\abserr{CodyWaite})(1+2^{-53})$
# When squared we get $y_h\otimes y_h = (y+\abserr{CodyWaite})^2(1+2^{-53})^3 $

epsy2 := evalf(((1+eps_ArgRed)**2) *  (1+2**(-53))**2 - 1):

############### Computing Ts

polySin:=  poly_exact (convert( series(sin(x), x=0, 9), polynom),x):
polyTs :=  expand(polySin/x-1):
delta_approx_Ts := numapprox[infnorm](polyTs - (sin(x)/x-1), x=-ymax..ymax):
maxTs := numapprox[infnorm](polyTs, x=-ymax..ymax):
polyTs2 :=  subs(x=sqrt(y),polyTs):

errlist:=errlist_quickphase_horner(degree(polyTs2),0,0, epsy2 , 0):
(eps_rounding_Ts, delta_rounding_Ts, minTs, maxTs):=
	compute_horner_rounding_error(polyTs2,y,y2max, errlist, true):


############### Computing Tc

if(1+1=3) then # We get a better polynomial below ;
 polyCos:= poly_exact(convert( series(cos(x), x=0, 8), polynom),x);
 polyTc:=expand(polyCos-1);
 deltaTc := numapprox[infnorm](polyTc -  (cos(x)-1), x=-ymax..ymax):
 maxTc := numapprox[infnorm](polyTc, x=-ymax..ymax):
 log2(deltaTc);
 log2(maxTc);
end if;

polyCos:= poly_exact(subs(y=x^2, numapprox[minimax]((cos(sqrt(y))), y=2^(-2048)..y2max, [3,0])), x);
polyTc:=polyCos - 1;
delta_approx_Tc:= numapprox[infnorm](polyCos -  cos(x), x=-ymax..ymax):
maxTc := numapprox[infnorm](polyTc, x=-ymax..ymax):
polyTc2 :=  subs(x=sqrt(y),polyTc):

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


# Output

filename:="TEMPTRIG/trigo.h":
fd:=fopen(filename, WRITE, TEXT):

fprintf(fd, "#include \"crlibm.h\"\n#include \"crlibm_private.h\"\n"):

fprintf(fd, "\n/*File generated by maple/coef_sine.mw*/\n"):


fprintf(fd, "#define XMAX_RETURN_X_FOR_SIN 0x%s\n", ieeehexa(xmax_return_x_for_sin)[1]):
fprintf(fd, "#define XMAX_SIN_FAST 0x%s\n", ieeehexa(xmaxSinFast)[1]):
fprintf(fd, "#define XMAX_RETURN_1_FOR_COS 0x%s\n", ieeehexa(xmax_return_1_for_cos)[1]):
fprintf(fd, "#define XMAX_COS_FAST 0x%s\n", ieeehexa(xmaxSinFast)[1]):
fprintf(fd, "\n"):

fprintf(fd, "#define RN_CST_SINFAST %f \n", rnconstantSinFast);
fprintf(fd, "#define RN_CST_COSFAST %f \n", rnconstantCosFast);
fprintf(fd, "\n"):

fprintf(fd, "#define INV_PIO256 %1.50f \n", 1/C);
fprintf(fd, "\n"):

fprintf(fd, "#define XMAX_CODY_WAITE_2 0x%s\n", ieeehexa(XMAX_CODY_WAITE_2)[1]):
fprintf(fd, "#define XMAX_CODY_WAITE_3 0x%s\n", ieeehexa(XMAX_CODY_WAITE_3)[1]):
fprintf(fd, "#define XMAX_DDRR 0x%s\n", ieeehexa(XMAX_DDRR)[1]):
fprintf(fd, "\n"):

fprintf(fd, "#define CW2_CH %1.50e\n", CW2_CH):
fprintf(fd, "#define CW2_MCL %1.50e\n", CW2_MCL):
fprintf(fd, "\n"):

fprintf(fd, "#define CW3_CH %1.50e\n", CW3_CH):
fprintf(fd, "#define CW3_CM %1.50e\n", CW3_CM):
fprintf(fd, "#define CW3_MCL %1.50e\n", CW3_MCL):
fprintf(fd, "\n"):

fprintf(fd, "#define RRK0_MCH %1.50e\n", RRK0_MCH):
fprintf(fd, "#define RRK0_MCM %1.50e\n", RRK0_MCM):
fprintf(fd, "#define RRK0_CL %1.50e\n", RRK0_CL):
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
  printendian(fd, coeff(polySin,x,3), isbig):
  fprintf(fd, ";\n"):
  fprintf(fd, "static db_number const s5 = "):
  printendian(fd, coeff(polySin,x,5), isbig):
  fprintf(fd, ";\n"):
  fprintf(fd, "static db_number const s7 = "):
  printendian(fd, coeff(polySinFast,x,7), isbig):
  fprintf(fd, ";\n"):
  fprintf(fd, "static db_number const s9 = "):
  printendian(fd, coeff(polySinFast,x,9), isbig):
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
  fprintf(fd, ";\n"):
  fprintf(fd, "static db_number const c8 = "):
  printendian(fd, coeff(polyCosFast,x,8), isbig):
  fprintf(fd, ";\n"):

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
Poly_P := series(sin(sqrt(x))/(x^(3/2))-1/x, x=0, 40):
Poly_Q := convert(Poly_P,polynom):
Poly_cheb := chebpade(Poly_Q, x=0..evalf(Pi/4), [17,0]);
Poly_Res := x + x^3 * subs(x=x^2, Poly_cheb);
log(infnorm( 1 - (Poly_Res)/sin(x),x=0..evalf(Pi/4), err))/log(2.);

fi: