
Digits := 100:
interface(quiet=true):
read "common-procedures.mpl":


mkdir("TEMPTRIG"):



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

xmax_return_x_for_tan := 2^(-27):




########################################################
# Case 2 : simple polynomial approximation
########################################################

# We want to use the same polynomial in case 2 and 3.
# So see after arg red

#################################################
#   Case 3 : Argument reduction
#################################################



#################################################
# CODY and WAITE  Argument reduction


C := Pi/256:
invC:= nearest(1/C):
reminvC := evalf(1/C - invC):
expC:=ieeedouble(C)[2]:
epsinvC := abs(reminvC*C):

# There are three sets of constants :
#  - split redC into two constants, for small values when we are concerned with absolute error
#  - split redC into three constants, for larger values
#  - split redC into three doubles, for the cases when we need
#     good relative precision on the result and fear cancellation




# Fastest reduction using two-part Cody and Waite (up to |k|=2^22)

bitsCh_0:=34:

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
delta_repr_C_cw2   := abs(C-Ch-Cl):
delta_round_cw2    := kmax_cw2* 1/2 * ulp(Cl) :
delta_cody_waite_2 := kmax_cw2 * delta_repr_C_cw2 + delta_round_cw2:
# This is the delta on y, the reduced argument

#log2(%);



# Slower reduction using three-part Cody and Waite, up to |k|=2^31

bitsCh_0:=23: # 22 or 23
Ch:= round(evalf(  C * 2^(bitsCh_0-expC-1))) / (2^(bitsCh_0-expC-1)):
# recompute bitsCh in case we are lucky
bitsCh:=1+log2(op(2,ieeedouble(Ch)[3])) :  # this means the log of the denominator

r := C-Ch:
Cmed := round(evalf(  r * 2^(2*bitsCh-expC-1))) / (2^(2*bitsCh-expC-1)):
bitsCmed:=1+log2(op(2,ieeedouble(Cmed)[3])) :

Cl:=nearest(C - Ch - Cmed):

kmax_cw3 := 2^min(53-bitsCh, 53-bitsCmed, 31):# Otherwise we have integer overflow



# The constants to move to the .h file
RR_CW3_CH  := Ch:
RR_CW3_CM  := Cmed:
RR_CW3_MCL := -Cl:
XMAX_CODY_WAITE_3 := nearest(kmax_cw3*C):

# The error in this case (we need absolute error)
delta_repr_C_cw3   := abs(C - Ch - Cmed - Cl):
delta_round_cw3    := kmax_cw3 * 1/2 * ulp(Cl) :
delta_cody_waite_3 := kmax_cw3 * delta_repr_C_cw3 + delta_round_cw3:
# This is the delta on y, the reduced argument

#log2(%);





# Third range reduction, using double-double arithmetic, for |k| up to 2^51-1

# This max int value can be produced by DOUBLE2LONGINT
kmax:=2^46-1:
XMAX_DDRR:=nearest(kmax*C):

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

print("delta_ArgRed to move to the .gappa file = ",  evalf(delta_ArgRed)):
#log2(delta_ArgRed);





# Now we use the above absolute error when k mod 256 <> 0
# otherwise we need to worry about relative accuracy of the result.

# First, what is the worst case for cancellation ?

emax := ieeedouble(XMAX_DDRR)[2] +1 :
# above emax, we will use Payne and Hanek so we do not worry

(wcn, wce, wceps) := WorstCaseForAdditiveRangeReduction(2,53,-8, emax, C):
wcx := wcn * 2^wce:
wck := round(wcx/C):
wcy := wcx - wck*C:

#log2(wcy);   # y > 2^(-67);

# In these cases we use the double-double range reduction, for |k|<kmax_cw3
# and the relative precision in the worst case is for wcy

delta_round := kmax_cw3 * 1/2 * ulp(Cl)      # for   kd*RR_DD_CL
              + kmax_cw3 * ulp(Cl) :         # for the subtraction

delta_RR_DD :=  kmax_cw3 * delta_repr_C + delta_round:

eps_ArgRed := (1+delta_RR_DD/wcy)*(1+2^(-100)) -1:

#log2(eps_ArgRed);

# In all the other cases we use Payne and Hanek, and eps_ArgRed is
# much smaller, so this is the max.



###########
# Polynomials for do_sine and do_cos, and for the case 2
degreeSin := 8:
degreeCos := 7:

maxepsk := (1+epsinvC)*(1+2^(-53))-1:

ymaxCase3  := evalf(Pi/512 + XMAX_DDRR*maxepsk):
print("ymaxCase3 to move to the .gappa file = ",  ymaxCase3):


y2maxCase3 := ymaxCase3^2:
# These are the parameters that can be varied to fine-tune performance
#   (they should always be larger than Pi/512
xmaxCosCase2   := Pi/256:
xmaxSinCase2   := Pi/256:


x2maxSinCase2:= xmaxSinCase2^2:
x2maxCosCase2:= xmaxCosCase2^2:

# We had the difficulty here to find minimax polynomials which are good for case 2
# as well as for case 3. A simple solution was to set xmaxCosCase2 =ymaxCase3...
# However we found another answer: in the future we intend to use these polynomials for second
# step, too. Therefore, no minimax, only Taylor.

polySin:=  poly_exact(convert( series(sin(x), x=0, degreeSin+1), polynom)):
polyTs := expand(polySin/x-1):
polyTs2 := subs(x=sqrt(y), polyTs):

polyCos  := poly_exact (convert( series(cos(x), x=0, degreeCos+1), polynom)):
polyTc2 := subs(x=sqrt(y), polyCos - 1):

eps_approx_Sin_Case2 := numapprox[infnorm]((x*polyTs+x -sin(x))/sin(x), x=0..xmaxSinCase2):
eps_approx_Sin_Case3 := numapprox[infnorm]((x*polyTs +x -sin(x))/sin(x), x=0..ymaxCase3):

delta_approx_Sin_Case2 := numapprox[infnorm]((x*polyTs+x -sin(x)), x=0..xmaxSinCase2):
delta_approx_Sin_Case3 := numapprox[infnorm]((x*polyTs +x -sin(x)), x=0..ymaxCase3):

delta_approx_Cos_Case2:= numapprox[infnorm](polyCos -  cos(x), x=0..xmaxCosCase2):
delta_approx_Cos_Case3:= numapprox[infnorm](polyCos -  cos(x), x=0..ymaxCase3):

print("delta_approx_Sin_Case3 to move to the .gappa file = ",  delta_approx_Sin_Case3):
print("delta_approx_Cos_Case3 to move to the .gappa file = ",  delta_approx_Cos_Case3):


########################## Case 2 for sine  ###########################



# evaluate this polynomial in double. The error on x*x is at most half an ulp
errlist:=errlist_quickphase_horner(degree(polyTs2),0,0,2^(-53), 0):
(eps_rounding_Ts, delta_rounding_Ts, minTs, maxTs):=
	compute_horner_rounding_error(polyTs2,y,x2maxSinCase2, errlist, true):

eps_poly_Ts_Case2 := numapprox[infnorm]((x*polyTs)/(sin(x)-x) -1, x=0..xmaxSinCase2):
maxeps2 := (1+eps_poly_Ts_Case2)*(1+eps_rounding_Ts)*(1+2^(-53))-1:

maxepstotalSinCase2 := maxeps2 * numapprox[infnorm](1-x/sin(x), x=0..xmaxSinCase2):
rnconstantSinCase2 := evalf(compute_rn_constant(maxepstotalSinCase2)):





##################################### Case2 cos ###########################

# evaluate this polynomial in double. The error on x*x is at most half an ulp
errlist        := errlist_quickphase_horner(degree(polyTc2),0,0,2**(-53), 0):
(eps_rounding_Tc, delta_rounding_Tc, minTc, maxTc):=
              compute_horner_rounding_error(polyTc2,y,x2maxCosCase2, errlist, true):

# Then we have an Add12 which is exact. The result is greater then cos(xmaxCosCase2):
miny := cos(xmaxCosCase2):
maxepstotalCosCase2 :=  (delta_rounding_Tc + delta_approx_Cos_Case2) / miny :
#log2(%);
rnconstantCosCase2 := evalf(compute_rn_constant(maxepstotalCosCase2)):




######################## Case2 Tangent #########################
#
# Compute the Taylor series
with(orthopoly):
degreeTanCase2 := 12:
xmaxTanCase2   := 2**(-4):
xminTanCase2   := 2**(-30):

Poly_P := convert(series(tan(sqrt(x))/(x^(3/2))-1/x, x=0, degreeTanCase2*4),polynom):
Poly_cheb := numapprox[chebpade](Poly_P, x=xminTanCase2..xmaxTanCase2^2, [degreeTanCase2/2-2,0]):
polyTanCase2 :=  poly_exact2(expand(x + x^3 * subs(x=x^2, Poly_cheb)), 4):

maxepsApproxTanCase2:=numapprox[infnorm](1 - polyTanCase2 / tan(x), x=xminTanCase2..xmaxTanCase2):

# Now we pass these values to Gappa

filename:="TEMPTRIG/TanCase2.sed":
fd:=fopen(filename, WRITE, TEXT):
  t3h, t3l := hi_lo(coeff(polyTanCase2,x,3)):
  fprintf(fd, " s/_t3h/%1.40e/g\n", t3h):
  fprintf(fd, " s/_t3l/%1.40e/g\n", t3l):
  fprintf(fd, " s/_t5/%1.40e/g\n",  coeff(polyTanCase2,x,5)):
  fprintf(fd, " s/_t7/%1.40e/g\n",  coeff(polyTanCase2,x,7)):
  fprintf(fd, " s/_t9/%1.40e/g\n",  coeff(polyTanCase2,x,9)):
  fprintf(fd, " s/_t11/%1.40e/g\n", coeff(polyTanCase2,x,11)):
  fprintf(fd, " s/_XMAX/%1.40e/g\n", xmaxTanCase2):
  fprintf(fd, " s/_MAXEPSAPPROX/%1.40e/g\n", maxepsApproxTanCase2*1.00001):
fclose(fd);

printf("Now you should use  \n    sed -f TEMPTRIG/TanCase2.sed trigoTanCase2.gappa | gappa  > /dev/null \n");



maxepstotalTanCase2:=4.59602e-19:  # Cut from Gappa output

log2(maxepstotalTanCase2): # almost 61 bits


#rnconstantTanCase22 := evalf(compute_rn_constant(maxepstotalTanCase22)):





###############################################################################
#   Computing errors for Case3 : now we have an error due to arg red


#  This value has been computed by Gappa (using all the previous)
maxepstotalSinCosCase3:=2**(-66):
rnconstantSinCosCase3 := evalf(compute_rn_constant(maxepstotalSinCosCase3)):


# The error of sin, the error of cos, then the error of Div22
maxepstotalTanCase3:= 2.1*maxepstotalSinCosCase3:
rnconstantTanCase3 := evalf(compute_rn_constant(maxepstotalTanCase3)):





##############################################
## Compute constants for SCS arg red
oldDigits:=Digits:
Digits:=1000:
# for 2/Pi:
n:=round(2^(30*48)*evalf(2/Pi)):
digitlist:=[]:
for i from 1 to 48 do
  r:=n mod (2^30):
  n:=floor(n/(2^30)):
  hexstring:= convert(convert(r,hex),string):
  digitlist:=[hexstring, op(digitlist)]:
end:
digitlist:

# for 256/Pi:
n:=round(2^(30*47)*evalf(256/Pi)):
digitlist:=[]:
for i from 1 to 48 do
  r:=n mod (2^30):
  n:=floor(n/(2^30)):
  hexstring:= convert(convert(r,hex),string):
  digitlist:=[hexstring, op(digitlist)]:
end:
digitlist:
Digits:=oldDigits:



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
fprintf(fd, "#define ONE_ROUNDED_DOWN  %1.25e \n", one_rounded_down):
fprintf(fd, "\n"):

fprintf(fd, "#define EPS_SIN_CASE2     %1.25e \n", maxepstotalSinCase2):
fprintf(fd, "#define RN_CST_SIN_CASE2  %1.25f \n", rnconstantSinCase2):
fprintf(fd, "\n"):
fprintf(fd, "#define EPS_COS_CASE2     %1.25e \n", maxepstotalCosCase2):
fprintf(fd, "#define RN_CST_COS_CASE2  %1.25f \n", rnconstantCosCase2):
fprintf(fd, "\n"):
fprintf(fd, "#define EPS_SINCOS_CASE3     %1.25e \n", maxepstotalSinCosCase3):
fprintf(fd, "#define RN_CST_SINCOS_CASE3  %1.25f \n", rnconstantSinCosCase3):
fprintf(fd, "\n"):
fprintf(fd, "#define EPS_TAN_CASE2     %1.25e \n", maxepstotalTanCase2):
fprintf(fd, "#define EPS_TAN_CASE3     %1.25e \n", maxepstotalTanCase3):
fprintf(fd, "#define RN_CST_TAN_CASE3  %1.25f \n", rnconstantTanCase3):

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

printf("\n\n************ DONE TEMPTRIG/trigo_fast.h ************\n Copy it to the crlibm source directory.\n\n");


# The Gappa files in TEMPTRIG
for i from 1 to SinCosSize/2 do
    filename:=cat("TEMPTRIG/SinACosA_",i,".gappa"):
    fd:=fopen(filename, WRITE, TEXT):
    s:=hi_lo(sin(i*Pi/(2*SinCosSize))):
    c:=hi_lo(cos(i*Pi/(2*SinCosSize))):
    fprintf(fd, "cah=%1.50e ;\n", c[1]):
    fprintf(fd, "cal=%1.50e ;\n", c[2]):
    fprintf(fd, "sah=%1.50e ;\n", s[1]):
    fprintf(fd, "sal=%1.50e ;\n", s[2]):
    fclose(fd):
od:


printf("************ DONE TEMPTRIG/*.gappa ************\n");

# A shell script to use them
filename:="trigo_test.sh":
fd:=fopen(filename, WRITE, TEXT):
fprintf(fd, "#!/bin/sh\n"):
fprintf(fd, "for file in TEMPTRIG/SinACosA*  \n"):
fprintf(fd, "do\n"):
fprintf(fd, "  echo $file:\n"):
fprintf(fd, "  cat  $file  trigoSinCosCase3.gappa | ~/gappa/src/gappa > /dev/null\n"):
fprintf(fd, "  echo\n"):
fprintf(fd, "done\n"):
fclose(fd):

printf("************ DONE trigo_test.sh ************"):
printf("Now you should run"):
printf(" sh trigo_test.sh 2>TEMPTRIG/Gappa.out"):
printf("Then"):
printf(" grep '{' TEMPTRIG/Gappa.out| sed  -e 's/}.*{/\n/g'  -e 's/^.*{-//g' -e 's/}]//g'  | sort | head "):

printf("If the first number is smaller than 2^(-66) then everything is OK and the rounding constants in TEMPTRIG/trigo_fast.h are proven upper bounds."):



#################################################
# Stuff for the proof of the accurate phase,
#
if(1+1=3) then
 polysin:=convert(series(sin(x), x=0, 32), polynom):
 polycos:=convert(series(cos(x), x=0, 32), polynom):
 log[2](numapprox[infnorm](1 - polysin/sin(x), x=0..Pi/4));
 log[2](numapprox[infnorm](1 - polysin/sin(x), x=0..2^(-17)));

 log[2](numapprox[infnorm](1 - polycos/cos(x), x=0..Pi/4));
 log[2](numapprox[infnorm](1 - polycos/cos(x), x=0..2^(-18)));

 polytan:=convert(series(tan(x), x=0, 128), polynom):
 log[2](numapprox[infnorm](  (tan(x) - polytan)/tan(x), x=0..Pi/4));
 log[2](numapprox[infnorm](  (tan(x) - polytan)/tan(x), x=0..2^(-18)));



 ftany := tan(sqrt(y))/sqrt(y^3)-1/y;
 ymin:=2^(-1075);
 ymax:=Pi^2/16;

 # to avoid pb of singularity in the minimax, work on a Taylor approx
 pTanTaylor := convert(series(ftany, y=0, 150), polynom):
 maxepsApprox1 := abs(evalf(eval(1-pTanTaylor/ftany, y=ymax))); # the max error is there
 log2(maxepsApprox1);

 pErrTan := convert(series( y^(3/2)/tan(sqrt(y)), y=0, 3), polynom);
 ptanscs_y := numapprox[minimax](pTanTaylor, y=ymin..ymax,[5,0], pErrTan,'err'    );
log2(err);
abs(evalf(eval(1-convert(series(ftany, y=0, 5), polynom)/ftany, y=ymax)))
 ptanscs_y := numapprox[minimax](ftany,
                                 y=ymin..ymax,
                                 [4,0],
                                 pErrTan,
                                 'err'    );
log2(err);

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

Poly_P    := series(tan(sqrt(x))/(x^(3/2))-1/x, x=0, 40):
Poly_Q    := convert(Poly_P,polynom):
Poly_cheb := chebpade(Poly_Q, x=0..evalf(Pi/4), [17,0]);
Poly_Res  := x + x^3 * subs(x=x^2, Poly_cheb);
log2(infnorm( 1 - (Poly_Res)/tan(x),x=0..evalf(Pi/4)));

fi:
