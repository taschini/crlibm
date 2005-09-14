#######################################################################
# This file is part of the crlibm library, and is distributed under
# the  LGPL.
# To use:
# restart; read "log-de.mpl";
Digits := 120:

interface(quiet=true):

read "common-procedures.mpl":
read "triple-double.mpl":
mkdir("TEMPLOG"):


# We want log2h + log2m + log2l + delta = log(2) such that
# log2h and log2m have at least 11 trailing zeros
# in order to have an exact multiplication with E, which is lower than 1024 in
# magnitude
# The resting accuracy is enough for both quick and accurate phases.

log2acc := log(2):
log2h := round(log2acc * 2**(floor(-log[2](abs(log2acc))) + (53 - 11))) / 2**(floor(-log[2](abs(log2acc))) + (53 - 11)):
log2m := round((log2acc - log2h) * 2**(floor(-log[2](abs((log2acc - log2h)))) + (53 - 11))) / 2**(floor(-log[2](abs((log2acc - log2h)))) + (53 - 11)):
log2l := log2acc - (log2h + log2m):


L := 7: # number of bits used to address the table

MAXINDEX    := round(2^L * (sqrt(2)-1)):

for i from 0 to MAXINDEX-1 do
    center[i] := 1 + i*2^(-L): # center[i] in [1, 2[
    t :=  evalf(1/center[i]):
    r[i] := round(t * 2**(floor(-log[2](abs(t))) + 23)) / 2**(floor(-log[2](abs(t))) + 23):
    (logih[i], logim[i], logil[i]) := hi_mi_lo(evalf(-log(r[i]))):
od:
for i from MAXINDEX to 2^L do
    # y has been divided by two, center[i] in [0.5, 1[
    center[i]:=(1 + i*2^(-L)) / 2:
    t :=  evalf(1/center[i]):
    r[i] := round(t * 2**(floor(-log[2](abs(t))) + 23)) / 2**(floor(-log[2](abs(t))) + 23):
    (logih[i], logim[i], logil[i]) := hi_mi_lo(evalf(-log(r[i]))):
od:




#Computation of ZMax.
for i from 0 to MAXINDEX-1 do
    __x := center[i] + 2^(-L-1) :
    zmax[i] := (__x*r[i]-1) :
    __x := center[i] - 2^(-L-1) :
    zmin[i] := (__x*r[i]-1) :
od:
for i from MAXINDEX to 2^L do
    __x := center[i] + 2^(-L-2) :
    zmax[i] := (__x*r[i]-1) :
    __x := center[i] - 2^(-L-2) :
    zmin[i] := (__x*r[i]-1) :
od:

zmaxmax:=0:
zminmin:=0:
for i from 0 to 2^L do
    tabulated_value := logih[i] + logim[i] + logil[i]:
    poly_approx_min := evalf(log(1+zmin[i])):
    poly_approx_max := evalf(log(1+zmax[i])):

    # Test if we have a case where we cancellate a lot
    # i.e. the polynomial approximation value is greater than half the tabulated value
    # the tabulated value is not exactly zero and we are of opposite sign

    if ((abs(poly_approx_min) > 0.75*abs(tabulated_value)) and (tabulated_value <> 0.0) and (poly_approx_min * tabulated_value < 0)) then
	printf("Polynomial approximation is greater in magnitude in zmin[%d] than half the tabluated value\n",i):

	printf("The tabulated value is %1.50e\n",tabulated_value):
	if (tabulated_value <> 0.0) then printf("i.e. the value has the exponent %d\n",floor(log2(abs(tabulated_value)))) fi:
	printf("The value of polynomial in zmin[%d] is %1.50e\n",i,poly_approx_min):
        if (poly_approx_min <> 0.0) then printf("i.e. the value has the exponent %d\n",floor(log2(abs(poly_approx_min)))) fi:

	summe := poly_approx_min + tabulated_value:
	printf("The exponent of the sum of both is %d\n",floor(log2(abs(summe)))):


    fi:

    if ((abs(poly_approx_max) > 0.75*abs(tabulated_value)) and (tabulated_value <> 0.0) and (poly_approx_max * tabulated_value <0)) then
	printf("Polynomial approximation is greater in magnitude in zmax[%d] than half the tabluated value\n",i):

	printf("The tabulated value is %1.50e\n",tabulated_value):
	if (tabulated_value <> 0.0) then printf("i.e. the value has the exponent %d\n",floor(log2(abs(tabulated_value)))) fi:
	printf("The value of polynomial in zmax[%d] is %1.50e\n",i,poly_approx_max):
        if (poly_approx_max <> 0.0) then printf("i.e. the value has the exponent %d\n",floor(log2(abs(poly_approx_max)))) fi:

	summe := poly_approx_max + tabulated_value:
	printf("The exponent of the sum of both is %d\n",floor(log2(abs(summe)))):


    fi:


    if zmax[i] > zmaxmax then zmaxmax := zmax[i]: fi:
    if zmin[i] < zminmin then zminmin := zmin[i]: fi:
od:
printf("zminmin = -2^(%2f)   zmaxmax = 2^(%2f)\n", log2(-zminmin), log2(zmaxmax) ) :

PolyDegreeQuick:=7:

printf("   degree of the polynomial used in the quick phase is %d\n",PolyDegreeQuick);

DDNumberQuick:=3:

printf("   number of double doubles used for the coefficients is %d\n",DDNumberQuick);

#Keep -zmaxmax..zmaxmax to keep c1=1, which is useful in the proof
#and constrain the first two coefficients to 1 and -1/2 in order to save up a full multiplication and a rounding error
polyQuick:= poly_exact2(x*(1+x*(-0.5+x*(numapprox[minimax]((((log(1+x)/x)-1)/x+0.5)/x,
			x=-zmaxmax..zmaxmax,  [PolyDegreeQuick-3,0], 1 ,  'deltaApprox')))), DDNumberQuick):

#Try to verify the bound for using double double arithmetic.
#Idea: compare the maximum absolute value of the polynomial in zmaxmax (the polynomial has its maxima at the limits)
#with the maximum value of the first term which is calculated in double precision only

p := unapply(polyQuick,x):
printf("   using only %d double doubles should be fine since the hi z ulp should affect the result starting from bit %f\n",
	DDNumberQuick,evalf(53 + log2(p(zmaxmax)) - log2(zmaxmax^(DDNumberQuick)),5)):


epsilonApproxQuick := numapprox[infnorm]( 1-polyQuick/log(1+x), x=zminmin..zmaxmax):
printf("   approximation rel error for the quick phase is 2^(%2f)\n", log2(epsilonApproxQuick) ) :
deltaApproxQuick := numapprox[infnorm]( polyQuick-log(1+x), x=zminmin..zmaxmax):
printf("   approximation abs error for the quick phase is 2^(%2f)\n", log2(deltaApproxQuick) ) :


PolyDegreeAccurate:=14:

printf("   degree of the polynomial used in the accurate phase is %d\n",PolyDegreeAccurate);

DDNumberAccu:=7:
TDNumberAccu:=3:

printf("   number of triple doubles used for the coefficients is %d\n",TDNumberAccu);
printf("   number of double doubles used for the coefficients is %d\n",DDNumberAccu);


#Keep -zmaxmax..zmaxmax to keep c1=1, which is useful in the proof
polyAccurate:= poly_exact32(x*(1+x*(-0.5+x*(numapprox[minimax]((((log(1+x)/x)-1)/x+0.5)/x,
				x=-zmaxmax..zmaxmax,  [PolyDegreeAccurate-3,0], 1 ,  'deltaApprox')))),
				TDNumberAccu, DDNumberAccu):

#Try to verify the bound for using double double arithmetic.
#Idea: compare the maximum absolute value of the polynomial in zmaxmax (the polynomial has its maxima at the limits)
#with the maximum value of the first term which is calculated in double precision only

pp := unapply(polyAccurate,x):
printf("   using only %d triple doubles should be fine since the mi z ulp should affect the result starting from bit %f\n",
	TDNumberAccu,evalf(106 + log2(p(zmaxmax)) - log2(zmaxmax^(TDNumberAccu)),5)):
printf("   using only %d double doubles should be fine since the hi z ulp should affect the result starting from bit %f\n",
	DDNumberAccu,evalf(53 + log2(p(zmaxmax)) - log2(zmaxmax^(TDNumberAccu + DDNumberAccu)),5)):


epsilonApproxAccurate := numapprox[infnorm]( 1-polyAccurate/log(1+x), x=zminmin..zmaxmax):
printf("   approximation rel error for the accurate phase is 2^(%2f)\n", log2(epsilonApproxAccurate) ) :
deltaApproxAccurate := numapprox[infnorm]( polyAccurate-log(1+x), x=zminmin..zmaxmax):
printf("   approximation abs error for the quick phase is 2^(%2f)\n", log2(deltaApproxAccurate) ) :



#-------------------------------------------------------------------
# Output


filename:="TEMPLOG/log-td.h":
fd:=fopen(filename, WRITE, TEXT):

fprintf(fd, "#include \"crlibm.h\"\n#include \"crlibm_private.h\"\n"):

fprintf(fd, "\n/*File generated by maple/log-td.mpl*/\n"):

fprintf(fd, "\n\#define L %d\n\n",L):
fprintf(fd, "\#define MAXINDEX %d\n\n", MAXINDEX):
fprintf(fd, "\#define INDEXMASK %d\n", 2^L-1):
fprintf(fd, "\#define two52 %1.50e\n", 2^(52)):
fprintf(fd, "\#define log2h %1.50e\n", log2h):
fprintf(fd, "\#define log2m %1.50e\n", log2m):
fprintf(fd, "\#define log2l %1.50e\n", log2l):

epsilon_quick_1 := 2^(-62.5): # The Gappa proof will show this bound
epsilon_quick_2 := 2^(-62.5): # The Gappa proof will show this bound
fprintf(fd, "\#define ROUNDCST1 %1.50e\n", compute_rn_constant(epsilon_quick_1)):
fprintf(fd, "\#define ROUNDCST2 %1.50e\n", compute_rn_constant(epsilon_quick_2)):
fprintf(fd, "\#define RDROUNDCST1 %1.50e\n", epsilon_quick_1):
fprintf(fd, "\#define RDROUNDCST2 %1.50e\n", epsilon_quick_2):


fprintf(fd, "\n\n"):


# Print the defines for the define statements

for i from 3 to PolyDegreeQuick do
	fprintf(fd, "\#define c%d %1.50e\n",i,coeff(polyQuick,x,i)):
od:

fprintf(fd, "\n\n"):

for i from 3 to (DDNumberAccu + TDNumberAccu -1) do
	(hi,lo) := hi_lo(coeff(polyAccurate,x,i)):
	fprintf(fd, "\#define accPolyC%dh %1.50e\n",i,hi):
	fprintf(fd, "\#define accPolyC%dl %1.50e\n",i,lo):
od:

for i from (DDNumberAccu + TDNumberAccu) to PolyDegreeAccurate do
	fprintf(fd, "\#define accPolyC%d %1.50e\n",i,coeff(polyAccurate,x,i)):
od:

fprintf(fd, "\n\n"):

# Print the table
fprintf(fd, "typedef struct rri_tag {float ri; double logih; double logim; double logil;} rri;  \n"):
fprintf(fd, "static const rri argredtable[%d] = {\n", 2^L):
for i from 0 to 2^L-1 do
      fprintf(fd, "  { \n"):
      fprintf(fd, "    %1.50e,   /* r[%d] */ \n", r[i], i):
      fprintf(fd, "    %1.50e, /* logih[%d] */ \n", logih[i], i):
      fprintf(fd, "    %1.50e, /* logim[%d] */ \n", logim[i], i):
      fprintf(fd, "    %1.50e, /* logil[%d] */ \n", logil[i], i):
      fprintf(fd, "  } "):
      if(i<2^L-1) then  fprintf(fd, ", \n"): fi
od:
fprintf(fd, "}; \n \n"):

fclose(fd):

for j from 0 to 2^L-1 do
    filename:=cat("TEMPLOG/log-td_",j,".sed"):
    fd:=fopen(filename, WRITE, TEXT):
    fprintf(fd, " s/_log2h/%1.50e/g\n", log2h):
    fprintf(fd, " s/_log2m/%1.50e/g\n", log2m):
    fprintf(fd, " s/_log2l/%1.50e/g\n", log2l):
    fprintf(fd, " s/_logih/%1.50e/g\n", logih[j]):
    fprintf(fd, " s/_logim/%1.50e/g\n", logim[j]):
    fprintf(fd, " s/_logil/%1.50e/g\n", logil[j]):
    fprintf(fd, " s/_zmin/%1.50e/g\n", zmin[j]):
    fprintf(fd, " s/_zmax/%1.50e/g\n", zmax[j]):
    for i from 3 to PolyDegreeQuick do
        fprintf(fd, " s/_c%d/%1.50e/g\n", i, coeff(polyQuick,x,i)):
    od:
    fprintf(fd, " s/_epsilonApproxQuick/%1.50e/g\n", epsilonApproxQuick):
    fclose(fd):
  od:

for j from 0 to 2^L-1 do
    filename:=cat("TEMPLOG/log-td-accurate_",j,".sed"):
    fd:=fopen(filename, WRITE, TEXT):
    fprintf(fd, " s/_log2h/%1.50e/g\n", log2h):
    fprintf(fd, " s/_log2m/%1.50e/g\n", log2m):
    fprintf(fd, " s/_log2l/%1.50e/g\n", log2l):
    fprintf(fd, " s/_logih/%1.50e/g\n", logih[j]):
    fprintf(fd, " s/_logim/%1.50e/g\n", logim[j]):
    fprintf(fd, " s/_logil/%1.50e/g\n", logil[j]):
    fprintf(fd, " s/_zmin/%1.50e/g\n", zmin[j]):
    fprintf(fd, " s/_zmax/%1.50e/g\n", zmax[j]):
    for i from 3 to (DDNumberAccu + TDNumberAccu -1) do
	(hi,lo) := hi_lo(coeff(polyAccurate,x,i)):
        fprintf(fd, " s/_accPolyC%dh/%1.50e/g\n", i, hi):
        fprintf(fd, " s/_accPolyC%dl/%1.50e/g\n", i, lo):
    od:
    for i from (DDNumberAccu + TDNumberAccu) to PolyDegreeAccurate do
        fprintf(fd, " s/_accPolyC%d/%1.50e/g\n", i, coeff(polyAccurate,x,i)):
    od:
    fprintf(fd, " s/_epsilonApproxAccurate/%1.50e/g\n", epsilonApproxAccurate):
    fclose(fd):
  od:

# A shell script to use them
filename:="run-log-td-proof.sh":
fd:=fopen(filename, WRITE, TEXT):
fprintf(fd, "#!/bin/sh\n"):
fprintf(fd, "# You probably need to edit the path to the gappa executable\n"):
fprintf(fd, "GAPPA=~/sangria/gappa-0.4.7/src/gappa\n"):
fprintf(fd, "# Test all the possible table value for E=1\n"):
fprintf(fd, "for num in `seq 0 %d`; do\n", 2^L-1):
fprintf(fd, "  echo $num, E=1:\n"):
fprintf(fd, "  sed -f ./TEMPLOG/log-td_$num.sed log-td.gappa | $GAPPA > /dev/null\n"):
fprintf(fd, "  echo\n"):
fprintf(fd, "done\n"):
fprintf(fd, "# For the case E=0 we first handle the cases 0 and %d using log-td-E0-logir0.gappa\n", 2^L):
fprintf(fd, "echo 0 and %d, E=0:\n", 2^L):
fprintf(fd, "sed -f log-td_0.sed log-td-E0-logir0.gappa | $GAPPA > /dev/null\n"):
fprintf(fd, "# then the other cases where logirh <> 0\n"):
fprintf(fd, "for num in `seq 1 %d`; do\n", 2^L-1):
fprintf(fd, "  echo $num, E=0:\n"):
fprintf(fd, "  sed -f ./TEMPLOG/log-td_$num.sed log-td-E0.gappa | $GAPPA > /dev/null\n"):
fprintf(fd, "  echo\n"):
fprintf(fd, "done\n"):
fprintf(fd, "# Accurate phase: Test all the possible table value for E=1\n"):
fprintf(fd, "for num in `seq 0 %d`; do\n", 2^L-1):
fprintf(fd, "  echo Accurate phase: $num, E=1:\n"):
fprintf(fd, "  sed -f ./TEMPLOG/log-td-accurate_$num.sed log-td-accurate.gappa | $GAPPA > /dev/null\n"):
fprintf(fd, "  echo\n"):
fprintf(fd, "done\n"):
fprintf(fd, "# Accurate phase: For the case E=0 we first handle the cases 0 and %d using log-td-accurate-E0-logir0.gappa\n", 2^L):
fprintf(fd, "echo 0 and %d, E=0:\n", 2^L):
fprintf(fd, "sed -f ./TEMPLOG/log-td-accurate_0.sed log-td-accurate-E0-logir0.gappa | $GAPPA > /dev/null\n"):
fprintf(fd, "# Accurate phase: then the other cases where logirh <> 0\n"):
fprintf(fd, "for num in `seq 1 %d`; do\n", 2^L-1):
fprintf(fd, "  echo $num, E=0:\n"):
fprintf(fd, "  sed -f ./TEMPLOG/log-td-accurate_$num.sed log-td-accurate-E0.gappa | $GAPPA > /dev/null\n"):
fprintf(fd, "  echo\n"):
fprintf(fd, "done\n"):
fclose(fd):

printf("----DONE---\n") :


