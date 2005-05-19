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
# in order to have an exact multiplication with E, which is lower that 1024 in
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
    _x := center[i] + 2^(-L-1) :
    zmax[i] := (_x*r[i]-1) :
    _x := center[i] - 2^(-L-1) :
    zmin[i] := (_x*r[i]-1) :
od:
for i from MAXINDEX to 2^L do
    _x := center[i] + 2^(-L-2) :
    zmax[i] := (_x*r[i]-1) :
    _x := center[i] - 2^(-L-2) :
    zmin[i] := (_x*r[i]-1) :
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
# + force les premiers deux coefficients a 1 et 1/2 pour economiser une multiplication et 
# une erreur de multiplication
polyQuick:= poly_exact2(x*(1+x*(-0.5+x*(numapprox[minimax]((((log(1+x)/x)-1)/x+0.5)/x,  
			x=-zmaxmax..zmaxmax,  [PolyDegreeQuick-3,0], 1 ,  'deltaApprox')))), DDNumberQuick):

#On essaie de verfier la borne pour l'utilisation de la double double.
#Idee: comparer la valeur maximale absolue du polynome en zmaxmax (le polynome a ces maxima aux bords) 
#avec la valeur maximale du premier terme seulement calcule en double precision
#Est-ce bien correct?

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

#On essaie de verfier la borne pour l'utilisation de la double double.
#Idee: comparer la valeur maximale absolue du polynome en zmaxmax (le polynome a ces maxima aux bords) 
#avec la valeur maximale du premier terme seulement calcule en double precision
#Est-ce bien correct?

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
fprintf(fd, "\#define two52 %1.50eL\n", 2^(52)):
fprintf(fd, "\#define log2h %1.50eL\n", log2h):
fprintf(fd, "\#define log2m %1.50eL\n", log2m):
fprintf(fd, "\#define log2l %1.50eL\n", log2l):

epsilon_quick_1 := 2^(-62.5): # ATTENTION: c est completement pifometrique pour l instant
epsilon_quick_2 := 2^(-62.5): # ATTENTION: c est completement pifometrique pour l instant
fprintf(fd, "\#define ROUNDCST1 %1.50eL\n", compute_rn_constant(epsilon_quick_1)):   
fprintf(fd, "\#define ROUNDCST2 %1.50eL\n", compute_rn_constant(epsilon_quick_2)):   
fprintf(fd, "\#define RDROUNDCST1 %1.50eL\n", epsilon_quick_1):   
fprintf(fd, "\#define RDROUNDCST2 %1.50eL\n", epsilon_quick_2):   


fprintf(fd, "\n\n"):


# Print the defines for the define statements 

for i from 3 to PolyDegreeQuick do
	fprintf(fd, "\#define c%d %1.50eL\n",i,coeff(polyQuick,x,i)):
od:

fprintf(fd, "\n\n"):

for i from 3 to (DDNumberAccu + TDNumberAccu -1) do
	(hi,lo) := hi_lo(coeff(polyAccurate,x,i)):
	fprintf(fd, "\#define accPolyC%dh %1.50eL\n",i,hi):
	fprintf(fd, "\#define accPolyC%dl %1.50eL\n",i,lo):
od:

for i from (DDNumberAccu + TDNumberAccu) to PolyDegreeAccurate do
	fprintf(fd, "\#define accPolyC%d %1.50eL\n",i,coeff(polyAccurate,x,i)):
od:

fprintf(fd, "\n\n"):

# Print the table
fprintf(fd, "typedef struct rri_tag {float ri; double logih; double logim; double logil;} rri;  \n"):
fprintf(fd, "static const rri argredtable[%d] = {\n", 2^L):
for i from 0 to 2^L-1 do
      fprintf(fd, "  { \n"):      
      fprintf(fd, "    %1.50eL,   /* r[%d] */ \n", r[i], i):
      fprintf(fd, "    %1.50eL, /* logih[%d] */ \n", logih[i], i):
      fprintf(fd, "    %1.50eL, /* logim[%d] */ \n", logim[i], i):
      fprintf(fd, "    %1.50eL, /* logil[%d] */ \n", logil[i], i):
      fprintf(fd, "  } "):
      if(i<2^L-1) then  fprintf(fd, ", \n"): fi
od:
fprintf(fd, "}; \n \n"):

fclose(fd):

printf("----DONE---\n") :


