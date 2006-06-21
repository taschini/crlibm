#######################################################################
# This file is part of the crlibm library, and is distributed under
# the  LGPL.
# To use:
# restart; read "exp-td.mpl";
Digits := 120:

interface(quiet=true):

read "common-procedures.mpl":
read "triple-double.mpl":
mkdir("TEMPEXPM1"):


printPolynomialIntoFile := proc(fd,s,p) 
local i, hi, mi, lo:
for i from 0 to degree(p(x),x) do
	(hi,mi,lo) := hi_mi_lo(coeff(p(x),x,i)):
	if ((abs(hi) = 1.0) and (mi = 0) and (lo = 0)) then 
		printf(
		"Coefficient %d of the polynomial is exactly %f and will not be stored in the table\n",i,hi): 
	else 
	if ((abs(hi) = 0.5) and (mi = 0) and (lo = 0)) then 
		printf(
		"Coefficient %d of the polynomial is exactly %f and will not be stored in the table\n",i,hi): 
	else 
	if (hi <> 0) then
 		fprintf(fd,"#define %s%dh %1.50e\n",s,i,hi):
	end if:
	if (mi <> 0) then
	 	fprintf(fd,"#define %s%dm %1.50e\n",s,i,mi):
	end if:
	if (lo <> 0) then
	 	fprintf(fd,"#define %s%dl %1.50e\n",s,i,lo):
	end if:
	end if:
	end if:
od:
end proc:


# First, we compute special values 

ReturnXBound := convert((ieeehexa(2^(-54)))[1],decimal,hex):

Largest := 2^(1023) * ((2^(53) - 1) / 2^(52)):

Smallest := 2^(-1023) * 1 * 2^(-51):

OverflowBound := nearest(log(Largest + 1)):

MinusOneBound := nearest(log(2^(-54))):

SimpleOverflowBound := convert(ieeehexa(OverflowBound)[1],decimal,hex):

DirectIntervalBound := convert((ieeehexa(0.25))[1],decimal,hex):

MinusOnePlusOneUlp := -1 + 2^(-53): # Attention: it's 2^(-53) because we are at a binade boundary

# Second, we have the computation of the values for the direct interval

# The function, that we approximate is 

directF := unapply(exp(x) - 1,x):

# The domain is 

directA := -2^(-5):
directB := 2^(-5):

# The polynomials are

quickDirectpoly := X -> X+1/2*X^2+(3360682229480701/1180591620717411303424*X^6+3660136839517697/147573952589676412928*X^5+7320130809407439/36893488147419103232*X^4+3202559734508631/2305843009213693952*X^3+4803839602572223/576460752303423488*X^2+6004799503160665/144115188075855872*X+6004799503160661/36028797018963968)*X^3:

accuDirectpoly := X -> X+1/2*X^2+(3786738884990361/4951760157141521099596496896*X^12+7100145222887513/618970019642690137449562112*X^11+6212541673969101/38685626227668133590597632*X^10+5047690109993399/2417851639229258349412352*X^9+3785767582868083/151115727451828646838272*X^8+5205430426443615/18889465931478580854784*X^7+29303968161043118891149009244865/10633823966279326983230456482242756608*X^6+65933928362347017505024866986963/2658455991569831745807614120560689152*X^5+65933928362347017505149159875899/332306998946228968225951765070086144*X^4+28846093658526820158502757550845/20769187434139310514121985316880384*X^3+21634570243895115118877068038417/2596148429267413814265248164610048*X^2+27043212804868893898596335048021/649037107316853453566312041152512*X+243583606221817153033947472119380503276473908509/1461501637330902918203684832716283019655932542976)*X^3:


# Truncate the quick phase direct interval polynomial to degree specialDegree 
# for special interval |x| <= specialBound (speed-up)

specialDegree := 5:
specialBound := 2^(-12):

specialPoly := unapply(sum(coeff(quickDirectpoly(x),x,i) * x^i,i=0..specialDegree),x):

printf("Special polynomial is the direct polynomial truncated to degree %d used in |x| < 2^(%f)\n",
	specialDegree, evalf(log[2](specialBound))):

# Compute the relative errors

errDirectQuick := numapprox[infnorm](quickDirectpoly(x)/directF(x) -1,x=directA..directB):
errDirectAccu := numapprox[infnorm](accuDirectpoly(x)/directF(x) -1,x=directA..directB):

errSpecialPoly := numapprox[infnorm](specialPoly(x)/directF(x) -1,x=-specialBound..specialBound):

errDirectAccuSpecial := numapprox[infnorm](accuDirectpoly(x)/directF(x) -1,x=2^(-12)..2^(-12)):

printf("The relative approximation error of the direct interval quick polynomial is 2^(%f)\n",
	evalf(log[2](abs(errDirectQuick)))):
printf("The relative approximation error of the direct interval accurate polynomial is 2^(%f)\n",
	evalf(log[2](abs(errDirectAccu)))):
printf("The relative approximation error of the special interval special polynomial is 2^(%f)\n",
	evalf(log[2](abs(errSpecialPoly)))):
printf("The relative approximation error of the direct interval accurate polynomial in special domain is 2^(%f)\n",
	evalf(log[2](abs(errDirectAccuSpecial)))):



# Third, we have the computation of the values for the common interval

# The function, that we approximate is 

commonF := unapply(exp(x),x):

# The domain is 

commonA := -log(2)*2^(-12) * (1/2 + 2^(-19)):
commonB := log(2)*2^(-12) * (1/2 + 2^(-19)):


quickCommonpoly := X -> 1+X+1/2*X^2+(6004799504593679/144115188075855872*X+6004799504235425/36028797018963968)*X^3:

accuCommonpoly := X -> 1+X+1/2*X^2+(3660068549402285/18446744073709551616*X^4+6405119471061623/4611686018427387904*X^3+4803839602528529/576460752303423488*X^2+54086425609737787796676993069745/1298074214633706907132624082305024*X+54086425609737787797192670135537/324518553658426726783156020576256)*X^3:

# Compute the relative errors

errCommonQuick := numapprox[infnorm](quickCommonpoly(x)/commonF(x) -1,x=commonA..commonB):
errCommonAccu := numapprox[infnorm](accuCommonpoly(x)/commonF(x) -1,x=commonA..commonB):

printf("The relative approximation error of the common interval quick polynomial is 2^(%f)\n",
	evalf(log[2](abs(errCommonQuick)))):
printf("The relative approximation error of the common interval accurate polynomial is 2^(%f)\n",
	evalf(log[2](abs(errCommonAccu)))):

epsilonApproxRmAccurate := numapprox[infnorm]( ((1+x)/(exp(x)))-1, x=commonA*2^(-52)..commonB*2^(-52)):
epsilonApproxRlAccurate := numapprox[infnorm]( ((1+x)/(exp(x)))-1, x=commonA*2^(-105)..commonB*2^(-105)):

printf("The approximation rel error for approximating exp(rm) by 1 + rm is 2^(%2f)\n", 
	log2(abs(epsilonApproxRmAccurate))):
printf("The approximation rel error for approximating exp(rl) by 1 + rl is 2^(%2f)\n", 
	log2(abs(epsilonApproxRlAccurate))):



# Compute the constants for argument reduction and the tables in the common path

MsLog2Div2L := evalf(-log(2)/(2^(12))):

msLog2Div2Lh, msLog2Div2Lm, msLog2Div2Ll := hi_mi_lo(MsLog2Div2L):

epsMsLog2Div2L := evalf(abs(((msLog2Div2Lh + msLog2Div2Lm + msLog2Div2Ll) - MsLog2Div2L)/MsLog2Div2L)):
epsDDMsLog2Div2L := evalf(abs(((msLog2Div2Lh + msLog2Div2Lm) - MsLog2Div2L)/MsLog2Div2L)):

printf("The error made by storing MsLog2Div2L as a double-double is 2^(%f)\n",log[2](epsDDMsLog2Div2L)):
printf("The error made by storing MsLog2Div2L as a triple-double is 2^(%f)\n",log[2](epsMsLog2Div2L)):

gap := -floor(-log[2](abs(msLog2Div2Lm/msLog2Div2Lh))):

printf("Information: |msLog2Div2Lm| <= 2^(%f) * |msLog2Div2Lh|\n",gap):


log2InvMult2L := nearest(2^(12) / (log(2))):

shiftConst := 2^(52) + 2^(51):

indexmask1 := 2^((12)/2) - 1:
indexmask2 := indexmask1 * 2^((12)/2):

for i from 0 to 2^(12/2) - 1 do
	twoPowerIndex1hi[i], twoPowerIndex1mi[i], twoPowerIndex1lo[i] := hi_mi_lo(evalf(2^(i/(2^12)))):
	twoPowerIndex2hi[i], twoPowerIndex2mi[i], twoPowerIndex2lo[i] := hi_mi_lo(evalf(2^(i/(2^(12/2))))):
od: 

# Estimate the error of the two quick phases 

# ATTENTION: C EST PIFOMETRIQUE POUR L INSTANT

epsQuickDirectOverall := 2^(-62):
epsQuickCommonOverall := 2^(-62):



# Write the tables

printf("Write tables...\n"):

filename:="TEMPEXPM1/expm1.h":
fd:=fopen(filename, WRITE, TEXT):

fprintf(fd, "#include \"crlibm.h\"\n#include \"crlibm_private.h\"\n"):

fprintf(fd, "\n/* File generated by maple/expm1.mpl */\n"):

fprintf(fd, "\#define log2InvMult2L %1.50e\n",log2InvMult2L):
fprintf(fd, "\#define msLog2Div2Lh %1.50e\n",msLog2Div2Lh):
fprintf(fd, "\#define msLog2Div2Lm %1.50e\n",msLog2Div2Lm):
fprintf(fd, "\#define msLog2Div2Ll %1.50e\n",msLog2Div2Ll):
fprintf(fd, "\#define shiftConst %1.50e\n",shiftConst):
fprintf(fd, "\#define INDEXMASK1 0x%08x\n",indexmask1):
fprintf(fd, "\#define INDEXMASK2 0x%08x\n",indexmask2):
fprintf(fd, "\#define RETURNXBOUND 0x%08x\n",ReturnXBound):
fprintf(fd, "\#define OVERFLOWBOUND %1.50e\n",OverflowBound):
fprintf(fd, "\#define LARGEST %1.50e\n",Largest): 
fprintf(fd, "\#define SMALLEST %1.50e\n",Smallest): 
fprintf(fd, "\#define MINUSONEBOUND %1.50e\n",MinusOneBound):
fprintf(fd, "\#define SIMPLEOVERFLOWBOUND 0x%08x\n",SimpleOverflowBound):
fprintf(fd, "\#define DIRECTINTERVALBOUND 0x%08x\n",DirectIntervalBound):
fprintf(fd, "\#define SPECIALINTERVALBOUND 0x%08x\n",convert((ieeehexa(specialBound))[1],decimal,hex)):
fprintf(fd, "\#define ROUNDCSTDIRECTRN %1.50e\n",compute_rn_constant(epsQuickDirectOverall)):
fprintf(fd, "\#define ROUNDCSTDIRECTRD %1.50e\n",epsQuickDirectOverall):
fprintf(fd, "\#define ROUNDCSTCOMMONRN %1.50e\n",compute_rn_constant(epsQuickCommonOverall)):
fprintf(fd, "\#define ROUNDCSTCOMMONRD %1.50e\n",epsQuickCommonOverall):
fprintf(fd, "\#define MINUSONEPLUSONEULP %1.50e\n",MinusOnePlusOneUlp):



fprintf(fd,"\n\n"):

printPolynomialIntoFile(fd,"quickDirectpolyC",quickDirectpoly):
fprintf(fd,"\n"):
printPolynomialIntoFile(fd,"accuDirectpolyC",accuDirectpoly):
fprintf(fd,"\n"):
printPolynomialIntoFile(fd,"quickCommonpolyC",quickCommonpoly):
fprintf(fd,"\n"):
printPolynomialIntoFile(fd,"accuCommonpolyC",accuCommonpoly):


fprintf(fd,"\n\n"):

# Print the tables
fprintf(fd, "typedef struct tPi_t_tag {double hi; double mi; double lo;} tPi_t;  \n"):
fprintf(fd, "static const tPi_t twoPowerIndex1[%d] = {\n", 2^(12/2)):
for i from 0 to 2^(12/2)-1 do
      fprintf(fd, "  { \n"):      
      fprintf(fd, "    %1.50e, /* twoPowerIndex1hi[%d] */ \n", twoPowerIndex1hi[i], i):
      fprintf(fd, "    %1.50e, /* twoPowerIndex1mi[%d] */ \n", twoPowerIndex1mi[i], i):
      fprintf(fd, "    %1.50e, /* twoPowerIndex1lo[%d] */ \n", twoPowerIndex1lo[i], i):
      fprintf(fd, "  } "):
      if(i<2^(12/2)-1) then  fprintf(fd, ", \n"): fi
od:
fprintf(fd, "}; \n \n"):
fprintf(fd, "static const tPi_t twoPowerIndex2[%d] = {\n", 2^(12/2)):
for i from 0 to 2^(12/2)-1 do
      fprintf(fd, "  { \n"):      
      fprintf(fd, "    %1.50e, /* twoPowerIndex2hi[%d] */ \n", twoPowerIndex2hi[i], i):
      fprintf(fd, "    %1.50e, /* twoPowerIndex2mi[%d] */ \n", twoPowerIndex2mi[i], i):
      fprintf(fd, "    %1.50e, /* twoPowerIndex2lo[%d] */ \n", twoPowerIndex2lo[i], i):
      fprintf(fd, "  } "):
      if(i<2^(12/2)-1) then  fprintf(fd, ", \n"): fi
od:
fprintf(fd, "}; \n \n"):

fprintf(fd, "\n\n"):

fclose(fd):

printf("               ...done\n"):