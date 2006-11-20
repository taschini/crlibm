#######################################################################
# This file is part of the crlibm library, and is distributed under
# the  LGPL.
# To use:
# restart; read "exp-td.mpl";
Digits := 145:

interface(quiet=true):

read "common-procedures.mpl":
read "triple-double.mpl":
read "gal.mpl":
mkdir("TEMPASIN"):

with(orthopoly,T):


truncPoly := proc(p,k) 
local i, q:
q := 0:
convert(q,polynom):
for i from 0 to min(degree(p,x),k) do
    q := q + coeff(p,x,i) * x^i:
od:
return (q):
end:

intervals := 10:

for i from 1 to intervals do
	epsAccurate[i] := 0.5:
	epsQuick[i] := 0.5:
	polyAccurate[i] := 1:
	polyQuick[i] := 1:
od:
epsAccurateSpecial := 0.5:
epsQuickExtra := 0.5:

lowestIntervalMax := 0.185:
highestIntervalMin := 0.78:

polyAccurateTDCoeffsLowest := 13:
polyAccurateDDCoeffsLowest := 12:
polyAccurateDCoeffsLowest := 16:

polyAccurateDegreeLowest := polyAccurateTDCoeffsLowest + polyAccurateDDCoeffsLowest + polyAccurateDCoeffsLowest:

polyQuickDDCoeffsLowest := 12:
polyQuickDCoeffsLowest := 10:

polyQuickDegreeLowest := polyQuickDDCoeffsLowest + polyQuickDCoeffsLowest:

polyAccurateTDCoeffsMiddle := 7:
polyAccurateDDCoeffsMiddle := 9:
polyAccurateDCoeffsMiddle := 19:

polyAccurateDegreeMiddle := polyAccurateTDCoeffsMiddle + polyAccurateDDCoeffsMiddle + polyAccurateDCoeffsMiddle:

polyQuickDDCoeffsMiddle := 7:
polyQuickDCoeffsMiddle := 7:

polyQuickDegreeMiddle := polyQuickDDCoeffsMiddle + polyQuickDCoeffsMiddle:


polyAccurateTDCoeffsHighest := 9:
polyAccurateDDCoeffsHighest := 9:
polyAccurateDCoeffsHighest := 11:

polyAccurateDegreeHighest := polyAccurateTDCoeffsHighest + polyAccurateDDCoeffsHighest + polyAccurateDCoeffsHighest:

polyQuickDDCoeffsHighest := 9:
polyQuickDCoeffsHighest := 9:

polyQuickDegreeHighest := polyQuickDDCoeffsHighest + polyQuickDCoeffsHighest:

extrabound := hexa2ieee(["3F500000","00000000"]):
polyQuickDegreeExtra := 5:

bound[0] := 0:
b := nearest(lowestIntervalMax):
he := ieeehexa(b):
bound[1] := hexa2ieee([he[1],"00000000"]):
b := nearest(highestIntervalMin):
he := ieeehexa(b):
bound[intervals - 1] := hexa2ieee([he[1],"00000000"]):
bound[intervals] := 1:
linearWidth := (highestIntervalMin - lowestIntervalMax) / (intervals - 2):
scaleSum := 0:
scale[2] := 1.5:
scale[3] := 1.35:
scale[4] := 1.18:
scale[5] := 1.00:
scale[6] := 0.915:
scale[7] := 0.74:
scale[8] := 0.595:
scale[9] := 0.50:
for i from 2 to (intervals-1) do
	scaleSum := scaleSum + scale[i]:
od:
for i from 2 to (intervals-1) do
	scale[i] := scale[i] / scaleSum;
od:
current := lowestIntervalMax:
for i from 2 to (intervals-2) do
	b := nearest(current + (highestIntervalMin - lowestIntervalMax) * scale[i]):
	current := b:
	he := ieeehexa(b):
	bound[i] := hexa2ieee([he[1],"00000000"]):
od:
printf("Using %d intervals with bounds:\n",intervals):
for i from 0 to intervals do
	printf("bound[%d] = %1.30e\n",i,bound[i]):
od:
printf("Using an extra bound for truncating the quick poly to deg. %d in interval #0 for small args up to %f (2^(%f))\n",
	polyQuickDegreeExtra,extrabound,log[2](abs(extrabound)));

printf("Computing Gal's accurate table values for interval midpoints\n"):
for i from 2 to (intervals-1) do
	m := nearest((bound[i] + bound[i-1]) / 2):
	mhe := ieeehexa(m):
	printf("Interval %d: accurate table research start value %f (%s%s)\n",i,m,mhe[1],mhe[2]):
	midpointFloat[i] := galDoubleToDoubleDouble(nearest(m),arcsin,2^(-121),2^(15)):
od:


printf("Using the following floating point midpoints for intervals 2-%d:\n",intervals-2):
for i from 2 to (intervals-1) do
	mhe := ieeehexa(midpointFloat[i]):
	printf("midpointFloat[%d] = %f (%s%s)\n",i,midpointFloat[i],mhe[1],mhe[2]):
od:

printf("The reduced argument z is therefore bounded by:\n"):
printf("Interval 1: |z| < 2^(%f)\n",
	log[2](abs(bound[1]))):
for i from 2 to (intervals-1) do
	printf("Interval %d: |z| < 2^(%f)\n",i,
	log[2](max(abs(midpointFloat[i] - bound[i-1]),abs(midpointFloat[i] - bound[i])))):
od:
printf("Interval %d: |z| < 2^(%f)\n",intervals,
	log[2](abs(1 - bound[intervals-1]))):



printf("Using a %d degree polynomial for lowest interval (1) (accurate phase)\n",
	polyAccurateDegreeLowest):
printf("with %d triple-double, %d double-double and %d double coefficients\n",
	polyAccurateTDCoeffsLowest, polyAccurateDDCoeffsLowest, polyAccurateDCoeffsLowest):
printf("Using a %d degree polynomial for lowest interval (1) (quick phase)\n",
	polyQuickDegreeLowest):
printf("with %d double-double and %d double coefficients\n",
	polyQuickDDCoeffsLowest, polyQuickDCoeffsLowest):
printf("Using a %d degree polynomial for middle intervals (2-%d) (accurate phase)\n",
	polyAccurateDegreeMiddle,intervals-2):
printf("with %d triple-double, %d double-double and %d double coefficients\n",
	polyAccurateTDCoeffsMiddle, polyAccurateDDCoeffsMiddle, polyAccurateDCoeffsMiddle):
printf("Using a %d degree polynomial for middle intervals (2-%d) (quick phase)\n",
	polyQuickDegreeMiddle,intervals-2):
printf("with %d double-double and %d double coefficients\n",
	polyQuickDDCoeffsMiddle, polyQuickDCoeffsMiddle):
printf("Using a %d degree polynomial for highest interval (%d) (accurate phase)\n",
	polyAccurateDegreeHighest,intervals):
printf("with %d triple-double, %d double-double and %d double coefficients\n",
	polyAccurateTDCoeffsHighest, polyAccurateDDCoeffsHighest, polyAccurateDCoeffsHighest):
printf("Using a %d degree polynomial for highest interval (%d) (quick phase)\n",
	polyQuickDegreeHighest,intervals):
printf("with %d double-double and %d double coefficients\n",
	polyQuickDDCoeffsHighest, polyQuickDCoeffsHighest):

if true then 

printf("Computing polynomials for interval 1 ([%f;%f])\n",bound[0],bound[1]):

f := unapply(convert(series((arcsin(sqrt(x))/sqrt(x) - 1)/x,x=0,300),polynom),x):

p := unapply(eval(numapprox[chebyshev](f(x),x=(bound[0])^2..(bound[1])^2,2^(-127))),x):

polyAccuExact[1] := truncPoly(subs(X=x^2,p(X))*x^3+x,polyAccurateDegreeLowest):

polyAccurate[1] := poly_exact32(polyAccuExact[1],polyAccurateTDCoeffsLowest,polyAccurateDDCoeffsLowest):


epsAccurate[1] := numapprox[infnorm]((polyAccurate[1]/arcsin(x))-1, x=bound[0]..bound[1]):
epsAccurateSpecial := numapprox[infnorm]((polyAccurate[1]/arcsin(x))-1, x=bound[0]..evalf(sin(2^(-18)))):


polyQuickExact[1] := truncPoly(polyAccuExact[1],polyQuickDegreeLowest):
polyQuick[1] := poly_exact2(polyQuickExact[1],polyQuickDDCoeffsLowest):

epsQuick[1] := numapprox[infnorm]((polyQuick[1]/arcsin(x))-1, x=bound[0]..bound[1]):

polyQuickExtraExact := truncPoly(polyAccuExact[1],polyQuickDegreeExtra):
polyQuickExtra := poly_exact2(polyQuickExtraExact,polyQuickDegreeExtra):
epsQuickExtra := numapprox[infnorm]((polyQuickExtra/arcsin(x))-1, x=bound[0]..extrabound):

end if:

for i from 2 to (intervals-1) do

printf("Computing polynomials for interval %d ([%f;%f])\n",i,bound[i-1],bound[i]):
printf("Reduced argument z will be in interval [%1.8e;%1.8e] ([-2^%f,2^%f]),\n",
	bound[i-1]-midpointFloat[i],bound[i]-midpointFloat[i],
	log[2](abs(bound[i-1]-midpointFloat[i])),log[2](abs(bound[i]-midpointFloat[i]))):

f := unapply(arcsin(x+midpointFloat[i]),x):

fhelp := unapply(convert(series((f(x)-arcsin(midpointFloat[i]))/x,x=0,polyAccurateDegreeMiddle*3),polynom),x):

polyAccuExact[i] := numapprox[minimax](fhelp(x),
				x=bound[i-1]-midpointFloat[i]..bound[i]-midpointFloat[i],
				[polyAccurateDegreeMiddle,0],1,'err')*x + nearestDD(evalf(arcsin(midpointFloat[i]))):

polyAccurate[i] := poly_exact32(polyAccuExact[i],polyAccurateTDCoeffsMiddle,polyAccurateDDCoeffsMiddle):

epsAccurate[i] := numapprox[infnorm]((polyAccurate[i]/f(x))-1, 
					x=bound[i-1]-midpointFloat[i]..bound[i]-midpointFloat[i]):

polyQuickExact[i] := truncPoly(polyAccuExact[i],polyQuickDegreeMiddle):
polyQuick[i] := poly_exact2(polyQuickExact[i],polyQuickDDCoeffsMiddle):

epsQuick[i] := numapprox[infnorm]((polyQuick[i]/f(x))-1, 
			 	   x=bound[i-1]-midpointFloat[i]..bound[i]-midpointFloat[i]):


od:


if true then 

printf("Computing polynomials for interval %d ([%f;%f])\n",intervals,bound[intervals-1],bound[intervals]):


g := unapply(((arcsin(1 - x) - Pi/2)/sqrt(2*x)),x):
f := unapply(convert(series(((g(x)+1)/x),x=0,polyAccurateDegreeHighest*4),polynom),x):


polyAccuExact[intervals] := numapprox[minimax](f(x),x=(1-bound[intervals]+2^(-53))..(1-bound[intervals-1]),
						[polyAccurateDegreeHighest-1,0],1,'err')*x-1:

polyAccurate[intervals] := poly_exact32(polyAccuExact[intervals],polyAccurateTDCoeffsHighest,polyAccurateDDCoeffsHighest):


epsAccurate[intervals] := numapprox[infnorm](((unapply(polyAccurate[intervals],x)(x)/g(x))-1), 
					x=(1-bound[intervals]+2^(-53))..(1-bound[intervals-1])):


polyQuickExact[intervals] := truncPoly(polyAccuExact[intervals],polyQuickDegreeHighest):
polyQuick[intervals] := poly_exact2(polyQuickExact[intervals],polyQuickDDCoeffsHighest):

epsQuick[intervals] := numapprox[infnorm](((unapply(polyQuick[intervals],x)(x)/g(x))-1), 
					x=(1-bound[intervals]+2^(-53))..(1-bound[intervals-1])):

printf("Checking if the polynomial for interval %d is exactly -1 in z = %f...\n",intervals,1-bound[intervals]);

if (unapply(polyAccurate[intervals],x)(1-bound[intervals]) = -1) then
	printf("  Check passed!\n"):
else
	printf("  Check failed!\n"):
end if:

end if:



for i from 1 to intervals do 
	printf("Relative error for accurate phase polynomial in interval %d ([%f;%f]) is 2^(%f)\n",
       		i,bound[i-1],bound[i],log[2](abs(epsAccurate[i]))):
	printf("Relative error for quick phase polynomial in interval %d ([%f;%f]) is 2^(%f)\n",
       		i,bound[i-1],bound[i],log[2](abs(epsQuick[i]))):
od:
printf("Relative error for accurate phase polynomial #1 in special interval [0;sin(2^(-18))]) is 2^(%f)\n",
       log[2](abs(epsAccurateSpecial))):

printf("Relative error for quick phase extra case truncated polynomial #1 in special interval [0;%f]) is 2^(%f)\n",
       extrabound,log[2](abs(epsQuickExtra))):


(PiHalfH, PiHalfM, PiHalfL) := hi_mi_lo(evalf(Pi/2)):

epsPiDD := evalf(((PiHalfH + PiHalfM) - Pi/2)/(Pi/2)):
epsPiTD := evalf(((PiHalfH + PiHalfM + PiHalfL) - Pi/2) / (Pi/2)):

printf("Relative error for storing Pi/2 as a double-double is 2^(%f)\n",evalf(log[2](abs(epsPiDD)))):
printf("Relative error for storing Pi/2 as a triple-double is 2^(%f)\n",evalf(log[2](abs(epsPiTD)))):

# Ce qui suit est pifometrique et doit etre prouve en Gappa ensuite

arithmeticalErrorQuick[1] := 2^(-61):
arithmeticalErrorQuick[2] := 2^(-61):
arithmeticalErrorQuick[3] := 2^(-61):
arithmeticalErrorQuick[4] := 2^(-61):
arithmeticalErrorQuick[5] := 2^(-61):
arithmeticalErrorQuick[6] := 2^(-61):
arithmeticalErrorQuick[7] := 2^(-61):
arithmeticalErrorQuick[8] := 2^(-61):
arithmeticalErrorQuick[9] := 2^(-61):
arithmeticalErrorQuick[10] := 2^(-68):

arithmeticalErrorQuickExtra := 2^(-80):

for i from 1 to intervals do
	estimatedOverallEpsQuick[i] := abs(epsQuick[i]) + 
				       abs(arithmeticalErrorQuick[i]) + 
                                       abs(epsQuick[i]) * abs(arithmeticalErrorQuick[i]):
	printf("Relative quick phase overall error bound to show in Gappa for interval %d ([%f;%f]) is 2^(%f)\n",
		i,bound[i-1],bound[i],log[2](abs(estimatedOverallEpsQuick[i]))):
od:

estimatedOverallEpsQuickExtra := abs(epsQuickExtra) + 
				 abs(arithmeticalErrorQuickExtra) + 
                                 abs(epsQuickExtra) * abs(arithmeticalErrorQuickExtra):
printf("Relative quick phase overall error bound to show for extra truncted poly in interval [0;%f]) is 2^(%f)\n",
	extrabound,log[2](abs(estimatedOverallEpsQuickExtra))):


printf("Write tables...\n"):

filename:="TEMPASIN/asin-td.h":
fd:=fopen(filename, WRITE, TEXT):

fprintf(fd, "#include \"crlibm.h\"\n#include \"crlibm_private.h\"\n"):

fprintf(fd, "\n/* File generated by maple/asin-td.mpl */\n\n"):

printf("Write table with bounds\n"):

fprintf(fd, "/* High order words of interval bounds (low order word is 0) */\n"):
for i from 1 to (intervals-1) do 
	heb := ieeehexa(bound[i]):
	fprintf(fd, "\#define BOUND%d 0x%s\n",i,heb[1]):
od:

heb := ieeehexa(extrabound):
fprintf(fd, "\#define EXTRABOUND 0x%s\n",heb[1]):

printf("Write additional constants\n"):

fprintf(fd, "\n\n/* Pi/2 as a triple-double*/\n"):
fprintf(fd, "\#define PIHALFH %1.50e\n",PiHalfH):
fprintf(fd, "\#define PIHALFM %1.50e\n",PiHalfM):
fprintf(fd, "\#define PIHALFL %1.50e\n",PiHalfL):

printf("Write table with midpoints and polynomial coefficients\n"):
k := 0:
for i from 0 to polyAccurateDegreeLowest do
	(hi,mi,lo) := hi_mi_lo(coeff(polyAccurate[1],x,i)):
	if ((abs(hi) = 1.0) and (mi = 0) and (lo = 0)) then 
		printf(
		"Coefficient %d of interval 1 polynomial is exactly %f and will not be stored in the table\n",i,hi): 
	else 
	g := 0;
	if (hi <> 0) then
		if (i <= polyQuickDegreeLowest) then
			tbl[k] := sprintf("%1.50e, \t/* %d, accPolyLowC%dh, quickPolyLowC%dh */",hi,k,i,i):
		else
			tbl[k] := sprintf("%1.50e, \t/* %d, accPolyLowC%dh */",hi,k,i):
		end if:
		k := k + 1:
	end if:
	if (mi <> 0) then
		if (i <= polyQuickDDCoeffsLowest) then
			tbl[k] := sprintf("%1.50e, \t/* %d, accPolyLowC%dm, quickPolyLowC%dl */",mi,k,i,i):
		else
			tbl[k] := sprintf("%1.50e, \t/* %d, accPolyLowC%dm */",mi,k,i):
		end if:
		k := k + 1:
	end if:
	if (lo <> 0) then
		tbl[k] := sprintf("%1.50e, \t/* %d, accPolyLowC%dl */",lo,k,i):
		k := k + 1:
	end if:
	end if:
od:
tbl[k] := sprintf("%1.50e, \t/* %d, RN rounding constant quick poly low*/",
		compute_rn_constant(estimatedOverallEpsQuick[1]),k):
k := k + 1:
tbl[k] := sprintf("%1.50e, \t/* %d, RD rounding constant quick poly low*/",
		estimatedOverallEpsQuick[1],k):
k := k + 1:
printf("Table for interval 1 written\n"):
for l from 2 to (intervals-1) do
	tbl[k] := sprintf("%1.50e, \t/* %d, midpoint in interval %d*/",midpointFloat[l],k,l):
	tblidx[l] := k;
	k := k + 1;
	for i from 0 to polyAccurateDegreeMiddle do
		(hi,mi,lo) := hi_mi_lo(coeff(polyAccurate[l],x,i)):
		if ((abs(hi) = 1.0) and (mi = 0) and (lo = 0)) then 
			printf(
		"Coefficient %d of interval %d polynomial is exactly %f and will not be stored in the table\n",i,l,hi): 
		else 
		g := 0;
		if (hi <> 0) then
			if (i <= polyQuickDegreeMiddle) then
				tbl[k] := sprintf(
					"%1.50e, \t/* %d, accPolyMid%dC%dh, quickPolyMid%dC%dh */",hi,k,l,i,l,i):
			else
				tbl[k] := sprintf("%1.50e, \t/* %d, accPolyMid%dC%dh */",hi,k,l,i):
			end if:
			k := k + 1:
		end if:
		if (mi <> 0) then
			if (i <= polyQuickDDCoeffsMiddle) then
				tbl[k] := sprintf(
					"%1.50e, \t/* %d, accPolyMid%dC%dm, quickPolyMid%dC%dl */",mi,k,l,i,l,i):
			else
				tbl[k] := sprintf("%1.50e, \t/* %d, accPolyMid%dC%dm */",mi,k,l,i):
			end if:
			k := k + 1:
		end if:
		if (lo <> 0) then
			tbl[k] := sprintf("%1.50e, \t/* %d, accPolyMid%dC%dl */",lo,k,l,i):
			k := k + 1:
		end if:
		end if:
	od:
	tbl[k] := sprintf("%1.50e, \t/* %d, RN rounding constant quick poly middle %d*/",
			compute_rn_constant(estimatedOverallEpsQuick[l]),k,l):
	k := k + 1:
	tbl[k] := sprintf("%1.50e, \t/* %d, RD rounding constant quick poly middle %d*/",
			estimatedOverallEpsQuick[l],k,l):
	k := k + 1:
	printf("Table for interval %d written\n",l):
od:
tblidx[intervals] := k:
for i from 0 to polyAccurateDegreeHighest do
	(hi,mi,lo) := hi_mi_lo(coeff(polyAccurate[intervals],x,i)):
	if ((abs(hi) = 1.0) and (mi = 0) and (lo = 0)) then 
		printf(
	"Coefficient %d of interval %d polynomial is exactly %f and will not be stored in the table\n",i,intervals,hi): 
	else 
	g := 0;	
	if (hi <> 0) then
		if (i <= polyQuickDegreeHighest) then
			tbl[k] := sprintf("%1.50e, \t/* %d, accPolyHighC%dh, quickPolyHighC%dh */",hi,k,i,i):
		else
			tbl[k] := sprintf("%1.50e, \t/* %d, accPolyHighC%dh */",hi,k,i):
		end if:
		k := k + 1:
	end if:
	if (mi <> 0) then
		if (i <= polyQuickDDCoeffsHighest) then
			tbl[k] := sprintf("%1.50e, \t/* %d, accPolyHighC%dm, quickPolyHighC%dl */",mi,k,i,i):
		else
			tbl[k] := sprintf("%1.50e, \t/* %d, accPolyHighC%dm */",mi,k,i):
		end if:
		k := k + 1:
	end if:
	if (lo <> 0) then
		tbl[k] := sprintf("%1.50e, \t/* %d, accPolyHighC%dl */",lo,k,i):
		k := k + 1:
	end if:
	end if:
od:
tbl[k] := sprintf("%1.50e, \t/* %d, RN rounding constant quick poly high*/",
		compute_rn_constant(estimatedOverallEpsQuick[intervals]),k):
k := k + 1:
tbl[k] := sprintf("%1.50e, \t/* %d, RD rounding constant quick poly high*/",
 		estimatedOverallEpsQuick[intervals],k):
k := k + 1:
printf("Table for interval %d written\n",intervals):
tbllen := k:
printf("The whole table has %d entries, so uses %d bytes of memory\n",tbllen,tbllen*8):
fprintf(fd,"\n\n/* Indices to the following table */\n"):
for i from 2 to intervals do
	fprintf(fd,"\#define TBLIDX%d %d\n",i,tblidx[i]):
od:
fprintf(fd, "\n\n/* Table with midpoints and polynomial coefficients */\n"):
fprintf(fd, "static const double tbl[%d] = {\n",tbllen):
for i from 0 to (tbllen - 1) do
	fprintf(fd, "%s\n",tbl[i]):
od:
fprintf(fd, "};\n\n"):

fclose(fd):

printf("----DONE---\n"):


