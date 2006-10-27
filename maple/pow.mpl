#######################################################################
# This file is part of the crlibm library, and is distributed under
# the  LGPL.
# To use:
# restart; read "exp-td.mpl";
Digits := 120:

interface(quiet=true):

read "common-procedures.mpl":
read "triple-double.mpl":

mkdir("TEMPPOW"):

# Some constants for special cases tests 

two52 := 2^(52):
two53 := 2^(53):
twoM53 := 2^(-53):
twoM54 := 2^(-54):
twoM1021 := 2^(-1021):
two1021 := 2^(1021):
twoM1000 := 2^(-1000):
two1000 := 2^(1000):
two999 := 2^(999):
two11 := 2^(11):
two64 := 2^(64):
twoM64 := 2^(-64):

scale := 2^(12):
rescale := 2^(-12):
shiftConst := 2^(52) + 2^(51):

largest := 2^(1023) * ((2^(53) - 1) / 2^(52)):
smallest := 2^(-1023) * 1 * 2^(-51):

# Logarithm log2_12 for power

L := 7: # number of bits used to address the table

MAXINDEX    := round(2^L * (sqrt(2)-1)):

for i from 0 to MAXINDEX-1 do
    center[i] := 1 + i*2^(-L): # center[i] in [1, 2[
    t :=  evalf(1/center[i]):
    r[i] := round(t * 2**(floor(-log[2](abs(t))) + 23)) / 2**(floor(-log[2](abs(t))) + 23):
    (logih[i], logil[i], logill[i]) := hi_mi_lo(evalf(-log[2](r[i]))):
od:
for i from MAXINDEX to 2^L do
    # y has been divided by two, center[i] in [0.5, 1[
    center[i]:=(1 + i*2^(-L)) / 2:
    t :=  evalf(1/center[i]):
    r[i] := round(t * 2**(floor(-log[2](abs(t))) + 23)) / 2**(floor(-log[2](abs(t))) + 23):
    (logih[i], logil[i], logill[i]) := hi_mi_lo(evalf(-log[2](r[i]))):
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
    if zmax[i] > zmaxmax then zmaxmax := zmax[i]: fi:
    if zmin[i] < zminmin then zminmin := zmin[i]: fi:
od:
printf("zminmin = -2^(%2f)   zmaxmax = 2^(%2f)\n", log2(-zminmin), log2(zmaxmax) ):

pLogExact := unapply(-6497523352912697/36028797018963968*X^8+464105351149111/2251799813685248*X^7-4331547231022885/18014398509481984*X^6+324866042375467/1125899906842624*X^5-6497320848556805/18014398509481984*X^4+8460053188225/17592186044416*X^3-58522663504933901518104329421789/81129638414606681695789005144064*X^2+29261331752466950759073917813481/20282409603651670423947251286016*X,X):

pLog := poly_exact2(pLogExact(x),3):

epsLog := numapprox[infnorm](pLog/log[2](1+x)-1,x=zminmin..zmaxmax):

printf("Relative error of pLog w.r.t log2(1+x) is 2^(%f)\n",log[2](abs(epsLog))):


# Exponential exp2_12 for power 


for i from 0 to 2^6 - 1 do
	twoPowerIndex1hi[i], twoPowerIndex1lo[i], twoPowerIndex1lolo[i] := hi_mi_lo(evalf(2^(i/(2^6)))):
	twoPowerIndex2hi[i], twoPowerIndex2lo[i], twoPowerIndex2lolo[i] := hi_mi_lo(evalf(2^(i/(2^(12))))):
od: 

pExpExact := unapply(2772236920359585/288230376151711744*X^4+1999746284130149/36028797018963968*X^3+8655072057804175/36028797018963968*X^2+6243314768165359/9007199254740992*X+1,X):

pExp := poly_exact(pExpExact(x)):

epsExp := numapprox[infnorm]((pExp/2^x)-1,x=-2^(-12)..2^(-12)):

printf("Relative error of pExp w.r.t 2^x is 2^(%f)\n",log[2](abs(epsExp))):

log2 := nearest(log(2)):


# Exponential exp2_33 for power

# Polynomial for approximating 2^x - 1 in x=-2^(-42)..2^(-42)

pExpXM := unapply(poly_exact2(x * ((6243314768165359 * 2^(-53) + 29397410857115 * 2^(-100)) + (x * 8655072057804175 * 2^(-55))),2),x):

epsExpXM := numapprox[infnorm](pExpXM(x)/(2^x-1)-1,x=-2^(-42)..2^(-42)):

printf("The relative error of pExpXM w.r.t. 2^x - 1 is 2^(%f)\n",log[2](abs(epsExpXM))):


# Polynomial for approximating 2^x - 1 in x=-2^(-95)..2^(-95)

pExpXL := unapply(poly_exact(x * 6243314768165359 * 2^(-53)), x):

epsExpXL := numapprox[infnorm](pExpXL(x)/(2^x-1)-1,x=-2^(-95)..2^(-95)):

printf("The relative error of pExpXL w.r.t. 2^x - 1 is 2^(%f)\n",log[2](abs(epsExpXL))):

# Polynomial for approximating 2^x in x=-2^(-12)..2^(-12)

pExpXH := unapply(poly_exact32(x * (506517869649829535567849302923399275789356375957 * 2^(-159) + (x * (38978979294391673005692521213079 * 2^(-107) + (x * (36024226132016099441525232746301 * 2^(-109) + (x * (3121261346907607936312652866425 * 2^(-108) + (x * (55385433661433492776134419224183 * 2^(-115) + (x * (5682899659966205 * 2^(-65) + (x * (4501812434047971 * 2^(-68) + (x * 6240991224781291 * 2^(-72))))))))))))))),2,4),x):

epsExpXH := numapprox[infnorm](pExpXH(x)/(2^x-1)-1,x=-2^(-12)..2^(-12)):

printf("The relative error of pExpXH w.r.t. 2^x - 1 is 2^(%f)\n",log[2](abs(epsExpXH))):


# Logarithm log2_13 for power

pLog13 := unapply(poly_exact32(x * (65890661388387311068680317907364672336343222485 * 2^(-155) + (x * ((-263562645553549244274721271629458689351564598207 * 2^(-158)) + (x * (39015109003289267678766993386435 * 2^(-106) + (x * ((-29261331752466950759075245039823 * 2^(-106)) + (x * (11704532700986780303630098000775 * 2^(-105) + (x * ((-78030218006578535357533995772145 * 2^(-108)) + (x * (66883044005638744592219245028355 * 2^(-108) + (x * ((-14630665876233475378839909031845 * 2^(-106)) + (x * (721924538728533 * 2^(-52) + (x * ((-324866042427843 * 2^(-51)) + (x * (4725324253587977 * 2^(-55) + (x * ((-8663094456247059 * 2^(-56)) + (x * (7996693900181249 * 2^(-56) + (x * ((-7425880934468497 * 2^(-56)) + (x * 7091529758988931 * 2^(-56))))))))))))))))))))))))))))),3,6),x):

epsLog13 := numapprox[infnorm](pLog13(x)/(log[2](1+x))-1,x=-2^(-8)..2^(-8)):

printf("The relative error of pLog13 w.r.t. log2(1+x) is 2^(%f)\n",log[2](abs(epsLog13))):


# exp2_30bits for exactness test 

coeff_0 :=  9.99999999947486895024439945700578391551971435546875000000000000000000000000000000e-01:
coeff_1 :=  6.93147180274189311788290979166049510240554809570312500000000000000000000000000000e-01:
coeff_2 :=  2.40226513201275415632096610352164134383201599121093750000000000000000000000000000e-01:
coeff_3 :=  5.55041194035996443556513213479775004088878631591796875000000000000000000000000000e-02:
coeff_4 :=  9.61801251055323207228564541537707555107772350311279296875000000000000000000000000e-03:
coeff_5 :=  1.33325640280455024258565721595459763193503022193908691406250000000000000000000000e-03:
coeff_6 :=  1.54736006782907911617439000728779774362919852137565612792968750000000000000000000e-04:
coeff_7 :=  1.55294506644329091183710789270122631933190859854221343994140625000000000000000000e-05:

pExp2Exact := unapply(coeff_0 + X * (coeff_1 + X * (coeff_2 + X * (coeff_3 + X * (coeff_4 + X * (coeff_5 + X * (coeff_6 + X * coeff_7)))))),X):

pExp2 := poly_exact(pExp2Exact(x)):

epsExp2 := numapprox[infnorm]((pExp2/2^x)-1,x=-0.5..0.5):

printf("Relative error of pExp2 w.r.t 2^x is 2^(%f)\n",log[2](abs(epsExp2))):

# Overall accuracy estimate

# C'est pifometrique
epsOverall := 2^(-62):

bi := floor(-log[2](abs(epsOverall))):
approxBoundFactor := 2^(-(bi - 54)):

epsOverallAccurate := 2^(-120):

bi2 := floor(-log[2](abs(epsOverallAccurate))) - 1:
approxBoundFactorAccurate := 2^(-(bi2 - 54)):


# Print out of the .h file 

filename:="TEMPPOW/pow.h":
fd:=fopen(filename, WRITE, TEXT):

fprintf(fd, "#include \"crlibm.h\"\n#include \"crlibm_private.h\"\n"):

fprintf(fd, "\n/*File generated by maple/pow.mpl*/\n"):

fprintf(fd, "\#define APPROXBOUNDFACTOR %1.50e\n", approxBoundFactor):   
fprintf(fd, "\#define APPROXBOUNDFACTORACCURATE %1.50e\n", approxBoundFactorAccurate):   
fprintf(fd, "\#define TWO52 %1.50e\n", two52):   
fprintf(fd, "\#define TWO53 %1.50e\n", two53):   
fprintf(fd, "\#define TWO11 %1.50e\n", two11):   
fprintf(fd, "\#define TWOM53 %1.50e\n", twoM53):   
fprintf(fd, "\#define TWOM54 %1.50e\n", twoM54):   
fprintf(fd, "\#define TWOM1021 %1.50e\n", twoM1021):   
fprintf(fd, "\#define TWO1021 %1.50e\n", two1021):   
fprintf(fd, "\#define TWOM1000 %1.50e\n", twoM1000):   
fprintf(fd, "\#define TWO999 %1.50e\n", two999):   
fprintf(fd, "\#define TWO1000 %1.50e\n\n", two1000):   
fprintf(fd, "\#define TWO64 %1.50e\n", two64):   
fprintf(fd, "\#define TWOM64 %1.50e\n", twoM64):   
fprintf(fd, "\#define SCALE %1.50e\n", scale):   
fprintf(fd, "\#define RESCALE %1.50e\n", rescale):   
fprintf(fd, "\#define SHIFTCONSTANT %1.50e\n", shiftConst):   
fprintf(fd, "\#define LARGEST %1.50e\n",largest):
fprintf(fd, "\#define SMALLEST %1.50e\n\n",smallest):


(log2coeff1dh,log2coeff1dl) := hi_lo(coeff(pLog,x,1)):
(log2coeff2dh,log2coeff2dl) := hi_lo(coeff(pLog,x,2)):
fprintf(fd, "\#define log2coeff1h %1.50e\n",log2coeff1dh):	
fprintf(fd, "\#define log2coeff1l %1.50e\n",log2coeff1dl):	
fprintf(fd, "\#define log2coeff2h %1.50e\n",log2coeff2dh):	
fprintf(fd, "\#define log2coeff2l %1.50e\n",log2coeff2dl):	
for i from 3 to 8 do
	fprintf(fd, "\#define log2coeff%d %1.50e\n",i,coeff(pLog,x,i)):
od:
fprintf(fd,"\n"):

for i from 1 to 4 do
	fprintf(fd, "\#define exp2coeff%d %1.50e\n",i,coeff(pExp,x,i)):
od:
fprintf(fd,"\n"):

for i from 0 to 7 do 
	fprintf(fd, "\#define exp2InaccuCoeff%d %1.50e\n",i,coeff(pExp2,x,i)):
od:
fprintf(fd,"\n"):

(exp2XMcoeff1dh, exp2XMcoeff1dl) := hi_lo(coeff(pExpXM(x),x,1)):
exp2XMcoeff2dh := nearest(coeff(pExpXM(x),x,2)):
fprintf(fd, "\#define exp2XMcoeff1h %1.50e\n",exp2XMcoeff1dh):	
fprintf(fd, "\#define exp2XMcoeff1l %1.50e\n",exp2XMcoeff1dl):	
fprintf(fd, "\#define exp2XMcoeff2h %1.50e\n\n",exp2XMcoeff2dh):	

exp2XLcoeff1dh := nearest(coeff(pExpXL(x),x,1)):
fprintf(fd, "\#define exp2XLcoeff1h %1.50e\n\n",exp2XLcoeff1dh):	

(exp2XHcoeff1dh,exp2XHcoeff1dm,exp2XHcoeff1dl) := hi_mi_lo(coeff(pExpXH(x),x,1)):
fprintf(fd, "\#define exp2XHcoeff1h %1.50e\n",exp2XHcoeff1dh):	
fprintf(fd, "\#define exp2XHcoeff1m %1.50e\n",exp2XHcoeff1dm):	
fprintf(fd, "\#define exp2XHcoeff1l %1.50e\n",exp2XHcoeff1dl):	

for i from 2 to 5 do
	(exp2XHcoeffidh,exp2XHcoeffidl) := hi_lo(coeff(pExpXH(x),x,i)):
	fprintf(fd, "\#define exp2XHcoeff%dh %1.50e\n",i,exp2XHcoeffidh):	
	fprintf(fd, "\#define exp2XHcoeff%dm %1.50e\n",i,exp2XHcoeffidl):	
od:

for i from 6 to 8 do 
	fprintf(fd, "\#define exp2XHcoeff%dh %1.50e\n",i,nearest(coeff(pExpXH(x),x,i))):		
od:

fprintf(fd,"\n"):

(log213coeff1dh,log213coeff1dm,log213coeff1dl) := hi_mi_lo(coeff(pLog13(x),x,1)):
fprintf(fd, "\#define log213coeff1h %1.50e\n",log213coeff1dh):	
fprintf(fd, "\#define log213coeff1m %1.50e\n",log213coeff1dm):	
fprintf(fd, "\#define log213coeff1l %1.50e\n",log213coeff1dl):	

(log213coeff2dh,log213coeff2dm,log213coeff2dl) := hi_mi_lo(coeff(pLog13(x),x,2)):
fprintf(fd, "\#define log213coeff2h %1.50e\n",log213coeff2dh):	
fprintf(fd, "\#define log213coeff2m %1.50e\n",log213coeff2dm):	
fprintf(fd, "\#define log213coeff2l %1.50e\n",log213coeff2dl):	

for i from 3 to 8 do
	(log213coeffidh, log213coeffidl) := hi_lo(coeff(pLog13(x),x,i)):
	fprintf(fd, "\#define log213coeff%dh %1.50e\n",i,log213coeffidh):	
	fprintf(fd, "\#define log213coeff%dm %1.50e\n",i,log213coeffidl):	
od:

for i from 9 to 15 do 
	fprintf(fd, "\#define log213coeff%dh %1.50e\n",i,nearest(coeff(pLog13(x),x,i))):		
od:


fprintf(fd, "\n\#define LOG2 %1.50e\n\n", log2):   

fprintf(fd, "typedef struct rri_tag {float ri; double logih; double logil; double logill;} rri;  \n"):
fprintf(fd, "static const rri argredtable[%d] = {\n", 2^L):
for i from 0 to 2^L-1 do
      fprintf(fd, "  { \n"):
      fprintf(fd, "    %1.50e,   /* r[%d] */ \n", r[i], i):
      fprintf(fd, "    %1.50e, /* logih[%d] */ \n", logih[i], i):
      fprintf(fd, "    %1.50e, /* logil[%d] */ \n", logil[i], i):
      fprintf(fd, "    %1.50e, /* logill[%d] */ \n", logill[i], i):
      fprintf(fd, "  } "):
      if(i<2^L-1) then  fprintf(fd, ", \n"): fi
od:
fprintf(fd, "}; \n \n"):

fprintf(fd, "typedef struct tPi_t_tag {double hi; double lo; double lolo;} tPi_t;  \n"):
fprintf(fd, "static const tPi_t twoPowerIndex1[%d] = {\n", 2^(6)):
for i from 0 to 2^(6)-1 do
      fprintf(fd, "  { \n"):      
      fprintf(fd, "    %1.50e, /* twoPowerIndex1hi[%d] */ \n", twoPowerIndex1hi[i], i):
      fprintf(fd, "    %1.50e, /* twoPowerIndex1lo[%d] */ \n", twoPowerIndex1lo[i], i):
      fprintf(fd, "    %1.50e, /* twoPowerIndex1lolo[%d] */ \n", twoPowerIndex1lolo[i], i):
      fprintf(fd, "  } "):
      if(i<2^(6)-1) then  fprintf(fd, ", \n"): fi
od:
fprintf(fd, "}; \n \n"):
fprintf(fd, "static const tPi_t twoPowerIndex2[%d] = {\n", 2^(6)):
for i from 0 to 2^(6)-1 do
      fprintf(fd, "  { \n"):      
      fprintf(fd, "    %1.50e, /* twoPowerIndex2hi[%d] */ \n", twoPowerIndex2hi[i], i):
      fprintf(fd, "    %1.50e, /* twoPowerIndex2lo[%d] */ \n", twoPowerIndex2lo[i], i):
      fprintf(fd, "    %1.50e, /* twoPowerIndex2lolo[%d] */ \n", twoPowerIndex2lolo[i], i):
      fprintf(fd, "  } "):
      if(i<2^(6)-1) then  fprintf(fd, ", \n"): fi
od:
fprintf(fd, "}; \n \n"):





fclose(fd):

printf("--------- DONE -----------\n");