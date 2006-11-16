#######################################################################
# This file is part of the crlibm library, and is distributed under
# the  LGPL.
# To use:
# restart; read "log-de.mpl";
#TODO output scripts for the Gappa proof out for Paterson/Stockmayer
Digits := 100:

interface(quiet=true):

read "common-procedures.mpl":
read "double-extended.mpl":
mkdir("TEMPLOG"):


log2h,log2l := hiloExt(log(2)):


L := 7: # number of bits used to address the table

MAXINDEX    := round(2^L * (sqrt(2)-1)):

for i from 0 to MAXINDEX-1 do
    center[i] := 1 + i*2^(-L): # center[i] in [1, 2[
    # We want it to fit on 11 bits of mantissa
    r[i] :=  round(evalf(  (1/center[i]) * 2^(11)) ) / 2^(11) ;

od:
for i from MAXINDEX to 2^L do
    # y has been divided by two, center[i] in [0.5, 1[
    center[i]:=(1 + i*2^(-L)) / 2:
    # We want it to fit on 11 bits of mantissa,
    r[i] :=  round(evalf(  (1/center[i]) * 2^(10)) ) / 2^(10) ;
od:

# Note that we go up to 2^L although the case 2^L is wrapped to zero
# in the C code. It could be important for zmax (but it turns out not).

for i from 0 to 2^L do
    (logirh[i], logirl[i]) := hiloExt(-log(r[i])):
od:

#Computation of ZMax.
for i from 0 to MAXINDEX-1 do
    t_x := center[i] + 2^(-L-1) :
    zmax[i] := (t_x*r[i]-1) :
    t_x := center[i] - 2^(-L-1) :
    zmin[i] := (t_x*r[i]-1) :
    zabsmax[i] := max(abs(zmin[i]), abs(zmax[i])):
od:
for i from MAXINDEX to 2^L do
    t_x := center[i] + 2^(-L-2) :
    zmax[i] := (t_x*r[i]-1) :
    t_x := center[i] - 2^(-L-2) :
    zmin[i] := (t_x*r[i]-1) :
    zabsmax[i] := max(abs(zmin[i]), abs(zmax[i])):
od:

zmaxmax:=0:
zminmin:=0:
for i from 0 to 2^L do
    if zmax[i] > zmaxmax then zmaxmax := zmax[i]: fi:
    if zmin[i] < zminmin then zminmin := zmin[i]: fi:
od:
printf("zminmin = -2^(%2f)   zmaxmax = 2^(%2f)\n", log2(-zminmin), log2(zmaxmax) ) :


PolyDegreeQuick:=7:

#Keep -zmaxmax..zmaxmax to keep c1=1, which is useful in the proof

if(1+1=3) then # The polynomial used to be computed here,
polyQuick:= polyExact2Ext(x  * numapprox[minimax](  log(1+x)/x,  x=-zmaxmax..zmaxmax,  [PolyDegreeQuick-1,0], 1 ,  'deltaApprox'), 0):

epsilonApproxQuick := numapprox[infnorm]( 1-polyQuick/log(1+x), x=zminmin..zmaxmax):
deltaApproxQuick := numapprox[infnorm]( polyQuick-log(1+x), x=zminmin..zmaxmax):
else
# magical polynomial given by Sylvain's tool
polyQuick := x * (1 + (x * ((-0.50000000000000000000000000000000000000000000000000000000000000e0) + (x * (0.33333333333333342585191871876304503530263900756835937500000000e0 + (x * ((-0.24999999998423708125194764306797878816723823547363281250000000e0) + (x * (0.19999999996748479835773082413652446120977401733398437500000000e0 + (x * ((-0.16666957260954737285452154083031928166747093200683593750000000e0) + (x * 0.14286056338555042088955815415829420089721679687500000000000000e0)))))))))))):

# validated infinite norms given by Christoph's tool
epsilonApproxQuick := 2^(-64.119667587221):
deltaApproxQuick := 2^(-72.230387592107):
fi:
printf(" rel approximation error for the quick phase is 2^(%2f)\n", log2(epsilonApproxQuick) ) :
printf(" abs approximation error for the quick phase is 2^(%2f)\n", log2(deltaApproxQuick) ) :

# For the Paterson/Stockmayer method
polyQuick0:= x  * numapprox[minimax](  log(1+x)/x,  x=-zmaxmax..zmaxmax,  [PolyDegreeQuick-1,0], 1 ,  'deltaApprox'):

## Création des nouveaux coefficients
polyQuick0:= x  * numapprox[minimax](  log(1+x)/x,  x=-zmaxmax..zmaxmax,  [PolyDegreeQuick-1,0], 1 ,  'deltaApprox'):

## Création des nouveaux coefficients
a7 := convert(coeff(polyQuick0,x,7),rational):
a6 := convert(coeff(polyQuick0,x,6)/a7,rational):
a5 := convert(coeff(polyQuick0,x,5)/a7,rational):
a4 := convert(coeff(polyQuick0,x,4)/a7,rational):

a3 := convert(coeff(polyQuick0,x,3),rational):
a2 := convert(coeff(polyQuick0,x,2),rational):
a1 := convert(coeff(polyQuick0,x,1),rational):
a0 := convert(coeff(polyQuick0,x,0),rational):

alpha0 := a5 - 1:
alpha1 := a6:
alpha2 := a4 - alpha0 * a6:

##
a7_hi,a7_lo         := hiloExt(a7):
a3_hi,a3_lo         := hiloExt(a3):
a2_hi,a2_lo         := hiloExt(a2):
a0_hi,a0_lo         := hiloExt(a0):

c_hi,c_lo         := hiloExt(alpha0):
alpha_hi,alpha_lo := hiloExt(alpha1):
beta_hi,beta_lo   := hiloExt(alpha2):

#coef_c := c_hi + c_lo:
#coef_alpha := alpha_hi + alpha_lo:
#coef_beta := beta_hi + beta_lo:
#coef_c7 := a7_hi + a7_lo:
#coef_c3 := a3_hi + a3_lo:
#coef_c2 := a2_hi + a2_lo:
#coef_c0 := a0_hi + a0_lo:

coef_c := c_hi:
coef_alpha := alpha_hi:
coef_beta := beta_hi:
coef_c7 := a7_hi:
coef_c3 := a3_hi:
coef_c2 := a2_hi:
coef_c0 := a0_hi:

polyQuickPat := ((x^2 + coef_c)*(x + coef_alpha) + (x + coef_beta)) * (coef_c7*x^4) + ((coef_c3 * x + coef_c2)*x^2 + (x)):

epsilonApproxQuickPat := numapprox[infnorm]( 1-polyQuickPat/log(1+x), x=zminmin..zmaxmax):
printf("   approximation rel error for Paterson/Stockmayer in the quick phase is 2^(%2f)\n", log2(epsilonApproxQuickPat) ) :
deltaApproxQuickPat := numapprox[infnorm]( polyQuickPat-log(1+x), x=zminmin..zmaxmax):
printf("   approximation abs error for Paterson/Stockmayer in the quick phase is 2^(%2f)\n", log2(deltaApproxQuickPat) ) :



PolyDegreeAccurate:=14:

printf("Computing the polynomial for accurate phase (this may take some time...)\n"):
pe:= x  * numapprox[minimax](  log(1+x)/x,  x=-zmaxmax..zmaxmax,  [PolyDegreeAccurate-1,0], 1 ,  'deltaApprox'):


MaxDegreeDDE:=8:  #

polyAccurate := polyExact2Ext(pe, MaxDegreeDDE):
deltaApproxAccurate := numapprox[infnorm](polyAccurate-log(1+x), x=-zmaxmax..zmaxmax):
epsilonApproxAccurate := numapprox[infnorm]( 1-polyAccurate/log(1+x), x=-zmaxmax..zmaxmax):
printf("   approximation error for the accurate phase is 2^(%2f)\n", log2(epsilonApproxAccurate) ) :




filename:="TEMPLOG/log-de.h":
fd:=fopen(filename, WRITE, TEXT):

fprintf(fd, "/*File generated by maple/log-de.mpl*/\n\n"):
# some constants
fprintf(fd, "#define L        %d\n", L):
fprintf(fd, "#define MAXINDEX %d\n", MAXINDEX):
fprintf(fd, "#define INDEXMASK %d\n", 2^L-1):
fprintf(fd, "static const double two64 = %1.25e ;\n", evalf(2^64)):

fprintf(fd, "#if 1\n#define ESTRIN\n#else\n#define PATERSON\n#endif\n"):

fprintf(fd, "#if defined(CRLIBM_TYPECPU_X86) || defined(CRLIBM_TYPECPU_AMD64)\n\n"):

fprintf(fd, "#ifdef PATERSON\n"):
  #  polynomial for quick phase
  fprintf(fd,"static const long double c7 = %1.50e;\n",coef_c7):
  fprintf(fd,"static const long double c3 = %1.50e;\n",coef_c3):
  fprintf(fd,"static const long double c2 = %1.50e;\n",coef_c2):
  fprintf(fd,"static const long double ps_alpha = %1.50e;\n",coef_alpha):
  fprintf(fd,"static const long double ps_beta = %1.50e;\n",coef_beta):
  fprintf(fd,"static const long double ps_c = %1.50e;\n",coef_c):
fprintf(fd, "#endif/* PATERSON*/\n\n"):

fprintf(fd, "#ifdef ESTRIN\n"):
 fprintf(fd, "  /* Coefficients are read directly from the array thanks to the following macros */ \n"):
  for i from PolyDegreeQuick to 2 by -1 do
    fprintf(fd, "#define c%d  c[%d]\n", i, PolyDegreeQuick-i):
  od:
   fprintf(fd, "static const double c[%d] =  {\n",PolyDegreeQuick):
   for i from PolyDegreeQuick to 2 by -1 do
     fprintf(fd, "   /* c%d  = */ %1.50eL, \n", i, coeff(polyQuick,x,i)):
   od:
  fprintf(fd, "}; \n \n"):
fprintf(fd, "#endif/* ESTRIN */\n\n"):



  for i from PolyDegreeAccurate to 1 by -1 do
    fprintf(fd, "#define c%dh  ch[%d]\n", i, PolyDegreeAccurate-i):
  od:

  for i from MaxDegreeDDE-1 to 1 by -1 do
    fprintf(fd, "#define c%dl  cl[%d]\n", i, MaxDegreeDDE-1-i):
  od:
  fprintf(fd, "#define PREFETCH_POLY_ACCURATE \n"):
  fprintf(fd, "\n#else /* not(defined(CRLIBM_TYPECPU_X86) || defined(CRLIBM_TYPECPU_AMD64)),\n   assuming Itanium, otherwise we shouldn't be there */ \n\n"):
  fprintf(fd, "#define PREFETCH_POLY_QUICK "):
  for i from PolyDegreeQuick to 2 by -1 do
    fprintf(fd, "c%d=c[%d]; ", i, PolyDegreeQuick-i):
  od:
  fprintf(fd, "\n#define PREFETCH_POLY_ACCURATE "):
  for i from PolyDegreeAccurate to 1 by -1 do
    fprintf(fd, "c%dh=ch[%d]; ", i, PolyDegreeAccurate-i):
    if i mod 4 =0 then  fprintf(fd, "\\\n        "): fi:
  od:
  fprintf(fd, "\\\n        "):
  for i from MaxDegreeDDE-1 to 1 by -1 do
    fprintf(fd, "c%dl=cl[%d]; ", i, MaxDegreeDDE-1-i):
  od:

  fprintf(fd, "\n#endif /* defined(CRLIBM_TYPECPU_X86) || defined(CRLIBM_TYPECPU_AMD64) */ \n\n"):

  # Various constants
  fprintf(fd, "static const long double log2H  = %1.50eL ;\n", log2h):
  fprintf(fd, "static const long double log2L  = %1.50eL ;\n", log2l):

  # The polynomials
  #  polynomial for quick phase
#  for i from PolyDegreeQuick to 1 by -1 do
#    fprintf(fd, "static const long double c%d =    %1.50eL ;\n", i, coeff(polyQuick,x,i)):
#  od:

  #  polynomial for accurate phase
  fprintf(fd, "static const long double ch[%d] =  {\n",PolyDegreeAccurate):
   for i from PolyDegreeAccurate to 1 by -1 do
     (ch, cl) := hiloExt(coeff(polyAccurate,x,i)):
     fprintf(fd, "   /* ch%d  = */ %1.50eL, \n", i, ch):
   od:
  fprintf(fd, "}; \n \n"):

  fprintf(fd, "static const long double cl[%d] =  {\n", MaxDegreeDDE):
  for i from MaxDegreeDDE-1 to 1 by -1 do
     (ch, cl) := hiloExt(coeff(polyAccurate,x,i)):
     fprintf(fd, "   /* cl%d  = */ %1.50eL, \n", i, cl):
   od:
  fprintf(fd, "}; \n \n"):


  # The tables


  fprintf(fd, "typedef struct rri_tag {float r; long double logirh;  long double logirl; } rri ;  \n"):
  fprintf(fd, "static const rri argredtable[%d] = {\n", 2^L):
  for i from 0 to 2^L-1 do
      fprintf(fd, "  { \n"):
      fprintf(fd, "    %1.50eL,   /* r[%d] */ \n", r[i], i):
      fprintf(fd, "    %1.50eL, /* logirh[%d] */ \n", logirh[i], i):
      fprintf(fd, "    %1.50eL, /* logirl[%d] */ \n", logirl[i], i):
      fprintf(fd, "  } "):
      if(i<2^L-1) then  fprintf(fd, ", \n"): fi
  od:
fprintf(fd, "}; \n \n"):

fclose(fd):

printf("\n************ DONE TEMPLOG/log-de.h ************\n");




filename:="TEMPLOG/polynomials.sed":
fd:=fopen(filename, WRITE, TEXT):
for i from PolyDegreeAccurate to 1 by -1 do
    (ch, cl) := hiloExt(coeff(polyAccurate,x,i)):
    fprintf(fd, " s/_c%dh/%1.40e/g\n", i, ch):
    if(i<MaxDegreeDDE) then
        fprintf(fd, " s/_c%dl/%1.40e/g\n", i, cl):
    fi:
od:
for i from PolyDegreeQuick to 1 by -1 do
    fprintf(fd, " s/_c%d/%1.40e/g\n", i, coeff(polyQuick,x,i)):
od:
fprintf(fd, " s/_deltaApproxQuick/%1.40e/g\n", deltaApproxQuick):
fprintf(fd, " s/_epsilonApproxQuick/%1.40e/g\n", epsilonApproxQuick):
fprintf(fd, " s/_deltaApproxAccurate/%1.40e/g\n", deltaApproxAccurate):
fprintf(fd, " s/_epsilonApproxAccurate/%1.40e/g\n", epsilonApproxAccurate):
fprintf(fd, " s/_log2h/%1.40e/g\n", log2h):
fprintf(fd, " s/_log2l/%1.40e/g\n", log2l):
fclose(fd):

printf("\n************ DONE TEMPLOG/polynomials.sed ************\n");

for j from 0 to 2^L-1 do
    filename:=cat("TEMPLOG/log-de_",j,".sed"):
    fd:=fopen(filename, WRITE, TEXT):
    fprintf(fd, " s/_rval/%1.40e/g\n",     r[j]):
    fprintf(fd, " s/_logirh/%1.40e/g\n", logirh[j]):
    fprintf(fd, " s/_logirl/%1.40e/g\n", logirl[j]):
    fprintf(fd, " s/_zabsmax/%1.40e/g\n", zabsmax[j]):
    fclose(fd):
  od:



printf("\n************ DONE TEMPLOG/log-de*.sed ************\n");

# A shell script to use them
filename:="../gappa/log-de/run-log-de-proof.sh":
fd:=fopen(filename, WRITE, TEXT):
fprintf(fd, "#!/bin/sh\n"):
fprintf(fd, "# You probably need to edit the path to the gappa executable\n"):
fprintf(fd, "GAPPA=~/local/src/gappa/src/gappa\n"):
fprintf(fd, "CRLIBMDIR=../..\n"):

fprintf(fd, "echo Accurate phase, case E=0, index=0: 1>&2\n"):
fprintf(fd, "sed  -f $CRLIBMDIR/maple/TEMPLOG/polynomials.sed  -f $CRLIBMDIR/maple/TEMPLOG/log-de_0.sed $CRLIBMDIR/gappa/log-de/log-de-acc-index0-E0.gappa | $GAPPA \n"):

fprintf(fd, "echo Accurate phase, case E!=0, index=0 1>&2\n"):
fprintf(fd, "sed  -f $CRLIBMDIR/maple/TEMPLOG/polynomials.sed  -f $CRLIBMDIR/maple/TEMPLOG/log-de_0.sed $CRLIBMDIR/gappa/log-de/log-de-acc-index0-E1N.gappa | $GAPPA \n"):

fprintf(fd, "for num in `seq 1 %d`; do\n", 2^L-1):
fprintf(fd, "  echo Accurate phase, case E=0, index=$num 1>&2\n"):
fprintf(fd, "  sed -f $CRLIBMDIR/maple/TEMPLOG/polynomials.sed  -f $CRLIBMDIR/maple/TEMPLOG/log-de_$num.sed $CRLIBMDIR/gappa/log-de/log-de-acc-index1N-E0.gappa | $GAPPA \n"):
fprintf(fd, "  echo 1>&2\n"):
fprintf(fd, "done\n"):

fprintf(fd, "for num in `seq 1 %d`; do\n", 2^L-1):
fprintf(fd, "  echo Accurate phase, case E!=0, index = $num 1>&2 \n"):
fprintf(fd, "  sed -f $CRLIBMDIR/maple/TEMPLOG/polynomials.sed  -f $CRLIBMDIR/maple/TEMPLOG/log-de_$num.sed $CRLIBMDIR/gappa/log-de/log-de-acc-index1N-E1N.gappa | $GAPPA \n"):
fprintf(fd, "  echo 1>&2\n"):
fprintf(fd, "done\n\n"):


fprintf(fd, "echo Quick phase, case E=0, index=0  1>&2\n"):
fprintf(fd, "sed  -f $CRLIBMDIR/maple/TEMPLOG/polynomials.sed  -f $CRLIBMDIR/maple/TEMPLOG/log-de_0.sed $CRLIBMDIR/gappa/log-de/log-de-index0-E0.gappa | $GAPPA \n"):
fprintf(fd, "  echo 1>&2\n"):

fprintf(fd, "echo Quick phase, case E!=0, index=0 1>&2 \n"):
fprintf(fd, "sed  -f $CRLIBMDIR/maple/TEMPLOG/polynomials.sed  -f $CRLIBMDIR/maple/TEMPLOG/log-de_0.sed $CRLIBMDIR/gappa/log-de/log-de-index0-E1N.gappa | $GAPPA \n"):
fprintf(fd, "  echo 1>&2\n"):


fprintf(fd, "for num in `seq 1 %d`; do\n", 2^L-1):
fprintf(fd, "  echo Quick phase, for all E,  index=$num 1>&2\n"):
fprintf(fd, "  sed  -f $CRLIBMDIR/maple/TEMPLOG/polynomials.sed  -f $CRLIBMDIR/maple/TEMPLOG/log-de_$num.sed $CRLIBMDIR/gappa/log-de/log-de-index1N-E0N.gappa | $GAPPA \n"):
fprintf(fd, "  echo 1>&2\n"):
fprintf(fd, "done\n"):

fclose(fd):

printf("\n************ DONE TEMPLOG/run-log-de-proof.sh ************\n"):
printf("Now you should run \n"):
printf(" sh ../../gappa/log-de/run-log-de-proof.sh  2> ../../maple/TEMPLOG/Gappa.out\n"):
printf("  (You probably need to edit the path to the gappa executable within run-log-de-proof.sh)\n"):
printf("Then look at TEMPLOG/Gappa.out. It shouldn't contain 'some enclosures were not satisfied'.\n If it does, report it !\n"):


printf("----DONE---\n") :

