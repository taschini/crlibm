#######################################################################
# To use:
# restart; read "log.mpl";
Digits := 100:

interface(quiet=true):

read "common-procedures.mpl":
mkdir("TEMPLOG"):
deg:=12:
# Error, (in mkdir) directory exists and is not empty


#######################################################################
# First, quick phase

# Approach for the first function:
# log(x) = E*log(2) + P(R)
# where P(R) is defined over such an interval [11/16; 12/16] ... [21/16; 23/16]
# At all, there is 8 intervals, corresponding to 8 polynomials.

# Computation of log(2) as a double-double
# We want to compute E.log(2) where |E|<1024, and get the result as a double-double.
# We store two double constants ln2hi and ln2lo, such that ln2hi has at most 43 non-zero bits in its mantissa,
# such that
# ln2hi+ln2lo = log(2) (1+2**-101) (see below log2(eps_ln2_cst) )
# and the multiplication of ln2hi by E is exact. Then we use a Fast2Sum to compute E.log2 as a double-double.

ln2hi:= round( ln(2.) * (2**43))    /    (2**43) :
r:=evalf(ln(2) - ln2hi):
bits_log2hi:=1+log2(op(2,ieeedouble(ln2hi)[3])):
ln2lo:=nearest(ln(2.) - ln2hi):
eps_ln2_cst:=abs( (ln2hi+ln2lo-ln(2.)) / ln(2.)):
log2(eps_ln2_cst):



# Performance parameters

# Try adding one or removing one to bitsFastPath. This value is optimal for the Pentium IV
bitsFastPath:=7:

#idem for MedPath. Look at the proof for explanations.
bitsMedPath:=4:

minFastPath:=(2^bitsFastPath)+0.5: # such that   |Elog(2) + p(z)| > 2^bitsFastPath   where |p(z)| < 0.4
EminFastPath:=floor(minFastPath/log(2)): # such that |E| > EminFastPath   => |Elog(2) + p(z)| > 2^bitsFastPath
evalf((EminFastPath+1)*log(2)-0.4): # just checking

minMedPath:=(2^bitsMedPath)+0.5: # such that   |Elog(2) + p(z)| > 2^bitsMedPath   where |p(z)| < 0.4
EminMedPath:=floor(minMedPath/log(2)): # such that |E| > EminMedPath   => |Elog(2) + p(z)| > 2^bitsMedPath
evalf((EminMedPath+1)*log(2)-0.4): # just checking



# Interval parameters

ymin := (1-2^(-50)) * sqrt(2)/2: # sqrt(2)/2 with margin for the FP comparison and for the computation of xmax
ymax := (1+2^(-50)) * sqrt(2): # sqrt(2)/2 with margin for the FP comparison and for the computation of xmax

Ilist:=[[ymin,12/16],[12/16,13/16],[13/16,14/16],[14/16,15/16],[15/16,17/16],[17/16,19/16],[19/16,21/16],[21/16,ymax]]:
Xmax:=[(12/16-ymin)/2, 2**(-5), 2**(-5), 2**(-5), 2**(-4), 2**(-4), 2**(-4), (ymax - 21/16)/2]:
midI:=[nearest((ymin+12/16)/2), 25/32, 27/32, 29/32, 1,18/16, 20/16, nearest((21/16+ymax)/2)]:


# Computation of the polynomials

# This is the function we use on all the intervals

poly_log_2 := proc(i,deg)
 local p, pe, e, pt, delta_approx, epsilon_approx, delta_rounding, minptr, err,
       maxptr, maxpt, epsilon, miny, epsilonZero, epsilonEgtMed, ptFastPath,epsilonFastPath,
       minptrFastPath, maxptrFastPath, errlist, epsilon_lastmult, maxEln2_lo,
	s1, p1, maxres, c0h, c0l, maxP_hi, maxP_lo,delta_rounding_s1,maxp1,delta_rounding_p1:

  printf("Interval %d  : ", i):
  if(i=5) # special case around 1, we impose that coeff0 = 0
  then
    pe:=x * numapprox[minimax](  (log(1+x)/x),  x=-Xmax[i]..Xmax[i],  [deg-1,0], 1 ,  'delta_approx'):
    pt := poly_exact2(pe, 2):
    delta_approx := numapprox[infnorm](pt-log(1+x), x=-Xmax[i]..Xmax[i]):
    epsilon_approx := numapprox[infnorm]( 1-pt/log(1+x), x=-Xmax[i]..Xmax[i]):
    maxpt:= numapprox[infnorm]( pt, x=-Xmax[i]..Xmax[i]):
  else
    pt, delta_approx, epsilon_approx,  maxpt := poly_trunc_f2d_2( log(x+midI[i]) ,  deg, -Xmax[i], Xmax[i] ):
  end if:

  printf(" delta_approx = %3.2f, epsilon_approx = %3.2f, maxpt = %3.2f\n", -log2(delta_approx), -log2(epsilon_approx), maxpt):

  #delta for E=0
  errlist := errlist_quickphase_horner(deg, 2, 2, 0,0): # two dd adds, two dd mul, no error on x
  epsilon_lastmult, delta_rounding, minptr, maxptr := compute_horner_rounding_error(pt, x, Xmax[i], errlist, true):
  if i=5 then
    epsilonZero := (1+epsilon_approx)  * (1 + epsilon_lastmult) - 1 :
  else
    epsilonZero := (delta_approx + delta_rounding)/minptr:
  fi:

  printf(" epsilonZero =  %3.2f ", -log2(epsilonZero)):


  # For 0<E<= EminFastPath, the Horner evaluation ends in double-double
  errlist := errlist_quickphase_horner(deg, 2, 1, 0,0): # two dd adds, one dd mul, no error on x
  epsilon_lastmult, delta_rounding, minptr, maxptr := compute_horner_rounding_error(pt, x, Xmax[i], errlist, true):

  # epsilon for 0<E<= EminMedPath

  miny:=nearest(log(2.)) - maxptr:       #  worst-case miny
  # and add to the absolute error that of computing Elog2 : 2**(-90)
  epsilon:=evalf(  (delta_approx + delta_rounding + 2**(-90) )/miny  )  :
  printf(" epsilon =  %3.2f ", -log2(epsilon)):

  #delta for EminMedPath < E <= EminFastPath

  miny := (EminMedPath+1)*log(2.) - maxptr:
  epsilonEgtMed := (delta_approx + delta_rounding + 2**(-90) ) /miny:
  printf(" epsilonEgtMed =  %3.2f ", -log2(epsilonEgtMed)):

  #epsilon for the fast path : in this case the polynomial is rounded to double, and evaluated only in double
  ptFastPath := poly_exact2(pt,1): # only coeff of degree 0 in double-double
  delta_approx :=  numapprox[infnorm](ptFastPath-log(x+midI[i]), x=-Xmax[i]..Xmax[i]):
  s1 := expand((ptFastPath - coeff(ptFastPath,x,0)) / x): # the polynomial computed in double
  errlist := errlist_quickphase_horner(deg-1, 0, 0, 0,0): # no dd adds, no dd mul, no error on x
  epsilon_lastmult, delta_rounding_s1, minptr, maxptrFastPath := compute_horner_rounding_error(s1, x, Xmax[i], errlist, true):
  p1 := s1*x :
  maxp1:=numapprox[infnorm](p1, x=-Xmax[i]..Xmax[i]):
  delta_rounding_p1 := delta_rounding_s1*Xmax[i] + 0.5*ulp(maxp1):   # the last mult by z
  c0h,c0l := hi_lo(coeff(ptFastPath, x, 0)):
  maxP_hi := 1075*log(2) + c0h:
  maxP_lo := maxP_hi*2^(-53):
  maxEln2_lo := maxP_lo:
  # the delta is that of the second argument of the last Add12.
  delta_rounding := delta_rounding_p1
                      +  0.5*ulp(maxp1 + c0l + maxEln2_lo + maxP_lo)
                      +  0.5*ulp(c0l + maxEln2_lo + maxP_lo)     # these two last terms are zero in the case i=5
                      +  0.5*ulp(maxEln2_lo + maxP_lo):          #  but it doesn't change much
  miny := (EminFastPath+1)*log(2.) - maxp1 :
  epsilonFastPath := (delta_approx + delta_rounding + 2**(-90) ) / miny :
  printf(" epsilonFastPath = %3.2f\n", -log2(epsilonFastPath) ):
  [pt, epsilon_approx, max(maxpt,maxptr), epsilonZero, epsilon, epsilonEgtMed, epsilonFastPath]
end proc:






printf("Calcul de PolyList\n"):

PolyList:=[
  poly_log_2(1,deg),
  poly_log_2(2,deg),
  poly_log_2(3,deg),
  poly_log_2(4,deg),
  poly_log_2(5,deg),
  poly_log_2(6,deg),
  poly_log_2(7,deg),
  poly_log_2(8,deg)
]:



save PolyList, "TEMPLOG/PolyList.m":



# Computation of constants for RN test
# We have in PolyList epsilon the total error for polynomial approximation, build tabrndcst[] the table of "e" needed for round to nearest

maxepsilonEZero:=0:
maxepsilon:=0:
maxepsilonEgtMed:=0:
maxepsilonFastPath:=0:

for i from 1 to 8 do
  epsilonEZero:=evalf(PolyList[i][4]):
  epsilon:=evalf(PolyList[i][5]):
  epsilonEgtMed:=evalf(PolyList[i][6]):
  epsilonFastPath:=evalf(PolyList[i][7]):
  if epsilonEZero > maxepsilonEZero then maxepsilonEZero := epsilonEZero : fi :
  if epsilon > maxepsilon then maxepsilon := epsilon : fi :
  if epsilonEgtMed > maxepsilonEgtMed then maxepsilonEgtMed := epsilonEgtMed : fi :
  if epsilonFastPath > maxepsilonFastPath then maxepsilonFastPath := epsilonFastPath : fi :
od:
rncstEZero := evalf(compute_rn_constant(maxepsilonEZero)):
rncst := evalf(compute_rn_constant(maxepsilon)):
rncstEgtMed := evalf(compute_rn_constant(maxepsilonEgtMed)):
rncstFastPath := evalf(compute_rn_constant(maxepsilonFastPath)):







#-------------------------------------------------------------------
# Output


filename:="TEMPLOG/log_fast.h":
fd:=fopen(filename, WRITE, TEXT):

fprintf(fd, "#include \"crlibm.h\"\n#include \"crlibm_private.h\"\n"):

fprintf(fd, "\n/*File generated by maple/log.mpl*/\n"):

fprintf(fd, "\n\#define SQRT_2 1.4142135623730950489e0 \n\n"):

fprintf(fd, "#define DEGREE %d\n\n", deg):

fprintf(fd,"#define EMIN_MEDIUMPATH %4.3f \n",EminMedPath ):
fprintf(fd,"#define EMIN_FASTPATH %4.3f \n",EminFastPath ):

fprintf(fd,"/* Constants for rounding  */\n"):

fprintf(fd,"static const double epsilon[4] =\n{\n"):
fprintf(fd,"  /* Case E=0  */\n"):
fprintf(fd,"  %1.50e, \n",maxepsilonEZero):
fprintf(fd,"  /* Middle case   */\n"):
fprintf(fd,"  %1.50e, \n",maxepsilon):
fprintf(fd,"  /* Case E>EminMedPath  */\n"):
fprintf(fd,"  %1.50e, \n",maxepsilonEgtMed):
fprintf(fd,"  /* And for the fast path  */\n"):
fprintf(fd,"  %1.50e \n",maxepsilonFastPath):
fprintf(fd,"\n};\n\n"):

fprintf(fd,"static const double rncst[4] =\n{\n"):
fprintf(fd,"  /* Case E=0  */\n"):
fprintf(fd,"  %1.50e, \n",rncstEZero):
fprintf(fd,"  /* Middle case   */\n"):
fprintf(fd,"  %1.50e, \n",rncst):
fprintf(fd,"  /* Case E>EminMedPath  */\n"):
fprintf(fd,"  %1.50e, \n",rncstEgtMed):
fprintf(fd,"  /* And for the fast path  */\n"):
fprintf(fd,"  %1.50e \n",rncstFastPath):
fprintf(fd,"};\n\n"):


# THE POLYNOMIALS

npol:=8:
n:=2:  # the degree from which we want two doubles

fprintf(fd,"#ifdef WORDS_BIGENDIAN\n"):
for isbig from 1 to 0 by -1 do

  if(isbig=0) then
    fprintf(fd,"#else\n"):
  fi:

  # Various constants
  fprintf(fd, "static const db_number ln2hi = "):
  printendian(fd,ln2hi,isbig):
  fprintf(fd, ";\n"):

  fprintf(fd, "static const db_number ln2lo = "):
  printendian(fd,ln2lo,isbig):
  fprintf(fd, ";\n"):

  fprintf(fd, "static const db_number two52 = "):
  printendian(fd, nearest(2.**52),isbig):
  fprintf(fd, ";\n"):


  # Write middle table

  fprintf(fd, "static const db_number middle[%d] =\n{\n",npol):
  for i from 1 to 7 do
    printendian(fd,midI[i],isbig):
    fprintf(fd," ,\n"):
  od:
  printendian(fd,midI[8],isbig):
  fprintf(fd, "\n};\n\n"):

  #Write Poly_h

  fprintf(fd,"static const db_number Poly_h[%d][%d] =\n{\n", npol, deg+1):
  for k from 1 to 8 do
    fprintf(fd," /* polynomial %d */\n",k):
    fprintf(fd,"{\n"):
    P:=PolyList[k][1]:
    for j from 0 to deg do
      coef:=hi_lo(coeff(P,x,j)):
      printendian(fd, coef[1], isbig):
      if j<deg then fprintf(fd," ,\n") fi:
    od:
    fprintf(fd, "\n}"):
    if k<8 then fprintf(fd, ",\n") fi:
  od:

  fprintf(fd, "\n};\n\n"):

  #Write Poly_fast_l

  fprintf(fd,"static const db_number Poly_l[%d][%d] =\n{\n",8,n):
  for k from 1 to 8 do
    fprintf(fd," /* polynomial %d */\n",k):
    fprintf(fd,"{\n"):
    P:=PolyList[k][1]:
    for j from 0 to n-1 do
      coef:=hi_lo(coeff(P,x,j)):
      printendian(fd, coef[2], isbig):
      if j<deg then fprintf(fd," ,\n") fi:
    od:
    fprintf(fd, "\n}"):
    if k<8 then fprintf(fd, ",\n") fi:
  od:
  fprintf(fd,"\n};\n\n"):

od:
fprintf(fd,"#endif\n\n\n"):



fclose(fd):

# Output of latex for the report

printf("Polynomial & Relative Approximation Error \\\\\n"):
   for k from 1 to 8 do
     P:=PolyList[k]:
     printf("P[%d] &  %2.2f \\\\ \n", k, -log2(P[2])):
    od:


   printf("Polynomial & Maxp  & Epsilon0 & Epsilon & EpsilonMed & EpsilonFast \\\\\n"):
   for k from 1 to 8 do
     P:=PolyList[k]:
     printf("P[%d] & %2.3f & %2.2f & %2.2f& %2.2f  & %2.2f\\\\ \n",k, P[3]+0.0005, -log2(P[4]), -log2(P[5]), -log2(P[6]), -log2(P[7])):
    od:






######################################################################
# Second, accurate phase in SCS
# Approach for the scs function:


#  Who wrote it ? Please comment and clean up.

#  x=2^e (1+f)
#  log(x) = e.log(2) + log(1+f)

#   log(1+f) = log(w) + log(1+(1+f-w)/w)
if(1+1=3) then
W := 2^5:
Poly_P := series((ln(1+x)/ln(2.))/x, x=0, 45):
Poly_Q := convert(Poly_P,polynom):
Poly_cheb := chebpade(Poly_Q, x=-1/W..1/W, [19,0]):
Poly_Res  := sort(expand(x * Poly_cheb)):
log2(numapprox[infnorm]((1-(expand(Poly_Res))/(ln(1+x)/ln(2.))), x=-1/W..1/W)):
# This procedure gives under a polynom form value needed by tabular range reduction for the logarithm.
# DON'T forget the zero !!!!

Digits       := 50:
table_ti     := 0:
table_inv_wi := 0:
N            := -4:
break_point  := sqrt(2.):
start_nb     := round((break_point/2) * 2.^(-N))*2.^(N):
nb_turn      := round((break_point/2) * 2.^(-N)):

for i from 0 to nb_turn do
   table_ti := table_ti * x + log(start_nb):
   table_inv_wi := table_inv_wi *x + 1./ start_nb:
   start_nb := start_nb + 2.^(N):
od:

fi:

# Tests and scratch






