restart:
Digits := 100:
with (numapprox):with(orthopoly):
read "common-procedures.mpl";
mkdir("TEMPATAN");

#First, quick phase
#Approach for the first function:
#atan(x) = atan(b(i)) + atan ( (x-bi) / (1+x*b(i) ) )
#We try to have 64 correct bits



#Constants


#if (x >= xmax) : return arctan(x)/2)=neareast(Pi/2)

log[2](evalf(min(solve(arctan(x)=nearest(Pi/2),x)))):
xmax:=2^54:


#halfPi=Pi/2 (remark that maxvalue < Pi/2 so that atan(x) < Pi/2

halfPi:=nearest(Pi/2):
halfPi_to_plus_infinity:=nearest(nearest(Pi/2)*(1+2^(-53))):


# invPi, for the atanpi function

invpih,invpil := hi_lo(1/Pi);       



# if (xmin<xmin) : return arctan(x)=x

xmin=2^(-27);

ieeehexa(halfPi); ieeehexa( halfPi_to_plus_infinity ); ieeehexa(2^(-27));ieeehexa(2^(54));ieeehexa(2^(-10));


#Interval parameters
#Some constants :
e:=2^(-6.3):
marge := 2^(-30):

#The intervals are (a[i],a[i+1]). the "middle" is b(i) :
maxx := bi -> simplify( solve( (x-bi) / (1+x*bi) = e ,x) ):
minx := bi -> simplify(solve( (x-bi) /(1+x*bi)=-e ,x)):

nextbi := proc (x) evalf( max(solve( minx(bi) = x ,bi) )*(1-marge)); end:

allbi := proc (n)
  local xi,nbi,x,i,j;
  global b,a, nb_of_ai, nb_of_bi, value_for_dicho;
  x := e;
  nbi := 0;
  i := 0;
  while(i<n and nbi < 1/e and nbi >= 0) do
     nbi := nearest ( nextbi(x) );
     b[i] := evalf( nbi );
     a[i] := x;
     x   := evalf(maxx(nbi));
     i:=i+1;
  od;
  j:=0;
  while ( 2^j < i ) do j:=j+1 od:
  nb_of_ai := i;
  nb_of_bi := i;
  b[i-1] := nearest(1/e+4):
  value_for_dicho := 2^j;
  return i,b[0],b[i-1];
end:

allbi(100);




#The polynome and its size ( 4 because 4 terms : x^2/3, ...)

Q9:= poly_exact(x-1/3*x^3+1/5*x^5 -1/7*x^7+1/9*x^9):
deg_poly := 4:
coef_poly[0]:=nearest (-1/3):
coef_poly[1]:=nearest (1/5):
coef_poly[2]:=nearest (-1/7):
coef_poly[3]:=nearest (1/9):

Q := x^2 * coef_poly[0] + x^4*coef_poly[1] + x^6*coef_poly[2] + x^8 * coef_poly[3]:

Qprime := poly_exact( x * coef_poly[0] + x^2*coef_poly[1] + x^3*coef_poly[2] + x^4 * coef_poly[3]):

log[2] (infnorm( (arctan(x)-x*(1+Q) )/x,x=2^(-53)..e));
# we can note that arctan(x) = sum((-1)^i*x^(2*i+1)/(2*i+1), i = (0 .. infinity))
# so ( arctan(x) - Q9(x) ) / x = sum((-1)^i*x^(2*i)/(2*i+1), i = (5 .. infinity))
# so abs ( (arctan(x) -Q9(x) )/x )  <= x^10 * sum((-1)^(i+1)*x^(2*i)/(2*i+1+10), i = (0 .. infinity))
#                                   <= x^10 * 1/11
#                                   <= 2^(-66) if x<2^(-6.3)









#Conputation of the RN constant


#Epsilon = relative error
#Delta = absolute error
#1st : Error about reduction:

EpsilonXminusBi := 2^(-105):

# 1 + x*b[i] :
Epsilon1A22xbi :=  2^(-104):


# Xred = (x-b[i]) / ( 1 + x*b[i] )

EpsilonXred := EpsilonXminusBi+Epsilon1A22xbi + 2^(-105):

DeltaXred := EpsilonXred * e:

log2(EpsilonXred);

#Polynomial evaluation :
#atan ~= x - x^3/3 + x^5/5 - x^7/7 + x^9/9
#        ~= x * ( 1 + Q(x^2) )
#    = (xhi+xlo) * ( 1 + Q)
#        ~= xhi + (xhi*Q + xlo)

#  Warning : the approx about x !!

#calc of Q:

EpsilonXred2 := (1+2^(-53) + 2^(-105) )^2 * (1+2^(-53)) -1 :

log2(EpsilonXred2);

errlist:=errlist_quickphase_horner( degree( x^4 ),0,0,EpsilonXred2, 2**(-53)):

rounding_error1:= compute_horner_rounding_error(Qprime,x,e,errlist,true):

EpsilonQ := rounding_error1[1]:
log2(EpsilonQ);

# Since q will be multiplied in x(1+q), its absolute error term is :

Qmax:=infnorm(Qprime,x=-e^2..e^2):

DeltaQ := e*Qmax*EpsilonQ;
log2(DeltaQ);

# calc of xi*Q + xlo + atan(b[i])lo

#First the truncation of xredlo*q : at most half an ulp of xredhi*q
DeltaTruncXredloQ :=  2^(-53) * e * Qmax:
log2(DeltaTruncXredloQ);
#then the sum of atanbilo and xredlo : the largest term is atanbilo
DeltaAdd1  := 0.5*ulp(2^(-53)*(Pi/2 + e)):
log2(DeltaAdd1);

#then compute xredhi*q
DeltaXredhiQ := 0.5*ulp(e*Qmax):
log2(DeltaXredhiQ);

# and the second addition
DeltaAdd2 :=  0.5*ulp(e*Qmax +  2^(-53)*(Pi/2 + e) ):
log2(DeltaAdd2);

DeltaAtanlolo := DeltaTruncXredloQ + DeltaAdd1 + DeltaXredhiQ + DeltaAdd2:
log2(DeltaAtanlolo);

# In the second reconstruction there is one more FP add which is
# aligned to the previous:
DeltaReconst := DeltaAtanlolo + DeltaAdd2 :
log2(DeltaReconst);



# error due to the approx of atan by P : < e^11/11
Deltaapprox := infnorm( (arctan(x)-x*(1+Q)), x=0..e):

log[2](Deltaapprox);





DeltaTotal := Deltaapprox + DeltaQ + DeltaReconst :
log2(DeltaTotal);
EpsilonTotal := evalf( DeltaTotal / arctan(xx) ):

EpsilonFinal := infnorm( EpsilonTotal,xx=2^(-6.3)..1,'xmaxmax'):
EpsilonFinal_i_10 := infnorm( EpsilonTotal,xx=a[10]..2^60,'xmaxmax'):
log2( EpsilonFinal );
log2( EpsilonFinal_i_10 );


# Computation of e and e_i_30 (which is e when i>30 (ie x>1)
E := evalf(compute_rn_constant(EpsilonFinal));
E_i_10 := evalf(compute_rn_constant(EpsilonFinal_i_10));


#An other Constant : when there is no reduction.
#We compute the polynom and then atan(x) = x + x.Q

Epsilonx2 := 2^(-53):
errlist:=errlist_quickphase_horner( degree( Qprime ),0,0,Epsilonx2, 2**(-53)):
rounding_error:= compute_horner_rounding_error(Qprime,x,e,errlist,true):

epsilon_Q := rounding_error[1]:
log[2](epsilon_Q);

epsilon_x_Q := 2^(-52):

delta_final := 2^(-105)*x + epsilon_Q*x^3 + epsilon_x_Q*x^3 - arctan(x) + x*(1+Q):
epsilon_final := infnorm(delta_final/x,x=2^(-10)..e):
epsilon_final_m10 := infnorm(delta_final/x,x=2^(-24)..2^(-10)):

log[2](epsilon_final);
E_no_reduction := evalf(compute_rn_constant(epsilon_final));
E_no_reduction_m10 := evalf(compute_rn_constant(epsilon_final_m10));



# Output


#test
filename:="TEMPATAN/atan_fast.h":

fd:=fopen(filename, WRITE, TEXT):
fprintf(fd, "\n/*File generated by maple/atan.mpl */\n"):

fprintf(fd, "#ifndef _ATAN_FAST_H\n#define _ATAN_FAST_H\n\n");

fprintf(fd, "#include \"crlibm.h\"\n#include \"crlibm_private.h\"\n"):

fprintf(fd,"#ifdef WORDS_BIGENDIAN\n"):
fprintf(fd, "static const db_number HALFPI = {{0x3FF921FB,0x54442D18}};"):
fprintf(fd, "static const db_number HALFPI_TO_PLUS_INFINITY = {{0x3FF921FB,0x54442D19}};"):
fprintf(fd,"\n#else\n"):
fprintf(fd, "static const db_number HALFPI = {{0x54442D18,0x3FF921FB}};"):
fprintf(fd, "static const db_number HALFPI_TO_PLUS_INFINITY = {{0x54442D19,0x3FF921FB}};"):
fprintf(fd,"\n"):
fprintf(fd,"#endif\n"):

fprintf(fd, "#define MIN_REDUCTION_NEEDED %1.50f\n",e):

fprintf(fd, "#define INVPIH %1.50f\n", invpih):
fprintf(fd, "#define INVPIL %1.50f\n", invpil):

fprintf(fd,"#define nb_of_ai %d\n",nb_of_ai):
fprintf(fd,"#define nb_of_bi %d\n",nb_of_bi):

# we multiply by 2^(-20) in order to prevent from rounding
fprintf(fd, "static const double rncst[4] ={\n"):
fprintf(fd, "%1.50f , /* i<10 */ \n",E*(1+2^(-20))) :
fprintf(fd, "%1.50f , /* i>10 */ \n",E_i_10*(1+2^(-20))) :
fprintf(fd, "%1.50f , /* e > 2^-10 */ \n",E_no_reduction*(1+2^(-20))) :
fprintf(fd, "%1.50f , /* e < 2^-10 */ \n };\n",E_no_reduction_m10*(1+2^(-20))) :

fprintf(fd, "static const double epsilon[4] ={\n"):
fprintf(fd, "%1.50e ,\n",EpsilonFinal *(1+2^(-20))) :
fprintf(fd, "%1.50e ,\n",EpsilonFinal_i_10*(1+2^(-20))) :
fprintf(fd, "%1.50e ,\n",epsilon_final*(1+2^(-20))) :
fprintf(fd, "%1.50e ,\n };\n",epsilon_final_m10*(1+2^(-20))) :

fprintf(fd, "#define DEGREE %d\n", deg_poly):
fprintf(fd, "static double const coef_poly[%d] = \n{\n",deg_poly):

for i from deg_poly-1 to 0 by -1 do
    fprintf(fd, "/* coef for degree %d */ %1.50f, \n" , 2*i+3, coef_poly[i]):
od:
fprintf(fd," }; \n");

fprintf(fd,"#define A 0\n"):
fprintf(fd,"#define B 1\n"):
fprintf(fd,"#define ATAN_BHI 2\n"):
fprintf(fd,"#define ATAN_BLO 3\n"):

fprintf(fd,"\n#ifdef WORDS_BIGENDIAN\n "):

for isbig from 1 to 0 by -1 do

    if(isbig=0) then
        fprintf(fd,"#else\n");
    fi:

    if(not (nb_of_ai = nb_of_bi)) then
        printf("Warning : nb_of_ai != nb_of_bi, this should not work");
    fi:

    fprintf(fd,"\n/* limits of the intervals [a[i],b[i]] */\n");
    fprintf(fd, "static db_number const arctan_table[%d][4] = \n{\n" , nb_of_ai );

   for i from 0 to nb_of_ai - 1 do
        fprintf(fd, "{\n/*a[%d]     */ ",i);
        printendian(fd,a[i],isbig);
        fprintf(fd," ,\n/*b[%d] :   */ ",i):
        printendian(fd,b[i],isbig):
        fprintf(fd," ,\n/*atan_b[%d]*/ ",i):
        printendian(fd, arctan(b[i]) ,isbig):
        fprintf(fd,",");
        printendian(fd, arctan(b[i])-nearest(arctan(b[i])) ,isbig);
        fprintf(fd," ,\n}\n,");
    od;
    fprintf(fd,"\n};\n");
od:

fprintf(fd,"\n#endif\n\n"):

fprintf(fd, "#endif /* def _ATAN_FAST_H */\n");

fclose(fd);



# Scs phase

#Now we want a precision of 130 bits in order to have correct rounding in all cases.
evalf(log[2](e^21/21));
#Test about Rounding :
#error on the polynom :
P_scs := convert(series(arctan(x),x=0,21),polynom);
evalf(log[2](e^21/21));
#Intervals parameters
#We choose the same e and the same intervals as in the quick phase.

#We have the same b[i]s and a[i]s but we save atan[b[i]] on 4 doubles

# Output
filename:="TEMPATAN/atan_accurate.h":
fd1:=fopen(filename, WRITE, TEXT):

fprintf(fd1, "#include \"crlibm.h\"\n#include \"crlibm_private.h\" \n#include \"atan_fast.h\"\n"):
fprintf(fd1, "\n/*File generated by maple/atan.mpl */\n"):

Write_SCS_poly(fd1, "constant_poly", P_scs);

fprintf(fd1,"#define constant_poly_ptr (scs_ptr)&constant_poly\n");

# The 1/Pi SCS constant
  fprintf(fd, "static const scs InvPiSCS=\n"):
  WriteSCS(fd, evalf(1/Pi)):
  fprintf(fd, ";\n#define InvPiSCS_ptr  (scs_ptr)(& InvPiSCS)\n\n"):


fprintf(fd1,"#ifdef WORDS_BIGENDIAN\n\n"):
for isbig from 1 to 0 by -1 do
   if isbig = 0 then
      fprintf(fd1, "\n#else\n"):
   fi:

   fprintf(fd1,"static const db_number atan_blolo[%d] = {\n",nb_of_ai);
   for i from 0 to nb_of_ai-1 do
      fprintf(fd1,"/* %d */  ",i):
      temp0 := nearest(arctan(b[i]));
      temp := nearest( arctan(b[i]) - nearest(arctan(b[i])));
      temp1 := nearest( arctan(b[i]) - temp -temp0);
      printendian(fd1,temp1,isbig):
      fprintf(fd1,", \n");
   od:

   fprintf(fd1,"};\n"):

od:
fprintf(fd1,"\n#endif  /* WORDS_BIGENDIAN */ \n");

fclose(fd1);
