restart:
Digits := 100:
 interface(quiet=true);
read "Common_maple_procedures.m";
read "temp_common.mpl";


mkdir("TEMPSIN");
# Error, (in mkdir) directory exists and is not empty

#################################
#Comments and todos : 

#don't forget the error on k in Cody and Waite



#########################################################
# The polynomials 

#########################################################
# Small arguments

# When do we return x ?
xmax_return_x_for_sin := 1.49e-8;  # TODO Demander à Muller comment il l'a calculé


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
errlist:=errlist_quickphase_horner(degree(polySinFast2),0,0,2**(-53), 2**(-70)):
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








#################################################
# CODY and WAITE  Argument reduction

# This parameter may be adjusted. 
bits_pio256hi0:=37:

invpio256:= nearest(256/Pi);
reminvpio2 := evalf(256/Pi - invpio256);

pio256hi:= round(evalf(  (Pi/256) * 2^(bits_pio256hi0))) / (2^(bits_pio256hi0));
r:=evalf(Pi/256 - pio256hi):
bits_pio256hi:=1+log2(op(2,ieeedouble(pio256hi)[3])) ;  # this means the log of the denominator

pio256lo:=nearest(Pi/256 - pio256hi);

# Cody and Waite argument reduction will work for k<kmax
kmax:=2^(53-bits_pio256hi);
xmax_cody_waite := kmax*Pi/256;
# Doing a 3-part Cody and Waite is needed in some rare cases
pio256med1:= round(r * (2^(2*bits_pio256hi))) / (2^(2*bits_pio256hi));
bits_pio256med1:=1+log2(op(2,ieeedouble(pio256med1)[3]));

r:=evalf(Pi/256 - pio256hi - pio256med1):
pio256med2:= round(r * (2^(3*bits_pio256hi))) / (2^(3*bits_pio256hi));
bits_pio256med2:=1+log2(op(2,ieeedouble(pio256med2)[3]));
r:=evalf(Pi/256 - pio256hi - pio256med1 - pio256med2):
pio256lo2:=nearest(r);

log2(Pi/256 -  pio256hi - pio256med1  - pio256med2 - pio256lo2);


# the only error in the Add12 that computes the range reduction is the following
rel_error_2nd_mul_CodyWaite:=2**(-53-bits_pio256hi);

# which may be multiplied by kmax 
# TODO it may be worth splitting the interval of k to optimize the rounding constants

abs_error_CodyWaite:= kmax * rel_error_2nd_mul_CodyWaite; # = 2**(1-2*bits_pio256hi)

# je ne comprends pas cette ligne
miny_accurate:=2**(1-2*bits_pio256hi+53+abs_error_CodyWaite);

rel_error_CodyWaite:= abs_error_CodyWaite/miny_accurate;

log2(rel_error_CodyWaite);


 #TODO ca peut pas être vrai, cela
rel_error_CodyWaite := 2**(-70);




###########################################################################################
#            Payne and Hanek argument reduction






################################## Polynomials for sine and cos 

ymax:=Pi/512;
y2max:= ymax**2; 

polySin:=  poly_exact (convert( series(sin(x), x=0, 9), polynom),x);
polyTs :=  expand(polySin/x-1);
deltaTs := numapprox[infnorm](polyTs - (sin(x)/x-1), x=-ymax..ymax):
maxTs := numapprox[infnorm](polyTs, x=-ymax..ymax):
log2(deltaTs); 
log2(maxTs);



polyCos:= poly_exact(subs(y=x^2, numapprox[minimax]((cos(sqrt(y))), y=2^(-2048)..y2max, [3,0])), x);
deltaCos:= numapprox[infnorm](polyCos -  cos(x), x=-ymax..ymax):
log2(deltaCos);
polyTc:=polyCos - 1;
deltaTc := numapprox[infnorm](polyTc -  (cos(x)-1), x=-ymax..ymax):
maxTc := numapprox[infnorm](polyTc, x=-ymax..ymax):
log2(deltaTc);
log2(maxTc);


if(1+1=3) then
# We get a better polynomial above
 polyCos:= poly_exact(convert( series(cos(x), x=0, 8), polynom),x);
 polyTc:=expand(polyCos-1);
 deltaTc := numapprox[infnorm](polyTc -  (cos(x)-1), x=-ymax..ymax):
 maxTc := numapprox[infnorm](polyTc, x=-ymax..ymax):
 log2(deltaTc);
 log2(maxTc);
end if;



# TODO add the error due to rounding k

# evaluate this polynomial in double. 

# The error on yh*yh 
# we had $y_h+y_l = y + \abserr{CodyWaite}$. 
# Now we take only $y_h$ :  $y_h = (y_h + y_l)(1+2^{-53}) = (y+\abserr{CodyWaite})(1+2^{-53})$
# When squared we get $y_h\otimes y_h = (y+\abserr{CodyWaite})^2(1+2^{-53})^3 $

epsy2 := ((1+rel_error_CodyWaite)**2) *  (1+2**(-53))**2 - 1;
 
log2(epsy2);

errlist:=errlist_quickphase_horner(degree(polyDoSin2),0,0, epsy2 , 0):
rounding_error_ts:=compute_horner_rounding_error(polyDoSin2,y,y2max, errlist, true):
delta_ts := rounding_error_ts[3]:
log2(delta_ts);




#################################################

# hi_lo takes an arbitrary precision number x and returns two doubles such that:
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






# Output

filename:="TEMPSIN/sine_fast.h":
fd:=fopen(filename, WRITE, TEXT):

fprintf(fd, "#include \"crlibm.h\"\n#include \"crlibm_private.h\"\n"):

fprintf(fd, "\n/*File generated by maple/coef_sine.mw*/\n"):

fprintf(fd, "#define XMAX_CODY_WAITE 0x%s\n", ieeehexa(xmax_cody_waite)[1]):
fprintf(fd, "#define XMAX_RETURN_X_FOR_SIN 0x%s\n", ieeehexa(xmax_return_x_for_sin)[1]):
fprintf(fd, "#define XMAX_SIN_FAST 0x%s\n", ieeehexa(xmaxSinFast)[1]):

fprintf(fd, "#define RND_CST_SINFAST %f \n", rnconstantSinFast);
#compute_rn_constant(2*rel_error_CodyWaite)):

fprintf(fd,"#ifdef WORDS_BIGENDIAN\n"):
for isbig from 1 to 0 by -1 do

  if(isbig=0) then
    fprintf(fd,"#else\n"):
  fi;

  fprintf(fd, "static db_number const pio256hi = "):
  printendian(fd, pio256hi, isbig):
  fprintf(fd, " ;\n"):

  fprintf(fd, "static db_number const mpio256lo = "):
  printendian(fd, -pio256lo, isbig):
  fprintf(fd, " ;\n"):

  fprintf(fd, "static db_number const mpio256med1 = "):
  printendian(fd, -pio256med1, isbig):
  fprintf(fd, " ;\n"):

  fprintf(fd, "static db_number const mpio256med2 = "):
  printendian(fd, -pio256med2, isbig):
  fprintf(fd, " ;\n"):

  fprintf(fd, "static db_number const mpio256lo2 = "):
  printendian(fd, -pio256lo2, isbig):
  fprintf(fd, " ;\n"):

  fprintf(fd, "static db_number const invpio256 = "):
  printendian(fd, invpio256, isbig):
  fprintf(fd, " ;\n\n"):

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

  fprintf(fd, "\n\n"):

  #The sincos table

  fprintf(fd, "\n/*  sine and cos of kPi/256 in double-double */\n"):
  SinCosSize:= 128;
  fprintf(fd, "static db_number const sincosTable[%d] =\n{\n",  4*(SinCosSize/2+1)):
  for i from 0 to SinCosSize/2 do
    if(1+1=2) then # normal tables
      s:=hi_lo(sin(i*Pi/(2*SinCosSize)));
      c:=hi_lo(cos(i*Pi/(2*SinCosSize)));
      printendian(fd,s[1],isbig);
      fprintf(fd," ,\n"):
      printendian(fd,s[2],isbig);
      fprintf(fd," ,\n"):
      printendian(fd,c[1],isbig);
      fprintf(fd," ,\n"):
      printendian(fd,c[2],isbig);
    else # tables with half-full MSdigit
      s:=hi_lo_sincos(sin(i*Pi/(2*SinCosSize)));
      c:=hi_lo_sincos(cos(i*Pi/(2*SinCosSize)));
      printendian(fd,s[1],isbig);
      fprintf(fd," ,\n"):
      printendian(fd,s[2],isbig);
      fprintf(fd," ,\n"):
      printendian(fd,c[1],isbig);
      fprintf(fd," ,\n"):
      printendian(fd,c[2],isbig);
    end if:

    if i<SinCosSize-1 then fprintf(fd," ,\n"): fi:
  od:
  fprintf(fd, "\n};\n\n"):


od:
fprintf(fd,"#endif /* WORDS_BIGENDIAN */\n\n\n"):

fclose(fd): 
