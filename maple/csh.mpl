restart:
Digits := 150;
with (numapprox):with(orthopoly):
interface(quiet=true);
read "common-procedures.mpl";
mkdir("TEMPCSH");


#######################################################################
# Some values to go to the .testdata
# What is the first value that rounds to nearest to +inf ?
ieeehexa(arcsinh(hexa2ieee(["7fefffff","ffffffff"])));

# What is the first value that rounds to +inf to +inf ?


######################################################################
#First, some constants variables (used in the polynomial evaluations)
n_double_ch := 11: # max degree for cosh's polynomial evaluation
n_double_sh := 8:  # max degree for sinh's polynomial evaluation
b_max := 2**(-9.): # max absolute input for polynomial evaluation

######################################################################
#some constants ...
inv_ln_2 := 1/ln(2.):
ln2_hi := hexa2ieee(["3FE62E42", "FEFA3800"]):
ln2_lo := nearest(ln(2.)-ln2_hi):
two_43_44 := 2^43 + 2^44:
bias := convert(op(2,ieeehexa(ln(2.)/2.+two_43_44)),'decimal','hex'): #to get maximum index for the second range reduction ...


######################################################################
#Bounds of our evaluation
max_input_ch := arccosh(hexa2ieee(["7FEFFFFF","FFFFFFFF"])):
k_max_ch := ceil(max_input_ch / ln(2)):
max_input_sh := arcsinh(hexa2ieee(["7FEFFFFF","FFFFFFFF"])):
k_max_sh := ceil(max_input_sh / ln(2)):
#towards +inf, we have cosh(x) = sinh(x) since exp(-x) is too small
#so, we are sure of k_max_sh == k_max_ch
k_max := k_max_ch;

######################################################################
# When can we ignore exp(-x) in front of exp(x) in the first step ?
# We want the same error as in the general case
k_max_csh_approx_exp := 35:
tempxmax:=(k_max_csh_approx_exp-1)*log(2):
eps_csh_approx_exp := exp(-tempxmax)/exp(tempxmax):
log2(%);  #

# When can we ignore exp(-x) in front of exp(x) in the second step ?
# The worst case for exp for large arguments requires 115 bits
k_max_csh_approx_exp_2 := 65:
tempxmax:=(k_max_csh_approx_exp_2-1)*log(2):
eps_csh_approx_exp_2 := exp(-tempxmax)/exp(tempxmax):
log2(%); # 118


######################################################################
#The Taylor polynoms
poly_ch :=series(cosh(x),x,n_double_ch):
poly_ch := convert(poly_ch,polynom)-1;
poly_sh :=series(sinh(x),x,n_double_sh):
poly_sh := (convert(poly_sh,polynom))/x-1;






####################################################################
# secondary functions

size_of_table := convert(op(2,ieeehexa(two_43_44+ln(2.)/2.)),decimal,hex);
#returns the float which follow immediately the input
next_float := proc(value)
	local hex1,hex2,hexcat,result;
	hex1:= op(1, value):
	hex2:= op(2, value):
	hexcat:= cat(hex1, hex2);
	result := convert(convert(convert(hexcat,decimal,hex)+1+2**64,hex),string);
	result := [substring(result,2..9), substring(result,10..18)];
end:
#compute the errors done in tabulated values for cosh
delta_table_cosh_func := proc()
	local result, i, value, temp, tmp, maxi;
	value := ieeehexa(two_43_44-ln(2.)/2.);
	result := 0;
	maxi := 0;
	for i from -size_of_table to size_of_table do
		tmp := cosh(hexa2ieee(value)-two_43_44):
		temp := nearest(tmp):
		result := max(result, abs(tmp - temp - nearest(tmp-temp)));
		maxi := max(maxi, abs(temp + nearest(tmp-temp)));
		value:=next_float(value):
	od:
	result,maxi;
end:
#compute the error done in tabulated values for sinh
delta_table_sinh_func := proc()
	local result, i, value, temp, tmp, maxi;
	value := ieeehexa(two_43_44-ln(2.)/2.);
	result := 0;
	maxi := 0;
	for i from -size_of_table to size_of_table do
		tmp := sinh(hexa2ieee(value)-two_43_44):
		temp := nearest(tmp):
		result := max(result, abs(tmp - temp - nearest(tmp-temp)));
		maxi := max(maxi, abs(temp + nearest(tmp-temp)));
		value:=next_float(value):
	od:
	result,maxi;
end:


#return the error on x_hi * y_lo
Mul11_Error := proc(x,err_x, y, err_y)
(2^(-53) * y + err_y) * (err_x + 1/2*ulp(x)) + x * err_y;
end:


#return the error on (x_hi * y_hi)_lo + x_lo * y_hi + x_hi * y_lo
Mul43_Error := proc(x,err_x, y, err_y)
1/2*ulp(3*2^(-53)* x * y) + 1/2*ulp(2*2^(-53) * x * y) + Mul11_Error(x, err_x, y, err_y) + Mul11_Error(y, err_y, x, err_x);
end:



cosh_0_35 := proc()
local k, delta, maxi,epsilon, mini,cosh_0_35_max, delta_cosh_0_35:
epsilon := 0; delta := 0:
for k from -35 to -1 do
delta_cosh_0_35 := 1/2*2^(-53)*ulp(1/(2^k)*(cosh_max + sinh_max)) + 1/(2^k)*(delta_sinh + delta_cosh):
cosh_0_35_max := 1/(2^k)*(cosh_max + sinh_max):
delta_cosh_0_35 :=  1/2*2^(-53)*ulp(cosh_0_35_max + 2^k * sinh_max) + delta_cosh_0_35 + 2^k  * delta_sinh:
cosh_0_35_max := cosh_0_35_max + 2^k * sinh_max:
delta_cosh_0_35 :=  1/2*2^(-53)*ulp(cosh_0_35_max + 2^k * cosh_max) + delta_cosh_0_35 + 2^k  * delta_cosh:
cosh_0_35_max := cosh_0_35_max + 2^k * cosh_max:
maxi := max(evalf(cosh((k-1/2)*ln(2))),evalf(cosh((k+1/2)*ln(2)))):
mini := min(evalf(cosh((k-1/2)*ln(2))),evalf(cosh((k+1/2)*ln(2)))):

epsilon := max(epsilon, delta_cosh_0_35/mini):
delta := max(delta, delta_cosh_0_35):
od;
for k from 1 to 35 do
delta_cosh_0_35 := 1/2*2^(-53)*ulp(1/(2^k)*(cosh_max + sinh_max)) + 1/(2^k)*(delta_sinh + delta_cosh):
cosh_0_35_max := 1/(2^k)*(cosh_max + sinh_max):
delta_cosh_0_35 :=  1/2*2^(-53)*ulp(cosh_0_35_max + 2^k * sinh_max) + delta_cosh_0_35 + 2^k  * delta_sinh:
cosh_0_35_max := cosh_0_35_max + 2^k * sinh_max:
delta_cosh_0_35 :=  1/2*2^(-53)*ulp(cosh_0_35_max + 2^k * cosh_max) + delta_cosh_0_35 + 2^k  * delta_cosh:
cosh_0_35_max := cosh_0_35_max + 2^k * cosh_max:

maxi := max(evalf(cosh((k-1/2)*ln(2))),evalf(cosh((k+1/2)*ln(2)))):
mini := min(evalf(cosh((k-1/2)*ln(2))),evalf(cosh((k+1/2)*ln(2)))):
epsilon := max(epsilon, delta_cosh_0_35/mini):
delta := max(delta, delta_cosh_0_35):
od;
delta, epsilon;
end:


cosh_35_inf := proc()
local k, delta, maxi,epsilon, mini,cosh_0_35_max, delta_cosh_0_35:
epsilon := 0; delta := 0:
for k from 35 to 1025 do
delta_cosh_0_35 := 1/2*2^(-53)*ulp(2^k*(cosh_max + sinh_max)) + 2^k*(delta_sinh + delta_cosh) + 1/(2^k)*(cosh_max + sinh_max + delta_sinh + delta_cosh):
cosh_0_35_max := 2^k*(cosh_max + sinh_max):

maxi := max(evalf(cosh((k-1/2)*ln(2))),evalf(cosh((k+1/2)*ln(2)))):
mini := min(evalf(cosh((k-1/2)*ln(2))),evalf(cosh((k+1/2)*ln(2)))):
epsilon := max(epsilon, delta_cosh_0_35/mini):
delta := max(delta, delta_cosh_0_35):
od;
delta, epsilon;
end:
#############################################################








######################################################################
#now we can begin the proof

######################################################################
#first, we must compute the error created by the first range reduction
# CODY and WAITE  Argument reduction

ln2 := ln(2.);
invln2:= nearest(1/ln2);
reminvln2 := evalf(1/ln2 - invln2);
expln2:=ieeedouble(ln2)[2]:   #get the exponent of ln2 in its IEEE representation

bits_ln2_hi_0 := ceil(log2(k_max));
# 1/2 <= ln2/2^(expln2 + 1) <  1
ln2_hi := round(evalf(ln2 * 2^(52 - bits_ln2_hi_0 - expln2)))/2^(52 - bits_ln2_hi_0 - expln2);#this trick is used to get a truncated mantissa for ln2_hi
#ln2_hi is a now exactly representable in the IEEE format and bits_ln2_hi_0 last bits of its mantissa are set to 0
#and bits_ln2_hi_0 is exactly the max number of bits to represent k
#so the k * ln2_hi-product is exact :)
ln2_lo:=nearest(ln2 - ln2_hi):

# The error in this case (we need absolute error)
delta_repr_ln2 := abs(ln2 - ln2_hi - ln2_lo);
delta_round := evalf(1/2 * ulp(ln2_lo));
delta_cody_waite := k_max * (delta_repr_ln2 + delta_round);

delta_b := delta_repr_ln2 + delta_round + delta_cody_waite;

#we have 2 cases:
#  * k != 0 and we have an inexact range reduction
#  * k == 0 and we have no range reduction and then delta_range_reduc = 0

#the second range reduction is exact, so it doesn't introduce new error
#after this second range reduction, we have a argument <= 2^(-9)
#'mathematical' reductions:
# x = k * ln(2) + y
# y = a + b
#'true' reductions:
# x = k * (ln2_hi + ln2_lo) + (b_hi + b_lo) + table_index_float
#   with table_index_float = table_index * 2^(-8)

#we'll use the following mathematical formulaes :
#cosh(a + b) = cosh(a) * cosh(b) + sinh(a) * sinh(b)
#sinh(a + b) = sinh(a) * cosh(b) + sinh(b) * cosh(a)
#sinh(a) and cosh(a) are tabulated as double double
# and we use Taylor series to compute approximations of the sums

#computation of the absolute error in the tabulated values :
delta_ca := delta_table_cosh_func()[1];
delta_sa := delta_table_sinh_func()[1];
ca_max := delta_table_cosh_func()[2];
sa_max := delta_table_sinh_func()[2];


#now we must compute the error done in polynomial evaluation
#we use   cosh(b) = 1 + sum(b^(2*k)/(2*k!), k > 0)
#         sinh(b) = b * (1 + sum(b^(2*k)/(2*k+1!), k > 0))

#both used polynoms are even (x², x^4, x^6, ...) and we can use this fact
y_max := b_max ^ 2;
delta_y := 1/2*ulp(y_max) + delta_b^2;
# remove the first x and compute the polynomial of y = x²
poly_ch2 :=  subs(x=sqrt(y), expand(poly_ch));
poly_sh2 :=  subs(x=sqrt(y), expand(poly_sh));
errlist_cosh := errlist_quickphase_horner(degree(poly_ch2), 0, 0, 2^(-53), 2^(-70)):
errlist_sinh := errlist_quickphase_horner(degree(poly_sh2), 0, 0, 2^(-53), 2^(-70)):

#error between effective result and theorical polynomial result
rounding_error_tcb := compute_horner_rounding_error(poly_ch2, y, y_max, errlist_cosh, true);
rounding_error_tsb := compute_horner_rounding_error(poly_sh2, y, y_max, errlist_sinh, true);

#error between therical polynomial result and cosh value
approx_error_tcb := infnorm((cosh(x)-1-poly_ch)/(cosh(x)-1),x= -b_max..b_max);
approx_error_tsb := infnorm((sinh(x)/x-1-poly_sh)/(sinh(x)/x-1),x= -b_max..b_max);

delta_tcb := rounding_error_tcb[2] + approx_error_tcb;
delta_tsb := rounding_error_tsb[2] + approx_error_tsb;
tsb_max := rounding_error_tsb[4];
tcb_max := rounding_error_tcb[4];

#now we must do the first reconstruction, which correspond to cosh(a + b) = cosh(a) * cosh(b) + sinh(a) * sinh(b)
#first case : sinh(a) = 0 = sa (= sa_hi + sa_lo in the C code)
#             cosh(a) = 1 = ca (= ca_hi + ca_lo in the C code)
#there is no error on sa and caj
delta_cosh0 := delta_tcb;
cosh0_max := 1+tcb_max;
delta_sinh0 := delta_b:
sinh0_max := b_max:
delta_sinh0 := delta_sinh0 + (tsb_max + delta_tsb)*(b_max + delta_b + 1/2*ulp(b_max)) - b_max * tsb_max + 1/2*ulp(sinh0_max + tsb_max*b_max):
sinh0_max := sinh0_max + tsb_max*b_max:
delta_sinh0 := delta_sinh0 + 2^(-53)*1/2*ulp(sinh0_max + tsb_max*b_max);
sinh0_max := sinh0_max + tsb_max*b_max;

#second case : sinh(a) <> 0
#there is a delta_table_cosh and delta_table_sinh absolute error on ca and sa.

delta_cosh1 := delta_ca:
cosh1_max := 2^(-53)*ca_max:
delta_cosh1 := 1/2*ulp(cosh1_max + 3*2^(-53)*b_max * sa_max) + delta_cosh1 + Mul43_Error(b_max, delta_b, sa_max, delta_sa):
cosh1_max := cosh1_max + 3 * 2^(-53) * b_max * sa_max:
delta_cosh1 := 1/2*ulp(cosh1_max + b_max * sa_max * tsb_max) + delta_cosh1 + ((sa_max+1/2*ulp(sa_max)+delta_sa)*(b_max+ulp(b_max)+delta_b)*(tsb_max+delta_tsb)-tsb_max*sa_max*b_max):
cosh1_max := cosh1_max + b_max * sa_max * tsb_max:
delta_cosh1 := 1/2*ulp(cosh1_max + tcb_max*ca_max) + delta_cosh1 + ((ca_max + 1/2*ulp(ca_max) + delta_ca)*(tcb_max-delta_tcb)-ca_max*tcb_max):
cosh11_max := cosh1_max + tcb_max*ca_max:
delta_cosh1 := 1/2*ulp(cosh1_max + b_max * sa_max) + delta_cosh1:
cosh11_max := cosh1_max + b_max * sa_max:
delta_cosh1 := 1/2*2^(-53)*ulp(cosh1_max + ca_max) + delta_cosh1;
cosh11_max := cosh1_max + ca_max;

cosh_max := max(cosh0_max, cosh1_max);
delta_cosh := max(delta_cosh0, delta_cosh1);

delta_sinh1 := delta_sa:
sinh1_max := 2^(-53)*sa_max:
delta_sinh1 := 1/2*ulp(sinh1_max + 3*2^(-53)*ca_max * b_max) + delta_sinh1 + Mul43_Error(ca_max, delta_ca, b_max, delta_b):
sinh1_max := sinh1_max + 3*2^(-53)*ca_max * b_max:
delta_sinh1 := delta_sinh1 + 1/2*ulp(sinh1_max + sa_max * tcb_max):
sinh1_max := sinh1_max + sa_max * tcb_max:
delta_sinh1 := 1/2*ulp(sinh1_max + b_max * ca_max * tsb_max) + delta_sinh1 + ((ca_max+1/2*ulp(ca_max)+delta_ca)*(b_max+ulp(b_max)+delta_b)*(tsb_max+delta_tsb)-tsb_max*ca_max*b_max):
sinh1_max := sinh1_max + b_max * sa_max * tsb_max:
delta_sinh1 := delta_sinh1 + 1/2*2^(-53)*ulp(sinh1_max + b_max * ca_max):
sinh1_max := sinh1_max + b_max * ca_max:
delta_sinh1 := delta_sinh1 + 2^(-53)*1/2*ulp(sinh1_max + sa_max);
sinh1_max := sinh1_max + sa_max;

sinh_max := max(sinh1_max, sinh0_max);
delta_sinh := max(delta_sinh0, delta_sinh1);
#so we have the error done on cosh(y) and sinh(y)
#now we must compute the error done on the last reconstruction
#there are many cases
#we begin by 0 < |k| < 35
epsilon_cosh_0_35 := cosh_0_35()[2];
#|k| > 35
epsilon_cosh_35_inf := cosh_35_inf()[2];

#rounding constant
maxepsilon_csh := max(epsilon_cosh_35_inf, epsilon_cosh_0_35, delta_cosh):
round_cst_csh := evalf(compute_rn_constant(maxepsilon_csh));;



#################################################################################################"
#now some functions used to build the .h header file.
#################################################################################################"
IEEE2db_number_BE := proc(ieee_number)
local hex1, hex2, hexcat;
hexcat=ieeehexa(ieee_number);
hex1:= op(1, ieeehexa(ieee_number)):
hex2:= op(2, ieeehexa(ieee_number)):
cat(cat("{{0x"||hex1||",0x"||hex2||"}};     /*",sprintf("%.10e",ieee_number)),"*/  \n");
end:
IEEE2db_number_LE := proc(ieee_number)
local hex1, hex2, hexcat;
hexcat=ieeehexa(ieee_number);
hex1:= op(2, ieeehexa(ieee_number)):
hex2:= op(1, ieeehexa(ieee_number)):
cat(cat("{{0x"||hex1||",0x"||hex2||"}};     /*",sprintf("%.10e",ieee_number)),"*/  \n");
end:
IEEE2db_number := proc(ieee_number, big_little)
if (big_little = 1) then
IEEE2db_number_BE(ieee_number):
else
IEEE2db_number_LE(ieee_number):
fi;
end:
lo_part := proc(ieee_number)
if (ieee_number <> 0) then
ieee_number - nearest(ieee_number);
else
0;
end if;
end:
IEEE2db_db_number_sinh := proc(number,big_little)
local hexstring1,hexstring2, hex1,hex2,hex3,hex4;
hexstring1 := ieeehexa(number);
hexstring2 := ieeehexa(lo_part(number));
if (big_little = 1) then
hex1:= op(1, hexstring1):
hex2:= op(2, hexstring1):
hex3:= op(1, hexstring2):
hex4:= op(2, hexstring2):
else
hex1:= op(2, hexstring1):
hex2:= op(1, hexstring1):
hex3:= op(2, hexstring2):
hex4:= op(1, hexstring2):
fi;
"	  {{0x"||hex1||", 0x"||hex2||"}}, {{0x"||hex3||", 0x"||hex4||"}}},\n";
end:

IEEE2db_db_number_cosh := proc(number,big_little)
local hexstring1,hexstring2, hex1,hex2,hex3,hex4;
hexstring1 := ieeehexa(number);
hexstring2 := ieeehexa(lo_part(number));
if (big_little = 1) then
hex1:= op(1, hexstring1):
hex2:= op(2, hexstring1):
hex3:= op(1, hexstring2):
hex4:= op(2, hexstring2):
else
hex1:= op(2, hexstring1):
hex2:= op(1, hexstring1):
hex3:= op(2, hexstring2):
hex4:= op(1, hexstring2):
fi;
"	{{{0x"||hex1||", 0x"||hex2||"}}, {{0x"||hex3||", 0x"||hex4||"}},\n";
end:

#####################################################################################



#now, we can produce the header file !
round_cst_cosh := 1.0020:
round_cst_sinh := round_cst_cosh:
filename := "TEMPCSH/csh_fast.h":
fd := fopen(filename, WRITE, TEXT):
fprintf(fd, "\n /* File generated by maple/csh.mpl */ \n"):
fprintf(fd, "\n"):
fprintf(fd, " static double maxepsilon_csh = %1.30e ;\n", maxepsilon_csh):
fprintf(fd, " static double round_cst_csh  = %1.30e ;\n", round_cst_csh):
fprintf(fd, "\n"):


fprintf(fd, "#ifdef WORDS_BIGENDIAN  \n"):

for big_little from 1 to 2 do
	if (big_little = 2) then
		fprintf(fd,                "#else  \n"):
	fi:
	fprintf(fd, cat(  "  static db_number const inv_ln_2 =     ",IEEE2db_number(inv_ln_2,big_little))):
	fprintf(fd, cat(  "  static db_number const ln2_hi =       ",IEEE2db_number(ln2_hi,big_little))):
	fprintf(fd, cat(  "  static db_number const ln2_lo =       ",IEEE2db_number(ln2_lo,big_little))):
	fprintf(fd, cat(  "  static db_number const two_43_44 =    ", IEEE2db_number(two_43_44,big_little))):
	fprintf(fd, cat(  "  static db_number const two_minus_30 = ", IEEE2db_number(2**(-40),big_little))):
	fprintf(fd,       "  static int const bias = %d ;\n", bias):
	fprintf(fd, "\n");

	fprintf(fd,"/* some bounds */ \n"):
	fprintf(fd, cat(  "  static db_number const max_input_csh =  ",IEEE2db_number(max_input_ch,big_little))):

    fprintf(fd, "\n"):
	fprintf(fd, cat(cat("  static const db_number cosh_sinh_table[",convert(2*size_of_table+1,string)),"][4] = { \n"));
	vvalue := ieeehexa(two_43_44-ln(2.)/2.);
	for i from -size_of_table to size_of_table do
		fprintf(fd, IEEE2db_db_number_cosh(cosh(hexa2ieee(vvalue)-two_43_44),big_little));
		fprintf(fd, IEEE2db_db_number_sinh(sinh(hexa2ieee(vvalue)-two_43_44),big_little));
		vvalue:=next_float(vvalue):
	od:
	fprintf(fd,"}; \n");

	fprintf(fd,"/* the coefficients for the cosh-approximations */ \n"):
	for i from 1 to (n_double_ch/2) do
		fprintf(fd, cat(cat(cat(  "  static const db_number c",convert(2*(i-1),string))," =   "),IEEE2db_number(coeff(poly_ch+1,x,2*(i-1)),big_little))):
	od:

	fprintf(fd,"/* the coefficients for the sinh-approximations */\n"):
	for i from 1 to (n_double_sh/2) do
		fprintf(fd, cat(cat(cat(  "  static  const db_number s",convert(2*i-1,string))," =   "),IEEE2db_number(coeff(poly_sh+1,x,2*(i-1)),big_little))):
	od:
od:
fprintf(fd, "#endif  \n"):

fclose(fd):
