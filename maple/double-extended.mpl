
#---------------------------------------------------------------------
# ieeedouble converts a number to IEEE double extended format.
# returns sign (-1 or 1), exponent between -16383 and 16383, mantissa as a fraction between 0.5 and 1.

# TODO : use JMM procedure; check subnormals etc
ieeedoubleExt:=proc(xx)
local x, sign, logabsx, exponent, mantissa, infmantissa;
    x:=evalf(xx):
    if (x=0) then
        sign,exponent,mantissa := 0,0,0;
    else
        if (x<0) then sign:=-1:
        else sign:=1:
        fi:
        exponent := floor(log2(sign*x));
        if (exponent>16383) then mantissa:=infinity: exponent:=16383:
        elif (exponent< -16382) then
            # denorm
            exponent := -16383
        fi:
        infmantissa := sign*x*2^(63-exponent);
        if frac(infmantissa) <> 0.5 then mantissa := round(infmantissa)
        else
            mantissa := floor(infmantissa);
            if type(mantissa,odd) then mantissa := mantissa+1 fi;
        fi;
        mantissa := mantissa*2^(-63);
    fi;
    sign,exponent,mantissa;
end:

nearestExt := proc(x)
    local sign, exponent, mantissa;

    sign, exponent, mantissa := ieeedoubleExt(x);
    sign*mantissa*2^(exponent);
end:



ieeehexaExt:= proc(x)
local resultat, sign, exponent, mantissa, t;

    if(x=0) then resultat:=["0000","0000","0000","0000","0000"];
    elif(x=-0) then resultat:=["8000","0000","0000","0000","0000"];
    else
        sign,exponent,mantissa := ieeedoubleExt(x);
        t := 2**80 +  (exponent+16383)*2^64 + mantissa*2^63;
        if (sign=-1) then
            t := t + 2**79;
        fi:
   t := convert(t, hex);
   t:=convert(t, string):

   resultat:=[substring(t, 2..5),
              substring(t, 6..9 ),  substring(t, 10..13 ),
              substring(t, 14..17 ),  substring(t, 18..21 )];

    end if:
    resultat;
end proc:



printDoubleAsShort:=proc(x)
    local ss;
    ss:=ieeehexa(x);
    cat( "DOUBLE_HEX(",
         substring(ss[1], 1..4),  ", " ,
         substring(ss[1], 5..8),  ", " ,
         substring(ss[2], 1..4),  ", " ,
         substring(ss[2], 5..8)) ;
end proc:

printDoubleExtAsShort:=proc(x)
    local ss;
    ss:=ieeehexaExt(x);
    cat( "LDOUBLE_HEX(", ss[1],  ", ", ss[2],  ", ",ss[3],  ", ", ss[4],  ", ", ss[5],  ")");
end proc:

printDoubleAsULL:=proc(x)
    local ss;
    ss:=ieeehexa(x);
    cat( "ULL(", ss[1], ss[2],  ")");
end proc:

printDoubleAsHexInt:=proc(x)
    local ss;
    ss:=ieeehexa(x);
    cat(ss[1], ss[2]);
end proc:




#---------------------------------------------------------------------
# hi_lo takes an arbitrary precision number x and returns two doubles such that:
# x ~ x_hi + x_lo
hiloExt:= proc(x)
local x_hi, x_lo, res:
    x_hi:= nearestExt(evalf(x)):
    res:=x-x_hi:
    if (res = 0) then
        x_lo:=0:
    else
        x_lo:=nearestExt(evalf(res)):
    end if;
    x_hi,x_lo;
end:


#---------------------------------------------------------------------
# Like poly_exact, but the n first coefficients are exactly representable as the sum of two doubles.
# (to actually get the two doubles, use procedure hi_lo)

polyExact2Ext:=proc(P,n)
local deg,i, coef, coef_hi, coef_lo, Q:
    Q:= 0:
    convert(Q, polynom):
    deg:=degree(P,x):
    for i from 0 to deg do
        coef :=coeff(P,x,i):
        coef_hi, coef_lo:=hiloExt(coef):
        Q:= Q + coef_hi*x^i:
        if(i<n) then
            Q := Q + coef_lo*x^i:
        fi:
    od:
    return(Q);
end:

printPolyExt := proc(fd,P,n, name_of_poly)
local deg,i, coef, coef_hi, coef_lo;
    convert(Q, polynom):
    deg:=degree(P,x):
    fprintf(fd, " static const long double %s[%d][2] = {\n", name_of_poly, deg+1);
    for i from 0 to deg do
        coef :=coeff(P,x,i):
        coef_hi, coef_lo:=hiloExt(coef):

        fprintf(fd,"{ %1.50eL, ",coef_hi);

        if(i<n) then
            fprintf(fd," %1.50eL},\n",coef_lo);
        else
            fprintf(fd,"0},\n");
        fi:
    od:
    fprintf(fd,"}; \n");
end:
