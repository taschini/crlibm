
#####################################################################

# Stuff for SCS
# All the fllowing may be outdated

# Global parameters
# Don´t forget to set all the parameters here
SCS_NB_WORDS := 8:
SCS_NB_BITS  := 30:

# GetSCS_real
# This procedure convert a decimal number into it SCS representation.
#        x : input number to convert into it SCS representation
real_to_SCS := proc(x)
        local exception, index, sign, mantissa, nb, i;

            if x <> 0 then
                exception := 1;
                if x > 0 then
                    sign  := 1;
                    nb    := x;
                elif x < 0 then
                    sign := -1;
                    nb   := -x;
                end if;

                index := 0;

                if nb >= 1 then
                    for i from 0 while nb > (2^(SCS_NB_BITS+1)-1) do
                        index := index+1;
                        nb    := nb * 2^(-SCS_NB_BITS);
                    end do;
                else
                    for i from 0 while nb < 1 do
                        index := index-1;
                        nb    := nb * 2^(SCS_NB_BITS);
                    end do;
                end if;

                for i from 0 by 1 to (SCS_NB_WORDS-1) do
                    mantissa[i] := trunc(nb);
                    nb          := (nb - mantissa[i]) * 2^(SCS_NB_BITS);
                end do;
            else
                for i from 0 by 1 to (SCS_NB_WORDS-1) do
                    mantissa[i] := 0;
                end do;

                index     := 1;
                exception := x;
                sign      := 1;
            end if;
            mantissa[SCS_NB_WORDS]   := exception;
            mantissa[SCS_NB_WORDS+1] := index;
            mantissa[SCS_NB_WORDS+2] := sign;

            return mantissa;
        end proc:




# SetSCS_real
# Convert an SCS number into a rational number

     SCS_to_real := proc(tab)
            local res, i;

                if (tab[SCS_NB_WORDS] <> 1) then
                    return tab[SCS_NB_WORDS];
                end if;

                res := 0;
                for i from (SCS_NB_WORDS-1) by -1 while i>=0 do
                    res := 2^(-SCS_NB_BITS)*res + tab[i]
                end do;

                res := tab[SCS_NB_WORDS+2]*(res * 2.^(SCS_NB_BITS * tab[SCS_NB_WORDS+1]));

                return res;

            end proc:


# WriteSCS
# Write Into file fd the SCSS number stored into the table tab where
# tab[0..(SCS_NB_WORDS-1)] store the mantissa
# tab[SCS_NB_WORDS] store the exception
# tab[SCS_NB_WORDS+1] store the index
# tab[SCS_NB_WORDS+2] store the sign


         WriteSCS_from_table := proc(fd, tab)
                       local i;

                           fprintf(fd,"{{");

                           fprintf(fd,"0x%+0.8x, ", tab[0]);
                           for i from 1 by 1 to (SCS_NB_WORDS-2) do
                               fprintf(fd,"0x%+0.8x, ", tab[i]);
                               if (i mod 4 = 3) then
                                   fprintf(fd,"\n");
                               fi;
                           end do;
                           fprintf(fd,"0x%+0.8x},\n", tab[SCS_NB_WORDS-1]);
                           if (tab[SCS_NB_WORDS]=1) then
                               fprintf(fd,"DB_ONE, %3d, %3d ", tab[SCS_NB_WORDS+1], tab[SCS_NB_WORDS+2]);
                           else
                               # the only other possible value is 0 so ...
                               fprintf(fd,"{0x00000000, 0x00000000}, %3d, %3d ", tab[SCS_NB_WORDS+1], tab[SCS_NB_WORDS+2]);
                           end if;

                           fprintf(fd, "} \n");
                       end proc:
# GetSCS_poly

                  get_nb_terms := proc(poly)
                         local i, deg_poly;

                             deg_poly := degree(poly);
                             for i from deg_poly by -1 while i>=0 do
                                 if coeff(poly, x, i)=0 then
                                     deg_poly := deg_poly-1;
                                 end if;
                             end do;

                             return deg_poly;
                         end proc:

                  WriteSCS := proc (fd,x)
                      WriteSCS_from_table(fd , real_to_SCS (x));
end:


# Convert each coefficient of a polynom into it SCSS representation


#poly : input polynom
#file : name of the file where to put the result
Write_SCS_poly := proc(poly, fd)
local i, deg;
           #fclose(fd);
           try
           finally
               fprintf(fd,"static const scs constant_poly[%d]=\n",get_nb_terms(poly)+1);
               deg := degree(poly);

               fprintf(fd,"/* ~%1.50e */ \n{", coeff(poly, x, deg));
               WriteSCS(fd, coeff(poly, x, deg));
               for i from (deg-1) by (-1) while i>=0 do
                   if (coeff(poly, x, i)<>0) then
                       fprintf(fd,",\n/* ~%1.50e */ \n", coeff(poly, x, i));
                       WriteSCS(fd, coeff(poly, x, i), 0);
                   end if;
               end do;
               fprintf(fd,"};\n");
            end try;
       end proc:


# hexa_scs

# Une dernière procédure pour écrire les nombres SCS (floatant exacts) en hexadécimal dans un fichier ".h"

# hexa_scs := proc (x)
#   local i, elem, mantisse, exposant, resultat;
#   global ulp, ulp_inv;
#   if (x :: list) then
#     [seq (hexa(elem), elem=x)];
#   else
#     Digits := 60;
#     if (x = 0) then
#       resultat := 0;
#     else
#       exposant := floor(log[2](abs(evalf(x))));
#       if (exposant < -126) then exposant := -126; fi;
#       resultat := round (ulp_inv * (abs(x) * 2^(-exposant) -1)) + ulp_inv / 2 * (exposant + 2^(lgex-1) - 1);
#       if (x < 0) then resultat := resultat + ulp_inv * 2^(lgex-1); fi;
#     fi;
#   fi;
#   convert(resultat, hex);
# end:
