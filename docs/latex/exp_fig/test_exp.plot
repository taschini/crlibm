set terminal postscript eps dashed "Times-Roman" 24

set xlabel "Nb d'appel"
set ylabel "Temps"
set logscale x


###
# Graphique de test sur l'exponentielle PENTIUM III
set output "fig_pentium_III.ps"
plot "res_pentiumIII.data" using 1:2 title "libm" with linespoint, "res_pentiumIII.data" using 1:3 title "crlibm" with linespoint lt 3, "res_pentiumIII.data" using 1:4 title "libultim" with linespoint lt 6 




###
# Graphique de test sur l'exponentielle XEON
set output "fig_xeon.ps"
plot "res_xeon.data" using 1:2 title "libm" with linespoint, "res_xeon.data" using 1:3 title "crlibm" with linespoint lt 3, "res_xeon.data" using 1:4 title "libultim" with linespoint lt 6 



###
# Graphique de test sur l'exponentielle ITANIUM
set output "fig_itanium.ps"
plot "res_itanium.data" using 1:2 title "libm" with linespoint, "res_itanium.data" using 1:3 title "crlibm" with linespoint lt 3, "res_itanium.data" using 1:4 title "libultim" with linespoint lt 6


###
# Graphique de test sur l'exponentielle SUN
set output "fig_sun.ps"
plot "res_sun.data" using 1:2 title "libm" with linespoint, "res_sun.data" using 1:3 title "crlibm" with linespoint lt 3, "res_sun.data" using 1:4 title "libultim" with linespoint lt 6 



###
# Graphique de test sur l'exponentielle PowerPC
set output "fig_ppc.ps"
plot "res_ppc.data" using 1:2 title "libm" with linespoint, "res_ppc.data" using 1:3 title "crlibm" with linespoint lt 3, "res_ppc.data" using 1:4 title "libultim" with linespoint lt 6 

