set terminal postscript eps dashed "Times-Roman" 24

set xlabel "x"
set ylabel "exp(x)"
set xrange [-5:5]
set yrange [-1:20]
set xzeroaxis lt 0 lw 0
set yzeroaxis lt 0 lw 0

set output "fig_exp.ps"
plot exp(x) notitle