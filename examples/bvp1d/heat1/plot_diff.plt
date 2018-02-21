unset key
set term wxt 0
set logscale y
set ylabel "r(t)"
set xlabel "t"
plot "output_diff.txt" using 1:3 with lines