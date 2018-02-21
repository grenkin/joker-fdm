unset key
set term wxt 0
set ylabel "phi"
set xlabel "t"
plot "output_mid.txt" using 1:3 with lines, "output_mid.txt" using 1:4 with lines
set term wxt 1
set ylabel "theta"
set xlabel "t"
plot "output_mid.txt" using 1:2 with lines