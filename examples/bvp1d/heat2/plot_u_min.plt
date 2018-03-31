set term postscript eps enhanced font ', 18' 
set output "monot_u_min.eps"
set ylabel "{/Helvetica-Italic u}({/Helvetica-Italic L}, {/Helvetica-Italic t})" norotate font ',20'
set xlabel '{/Helvetica-Italic t}' font ',20'
set yrange [0:0.51]
unset key
plot "output_min.txt" using 1:3 with lines linewidth 3