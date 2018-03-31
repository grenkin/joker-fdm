set term postscript eps enhanced font ', 16' 
set output "monot_theta.eps"
set ylabel "{/Symbol-Oblique q}" norotate font ',20' offset 1.5,0
set xlabel '{/Helvetica-Italic x}' font ',20' offset 0,0.5
set key bottom right spacing 1.2
plot "output_theta_phi_min.txt" using 1:2 with lines lw 2 lt 2 title "min", \
  "output_theta_phi_max.txt" using 1:2 with lines lw 2 lt 1 title "max"