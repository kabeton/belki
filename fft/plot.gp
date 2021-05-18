set term png size 1300,600
set output "gr.png"
set multiplot layout 1, 2
set grid x y
show grid
plot "rfft.dat" u 1:2 t 'transformed' w lines, "rfft.dat" u 1:3 t 'initial' w lines
set grid x y
show grid
set xrange [-5:5]
plot "fft.dat" u 1:2 w lines t 'before filter', "fftt.dat" u 1:2 w lines t 'after filter'
unset multiplot
pause mouse close

