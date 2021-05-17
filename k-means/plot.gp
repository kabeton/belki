stats "clusters.dat" u 1:2 nooutput
blocks = STATS_blocks
set grid x y
show grid

plot for [i=0:blocks-2] "clusters.dat" index i u 1:2 pt 7 ps 2 t columnheader(1)

pause mouse close
