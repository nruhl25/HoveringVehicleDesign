set term gif animate delay 4
set output 'animate.gif'

stats 'animate.txt' name 'A'
set xrange [A_min_x:A_max_x]
set yrange [A_min_y:A_max_y]
set datafile separator whitespace
do for [i=1:int(A_blocks-1)] { plot 'animate.txt' index i}