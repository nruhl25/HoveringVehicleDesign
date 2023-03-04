set term gif animate delay 20
set output 'animate.gif'
vtip(n) = word("400.00 413.88 427.77 441.66 455.55 469.44 483.33 497.22 511.11 525.00",n)

stats 'animate.txt' name 'A'
set xrange [A_min_x:A_max_x]
set yrange [A_min_y:A_max_y]
set grid
set ylabel 'Inflow ratio'
set xlabel 'Freestream velocity'
set datafile separator whitespace
do for [i=1:int(A_blocks-1)] { plot 'animate.txt' index i w l title 'v_{tip}= '.vtip(i).' ft/sec'}