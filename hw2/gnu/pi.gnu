plot "inflow_alpha_n10.txt" u 1:($3*0.0018) w l lc 1 title "{/Symbol a}=-10 deg", "inflow_alpha0.txt" u 1:($3*0.0018) w l lc 2 title "{/Symbol a}=0 deg", "inflow_alpha10.txt" u 1:($3*0.0018) w l lc 3 title "{/Symbol a}=10 deg"
set grid
set ylabel "Induced power (hp)"
set xlabel "Free-stream Velocity (ft/sec)"