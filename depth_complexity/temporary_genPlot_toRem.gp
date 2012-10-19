set terminal pdf size 10, 6
set output "elephav.obj_20_r.pdf"
set style fill solid 1
set boxwidth 0.8 relative
set xtics 1
set ytics 1
set yrange [0:12]
set key off
set title "elephav.obj_20_r MDC Behavior"
set xlabel "Discretization Steps" 
set ylabel "Maximum Depth Complexity"
