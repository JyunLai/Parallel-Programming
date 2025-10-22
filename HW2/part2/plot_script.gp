set terminal pngcairo size 800,600 font "Verdana,12"

set output 'speedup_graph.png'

set title "Speedup vs. Number of Threads (Block Decomposition)"

set xlabel "Number of Threads"
set ylabel "Speedup"

set xrange [0:7]
set yrange [0:7]

set grid

plot "speedup_data.txt" using 1:2 with linespoints title "Actual Speedup"
