#!/bin/sh
gnuplot << EOF
set term png size 1600, 1300
set output "print_files/${2}.png"
set origin 0,0
set title "$1"
set grid
set style line 1 lt 1 lc rgb "magenta"      lw 1 pt 4 ps 0 pi 10
set style line 2 lt 1 lc rgb "blue"         lw 1 pt 1 ps 0 pi 10
set style line 3 lt 1 lc rgb "green"        lw 1 pt 2 ps 0 pi 10
set style line 4 lt 1 lc rgb "yellow"       lw 1 pt 3 ps 0 pi 10
set style line 5 lt 1 lc rgb "brown"        lw 1 pt 5 ps 0 pi 10
set style line 6 lt 1 lc rgb "black"        lw 1 pt 6 ps 0 pi 10
set style line 7 lt 1 lc rgb "purple"       lw 1 pt 7 ps 0 pi 10
set style line 8 lt 1 lc rgb "red"       lw 1 pt 8 ps 0 pi 10
set style line 9 lt 1 lc rgb "grey80"       lw 1 pt 9 ps 0 pi 10
set style line 10 lt 1 lc rgb "violet"      lw 1 pt 10 ps 0 pi 10
set style line 11 lt 1 lc rgb "black"       lw 1 pt 11 ps 1 pi 10
set style line 12 lt 1 lc rgb "yellow"         lw 1 pt 4  ps 1 pi 10
set y2tics
set y2label "In %"
set ylabel "In Energy units"
#set y2range [0:0.014]
plot "${2}.dat" using 1:2 title "eKin" with lp ls 1, "${2}.dat" using 1:3 title "eCOL" with lp ls 2, "${2}.dat" using 1:4 title "ePG" with lp ls 3,"${2}.dat" using 1:5 title "DEH" with lp ls 4, "${2}.dat" using 1:6 title "IKE" with lp ls 5, "${2}.dat" using 1:7 title "LHS" with lp ls 6, "${2}.dat" using 1:(1*\$8) title "errA" with lp ls 7, "${2}.dat" using 1:9 axes x1y2 title "errR %" with lp ls 8
#plot [0:200000] [0:35] "${2}.dat" using 1:2 title "eKin" with lp ls 1, "${2}.dat" using 1:3 title "eCOL" with lp ls 2, "${2}.dat" using 1:4 title "ePG" with lp ls 3,"${2}.dat" using 1:5 title "DEH" with lp ls 4, "${2}.dat" using 1:6 title "IKE" with lp ls 5, "${2}.dat" using 1:7 title "LHS" with lp ls 6, "${2}.dat" using 1:8 axes x1y2 title "errA" with lp ls 7, "${2}.dat" using 1:9 axes x1y2 title "errR %" with lp ls 8
EOF
