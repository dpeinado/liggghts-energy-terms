#!/bin/sh
gnuplot << EOF
set term png size 1600, 1400
set output "print_files/${2}.png"
set origin 0,0
set multiplot layout 2,2 scale 1,1  title "$1"

set style line 1 lt 1 lc rgb "magenta"      lw 2 pt 4 ps 0 pi 30
set style line 2 lt 1 lc rgb "blue"         lw 2 pt 1 ps 0 pi 30
set style line 3 lt 1 lc rgb "green"        lw 2 pt 2 ps 0 pi 30
set style line 4 lt 1 lc rgb "yellow"       lw 2 pt 3 ps 2 pi 30
set style line 5 lt 1 lc rgb "brown"        lw 2 pt 5 ps 0 pi 30
set style line 6 lt 1 lc rgb "black"        lw 2 pt 6 ps 1 pi 20
set style line 7 lt 1 lc rgb "purple"       lw 2 pt 7 ps 1 pi 30
set style line 8 lt 1 lc rgb "yellow"       lw 2 pt 8 ps 1 pi 30
set style line 9 lt 1 lc rgb "grey80"       lw 2 pt 9 ps 0 pi 30
set style line 10 lt 1 lc rgb "violet"      lw 2 pt 10 ps 0 pi 30
set style line 11 lt 1 lc rgb "black"       lw 2 pt 11 ps 1 pi 30
set style line 12 lt 1 lc rgb "red"         lw 2 pt 4  ps 1 pi 30
set grid 
plot "files/${2}0.dat" using 1:2 title "kE" with lp ls 1, "files/${2}0.dat" using 1:3 title "kR" with lp ls 2,"files/${2}0.dat" using 1:4 title "peN" with lp ls 3,"files/${2}0.dat" using 1:5 title "peT" with lp ls 4, "files/${2}0.dat" using 1:6 title "DeN" with lp ls 5,"files/${2}0.dat" using 1:7 title "DeTV" with lp ls 6,"files/${2}0.dat" using 1:8 title "DeTF" with lp ls 12,"files/${2}0.dat" using 1:9 title "TWF" with lp ls 7,"files/${2}0.dat" using 1:10 title "DEH" with lp ls 8, "files/${2}0.dat" using 1:11 title "ECons" with lp ls 9, "files/${2}0.dat" using 1:12 title "ETot" with lp ls 10
plot "files/${2}1.dat" using 1:2 title "kE" with lp ls 1, "files/${2}1.dat" using 1:3 title "kR" with lp ls 2,"files/${2}1.dat" using 1:4 title "peN" with lp ls 3,"files/${2}1.dat" using 1:5 title "peT" with lp ls 4, "files/${2}1.dat" using 1:6 title "DeN" with lp ls 5,"files/${2}1.dat" using 1:7 title "DeTV" with lp ls 6,"files/${2}1.dat" using 1:8 title "DeTF" with lp ls 12,"files/${2}1.dat" using 1:9 title "TWF" with lp ls 7,"files/${2}1.dat" using 1:10 title "DEH" with lp ls 8, "files/${2}1.dat" using 1:11 title "ECons" with lp ls 9, "files/${2}1.dat" using 1:12 title "ETot" with lp ls 10
plot "files/${2}2.dat" using 1:2 title "kE" with lp ls 1, "files/${2}2.dat" using 1:3 title "kR" with lp ls 2,"files/${2}2.dat" using 1:4 title "peN" with lp ls 3,"files/${2}2.dat" using 1:5 title "peT" with lp ls 4, "files/${2}2.dat" using 1:6 title "DeN" with lp ls 5,"files/${2}2.dat" using 1:7 title "DeTV" with lp ls 6,"files/${2}2.dat" using 1:8 title "DeTF" with lp ls 12,"files/${2}2.dat" using 1:9 title "TWF" with lp ls 7,"files/${2}2.dat" using 1:10 title "DEH" with lp ls 8, "files/${2}2.dat" using 1:11 title "ECons" with lp ls 9, "files/${2}2.dat" using 1:12 title "ETot" with lp ls 10
EOF
