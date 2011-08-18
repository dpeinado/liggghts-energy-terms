#!/bin/bash
pitos="plot "
a=1
for name in files/plot_*.dat
do
    pitos="${pitos} \"${name}\" with lp ls $a,"
    a=`calc $a+1`
done
longitud=`calc "${#pitos}-1"`
gnuplot<<EOF
set term png size 1024,780
set output "${1}"
set grid
set key on left
set style line 1 lt 1 lc rgb "red" 	         lw 1 pt 1  ps 2 pi 1
set style line 2 lt 1 lc rgb "blue"          lw 1 pt 2  ps 2 pi 1
set style line 3 lt 1 lc rgb "green"         lw 1 pt 3  ps 2 pi 1
set style line 4 lt 1 lc rgb "black"         lw 1 pt 4  ps 2 pi 1
set style line 5 lt 1 lc rgb "yellow"        lw 1 pt 5  ps 2 pi 1
set style line 6 lt 1 lc rgb "red"  		     lw 1 pt 6  ps 2 pi 1
set style line 7 lt 1 lc rgb "blue"     	   lw 1 pt 7  ps 2 pi 1
set style line 8 lt 1 lc rgb "green"    	   lw 1 pt 8  ps 2 pi 1
set style line 9 lt 1 lc rgb "black"      	 lw 1 pt 9  ps 2 pi 1
set style line 10 lt 1 lc rgb "yellow"       lw 1 pt 10 ps 2 pi 1
${pitos:0:${longitud}}
pause -1
EOF
