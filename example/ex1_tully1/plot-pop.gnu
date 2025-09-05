#!/usr/bin/gnuplot
set encoding iso_8859_1
set key right top

set terminal pdf size 4in,3in enhanced color font 'Helvetica,16'
set border lw 2.5
set tics scale 1.2

set output 'fig-tully1-pop.pdf'
set key  samplen 1.0 spacing 1.3 font "Helvetica, 14"
set multiplot layout 1,1
set xlabel 'Time (a.u.)'
set ylabel 'Population'

set mxtics 2
set mytics 2
set ytics 0.2
set xtics 500
set yr[-0.01:1.01]

plot 'pop_diabat3.out' u 1:2 w l lw 6 t'Diabat State 1',\
    'pop_diabat3.out' u 1:3 w l lw 6 t'Diabat State 2'

