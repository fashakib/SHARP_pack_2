#!/usr/bin/gnuplot
set encoding iso_8859_1
set key right top

set terminal pdf size 4in,3in enhanced color font 'Helvetica,16'
set border lw 2.5
set tics scale 1.2

set output 'fig-tully3-branch.pdf'
set key  samplen 1.0 spacing 1.3 font "Helvetica, 14"
set multiplot layout 1,1
set xlabel 'k (a.u.)'
set ylabel 'P(k)'

set mytics 2
set ytics 0.2
set xtics 5
set yr[0:1]
ps=0.7

plot 'pop_branch_all.out' u 1:2 w lp lc 2 lw 2 dt 2 pt 7 ps ps t'T1',\
     'pop_branch_all.out' u 1:3 w lp lc 6 lw 2 dt 2 pt 7 ps ps t'T2',\
     'pop_branch_all.out' u 1:4 w lp lc 7 lw 2 dt 2 pt 7 ps ps  t'R1',\
     'pop_branch_all.out' u 1:5 w lp lc 4 lw 2 dt 2 pt 7 ps ps t'R2'

