#!/usr/bin/gnuplot
set encoding iso_8859_1
set key right center
set terminal pdf size 4in,2.8in enhanced color font 'Helvetica,16'
set border lw 2.5
set tics scale 1.2
set ylabel "{/Helvetica=22 Population}"
set xlabel "{/Helvetica=22 time (ps)}"
set mytics 2
set mxtics 2
set grid xtics ytics

e1=0;e2=8.0;
temp=500.0  #in Kevlin
kbT=0.0083*temp  

p1(x)=exp(-e1/kbT)/(exp(-e1/kbT)+exp(-e2/kbT));
p2(x)=exp(-e2/kbT)/(exp(-e1/kbT)+exp(-e2/kbT));

outfile = 'fig-population.pdf'
set output outfile
set key opaque samplen 1.0 spacing 1.3 font "Helvetica, 14"
plot 'pop_adiabat1.out' u (column(1)/41340):2 w l lc 3 t'FSSH: P_1',\
     'pop_adiabat1.out' u (column(1)/41340):3 w l lc 4 t'FSSH: P_2',\
     p1(x) w l lc 6 lw 5 t'Exact: P_1',\
     p2(x) w l lc 7 lw 5 t'Exact: P_2'


