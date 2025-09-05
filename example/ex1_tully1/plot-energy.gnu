#!/usr/bin/gnuplot
set encoding iso_8859_1
set key right top

set terminal pdf size 4in,4.5in enhanced color font 'Helvetica,18'
set border lw 2.5
set tics scale 2

set output 'fig-energy.pdf'
NOXTICS = "set xtics (''-10,''-5,''0,''5,''10); \
          unset xlabel"
XTICS="set xtics -10,5,10;\
       set xlabel '{/Helvetica=22 R (a.u.)}'"

set key  samplen 1.0 spacing 1.3 font "Helvetica, 15"
set multiplot layout 2,1

@NOXTICS
set ylabel "{/Helvetica=22 H_{dia} (a.u.)}"
set yr[-.015:0.02]
set ytics 0.01
plot 'energy_surface.out' u 1:2 w l lw 4 t'V_{11}',\
     'energy_surface.out' u 1:5 w l lw 4 t'V_{22}',\
     'energy_surface.out' u 1:3 w l lw 4 dt 2 t'V_{12}',\

@XTICS
set ylabel "{/Helvetica=22 H_{adia} (a.u.)}"
set yr[-.015:0.02]
set ytics 0.01

plot 'energy_surface.out' u 1:6 w l lw 4 t'{/Symbol e}_{1}',\
     'energy_surface.out' u 1:7 w l lw 4 t'{/Symbol e}_{2}',\
     'energy_surface.out' u 1:($10/100) w l lw 4 dt 2 t'd_{12}/100',\
