

set title "Gaussian Process Generation"
set xlabel "Time"
set ylabel "Magnitude"

#This file is called from within the directory of the makefile and so it is 
#important to reference the .dat files in the way that this file has them

set terminal png size 1920, 1080
set output "../plots/partial.png"
plot '../plotData/pgen.dat' u 1:2 w l lt rgb "black" title "Gaussian Process", \
     '../plotData/opartial.dat' u 1:2 pt 5 lt rgb "black" title "Rolling Window"