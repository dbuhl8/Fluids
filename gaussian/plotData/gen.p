
set title "Gaussian Process Generation"

set xlabel "Time"
set ylabel "Magnitude"

set terminal png size 1600, 900
set output "../plots/gen.png"
# This file is called within the makefile in a different director
# and will need to navigate to the plotData directory
plot "../plotData/ogen.dat" u 1:2 w l title "GP"




