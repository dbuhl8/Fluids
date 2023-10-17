

set title "Gaussian Process Generation"
set xlabel "Time"
set ylabel "Magnitude"

plot 'plotData/pgen.dat' u 1:2 w l lt rgb "black" title "Neg Exp", \
     'plotData/opartial.dat' u 1:2 pt 5 lt rgb "black" title "Neg Exp Rolling Window", \
     'plotData/pgen.dat' u 1:3 w l lt rgb "red"  title "Sin", \
     'plotData/opartial.dat' u 1:3 pt 5 lt rgb "red" title "Sin Rolling Window", \
     'plotData/pgen.dat' u 1:4 w l lt rgb "blue" title "Cos", \
     'plotData/opartial.dat' u 1:4 pt 5 lt rgb "blue" title "Cos Rolling Window", \
     'plotData/pgen.dat' u 1:5 w l lt rgb "forest-green" title "Abs Val", \
     'plotData/opartial.dat' u 1:5 pt 5 lt rgb "dark-green" title "Abs Val Rolling Window"
