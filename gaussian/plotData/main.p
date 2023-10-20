

set title "Gaussian Process Regression"
set xlabel "Time"
set ylabel "Amplitude"

set terminal png size 1920, 1080
set output "../plots/main.png"

plot '../plotData/regressiondata.dat' u 1:2 lt rgb "red" w l title "Gaussian Regression",\
     'plotData/regressiondata.dat' u 1:3 lt rgb "black" w l title "Actual - Sin", \
     'plotData/traindata.dat' u 1:2 title "Training Data"
