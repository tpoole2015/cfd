# see http://www.gnuplotting.org/plotting-data/ for gnuplot commands
set macros
POS = "at graph 0.05,0.9 font ',8'"

set terminal png
set output 'linear_advection.png'

set xrange [0:1]
set yrange [-0.2:1.2]

unset key

set terminal
set multiplot layout 4,3 rowsfirst

do for [beta in "0 0.1 0.5 1"] {
    do for [kMax in "10 50 100"] {
        myLabel = "kMax=".kMax.",beta=".beta
        myFile = "kMax=".kMax."+beta=".beta.".dat"
        set label 1 myLabel @POS
        plot myFile using 1:2 with lines
    }
}

unset multiplot
