#!/bin/bash
N=200
M=200
h=0.001
a=1

betaValues="0 0.1 0.5 1"
kMaxValues="10 50 100"
for beta in $betaValues; do
    for kMax in $kMaxValues; do
        printf "$N\n$M\n$h\n$a\n$beta\n$kMax\nkMax=$kMax+beta=$beta.dat" | ./la
    done
done

gnuplot cmds.gnuplot
