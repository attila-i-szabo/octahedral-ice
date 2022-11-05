#!/bin/bash
# Cmdline arguments:
#   1) index of point in parameter space
#   2) output root directory

J2s=(`seq 0.00 0.05 1.20`)
J2ps=(`seq -0.30 0.05 0.20`)

J2=${J2s[$1/11]}
J2p=${J2ps[$1%11]}

dir="$2"/"$J2"_"$J2p"
mkdir $dir
c++/dipolar_correlators.out 16 0.0 -$J2 $J2p 512 4096 48 single 321 Tpoints/singlespin $dir 8 > "$dir"/log
