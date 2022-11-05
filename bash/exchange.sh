#!/bin/bash
# Cmdline arguments:
#   1) J2
#   2) index of independent run
#   3) output root directory

J2=$1
n=$2
dir="$3"/"$J2"/"$n"
mkdir $dir
c++/exchange_correlators.out 16 1.0 -$J2 0.0 256 4096 32 loop 23"$n"1 Tpoints/exchange $dir > "$dir"/log & 
