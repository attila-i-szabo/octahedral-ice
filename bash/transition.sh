#!/bin/bash
# Cmdline arguments:
#   1) name of transition studied
#   2) J2
#   3) J2'

mkdir data/transition/"$1"
c++/dipolar_correlators.out 16 0.0 -$2 $3 128 512 48 loop 345 Tpoints/transition/"$1" data/transition/"$1" > data/transition/"$1"/log
