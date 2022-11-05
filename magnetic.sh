#!/bin/bash
# Generates structure factors and other magnetic observables from data produced by generate.sh

# Figure 3b

for J2 in `seq 0.00 0.05 1.20`;
do
    for J2p in `seq -0.30 0.05 0.20`;
    do
        python3 python/symmetrise.py data/loop/"$J2"_"$J2p"/T69.corr 16
    done
done

cd data/loop/
python3 ../../python/OP_grid.py
cd -

# Neutron patterns

## Near-neighbour

for J2 in 0.95 1.00 1.05;
do
    python3 python/average.py data/exchange/"$J2"/*/T99.corr data/exchange/"$J2"/T99.corr
    python3 python/symmetrise.py data/exchange/"$J2"/T99.corr 16
    python3 python/neutron.py data/exchange/"$J2"/T99.corr.symm 16 3
    python3 python/powder.py data/exchange/"$J2"/T99.corr 16 4
done

python3 python/symmetrise.py data/exchange/otsuka/T90.corr 64
python3 python/neutron.py data/exchange/otsuka/T90.corr.symm 64 3

## Dipolar, single-spin-flip

for dir in "0.40_-0.05" "0.80_-0.05";
do
    python3 python/symmetrise.py data/singlespin/"$dir"/T24.corr 16
    python3 python/neutron.py data/singlespin/"$dir"/T24.corr.symm 16 3
    python3 python/powder.py data/singlespin/"$dir"/T24.corr 16 4
done

## Dipolar, loop-update, also order parameters

for J2 in 0.00 1.00;
do
    for J2p in -0.30 0.20;
    do
        for T in `seq 0 69`;
        do
            python3 python/symmetrise.py data/loop/"$J2"_"$J2p"/T"$T".corr 16
        done
        python3 python/neutron.py data/loop/"$J2"_"$J2p"/T69.corr.symm 16 3
        python3 python/powder.py data/loop/"$J2"_"$J2p"/T69.corr 16 4 order
        python3 OP.py data/loop/"$J2"_"$J2p" Tpoints/loop 16
    done
done

# Order parameters near transition

for dir in `ls data/transition`;
do
    for corr in `ls data/transition/"$dir"/*.corr`;
    do
        python3 python/symmetrise.py $corr 16
    done
    python3 OP.py data/transition/"$dir" Tpoints/transition/"$dir" 16
done