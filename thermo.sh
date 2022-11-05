#!/bin/bash
# Generates thermodynamic observables from data produced by generate.sh


# Figure 3a
for J2 in `seq 0.00 0.05 1.20`;
do
    for J2p in `seq -0.30 0.05 0.20`;
    do
        python3 python/heatcap_combi.py "$J2"_"$J2p"
    done
done

cd data/singlespin/
python3 ../../python/entropy_grid.py 
cd -

# Figure 2
for J2 in 0.95 1.00 1.05;
do
    cat data/exchange/"$J2"/*/energy > data/exchange/"$J2"/energy
    python3 python/heatcap.py data/exchange/"$J2" Tpoints/exchange 16 256
done