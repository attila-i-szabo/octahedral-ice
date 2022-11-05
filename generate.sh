#!/bin/bash
# Generates the data used in the paper
# Assumes there is a directory data/ for keeping the output files

# Loop-update simulations of dipolar ice
mkdir data/loop
for i in `seq 0 274`;
do
    bash/loop.sh $i data/loop
done

# Single-spin-flip simulations for low-temperature dipolar ice
mkdir data/singlespin
for i in `seq 0 274`;
do
    bash/singlespin.sh $i data/singlespin
done

# Loop-update simulations of near-neighbour ice
mkdir data/exchange
for J2 in 0.95 1.00 1.05;
do
    mkdir data/exchange/"$J2"
    for n in `seq 8`;
    do
        bash/exchange.sh $J2 $n data/exchange/
    done
done

# Loop-string simulation of idealised ice
mkdir data/exchange/otsuka
c++/otsuka_correlators.out 64 64 8192 123 Tpoints/otsuka data/exchange/otsuka > data/exchange/otsuka/log

# Detailed samples for calculating autocorrelation times
mkdir data/autocorr
for type in loop singlespin;
do
    mkdir data/autocorr/"$type"
    c++/all_samples.out 16 0.0 -0.8 -0.05 1024 16384 1 $type 543 Tpoints/autocorr data/autocorr/"$type"
done

# Fine-grained simulations near transitions
mkdir data/transition
bash/transition.sh 3k 1.0 0.2
bash/transition.sh FM 1.0 -0.3
bash/transition.sh AFM2k_high 0.0 0.2
bash/transition.sh AFM2k_low 0.0 0.2
bash/transition.sh FM2k 0.0 -0.3
