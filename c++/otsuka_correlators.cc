/* A program to measure correlation functions in idealised NN octahedral ice.
 *
 * Uses single-spin-flip dynamics at a number of decreasing temperature points.
 *
 * Cmdline arguments:
 *     1) system size
 *     2) number of burn-in samples per T-point
 *     3) number of samples per T-point
 *     4) random seed used for Monte-Carlo sampling
 *     5) input file listing a temperature point (in units of J) on each line
 *     6) output filename root
 *        filenames are of the form root_T[#T-point].corr 
 *
 * Copyright (C) 2022 Attila Szabó <attila.szabo@physics.ox.ac.uk>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License,
 * version 2, as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

#include <vec3.hh>
#include <spinice.hh>
#include <ewald_generator.hh>
#include <ewald_rs.hh>
#include <spin_correlator.hh>
#include <basic_stat.hh>
#include <octa.hh>

#include <cmath>
#include <cstdio>
#include <cstdlib>

using namespace std;

void print_line(const basic_stat& bs) {
    printf("%11.2f %7d %7d %7d", bs.E, bs.M[0], bs.M[1], bs.M[2]);
    for (int i = 0; i < 7; ++i)
        printf(" %6d", bs.n_q[i][0]);
    for (int i = 0; i < 7; ++i)
        printf(" %6d", bs.n_q[i][1]);
    printf("\n");
    fflush(stdout);
}

int main (int argc, char** argv) {
    unsigned
        L = atoi(argv[1]),
        burnin = atoi(argv[2]),
        sample = atoi(argv[3]),
        seed = atoi(argv[4]);

    // set up MC simulation
    ewald_generator gen(L,L,L, ewald_generator::nn_type::nn);
    gen.nn(ewald_generator::nn_type::nnn_octa, -1.0, true);

    spin_correlator sc(L,L,L, true);
    sc.load_pot(gen);

    spinice s (L,L,L, NULL, &sc);
    s.add_J1(1.0);
    s.add_J2(-1.0);
    s.seed(seed);
    s.randomize();

    // run the simulation
    FILE* Tpoints = fopen(argv[5], "r");
    int nT;
    fscanf(Tpoints, "%d", &nT);
    for (int iT = 0; iT < nT; ++iT) {
        double T;
        fscanf(Tpoints, "%lf", &T);
        printf("# Starting T = %.5f\n#", T);
        fflush(stdout);

        for (int i = 0; i < burnin; ++i) {
            s.otsuka(T, false);
            printf("."); fflush(stdout);
        }
        printf("\n");
        s.reset_corr();

        long n_success = 0;
        double n_total = 3.0*L*L*L * sample;
        basic_stat bs;
        for (int i = 0; i < sample; ++i) {
            n_success += s.otsuka(T, true, &bs);
            print_line(bs);
        }
        char outfn_buf[256];
        sprintf(outfn_buf, "%s/T%d.corr", argv[6], iT);
        s.save_corr(outfn_buf);
        printf("# %10.6f %12.10f # T, Acceptance\n\n", T, n_success/n_total);
    }
    fclose(Tpoints);

    return 0;
}
