/* A program to generate raw samples of dipolar octahedral ice.
 *
 * Uses single-spin-flip dynamics at a number of decreasing temperature points.
 *
 * Cmdline arguments:
 *     1)   system size
 *     2-4) J1, J2, J2'
 *     5)   number of burn-in samples
 *     6)   number of samples per T-point
 *     7)   number of annealign cycles
 *     8)   whether it should use (l)oop-update or (s)ingle-spin-flip dynamics
 *     9)   random seed used for Monte-Carlo sampling
 *     10)  input file listing a temperature point (in units of D) on each line
 *     11)  output filename root
 *          filenames are of the form root/T[#T-point]_[#cycle].samples 
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
#include <octa.hh>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <vector>

using namespace std;

int main (int argc, char** argv) {
    unsigned
        L = atoi(argv[1]),
        burnin = atoi(argv[5]),
        sample = atoi(argv[6]),
        cycle = atoi(argv[7]),
        seed = atoi(argv[9]);

    double
        J1 = atof(argv[2]),
        J2 = atof(argv[3]),
        J2p = atof(argv[4]);

    bool do_loop_update = argv[8][0] == 'l';

    printf("# J1:  %.3f\n# J2:  %.3f\n# J2': %.3f\n", J1, J2, J2p);

    // set up MC simulation, no demag factor
    ewald_generator gen(L,L,L, 0.0);
    gen.nn(ewald_generator::nn_type::nn, J1, true);
    gen.nn(ewald_generator::nn_type::nnn_octa, J2, true);
    gen.nn(ewald_generator::nn_type::nnn_perp, J2p, true);

    ewald_rs ew(L,L,L, gen);

    // get the T-points
    FILE* Tfile = fopen(argv[10], "r");
    int nT;
    fscanf(Tfile, "%d", &nT);
    double* Tpoints = new double[nT];
    for (int iT = 0; iT < nT; ++iT) 
        fscanf(Tfile, "%lf", &Tpoints[iT]);
    fclose(Tfile);
    
    spinice s (L,L,L, &ew);
    s.seed(seed);
    
    for (int iC = 0; iC < cycle; ++iC) {
        s.randomize();
        for (int iT = 0; iT < nT; ++iT) {
            double T = Tpoints[iT];
            printf("# Starting T = %.5f\n", T);
            auto now = std::time(nullptr);
            printf("# %s", std::ctime(&now));
            fflush(stdout);
            
            char outfn_buf[256];
            if (cycle > 1)
                sprintf(outfn_buf, "%s/T%d_%d.samples", argv[11], iT, iC);
            else
                sprintf(outfn_buf, "%s/T%d.samples", argv[11], iT);
            FILE* f = fopen(outfn_buf, "wb");
            
            for (int i = 0; i < sample; ++i) {
                if (do_loop_update)
                    s.loop_update(T, false);
                s.montecarlo(T, false);
                s.save_config(f);
                // fflush(f);
            }
            
            fclose(f);
        }
    }
    
    delete[] Tpoints;

    return 0;
}
