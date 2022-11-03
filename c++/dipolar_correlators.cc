/* A program to measure correlation functions in dipolar octahedral ice.
 *
 * Uses single-spin-flip and (optionally) loop-update dynamics at a number of
 * decreasing temperature points.
 *
 * Cmdline arguments:
 *     1)   system size
 *     2-4) J1, J2, J2'
 *     5)   number of burn-in samples per T-point
 *     6)   number of samples per T-point
 *     7)   number of annealing cycles
 *     8)   algorithm to be used: (l)oop-update or (s)ingle-spin-flip
 *     9)   random seed used for Monte-Carlo sampling
 *     10)  input file listing a temperature point (in units of D) on each line
 *     11)  output filename root
 *          filenames are of the form root/T[#T-point].corr
 *     12)  frequency of recording samples (optional, default: 1)
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
#include <ctime>
#include <vector>

using namespace std;

void print_line(const basic_stat& bs) {
    printf("%11.6f %5d %5d %5d", bs.E, bs.M[0], bs.M[1], bs.M[2]);
    for (int i = 0; i < 7; ++i)
        printf(" %5d", bs.n_q[i][0]);
    for (int i = 0; i < 7; ++i)
        printf(" %5d", bs.n_q[i][1]);
    printf("\n");
    fflush(stdout);
}

void save_stat (const stat_0pi& bs, FILE* Efile, FILE* Mfile) {
    fwrite(&(bs.E), sizeof(double), 1, Efile);
    fwrite(bs.M, sizeof(int16_t), 24, Mfile);
}

int main (int argc, char** argv) {
    unsigned
        L = atoi(argv[1]),
        burnin = atoi(argv[5]),
        sample = atoi(argv[6]),
        cycle = atoi(argv[7]),
        seed = atoi(argv[9]),
        freq = 1;

    if (argc > 12)
        freq = atoi(argv[12]);

    double
        J1 = atof(argv[2]),
        J2 = atof(argv[3]),
        J2p = atof(argv[4]);

    bool consider_loop_update = argv[8][0] == 'l'; 

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
    spin_correlator** scs = new spin_correlator*[nT];
    for (int iT = 0; iT < nT; ++iT) {
        fscanf(Tfile, "%lf", &Tpoints[iT]);
        scs[iT] = new spin_correlator(L,L,L, true);
        scs[iT] -> load_pot(gen);
    }
    fclose(Tfile);
    
    spinice s (L,L,L, &ew, scs[0]);
    s.seed(seed);

    char outfn_buf[256];
    sprintf(outfn_buf, "%s/samples", argv[11]);
    FILE* f_samples = fopen(outfn_buf, "wb");
    sprintf(outfn_buf, "%s/energy", argv[11]);
    FILE* f_energy = fopen(outfn_buf, "wb");
    sprintf(outfn_buf, "%s/magnet", argv[11]);
    FILE* f_magnet = fopen(outfn_buf, "wb");

    for (int iC = 0; iC < cycle; ++iC) {
        s.randomize();
        for (int iT = 0; iT < nT; ++iT) {
            double T = Tpoints[iT];
            printf("# Starting T = %.5f\n", T);
            auto now = std::time(nullptr);
            printf("# %s#", std::ctime(&now));
            fflush(stdout);
            
            s.set_magnetic(scs[iT]);
        
            bool do_loop_update = (consider_loop_update ? (T < 1.25) : false);

            for (int i = 0; i < burnin; ++i) {
                if (do_loop_update)
                    s.loop_update(T, false);
                s.montecarlo(T, false);
                printf("."); fflush(stdout);
            }
            printf("\n");

            long n_success = 0;
            double n_total = 3.0*L*L*L * sample;
            stat_0pi bs;
            long lu_stat[4] = {0,0,0,0};
            for (int i = 0; i < sample; ++i) {
                if (do_loop_update) {
                    auto stat = s.loop_update(T, false);
                    for (int j = 0; j < 4; ++j)
                        lu_stat[j] += stat[j];
                }
                bool record = (i % freq) == 0;
                n_success += s.montecarlo(T, record);
                if (record) {
                    s.MC_record(bs);
                    save_stat(bs, f_energy, f_magnet);
                }
            }
            if (do_loop_update)
                printf("%10.6f %12.10f %12.10f %5.2f # T, Acceptance, Loop update acceptance, Loop length\n", 
                       T, n_success/n_total, lu_stat[0]/double(lu_stat[1]),
                       lu_stat[1]/double(lu_stat[3]));
            else
                printf("%10.6f %12.10f # T, Acceptance\n", T, n_success/n_total);
                
            s.save_config(f_samples);
            fflush(f_samples);
            fflush(f_energy);
            fflush(f_magnet);
        }
    }
    
    fclose(f_samples);
    fclose(f_energy);
    fclose(f_magnet);
    
    for (int iT = 0; iT < nT; ++iT) {
        sprintf(outfn_buf, "%s/T%d.corr", argv[10], iT);
        FILE* f = fopen(outfn_buf, "wb");
        scs[iT] -> save_corr(f);
        fclose(f);
        delete scs[iT];
    }
    
    delete[] Tpoints;
    delete[] scs;

    return 0;
}
