/* Class basic_stat
 * Collects basic statistics of the spin ice sample, namely magnetisation and
 * number of monopoles of different charges and sublattices
 *
 * Created on 14/08/2019
 * Copyright (C) 2019, 2021 Attila Szabó <attila.szabo@physics.ox.ac.uk>
 */

#ifndef basic_stat_hh
#define basic_stat_hh

#include <vec3.hh>
#include <cstring>
#include <cstdio>
#include <cstdint>

struct basic_stat {
    // monopole count; charge q, sublattice x stored in n_q[q+3][x]
    unsigned n_q[7][2];
    
    vec3_int M; // magnetisation
    double E; // energy
    
    inline void clear() {
        M = vec3_int();
        memset(n_q, 0, 14*sizeof(unsigned));
    }

    inline unsigned& Q(int q, int x) {return n_q[q+3][x];}

    // Calculates sum of q^2 for all tetrahedra
    inline int SJ () {
        return 
            (n_q[0][0] + n_q[0][1] + n_q[6][0] + n_q[6][1]) * 9 +
            (n_q[1][0] + n_q[1][1] + n_q[5][0] + n_q[5][1]) * 4 +
            (n_q[2][0] + n_q[2][1] + n_q[4][0] + n_q[4][1]);
    }
    
    // Calculates sum of all spins (that is, 2 * sum of q on A tetrahedra)
    inline int Sh () {
        return 
            (-int(n_q[0][0]) + int(n_q[6][0])) * 6 +
            (-int(n_q[1][0]) + int(n_q[5][0])) * 4 +
            (-int(n_q[2][0]) + int(n_q[4][0])) * 2;
    }
};

struct stat_0pi {
    // monopole count; charge q, sublattice x stored in n_q[q+3][x]
    uint16_t n_q[7][2];
    
    // magnetisation, separated by even/odd x,y,z coordinates
    // spin at (x,y,z) is stored in M[floor(x)%2][floor(y)%2][floor(z)%2][sl]
    int16_t M[2][2][2][3];
    double E; // energy
    
    inline void clear() {
        memset(n_q, 0, 14*sizeof(uint16_t));
        memset(M, 0, 24*sizeof(int16_t));
    }

    inline uint16_t& Q(int q, int x) {return n_q[q+3][x];}

    // Calculates sum of q^2 for all tetrahedra
    inline int SJ () {
        return 
            (n_q[0][0] + n_q[0][1] + n_q[6][0] + n_q[6][1]) * 9 +
            (n_q[1][0] + n_q[1][1] + n_q[5][0] + n_q[5][1]) * 4 +
            (n_q[2][0] + n_q[2][1] + n_q[4][0] + n_q[4][1]);
    }
    
    // Calculates sum of all spins (that is, 2 * sum of q on A tetrahedra)
    inline int Sh () {
        return 
            (-int(n_q[0][0]) + int(n_q[6][0])) * 6 +
            (-int(n_q[1][0]) + int(n_q[5][0])) * 4 +
            (-int(n_q[2][0]) + int(n_q[4][0])) * 2;
    }
};

#endif
