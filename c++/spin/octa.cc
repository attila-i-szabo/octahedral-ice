/* Implementation of functions in class octa 
 *
 * Created on 16/11/2017
 * Copyright (C) 2017, 2018, 2021 Attila Szabó <attila.szabo@physics.ox.ac.uk>
 */

#include "octa.hh"
#include "spin.hh"
#include <vec3.hh>
#include <random>
#include <algorithm>
#include <vector>
#include <utility>

static octa::sublattice pos2sl(const vec3_int& pos) {
    if ((pos[0]+pos[1]+pos[2]) % 2)
        return octa::sublattice::B;
    else
        return octa::sublattice::A;
}

octa::octa(const vec3_int& pos):
    m_pos(pos),
    m_sl(pos2sl(pos/2))
{}

void octa::reg(spin* s, unsigned index) {
    neigh[index] = s;
}

int octa::sublat_01() const {
    switch (m_sl) {
    case octa::sublattice::A: 
        return 0;
    case octa::sublattice::B:
        return 1;
    }
}

int octa::charge() const {
    /* Ising spin is +1 if it points to +x/y/z
     * Sign convention makes them +1 if they point out of this octa */
    int q = (neigh[0] -> ising() +
             neigh[1] -> ising() +
             neigh[2] -> ising() -
             neigh[3] -> ising() -
             neigh[4] -> ising() -
             neigh[5] -> ising() );
    // Divide by 2 to get +/-1 for single monopole
    return q/2;
}


static int compatible[4][4] = {
    {1, 9, 18, 6},
    {1, 8, 12, 0},
    {1, 5, 0, 0},
    {1, 0, 0, 0}
};

// inline static int abs(int x) {return x < 0 ? -x : x;}

void octa::pair_up(const double* w, std::mt19937& random) {
    // clean up previous pairing
    for (int i = 0; i < 6; ++i)
        pairs[i] = NULL;

    // separate in and out spins, shuffle them
    std::vector<int> ins, outs;
    for (int i = 0; i < 6; ++i)
        if (neigh[i]->ising() == (i < 3 ? 1 : -1))
            ins.push_back(i);
        else
            outs.push_back(i);
    shuffle(ins.begin(), ins.end(), random);
    shuffle(outs.begin(), outs.end(), random);

    // divine the number of pairings needed
    int q = abs(charge());
    std::discrete_distribution<int> dist_pair {
        compatible[q][0]*w[0],
        compatible[q][1]*w[1],
        compatible[q][2]*w[2],
        compatible[q][3]*w[3]
    };
    int n_pair = dist_pair(random);

    for (int i = 0; i < n_pair; ++i) {
        pairs[ins[i]] = neigh[outs[i]];
        pairs[outs[i]] = neigh[ins[i]];
    }
}

spin* octa::get_pair(const spin* s) const {
    for (int i = 0; i < 6; ++i)
        if (neigh[i] == s)
            return pairs[i];
    throw "The given spin doesn't belong to this octahedron";
}

bool octa::get_next (std::mt19937& random) {
    // list in spins, shuffle, take first as random choice
    std::vector<int> ins;
    for (int i = 0; i < 6; ++i)
        if (neigh[i]->ising() == (i < 3 ? 1 : -1))
            ins.push_back(i);
    if (ins.size() == 0)
        return false;
        
    shuffle(ins.begin(), ins.end(), random);
    
    spin* next_spin = neigh[ins[0]];    
    octa* next_octa = next_spin -> other(this);
    next = std::make_pair(next_spin, next_octa);
    return true;
}
