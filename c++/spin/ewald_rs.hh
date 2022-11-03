/* Class for providing real-space Ewald sum terms
 *
 * Created on 06/07/2021
 * Copyright (C) 2021 Attila Szabó <attila.szabo@physics.ox.ac.uk>
 */

#ifndef ewald_rs_hh
#define ewald_rs_hh

#include <vec3.hh>
#include <misc.hh>
#include <cstdlib>
#include <cstdio>

class ewald_generator;

class ewald_rs {
protected:
    // Number of reciprocal lattice points per dimension and in total
    const int Na, Nb, Nc;
    const size_t N;
    
    /* Array of potential terms
     * Stored as one big array of doubles, indexed as [i, x, y, z, j], where
     * i and j index sublattices within the cubic unit cell.
     * potential[i,x,y,z,j] gives the pot. between spins at R_i and r+R_j */
    double *potential;

    // Indexing of arrays in a row-major order
    inline size_t index(int x, int y, int z, int i, int j) const {
        x = mod(x,Na);
        y = mod(y,Nb);
        z = mod(z,Nc);
        return (((size_t(i) * Na + x) * Nb + y) * Nc + z) * 3 + j;
    }
    inline size_t index(vec3_int v, int i, int j) const {
      return index(v[0],v[1],v[2],i,j);
    }

    inline double& operator()(int x, int y, int z, int i, int j) {
        return potential[index(x,y,z,i,j)];
    }

public:
    ewald_rs(size_t na, size_t nb, size_t nc):
        Na(na), Nb(nb), Nc(nc), N(na*nb*nc)
    {
        potential = new double[9*N];
    }

    ewald_rs(size_t na, size_t nb, size_t nc, const ewald_generator& gen):
        Na(na), Nb(nb), Nc(nc), N(na*nb*nc)
    {
        potential = new double[9*N];
        load(gen);
    }

    virtual ~ewald_rs() {delete[] potential;}

    // Load potential terms from an ewald_generator object
    void load(const ewald_generator& gen);

    // Returns potential term
    inline double operator()(vec3_int v, int i, int j) const {
        return potential[index(v,i,j)];
    }

    // Returns pointer to potential term; useful for speeding up bulk access
    inline const double* ptr(int x, int y, int i, int z = 0, int j = 0) const {
        return &(potential[index(x,y,z,i,j)]);
    }
};
#endif
