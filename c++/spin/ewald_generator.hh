/* Class to evaluate dipolar Ewald sums and return results in real space
 *
 * Created on 13/08/2019
 * Copyright (C) 2019, 2021 Attila Szabó <attila.szabo@physics.ox.ac.uk>
 */

#ifndef ewald_generator_hh
#define ewald_generator_hh

#include <vec3.hh>
#include <misc.hh>
#include <cstdio>

class ewald_generator{
    // Number of reciprocal lattice points per dimension and in total
    const int Na, Nb, Nc;
    const size_t N;

    /* Array of potential terms
     * Stored as 9 arrays of doubles, each storing coupling between a given pair
     * of sublattices
     * potential[i][j][x,y,z] gives the pot. between spins at R_i and (x,y,z)+R_j */
    double* potential[3][3];

    // Indexing of arrays in a row-major order
    inline size_t index(int x, int y, int z) const {
        x = mod(x,Na);
        y = mod(y,Nb);
        z = mod(z,Nc);
        return (size_t(x) * Nb + y) * Nc + z;
    }

    inline double& operator()(int x, int y, int z, int i, int j) {
        return potential[i][j][index(x,y,z)];
    }

    // Indexing FFTW arrays in row-major order
    inline size_t index_fft(int x, int y, int z) const {
        x = mod(x,2*Na);
        y = mod(y,2*Nb);
        z = mod(z,2*Nc);
        return (size_t(x) * 2*Nb + y) * 2*Nc + z;
    }
    inline size_t index_fft(vec3_int v) const {
        return index_fft(v[0],v[1],v[2]);
    }

public:
    enum class nn_type {nn = 0, nnn_octa = 1, nnn_perp = 2};

    ewald_generator(size_t Na, size_t Nb, size_t Nc, double demag = 1.0/3.0);
    ewald_generator(size_t Na, size_t Nb, size_t Nc, nn_type what);
    ewald_generator(size_t Na, size_t Nb, size_t Nc, std::FILE* f);
    ~ewald_generator();

    // Computes dipolar potential terms for a given demag. factor
    void dipolar(double demag = 1.0/3.0);

    // Generates nearest-neighbour terms
    void nn(nn_type what, double J = 1.0, bool cumulative = false);

    /* Load potential terms from file
     * Format: binary, arrays in potential listed one after the other */
    void load(std::FILE* f);

    // Save potential terms to file
    void save(std::FILE* f) const;

    // Returns single potential term
    inline double operator()(int x, int y, int z, int i, int j) const {
        return potential[i][j][index(x,y,z)];
    }

    // Read-only bulk access to arrays
    inline const double* ptr(int i, int j) const {
        return potential[i][j];
    }
};

#endif
