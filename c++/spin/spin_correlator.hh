/* Class spin_correlator: calculates spin correlators in reciprocal space
 *
 * Ising spins are Fourier transformed sublattice by sublattice and correlators
 * are computed for each pair of sublattices using FFTW
 *         
 * Created on 01/12/2017
 * Copyright (C) 2017-8, 2020-1 Attila Szabó <attila.szabo@physics.ox.ac.uk>
 */

#ifndef spin_correlator_hh
#define spin_correlator_hh

#include <complex>
#include <fftw3.h>

class ewald_generator;

class spin_correlator {
    // Number of reciprocal lattice points per dimension and in total
    const size_t Na, Nb, Nc;
    const size_t N;

    // Number of stored correlators
    size_t n_corr;

    /* Arrays of Fourier transforms */
    
    // Arrays to keep the spin FT
    std::complex<double>* spin[3];
    
    // Collects Ising spin correlators
    std::complex<double>* correlator[3][3];

    // Arrays to keep potential term FTs
    std::complex<double>* potential[3][3];

    /* FT workhorse array */    
    std::complex<double>* fft;
    fftw_complex*         fft_w;

    // FFTW plan for normal Fourier transforms
    fftw_plan plan;

    // Indexing of arrays in a row-major order as required by FFTW
    size_t index(int x, int y, int z) const;

public:
    // Whether the object can handle potential terms
    const bool has_potential;

    // Constructor
    spin_correlator(size_t Na, size_t Nb, size_t Nc, bool has_potential = false);
    // Destructor
    ~spin_correlator();

    // Set real space input, indexing by dimension
    inline std::complex<double>& operator()(int x, int y, int z) {
        return fft[index(x,y,z)];
    }
    
    inline const
    std::complex<double>& operator()(int x, int y, int z) const {
        return fft[index(x,y,z)];
    }

    // Returns number of stored correlators
    size_t count() const {return n_corr;}

    // Increase correlator count
    void incr_correlator(size_t n) {n_corr += n;}

    // Zero out input array
    void clear_input();

    // Zero out correlator arrays
    void reset();

    // Set Ising spins of a given sublattice from the real space input
    void set_ising(size_t i);

    // Add current spin set-up to correlator tally
    void add_correlator(size_t n_sample);

    // Save correlators into binary file
    void save_corr(std::FILE* f) const;

    // Load potential terms from an ewald_generator object
    void load_pot(const ewald_generator& gen); 

    // Save FT'd potential terms into binary file
    void save_pot(std::FILE* f) const;

    // Get Ewald summed energy
    double energy() const;
};

#endif
