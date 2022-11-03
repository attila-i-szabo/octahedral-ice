/* Implementation of class spin_correlator.
 *
 * Created on 01/12/2017
 * Copyright (C) 2017-8, 2020-1 Attila Szabó <attila.szabo@physics.ox.ac.uk>
 */

#include "spin_correlator.hh"
#include "ewald_generator.hh"
#include <misc.hh>
#include <complex>
#include <fftw3.h>
#include <cstdio>
// Memory operations like memset
#include <cstring>

using namespace std;

spin_correlator::spin_correlator (size_t Na, size_t Nb, size_t Nc, bool pot):
    Na(Na),
    Nb(Nb),
    Nc(Nc),
    N(Na*Nb*Nc),
    has_potential(pot),
    n_corr(0)
{
    // Creating various arrays, all of size N
    // Zeroing out correlator arrays
    for(size_t i = 0; i < 3; i++) {
        spin[i] = new complex<double>[N];
        for(size_t j = 0; j < 3; j++) {
            correlator[i][j] = new complex<double>[N];
            memset(correlator[i][j], 0, N*sizeof(complex<double>));
            if (has_potential)
                potential[i][j] = new complex<double>[N];
        }
    }

    // Creating Fourier transform arrays using FFTW methods
    // and casting them onto complex<double>
    fft_w = fftw_alloc_complex(N);
    fft = (complex<double>*) fft_w;

    // Creating FFTW plan
    plan  = fftw_plan_dft_3d(Na, Nb, Nc, fft_w, fft_w,
                             FFTW_FORWARD, FFTW_MEASURE);
}

spin_correlator::~spin_correlator() {
    // Removing various arrays
    for(size_t i = 0; i < 3; i++) {
        delete[] spin[i];
        for(size_t j = 0; j < 3; j++) {
            delete[] correlator[i][j];
            if (has_potential)
                delete[] potential[i][j];
        }
    }

    // Removing Fourier transform arrays
    // NB complex<double> pointers are but typecasts, no need to delete[] them
    fftw_free(fft_w);

    // Removing FFTW plan
    fftw_destroy_plan(plan);
}

// Indexing of arrays in a row-major order as required by FFTW
size_t spin_correlator::index(int x, int y, int z) const {
    // Sanitise input, PBC
    x = mod(x, Na);
    y = mod(y, Nb);
    z = mod(z, Nc);
    return (x*Nb + y)*Nc + z;
}

void spin_correlator::clear_input() {
    memset(fft, 0, N*sizeof(complex<double>));
}

void spin_correlator::reset() {
    n_corr = 0; // functionality of part class reset
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            memset(correlator[i][j], 0, N*sizeof(complex<double>));
}

void spin_correlator::set_ising(size_t i) {
    fftw_execute(plan);
    memcpy(spin[i], fft, N*sizeof(complex<double>));
}

/* We simply need to calculate sigma_i(q) sigma_j*(q) for all q,i,j
 * and add them to whatever already is in correlator */
void spin_correlator::add_correlator(size_t n_sample) {
    for (size_t i=0; i<3; i++)
        for (size_t j=0; j<3; j++)
            for (size_t q=0; q<N; q++)
                correlator[i][j][q] += spin[i][q] * conj(spin[j][q]);
    n_corr += n_sample;
}

// Save correlators into binary file
void spin_correlator::save_corr(FILE* f) const {
    // array in which to divide total correlators by n_corr
    complex<double>* div = new complex<double>[N];
    double rec = 1.0 / n_corr;
    
    // for each array, divide each q-point by n_corr, then output
    for(size_t i=0; i<3; ++i)
        for(size_t j=0; j<3; ++j) {
            for(size_t q=0; q < N; ++q)
                div[q] = correlator[i][j][q] * rec;
            fwrite(div, sizeof(complex<double>), N, f);
        }

    // delete tmp array
    delete[] div;
}

void spin_correlator::load_pot(const ewald_generator& gen) {
    if (!has_potential)
        throw "Cannot store potential terms";
    for (int i = 0; i < 3; ++i) 
        for (int j = 0; j < 3; ++j) {
            const double* source = gen.ptr(i,j);
            for (int q = 0; q < N; ++q)
                fft[q] = source[q];
            fftw_execute(plan);
            memcpy(potential[i][j], fft, N*sizeof(complex<double>));
        }
}

void spin_correlator::save_pot(FILE* f) const {
    for(size_t i=0; i<3; ++i)
        for(size_t j=0; j<3; ++j) 
            fwrite(potential[i][j], sizeof(complex<double>), N, f);
}

double spin_correlator::energy() const {
    if (!has_potential)
        throw "Cannot compute energy without potential terms";
    double result = 0.0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int q = 0; q < N; ++q)
                result += real(potential[i][j][q] * spin[i][q] * conj(spin[j][q]));
    // factor of N for FFT, factor of 2 for double counting
    return result/(2*N);
}
