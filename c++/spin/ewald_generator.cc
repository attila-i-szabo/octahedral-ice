/* Implementation of class ewald_generator
 *
 * Created on 13/08/2019
 * Copyright (C) 2019, 2021 Attila Szabó <attila.szabo@physics.ox.ac.uk>
 */

#include "ewald_generator.hh"
#include "spin.hh"
#include <vec3.hh>
#include <misc.hh>
#include <ewald.hh>
#include <complex>
#include <fftw3.h>
#include <cmath>
#include <cstdio>
#include <cstring>

using namespace std;
using namespace ewald_def;

#define SQRT2 1.414213562373095048801688724209698078569671875376948073176
#define PI 3.141592653589793238462643383279502884197169399375105820974

// same as spin::directions
static const vec3_int directions[3] = {
    vec3_int(1, 0, 0),
    vec3_int(0, 1, 0),
    vec3_int(0, 0, 1)};

ewald_generator::ewald_generator(size_t na, size_t nb, size_t nc, double demag):
    Na(na), Nb(nb), Nc(nc), N(na*nb*nc)
{
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            potential[i][j] = new double[N];
    dipolar(demag);
}

ewald_generator::ewald_generator(size_t na, size_t nb, size_t nc, 
                                 ewald_generator::nn_type what):
    Na(na), Nb(nb), Nc(nc), N(na*nb*nc)
{
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            potential[i][j] = new double[N];
    nn(what);
}

ewald_generator::ewald_generator(size_t na, size_t nb, size_t nc, FILE* f):
    Na(na), Nb(nb), Nc(nc), N(na*nb*nc)
{
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            potential[i][j] = new double[N];
    load(f);
}

ewald_generator::~ewald_generator() {
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            delete[] potential[i][j];
}

void ewald_generator::dipolar(double demag) {
    /* Fourier transform array and plans to evaluate reciprocal space sums
     * Real space lattice unit is a_0/2 */
    fftw_complex* fft_w = fftw_alloc_complex(3*8*N);
    vec3_cplx* fft = (vec3_cplx*) fft_w;
    int n[] = {2*Na,2*Nb,2*Nc};
    fftw_plan plan = fftw_plan_many_dft(3,        // rank
                    n,        // dimensions
                    3,        // howmany,
                    fft_w, n, // in, inembed
                    3, 1,     // istride, idist
                    fft_w, n, // out, onembed
                    3, 1,     // ostride, odist
                    FFTW_FORWARD,  // sign = -1
                    FFTW_ESTIMATE); // flags

    /* Constant prefactors of reciprocal and real space sums and surface term
     * See Ewald summation literature (and "Ewald summation" in Notion)
     * Results returned in units of D for general usability */
    const double RECIP = PI * SQRT2 / N;
    const double REAL = 1.0 / (2.0 * SQRT2);
    const double SURFACE = PI * SQRT2 / N * demag;

    // Evaluate and store `fields' due to different directions of spins
    for (int i = 0; i < 3; ++i) {
        vec3_int ei = directions[i];

        // Fill out array in reciprocal space
        for (int x = -Na; x < Na; ++x)
            for (int y = -Nb; y < Nb; ++y)
                for (int z = -Nc; z < Nc; ++z) {
                    vec3 k( 2.0 * PI * x / Na,
                            2.0 * PI * y / Nb,
                            2.0 * PI * z / Nc );
                    fft[index_fft(x,y,z)] =
                        (x||y||z ?                         
                         RECIP * (ei%k)/k.len2() * exp(-0.5*A2*k.len2()) * k
                         : SURFACE * ei );
                }
        // Transform it to obtain `fields' in real space
        fftw_execute(plan);

        // Calculate reciprocal & real space term by sublattice
        for (int j = 0; j < 3; ++j) {
            vec3_int ej = directions[j];

            // Reciprocal space
            for (int x = 0; x < Na; ++x)
              for (int y = 0; y < Nb; ++y)
                for (int z = 0; z < Nc; ++z) {
                    // Position in FFT array
                    vec3_int rf = (2*vec3_int(x,y,z) + ej - ei);
                    (*this)(x,y,z,i,j) = fft[index_fft(rf)].real() % ej;
                }

            // Real space - (2L)^3 cubic unit cells
            for (int x = -L; x < L; ++x)
              for (int y = -L; y < L; ++y)
                for (int z = -L; z < L; ++z) {
                    // Physical separation vector
                    vec3 r = vec3(x,y,z) + (ej-ei)/2.0;
                    double l = r.len();

                    if (l < 1e-3) {
                        // Self-interaction is set to zero
                        (*this)(x,y,z,i,j) = 0.0;
                    } else {
                        // Otherwise add real space term
                        (*this)(x,y,z,i,j) +=
                            REAL * ((ei%ej) * B(l) - (ei%r)*(ej%r) * C(l));
                    }
                }
        }
    }

    // Ensure there is no self-interaction
    for (int i = 0; i < 3; ++i)
        (*this)(0,0,0,i,i) = 0.0;

    fftw_free(fft_w);
    fftw_destroy_plan(plan);
}

void ewald_generator::nn(ewald_generator::nn_type what, double J, bool cumul) {
    if (!cumul)
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                memset(potential[i][j], 0, N*sizeof(double));

    switch(what) {
    case ewald_generator::nn_type::nn:
        (*this)(0, 0,0,0,1) += J;
        (*this)(1, 0,0,0,1) -= J;
        (*this)(0,-1,0,0,1) -= J;
        (*this)(1,-1,0,0,1) += J;

        (*this)( 0,0,0,1,0) += J;
        (*this)(-1,0,0,1,0) -= J;
        (*this)( 0,1,0,1,0) -= J;
        (*this)(-1,1,0,1,0) += J;

        (*this)(0,0, 0,1,2) += J;
        (*this)(0,1, 0,1,2) -= J;
        (*this)(0,0,-1,1,2) -= J;
        (*this)(0,1,-1,1,2) += J;

        (*this)(0, 0,0,2,1) += J;
        (*this)(0,-1,0,2,1) -= J;
        (*this)(0, 0,1,2,1) -= J;
        (*this)(0,-1,1,2,1) += J;

        (*this)( 0,0,0,2,0) += J;
        (*this)( 0,0,1,2,0) -= J;
        (*this)(-1,0,0,2,0) -= J;
        (*this)(-1,0,1,2,0) += J;

        (*this)(0,0, 0,0,2) += J;
        (*this)(0,0,-1,0,2) -= J;
        (*this)(1,0, 0,0,2) -= J;
        (*this)(1,0,-1,0,2) += J;
        break;
    case ewald_generator::nn_type::nnn_octa:
        (*this)( 1,0,0,0,0) += J;
        (*this)(-1,0,0,0,0) += J;
        (*this)(0, 1,0,1,1) += J;
        (*this)(0,-1,0,1,1) += J;
        (*this)(0,0, 1,2,2) += J;
        (*this)(0,0,-1,2,2) += J;
        break;
    case ewald_generator::nn_type::nnn_perp:
        (*this)(0, 1,0,0,0) += J;
        (*this)(0,-1,0,0,0) += J;
        (*this)(0,0, 1,0,0) += J;
        (*this)(0,0,-1,0,0) += J;

        (*this)( 1,0,0,1,1) += J;
        (*this)(-1,0,0,1,1) += J;
        (*this)(0,0, 1,1,1) += J;
        (*this)(0,0,-1,1,1) += J;

        (*this)( 1,0,0,2,2) += J;
        (*this)(-1,0,0,2,2) += J;
        (*this)(0, 1,0,2,2) += J;
        (*this)(0,-1,0,2,2) += J;
        break;
    default:
        throw "Invalid NN type";
    }
}

void ewald_generator::load(std::FILE* f) {
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            fread(potential[i][j], sizeof(double), N, f);
}

void ewald_generator::save(std::FILE* f) const {
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            fwrite(potential[i][j], sizeof(double), N, f);
}
    
