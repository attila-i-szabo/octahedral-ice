/* Class spinice
 *
 * Implementation of spinice class using single spin flip dynamics
 *
 * Created on 13/08/2019
 * Copyright (C) 2019,2021 Attila Szabó <attila.szabo@physics.ox.ac.uk>
 */

#ifndef spinice_hh
#define spinice_hh

#include <vec3.hh>
#include <random>
#include <vector>
#include <cstdio>

class spin;
class octa;
class basic_stat;
class stat_0pi;
class spin_correlator;
class ewald_rs;

class spinice {
    // Number of unit cells along different directions
    const int na, nb, nc;

    // Number of spins
    const size_t N;

    // Number of lattice sites / octahedra
    const size_t Nl;

    // Spins, Ising components stored for speed
    spin** spins;
    double* spin_signs;

    // Octahedra
    octa** octas;
    
private:
    // Stores Ewald summed potential terms in real space
    ewald_rs* ewald;

    // Correlators of spins in reciprocal space
    bool new_magnetic;
    spin_correlator* magnetic;

    // Random number engine
    std::mt19937 random;
    // Uniform random numbers between [0,1]
    std::uniform_real_distribution<double> dist;
    
    //---------- ROUTINES TO ASSEMBLE STRUCTURE ------------------------------
    // NB basic unit of distance is a_0/2

    // Returns index of cell corresponding to reduced position (units of a_0)
    int cell_at (vec3_int pos) const;
    
    /* Returns spin pointer corresponding to given position,
     * provided there is one there */
    spin* spin_at (vec3_int pos) const;

    /* Returns octaahedron pointer at given position,
     * provided there is one there */
    octa* octa_at (vec3_int pos) const;

    // Returns valid cubic lattice points inside the box
    vec3_int lattice (size_t n) const;

    // Constructs the full spin ice structure
    void init();

    /* Adds a set of near-neighbour couplings
     * spin[r,sl1] and spin[r+dr,sl2] get coupled with strength J for all r 
     * coupling registered both ways */
    void add_coupling(vec3_int dr, int sl1, int sl2, double J);

    //---------- PRIVATE IMPLEMENTATION OF MC --------------------------------
    // Evaluates the local "Ising field" for a given spin
    double field (size_t i) const;

    // Record stuff
    void MC_record(bool record, basic_stat* bs);

    //---------- PUBLIC INTERFACE --------------------------------------------
public:
    // External magnetic field times atomic moment (viz. in units of energy)
    vec3 B;

    // preloaded/precomputed Ewald sum terms, pass NULL for pure NN ice
    // spin_correlator object created on the fly if not supplied
    spinice (size_t na, size_t nb, size_t nc, ewald_rs* ew, 
             spin_correlator* sc = NULL);
    ~spinice();

    // add various near-neighbour couplings
    // +J favours spin ice, -J favours AIAO
    void add_J1 (double J);
    // +J is AFM, -J is FM
    void add_J2 (double J);
    void add_J2p (double J);

    // Seed random number generator from preset list of seeds
    void seed (unsigned no);
    // Resets default orientations of spins
    void reset_spins();
    // Randomises spins
    void randomize();

    // single spin flip MC step, returns number of successful flips
    int montecarlo (double T, bool record, basic_stat* bs = NULL);

    // MC step using Otsuka's algorithm; assumes model is ideal NN ice
    int otsuka (double T, bool record, basic_stat* bs = NULL);
    
    // loop-update MC step, returns number of flipped spins and attempts
    std::vector<int> loop_update(double T, bool record, basic_stat* bs = NULL);
    
    // Record stuff, split by even/odd coordinates
    void MC_record (stat_0pi& stat);

    // Number of spin sites
    size_t n_spin() const {return N;}
    /* Read-only pointer to spin
     * Arg:
     *    n: position of spin in array spins */
    const spin* get_spin(size_t n) const {return spins[n];}
    /* Read-only pointer to spin
     * Arg:
     *    pos: position of spin [a_fcc/2]
     * Throws exception if there's no spin at given point */
    const spin* get_spin(const vec3_int& pos) const {return spin_at(pos);}

    // Number of octahedra
    size_t n_octa() const {return N/3;}
    /* Read-only pointer to octahedron
     * Arg:
     *    n: position of octahedron in array tetras */
    const octa* get_octa(size_t n) const {return octas[n];}
    /* Read-only pointer to octahedron
     * Arg:
     *    pos: position of octahedron [a_fcc/2]
     * Throws exception if there's no octahedron at given point */
    const octa* get_octa(const vec3_int& pos) const {return octa_at(pos);}
    
    // Change correlator object
    void set_magnetic(spin_correlator* m);

    /* Saves correlators to given fliename */
    void save_corr(const char* s) const;
    
    // Save current configuration as a bit string
    void save_config(std::FILE* f) const;

    // Reset correlators
    void reset_corr();
    
    // Calculate energy 
    // parameter indicates whether magnetic object already contains the spins
    double energy(bool filled = false);
};

#endif
