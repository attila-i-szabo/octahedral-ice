/* Class octa
 * Describes an octahedron centre and thus a monopole
 *
 * Created on 16/11/2017
 * Copyright (C) 2017, 2018, 2021 Attila Szabó <attila.szabo@physics.ox.ac.uk>
 */

#ifndef octa_hh
#define octa_hh

#include <vec3.hh>
#include <random>
#include <utility>

class spin;

class octa {
public:
    // Possible sublattices of octahedron centres
    enum class sublattice {A = 1, B = -1};
    
private:
    // Neighbouring spins
    spin* neigh[6];

    // Sublattice octahedron belongs to
    const sublattice m_sl;
    
    // Position (in units of a_0/2, must be all even)
    const vec3_int m_pos;

public:
    // Construct from coordinates in primitive unit cell
    // Sublattice is determined by this position
    octa(const vec3_int& pos);

    // Register with a spin
    void reg(spin* s, unsigned index);

    // Get sublattice index
    int sublat() const {return (int)m_sl;}
    sublattice sublat_enum() const {return m_sl;}
    int sublat_01() const;

    // Get position
    vec3_int pos() const {return m_pos;}

    // Get corresponding cubic lattice site
    vec3_int cubic() const {return m_pos/2;}

    // Get a neighbour
    spin* neighbour(size_t n) const {return neigh[n];}

    //---------- PHYSICAL PARAMETERS -----------------------------------------
    
    // Monopole charge
    int charge() const;

    //---------- OTSUKA'S ALGORITHM ------------------------------------------

    // Makes the portion of the graph inside this octahedron
    void pair_up(const double* w, std::mt19937& random);

    // Returns the pair of the given spin
    spin* get_pair(const spin* s) const;
    
    //---------- LOOP UPDATE ALGORITHM ---------------------------------------
    // Next spin & octahedron in the chain
    std::pair<spin*, octa*> next;
    
    // Fills out `next` randomly from outgoing edges
    bool get_next(std::mt19937& random);

private:
    // pairs of each spin in Otsuka's algorithm
    spin* pairs[6];
};

#endif
