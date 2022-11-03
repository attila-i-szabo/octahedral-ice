/* Class spin
 * Describes an Ising spin and will allow for extensions to realistic spins
 * like dipoles or solenoids.
 *
 * Created on 26/10/2017
 * Copyright (C) 2017-8, 2020-1 Attila Szabó <attila.szabo@physics.ox.ac.uk>
 */

#ifndef spin_hh
#define spin_hh

#include <vec3.hh>
#include <vector>
#include <utility>

class octa;

class spin {
public:
    // Possible sublattices of spins
    enum class sublattice {x = 0, y = 1, z = 2};

    // Directions of spin vectors
    static const vec3_int directions[3];

    // Each entry contains a near neighbour and the Ising interaction strength
    // between the two
    std::vector<std::pair<spin*,double>> neighbours;
    
protected:
    // Sign of the Ising spin relative to standard orientation
    int m_sign;

    // Position (in units of a_0/2)
    const vec3_int m_pos;

    // Tetrahedron centres connected to this spin
    octa* m_A;
    octa* m_B;

    // Sublattice index
    const sublattice m_sl;
    
public:    
    //---------- CONSTRUCTION, FITTING INTO STRUCTURE -------------------------
    // index in `spins` array of `spinice`
    const unsigned index;
    
    // Can only be constructed with position in units of a_0/2
    // sublattice is computed on the fly
    spin(const vec3_int& pos, unsigned index);
    
    // Virtual destructor to make sure descendants are destroyed
    virtual ~spin() {}

    // Register spin with tetrahedron centres
    void reg(octa* t);
    // Registered tetrahedron centres
    octa* A() const {return m_A;}
    octa* B() const {return m_B;}
    octa* other (const octa* one) const;  

    // Register a neighbour
    void reg(spin* s, double J) {
        neighbours.push_back(std::make_pair(s,J));
    }

    // Get sublattice index
    int sublat() const {return (int)m_sl;}
    
    // Get position
    vec3_int pos() const {return m_pos;}

    // Get cubic lattice reference site
    vec3_int cubic() const {
        return (m_pos - directions[sublat()]) / 2;
    }

    //---------- BASIC SPIN OPERATIONS ---------------------------------------

    // Kill the spin (i.e. set it to 0 forever)
    void kill() {m_sign = 0;}

    // Flip the spin (fine even if the spin is killed)
    void flip() {m_sign *= -1;}

    // Set the Ising spin (may only do if the spin is alive)
    void set(int s) {
        if (m_sign) m_sign = s;
    }

    // Get the Ising spin
    int ising() const {return m_sign;}

    // Principal direction
    vec3_int dir_plus() const {return directions[sublat()];}

    // Direction of spin
    vec3_int dir() const {return m_sign * directions[sublat()];}

    // Spin as a unit vector (note directions contains unit vectors)
    vec3 unit() const {return 1.0 * dir();}

    // Net interaction due to near neighbours
    double field() const;

    // Need this marker for Otsuka algorithm
    bool visited;
};

#endif
