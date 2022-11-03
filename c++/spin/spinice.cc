/* Implementation of class spinice
 *
 * Created on 13/08/2019
 * Copyright (C) 2019,2021 Attila Szabó <attila.szabo@physics.ox.ac.uk>
 */

#include "spinice.hh"
#include "spin.hh"
#include "octa.hh"
#include "basic_stat.hh"
#include "ewald_rs.hh"
#include "spin_correlator.hh"
#include <vec3.hh>
#include <misc.hh>
#include <cmath>
#include <string>
#include <cstdio>
#include <vector>
// Random numbers
#include <random>
#include <ctime>

using namespace std;


//-------------- ASSEMBLING STRUCTURE -----------------------------------------

/* Generate lattice points inside rectangular box
 * Order: pass unit cells by increasing x, y, and z in precedence
 */
vec3_int spinice::lattice(size_t n) const {
    // Separate various indices
    int z      = n % nc; n /= nc;
    int y      = n % nb; n /= nb;
    int x      = n % na; // hopefully this is a trivial step

    vec3_int v(x,y,z);
    return 2*v;
}

int spinice::cell_at(vec3_int vec) const {
    int x = mod(vec[0], na);
    int y = mod(vec[1], nb);
    int z = mod(vec[2], nc);
    return (x*nb+y)*nc+z;
}

/* Returns the spin sitting at some point or throws an error if there is
 * nothing there.
 * Order in array spins: lattice(n)+e[i] is at spins[3*n+i] */
spin* spinice::spin_at(vec3_int vec) const {
    int sl = -1;
    for (int i=0; i<3; i++) 
        if (alldiv(vec - spin::directions[i], 2)) {
            sl = i;
            break;
        }
    if (sl<0) throw "Invalid position of spin";
    vec -= spin::directions[sl];
    vec /= 2;

    return spins[3*cell_at(vec)+sl];
}

/* Returns the spin sitting at some point or throws an error if there is
 * nothing there.
 * Order in array spins: lattice(n) is at octas[n] */
octa* spinice::octa_at(vec3_int vec) const {
    if (!(alldiv(vec, 2))) throw "Invalid position of octahedron";

    return octas[cell_at(vec/2)];
}

void spinice::add_coupling(vec3_int dr, int sl1, int sl2, double J) {
    for (int ic = 0; ic < Nl; ++ic) {
        vec3_int r = lattice(ic)/2;
        spin* s1 = spins[3*ic + sl1];
        spin* s2 = spins[3*cell_at(r+dr) + sl2];
        s1->reg(s2,J);
        s2->reg(s1,J);
    }
}

void spinice::add_J1(double J) {
    add_coupling(vec3_int(0, 0,0),0,1,+J);
    add_coupling(vec3_int(1, 0,0),0,1,-J);
    add_coupling(vec3_int(0,-1,0),0,1,-J);
    add_coupling(vec3_int(1,-1,0),0,1,+J);

    add_coupling(vec3_int(0,0, 0),1,2,+J);
    add_coupling(vec3_int(0,1, 0),1,2,-J);
    add_coupling(vec3_int(0,0,-1),1,2,-J);
    add_coupling(vec3_int(0,1,-1),1,2,+J);

    add_coupling(vec3_int( 0,0,0),2,0,+J);
    add_coupling(vec3_int( 0,0,1),2,0,-J);
    add_coupling(vec3_int(-1,0,0),2,0,-J);
    add_coupling(vec3_int(-1,0,1),2,0,+J);
}

void spinice::add_J2(double J) {
    add_coupling(vec3_int(1,0,0),0,0,J);
    add_coupling(vec3_int(0,1,0),1,1,J);
    add_coupling(vec3_int(0,0,1),2,2,J);
}

void spinice::add_J2p(double J) {
    add_coupling(vec3_int(0,1,0),0,0,J);
    add_coupling(vec3_int(0,0,1),0,0,J);
    add_coupling(vec3_int(1,0,0),1,1,J);
    add_coupling(vec3_int(0,0,1),1,1,J);
    add_coupling(vec3_int(1,0,0),2,2,J);
    add_coupling(vec3_int(0,1,0),2,2,J);
}

//-------------- CONSTRUCTOR, DESTRUCTOR --------------------------------------

void spinice::init() {
    spins = new spin*[N];
    spin_signs = new double[N];
    octas = new octa*[Nl];

    // Creating spins and tetrahedra
    for(size_t ic = 0; ic < Nl; ic++) {
        for(int ssl = 0; ssl < 3; ssl++) {
            spins[3*ic+ssl] = new spin(lattice(ic) + spin::directions[ssl],
                                       3*ic+ssl);
            spin_signs[3*ic+ssl] = spins[3*ic+ssl] -> ising();
        }
        octas[ic] = new octa(lattice(ic));
    }

    // Registering tetrahedra with neighbouring spins
    for(size_t i = 0; i < Nl; i++) {
        octa* t = octas[i];
        for(int j = 0; j < 3; j++) {
            vec3_int r = t->pos() + spin::directions[j];
            spin* s = spin_at(r);
            s->reg(t);
            t->reg(s,j);

            r = t->pos() - spin::directions[j];
            s = spin_at(r);
            s->reg(t);
            t->reg(s,j+3);
        }
    }

    // Creates magnetic correlator objects as necessary
    new_magnetic = (magnetic == NULL);
    if (new_magnetic)
        magnetic = new spin_correlator(na, nb, nc);
}


spinice::spinice (size_t na, size_t nb, size_t nc, ewald_rs* ew, spin_correlator* sc):
    na(na),
    nb(nb),
    nc(nc),
    N(3lu*na*nb*nc),
    Nl(1lu*na*nb*nc),
    ewald(ew),
    magnetic(sc)
{
    init();
}

spinice::~spinice() {
    for (size_t i = 0; i < N; i++) 
        delete spins[i];
    delete[] spins;
    delete[] spin_signs;

    for (size_t i = 0; i < Nl; i++)
        delete octas[i];
    delete[] octas;

    if (new_magnetic)
        delete magnetic;
}

void spinice::seed(unsigned no) {
    random.seed(no);
}

// Reset spins into default direction: everything is 1
void spinice::reset_spins() {
    for (size_t i = 0; i < N; ++i) {
        spins[i]->set(1);
        spin_signs[i] = 1;
    }
}

// Randomises spins
void spinice::randomize() {
    for (size_t i = 0; i < N; ++i) {
        spin_signs[i] = (dist(random) < 0.5 ? +1 : -1);
        spins[i] -> set(spin_signs[i]);
    }
}

//-------------- MONTE CARLO --------------------------------------------------

double spinice::field (size_t i) const {
    spin* s = spins[i];
    double retval = B % s->dir_plus() + s->field();

    if (!ewald) // no dipolar interactions
        return retval;

    // Dipolar interactions, if tracked
    int sl_i = i % 3; // index of the cubic sublattice spin i belongs to
    vec3_int cubic_i = s->cubic();
    int xi = cubic_i[0],
        yi = cubic_i[1],
        zi = cubic_i[2];
    size_t j = 0;
    for (int x = 0; x < na; ++x)
        for (int y = 0; y < nb; ++y) {
            // exploit that the potential terms for a row of constant x, y are
            // contiguous in memory
            // however, the sum is split up once as ewald terms are cyclically
            // shifted from the spins
	  const double* pot = ewald->ptr(x-xi, y-yi, sl_i, /*z =*/ nc-zi);
            for (int jj = 0; jj < 3*zi; ++jj)
                retval += spin_signs[j+jj] * pot[jj];
            j += 3*zi;

            pot = ewald->ptr(x-xi, y-yi, sl_i);
            for (int jj = 0; jj < 3*(nc-zi); ++jj)
                retval += spin_signs[j+jj] * pot[jj];
            j += 3*(nc-zi);
        }

    return retval;
}

void spinice::MC_record(bool record, basic_stat* bs) {
    // Store spin correlators if needed
    if (record) {
        for (int i = 0; i < 3; ++i) {
            magnetic->clear_input();
            for (size_t n = i; n < N; n += 3) {
                vec3_int v = spins[n]->cubic();
                (*magnetic)(v[0], v[1], v[2]) = spins[n]->ising();
            }
            magnetic->set_ising(i);
        }
        magnetic->add_correlator(N);
    }

    // Record basic stats if object is provided
    if (bs) {
        bs->clear();
        for (size_t i = 0; i < N; ++i)
            bs->M += spins[i]->dir();
        for (size_t i = 0; i < Nl; ++i)
            ++(bs->Q(octas[i]->charge(), octas[i]->sublat_01()));
    }

    // Record energy if equipped to do so
    if (record && bs) {
        bs->E = energy(true);
    }
}

void spinice::MC_record (stat_0pi& bs) {
    bs.clear();
    
    size_t i = 0;
    for (int x = 0; x < na; ++x) 
        for (int y = 0; y < nb; ++y)
            for (int z = 0; z < nc; ++z)
                for (int sl = 0; sl < 3; ++sl) {
                    bs.M[x%2][y%2][z%2][sl] += spin_signs[i];
                    ++i;
                }
                
    for (size_t i = 0; i < Nl; ++i)
        ++(bs.Q(octas[i]->charge(), octas[i]->sublat_01()));
        
    // Record energy 
    bs.E = energy(true);
}

int spinice::montecarlo (double T, bool record, basic_stat* bs) {
    int success = 0;
    for (size_t i = 0; i < N; ++i) {
        double f = field(i);

        if ( f * (spins[i]->ising()) > 0 ) {
            spins[i] -> flip();
            spin_signs[i] *= -1;
            success++;
        } else if ( dist(random) < exp(2*f*(spins[i]->ising())/T) ) {
            spins[i] -> flip();
            spin_signs[i] *= -1;
            success++;
        }
    }

    MC_record(record, bs);

    return success;
}

vector<int> spinice::loop_update (double T, bool record, basic_stat* bs) {
    vector<int> stat = {0,0,0,0};
    std::uniform_int_distribution<> octa_dist(0, Nl-1);
    
    while (stat[1] < N) {
        got_stuck:
        // wipe previous loops
        for (int i = 0; i < Nl; ++i)
            octas[i] -> next = make_pair((spin*)NULL, (octa*)NULL);
        
        // build a chain until a closed loop forms in it
        octa* current = octas[octa_dist(random)];
        while (current -> next.first == NULL) {
            if (!current -> get_next(random))
                goto got_stuck; // `continue` the outer loop
            current = current -> next.second;
        }
        
        // traverse the loop, calculate energy change, flip spins 
        octa* loop_start = current;
        int loop_length = 0;
        double deltaE = 0.0;
        do {
            spin* s = current -> next.first;
            deltaE -= 2 * s->ising() * field(s->index);
            s -> flip();
            spin_signs[s->index] *= -1;
            current = current -> next.second;
            ++loop_length;
        } while (current != loop_start);
        stat[1] += loop_length;
        stat[3]++;
        
        if ( (deltaE > 0) && (dist(random) > exp(-deltaE/T)) ) {
            // unsuccessful, flip everything back
            do {
                spin* s = current -> next.first;
                s -> flip();
                spin_signs[s->index] *= -1;
                current = current -> next.second;
            } while (current != loop_start);
        } else {
            // record that it was successful
            stat[0] += loop_length;
            stat[2]++;
        }
    }
    
    MC_record(record, bs);
    return stat;
}
        

/*
vector<int> spinice::loop_update (double T, bool record, basic_stat* bs) {
    vector<int> stat = {0,0,0,0};
    double E = energy();
    std::uniform_int_distribution<> octa_dist(0, Nl-1);
    
    while (stat[1] < N) {
        // wipe previous loops
        for (int i = 0; i < Nl; ++i)
            octas[i] -> next = make_pair((spin*)NULL, (octa*)NULL);
        
        // build a chain until a closed loop forms in it
        octa* current = octas[octa_dist(random)];
        while (current -> next.first == NULL) {
            current -> get_next(random);
            current = current -> next.second;
        }
        
        // traverse the loop and flip spins *only in spin objects*
        octa* loop_start = current;
        int loop_length = 0;
        do {
            current -> next.first -> flip();
            current = current -> next.second;
            ++loop_length;
        } while (current != loop_start);
        stat[1] += loop_length;
        stat[3]++;
        
        // calculate the energy of the flipped configuration
        double E_new = energy();
        if ( (E_new > E) && (dist(random) > exp((E - E_new)/T)) ) {
            // unsuccessful, flip everything back
            do {
                current -> next.first -> flip();
                current = current -> next.second;
            } while (current != loop_start);
        } else {
            // successful, flip spin_signs and update E
            for (int i = 0; i < N; ++i)
                spin_signs[i] = spins[i] -> ising();
            E = E_new;
            stat[0] += loop_length;
            stat[2]++;
        }
    }
    
    MC_record(record, bs);
    return stat;
}
*/

int spinice::otsuka (double T, bool record, basic_stat* bs) {
    // generate graph
    double w[4];
    w[0] = exp(-18/T);
    w[1] = (exp(-8/T) - w[0])/5;
    w[2] = (exp(-2/T) - w[0] - 8*w[1])/12;
    w[3] = (1 - w[0] - 9*w[1] - 18*w[2])/6;
    for (size_t i = 0; i < Nl; ++i)
        octas[i]->pair_up(w, random);

    // traverse graph and decide what to flip
    int success = 0;
    std::bernoulli_distribution dist_flip(0.5);
    for (size_t i = 0; i < N; ++i)
        spins[i] -> visited = false;
    for (size_t i = 0; i < N; ++i) {
        if (spins[i] -> visited) continue;

        bool flip = dist_flip(random);

        spins[i]->visited = true;
        if (flip) {
            spins[i]->flip();
            ++success;
        }

        // traverse string towards A octahedron
        spin* s = spins[i];
        octa* o; 
        bool is_A = true;
        while (s) {
            if (is_A)
                o = s->A();
            else
                o = s->B();
            is_A = !is_A;
            s = o->get_pair(s);

            // check whether we got to the end of the string
            if (!s) break;
            // check whether we got around a closed loop
            if (s->visited) break;

            s->visited = true;
            if (flip) {
                s->flip();
                ++success;
            }
        }

        // traverse string towards B octahedron
        s = spins[i];
        is_A = false;
        while (s) {
            if (is_A)
                o = s->A();
            else
                o = s->B();
            is_A = !is_A;
            s = o->get_pair(s);

            // check whether we got to the end of the string
            if (!s) break;
            // check whether we got around a closed loop
            if (s->visited) break;

            s->visited = true;
            if (flip) {
                s->flip();
                ++success;
            }
        }
    }

    // ensure spin_signs matches with spins
    for (size_t i = 0; i < N; ++i)
        spin_signs[i] = spins[i]->ising();

    MC_record(record, bs);

    return success;
} 

double spinice::energy (bool filled) {
    double E = 0;
    // There are long-range terms
    if (ewald) {
        // fill out magnetic if necessary
        if (!filled) {
            for (int i = 0; i < 3; ++i) {
                magnetic->clear_input();
                for (size_t n = i; n < N; n += 3) {
                    vec3_int v = spins[n]->cubic();
                    (*magnetic)(v[0], v[1], v[2]) = spins[n]->ising();
                }
                magnetic->set_ising(i);
            }
        }
        E = magnetic->energy();
    }
    // Short-range terms; divide by 2 to avoid double counting
    for (size_t i = 0; i < N; ++i)
        E += spins[i]->ising() * spins[i]->field() / 2;
    return E;
}

//---------- SAVE TO FILE -----------------------------------------------------

void spinice::set_magnetic (spin_correlator* m) {
    if (new_magnetic)
        delete magnetic;
    new_magnetic = false;
    magnetic = m;
}

void spinice::save_corr (const char* s) const {
    FILE* f = fopen(s, "wb");
    magnetic -> save_corr(f);
    fclose(f);
}

void spinice::save_config (FILE* f) const {
    unsigned char* bits = new unsigned char[(N+7)/8] ();
    for (int i = 0; i < N; ++i)
        bits[i/8] |= (spin_signs[i] > 0 ? 1u : 0u) << (i%8);
    fwrite(bits, 1, (N+7)/8, f);
}

void spinice::reset_corr() {
    magnetic -> reset();
}
