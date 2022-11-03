/* Implementation of functions in class spin 
 *
 * Created on 29/11/2017
 * Copyright (C) 2017, 2018, 2021 Attila Szabó <attila.szabo@physics.ox.ac.uk>
 */

#include "spin.hh"
#include "octa.hh"
#include <vec3.hh>
#include <vector>
#include <utility>

const vec3_int spin::directions[3] = {
    vec3_int(1, 0, 0),
    vec3_int(0, 1, 0),
    vec3_int(0, 0, 1)};

static spin::sublattice pos2sl(const vec3_int& pos) {
    if (pos[0] % 2) 
        return spin::sublattice::x;
    else if (pos[1] % 2) 
        return spin::sublattice::y;
    else
        return spin::sublattice::z;
}

spin::spin(const vec3_int& pos, unsigned index):
    m_pos(pos),
    m_sl(pos2sl(pos)),
    m_sign(1),
    index(index)
{}

void spin::reg(octa* t) {
    switch(t->sublat_enum()) {
    case octa::sublattice::A: m_A = t;
        break;
    case octa::sublattice::B: m_B = t;
        break;
    }
}

octa* spin::other (const octa* one) const {
    if (one == m_A)
        return m_B;
    else if (one == m_B)
        return m_A;
    else
        throw "Spin doesn't belong to given octahedron!";
}

double spin::field() const {
    double result = 0.0;
    for (auto& n: neighbours) 
        result += n.second * n.first->m_sign;
    return result;
}
