/* Implementation of class ewald_rs
 *
 * Created on 22/07/2021
 * Copyright (C) 2021 Attila Szabó <attila.szabo@physics.ox.ac.uk>
 */

#include "ewald_rs.hh"
#include "ewald_generator.hh"

void ewald_rs::load(const ewald_generator& gen) {
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int x = 0; x < Na; ++x)
                for (int y = 0; y < Nb; ++y)
                    for (int z = 0; z < Nc; ++z)
                        (*this)(x,y,z,i,j) = gen(x,y,z,i,j);
}