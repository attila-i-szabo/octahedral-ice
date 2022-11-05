''' Generates the combined order parameter plot in Fig. 3b.

Should be run in data/loop. 

Copyright (C) 2022 Attila Szabó <attila.szabo@physics.ox.ac.uk>

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License,
version 2, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
'''

import numpy as np
from sys import argv

L = 16

J2s = np.arange(25)/20
J2xs = np.arange(-6,5)/20

OP = np.zeros((J2s.size,J2xs.size))

for i,J2 in enumerate(J2s):
    for j,J2x in enumerate(J2xs):
        S = np.fromfile(f'{J2:.2f}_{J2x:.2f}/T69.corr.symm',dtype=complex)
        S = S.reshape(3,3,L,L,L)[:,:,L//2,L//2,0].real
        OP[i,j] = S[0,0] + S[1,1] - S[2,2]

OP.tofile('combi_order_parameter')
