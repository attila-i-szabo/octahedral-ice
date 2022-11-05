'''Generates order parameter plots as a function of temperature

Cmdline arguments:
1) directory with (symmetrised) correlator arrays
2) list of T-points
3) L

Output has five columns: T, Gamma, (π,π,π), (π,π,0) normal to spin, (π,π,0) in spin plane

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

L = int(argv[3])
π = L//2

# read T-points
T = np.loadtxt(argv[2])
N_Tp = int(T[0])
T = T[1:(N_Tp+1)]

OP = np.zeros((N_Tp,4))

for i in range(N_Tp):
    S = np.fromfile(f'{argv[1]}/T{i}.corr.symm', dtype=complex).reshape(3,3,L,L,L)[2,2].real
    OP[i,0] = S[0,0,0]
    OP[i,1] = S[π,π,π]
    OP[i,2] = S[π,π,0]
    OP[i,3] = S[0,π,π]

todos = np.column_stack((T, 3*OP/L**3))
np.savetxt(f'{argv[1]}/OP',todos[::-1])
