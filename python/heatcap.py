''' Computes thermodynamic data from raw energy measurements.

Cmdline arguments: 
1) working directory
2) file with T-points
3) system size
4) #independent runs

Output has seven columns: T, E, std(E)_1, std(E)_2, C_1, C_2, S
    index 1 means variance taken by independent run and averaged
    index 2 means variance taken across all runs

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

T = np.loadtxt(argv[2])
N_Tp = int(T[0])
T = T[1:(N_Tp+1)]

L = int(argv[3])
N = 3*L**3

N_runs = int(argv[4])

energy = np.fromfile(argv[1]+'/energy').reshape(N_runs, N_Tp, -1)

E = np.average(energy, axis=(0,2))
vE1 = np.average(np.var(energy, axis=2), axis=0)
C1 = vE1/T**2
vE2 = np.var(energy, axis=(0,2))
C2 = vE2/T**2

# More numerically stable near transitions and doesn't depend on equilibration
dS = (E[1:]-E[:-1]) / (T[1:] + T[:-1]) * 2
S = [0.0] + np.cumsum(dS).tolist()

# Add high-T expansion approximation to entropy
S = S + E[0]/2/T[0] + N*np.log(2)

# Return order: E, ΔE by run, ΔE overall, C by run, C overall, S
allE = np.column_stack((E,vE1**0.5,vE2**0.5,C1,C2,S)) / N

todos = np.column_stack((T,allE))

np.savetxt(argv[1]+'/heatcap', todos[::-1])
