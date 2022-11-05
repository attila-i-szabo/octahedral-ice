''' Processes thermodynamic data from loop-update and single-spin-flip simulations
of dipolar spin ice and generates entropy curves by merging the two datasets.

Command-line argument: J2_J2' (directory name)

`heatcap` files has seven columns: T, E, std(E)_1, std(E)_2, C_1, C_2, S
    index 1 means variance taken by independent run and averaged
    index 2 means variance taken across all runs

`entropy` files have four columns: T, E, S, C_2

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
N = 3*L**3

# parameters of the loop-update run
T1 = 10**(np.arange(39,-31,-1)/30)
N1 = 512 # number of samples per run

# parameters of the single-spin run
T2 = 10**(-np.arange(25)/40)
N2 = 512

def E_C(samples,T):
    E = np.average(samples, axis=(0,2))
    vE1 = np.average(np.var(samples, axis=2), axis=0)
    C1 = vE1/T**2
    vE2 = np.var(samples, axis=(0,2))
    C2 = vE2/T**2
    return E,vE1,vE2,C1,C2

def S(E,T):
    # More numerically stable near transitions and doesn't depend on equilibration
    dS = (E[1:]-E[:-1]) / (T[1:] + T[:-1]) * 2
    S = [0.0] + np.cumsum(dS).tolist()

    # Add high-T expansion approximation to entropy
    S = S + E[0]/2/T[0] + N*np.log(2)
    return S

def table(EC,S,T):
    E, vE1, vE2, C1, C2 = EC
    allE = np.column_stack((E,vE1**0.5,vE2**0.5,C1,C2,S)) / N
    todos = np.column_stack((T,allE))
    return todos[::-1]

# process loop-update data
samples = np.fromfile(f'data/loop/{argv[1]}/energy').reshape(-1, T1.size, N1)
EC1 = E_C(samples,T1)
S1 = S(EC1[0],T1)

# process single-spin-flip data
samples = np.fromfile(f'data/singlespin/{argv[1]}/energy').reshape(-1,T2.size,N2)
EC2 = E_C(samples,T2)
E_cum = np.concatenate((EC1[0][:39], [(EC1[0][39]+EC2[0][0])/2], EC2[0][1:]))
C_cum = np.concatenate((EC1[3][:39], [(EC1[3][39]+EC2[3][0])/2], EC2[3][1:]))
T_cum = np.concatenate((T1[:39], T2))
S_cum = S(E_cum, T_cum)
S2 = S_cum[39:]

# output
np.savetxt(f'data/loop/{argv[1]}/heatcap', table(EC1,S1,T1))
np.savetxt(f'data/singlespin/{argv[1]}/heatcap', table(EC2,S2,T2))
np.savetxt(f'data/singlespin/{argv[1]}/entropy',
           np.column_stack((T_cum, E_cum/N, S_cum/N, C_cum/N))[::-1])
