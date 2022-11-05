''' Computes heat capacity and statistics at wave vectors that are multiples of pi

Cmdline args:
1)   file of T-points
2)   number of independent runs
3)   linear system size

Should be run in the directory containing the raw data.

Outputs a table of
1)   T
2)   energy
3)   heat capacity (all runs mixed together)
4)   heat capacity (average of runs)
5-10)  details of FM dipolar OP
11-16) details of 3k dipolar OP
17-22) details of FM2k quadrupolar OP
23-28) details of AFM2k quadrupolar OP

"details" contain, in this order,
    <|M|>, M_rms, χ (all runs mixed), χ (average of runs), 
    <M_x^4+...>/<M^2>^2, <M^4>/<M^2>^2

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

Ts = np.loadtxt(argv[1])
nT = int(Ts[0])
Ts = Ts[1:(nT+1)]
nR = int(argv[2])
L = int(argv[3])
N = 3*L**3

energy_dat = np.fromfile('energy').reshape(nR,nT,-1)
energy = np.average(energy_dat, axis=(0,2))
C1 = np.var(energy_dat, axis=(0,2)) / Ts**2
C2 = np.average(np.var(energy_dat, axis=2), axis=0) / Ts**2

magnet_dat = np.fromfile('magnet', dtype=np.int16).reshape(nR,nT,-1,8,3)
ks = [[1,1],[1,-1]]
ks = np.kron(np.kron(ks,ks),ks)
x,y,z = np.eye(3)
OPs = [
    [np.outer(ks[0],x), np.outer(ks[0],y), np.outer(ks[0],z)], # FM dipolar
    [np.outer(ks[3],x), np.outer(ks[5],y), np.outer(ks[6],z)], # 3-k dipolar
    [
        np.outer(ks[3], y-z),
        np.outer(ks[5], z-x),
        np.outer(ks[6], x-y),
    ], #FM quadrupolar
    [np.outer(ks[7],x), np.outer(ks[7],y), np.outer(ks[7],z)], # AFM quadrupolar
]

# extract order parameters
M = np.einsum('rtskl,ockl->rtsoc',magnet_dat, OPs)
# ensure that AFM quadrupolar OP components sum to 0
M[:,:,:,-1] -= np.average(M[:,:,:,-1], axis=-1, keepdims=True)

M4 = np.sum(M**4, axis=-1)
M = np.sum(M**2, axis=-1)**0.5

# magnetisation measures
Mav = np.average(M, axis=(0,2))
Mrms = np.average(M**2, axis=(0,2))**0.5
# susceptibility
chi1 = np.var(M, axis=(0,2)) / Ts[:,None]
chi2 = np.average(np.var(M, axis=2), axis=0) / Ts[:,None]
# Binder cumulants
b1 = np.average(M4, axis=(0,2))/Mrms**4
b2 = np.average(M**4, axis=(0,2))/Mrms**4

todos = np.stack((Mav/N,Mrms/N,chi1/N,chi2/N,b1,b2)).transpose(1,2,0).reshape(len(Ts),-1)
todos = np.column_stack((Ts, energy/N, C1/N, C2/N, todos))

np.savetxt('thermo',todos[::-1])
