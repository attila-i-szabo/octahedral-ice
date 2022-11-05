''' Converts temperature series into maps of entropy in parameter space 
at a given temperature.

Must be run in data/loop or data/singlespin.

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

J2s = np.arange(25)/20
J2xs = np.arange(-6,5)/20

S = [] 

for i,J2 in enumerate(J2s):
    S.append([])
    for j,J2x in enumerate(J2xs):
        a = np.loadtxt(f'{J2:.2f}_{J2x:.2f}/heatcap')
        S[-1].append(a[:,6])

Ts = a[:,0]
S = np.array(S)

for i,T in enumerate(Ts):
    np.savetxt(f'entropy/{T:.3f}.txt', S[:,:,i], fmt='%.6f')
    S[:,:,i].tofile(f'entropy/{T:.3f}.bin')
