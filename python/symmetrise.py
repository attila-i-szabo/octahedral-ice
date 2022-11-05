''' Applies cubic point group symmetry to raw correlators.

Cmdline arguments:
1) input file
2) linear system size

Output is in same format as input, saved into [input filename].symm. 

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

L = int(argv[2])

def reflect(S, axis):
    '''Calculates correlators of a configuration reflected across coordinate 
    axis #axis.'''

    # nontrivial phase factors
    φr = -np.exp(1j * np.linspace(0,2*np.pi,num=L,endpoint=False))
    # should only appear in spin coordinate/sublattice #axis
    φ = np.ones((3,L),dtype=complex)
    φ[axis] = φr
    # and in spatial coordinate #axis
    ax = [1,1,1]
    ax[axis] = L
    φ = φ.reshape(3,1,*ax)
    φc = φ.conj().reshape(1,3,*ax)

    # reverse the appropriate axis
    S = np.roll(np.flip(S,axis=axis+2),1,axis=axis+2)

    return S*φ*φc

def axis_permute(S, *order):
    '''Reshuffles spin coordinates/sublattices and spatial coordinates.'''
    
    return np.transpose(S, (0,1, order[0]+2, order[1]+2, order[2]+2))[list(order)][:,list(order)]

S = np.fromfile(argv[1], dtype=complex).reshape(3,3,L,L,L)
S = S + reflect(S,0)
S = S + reflect(S,1)
S = S + reflect(S,2)
S = S + axis_permute(S, 1,0,2)
S = S + axis_permute(S, 1,2,0) + axis_permute(S, 2,0,1)
S /= 48
S.tofile(argv[1]+'.symm')
    
