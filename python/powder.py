'''Computes powder neutron-scattering patterns from raw correlators.

Cmdline arguments:
1) input correlator array
2) linear system size
3) largest desired wave vector (2π/a_0)
4) if "order", smearing of discrete k-points is halved.
    This makes Bragg peaks sharper, but produces artefacts for diffuse patterns.

Output is a text table of k(2π/a_0), powder intensity

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
from scipy.special import erf
from sys import argv
#import matplotlib.pyplot as plt

def correlator_tensor(S, n):
    '''
    S: array of correlators with shifted origins
        shape: (3,3,L,L,L)
    n: number of unit cells to be computed along each direction
    output shape: (3,3,2*n*L,2*n*L,2*n*L)
    '''
    S = np.tile(S,(2*n,2*n,2*n))
    # phase factors
    Lf = 2*n*L
    φ = np.exp(-0.5j * np.linspace(-2*np.pi*n, 2*np.pi*n, num=Lf, endpoint=False))
    ι = np.ones(Lf)
    # the Lf-length axis broadcasts on the corresponding wave-vector axis
    # if wanted to do 3D rather than (hhk), it would go as Lf,1,1; 1,Lf,1; 1,1,Lf
    φx = np.stack((φ,ι,ι)).reshape(3,1,Lf,1,1)
    φy = np.stack((ι,φ,ι)).reshape(3,1,1,Lf,1)
    φz = np.stack((ι,ι,φ)).reshape(3,1,1,1,Lf)
    φxc = φx.conj().reshape(1,3,Lf,1,1)
    φyc = φy.conj().reshape(1,3,1,Lf,1)
    φzc = φz.conj().reshape(1,3,1,1,Lf)

    S = S * φx*φy*φz * φxc*φyc*φzc
    print("Expanding S done")
    return S

def NS(S):
    '''takes extended correlators in (hkl) space and returns the full NS intensity'''
    Lf = S.shape[-1]
    # absolute magnitude of q doesn't matter, only its direction
    q = np.linspace(-1,1, num=Lf, endpoint=False)
    h,k,l = np.meshgrid(q,q,q,indexing='ij')
    q = np.stack((h,k,l))
    dipolar = np.einsum('ihkl,jhkl->ijhkl',q,q)/np.einsum('ihkl,ihkl->hkl',q,q)
    dipolar = np.eye(3).reshape(3,3,1,1,1) - dipolar
    dipolar[:,:,Lf//2,Lf//2,Lf//2] = np.eye(3) * 2/3

    return np.einsum('ijhkl,ijhkl->hkl',S,dipolar).real

def powder(S,n,order=False):
    '''calculates powder average from the output of NS'''
    N = S.shape[-1]
    L = N//(2*n)
    # Gaussian blur used to smooth out the discrete k-points
    # we measure q in units of 2pi/a
    if order:
        # for ordered patterns, smaller blur is OK
        σ = 0.5/(2*np.pi)**0.5/L
    else:
        σ = 1.2/(2*np.pi)**0.5/L
    # magnitude of each k-point
    q = np.linspace(-n, n, num=N, endpoint=False)
    h,k,l = np.meshgrid(q,q,q,indexing='ij',sparse=True)
    q = (h**2+k**2+l**2)**0.5
    # wave vector bins matched to L,n
    if order:
        kbin = (np.arange(5*L,10*N+1)-0.5) / (20*L)
    else:
        kbin = (np.arange(2*L,4*N+1)-0.5) / (8*L)
    kavg = (kbin[1:] + kbin[:-1])/2
    dk = kbin[1] - kbin[0]
    # portion of Gaussian falling inside each wave vector bin
    powder_ = np.zeros(kavg.shape)
    erfs = erf((kbin[0]-q) / (2**0.5*σ))
    powder_int = np.einsum('hkl,hkl', S, erfs)
    print(kavg.shape)
    for i,k in enumerate(kbin[1:]):
        print(i,end=',')
        erfs = erf((k-q) / (2**0.5*σ))
        new_powder_int = np.einsum('hkl,hkl',S,erfs)
        powder_[i] = new_powder_int - powder_int
        powder_int = new_powder_int
    powder_ = powder_ / (4 * np.pi * kavg**2 * dk * L**3)
    return kavg,powder_

if __name__ == "__main__":
    L = int(argv[2])
    n = int(argv[3])
    order = len(argv)>4 and argv[4] == 'order'
    S = np.fromfile(argv[1], dtype=complex).reshape(3,3,L,L,L)
    S = correlator_tensor(S,n)
    S = NS(S)
    print("full NS intensity done")
    k,pwdr = powder(S,n,order=order)
    #plt.plot(k,powder)
    #plt.show()
    np.savetxt(argv[1]+'.powder', np.column_stack((k,pwdr)))
