''' Generates neutron scattering pattern in the [hhl] plane.

Cmdline arguments: 
1) raw correlator tensor
2) linear system size (L)
3) highest desired h/l (n)

Output saved as binary arrays of size (2nL+1)*(2nL+1). 

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
import matplotlib.pyplot as plt

def correlator_tensor(S, n):
    '''
    S: (hhk) cut from array of correlators with shifted origins
        shape: (3,3,L,L)
    n: number of unit cells to be computed along h and k directions
    output shape: (3,3,2*n*L+1,2*n*L+1)
    '''
    Lf = 2*n*L+1
    S = np.tile(S,(2*n+1,2*n+1))[:,:,:Lf,:Lf]
    # phase factors
    φ = np.exp(-0.5j * np.linspace(-2*np.pi*n, 2*np.pi*n, num=Lf))
    ι = np.ones(Lf)
    # the Lf-length axis broadcasts on the corresponding wave-vector axis
    # if wanted to do 3D rather than (hhk), it would go as Lf,1,1; 1,Lf,1; 1,1,Lf
    φx = np.stack((φ,ι,ι)).reshape(3,1,Lf,1)
    φy = np.stack((ι,φ,ι)).reshape(3,1,Lf,1)
    φz = np.stack((ι,ι,φ)).reshape(3,1,1,Lf)
    φxc = φx.conj().reshape(1,3,Lf,1)
    φyc = φy.conj().reshape(1,3,Lf,1)
    φzc = φz.conj().reshape(1,3,1,Lf)

    S = S * φx*φy*φz * φxc*φyc*φzc
    return S

def NS(S):
    '''takes an extended correlator cut along the (hhk) axis and returns the
    unpolarised NS intensity'''
    Lf = S.shape[-1]
    # absolute magnitude of q doesn't matter, only its direction
    q = np.linspace(-1,1, num=Lf)
    h,k = np.meshgrid(q,q,indexing='ij')
    q = np.stack((h,h,k))
    dipolar = np.einsum('ihk,jhk->ijhk',q,q)/np.einsum('ihk,ihk->hk',q,q)
    dipolar = np.eye(3).reshape(3,3,1,1) - dipolar
    dipolar[:,:,Lf//2,Lf//2] = np.eye(3) * 2/3

    return np.einsum('ijhk,ijhk->hk',S,dipolar).real

def chiral(S):
    '''takes an extended correlator cut along the (hhk) axis and returns the
    chiral NS intensity'''
    Lf = S.shape[-1]
    # absolute magnitude of q doesn't matter, only its direction
    q = np.linspace(-1,1, num=Lf)
    h,k = np.meshgrid(q,q,indexing='ij')
    q = np.stack((h,h,k))
    q /= np.sum(q**2, axis=0,keepdims=True)**0.5
    q[:,Lf//2,Lf//2] = 0
    # construct "cross product" tensors
    chi = np.cross(np.eye(3).reshape(3,3,1,1),q.reshape(3,1,Lf,Lf),axis=0)

    return np.einsum('ijhk,ijhk->hk',S,chi).real
    
def NSF(S):
    '''takes an extended correlator cut along the (hhk) axis and returns the
    non-spin-flip NS intensity with [1,-1,0] neutrons'''
    Lf = S.shape[-1]
    parallel = np.array([1,-1,0, -1,1,0, 0,0,0]).reshape(3,3,1,1)/2
    result =  np.einsum('ijhk,ijhk->hk',S,parallel).real
    #result[Lf//2,Lf//2] = 0
    return result


if __name__ == "__main__":
    L = int(argv[2])
    n = int(argv[3])
    S = np.fromfile(argv[1], dtype=complex).reshape(3,3,L*L,L)[:,:,::(L+1),:]
    S = correlator_tensor(S,n)
    total = NS(S)
    print(np.max(total), np.sum(total))
    plt.imshow(total)
    plt.colorbar()
    plt.show()
    total.tofile(argv[1]+'.total')
    nsf = NSF(S)
    print(np.max(nsf), np.sum(nsf))
    plt.imshow(nsf)
    plt.colorbar()
    plt.show()
    nsf.tofile(argv[1]+'.nsf')
    sf = total-nsf
    print(np.max(sf), np.sum(sf))
    plt.imshow(sf)
    plt.colorbar()
    plt.show()
    sf.tofile(argv[1]+'.sf')
    '''
    chi=chiral(S)
    plt.imshow(chi)
    plt.colorbar()
    plt.show()
    chi.tofile(argv[1]+'.chi')
    '''
