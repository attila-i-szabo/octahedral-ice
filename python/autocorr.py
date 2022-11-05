''' Computes the autocorrelation function from raw samples.

Cmdline arguments: 
1) input filename root (without `.samples`)
2) linear system size 
3) largest time for which autocorrelation is required

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
from scipy.signal import correlate
from sys import argv

L = int(argv[2])
N = 3*L**3

# load file, unpack the bits
σ = np.fromfile(argv[1]+'.samples', dtype=np.uint8)
σ = np.unpackbits(σ, bitorder='little').reshape(-1, N)
σ = 2.0 * σ - 1.0

# we're going to limit memory usage with the "valid" mode in correlate
# this function operates on every dimension, so it will also automatically
# collate over every spin

tmax = int(argv[3])
corr = correlate(σ, σ[:-tmax], mode='valid')
corr /= corr[0,0]

np.savetxt(argv[1]+'.corr', corr)
