''' Averages raw binary files.

Cmdline arguments: all input files, output filename

Input files are expected to be binary files of `double`s of equal size. 

Copyright (C) 2022 Attila Szabó <attila.szabo@physics.ox.ac.uk>

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License,
version 2, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
'''

from sys import argv
import numpy as np

cumul = 0
n = 0
for fn in argv[1:-1]:
  cumul = cumul + np.fromfile(fn)
  n += 1
cumul /= n
cumul.tofile(argv[-1])
