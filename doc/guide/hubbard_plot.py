import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math

from triqs.gf import *
from h5 import *
from nrgljubljana_interface import MeshReFreqPts

def GF_to_nparrays(A):
    lx = np.array(list(A.mesh.values()))
    ly = np.array(-1.0/math.pi*A.data[:,0,0].imag)
    return lx, ly

with HDFArchive('solution.h5','r') as ar:
    # Expectation values
    print("<n>=",ar['expv']['n_d'])
    print("<n^2>=",ar['expv']['n_d^2'])

    # Spectral function (of local lattice GF)
    A_w = ar['Gloc']['imp']
    lx, ly = GF_to_nparrays(A_w)
    plt.plot(lx, ly)
    plt.show()
