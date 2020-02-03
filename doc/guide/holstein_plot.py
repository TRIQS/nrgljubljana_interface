import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from pytriqs.gf import *
from pytriqs.archive import *
from nrgljubljana_interface import MeshReFreqPts

def A_to_nparrays(A):
    lx = np.array(list(A.mesh.values()))
    ly = np.array(A.data[:,0,0].real)
    return lx, ly

with HDFArchive('holstein_solution.h5','r') as ar:
    A_w = ar['A_w']['imp']     # Spectral function
    lx, ly = A_to_nparrays(A_w)
    plt.plot(lx, ly)
    plt.show()
