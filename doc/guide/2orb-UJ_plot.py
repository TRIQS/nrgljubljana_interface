import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math

from pytriqs.gf import *
from pytriqs.archive import *
from nrgljubljana_interface import MeshReFreqPts

def A_to_nparrays(A):
    lx = np.array(list(A.mesh.values()))
    ly = np.array(A.data[:,0,0].real)
    return lx, ly

def save_A(fn, gf):
  f = open(fn, "w")
  for w in gf.mesh:
    z = gf[w][0,0]
    f.write("%f %f\n" % (w, -1.0/math.pi*z.imag))

with HDFArchive('2orb-UJ_solution.h5','r') as ar:
    # Expectation values
    print("<n>=",ar['expv']['n_d1'])
    print("<n^2>=",ar['expv']['n_d1^2'])

    save_A("A_w_NRG.dat", ar['G_w']['imp'])

    # Spectral function
    A_w = ar['A_w']['imp']
    lx, ly = A_to_nparrays(A_w)
    plt.plot(lx, ly)
    plt.show()
    plt.savefig("plot.pdf")
