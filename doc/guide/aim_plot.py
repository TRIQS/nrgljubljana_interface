import numpy as np
import matplotlib as mpl
#mpl.use('PDF')
import matplotlib.pyplot as plt

from pytriqs.gf import *
from pytriqs.archive import *
from nrgljubljana_interface import MeshReFreqPts

def A_to_nparrays(A):
  lx = []
  ly = []
  for w in A.mesh:
    lx.append(float(w))
    ly.append(A[w][0,0].real)
  lx = np.array(lx)
  ly = np.array(ly)
  return lx, ly

with HDFArchive('aim_solution.h5','r') as ar:
    # Expectation values
    print("<n>=",ar['expv']['n_d'])
    print("<n^2>=",ar['expv']['n_d^2'])

    # Spectral function
    A_w = ar['A_w']['imp']
    lx, ly = A_to_nparrays(A_w)

    plt.plot(lx, ly)
    plt.show()
#    plt.savefig("A_w.pdf")
