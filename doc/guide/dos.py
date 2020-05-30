from triqs.gf import *
from nrgljubljana_interface import MeshReFreqPts, hilbert_transform_refreq, hilbert_transform_elementwise
import numpy as np
import math

dosfn = lambda eps : 2.0/math.pi*np.sqrt(1-eps**2) if abs(eps)<1.0 else 0.0
lin_mesh = np.linspace(-1.0, 1.0, num=2001, endpoint=True)
m = MeshReFreqPts(lin_mesh)
dosA = Gf(mesh=m, target_shape=[])
for w in dosA.mesh:
  dosA[w] = dosfn(float(w))

f = open("dos.dat", "w");
for w in dosA.mesh:
    f.write("%f %f\n" % (float(w), dosA[w].real))
