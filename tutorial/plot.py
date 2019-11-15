from pytriqs.gf import *
from pytriqs.archive import *

import numpy
from matplotlib import pyplot as plt

with HDFArchive('Solver.h5','r') as arch:

    x = arch['S/A_w/imp/mesh/points']
    y = arch['S/A_w/imp/data']

    plt.plot(x, y[::,0,0])
    plt.show()
