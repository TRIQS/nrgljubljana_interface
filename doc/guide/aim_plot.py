from pytriqs.gf import *
from pytriqs.archive import *
from pytriqs.plot.mpl_interface import oplot
from nrgljubljana_interface import Solver, FlatNew

with HDFArchive('aim_solution.h5','r') as ar:
    print(ar['expv']['n_d^2'])
#    oplot(ar['A_w']['imp'], '-o', x_window = (-10,10))
