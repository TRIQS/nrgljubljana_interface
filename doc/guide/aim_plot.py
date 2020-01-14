import matplotlib as mpl
mpl.use('PDF')

from pytriqs.gf import *
from pytriqs.archive import *
from pytriqs.plot.mpl_interface import *
from nrgljubljana_interface import MeshReFreqPts

with HDFArchive('aim_solution.h5','r') as ar:
    # Expectation values
    print("<n>=",ar['expv']['n_d'])
    print("<n^2>=",ar['expv']['n_d^2'])

    a = ar['A_w']['imp']
    g = GfReFreq(indices=[0], window=(-2,2), n_points=1000, name='imp')
    for w in g.mesh:
      g[w] = a(w.value) # not implemented yet (Jan 2020)

    print("g_dens=", g.density()) # seems incorrect ??

    oplot(g, '-o', mode = 'R', name = "A")
    plt.show()
    plt.savefig("A_w.pdf")
