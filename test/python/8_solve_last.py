#!/usr/bin/env python2

import unittest

from nrgljubljana_interface import Solver, SemiCircular, MeshReFreqPts

from h5 import *
from triqs.utility.h5diff import h5diff
from triqs.utility import mpi
from triqs.utility.comparison_tests import *

class test_SIAM(unittest.TestCase):

    # Construct Parameters
    cp = {}
    cp["model"] = "SIAM"
    cp["symtype"] = "QS"
    cp["mesh_max"] = 1.0
    cp["mesh_min"] = 1e-3
    cp["mesh_ratio"] = 1.1

    # Set up the Solver
    S = Solver(**cp)

    # Solve Parameters
    sp = {}
    sp["T"] = 1e-3
    sp["Lambda"] = 4.0
    sp["Nz"] = 2
    sp["Tmin"] = 1e-4
    sp["keep"] = 50
    sp["keepenergy"] = 6.0

    # Model Parameters
    mp = {}
    mp["U1"] = 0.5
    mp["eps1"] = -0.24
    sp["model_parameters"] = mp

    # Low-level NRG Parameters
    np = {}
    np["bins"] = 50
    S.set_nrg_params(**np)

    # # Initialize hybridization function
    S.Delta_w['imp'] << 0.1 * SemiCircular(1.0)

    # Solve the impurity model
    S.solve(**sp)
    G_w_1 = S.G_w.copy()

    S.solve(**S.last_solve_params)
    G_w_2 = S.G_w.copy()

    assert_block_gfs_are_close(G_w_1, G_w_2, 1e-12)

if __name__ == '__main__':
    unittest.main()
