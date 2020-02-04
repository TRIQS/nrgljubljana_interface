#!/usr/bin/env python2

import unittest

from nrgljubljana_interface import Solver, SemiCircular, MeshReFreqPts

from pytriqs.archive import *
from pytriqs.utility.h5diff import h5diff
from pytriqs.utility import mpi
from pytriqs.utility.comparison_tests import *

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

    if mpi.is_master_node(): HDFArchive('Solver.h5', 'w')['S'] = S

    # Rerun the Solver and Compare
    S_old = HDFArchive('Solver.h5', 'r')['S']
    S_old.solve(**S_old.last_solve_params)
    assert_block_gfs_are_close(S_old.G_w, S.G_w, 1e-12)

if __name__ == '__main__':
    unittest.main()
