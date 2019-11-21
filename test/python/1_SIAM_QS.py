#!/usr/bin/env python2

import unittest

from nrgljubljana_interface import Solver, FlatNew, SemiCircularNew

from pytriqs.gf import *
from pytriqs.archive import *
from pytriqs.operators import *
from pytriqs.utility.h5diff import h5diff


class test_SIAM(unittest.TestCase):

    # Construct Parameters
    cp = {}
    cp["model"] = "SIAM"
    cp["symtype"] = "QS"
    cp["mesh_max"] = 1.0
    cp["mesh_min"] = 1e-6
    cp["mesh_ratio"] = 1.03

    # Set up the Solver
    S = Solver(**cp)

    # Solve Parameters
    sp = {}
    sp["Lambda"] = 2.0
    sp["Nz"] = 2
    sp["Tmin"] = 1e-7
    sp["keep"] = 100
    sp["keepenergy"] = 8.0

    # Model Parameters
    mp = {}
    mp["U1"] = 0.2
    mp["eps1"] = -0.08
    sp["model_parameters"] = mp

    # Low-level NRG Parameters
    #np = {}
    #np["bandrescale"] = 1
    #S.set_nrg_params(**np)

    # # Initialize hybridization function
    S.Delta_w['imp'] << 0.01 * SemiCircularNew(1.0)

    # Solve the impurity model
    S.solve(**sp)

    # # Store the Result
    with HDFArchive("1_SIAM_QS.out.h5", 'w') as arch:
        arch["A_w"] = S.A_w
        arch["G_w"] = S.G_w
        arch["F_w"] = S.F_w
        arch["Sigma_w"] = S.Sigma_w

    # Compare against reference result
    h5diff("1_SIAM_QS.out.h5", "1_SIAM_QS.ref.h5")

if __name__ == '__main__':
    unittest.main()
