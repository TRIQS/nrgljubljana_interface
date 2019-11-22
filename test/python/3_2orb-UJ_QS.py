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
    cp["model"] = "2orb-UJ"
    cp["symtype"] = "QS"
    cp["mesh_max"] = 1.0
    cp["mesh_min"] = 1e-3
    cp["mesh_ratio"] = 1.1

    # Set up the Solver
    S = Solver(**cp)

    # Solve Parameters
    sp = {}
    sp["Lambda"] = 8.0
    sp["Nz"] = 2
    sp["Tmin"] = 1e-4
    sp["keep"] = 50
    sp["keepenergy"] = 6.0

    # Model Parameters
    mp = {}
    mp["U1"] = 2.0
    mp["U2"] = 2.0
    mp["eps1"] = -0.8
    mp["eps2"] = -0.8
    mp["U12"] = 1.0
    mp["J12"] = 0.5
    sp["model_parameters"] = mp

    # Low-level NRG Parameters
    np = {}
    np["bins"] = 50
    S.set_nrg_params(**np)

    # # Initialize hybridization function
    S.Delta_w['imp'][0,0] << 0.5 * SemiCircularNew(1.0)
    S.Delta_w['imp'][1,1] << 0.5 * SemiCircularNew(1.0)
    # Out-of-diagonal Delta is zero

    # Solve the impurity model
    S.solve(**sp)

    # # Store the Result
    with HDFArchive("3_2orb-UJ_QS.out.h5", 'w') as arch:
        arch["A_w"] = S.A_w
        arch["G_w"] = S.G_w
        arch["F_w"] = S.F_w
        arch["Sigma_w"] = S.Sigma_w

    # Compare against reference result
    h5diff("3_2orb-UJ_QS.out.h5", "3_2orb-UJ_QS.ref.h5")

if __name__ == '__main__':
    unittest.main()
