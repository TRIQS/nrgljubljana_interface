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
    cp["symtype"] = "QSZ"
    cp["mesh_max"] = 1.0
    cp["mesh_min"] = 1e-14
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
    sp["ops"] = "n_d n_d^2 hop0 A_d sigma_d self_d"
    sp["specs"] = "n_d-n_d"
    sp["specd"] = "A_d-A_d self_d-A_d"
    sp["spect"] = "sigma_d-sigma_d"
    sp["dmnrg"] = True
    sp["fdm"] = True

    # Model Parameters
    mp = {}
    mp["U1"] = 0.2
    mp["eps1"] = -0.08
    mp["B1"] = 0.
    sp["model_parameters"] = mp

    # Low-level NRG Parameters
    #np = {}
    #np["bandrescale"] = 1
    #S.set_nrg_params(**np)

    # # Initialize hybridization function
    # for w in S.Delta_w.mesh: S.Delta_w[w] = 1 / w
    S.Delta_w['up'] << 0.01 * SemiCircularNew(1.0)
    S.Delta_w['dn'] << 0.02 * SemiCircularNew(1.0)

    # Solve the impurity model
    S.solve(**sp)

    # # Store the Result
    with HDFArchive("SIAM.out.h5", 'w') as arch:
        arch["S"] = S

    # # # Compare against reference result
    # # h5diff("SIAM.out.h5", "SIAM.ref.h5")

if __name__ == '__main__':
    unittest.main()
