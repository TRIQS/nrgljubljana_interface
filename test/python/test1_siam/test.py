#!/usr/bin/env python2

import unittest

from nrgljubljana_interface import Solver

from pytriqs.gf import *
from pytriqs.archive import *
from pytriqs.operators import *
from pytriqs.utility.h5diff import h5diff


class test_SIAM(unittest.TestCase):

    # Construct Parameters
    cp = {}
    cp["problem"] = "SIAM"
    cp["mesh_max"] = 1.0
    cp["mesh_min"] = 1e-4
    cp["mesh_ratio"] = 1.02

    # Set up the Solver
    S = Solver(**cp)

    # Solve Parameters
    sp = {}
    sp["Lambda"] = 2.0
    sp["Nz"] = 1
    sp["Tmin"] = 1e-5
    sp["keep"] = 1000
    sp["keepenergy"] = 8.0
    sp["ops"] = "n_d n_d^2 hop0 A_d sigma_d"
    sp["specs"] = "n_d-n_d"
    sp["specd"] = "A_d-A_d"
    sp["spect"] = "sigma_d-sigma_d"
    sp["dmnrg"] = True
    sp["fdm"] = True

    # Model Parameters
    mp = {}
    mp["U1"] = 0.1
    mp["eps1"] = -0.04
    sp["model_parameters"] = mp

    # Low-level NRG Parameters
    #np = {}
    #np["bandrescale"] = 1
    #S.set_nrg_params(**np)

    # Solve the impurity model
    S.solve(**sp)

    # Store the Result
    with HDFArchive("test.out.h5", 'w') as arch:
        arch["S"] = S

    # -------- Compare ---------
    # h5diff("hubbard.out.h5", "hubbard.ref.h5")


if __name__ == '__main__':
    unittest.main()