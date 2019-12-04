from pytriqs.gf import *
from pytriqs.operators import *
from pytriqs.archive import *
from nrgljubljana_interface import Solver, FlatNew

# Parameters
D, V, U = 1.0, 0.2, 4.0
e_f, beta = -U/2.0, 50.0
T = 1.0/beta

# Set up the Solver
S = Solver(model = "SIAM", symtype = "QS", mesh_max = 10.0, mesh_min = 1e-3, mesh_ratio = 1.01)

# Solve Parameters
sp = {}
sp["T"] = T
sp["Lambda"] = 2.0
sp["Nz"] = 2
sp["Tmin"] = 1e-4
sp["keep"] = 1000
sp["keepenergy"] = 10.0

# Model Parameters
mp = {}
mp["U1"] = U
mp["eps1"] = e_f
sp["model_parameters"] = mp

# Low-level NRG Parameters
np = {}
np["bandrescale"] = 1.0
S.set_nrg_params(**np)

# Initialize hybridization function
S.Delta_w['imp'] << V**2 * FlatNew(D)

# Solve the impurity model
S.solve(**sp)

# Store the Result
with HDFArchive("aim_solution.h5", 'w') as arch:
    arch["A_w"] = S.A_w
    arch["G_w"] = S.G_w
    arch["F_w"] = S.F_w
    arch["Sigma_w"] = S.Sigma_w
