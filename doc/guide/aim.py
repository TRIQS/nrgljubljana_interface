from nrgljubljana_interface import Solver, Flat
from h5 import *

# Parameters
D, V, U = 1.0, 0.25, 1.0
e_f = -U/2.0 # particle-hole symmetric case
T = 1e-5

# Set up the Solver
S = Solver(model = "SIAM", symtype = "QS", mesh_max = 2.0, mesh_min = 1e-5, mesh_ratio = 1.01)

# Solve Parameters
sp = { "T": T, "Lambda": 2.0, "Nz": 4, "Tmin": 1e-6, "keep": 2000, "keepenergy": 10.0, "bandrescale": 1.0 }

# Model Parameters
mp = { "U1": U, "eps1": e_f }
sp["model_parameters"] = mp

# Initialize hybridization function
S.Delta_w['imp'] << V**2 * Flat(D)

# Solve the impurity model
S.solve(**sp)

# Store the Result
with HDFArchive("aim_solution.h5", 'w') as arch:
    arch["A_w"] = S.A_w
    arch["G_w"] = S.G_w
    arch["F_w"] = S.F_w
    arch["Sigma_w"] = S.Sigma_w
    arch["expv"] = S.expv

print("<n>=", S.expv["n_d"])
