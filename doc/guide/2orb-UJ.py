from nrgljubljana_interface import Solver, Flat
from h5 import *
import math

# Parameters
D = 1.0 # bandwidth
U = 1.0
Gamma = 0.1
e_f = -0.4
U12 = 1.0
J12 = 0.0
T = 1e-3

# Set up the Solver
S = Solver(model = "2orb-UJ", symtype = "QS", mesh_max = 2.00676, mesh_min = 1e-5, mesh_ratio = 1.01)

# Solve Parameters
sp = { "T": T, "Lambda": 4.0, "Nz": 4, "Tmin": 1e-5, "keep": 4000, "keepenergy": 10.0, "bandrescale": 1.0 }

# Model Parameters
mp = { "U1": U, "U2": U, "eps1": e_f, "eps2": e_f, "U12": U12, "J12": J12 }
sp["model_parameters"] = mp

# Initialize hybridization function
S.Delta_w['imp'] << Gamma * (2.0/math.pi) * Flat(D)

# Solve the impurity model
S.solve(**sp)

# Store the Result
with HDFArchive("2orb-UJ_solution.h5", 'w') as arch:
    arch["A_w"] = S.A_w
    arch["G_w"] = S.G_w
    arch["F_w"] = S.F_w
    arch["Sigma_w"] = S.Sigma_w
    arch["expv"] = S.expv

print("<n>=", S.expv["n_d1"])
print("<n^2>=", S.expv["n_d1^2"])
print("<SS>=", S.expv["S_d1S_d2"])
print("<nn>=", S.expv["n_d1n_d2"])
