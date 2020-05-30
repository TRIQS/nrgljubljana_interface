from nrgljubljana_interface import Solver, Flat
from h5 import *

# Parameters
D, V, U = 1.0, 0.2, 0    # Zero e-e interaction
e_f = -U/2.0             # particle-hole symmetric case
T = 1e-5
omega = 0.2              # phonon frequency
g1 = 0.2                 # electron-phonon coupling strength
n1 = 1                   # n1 in g1*(a+a^dag)*(n-n1)

# Set up the Solver
S = Solver(model = "Holstein/Nph=10", symtype = "QS", mesh_max = 2.0, mesh_min = 1e-5, mesh_ratio = 1.01)

# Solve Parameters
sp = { "T": T, "Lambda": 2.0, "Nz": 4, "Tmin": 1e-6, "keep": 1000, "keepenergy": 10.0, "bandrescale": 1.0 }

# Model Parameters
mp = { "U1": U, "eps1": e_f, "omega": omega, "g1": g1, "n1": n1 }
sp["model_parameters"] = mp

# Initialize hybridization function
S.Delta_w['imp'] << V**2 * Flat(D)

# Solve the impurity model
S.solve(**sp)

# Store the Result
with HDFArchive("holstein_solution.h5", 'w') as arch:
    arch["A_w"] = S.A_w
    arch["G_w"] = S.G_w
    arch["F_w"] = S.F_w
    arch["Sigma_w"] = S.Sigma_w
    arch["expv"] = S.expv

print("<n>=", S.expv["n_d"])           # occupancy
print("<nph>=", S.expv["nph"])         # average phonon number
print("<displ^2>=", S.expv["displ^2"]) # oscillator displacement squared
