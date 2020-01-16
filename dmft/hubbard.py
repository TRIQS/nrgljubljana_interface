# DMFT(NRG) for Hubbard model on Bethe lattice
# RZ, Jan 2020

from pytriqs.gf import *
from pytriqs.operators import *
from pytriqs.archive import *
from nrgljubljana_interface import Solver # , Omega, MeshReFreqPts
import numpy as np
import math, os, warnings
from scipy import interpolate, integrate, special, optimize

# Parameters
U = 2.0
goal = 0.8 # target occupancy
T = 1e-4
beta = 1/T

# Set up the Solver
S = Solver(model = "SIAM", symtype = "QS", mesh_max = 10.0, mesh_min = 1e-5, mesh_ratio = 1.01)

# Solve Parameters
sp = { "T": T, "Lambda": 2.0, "Nz": 4, "Tmin": 1e-5, "keep": 4000, "keepenergy": 10.0 }

# Model Parameters
mp = { "U1": U }
sp["model_parameters"] = mp

def set_mu(new_mu):
    global mu, old_mu
    try:
      old_mu = mu
    except NameError:
      old_mu = None
    mu = new_mu
    mp["eps1"] = -mu
    sp["model_parameters"] = mp

set_mu(U/2.0) # initial approximaiton for chemical potential

# Low-level NRG Parameters (optional tweaks)
nrgp = {}
#nrgp["bandrescale"] = 1.0
S.set_nrg_params(**nrgp)

def store_result(fn, S):
  with HDFArchive(fn, 'w') as arch:
    arch["A_w"] = S.A_w
    arch["G_w"] = S.G_w
    arch["F_w"] = S.F_w
    arch["Sigma_w"] = S.Sigma_w
    arch["Delta_w"] = S.Delta_w
    arch["expv"] = S.expv

def gf_diff(a, b):
  lx = []
  ly = []
  for w in a.mesh:
    lx.append(float(w))
    ly.append((-1.0/math.pi)*a['imp'][w][0,0].imag - (-1.0/math.pi)*b['imp'][w][0,0].imag) # difference of spectral functions
    # TO DO: matrix case: tr(A-A^dag)
  lx = np.array(lx)
  ly = np.array(ly)
  f = interpolate.interp1d(lx, ly, kind='cubic', bounds_error=False, fill_value=0)
  maxw = S.constr_params["mesh_max"]
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    diff = integrate.quad(lambda x : f(x)**2, -maxw, maxw)
  return diff[0]

# Analytical expression for Hilbert transform of Bethe lattice DOS.
ht0 = lambda z: 2*(z-1j*np.sign(z.imag)*np.sqrt(1-z**2))
EPS = 1e-20
ht = lambda z: ht0(z.real+1j*(z.imag if z.imag>0.0 else EPS)) # Fix problems due to causality violation

def save_gf(fn, gf):
    f = open(fn, "w")
    for w in gf.mesh:
        z = gf[w]
        f.write("%f %f %f\n" % (w, z.real, z.imag))

# Initialize self-energy
for w in S.G_w.mesh:
  S.Sigma_w['imp'][w] = mu

def self_consistency(S):
  Gloc = S.G_w.copy() # Gloc is local GF
  Delta = S.Delta_w.copy()
  for w in Gloc.mesh:
    Gloc['imp'][w] = ht(w + mu - S.Sigma_w['imp'][w]) # "Hilbert-transform" step
  for w in Gloc.mesh:
    Delta['imp'][w] = w + mu - S.Sigma_w['imp'][w] - 1.0/Gloc['imp'][w]
#  Delta << Omega + mu + S.Sigma_w - inverse(Gloc) # not implemented yet
  return Gloc, Delta

# Initilize hybridisation
Gloc, Delta = self_consistency(S)
S.Delta_w << Delta

# Store copies for checking convergence and for mixing
Gloc_prev = Gloc.copy()
Delta_prev = Delta.copy()

verbose = 1 # store intermediate results in each iteration
stats = open("stats.dat", "w", buffering=1) # line buffered
miniter = 5 # minimum number of iterations (prevent premature exit from the loop)
maxiter = 50 # maximum number of iterations
eps_prev = 1e-5 # convergence criterium: diff between two consecutive local GFs
eps_loc_imp = 1e-5 # convergence criterium: diff between local and impurity GF
eps_occupancy = 1e-5 # convergence criterium: occupancy
mixing = 0.5 # alpha parameter in the simple mixing scheme

for itern in range(1,maxiter+1):
  if verbose:
    os.mkdir(str(itern)) # one subdirectory per iteration
    if itern == 1:
      save_gf(str(itern)+"/Sigma_in.dat", S.Sigma_w['imp'])
      save_gf(str(itern)+"/G_in.dat", S.G_w['imp'])
    save_gf(str(itern)+"/Delta.dat", S.Delta_w['imp'])

  S.solve(**sp) # Solve the impurity model

  if verbose:
    save_gf(str(itern)+"/Sigma_out.dat", S.Sigma_w['imp'])
    save_gf(str(itern)+"/G_out.dat", S.G_w['imp'])
    store_result(str(itern)+"/solution.h5", S)

  Gself = S.G_w.copy() # impurity GF (self-energy-trick improved)
  for w in S.G_w.mesh:
    Gself['imp'][w] = 1.0/(w + mu - S.Delta_w['imp'][w] - S.Sigma_w['imp'][w])

  Gloc, Delta = self_consistency(S)

  # Mixing
  Delta1 = Delta.copy()
  Delta2 = Delta_prev.copy()
  newDelta = mixing*Delta1 + (1-mixing)*Delta2
  S.Delta_w << newDelta
  Delta_prev = newDelta.copy()

  if verbose:
    save_gf(str(itern)+"/Gself.dat", Gself['imp'])
    save_gf(str(itern)+"/Gloc.dat", Gloc['imp'])

  diff1 = gf_diff(Gself, Gloc)
  diff2 = gf_diff(Gloc, Gloc_prev)
  Gloc_prev = Gloc.copy()

  def n_vs_mu(x):
    for w in S.G_w.mesh:
      Gself['imp'][w] = 1.0/(w + x - S.Delta_w['imp'][w] - S.Sigma_w['imp'][w])
    lx = []
    ly = []
    for w in Gself.mesh:
      lx.append(float(w))
      ly.append((-1.0/math.pi)*Gself['imp'][w][0,0].imag)
    lx = np.array(lx)
    ly = np.array(ly)
    f = interpolate.interp1d(lx, ly, kind='cubic', bounds_error=False, fill_value=0)
    maxw = S.constr_params["mesh_max"]
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      n = integrate.quad(lambda x : 2*f(x)*special.expit(-beta*x), -maxw, maxw)
    return n[0]

  occupancy1 = S.expv["n_d"]
  occupancy2 = n_vs_mu(mu)

  sol = optimize.root_scalar(lambda x : n_vs_mu(x)-goal, x0=mu, x1=mu-0.1)
  print("sol=", sol)
  set_mu(sol.root)

  stats.write("%i %f %f %f %f %f %f\n" % (itern, diff1, diff2, old_mu, occupancy1, occupancy2, mu))
  if (diff1 < eps_loc_imp and diff2 < eps_prev and abs(occupancy2-goal) < eps_occupancy and itern >= miniter):
    break

store_result("solution.h5", S) # to do: Gself, Gloc
