# DMFT(NRG) for Hubbard and Hubbard-like models (including Hubbard-Holstein model)
# RZ, Feb 2020

from triqs.gf import *
from triqs.operators import *
from h5 import *
from triqs.utility import mpi
from nrgljubljana_interface import Solver, MeshReFreqPts, hilbert_transform_refreq
import math, os, warnings
import numpy as np
from scipy import interpolate, integrate, special, optimize
from collections import OrderedDict

# Parameters
model = "SIAM"            # impurity model (from the template library)
dos = "Bethe"             # "Bethe" for a Bethe lattice (semicircular) DOS or a file name for reading tabulated DOS data from a file;
                          # format of the DOS file: lines with "eps rho(eps)" pairs; arbitrary bandwidth, but should be within the NRG mesh
Bethe_unit = 1.0          # half-bandwidth in the Bethe lattice DOS for dos="Bethe"; 1.0 for D=1 units, 2.0 for t=1 units
U = 2.0                   # Hubbard electron-electron repulsion
occupancy_goal = 0.8      # band occupancy (summed over spin)
T = 1e-4                  # temperature
B = None                  # magnetic field for spin-polarized calculation with U(1) symmetry;
                          # None = non-spin-polarized calculation with SU(2) symmetry
omega = None              # Holstein phonon frequency
g1 = None                 # electron-phonon coupling strength
n1 = None                 # shift n1 in the Holstein model definition
verbose = True            # show info messages during the iteration
verbose_nrg = False       # show detailed output from the NRG solver
store_steps = True        # store intermediate results to files (one subdirectory per iteration)
min_iter = 5              # minimum number of iterations (prevents premature exit from the DMFT loop)
max_iter = 50             # maximum number of iterations (signals convergence failure)
eps_prev = 1e-5           # convergence criterium: integrated diff between two consecutive local spectral functions
eps_loc_imp = 1e-5        # convergence criterium: integrated diff between local and impurity spectral function
eps_occupancy = 1e-4      # convergence criterium: difference from the target occupancy (occupancy_goal)
mixing_method = "broyden" # "linear" or "broyden" mixing; the quantity being mixed is the hybridisation function
occup_method = "adjust"   # "adjust" or None; adjust=shift mu so that the current GF meets the occupancy goal
alpha = 0.5               # mixing parameter from both linear and Broyden mixing
max_mu_adjust = 10        # number of cycles for adjusting the value of the chemical potential
mesh_max = 10.0           # logarithmic mesh: maximum frequency (should be large enough to contain the full spectrum)
mesh_min = 1e-5           # logarithmic mesh: minimum frequency (should be smaller than the temperature T)
mesh_ratio = 1.01         # logarithmic mesh: common ratio of the mesh (1.01 is usually a good choice)
Delta_min = 1e-5          # minimal value for the hybridisation function; if too large, it produces incorrect spectral distribution,
                          # if too small, it leads to discretization problems (1e-5 is usually a good and rather safe choice)
normalize_to_one = True   # scaling of the spectral function (relevant in the case of block structure and/or matrix structure of GF)
solution_filename = "solution.h5"
checkpoint_filename = "checkpoint.h5"

# Additional quantities of interest
observables = ["n_d", "n_d^2"]
if B is not None:
  observables.extend(["SZd"])
if "Holstein" in model:
  observables.extend(["nph", "displ", "displ^2"])

# Run-time global variables
itern = 0                                          # iteration counter
verbose = verbose and mpi.is_master_node()         # output is only produced by the master node
store_steps = store_steps and mpi.is_master_node() # output is only produced by the master node
symtype = ("QS" if B is None else "QSZ")

# Set up the Solver
S = Solver(model=model, symtype=symtype, mesh_max=mesh_max, mesh_min=mesh_min, mesh_ratio=mesh_ratio)
S.set_verbosity(verbose_nrg)

newG = lambda : S.G_w.copy()                             # Creates a BlockGf object of appropriate structure
nr_blocks = lambda bgf : len([bl for bl in bgf.indices]) # Returns the number of blocks in a BlockGf object
block_size = lambda bl : len(S.G_w[bl].indices[0])       # Matrix size of Green's functions in block 'bl'
identity = lambda bl : np.identity(block_size(bl))       # Returns the identity matrix in block 'bl'

# Solve Parameters
sp = { "T": T, "Lambda": 2.0, "Nz": 4, "Tmin": 1e-5, "keep": 10000, "keepenergy": 10.0 }

# Model Parameters
mp = { "U1": U, "B1": B, "omega": omega, "g1": g1, "n1": n1 }
for k,v in mp.items():
  if v is None: del mp[k] # Remove undefined parameters from the dictionary
sp["model_parameters"] = mp

# Update the chemical potential
def set_mu(new_mu):
    global mu
    mu = new_mu
    sp["model_parameters"]["eps1"] = -mu


# Low-level NRG Parameters (optional tweaks)
nrgp = {}
#nrgp["bandrescale"] = 5.0
S.set_nrg_params(**nrgp)


# Initialize a function ht0 for calculating the Hilbert transform of the DOS
if (dos == "Bethe"):
  ht1 = lambda z: 2*(z-1j*np.sign(z.imag)*np.sqrt(1-z**2)) # Analytical expression for Hilbert transform of Bethe lattice DOS
  global ht0
  ht0 = lambda z: ht1(z/Bethe_unit)
else:
  table = np.loadtxt(dos)
  global dosA
  dosA = Gf(mesh=MeshReFreqPts(table[:,0]), target_shape=[])
  for i, w in enumerate(dosA.mesh):
      dosA[w] = np.array([[ table[i,1] ]])
  ht0 = lambda z: hilbert_transform_refreq(dosA, z)

# Wrapper for ht0 that ensures that the imaginary part of Sigma is negative
EPS = 1e-20 # Cut-off parameter for fixing causality violation
ht = lambda z: ht0(z.real+1j*(z.imag if z.imag>0.0 else EPS)) # Fix problems due to causality violation

# Calculate a GF from hybridisation and self-energy
def calc_G(Delta, Sigma, mu):
  G = newG()
  for bl in G.indices:
    for w in G.mesh:
      G[bl][w] = np.linalg.inv( (w + mu)*identity(bl) - Delta[bl][w] - Sigma[bl][w] ) # !!!
  return G

# Index range of a GF
index_range = lambda G : range(len(G.indices[0]))

# Return an interpolation-object representation of a spectral function for GF G
def interp_A(G):
  lx = np.array(list(G.mesh.values()))
  ly = sum( sum( -(1.0/math.pi)*np.array(G[bl].data[:,i,i].imag) for i in index_range(G[bl]) ) # matrix trace
    for bl in G.indices )                                                                      # sum over blocks
  if normalize_to_one:
    nr = sum( sum( 1 for i in index_range(G[bl]) ) for bl in G.indices )                       # number of contributions
    ly = ly/nr                                                                                 # rescale
  return interpolate.interp1d(lx, ly, kind='cubic', bounds_error=False, fill_value=0)

# Calculate occupancy for given hybridisation, self-energy and chemical potential
def calc_occupancy(Delta, Sigma, mu):
  Gtrial = calc_G(Delta, Sigma, mu)
  f = interp_A(Gtrial)
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    n = integrate.quad(lambda x : 2*f(x)*special.expit(-x/T), -mesh_max, mesh_max)
  return n[0]

# Update mu towards reaching the occupancy goal
def update_mu(Delta, Sigma):
  sol = optimize.root_scalar(lambda x : calc_occupancy(Delta, Sigma, x)-occupancy_goal, x0=mu, x1=mu-0.1)
  set_mu(sol.root)

# Calculate the local lattice GF and the hybridisation function for the effective impurity model
# for a given self-energy function
def self_consistency(Sigma):
  Gloc = newG()
  for bl in Gloc.indices:
    for w in Gloc.mesh:
      for i in range(block_size(bl)):
        for j in range(block_size(bl)): # assuming square matrix
          if i == j:
            Gloc[bl][w][i,i] = ht(w + mu - Sigma[bl][w][i,i]) # Hilbert-transform
          else:
            assert abs(Sigma[bl][w][i,j])<1e-10, "This implementation only supports diagonal self-energy"
            Gloc[bl][w][i,j] = 0.0
  Glocinv = Gloc.inverse()
  Delta = newG()
  for bl in Delta.indices:
    for w in Delta.mesh:
      Delta[bl][w] = (w+mu)*identity(bl) - Sigma[bl][w] - Glocinv[bl][w] # !!!
  return Gloc, Delta

# Iteratively adjust mu, taking into account the self-consistency equation.
# Returns an improved estimate for the hybridisation function.
def adjust_mu(Delta_in, Sigma):
  old_mu = mu
  Delta = Delta_in.copy()
  for _ in range(max_mu_adjust):
    update_mu(Delta, Sigma)
    Gloc, Delta = self_consistency(Sigma)
  new_mu = mu
  if verbose: print("Adjusted mu from %f to %f\n" % (old_mu,new_mu))
  return Gloc, Delta

# Difference between two Green's functions evaluated as the integrated squared difference between the
# corresponding spectral functions.
def gf_diff(a, b):
  f_a = interp_A(a)
  f_b = interp_A(b)
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    diff = integrate.quad(lambda x : (f_a(x)-f_b(x))**2, -mesh_max, mesh_max)
  return diff[0]

# Save a Green's function to a tabulated ASCII file
def save_Gf(fn, gf):
  f = open(fn, "w")
  for w in gf.mesh:
    z = gf[w]
    f.write("%f %f %f\n" % (w, z.real, z.imag))

# Save all blocks for a block GF to tabulated ASCII files
def save_BlockGf(fn, bgf):
  for bl in bgf.indices:
    save_Gf(fn + "_" + bl + ".dat", bgf[bl])

# Save a spectral function (-1/Pi Im GF) to a tabulated ASCII file
def save_A(fn, gf):
  f = open(fn, "w")
  for w in gf.mesh:
    z = gf[w]
    f.write("%f %f\n" % (w, -1.0/math.pi*z.imag))

# Save spectral functions for all blocks of the block GF
def save_BlockA(fn, bgf):
  for bl in bgf.indices:
    save_A(fn + "_" + bl + ".dat", bgf[bl])

# Store the complete set of results
def store_result(fn, S):
  with HDFArchive(fn, 'w') as arch:
    arch["S"] = S
    # Global variables
    arch["Gself"] = Gself
    arch["Gloc"] = Gloc
    arch["mu"] = mu


    
    
    
    

# Load the minimal set of stored results for restarting the calculation
def load_Sigma_mu(fn):
  with HDFArchive(fn, 'r') as arch:
    Sigma = arch["S"].Sigma_w
    mu = arch["mu"]
  return Sigma, mu


# Initial Sigma and mu from a file
def restart_calculation(fn):
  if verbose: print("Starting from stored results in file %s\n" % fn)
  return load_Sigma_mu(fn) # Load data from an HDF5 archive

# Initial Sigma and mu when starting from scratch
def new_calculation():
  if verbose: print("Starting from scratch\n")
  Sigma = newG()
  for bl in Sigma.indices:
    for w in Sigma.mesh:
      Sigma[bl][w] = U*occupancy_goal/2.0  # Initialize self-energy with the Hartree shift
  mu = U/2.0 # initial approximaiton for the chemical potential
  return Sigma, mu

# Provides the initial Sigma and mu
def initial_Sigma_mu():
  if os.path.isfile(solution_filename):     # continue DMFT iteration after failed convergence
    return restart_calculation(solution_filename)
  elif os.path.isfile(checkpoint_filename): # continue DMFT iteration after interruption
    return restart_calculation(checkpoint_filename)
  else:                                     # start from scratch
    return new_calculation()

# Exception to raise when convergence is reached
class Converged(Exception):
  def __init__(self, message):
      self.message = message

# Exception to raise when convergence is not reached (e.g. maximum nr of iterations exceeded)
class FailedToConverge(Exception):
  def __init__(self, message):
      self.message = message

# Formatting of the header in stats.dat
def fmt_str_header(nr_val):
  str = "{:>5}" # itern
  for _ in range(nr_val-1): str += " {:>15}"
  return str + "\n"

# Formatting of the results in stats.dat
def fmt_str(nr_val):
  str = "{:>5}" # itern
  for _ in range(nr_val-1): str += " {:>15.8g}"
  return str + "\n"

# Adjust Im(Delta) so that the hybridisation strength is not too small for the NRG discretization
# scheme to break down.
def fix_hyb_function(Delta, Delta_min):
  Delta_fixed = Delta.copy()
  for bl in Delta.indices:
    for w in Delta.mesh:
      for n in range(block_size(bl)): # only diagonal parts
        r = Delta[bl][w][n,n].real
        i = Delta[bl][w][n,n].imag
        Delta_fixed[bl][w][n,n] = r + 1j*(i if i<-Delta_min else -Delta_min)
  # Possible improvement: re-adjust the real part so that the Kramers-Kronig relation is maintained
  return Delta_fixed

# Perform a DMFT step. Input is the hybridization function for solving the effective impurity model,
# output is a new hybridization function resulting from the application of the DMFT self-consistency equation.
def dmft_step(Delta_in):
  global itern
  itern += 1
  if verbose:
    print("Iteration %i min_iter=%i max_iter=%i\n" % (itern, min_iter, max_iter))
  Delta_in_fixed = fix_hyb_function(Delta_in, Delta_min)
  S.Delta_w << Delta_in_fixed

  S.solve(**sp) # Solve the impurity model
  global Gself, Gloc, Gloc_prev
  Gself = calc_G(Delta_in_fixed, S.Sigma_w, mu) # impurity GF ("self-energy-trick" improved)
  Gloc, Delta = self_consistency(S.Sigma_w)     # apply the DMFT self-consistency equation

  diff_loc_imp = gf_diff(Gself, Gloc)            # difference between impurity and local lattice GF
  diff_prev = gf_diff(Gloc, Gloc_prev)           # difference between two consecutively computed local latice GFs
  Gloc_prev = Gloc.copy()
  occupancy = calc_occupancy(Delta, S.Sigma_w, mu)
  diff_occupancy = abs(occupancy-occupancy_goal) # this difference is used as the measure of deviation

  stats = OrderedDict([("itern", itern), ("mu", mu), ("diff_loc_imp", diff_loc_imp), ("diff_prev", diff_prev),
                       ("diff_occupancy", diff_occupancy), ("occupancy", occupancy)])
  for i in observables:
    stats[i] = S.expv[i]
  header_string = fmt_str_header(len(stats)).format(*[i for i in stats.keys()])
  stats_string  = fmt_str(len(stats)).format(*[i for i in stats.values()])
  if mpi.is_master_node():
    if itern == 1: stats_file.write(header_string)
    stats_file.write(stats_string)
  if verbose: print("stats: %sstats: %s" % (header_string, stats_string))

  if store_steps:
    os.mkdir(str(itern)) # one subdirectory per iteration
    save_BlockGf(str(itern)+"/Delta", Delta_in_fixed)
    save_BlockGf(str(itern)+"/Sigma", S.Sigma_w) # self-energy
    save_BlockGf(str(itern)+"/G", Gloc)          # local lattice Green's function
    save_BlockA(str(itern)+"/A", Gloc)           # spectral function of local lattice GF
    store_result(str(itern)+"/"+solution_filename, S)

  if mpi.is_master_node():
    store_result(checkpoint_filename, S) # for checkpoint/restart functionality

  # Check for convergence. The only way to exit the DMFT loop is by generating exceptions.
  if (diff_loc_imp   < eps_loc_imp   and
      diff_prev      < eps_prev      and
      diff_occupancy < eps_occupancy and
      itern >= min_iter):
    raise Converged(stats_string)
  if (itern == max_iter):
    raise FailedToConverge(stats_string)

  if occup_method == "adjust":
    Gloc, Delta = adjust_mu(Delta, S.Sigma_w) # here we update mu to get closer to target occupancy

  return Delta

# DMFT driver routine with linear mixing
def solve_with_linear_mixing(Delta, alpha):
  Delta_in = Delta.copy()
  while True:
    Delta_out = dmft_step(Delta_in)
    newDelta = alpha*Delta_out + (1-alpha)*Delta_in
    Delta_in << newDelta

# Convert GF object to a linear nparray
def gf_to_nparray(gf):
  return gf.data.flatten()

# Stack elements from all blocks in a block BF
def bgf_to_nparray(bgf):
  return np.hstack((bgf[bl].data.flatten() for bl in bgf.indices))

# Convert a linear numpy array to a single GF object
def nparray_to_gf(a, gf):
  b = a.reshape(gf.data.shape)
  gf.data[:,:,:] = b[:,:,:]

# Extract blocks from linear numpy array to a block GF
def nparray_to_bgf(a):
  G = newG()
  split = np.split(a, nr_blocks(G)) # Here we assume all blocks are equally large
  for i, bl in enumerate(G.indices):
    nparray_to_gf(split[i], G[bl])
  return G

# DMFT driver routine with Broyden mixing: the goal is to find a root of the function F(Delta)=dmft_step(Delta)-Delta.
F = lambda Delta : dmft_step(Delta)-Delta
npF = lambda x : bgf_to_nparray(F(nparray_to_bgf(x)))
def solve_with_broyden_mixing(Delta, alpha):
  xin = bgf_to_nparray(Delta)
  optimize.broyden1(npF, xin, alpha=alpha, reduction_method="svd", max_rank=10, verbose=verbose, f_tol=1e-99) # Loop forever (small f_tol!)

def solve(Delta):
  if (mixing_method=="linear"):
      solve_with_linear_mixing(Delta, alpha)
  if (mixing_method=="broyden"):
      solve_with_broyden_mixing(Delta, alpha)

# Initialization
Sigma0, mu0 = initial_Sigma_mu()       # Get initial self-energy and chemical potential
set_mu(mu0)                            # Set global variable mu and update the impurity level (epsilon_d = -mu)
Gloc, Delta = self_consistency(Sigma0) # Initialize local GF Gloc and hybridization function Delta
Gloc, Delta = adjust_mu(Delta, Sigma0) # Adjust mu (and return updated Glocal and Delta as well)
Gloc_prev = Gloc.copy()                # Store copy for checking convergence

if mpi.is_master_node():
  global stats_file
  stats_file = open("stats.dat", "w", buffering=1) # line buffered

try:
  solve(Delta)

except Converged as c:
  if verbose: print("Converged: %s" % c.message)
  if mpi.is_master_node():
    store_result(solution_filename, S) # full converged results as an HDF5 file
    os.remove(checkpoint_filename)     # checkpoint file is no longer needed
    save_BlockA("A", Gloc)             # converged spectral function for quick plotting

except FailedToConverge as c:
  if verbose: print("Failed to converge: %s" % c.message) # ... but restart is possible from the checkpoint file

mpi.barrier() # Synchronized exit
