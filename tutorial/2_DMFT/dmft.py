# DMFT(NRG) for Hubbard and Hubbard-like models (including Hubbard-Holstein model)
# RZ, Feb 2020

from __future__ import print_function  # Only needed for Python 2
from pytriqs.gf import *
from pytriqs.operators import *
from pytriqs.archive import *
from pytriqs.utility import mpi
from nrgljubljana_interface import Solver, MeshReFreqPts, hilbert_transform_refreq
import math, os, warnings
import numpy as np
from scipy import interpolate, integrate, special, optimize
from collections import OrderedDict

verbose = True # show info messages

newG = lambda S : S.G_w.copy()                           # Creates a BlockGf object of appropriate structure for the Solver
nr_blocks = lambda bgf : len([bl for bl in bgf.indices]) # Returns the number of blocks in a BlockGf object
block_size = lambda G, bl : len(G[bl].indices[0])        # Matrix size of Green's functions in block 'bl'
identity = lambda G, bl : np.identity(block_size(G,bl))  # Returns the identity matrix in block 'bl'

# Calculate a GF from hybridisation and self-energy
def calc_G(Delta, Sigma, mu):
  G = Delta.copy() # copy structure
  for bl in G.indices:
    for w in G.mesh:
      G[bl][w] = np.linalg.inv( (w+mu)*identity(G, bl) - Delta[bl][w] - Sigma[bl][w] ) # !!!
  return G

# Index range of a GF
index_range = lambda G : range(len(G.indices[0]))

# Return an interpolation-object representation of a spectral function for GF G
def interp_A(G, normalize_to_one = True):
  lx = np.array(list(G.mesh.values()))
  ly = sum( sum( -(1.0/math.pi)*np.array(G[bl].data[:,i,i].imag) for i in index_range(G[bl]) ) # matrix trace
    for bl in G.indices )                                                                      # sum over blocks
  if normalize_to_one:
    nr = sum( sum( 1 for i in index_range(G[bl]) ) for bl in G.indices )                       # number of contributions
    ly = ly/nr                                                                                 # rescale
  return interpolate.interp1d(lx, ly, kind='cubic', bounds_error=False, fill_value=0)

# Calculate occupancy for given hybridisation, self-energy and chemical potential
def calc_occupancy(Delta, Sigma, mu, T):
  mesh_max = max([float(w) for w in Delta.mesh]) # get the extent of the mesh, [-mesh_max:mesh_max]
  Gtrial = calc_G(Delta, Sigma, mu)
  f = interp_A(Gtrial)
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    n = integrate.quad(lambda x : 2*f(x)*special.expit(-x/T), -mesh_max, mesh_max)
  return n[0]

# Calculate the local lattice GF and the hybridisation function for the effective impurity model
# for a given self-energy function, chemical potential and Hilbert transform
def self_consistency(Sigma, mu, ht):
  Gloc = Sigma.copy() # copy structure
  for bl in Gloc.indices:
    for w in Gloc.mesh:
      for i in range(block_size(Gloc, bl)):
        for j in range(block_size(Gloc, bl)): # assuming square matrix
          if i == j:
            Gloc[bl][w][i,i] = ht(w + mu - Sigma[bl][w][i,i]) # Hilbert-transform
          else:
            assert abs(Sigma[bl][w][i,j])<1e-10, "This implementation only supports diagonal self-energy"
            Gloc[bl][w][i,j] = 0.0
  Glocinv = Gloc.inverse()
  Delta = Sigma.copy() # copy structure
  for bl in Delta.indices:
    for w in Delta.mesh:
      Delta[bl][w] = (w+mu)*identity(Delta, bl) - Sigma[bl][w] - Glocinv[bl][w] # !!!
  return Gloc, Delta

# Update mu towards reaching the occupancy goal
def update_mu(Delta, Sigma, occupancy_goal, T, old_mu):
  return optimize.root_scalar(lambda x : calc_occupancy(Delta, Sigma, x, T)-occupancy_goal, x0=old_mu, x1=old_mu-0.1).root

# Iteratively adjust mu, taking into account the self-consistency equation.
# max_mu_adjust = number of cycles for adjusting the value of the chemical potential
# Returns an improved estimate for the hybridisation function.
def adjust_mu(Delta_in, Sigma, occupancy_goal, T, old_mu, ht, max_mu_adjust = 10):
  Delta = Delta_in.copy()
  mu = old_mu
  for _ in range(max_mu_adjust):
    mu = update_mu(Delta, Sigma, occupancy_goal, T, mu)
    Gloc, Delta = self_consistency(Sigma, mu, ht)
  new_mu = mu
  if verbose and mpi.is_master_node(): print("Adjusted mu from %.10f to %.10f" % (old_mu,new_mu))
  return Gloc, Delta, new_mu

# Difference between two Green's functions evaluated as the integrated squared difference between the
# corresponding spectral functions.
def gf_diff(a, b):
  mesh_max = max([float(w) for w in a.mesh]) # mesh interval is [-mesh_max:mesh_max]
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
    f.write("%.18f %.18f %.18f\n" % (w, z.real, z.imag))

# Save all blocks for a block GF to tabulated ASCII files
def save_BlockGf(fn, bgf):
  for bl in bgf.indices:
    save_Gf(fn + "_" + bl + ".dat", bgf[bl])

# Save a spectral function (-1/Pi Im GF) to a tabulated ASCII file
def save_A(fn, gf):
  f = open(fn, "w")
  for w in gf.mesh:
    z = gf[w]
    f.write("%.18f %.18f\n" % (w, -1.0/math.pi*z.imag))

# Save spectral functions for all blocks of the block GF
def save_BlockA(fn, bgf):
  for bl in bgf.indices:
    save_A(fn + "_" + bl + ".dat", bgf[bl])

# Load the minimal set of stored results for restarting the calculation
def load_Sigma_mu(fn):
  with HDFArchive(fn, 'r') as arch:
    Sigma = arch["S"].Sigma_w
    mu = arch["mu"]
  return Sigma, mu

# Initial Sigma and mu from a file
def restart_calculation(fn):
  if verbose and mpi.is_master_node(): print("Starting from stored results in file %s" % fn)
  return load_Sigma_mu(fn) # Load data from an HDF5 archive

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
  return str

# Formatting of the results in stats.dat
def fmt_str(nr_val):
  str = "{:>5}" # itern
  for _ in range(nr_val-1): str += " {:>15.8g}"
  return str

# Adjust Im(Delta) so that the hybridisation strength is not too small for the NRG discretization
# scheme to break down.
def fix_hyb_function(Delta, Delta_min):
  Delta_fixed = Delta.copy()
  for bl in Delta.indices:
    for w in Delta.mesh:
      for n in range(block_size(Delta, bl)): # only diagonal parts
        r = Delta[bl][w][n,n].real
        i = Delta[bl][w][n,n].imag
        Delta_fixed[bl][w][n,n] = r + 1j*(i if i<-Delta_min else -Delta_min)
  # Possible improvement: re-adjust the real part so that the Kramers-Kronig relation is maintained
  return Delta_fixed

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
def nparray_to_bgf(a, S):
  G = newG(S)
  split = np.split(a, nr_blocks(G)) # Here we assume all blocks are equally large
  for i, bl in enumerate(G.indices):
    nparray_to_gf(split[i], G[bl])
  return G

class DMFT_solver(object):
  itern = 0                 # iteration counter
  solution_filename = "solution.h5"
  checkpoint_filename = "checkpoint.h5"
  stats_filename = "stats.dat"

  # Provides the initial Sigma and mu
  def initial_Sigma_mu(self):
    if os.path.isfile(self.solution_filename):     # continue DMFT iteration after failed convergence
      return restart_calculation(self.solution_filename)
    elif os.path.isfile(self.checkpoint_filename): # continue DMFT iteration after interruption
      return restart_calculation(self.checkpoint_filename)
    else:                                     # start from scratch
      return self.new_calculation()

  def __init__(self, param, dmft_param):
    self.param = param
    self.dmft_param = dmft_param
    self.occupancy_goal = param["occupancy"]
    self.T = param["T"]
    self.observables = []                                                      # List of expectation values to compute
    self.cp = {}                                                               # Dictinary with creator parameters
    self.sp = { "T": self.T, "model_parameters": {} }                          # Dictionary with solver parameters
    self.mp = {}                                                               # Dictionary with model parameters
    self.nrgp = {}                                                             # Dictionary with low-level NRG Parameters (optional tweaks)
    # Open file with basic information (convergence criteria, occupancy, expectation values of operators) for monitoring the iteration process
    if mpi.is_master_node():
      self.stats_file = open(self.stats_filename, "w", buffering=1) # line buffered
    # Initialize a function ht0 for calculating the Hilbert transform of the DOS
    if (param["dos"] == "Bethe"):
      self.ht1 = lambda z: 2*(z-1j*np.sign(z.imag)*np.sqrt(1-z**2)) # Analytical expression for Hilbert transform of Bethe lattice DOS
      self.ht0 = lambda z: self.ht1(z/param["Bethe_unit"])
    else:
      table = np.loadtxt(param["dos"])
      self.dosA = Gf(mesh=MeshReFreqPts(table[:,0]), target_shape=[])
      for i, w in enumerate(self.dosA.mesh):
        self.dosA[w] = np.array([[ table[i,1] ]])
      self.ht0 = lambda z: hilbert_transform_refreq(self.dosA, z)

  # Performs the Hilbert transform with argument z, ensuring that the imaginary part is positive (to produce retarded GFs)
  def ht(self, z, EPS = 1e-20): 
    return self.ht0(z.real+1j*(z.imag if z.imag>0.0 else EPS)) # Fix problems due to causality violation

  # This may be called after all constructors (base and all derived classes) have done their job
  def setup_Impurity_Solver(self, mesh_param={}, sp = {}, nrgp = {}, verbose_nrg = False):
    self.cp.update(mesh_param) # dictionary cp is set up by model-specific constructor of the derived class, but we need to add mesh_parameters
    self.S = Solver(**self.cp)
    self.S.set_verbosity(verbose_nrg)
    self.sp.update(sp)     # add remaining solve parameters
    self.nrgp.update(nrgp) # optionally add low-level NRG parameters
    self.S.set_nrg_params(**self.nrgp)
    # Now we initialize the physical quantities
    Sigma, mu = self.initial_Sigma_mu()                                                 # Get initial self-energy and chemical potential
    Gloc, Delta = self_consistency(Sigma, mu, self.ht)                                  # Initialize local GF Gloc and hybridization function Delta
    Gloc, Delta, mu = adjust_mu(Delta, Sigma, self.occupancy_goal, self.T, mu, self.ht) # Adjust mu toward satisfying occupancy goal (and return updated Glocal and Delta as well)
    self.Gloc = Gloc.copy()                                                             # Initial approximation for local lattice GF
    self.Delta = Delta.copy()                                                           # Initial approximation for hybridization function
    self.set_mu(mu)                                                                     # Sets mu and updates the impurity level

  def set_model_parameters(self, mp):
    for k,v in mp.items():
      if v is None: del mp[k] # Remove undefined parameters from the dictionary
    self.sp["model_parameters"].update(mp) # Update, not an assignment!

  # Store the complete set of results
  def store_result(self, fn):
    with HDFArchive(fn, 'w') as arch:
      arch["S"] = self.S
      arch["Gloc"] = self.Gloc
      arch["Delta"] = self.Delta
      arch["mu"] = self.mu
      arch["Gself"] = self.Gself
      arch["occupancy"] = self.occupancy

  # Perform a DMFT step. Input is the hybridization function for solving the effective impurity model,
  # output is a new hybridization function resulting from the application of the DMFT self-consistency equation.
  def dmft_step(self, Delta_in):
    self.itern += 1
    if verbose and mpi.is_master_node():
      print("\nIteration %i min_iter=%i max_iter=%i" % (self.itern, self.dmft_param["min_iter"], self.dmft_param["max_iter"]))
    Delta_in_fixed = fix_hyb_function(Delta_in, self.dmft_param["Delta_min"])
    self.S.Delta_w << Delta_in_fixed

    self.S.solve(**self.sp) # Solve the impurity model

    self.Gself = calc_G(Delta_in_fixed, self.S.Sigma_w, self.mu)                 # impurity GF ("self-energy-trick" improved)
    Gloc_prev = self.Gloc.copy()                                                 # store copy of previous Gloc for checking convergence
    self.Gloc, self.Delta = self_consistency(self.S.Sigma_w, self.mu, self.ht)   # apply the DMFT self-consistency equation
    self.occupancy = calc_occupancy(self.Delta, self.S.Sigma_w, self.mu, self.T) # occupancy calculated from the local lattice GF

    diff_loc_imp = gf_diff(self.Gself, self.Gloc)            # difference between impurity and local lattice GF
    diff_prev = gf_diff(self.Gloc, Gloc_prev)                # difference between two consecutively computed local latice GFs
    diff_occupancy = abs(self.occupancy-self.occupancy_goal) # this difference is used as the measure of deviation

    stats = OrderedDict([("itern", self.itern), ("mu", self.mu), ("diff_loc_imp", diff_loc_imp), ("diff_prev", diff_prev),
                         ("diff_occupancy", diff_occupancy), ("occupancy", self.occupancy)])
    for i in self.observables:
      stats[i] = self.S.expv[i]
    header_string = fmt_str_header(len(stats)).format(*[i for i in stats.keys()])
    stats_string  = fmt_str(len(stats)).format(*[i for i in stats.values()])
    if mpi.is_master_node():
      if self.itern == 1: print(header_string, file = self.stats_file)
      print(stats_string, file = self.stats_file)
      if verbose: 
        print(header_string)
        print(stats_string)

    if self.dmft_param["store_steps"] and mpi.is_master_node():
      os.mkdir(str(self.itern)) # one subdirectory per iteration
      save_BlockGf(str(self.itern)+"/Delta", Delta_in_fixed)
      save_BlockGf(str(self.itern)+"/Sigma", self.S.Sigma_w) # self-energy
      save_BlockGf(str(self.itern)+"/G", self.Gloc)     # local lattice Green's function
      save_BlockA(str(self.itern)+"/A", self.Gloc)      # spectral function of local lattice GF
      self.store_result(str(self.itern)+"/"+self.solution_filename)

    if mpi.is_master_node():
      self.store_result(self.checkpoint_filename) # for checkpoint/restart functionality

    # Check for convergence. The only way to exit the DMFT loop is by generating exceptions.
    if (diff_loc_imp   < self.dmft_param["eps_loc_imp"]   and
        diff_prev      < self.dmft_param["eps_prev"]      and
        diff_occupancy < self.dmft_param["eps_occupancy"] and
        self.itern >= self.dmft_param["min_iter"]):
      if mpi.is_master_node():
        self.store_result(self.solution_filename) # full converged results as an HDF5 file
        os.remove(self.checkpoint_filename)       # checkpoint file is no longer needed
      raise Converged("%s\n%s" % (header_string, stats_string))
    if (self.itern == self.dmft_param["max_iter"]):
      raise FailedToConverge("%s\n%s" % (header_string, stats_string))

    if  self.dmft_param["occup_method"] == "adjust":
      self.Gloc, self.Delta, new_mu = adjust_mu(self.Delta, self.S.Sigma_w, self.occupancy_goal, self.T, self.mu, self.ht) # here we update mu to get closer to target occupancy
      self.set_mu(new_mu)

    return self.Delta

  # DMFT driver routine with linear mixing
  def solve_with_linear_mixing(self, alpha):
    Delta_in = self.Delta.copy()
    while True:
      Delta_out = self.dmft_step(Delta_in)
      newDelta = alpha*Delta_out + (1-alpha)*Delta_in
      Delta_in << newDelta

  # DMFT driver routine with Broyden mixing: the goal is to find a root of the function F(Delta)=dmft_step(Delta)-Delta.
  def solve_with_broyden_mixing(self, alpha):
    F = lambda Delta : self.dmft_step(Delta)-Delta
    npF = lambda x : bgf_to_nparray(F(nparray_to_bgf(x, self.S)))
    xin = bgf_to_nparray(self.Delta)
    optimize.broyden1(npF, xin, alpha=alpha, reduction_method="svd", max_rank=10, verbose=verbose and mpi.is_master_node(), f_tol=1e-99) # Loop forever (small f_tol!)

  def solve(self):
    if (self.dmft_param["mixing_method"] == "linear"):
      self.solve_with_linear_mixing(self.dmft_param["alpha"])
    if (self.dmft_param["mixing_method"] == "broyden"):
      self.solve_with_broyden_mixing(self.dmft_param["alpha"])

class Hubbard_solver(DMFT_solver):
  def __init__(self, param, dmft_param):
    super(Hubbard_solver, self).__init__(param, dmft_param)
    self.cp["model"] = "SIAM"
    self.cp["symtype"] = ("QS" if self.param["B"] is None else "QSZ")
    # Model parameters
    self.mp.update({ "U1": param["U"], "B1": param["B"] })
    self.set_model_parameters(self.mp)
    # Model observables
    self.observables.extend(["n_d", "n_d^2"])
    if param["B"] is not None:
      self.observables.extend(["SZd"])

  def set_mu(self, mu):
    self.mu = mu
    self.sp["model_parameters"]["eps1"] = -mu # important: update the impurity level!
      
  # Initial Sigma and mu for the Hubbard model when starting from scratch
  def new_calculation(self):
    if verbose and mpi.is_master_node(): print("Starting from scratch")
    Sigma = newG(self.S)
    for bl in Sigma.indices:
      for w in Sigma.mesh:
        Sigma[bl][w] = self.param["U"]*self.param["occupancy"]/2.0 # Initialize self-energy with the Hartree shift
    mu = self.param["U"]/2.0                                       # Initial approximation for the chemical potential
    return Sigma, mu

class Holstein_solver(Hubbard_solver): # we are deriving from Hubbard_solver !
  def __init__(self, param, dmft_param):
    super(Holstein_solver, self).__init__(param, dmft_param)
    self.cp["model"] = "Holstein/Np=10"
    # Model parameters
    self.mp.update({ "omega": param["omega"], "g1": param["g"], "n1": param["n"] }) # note: g1 -> g, n1 -> n
    self.set_model_parameters(self.mp)
    # Model observables
    self.observables.extend(["nph", "displ", "displ^2"])
