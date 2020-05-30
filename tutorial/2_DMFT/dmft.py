# DMFT(NRG) for Hubbard and Hubbard-like models (including Hubbard-Holstein model)
# RZ, Feb 2020

from __future__ import print_function  # Only needed for Python 2

import os
import warnings
from collections import OrderedDict
import timeit
import datetime
import math
import numpy as np
from scipy import interpolate, integrate, special, optimize

from triqs.gf import *
from triqs.operators import *
from h5 import *
from triqs.utility import mpi

from nrgljubljana_interface import Solver, MeshReFreqPts, hilbert_transform_refreq

be_verbose = True # show info messages

def newG(S): return S.G_w.copy()                           # Creates a BlockGf object of appropriate structure for the Solver
def nr_blocks(bgf): return len([bl for bl in bgf.indices]) # Returns the number of blocks in a BlockGf object
def block_size(G): return len(G.indices[0])                # Matrix size of Green's function G
def index_range(G): return range(block_size(G))            # Iterator over matrix indeces
def identity(G): return np.identity(block_size(G))         # Returns the identity matrix of appropriate dimension for Green's function G

# Calculate a GF from hybridisation and self-energy
def calc_G(Delta, Sigma, mu):
  G = Delta.copy() # copy structure
  for bl in G.indices:
    for w in G.mesh:
      G[bl][w] = np.linalg.inv( (w+mu)*identity(G[bl]) - Delta[bl][w] - Sigma[bl][w] ) # !!!
  return G

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
      for i in range(block_size(Gloc[bl])):
        for j in range(block_size(Gloc[bl])): # assuming square matrix
          if i == j:
            Gloc[bl][w][i,i] = ht(w + mu - Sigma[bl][w][i,i]) # Hilbert-transform
          else:
            assert abs(Sigma[bl][w][i,j])<1e-10, "This implementation only supports diagonal self-energy"
            Gloc[bl][w][i,j] = 0.0
  Glocinv = Gloc.inverse()
  Delta = Sigma.copy() # copy structure
  for bl in Delta.indices:
    for w in Delta.mesh:
      Delta[bl][w] = (w+mu)*identity(Delta[bl]) - Sigma[bl][w] - Glocinv[bl][w] # !!!
  return Gloc, Delta

def verbose():
  return be_verbose and mpi.is_master_node()

# Update mu towards reaching the occupancy goal
def update_mu(Delta, Sigma, occupancy_goal, T, old_mu):
  sol = optimize.root_scalar(lambda x : calc_occupancy(Delta, Sigma, x, T)-occupancy_goal, x0=old_mu-0.05, x1=old_mu+0.05)
  if not sol.converged:
    if verbose(): 
      print("root_scalar() failed to converged. Keeping the old value of mu.")
    return old_mu
  return sol.root

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
  with open(fn, "w") as f:
    for w in gf.mesh:
      z = gf[w]
      f.write("%.18f %.18f %.18f\n" % (w, z.real, z.imag))

# Save all blocks for a block GF to tabulated ASCII files
def save_BlockGf(fn, bgf):
  for bl in bgf.indices:
    save_Gf(fn + "_" + bl + ".dat", bgf[bl])

# Save separately re and im parts of a GF object
def save_Gf_re_im(fn, gf):
  with open("re" + fn, "w") as fr, open("im" + fn, "w") as fi:
    for w in gf.mesh:
      z = gf[w]
      fr.write("%.18f %.18f\n" % (w, z.real))
      fi.write("%.18f %.18f\n" % (w, z.imag))

# Save all blocks separately as re and im parts
def save_BlockGf_re_im(fn, bgf):
  for bl in bgf.indices:
    save_Gf_re_im(fn + "_" + bl + ".dat", bgf[bl])

# Save a spectral function (-1/Pi Im GF) to a tabulated ASCII file
def save_A(fn, gf):
  with open(fn, "w") as f:
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
  if verbose(): 
    print("Starting from stored results in file %s" % fn)
  return load_Sigma_mu(fn) # Load data from an HDF5 archive

# Exception to raise when convergence is reached
class Converged(Exception):
  def __init__(self, message):
      self.message = message

# Exception to raise when convergence is not reached (e.g. maximum nr of iterations exceeded)
class FailedToConverge(Exception):
  def __init__(self, message):
      self.message = message
      
# Exception to raise when the calculation is stopped through a flag file
class ForcedStop(Exception):
  def __init__(self):
    pass

# Formatting of the header in stats.dat
def fmt_str_header(nr_val):
  str = "{:>5} {:>9}" # itern, time
  for _ in range(nr_val-2): str += " {:>15}"
  return str

# Formatting of the results in stats.dat
def fmt_str(nr_val):
  str = "{:>5} {:>9}" # itern, time
  for _ in range(nr_val-2): str += " {:>15.10g}"
  return str

# Reconstruct the real part of a GF from its imaginary part using
# the Kramers-Kronig relation.
# Returns a new gf object, input gf is not modified.
def kk_real_from_imag(gf, eps = 1e-10):
  gf_in = gf.copy()
  gf_out = gf.copy()
  for i in range(block_size(gf)):
    for j in range(block_size(gf)):
      if i == j:
        dosA = Gf(mesh=gf.mesh, target_shape=[])
        for w in gf.mesh:
          dosA[w] = np.array([[ (-1.0/math.pi)*gf_in[w][i,i].imag ]]) # dos = -1/pi Im[gf], ignore real part
        for w in gf.mesh:
          gf_out[w][i,i] = hilbert_transform_refreq(dosA, float(w) + eps*1j)
      else:
        gf_out[w][i,j] = 0
        assert abs(gf_in[w][i,j]) < 1e-16, "This implementation is only valid for diagonal matrix case"
  return gf_out

# Adjust Im(Delta) so that the hybridisation strength is not too small for the NRG discretization
# scheme to break down.
def fix_hyb_function(Delta, Delta_min):
  Delta_fixed = Delta.copy()
  for bl in Delta.indices:
    for w in Delta.mesh:
      for n in range(block_size(Delta[bl])): # only diagonal parts
        r = Delta[bl][w][n,n].real
        i = Delta[bl][w][n,n].imag
        Delta_fixed[bl][w][n,n] = r + 1j*(i if i<-Delta_min else -Delta_min)
    Delta_fixed[bl] << kk_real_from_imag(Delta_fixed[bl])
  return Delta_fixed

# Convert GF object to a linear nparray
def gf_to_nparray(gf):
  return gf.data.flatten().imag # only imag parts

# Stack elements from all blocks in a block BF
def bgf_to_nparray(bgf):
  return np.hstack((gf_to_nparray(bgf[bl]) for bl in bgf.indices))

# Convert a linear numpy array to a single GF object
def nparray_to_gf(a, gf):
  b = a.reshape(gf.data.shape)
  gf.data[:,:,:] = b[:,:,:]*1j # only imag parts (b is real)
  gf_fixed = kk_real_from_imag(gf)
  gf << gf_fixed

# Extract blocks from linear numpy array to a block GF
def nparray_to_bgf(a, S):
  G = newG(S)
  split = np.split(a, nr_blocks(G)) # Here we assume all blocks are equally large
  for i, bl in enumerate(G.indices):
    nparray_to_gf(split[i], G[bl])
  return G

# Pack Delta (and optionally a scalar element) to numpy array
def wrap(Delta, val = None):
  x = bgf_to_nparray(Delta)
  if val is not None:
    x = np.append(x, [val])
  return x

# Extract Delta (and optionally a scalar element) from numpy array
def unwrap(x, S, scalar = False):
  if scalar:
    val = x[-1].real
    x = np.delete(x, -1)
  else:
    val = None
  return (nparray_to_bgf(x, S), val) # a tuple

# Touch file
def touch(fname):
  if os.path.exists(fname):
    os.utime(fname, None)
  else:
    open(fname, 'a').close()

# Clamp a value to an interval
def clamp(minimum, x, maximum):
  return max(minimum, min(x, maximum))

class DMFT_solver(object):
  itern = 0                             # iteration counter
  okn = 0                               # number of iterations in which convergent criteria were satisfied
  solution_filename = "solution.h5"     # converged solution (all convergence criteria satisfied)
  checkpoint_filename = "checkpoint.h5" # current solution (for checkpoint/restart functionality)
  stats_filename = "stats.dat"          # some statistics about the DMFT iteration
  stop_flag_file = "STOP"               # the calculation can be aborted by creating this file
  converged_flag_file = "CONVERGED"     # the script touches this file upon reaching convergence

  # Provides the initial Sigma and mu
  def initial_Sigma_mu(self):
    if os.path.isfile(self.checkpoint_filename):  # continue DMFT iteration after interruption
      return restart_calculation(self.checkpoint_filename)
    elif os.path.isfile(self.solution_filename):  # continue DMFT iteration after failed convergence
      return restart_calculation(self.solution_filename)
    else:                                         # start from scratch
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

  # Model specific changes to the impurity model parameters need to be defined in derived classes
  def set_mu(self, mu):
    self.mu = mu

  # This may be called after all constructors (base and all derived classes) have done their job
  def setup_impurity_solver(self, mesh_param={}, sp = {}, nrgp = {}, verbose_nrg = False):
    self.cp.update(mesh_param) # dictionary cp is set up by model-specific constructor of the derived class, but we need to add mesh_parameters
    self.S = Solver(**self.cp)
    self.S.set_verbosity(verbose_nrg)
    self.sp.update(sp)     # add remaining solve parameters
    self.nrgp.update(nrgp) # optionally add low-level NRG parameters
    self.S.set_nrg_params(**self.nrgp)
    # Now we initialize the physical quantities
    Sigma, mu = self.initial_Sigma_mu()                # Get initial self-energy and chemical potential
    Gloc, Delta = self_consistency(Sigma, mu, self.ht) # Initialize local GF Gloc and hybridization function Delta
    Gloc, Delta, mu = self.adjust_mu(Delta, Sigma, mu) # Adjust mu toward satisfying occupancy goal (and return updated Glocal and Delta as well)
    self.Gloc = Gloc.copy()                            # Initial approximation for local lattice GF
    self.Delta = Delta.copy()                          # Initial approximation for hybridization function
    self.set_mu(mu)                                    # Sets mu and updates the impurity level

  # Iteratively adjust mu, taking into account the self-consistency equation.
  # Returns an improved estimate for the hybridisation function.
  def adjust_mu(self, Delta_in, Sigma, old_mu):
    Delta = Delta_in.copy()
    mu = old_mu
    max_mu_adjust = self.dmft_param.get("max_mu_adjust", 10) #  number of cycles for adjusting the value of the chemical potential
    for _ in range(max_mu_adjust):
      mu = update_mu(Delta, Sigma, self.occupancy_goal, self.T, mu)
      Gloc, Delta = self_consistency(Sigma, mu, self.ht)
    new_mu = mu
    max_mu_diff = self.dmft_param.get("max_mu_diff", None)   # maximum change in mu that is allowed
    if max_mu_diff is not None:
      new_mu = clamp(old_mu-max_mu_diff, new_mu, old_mu+max_mu_diff) 
    mu_alpha = self.dmft_param.get("mu_alpha", 1.0)          # mixing parameter for the chemical potential adjustment
    new_mu = mu_alpha*new_mu + (1.0-mu_alpha)*old_mu
    if verbose(): 
      print("Adjusted mu from %.10f to %.10f" % (old_mu,new_mu))
    return Gloc, Delta, new_mu

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
    if verbose():
      print("\nIteration %i min_iter=%i max_iter=%i" % (self.itern, self.dmft_param["min_iter"], self.dmft_param["max_iter"]))
    Delta_in_fixed = fix_hyb_function(Delta_in, self.dmft_param["Delta_min"])
    self.S.Delta_w << Delta_in_fixed

    t0 = timeit.default_timer()
    self.S.solve(**self.sp)                                                  # Solve the impurity model
    t1 = timeit.default_timer()
    dt = int(t1-t0) 
    solver_time = '{:02}:{:02}:{:02}'.format(dt//3600, dt%3600//60, dt%60)   # hh:mm:ss format

    self.Gself = calc_G(Delta_in_fixed, self.S.Sigma_w, self.mu)                 # impurity GF ("self-energy-trick" improved)
    Gloc_prev = self.Gloc.copy()                                                 # store copy of previous Gloc for checking convergence
    self.Gloc, self.Delta = self_consistency(self.S.Sigma_w, self.mu, self.ht)   # apply the DMFT self-consistency equation
    self.occupancy = calc_occupancy(self.Delta, self.S.Sigma_w, self.mu, self.T) # occupancy calculated from the local lattice GF

    diff_loc_imp = gf_diff(self.Gself, self.Gloc)            # difference between impurity and local lattice GF
    diff_prev = gf_diff(self.Gloc, Gloc_prev)                # difference between two consecutively computed local latice GFs
    diff_occupancy = abs(self.occupancy-self.occupancy_goal) # this difference is used as the measure of deviation

    stats = OrderedDict([("itern", self.itern), ("time", solver_time), ("mu", self.mu), ("diff_loc_imp", diff_loc_imp), ("diff_prev", diff_prev),
                         ("diff_occupancy", diff_occupancy), ("occupancy", self.occupancy)])
    for i in self.observables:
      stats[i] = self.S.expv[i]
    header_string = fmt_str_header(len(stats)).format(*[i for i in stats.keys()])
    stats_string  = fmt_str(len(stats)).format(*[i for i in stats.values()])
    if mpi.is_master_node():
      if self.itern == 1: 
        print(header_string, file = self.stats_file)
      print(stats_string, file = self.stats_file)
    if verbose():
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

    # Check for convergence
    if (diff_loc_imp   < self.dmft_param["eps_loc_imp"]   and
        diff_prev      < self.dmft_param["eps_prev"]      and
        diff_occupancy < self.dmft_param["eps_occupancy"]):
      self.okn += 1
    else:
      self.okn = 0

    # The only way to exit the DMFT loop is by generating exceptions.
    if (self.okn >= self.dmft_param.get("conv_iter", 1) and
        self.itern >= self.dmft_param["min_iter"]):
      if mpi.is_master_node():
        self.store_result(self.solution_filename) # full converged results as an HDF5 file
        try:
          os.remove(self.checkpoint_filename)     # checkpoint file is no longer needed
        except OSError:
          pass
        if self.converged_flag_file:
          touch(self.converged_flag_file)         # signal convergence through a file
      raise Converged("%s\n%s" % (header_string, stats_string))
    if (self.itern == self.dmft_param["max_iter"]):
      raise FailedToConverge("%s\n%s" % (header_string, stats_string))
    if self.stop_flag_file and os.path.isfile(self.stop_flag_file):
      raise ForcedStop()

    if self.dmft_param["occup_method"] == "adjust":
      self.Gloc, self.Delta, new_mu = self.adjust_mu(self.Delta, self.S.Sigma_w, self.mu) # here we update mu to get closer to target occupancy
      self.set_mu(new_mu)

    return self.Delta

  # DMFT driver routine with linear mixing
  def solve_with_linear_mixing(self, alpha, mixer_param):
    Delta_in = self.Delta.copy()
    while True:
      Delta_out = self.dmft_step(Delta_in)
      newDelta = alpha*Delta_out + (1-alpha)*Delta_in
      Delta_in << newDelta

  # DMFT driver routine with Broyden mixing: the goal is to find a root of the function F(Delta)=dmft_step(Delta)-Delta.
  def solve_with_broyden_mixing(self, alpha, broyden_param):
    if self.dmft_param["occup_method"] == "broyden":
      def F(Delta_in, mu):
        if verbose(): 
          print("(Broyden) Setting mu to %.18f" % mu)
        self.set_mu(mu)
        Delta_out = self.dmft_step(Delta_in)
        occup_factor = self.dmft_param.get("occup_factor", 1.0)
        return Delta_out-Delta_in, occup_factor*(self.occupancy_goal-self.occupancy)
      npF = lambda x : wrap(*F(*unwrap(x, self.S, scalar=True))) # unpack the tuples
      xin = wrap(self.Delta, self.mu)
    else:
      def F(Delta): return self.dmft_step(Delta)-Delta
      npF = lambda x : bgf_to_nparray(F(nparray_to_bgf(x, self.S)))
      xin = bgf_to_nparray(self.Delta)
    bp = { "reduction_method": "svd", "max_rank": 20, "f_tol": 1e-99 } # Loop forever (small f_tol!)
    bp.update(broyden_param)
    optimize.broyden1(npF, xin, alpha=alpha, verbose=verbose(), **bp)
    raise Converged("broyden1() converged")

  def solve(self):
    alpha = self.dmft_param["alpha"]
    if (self.dmft_param["mixing_method"] == "linear"):
      self.solve_with_linear_mixing(alpha)
    if (self.dmft_param["mixing_method"] == "broyden"):
      self.solve_with_broyden_mixing(alpha, self.dmft_param.get("broyden_param", {}))

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
    super(Hubbard_solver, self).set_mu(mu)
    self.sp["model_parameters"]["eps1"] = -mu # update the impurity level

  # Initial Sigma and mu for the Hubbard model when starting from scratch
  def new_calculation(self):
    if verbose(): 
      print("Starting from scratch")
    Sigma = newG(self.S)
    for bl in Sigma.indices:
      for w in Sigma.mesh:
        Sigma[bl][w] = identity(Sigma[bl]) * self.param["U"]*self.param["occupancy"]/2.0  # Initialize self-energy with the Hartree shift
    mu = self.param["U"]*self.param["occupancy"]/2.0                                      # Initial approximation for the chemical potential
    return Sigma, mu

class Holstein_solver(Hubbard_solver): # we are deriving from Hubbard_solver !
  def __init__(self, param, dmft_param):
    super(Holstein_solver, self).__init__(param, dmft_param)
    self.cp["model"] = "Holstein/Nph=10"
    # Model parameters
    self.mp.update({ "omega": param["omega"], "g1": param["g"], "n1": param["n"] }) # note: g1 -> g, n1 -> n
    self.set_model_parameters(self.mp)
    # Model observables
    self.observables.extend(["nph", "displ", "displ^2"])

  # Initial Sigma and mu for the Hubbard-Holstein model when starting from scratch
  def new_calculation(self):
    if verbose(): 
      print("Starting from scratch")
    Sigma = newG(self.S)
    for bl in Sigma.indices:
      for w in Sigma.mesh:
        Sigma[bl][w] = identity(Sigma[bl]) * (self.param["U"]-2*self.param["g"]**2/self.param["omega"])*self.param["occupancy"]/2.0 # Initialize self-energy with the Hartree shift
    mu = 0 
    return Sigma, mu
