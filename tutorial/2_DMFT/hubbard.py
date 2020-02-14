# Example of a DMFT(NRG) calculation
# RZ, Feb 2020

import dmft
from pytriqs.utility import mpi

param = {           # Physical parameters:
  "U": 2.0,         #  Hubbard electron-electron repulsion
  "occupancy": 0.8, #  band occupancy (summed over spin, 1.0 is half-filling)
  "T": 1e-4,        #  temperature
  "B": None,        #  magnetic field for spin-polarized calculation with U(1) symmetry; None = non-spin-polarized calculation with SU(2) symmetry
  "dos": "Bethe",   #  "Bethe" for a Bethe lattice (semicircular) DOS or a file name for reading tabulated DOS data from a file;
                    #     format of the DOS file: lines with "eps rho(eps)" pairs; arbitrary bandwidth, but should be within the NRG mesh
  "Bethe_unit": 1.0 #  half-bandwidth in the Bethe lattice DOS for dos="Bethe"; 1.0 for D=1 units, 2.0 for t=1 units
}

dmft_param = {                # DMFT loop parameters:
  "min_iter": 5,              #  minimum number of iterations (prevents premature exit from the DMFT loop)
  "max_iter": 50,             #  maximum number of iterations (signals convergence failure)
  "eps_prev": 1e-5,           #  convergence criterium: integrated diff between two consecutive local spectral functions
  "eps_loc_imp": 1e-5,        #  convergence criterium: integrated diff between local and impurity spectral function
  "eps_occupancy": 1e-4,      #  convergence criterium: difference from the target occupancy
  "mixing_method": "broyden", #  "linear" or "broyden" mixing; the quantity being mixed is the hybridisation function
  "alpha": 0.5,               #  mixing parameter from both linear and Broyden mixing
  "occup_method": "adjust",   #  "adjust" or None; adjust=shift mu so that the current GF meets the occupancy goal
  "Delta_min": 1e-5,          #  minimal value for the hybridisation function; if too large, it produces incorrect spectral distribution,
                              #     if too small, it leads to discretization problems (1e-5 is usually a good and rather safe choice)
  "store_steps": True         #  store intermediate results to files (one subdirectory per iteration)
}

mesh_param = {       # Logarithmic mesh parameters:
  "mesh_max": 10.0,  #  maximum frequency (should be large enough to contain the full spectrum)
  "mesh_min": 1e-5,  #  minimum frequency (should be smaller than the temperature T)
  "mesh_ratio": 1.01 #  common ratio of the mesh (1.01 is usually a good choice)
}

sol_param = {        # Impurity solver parameters:
  "Lambda": 2.0,     #  logarithmic discretization parameter (2.0 is a good choice)
  "Nz": 4,           #  number of interleaved discretization meshes (4 is a good choice for Lambda=2.0)
  "Tmin": 1e-5,      #  lowest temperature/energy scale considered (controls the length of the Wilson chain, choose Tmin<T)
  "keep": 10000,     #  maximum number of states kept in NRG truncation (ensure this number is high enough!)
  "keepenergy": 10.0 #  maximum energy of states kept in NRG truncation (10.0 is a good choice)
}

solver = dmft.Hubbard_solver(param, dmft_param)
solver.setup_Impurity_Solver(mesh_param, sol_param)

try:
  solver.solve()

except dmft.Converged as c:
  if mpi.is_master_node(): 
    print("\nConverged:\n%s" % c.message)
    dmft.save_BlockA("A", solver.Gloc) # converged spectral function for quick plotting

except dmft.FailedToConverge as c:
  if mpi.is_master_node():
    print("\nFailed to converge:\n%s" % c.message) # ... but restart is possible from the checkpoint file

mpi.barrier() # Synchronized exit
