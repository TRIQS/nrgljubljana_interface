/*******************************************************************************
 *
 * nrgljubljana_interface: A TRIQS interface to the nrgljubliana impurity solver
 *
 * Copyright (c) 2019 The Simons foundation
 *   authors: Rok Zitko, Nils Wentzell
 *
 * nrgljubljana_interface is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * nrgljubljana_interface is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * nrgljubljana_interface. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once
#include "./types.hpp"

namespace nrgljubljana_interface {

  /// The parameters for the solver construction
  struct constr_params_t {

    /// Path to the template library ("" = bundled library)
    std::string templatedir = ""; // QQQ: store this or not??

    /// Model considered (templated)
    std::string problem = "SIAM";

    /// Mesh maximum frequency
    double mesh_max = 10;

    /// Mesh minimum frequency
    double mesh_min = 1e-4;

    /// Common ratio of the geometric sequence
    double mesh_ratio = 1.05;

    /// Block structure of the gf
    //gf_struct_t gf_struct;

    /// Number of block indeces for the Green function
    [[nodiscard]] int n_blocks() const { return 1; }
    //return gf_struct.size(); }

    // Names of block indeces for the Green function
    [[nodiscard]] auto block_names() const {
      std::vector<std::string> v;
      //for (auto const &bl : gf_struct) v.push_back(bl.first);
      return v;
    }

    /// Write constr_params_t to hdf5
    friend void h5_write(triqs::h5::group h5group, std::string subgroup_name, constr_params_t const &cp);

    /// Read constr_params_t from hdf5
    friend void h5_read(triqs::h5::group h5group, std::string subgroup_name, constr_params_t &cp);
  };

  /// The parameters for the solve function
  struct solve_params_t {

    /// Logarithmic discretization parameter
    double Lambda = 2.0;

    /// Number of discretization meshes
    int Nz = 1;

    /// Lowest scale on the Wilson chain
    double Tmin = 1e-4;

    /// Maximum number of states to keep at each step
    size_t keep = 100;

    /// Cut-off energy for truncation
    double keepenergy = -1.0;

    /// Minimum number of states to keep at each step
    size_t keepmin = 0;

    /// Temperature, k_B T/D,
    double T = 0.001;

    /// Operators to be calculated
    std::string ops = "";

    /// Spectral functions (singlet ops) to compute
    std::string specs = "";

    /// Spectral functions (doublet ops) to compute
    std::string specd = "";

    /// Spectral functions (triplet ops) to compute
    std::string spect = "";

    /// Spectral functions (quadruplet ops) to compute
    std::string specq = "";

    /// Spectral functions (orbital triplet ops) to compute
    std::string specot = "";

    /// Conductance curves to compu<te
    std::string specgt = "";

    /// I_1 curves to compute
    std::string speci1t = "";

    /// I_2 curves to compute
    std::string speci2t = "";

    /// Susceptibilities to compute
    std::string specchit = "";

    /// 3-leg vertex functions to compute?
    std::string specv3 = "";

    /// Compute 3-leg vertex on matsubara/matsubara axis?
    bool v3mm = false;

    /// Perform DMNRG (density-matrix NRG) calculation
    bool dmnrg = false;

    /// Perform CFS (complete Fock space) calculation
    bool cfs = false;

    /// Perform FDM (full-density-matrix) calculation
    bool fdm = false;

    /// Calculate expectation values using FDM
    bool fdmexpv = false;

    /// DMNRG calculation on Matsubara axis
    bool dmnrgmats = false;

    /// FDM calculation on Matsubara axis
    bool fdmmats = false;

    /// Number of Matsubara points to collect
    size_t mats = 100;

    /// Width of logarithmic gaussian
    double alpha = 0.3;

    /// Parameter for Gaussian convolution step
    double gamma = 0.2;

    /// Model parameters
    std::map<std::string, double> model_parameters;

    /// Perform post processing
    bool post_process = true;

    /// Write constr_params_t to hdf5
    friend void h5_write(triqs::h5::group h5group, std::string subgroup_name, solve_params_t const &sp);

    /// Read constr_params_t from hdf5
    friend void h5_read(triqs::h5::group h5group, std::string subgroup_name, solve_params_t &sp);
  };

  /// NRG low-level parameters
  struct nrg_params_t {

    /// Band rescaling factor
    double bandrescale = -1.0;

    /// Number of sites in the star representation
    int mMAX = -1;

    /// Number of sites in the Wilson chain
    int Nmax = -1;

    /// Largest x in the discretization ODE solver
    double xmax = -1.0;

    /// Discretization scheme
    std::string discretization = "Z";

    /// Parameter z in the logarithmic discretization
    double z = 1.0;

    /// Spin-polarized Wilson chain
    bool polarized = false;

    /// 2x2 spin structure in Wilson chain
    bool pol2x2 = false;

    /// Channel-mixing terms in Wilson chain
    bool rungs = false;

    /// Tridiagonalisation approach
    std::string tri = "old";

    /// Precision for tridiagonalisation
    size_t preccpp = 2000;

    /// Eigensolver routine (dsyev|dsyevr|zheev|zheevr|default)
    std::string diag = "default";

    /// Ratio of eigenstates computed in partial diagonalisation
    double diagratio = 1.0;

    /// Minimal matrix size for dsyevr
    size_t dsyevrlimit = 100;

    /// Minimal matrix size for zheevr
    size_t zheevrlimit = 100;

    /// Restart calculation to achieve truncation goal?
    bool restart = true;

    /// Rescale factor for restart=true
    double restartfactor = 2.0;

    /// Additional states to keep in case of a near degeneracy
    double safeguard = 1e-5;

    /// Maximal number of additional states
    size_t safeguardmax = 200;

    /// Threshold value for eigenvalue splitting corrections
    double fixeps = 1e-15;

    /// Parameter \bar{\beta} for thermodynamics
    double betabar = 1.0;

    /// Parameter p for G(T) calculations
    double gtp = 0.7;

    /// Parameter p for chi(T) calculations
    double chitp = 1.0;

    /// Perform Costi-Hewson-Zlatic finite-T calculation
    bool finite = false;

    /// CFS greater correlation function
    bool cfsgt = false;

    /// CFS lesser correlation function
    bool cfsls = false;

    /// FDM greater correlation function?
    bool fdmgt = false;

    /// FDM lesser correlation function?
    bool fdmls = false;

    /// Iteration where we evaluate expv
    size_t fdmexpvn = 0;

    /// T>0 calculation on Matsubara axis
    bool finitemats = false;

    /// Compute density matrixes?
    bool dm = false;

    /// Broadening mesh maximum frequency
    //double broaden_max = 10;

    /// Broadening mesh minimum frequency
    //double broaden_min = -99.;

    /// Auto-tune broaden_min parameter
    double broaden_min_ratio = 3.0;

    /// Common ration of the geometric sequence
    //double broaden_ratio = 1.05;

    /// Smallest energy scale in the problem
    double omega0 = -1.0;

    /// omega0 = omega0_ratio x T
    double omega0_ratio = 1.0;

    /// Diagonalisation threads
    int diagth = 1;

    /// Interleaved diagonalization
    bool substeps = false;

    /// Recalculation strategy
    std::string strategy = "kept";

    /// Initial Wilson chain ops
    size_t Ninit = 0;

    /// Output imaginary parts of correlators?
    bool reim = false;

    /// Number of eigenvalues to dump
    size_t dumpannotated = 0;

    /// Dump in terms of absolute energies
    bool dumpabs = false;

    /// Dump using omega_N energy units
    bool dumpscaled = true;

    /// Dump with # digits of precision
    size_t dumpprecision = 8;

    /// Dump by grouping degenerate states
    bool dumpgroups = true;

    /// Energy tolerance for considering two states as degenerate
    double grouptol = 1e-6;

    /// Dump diagonal matrix elements
    size_t dumpdiagonal = 0;

    /// Save binned (unbroadened) data
    bool savebins = true; // !!

    /// Enable broadening of spectra
    bool broaden = true;

    /// Lower binning limit
    double emin = -1.0;

    /// Upper binning limit
    double emax = -1.0;

    /// bins/decade for spectral data
    size_t bins = 1000;

    /// Shift of the accumulation points for binning
    double accumulation = 0.0;

    /// Bin width for linear mesh
    double linstep = 0;

    /// Peak clipping at the end of the run
    double discard_trim = 1e-16;

    /// Peak clipping on the fly
    double discard_immediately = 1e-16;

    /// Energy window parameter for patching
    double goodE = 2.0;

    /// Do N/N+1 patching?
    bool NN1 = false;

    /// Use even iterations in N/N+2 patching
    bool NN2even = true;

    /// Average over even and odd N/N+2 spectra
    bool NN2avg = false;

    /// a in tanh[a(x-0.5)] window function
    double NNtanh = 0.0;

    /// Widht of columns in 'td'
    size_t width_td = 16;

    /// Width of columns in 'custom'
    size_t width_custom = 16;

    /// Precision of columns in 'td'
    size_t prec_td = 10;

    /// Precision of columns in 'custom'
    size_t prec_custom = 10;

    /// Precision of spectral function output
    size_t prec_xy = 10;

    /// Attempt restart?
    bool resume = false;

    /// list of tokens to define what to log
    std::string log = "";

    /// Log everything
    bool logall = false;

    /// Create DONE file?
    bool done = true;

    /// Perform calculations at 0-th iteration?
    bool calc0 = true;

    /// Keep all states in the last iteratio for DMNRG
    bool lastall = false;

    /// Override automatic lastall setting
    bool lastalloverride = false;

    /// Save detailed subspace info
    bool dumpsubspaces = false;

    /// Dump <f> matrix elements
    bool dump_f = false;

    /// Dump (all) energies to file?
    bool dumpenergies = false;

    /// # of eigenvalues to show for log=e
    size_t logenumber = 10;

    /// Stop calculation at some point?
    std::string stopafter = "";

    /// Stop iteration?
    int forcestop = -1;

    /// Remove temporary data files?
    bool removefiles = true;

    /// Do not output imaginary parts of expvs
    bool noimag = true;

    /// Check operator sumrules
    bool checksumrules = false;

    /// Test diag results
    bool checkdiag = false;

    /// Test tr(rho)=1
    bool checkrho = false;

    /// Write nrg_params_t to hdf5
    friend void h5_write(triqs::h5::group h5group, std::string subgroup_name, nrg_params_t const &sp);

    /// Read nrg_params_t from hdf5
    friend void h5_read(triqs::h5::group h5group, std::string subgroup_name, nrg_params_t &sp);
  };

} // namespace nrgljubljana_interface
