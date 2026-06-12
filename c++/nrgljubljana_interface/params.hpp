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

  /// Construction parameters for the NRGLjubljana solver.
  struct constr_params_t {

    /// Path to the template library (defaults to the bundled templates).
    std::string templatedir = NRGIF_TEMPLATE_DIR;

    /// Impurity model to solve (selects a template directory).
    std::string model = "SIAM";

    /// Symmetry type (NRGLjubljana symmetry code, e.g. QS, QSZ, ISO).
    std::string symtype = "QS";

    /**
     * @brief Resolve the template directory for the chosen model and symmetry.
     *
     * @details Returns ``templatedir/model/symtype``, with the ``$NRGIF_TEMPLATE_DIR``
     * environment variable overriding ``templatedir`` when set.
     *
     * @return Absolute path to the model/symmetry template directory.
     */
    [[nodiscard]] std::string get_model_dir() const {
      if (const char *env_tdir = std::getenv("NRGIF_TEMPLATE_DIR")) {
        return std::string{env_tdir} + "/" + model + "/" + symtype;
      } else {
        return templatedir + "/" + model + "/" + symtype;
      }
    }

    /// Maximum frequency of the logarithmic mesh.
    double mesh_max = 10;

    /// Minimum frequency of the logarithmic mesh.
    double mesh_min = 1e-4;

    /// Common ratio of the geometric (logarithmic) frequency mesh.
    double mesh_ratio = 1.05;

    /// Use a spin-polarized Wilson chain.
    bool polarized = false;

    /// Use a 2x2 spin structure in the Wilson chain.
    bool pol2x2 = false;

    /// Include channel-mixing terms in the Wilson chain.
    bool rungs = false;

    /// Operators whose expectation values are to be calculated.
    std::string ops = "";

    /// Spectral functions of singlet operators to compute.
    std::string specs = "";

    /// Spectral functions of doublet operators to compute.
    std::string specd = "";

    /// Spectral functions of triplet operators to compute.
    std::string spect = "";

    /// Spectral functions of quadruplet operators to compute.
    std::string specq = "";

    /// Spectral functions of orbital-triplet operators to compute.
    std::string specot = "";

    /// Susceptibilities to compute.
    std::string specchit = "";

    /// 3-leg vertex functions to compute.
    std::string specv3 = "";

    /// List of model parameters that need to be specified.
    std::string params = "";

    /// Write constr_params_t to HDF5.
    friend void h5_write(h5::group h5group, std::string subgroup_name, constr_params_t const &cp);

    /// Read constr_params_t from HDF5.
    friend void h5_read(h5::group h5group, std::string subgroup_name, constr_params_t &cp);
  };

  /// Parameters for the solve() method.
  struct solve_params_t {

    /// Logarithmic discretization parameter.
    double Lambda = 2.0;

    /// Number of discretization meshes (interleaved twist parameters z).
    int Nz = 1;

    /// Lowest energy scale on the Wilson chain.
    double Tmin = 1e-4;

    /// Maximum number of states to keep at each NRG step.
    size_t keep = 100;

    /// Cut-off energy for truncation.
    double keepenergy = -1.0;

    /// Minimum number of states to keep at each NRG step.
    size_t keepmin = 0;

    /// Temperature, \f$ k_B T / D \f$.
    double T = 0.001;

    /// Width of the logarithmic gaussian used for broadening.
    double alpha = 0.3;

    /// Parameter for the Gaussian convolution step.
    double gamma = 0.2;

    /// Method for calculating the dynamical quantities.
    std::string method = "fdm";

    /// Band rescaling factor (half-width of the support of the hybridisation function); set to mesh_max if negative.
    double bandrescale = -1.0;

    /// Model parameters (name to value map, e.g. U1, eps1).
    std::map<std::string, double> model_parameters;

    /// Write solve_params_t to HDF5.
    friend void h5_write(h5::group h5group, std::string subgroup_name, solve_params_t const &sp);

    /// Read solve_params_t from HDF5.
    friend void h5_read(h5::group h5group, std::string subgroup_name, solve_params_t &sp);
  };

  /// Low-level NRG parameters.
  struct nrg_params_t {

    /// Perform a DMNRG (density-matrix NRG) calculation.
    bool dmnrg = false;

    /// Perform a CFS (complete Fock space) calculation.
    bool cfs = false;

    /// Perform an FDM (full-density-matrix) calculation.
    bool fdm = true;

    /// Calculate expectation values using the FDM algorithm.
    bool fdmexpv = true;

    /// Perform the DMNRG calculation on the Matsubara axis.
    bool dmnrgmats = false;

    /// Perform the FDM calculation on the Matsubara axis.
    bool fdmmats = false;

    /// Number of Matsubara points to collect.
    size_t mats = 100;

    /// Conductance curves to compute.
    std::string specgt = "";

    /// \f$ I_1 \f$ curves to compute.
    std::string speci1t = "";

    /// \f$ I_2 \f$ curves to compute.
    std::string speci2t = "";

    /// Compute the 3-leg vertex on the Matsubara/Matsubara axis.
    bool v3mm = false;

    /// Number of sites in the star representation (\f$ -1 \f$: automatically determined).
    int mMAX = -1;

    /// Number of sites in the Wilson chain (\f$ -1 \f$: automatically determined).
    int Nmax = -1;

    /// Largest \f$ x \f$ in the discretization ODE solver (\f$ -1 \f$: automatically determined).
    double xmax = -1.0;

    /// Discretization scheme.
    std::string discretization = "Z";

    /// Parameter \f$ z \f$ (twist) in the logarithmic discretization.
    double z = 1.0;

    /// Tridiagonalisation approach.
    std::string tri = "old";

    /// Precision for tridiagonalisation.
    size_t preccpp = 2000;

    /// Eigensolver routine (dsyev|dsyevr|zheev|zheevr|default).
    std::string diag = "default";

    /// Ratio of eigenstates computed in partial diagonalisation.
    double diagratio = 1.0;

    /// Minimal matrix size for dsyevr.
    size_t dsyevrlimit = 100;

    /// Minimal matrix size for zheevr.
    size_t zheevrlimit = 100;

    /// Restart the calculation to achieve the truncation goal.
    bool restart = true;

    /// Rescale factor used when restart is true.
    double restartfactor = 2.0;

    /// Additional states to keep in case of a near degeneracy.
    double safeguard = 1e-5;

    /// Maximal number of additional states to keep.
    size_t safeguardmax = 200;

    /// Threshold value for eigenvalue splitting corrections.
    double fixeps = 1e-15;

    /// Parameter \f$ \bar{\beta} \f$ for thermodynamics.
    double betabar = 1.0;

    /// Parameter \f$ p \f$ for \f$ G(T) \f$ calculations.
    double gtp = 0.7;

    /// Parameter \f$ p \f$ for \f$ \chi(T) \f$ calculations.
    double chitp = 1.0;

    /// Perform a Costi-Hewson-Zlatic finite-T calculation.
    bool finite = false;

    /// Compute the CFS greater correlation function.
    bool cfsgt = false;

    /// Compute the CFS lesser correlation function.
    bool cfsls = false;

    /// Compute the FDM greater correlation function.
    bool fdmgt = false;

    /// Compute the FDM lesser correlation function.
    bool fdmls = false;

    /// Iteration at which the expectation values are evaluated.
    size_t fdmexpvn = 0;

    /// Perform a \f$ T > 0 \f$ calculation on the Matsubara axis.
    bool finitemats = false;

    /// Compute density matrices.
    bool dm = false;

    /// Auto-tune the `broaden_min` parameter.
    double broaden_min_ratio = 3.0;

    /// Smallest energy scale in the problem, \f$ \omega_0 \f$.
    double omega0 = -1.0;

    /// Sets \f$ \omega_0 = \mathtt{omega0\_ratio} \times T \f$.
    double omega0_ratio = 1.0;

    /// Number of diagonalisation threads.
    int diagth = 1;

    /// Use the interleaved diagonalization scheme.
    bool substeps = false;

    /// Recalculation strategy.
    std::string strategy = "kept";

    /// Number of initial Wilson chain operators.
    size_t Ninit = 0;

    /// Output the imaginary parts of the correlators.
    bool reim = false;

    /// Number of eigenvalues to dump.
    size_t dumpannotated = 0;

    /// Dump in terms of absolute energies.
    bool dumpabs = false;

    /// Dump using omega_N energy units.
    bool dumpscaled = true;

    /// Number of digits of precision used when dumping.
    size_t dumpprecision = 8;

    /// Dump by grouping degenerate states.
    bool dumpgroups = true;

    /// Energy tolerance for considering two states as degenerate.
    double grouptol = 1e-6;

    /// Dump diagonal matrix elements.
    size_t dumpdiagonal = 0;

    /// Save binned (unbroadened) data.
    bool savebins = true;

    /// Enable broadening of spectra.
    bool broaden = false;

    /// Lower binning limit.
    double emin = -1.0;

    /// Upper binning limit.
    double emax = -1.0;

    /// Number of bins per decade for spectral data.
    size_t bins = 1000;

    /// Shift of the accumulation points for binning.
    double accumulation = 0.0;

    /// Bin width for the linear mesh.
    double linstep = 0;

    /// Peak clipping at the end of the run.
    double discard_trim = 1e-16;

    /// Peak clipping on the fly.
    double discard_immediately = 1e-16;

    /// Energy window parameter for patching.
    double goodE = 2.0;

    /// Perform N/N+1 patching.
    bool NN1 = false;

    /// Use even iterations in N/N+2 patching.
    bool NN2even = true;

    /// Average over even and odd N/N+2 spectra.
    bool NN2avg = false;

    /// Parameter \f$ a \f$ in the \f$ \tanh[a(x-0.5)] \f$ window function.
    double NNtanh = 0.0;

    /// Width of columns in the 'td' output file.
    size_t width_td = 16;

    /// Width of columns in the 'custom' output file.
    size_t width_custom = 16;

    /// Precision of columns in the 'td' output file.
    size_t prec_td = 10;

    /// Precision of columns in the 'custom' output file.
    size_t prec_custom = 10;

    /// Precision of the spectral function output.
    size_t prec_xy = 10;

    /// Attempt to restart the calculation.
    bool resume = false;

    /// List of tokens defining what to log.
    std::string log = "";

    /// Log everything.
    bool logall = false;

    /// Create a DONE file.
    bool done = true;

    /// Perform calculations at the 0-th iteration.
    bool calc0 = true;

    /// Keep all states in the last iteration for DMNRG.
    bool lastall = false;

    /// Override the automatic lastall setting.
    bool lastalloverride = false;

    /// Save detailed subspace info.
    bool dumpsubspaces = false;

    /// Dump \f$ \langle f \rangle \f$ matrix elements.
    bool dump_f = false;

    /// Dump all energies to a file.
    bool dumpenergies = false;

    /// Number of eigenvalues to show for log=e.
    size_t logenumber = 10;

    /// Stop the calculation at a given point.
    std::string stopafter = "";

    /// Force stop at the given iteration (-1: disabled).
    int forcestop = -1;

    /// Remove temporary data files.
    bool removefiles = true;

    /// Do not output the imaginary parts of expectation values.
    bool noimag = true;

    /// Check operator sum rules.
    bool checksumrules = false;

    /// Test the diagonalisation results.
    bool checkdiag = false;

    /// Test that \f$ \mathrm{tr}(\rho) = 1 \f$.
    bool checkrho = false;

    /// Write nrg_params_t to HDF5.
    friend void h5_write(h5::group h5group, std::string subgroup_name, nrg_params_t const &sp);

    /// Read nrg_params_t from HDF5.
    friend void h5_read(h5::group h5group, std::string subgroup_name, nrg_params_t &sp);
  };

} // namespace nrgljubljana_interface
