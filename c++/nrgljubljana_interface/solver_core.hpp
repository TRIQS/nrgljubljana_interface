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
#include "./container_set.hpp"
#include "./params.hpp"
#include "./types.hpp"

#include <iostream>
#include <fstream>
#include <string>

#include <triqs/gfs/hilbert_transform.hpp>
#include <triqs/utility/macros.hpp>

namespace nrgljubljana_interface {

  /**
   * @brief TRIQS interface to the NRGLjubljana numerical renormalization group impurity solver.
   *
   * @details Given a hybridization function, drives the external NRGLjubljana solver to compute
   * spectral functions, Green's functions, the self-energy, susceptibilities and
   * thermodynamic/expectation values for quantum impurity models.
   */
  class solver_core : public container_set {

    private:
    // Mpi Communicator
    mpi::communicator world;

    // Return reference to container_set
    container_set &result_set() { return static_cast<container_set &>(*this); }
    [[nodiscard]] container_set const &result_set() const { return static_cast<container_set const &>(*this); }

    public:
    /**
     * @brief Construct an NRGLjubljana_interface solver.
     *
     * @param cp Construction parameters.
     */
    explicit solver_core(constr_params_t cp);

    // Delete assignement operator because of const members
    solver_core(solver_core const &s)            = default;
    solver_core(solver_core &&s)                 = default;
    solver_core &operator=(solver_core const &s) = delete;
    solver_core &operator=(solver_core &&s)      = default;
    ~solver_core()                               = default;

    /**
     * @brief Solve the impurity problem.
     *
     * @details Runs the full NRGLjubljana pipeline for the assigned hybridization function and
     * populates the result containers.
     *
     * @param solve_params Parameters specific to the NRGLjubljana run.
     */
    void solve(solve_params_t const &solve_params);

    /**
     * @brief Read the block structure of Green's function objects from a file.
     *
     * @param filename File describing the block (gf) structure.
     * @param mandatory If true, a missing file is an error.
     * @return The Green's function block structure.
     */
    gf_struct_t read_structure(const std::string &filename, bool mandatory);

    /**
     * @brief Create a temporary working directory for a series of NRG runs.
     *
     * @param tempdir_ Base directory in which to create the temporary directory.
     * @return Path to the created temporary directory.
     */
    std::string create_tempdir(const std::string &tempdir_);

    /// Write \f$ \Gamma = -\mathrm{Im}\,\Delta(\omega) \f$ to a file.
    void write_gamma();

    /**
     * @brief Prepare the input files for an individual NRG calculation. Called from solve().
     *
     * @param z Discretization twist parameter.
     * @param taskdir Directory in which to prepare the calculation.
     */
    void instantiate(double z, const std::string &taskdir);

    /**
     * @brief Perform an individual NRG calculation. Called from ``solve()``.
     *
     * @param taskdir Directory containing the prepared calculation.
     */
    void solve_one(const std::string &taskdir);

    /**
     * @brief Adjust the advanced (low-level) NRG parameters.
     *
     * @param nrg_params Low-level NRG parameters.
     */
    void set_nrg_params(nrg_params_t const &nrg_params);

    /// Establish good defaults for the low-level NRG parameters.
    void set_params();

    /**
     * @brief Check that all required model parameters have been defined.
     *
     * @param sp Solve parameters to validate.
     */
    void check_model_params(const solve_params_t &sp);

    /**
     * @brief Produce the param file for a given value of the twist parameter \f$ z \f$.
     *
     * @param z Discretization twist parameter.
     */
    void generate_param_file(double z);

    /// Parameters used for the solver construction.
    constr_params_t constr_params;

    /// Low-level NRG parameters.
    nrg_params_t nrg_params;

    /// Parameters used for the most recent solve process.
    std::optional<solve_params_t> last_solve_params;

    /// If true, detailed output from NRGLjubljana and its tools is sent to stdout.
    bool verbose = false;

    /// Keep the temporary directories after the calculation.
    bool keep_temp_dir = false;

    /// The Green's function structure object.
    gf_struct_t gf_struct;

    /// The hybridization function structure object.
    gf_struct_t Delta_struct;

    /// The susceptibility structure object.
    gf_struct_t chi_struct;

    /// Logarithmic real-frequency mesh.
    refreq_log log_mesh;

    /// The hybridization function on the real-frequency axis.
    g_w_t Delta_w;

    /**
     * @brief Read expectation values from the NRG output files.
     *
     * @param Nz Number of discretization twists.
     */
    void readexpv(int Nz);

    /**
     * @brief Read thermodynamic variables (FDM algorithm) from the NRG output files.
     *
     * @param Nz Number of discretization twists.
     */
    void readtdfdm(int Nz);

    /**
     * @brief Read a block Green's function from ``(im/re)name-block-ij.dat`` files.
     *
     * @param name Base name of the Green's function output files.
     * @param G_w Container that receives the Green's function.
     * @param _gf_struct Block structure of the Green's function.
     */
    C2PY_IGNORE void readGF(const std::string &name, std::optional<g_w_t> &G_w, gf_struct_t &_gf_struct);

    /**
     * @brief Read a block spectral function from ``name-block-ij.dat`` files.
     *
     * @details The spectral function is assumed to be purely real.
     *
     * @param name Base name of the spectral function output files.
     * @param A_w Container that receives the spectral function.
     * @param _gf_struct Block structure of the spectral function.
     */
    C2PY_IGNORE void readA(const std::string &name, std::optional<g_w_t> &A_w, gf_struct_t &_gf_struct);

    /// Read a scalar real-valued function name.dat
    // void readc(const std::string &name, std::optional<s_w_t> &s_w); // TO DO

    /// Suppress verbose output from the NRG solver.
    void be_quiet() { verbose = false; }

    /**
     * @brief Set the verbosity (see also ``be_quiet()``).
     *
     * @param v If true, enable verbose output.
     */
    void set_verbosity(bool v) { verbose = v; }

    /// HDF5 format tag for the solver object.
    static std::string hdf5_format() { return "NRGLJUBLJANA_INTERFACE_SolverCore"; }

    /// Write a solver object to an HDF5 file.
    friend void h5_write(h5::group h5group, std::string subgroup_name, solver_core const &s);

    /// Construct a solver object from an HDF5 file.
    C2PY_IGNORE
    static solver_core h5_read_construct(h5::group h5group, std::string subgroup_name);
  };

  /**
   * @brief Hilbert transform for refreq objects (scalar).
   *
   * @param gf Scalar real-frequency Green's function.
   * @param z Complex frequency at which to evaluate the transform.
   * @return Value of the Hilbert transform at \f$ z \f$.
   */
  std::complex<double> hilbert_transform_refreq(const c_w_cvt &gf, std::complex<double> z);

  /**
   * @brief Hilbert transform for refreq objects (matrix, elementwise).
   *
   * @param gf Matrix-valued real-frequency Green's function.
   * @param z Complex frequency at which to evaluate the transform.
   * @return Value of the elementwise Hilbert transform at \f$ z \f$.
   */
  matrix<std::complex<double>> hilbert_transform_elementwise(const m_w_cvt &gf, std::complex<double> z);

} // namespace nrgljubljana_interface
