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

namespace nrgljubljana_interface {

  /// The Solver class
  class solver_core : public container_set {

    private:
    // Mpi Communicator
    mpi::communicator world;

    // Return reference to container_set
    container_set &result_set() { return static_cast<container_set &>(*this); }
    [[nodiscard]] container_set const &result_set() const { return static_cast<container_set const &>(*this); }

    public:
    /**
     * Construct a NRGLJUBLJANA_INTERFACE solver
     *
     * @param construct_parameters Set of parameters specific to the NRGLJUBLJANA_INTERFACE solver
     */
    CPP2PY_ARG_AS_DICT
    explicit solver_core(constr_params_t cp);

    // Delete assignement operator because of const members
    solver_core(solver_core const &s) = default;
    solver_core(solver_core &&s)      = default;
    solver_core &operator=(solver_core const &s) = delete;
    solver_core &operator=(solver_core &&s) = default;
    ~solver_core() = default;

    /**
     * Solve method that performs NRGLJUBLJANA_INTERFACE calculation
     *
     * @param solve_params_t Set of parameters specific to the NRGLJUBLJANA_INTERFACE run
     */
    CPP2PY_ARG_AS_DICT
    void solve(solve_params_t const &solve_params);

    // Read the block structure of Green's function objects from file
    gf_struct_t read_structure(const std::string &filename, bool mandatory);

    // Create a temporary directory for a series of NRG runs
    std::string create_tempdir();

    // Perform an individual NRG calculation. Called from solve()
    void solve_one_z(double z, const std::string &taskdir);

    // Adjust the advanced NRG parameters
    CPP2PY_ARG_AS_DICT
    void set_nrg_params(nrg_params_t const &nrg_params);

    // Establish good defaults for nrg_params
    void set_params();

    // Produce param file for a given value of the twist parameter z.
    void generate_param_file(double z);

    // Struct containing the parameters relevant for the solver construction
    constr_params_t constr_params;

    /// Low-level NRG parameters
    nrg_params_t nrg_params;

    // Struct containing the parameters relevant for the solve process
    std::optional<solve_params_t> last_solve_params;

    /// The Green function structure object
    gf_struct_t gf_struct;

    /// The susceptibility structure object
    gf_struct_t chi_struct;

    /// Logarithmic mesh
    gf_mesh<refreq_pts> log_mesh;

    /// The hybridization function in real frequencies
    g_w_t Delta_w;

    /// Read expectation values
    void readexpv(int Nz);
    
    /// Read a block Green's function (im/re)name-block-ij.dat
    void readGF(const std::string &name, std::optional<g_w_t> &G_w, gf_struct_t &_gf_struct);

    /// Read a block spectral function name-block-ij.dat; here we assume that the
    /// spectral function is purely real.
    void readA(const std::string &name, std::optional<g_w_t> &A_w, gf_struct_t &_gf_struct);

    /// Read a scalar real-valued function name.dat
    // void readc(const std::string &name, std::optional<s_w_t> &s_w); TO DO

    static std::string hdf5_scheme() { return "NRGLJUBLJANA_INTERFACE_SolverCore"; }

    // Function that writes a solver object to hdf5 file
    friend void h5_write(triqs::h5::group h5group, std::string subgroup_name, solver_core const &s);

    // Function that constructs a solver object from an hdf5 file
    CPP2PY_IGNORE
    static solver_core h5_read_construct(triqs::h5::group h5group, std::string subgroup_name);
  };
} // namespace nrgljubljana_interface
