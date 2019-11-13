/*******************************************************************************
 *
 * nrgljubljana_interface: A TRIQS based impurity solver
 *
 * Copyright (c) 2019 The Simons foundation
 *   authors: Nils Wentzell
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

namespace nrgljubljana_interface {

  /// The Solver class
  class solver_core : public container_set {

    private:
    // Mpi Communicator
    mpi::communicator world;

    // Return reference to container_set
    container_set &result_set() { return static_cast<container_set &>(*this); }
    container_set const &result_set() const { return static_cast<container_set const &>(*this); }

    public:
    /**
     * Construct a NRGLJUBLJANA_INTERFACE solver
     *
     * @param construct_parameters Set of parameters specific to the NRGLJUBLJANA_INTERFACE solver
     */
    CPP2PY_ARG_AS_DICT
    solver_core(constr_params_t const &constr_params_);

    // Delete assignement operator because of const members
    solver_core(solver_core const &p) = default;
    solver_core(solver_core &&p)      = default;
    ~solver_core()                    = default;
    solver_core &operator=(solver_core const &p) = delete;
    solver_core &operator=(solver_core &&p) = default;

    /**
     * Solve method that performs NRGLJUBLJANA_INTERFACE calculation
     *
     * @param solve_params_t Set of parameters specific to the NRGLJUBLJANA_INTERFACE run
     */
    CPP2PY_ARG_AS_DICT
    void solve(solve_params_t const &solve_params);

    void solve_one_z(solve_params_t const &solve_params, double z);

    CPP2PY_ARG_AS_DICT
    void set_nrg_params(nrg_params_t const &nrg_params);

    void set_params(constr_params_t const &cp,
		    solve_params_t const &sp,
		    nrg_params_t &np);
	
    void generate_param_file(constr_params_t const &cp,
			     solve_params_t const &sp,
			     nrg_params_t const &np);

//    void run_single(all_solve_params_t const &solve_params);

    // Struct containing the parameters relevant for the solver construction
    constr_params_t constr_params;

    // Low-level NRG parameters
    nrg_params_t nrg_params;

    // Struct containing the parameters relevant for the solve process
    std::optional<solve_params_t> last_solve_params;

    static std::string hdf5_scheme() { return "NRGLJUBLJANA_INTERFACE_SolverCore"; }

    // Function that writes a solver object to hdf5 file
    friend void h5_write(triqs::h5::group h5group, std::string subgroup_name, solver_core const &s);

    // Function that constructs a solver object from an hdf5 file
    CPP2PY_IGNORE
    static solver_core h5_read_construct(triqs::h5::group h5group, std::string subgroup_name);
  };
} // namespace nrgljubljana_interface
