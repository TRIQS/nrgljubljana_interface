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

namespace nrgljubljana_interface {

  /// The Solver class
  class solver_core : public container_set {

    private:
    // Mpi Communicator
    mpi::communicator world;

    // Return reference to container_set
    container_set &result_set() { return static_cast<container_set &>(*this); }
    container_set const &result_set() const { return static_cast<container_set const &>(*this); }

    // Function to perform the post-processing steps
    void post_process(params_t const &p);

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

    // Struct containing the parameters relevant for the solver construction
    constr_params_t constr_params;

    // Struct containing the parameters relevant for the solve process
    std::optional<solve_params_t> last_solve_params;

    /// Noninteracting Green Function in Matsubara frequencies
    g_iw_t G0_iw;

    // Allow the user to retrigger post-processing with the last set of parameters
    void post_process() {
      if (not last_solve_params) TRIQS_RUNTIME_ERROR << "You need to run the solver once before you post-process";
      post_process({constr_params, *last_solve_params});
    }

    static std::string hdf5_scheme() { return "NRGLJUBLJANA_INTERFACE_SolverCore"; }

    // Function that writes a solver object to hdf5 file
    friend void h5_write(triqs::h5::group h5group, std::string subgroup_name, solver_core const &s) {
      auto grp = h5group.create_group(subgroup_name);
      h5_write_attribute(grp, "TRIQS_HDF5_data_scheme", solver_core::hdf5_scheme());
      h5_write_attribute(grp, "TRIQS_GIT_HASH", std::string(AS_STRING(TRIQS_GIT_HASH)));
      h5_write_attribute(grp, "NRGLJUBLJANA_INTERFACE_GIT_HASH", std::string(AS_STRING(NRGLJUBLJANA_INTERFACE_GIT_HASH)));
      h5_write(grp, "", s.result_set());
      h5_write(grp, "constr_params", s.constr_params);
      h5_write(grp, "last_solve_params", s.last_solve_params);
      h5_write(grp, "G0_iw", s.G0_iw);
    }

    // Function that constructs a solver object from an hdf5 file
    CPP2PY_IGNORE
    static solver_core h5_read_construct(triqs::h5::group h5group, std::string subgroup_name) {
      auto grp           = h5group.open_group(subgroup_name);
      auto constr_params = h5_read<constr_params_t>(grp, "constr_params");
      auto s             = solver_core{constr_params};
      h5_read(grp, "", s.result_set());
      h5_read(grp, "last_solve_params", s.last_solve_params);
      h5_read(grp, "G0_iw", s.G0_iw);
      return s;
    }
  };
} // namespace nrgljubljana_interface
