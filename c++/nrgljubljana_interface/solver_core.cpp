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
#include "./solver_core.hpp"
#include "./post_process.hpp"

#include <nrg-lib.h>

namespace nrgljubljana_interface {

  solver_core::solver_core(constr_params_t const &p) : constr_params(p) {
  }

  // -------------------------------------------------------------------------------

  void solver_core::solve(solve_params_t const &solve_params) {

    last_solve_params = solve_params;

    if (world.rank() == 0)
      std::cout << "\n"
                   "NRGLJUBLJANA_INTERFACE Solver\n";

    // Reset the results
    container_set::operator=(container_set{});

    // TODO Solve the impurity model
    set_workdir(".");
    run_nrg_master();
  }

  void solver_core::set_nrg_params(nrg_params_t const &nrg_params_) {
    nrg_params = nrg_params_;
  }

//  void solver_core::run_single(all_solve_params_t const &all_solve_params) {
//  }

  // -------------------------------------------------------------------------------

    // Function that writes a solver object to hdf5 file

  void h5_write(triqs::h5::group h5group, std::string subgroup_name, solver_core const &s) {
      auto grp = h5group.create_group(subgroup_name);
      h5_write_attribute(grp, "TRIQS_HDF5_data_scheme", solver_core::hdf5_scheme());
      h5_write_attribute(grp, "TRIQS_GIT_HASH", std::string(AS_STRING(TRIQS_GIT_HASH)));
      h5_write_attribute(grp, "NRGLJUBLJANA_INTERFACE_GIT_HASH", std::string(AS_STRING(NRGLJUBLJANA_INTERFACE_GIT_HASH)));
      h5_write(grp, "", s.result_set());
      h5_write(grp, "constr_params", s.constr_params);
      h5_write(grp, "last_solve_params", s.last_solve_params);
      h5_write(grp, "nrg_params", s.nrg_params);
    }

    // Function that constructs a solver object from an hdf5 file
    solver_core solver_core::h5_read_construct(triqs::h5::group h5group, std::string subgroup_name) {
      auto grp           = h5group.open_group(subgroup_name);
      auto constr_params = h5_read<constr_params_t>(grp, "constr_params");
      auto s             = solver_core{constr_params};
      h5_read(grp, "", s.result_set());
      h5_read(grp, "last_solve_params", s.last_solve_params);
      h5_read(grp, "nrg_params", s.nrg_params);
      return s;
    }
} // namespace nrgljubljana_interface
