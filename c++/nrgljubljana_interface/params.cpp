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
#include "./params.hpp"

namespace nrgljubljana_interface {

  void h5_write(triqs::h5::group h5group, std::string subgroup_name, constr_params_t const &cp) {
    auto grp = h5group.create_group(subgroup_name);
    h5_write(grp, "problem", cp.problem);
    h5_write(grp, "mesh_max", cp.mesh_max);
    h5_write(grp, "mesh_min", cp.mesh_min);
    h5_write(grp, "mesh_ratio", cp.mesh_ratio);
  }

  void h5_read(triqs::h5::group h5group, std::string subgroup_name, constr_params_t &cp) {
    auto grp = h5group.open_group(subgroup_name);
    h5_read(grp, "problem", cp.problem);
    h5_read(grp, "mesh_max", cp.mesh_max);
    h5_read(grp, "mesh_min", cp.mesh_min);
    h5_read(grp, "mesh_ratio", cp.mesh_ratio);
  }

  void h5_write(triqs::h5::group h5group, std::string subgroup_name, solve_params_t const &sp) {
    auto grp = h5group.create_group(subgroup_name);
    h5_write(grp, "Lambda", sp.Lambda);
    // to do
  }

  void h5_read(triqs::h5::group h5group, std::string subgroup_name, solve_params_t &sp) {
    auto grp = h5group.open_group(subgroup_name);
    h5_read(grp, "Lambda", sp.Lambda);
  }

  void h5_write(triqs::h5::group h5group, std::string subgroup_name, nrg_params_t const &np) {
    auto grp = h5group.create_group(subgroup_name);
    h5_write(grp, "bandrescale", np.bandrescale);
    h5_write(grp, "discretization", np.discretization);
  }

  void h5_read(triqs::h5::group h5group, std::string subgroup_name, nrg_params_t &np) {
    auto grp = h5group.open_group(subgroup_name);
    h5_read(grp, "discretization", np.discretization);
    h5_read(grp, "bandrescale", np.bandrescale);
  }

} // namespace nrgljubljana_interface
