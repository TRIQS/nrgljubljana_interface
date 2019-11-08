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
#include "./params.hpp"

namespace nrgljubljana_interface {

  void h5_write(triqs::h5::group h5group, std::string subgroup_name, constr_params_t const &cp) {
    auto grp = h5group.create_group(subgroup_name);
    h5_write(grp, "gf_struct", cp.gf_struct);
  }

  void h5_read(triqs::h5::group h5group, std::string subgroup_name, constr_params_t &cp) {
    auto grp = h5group.open_group(subgroup_name);
    h5_read(grp, "gf_struct", cp.gf_struct);
  }

  void h5_write(triqs::h5::group h5group, std::string subgroup_name, solve_params_t const &sp) {
    auto grp = h5group.create_group(subgroup_name);
    h5_write(grp, "post_process", sp.post_process);
  }

  void h5_read(triqs::h5::group h5group, std::string subgroup_name, solve_params_t &sp) {
    auto grp = h5group.open_group(subgroup_name);
    h5_read(grp, "post_process", sp.post_process);
  }

} // namespace nrgljubljana_interface
