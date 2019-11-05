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
#include "./types.hpp"

namespace nrgljubljana_interface {

  /// The collection of all output containers in solver_core
  struct container_set {

    /// Greens function in imaginary time
    g_tau_t G_tau;

    /// Greens function in Matsubara frequencies
    g_iw_t G_iw;

    /// Self-energy in Matsubara frequencies
    g_iw_t Sigma_iw;

    /// Function that writes all containers to hdf5 file
    friend void h5_write(triqs::h5::group h5group, std::string subgroup_name, container_set const &c) {
      auto grp = h5group.create_group(subgroup_name);
      h5_write(grp, "G_tau", c.G_tau);
      h5_write(grp, "G_iw", c.G_iw);
      h5_write(grp, "Sigma_iw", c.Sigma_iw);
    }

    /// Function that reads all containers from hdf5 file
    friend void h5_read(triqs::h5::group h5group, std::string subgroup_name, container_set &c) {
      auto grp = h5group.open_group(subgroup_name);
      h5_read(grp, "G_tau", c.G_tau);
      h5_read(grp, "G_iw", c.G_iw);
      h5_read(grp, "Sigma_iw", c.Sigma_iw);
    }
  };

} // namespace nrgljubljana_interface
