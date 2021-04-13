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

  /// The collection of all output containers in solver_core
  struct container_set {

    /// The spectral function
    std::optional<g_w_t> A_w;

    /// The spectral function of the auxiliary correlator F_w
    std::optional<g_w_t> B_w;

    /// The retarded Greens function
    std::optional<g_w_t> G_w;

    /// The auxiliary Green function F_w = Sigma_w * G_w
    std::optional<g_w_t> F_w;

    /// The retarded Self energy
    std::optional<g_w_t> Sigma_w;

    /// Expectation values of local impurity operators
    std::map<std::string, double> expv;

    /// Thermodynamic variables (FDM algorithm)
    std::map<std::string, double> tdfdm;

    /// Charge susceptibility
    std::optional<g_w_t> chi_NN_w;

    /// Spin susceptibility
    std::optional<g_w_t> chi_SS_w;

    /// Function that writes all containers to hdf5 file
    friend void h5_write(h5::group h5group, std::string subgroup_name, container_set const &c) {
      auto grp = h5group.create_group(subgroup_name);
      h5_write(grp, "A_w", c.A_w);
      h5_write(grp, "B_w", c.B_w);
      h5_write(grp, "G_w", c.G_w);
      h5_write(grp, "F_w", c.F_w);
      h5_write(grp, "Sigma_w", c.Sigma_w);
      h5_write(grp, "expv", c.expv);
      h5_write(grp, "tdfdm", c.tdfdm);
      h5_write(grp, "chi_NN_w", c.chi_NN_w);
      h5_write(grp, "chi_SS_w", c.chi_SS_w);
    }

    /// Function that reads all containers from hdf5 file
    friend void h5_read(h5::group h5group, std::string subgroup_name, container_set &c) {
      auto grp = h5group.open_group(subgroup_name);
      h5_read(grp, "A_w", c.A_w);
      h5_read(grp, "B_w", c.B_w);
      h5_read(grp, "G_w", c.G_w);
      h5_read(grp, "F_w", c.F_w);
      h5_read(grp, "Sigma_w", c.Sigma_w);
      h5_read(grp, "expv", c.expv);
      h5_read(grp, "tdfdm", c.tdfdm);
      h5_read(grp, "chi_NN_w", c.chi_NN_w);
      h5_read(grp, "chi_SS_w", c.chi_SS_w);
    }
  };

} // namespace nrgljubljana_interface
