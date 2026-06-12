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

  /// Collection of all output containers held by the solver.
  struct container_set {

    /// The spectral function \f$ A(\omega) \f$.
    std::optional<g_w_t> A_w;

    /// The spectral function \f$ B_l(\omega) \f$ of the auxiliary correlator \f$ F_l(\omega) \f$.
    std::optional<g_w_t> B_l_w;

    /// The spectral function \f$ B_r(\omega) \f$ of the auxiliary correlator \f$ F_r(\omega) \f$.
    std::optional<g_w_t> B_r_w;

    /// The spectral function \f$ C(\omega) \f$ of the auxiliary correlator \f$ I(\omega) \f$.
    std::optional<g_w_t> C_w;

    /// The retarded Green's function \f$ G(\omega) \f$.
    std::optional<g_w_t> G_w;

    /// The auxiliary Green's function \f$ F_l(\omega) = \Sigma(\omega)\, G(\omega) \f$.
    std::optional<g_w_t> F_l_w;

    /// The auxiliary Green's function \f$ F_r(\omega) = G(\omega)\, \Sigma(\omega) \f$.
    std::optional<g_w_t> F_r_w;

    /// The auxiliary Green's function \f$ I(\omega) \f$.
    std::optional<g_w_t> I_w;

    /// Constant Hartree shift to the self-energy, stored as a Green's function.
    std::optional<g_w_t> SigmaHartree_w;

    /// The retarded self-energy \f$ \Sigma(\omega) \f$ (computed from \f$ F_l \f$, \f$ F_r \f$, \f$ G \f$ and \f$ I \f$).
    std::optional<g_w_t> Sigma_w;

    /// Expectation values of local impurity operators.
    std::map<std::string, double> expv;

    /// Thermodynamic variables (FDM algorithm).
    std::map<std::string, double> tdfdm;

    /// Charge susceptibility \f$ \chi_{NN}(\omega) \f$.
    std::optional<g_w_t> chi_NN_w;

    /// Spin susceptibility \f$ \chi_{SS}(\omega) \f$.
    std::optional<g_w_t> chi_SS_w;

    /// Write all containers to an HDF5 file.
    friend void h5_write(h5::group h5group, std::string subgroup_name, container_set const &c) {
      auto grp = h5group.create_group(subgroup_name);
      h5_write(grp, "A_w", c.A_w);
      h5_write(grp, "B_l_w", c.B_l_w);
      h5_write(grp, "B_r_w", c.B_r_w);
      h5_write(grp, "C_w", c.C_w);
      h5_write(grp, "G_w", c.G_w);
      h5_write(grp, "F_l_w", c.F_l_w);
      h5_write(grp, "F_r_w", c.F_r_w);
      h5_write(grp, "I_w", c.I_w);
      h5_write(grp, "SigmaHartree_w", c.SigmaHartree_w);
      h5_write(grp, "Sigma_w", c.Sigma_w);
      h5_write(grp, "expv", c.expv);
      h5_write(grp, "tdfdm", c.tdfdm);
      h5_write(grp, "chi_NN_w", c.chi_NN_w);
      h5_write(grp, "chi_SS_w", c.chi_SS_w);
    }

    /// Read all containers from an HDF5 file.
    friend void h5_read(h5::group h5group, std::string subgroup_name, container_set &c) {
      auto grp = h5group.open_group(subgroup_name);
      h5_read(grp, "A_w", c.A_w);
      h5_read(grp, "B_l_w", c.B_l_w);
      h5_read(grp, "B_r_w", c.B_r_w);
      h5_read(grp, "C_w", c.C_w);
      h5_read(grp, "G_w", c.G_w);
      h5_read(grp, "F_l_w", c.F_l_w);
      h5_read(grp, "F_r_w", c.F_r_w);
      h5_read(grp, "I_w", c.I_w);
      h5_read(grp, "SigmaHartree_w", c.SigmaHartree_w);
      h5_read(grp, "Sigma_w", c.Sigma_w);
      h5_read(grp, "expv", c.expv);
      h5_read(grp, "tdfdm", c.tdfdm);
      h5_read(grp, "chi_NN_w", c.chi_NN_w);
      h5_read(grp, "chi_SS_w", c.chi_SS_w);
    }
  };

} // namespace nrgljubljana_interface
