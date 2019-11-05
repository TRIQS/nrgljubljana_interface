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

namespace nrgljubljana_interface {

  solver_core::solver_core(constr_params_t const &p) : constr_params(p) {

    // Initialize the non-interacting Green function
    G0_iw = block_gf<imfreq>{{p.beta, Fermion, p.n_iw}, p.gf_struct};

    // Initialize the result containers
    G_tau    = block_gf<imtime>{{p.beta, Fermion, p.n_tau}, p.gf_struct};
    G_iw     = G0_iw;
    Sigma_iw = G0_iw;
  }

  // -------------------------------------------------------------------------------

  void solver_core::solve(solve_params_t const &solve_params) {

    last_solve_params = solve_params;

    if (world.rank() == 0)
      std::cout << "\n"
                   "NRGLJUBLJANA_INTERFACE Solver\n";

    // Assert hermiticity of the given Weiss field
    if (!is_gf_hermitian(G0_iw)) TRIQS_RUNTIME_ERROR << "Please make sure that G0_iw fullfills the hermiticity relation G_ij[iw] = G_ji[-iw]*";

    // Merge constr_params and solve_params
    params_t params(constr_params, solve_params);

    // Reset the results
    container_set::operator=(container_set{});

    // TODO Solve the impurity model

    // Post Processing
    if (params.post_process) { post_process(params); }
  }

  // -------------------------------------------------------------------------------

  void solver_core::post_process(params_t const &p) {

    if (world.rank() == 0)
      std::cout << "\n"
                   "Post-processing ... \n";

    // TODO
  }

} // namespace nrgljubljana_interface
