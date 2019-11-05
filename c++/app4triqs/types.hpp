/*******************************************************************************
 *
 * app4triqs: A TRIQS based impurity solver
 *
 * Copyright (c) 2019 The Simons foundation
 *   authors: Nils Wentzell
 *
 * app4triqs is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * app4triqs is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * app4triqs. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include <triqs/gfs.hpp>
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/utility/macros.hpp>

#include <mpi/mpi.hpp>

#include <iostream>
#include <string>
#include <utility>

namespace app4triqs {

  using namespace std::complex_literals; // Complex Unity 1i
  using namespace triqs::gfs;
  using namespace triqs::arrays;
  using namespace triqs::operators;
  using namespace triqs::hilbert_space;
  using namespace triqs::utility;
  using namespace triqs::h5;

  using namespace itertools;

  /// The structure of the gf : block_idx -> pair of block_name and index list (int/string)
  using triqs::hilbert_space::gf_struct_t;

  /// Container type of one-particle Green and Vertex functions in imaginary times
  using g_tau_t = block_gf<imtime, matrix_valued>;

  /// A view to a g_tau_t
  using g_tau_vt = g_tau_t::view_type;

  /// A const_view to a g_tau_t
  using g_tau_cvt = g_tau_t::const_view_type;

  /// Container type of one-particle Green and Vertex functions in Matsubara frequencies
  using g_iw_t = block_gf<imfreq, matrix_valued>;

  /// A view to a g_iw_t
  using g_iw_vt = g_iw_t::view_type;

  /// A const_view to a g_iw_t
  using g_iw_cvt = g_iw_t::const_view_type;

  /// Type of the Monte-Carlo weight. Either double or dcomplex
  using mc_weight_t = g_tau_t::g_t::scalar_t;

  // Declare some placeholders for the rest of the code. Use anonymous namespace for proper linkage
  // in this code, all variables with trailing _ are placeholders by convention.
  namespace {
    triqs::clef::placeholder<0> i_;
    triqs::clef::placeholder<1> j_;
    triqs::clef::placeholder<2> k_;
    triqs::clef::placeholder<3> l_;
    triqs::clef::placeholder<4> iw_;
    triqs::clef::placeholder<5> iw1_;
    triqs::clef::placeholder<6> iw2_;
    triqs::clef::placeholder<7> iw3_;
    triqs::clef::placeholder<8> iw4_;
    triqs::clef::placeholder<9> t_;
    triqs::clef::placeholder<10> t1_;
    triqs::clef::placeholder<11> t2_;
    triqs::clef::placeholder<12> t3_;
    triqs::clef::placeholder<13> bl_;
    triqs::clef::placeholder<14> bl1_;
    triqs::clef::placeholder<15> bl2_;
    triqs::clef::placeholder_prime<0> iW_;
    triqs::clef::placeholder_prime<1> iwp_;
  } // anonymous namespace

} // namespace app4triqs
