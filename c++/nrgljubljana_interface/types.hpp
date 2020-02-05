/*******************************************************************************
 *
 * nrgljubljana_interface: A TRIQS interface to the nrgljubliana impurity solver
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

#include <triqs/gfs.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/utility/macros.hpp>

#include <triqs/gfs/meshes/refreq_pts.hpp>

#include <itertools/itertools.hpp>
#include <mpi/mpi.hpp>

#include <iostream>
#include <string>
#include <utility>

namespace nrgljubljana_interface {

  using namespace std::complex_literals; // Complex Unity 1i
  using namespace triqs::gfs;
  using namespace triqs::arrays;
  using namespace triqs::hilbert_space;
  using namespace triqs::utility;
  using namespace triqs::h5;

  using namespace itertools;

  ///// The structure of the gf : block_idx -> pair of block_name and index list (int/string)
  //using triqs::hilbert_space::gf_struct_t;

  /// Container type
  using g_w_t = block_gf<refreq_pts, matrix_valued>;

  /// A view to a g_w_t
  using g_w_vt = g_w_t::view_type;

  /// A const_view to a g_w_t
  using g_w_cvt = g_w_t::const_view_type;

  /// Container for scalar real-valued quantities with no block structure
  using s_w_t = gf<refreq_pts, scalar_real_valued>;

  /// A view to s_w_t
  using s_w_vt = s_w_t::view_type;

  /// A const_view to s_w_t
  using s_w_cvt = s_w_t::const_view_type;

  /// Container for scalar complex-valued quantities with no block structure
  using c_w_t = gf<refreq_pts, scalar_valued>;

  /// A view to c_w_t
  using c_w_vt = c_w_t::view_type;

  /// A const_view to c_w_t
  using c_w_cvt = c_w_t::const_view_type;

  /// Container for matrix complex-valued quantities with no block structure
  using m_w_t = gf<refreq_pts, matrix_valued>;

  /// A view to m_w_t
  using m_w_vt = m_w_t::view_type;

  /// A const_view to m_w_t
  using m_w_cvt = m_w_t::const_view_type;

  // Declare some placeholders for the rest of the code. Use anonymous namespace for proper linkage
  // in this code, all variables with trailing _ are placeholders by convention.
  namespace {
    triqs::clef::placeholder<0> i_;
    triqs::clef::placeholder<1> j_;
    triqs::clef::placeholder<2> w_;
  } // anonymous namespace

} // namespace nrgljubljana_interface
