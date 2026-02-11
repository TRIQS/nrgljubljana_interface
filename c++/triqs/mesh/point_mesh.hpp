// Copyright (c) 2019-2023 Simons Foundation
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You may obtain a copy of the License at
//     https://www.gnu.org/licenses/gpl-3.0.txt
//
// Authors: Rok Zitko, Nils Wentzell

/**
 * @file
 * @brief Provides a point mesh type with arbitrary, sorted mesh points.
 */

#pragma once

#include <triqs/mesh.hpp>
#include <triqs/mesh/utils.hpp>
#include <triqs/mesh/mesh_iterator.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <numeric>
#include <string>
#include <vector>

namespace triqs::mesh {

  /**
   * @ingroup triqs-meshes-other
   * @brief Point mesh type with arbitrary, sorted mesh points.
   *
   * @details A point mesh satisfies the triqs::mesh::Mesh and triqs::mesh::MeshWithValues concepts. It is defined by
   * a sorted vector of values \f$ \{x_0, x_1, \ldots, x_{N-1}\} \f$ with \f$ x_i < x_{i+1} \f$ and has the following
   * properties:
   *
   * - Each mesh point is identified by a unique index \f$ n \in \{0, 1, \ldots, N-1\} \f$.
   * - An index \f$ n \f$ is mapped to the corresponding data index \f$ d \f$ by the identity function \f$ d(n) = n \f$
   * and vice versa.
   * - An index \f$ n \f$ is mapped to its corresponding value \f$ x_n \f$ by the function \f$ x(n) = x_n \f$.
   * - An arbitrary value \f$ x \in [x_0, x_{N-1}] \f$ is mapped to the closest mesh point using binary search.
   *
   * @tparam Value Value type of the mesh points (must be totally ordered).
   */
  template <typename Value>
  class point_mesh {
    static_assert(std::totally_ordered<Value>);

    public:
    /// Value type.
    using value_t = Value;

    /// Index type.
    using index_t = long;

    /// Data index type.
    using data_index_t = long;

    /**
     * @brief %Mesh point of a triqs::mesh::point_mesh.
     *
     * @details It stores the index \f$ n \f$, the data index \f$ d \f$, the hash value of the parent mesh and the
     * value \f$ x_n \f$ of the mesh point.
     *
     * Arithmetic operations are defined for mesh points and scalars of the underlying value type. The operations are
     * performed between the value \f$ x_n \f$ of the mesh point and the given scalar.
     */
    class mesh_point_t {
      public:
      /// Parent mesh type.
      using mesh_t = point_mesh;

      /// Default constructor leaves the mesh point uninitialized.
      mesh_point_t() = default;

      /**
       * @brief Construct a mesh point with a given index \f$ n \f$, data index \f$ d \f$, hash value of the parent
       * mesh and value \f$ x_n \f$.
       *
       * @param n Index \f$ n \f$ of the mesh point.
       * @param d Data index \f$ d \f$ of the mesh point.
       * @param mhash Hash value of the parent mesh.
       * @param x Value \f$ x_n \f$ of the mesh point.
       */
      mesh_point_t(long n, long d, uint64_t mhash, value_t x) : index_(n), data_index_(d), mesh_hash_(mhash), value_(x) {}

      /// Get the index \f$ n \f$ of the mesh point.
      [[nodiscard]] long index() const noexcept { return index_; }

      /// Get the data index \f$ d \f$ of the mesh point.
      [[nodiscard]] long data_index() const noexcept { return data_index_; }

      /// Get the value \f$ x_n \f$ of the mesh point.
      [[nodiscard]] value_t value() const noexcept { return value_; }

      /// Get the hash value of the parent mesh.
      [[nodiscard]] uint64_t mesh_hash() const noexcept { return mesh_hash_; }

      /// Conversion to the value type of the parent mesh.
      operator value_t() const { return value_; }

#define IMPL_OP(OP)                                                                                                                                  \
  /** @brief Binary `OP` operation for a mesh_point_t and some type `U`. */                                                                          \
  template <typename U> friend auto operator OP(mesh_point_t const &mp, U &&y) { return mp.value() OP std::forward<U>(y); }                          \
  /** @brief Binary `OP` operation for some type `U` and a mesh_point_t. */                                                                          \
  template <typename U>                                                                                                                              \
    requires(not std::is_same_v<std::decay_t<U>, mesh_point_t>)                                                                                      \
  friend auto operator OP(U &&x, mesh_point_t const &mp) {                                                                                           \
    return std::forward<U>(x) OP mp.value();                                                                                                         \
  }
      IMPL_OP(+)
      IMPL_OP(-)
      IMPL_OP(*)
      IMPL_OP(/)
#undef IMPL_OP

      private:
      long index_         = 0;
      long data_index_    = 0;
      uint64_t mesh_hash_ = 0;
      value_t value_      = {};
    };

    /**
     * @brief Construct a point mesh from a sorted vector of points.
     *
     * @param pts Sorted vector of mesh point values.
     * @throws triqs::runtime_error if the vector is not sorted.
     */
    point_mesh(std::vector<value_t> pts)
       : pts_(std::move(pts)), mesh_hash_(hash(std::accumulate(pts_.begin(), pts_.end(), value_t{0}), static_cast<long>(pts_.size()))) {
      if (not std::is_sorted(pts_.begin(), pts_.end())) TRIQS_RUNTIME_ERROR << "Point mesh must be constructed with a sorted list of points";
    }

    /**
     * @brief Construct a point mesh from an initializer list of points.
     *
     * @param l Initializer list of mesh point values.
     */
    point_mesh(std::initializer_list<value_t> l) : point_mesh(std::vector<value_t>(l)) {}

    /// Default constructor creates an empty mesh.
    point_mesh() = default;

    /// Equal-to comparison operator compares the mesh points.
    bool operator==(point_mesh const &) const = default;

    /// Not-equal-to comparison operator compares the mesh points.
    bool operator!=(point_mesh const &) const = default;

    /**
     * @brief Check if an index \f$ n \f$ is valid.
     *
     * @param n Index \f$ n \f$ to check.
     * @return True if \f$ 0 \leq n < N \f$, false otherwise.
     */
    [[nodiscard]] bool is_index_valid(index_t n) const noexcept { return 0 <= n and n < size(); }

    /**
     * @brief Check if a value \f$ x \f$ is within the mesh range.
     *
     * @param x Value to check.
     * @return True if \f$ x_0 \leq x \leq x_{N-1} \f$, false otherwise.
     */
    [[nodiscard]] bool is_value_valid(value_t x) const noexcept { return pts_.front() <= x and x <= pts_.back(); }

    /**
     * @brief Map an index \f$ n \in \{0, 1, \ldots, N-1\} \f$ to its corresponding data index \f$ d(n) \f$.
     *
     * @param n Index \f$ n \f$ to map.
     * @return Data index \f$ d(n) = n \f$.
     */
    [[nodiscard]] data_index_t to_data_index(index_t n) const noexcept {
      EXPECTS(is_index_valid(n));
      return n;
    }

    /**
     * @brief Map a value \f$ x \f$ to the closest mesh point and return its data index.
     *
     * @param cmp triqs::mesh::closest_mesh_point_t containing the value \f$ x \f$ to map.
     * @return Data index of the closest mesh point.
     */
    [[nodiscard]] data_index_t to_data_index(closest_mesh_point_t<value_t> const &cmp) const noexcept {
      EXPECTS(is_value_valid(cmp.value));
      return to_data_index(to_index(cmp));
    }

    /**
     * @brief Map a data index \f$ d \in \{0, 1, \ldots, N-1\} \f$ to the corresponding index \f$ n(d) \f$.
     *
     * @param d Data index \f$ d \f$ to map.
     * @return Index \f$ n(d) = d \f$.
     */
    [[nodiscard]] index_t to_index(data_index_t d) const noexcept {
      EXPECTS(is_index_valid(d));
      return d;
    }

    /**
     * @brief Map a value \f$ x \f$ to the closest mesh point and return its index.
     *
     * @details Uses binary search to find the closest mesh point.
     *
     * @param cmp triqs::mesh::closest_mesh_point_t containing the value \f$ x \f$ to map.
     * @return Index of the closest mesh point.
     */
    [[nodiscard]] index_t to_index(closest_mesh_point_t<value_t> const &cmp) const noexcept {
      EXPECTS(is_value_valid(cmp.value));

      auto itr_r = std::lower_bound(pts_.begin(), pts_.end(), cmp.value);
      long i_r   = itr_r - pts_.begin();

      if (i_r == 0) { return 0; }
      if (i_r == size()) { return size() - 1; }

      long i_l = i_r - 1;
      if (std::abs(cmp.value - pts_[i_l]) < std::abs(cmp.value - pts_[i_r]))
        return i_l;
      else
        return i_r;
    }

    /**
     * @brief Subscript operator to access a mesh point by its data index \f$ d \in \{0, 1, \ldots, N-1\} \f$.
     *
     * @param d Data index \f$ d \f$ of the mesh point.
     * @return mesh_point_t with the index \f$ n(d) = d \f$, data index \f$ d \f$, hash value of the current mesh and
     * value \f$ x_d \f$.
     */
    [[nodiscard]] mesh_point_t operator[](long d) const noexcept { return (*this)(d); }

    /**
     * @brief Subscript operator to access a mesh point by a value \f$ x \f$ contained in a
     * triqs::mesh::closest_mesh_point_t.
     *
     * @param cmp triqs::mesh::closest_mesh_point_t containing the value \f$ x \f$.
     * @return mesh_point_t of the closest mesh point.
     */
    [[nodiscard]] mesh_point_t operator[](closest_mesh_point_t<value_t> const &cmp) const noexcept { return (*this)[this->to_data_index(cmp)]; }

    /**
     * @brief Function call operator to access a mesh point by its index \f$ n \in \{0, 1, \ldots, N-1\} \f$.
     *
     * @param n Index \f$ n \f$ of the mesh point.
     * @return mesh_point_t with the index \f$ n \f$, data index \f$ d(n) = n \f$, hash value of the current mesh and
     * value \f$ x_n \f$.
     */
    [[nodiscard]] mesh_point_t operator()(index_t n) const noexcept {
      EXPECTS(is_index_valid(n));
      return {n, n, mesh_hash_, to_value(n)};
    }

    /**
     * @brief Map an index \f$ n \in \{0, 1, \ldots, N-1\} \f$ to its corresponding value \f$ x_n \f$.
     *
     * @param n Index \f$ n \f$ to map.
     * @return Value of the mesh point \f$ x_n \f$.
     */
    [[nodiscard]] value_t to_value(index_t n) const noexcept {
      EXPECTS(is_index_valid(n));
      return pts_[n];
    }

    /// Get the hash value of the mesh.
    [[nodiscard]] uint64_t mesh_hash() const noexcept { return mesh_hash_; }

    /// Get the size \f$ N \f$ of the mesh, i.e. the number of mesh points.
    [[nodiscard]] long size() const noexcept { return static_cast<long>(pts_.size()); }

    /// Get the first index of the mesh, i.e. \f$ 0 \f$.
    [[nodiscard]] static constexpr long first_index() noexcept { return 0; }

    /// Get the last index of the mesh, i.e. \f$ N - 1 \f$.
    [[nodiscard]] long last_index() const noexcept { return size() - 1; }

    /// Get the vector of mesh point values.
    [[nodiscard]] std::vector<value_t> const &points() const noexcept { return pts_; }

    /// Get an iterator to the beginning of the mesh.
    [[nodiscard]] auto begin() const { return mesh_iterator<point_mesh>{.mesh_ptr = this, .data_index = 0}; }

    /// Get a const iterator to the beginning of the mesh.
    [[nodiscard]] auto cbegin() const { return begin(); }

    /// Get an iterator to the end of the mesh.
    [[nodiscard]] auto end() const { return mesh_iterator<point_mesh>{.mesh_ptr = this, .data_index = size()}; }

    /// Get a const iterator to the end of the mesh.
    [[nodiscard]] auto cend() const { return end(); }

    /**
     * @brief Write a triqs::mesh::point_mesh to a `std::ostream`.
     *
     * @param sout `std::ostream` object.
     * @param m %Mesh to be written.
     * @return Reference to `std::ostream` object.
     */
    friend std::ostream &operator<<(std::ostream &sout, point_mesh const &m) { return sout << "Point mesh of size " << m.size(); }

    /**
     * @brief Serialize the mesh to a generic archive.
     * @param ar Archive to serialize to.
     */
    void serialize(auto &ar) const { ar & pts_ & mesh_hash_; }

    /**
     * @brief Deserialize the mesh from a generic archive.
     * @param ar Archive to deserialize from.
     */
    void deserialize(auto &ar) { ar & pts_ & mesh_hash_; }

    /// Get the HDF5 format tag.
    [[nodiscard]] static std::string hdf5_format() { return "PointMesh"; }

    /**
     * @brief Write a triqs::mesh::point_mesh to HDF5.
     *
     * @param g `h5::group` to be written to.
     * @param name Name of the subgroup.
     * @param m %Mesh object to be written.
     */
    friend void h5_write(h5::group g, std::string const &name, point_mesh const &m) {
      h5::group gr = g.create_group(name);
      h5::write_hdf5_format(gr, m);
      h5::write(gr, "points", m.points());
    }

    /**
     * @brief Read a triqs::mesh::point_mesh from HDF5.
     *
     * @param g `h5::group` to be read from.
     * @param name Name of the subgroup.
     * @param m %Mesh object to be read into.
     */
    friend void h5_read(h5::group g, std::string const &name, point_mesh &m) {
      h5::group gr = g.open_group(name);
      h5::assert_hdf5_format(gr, m, true);
      auto pts = h5::read<std::vector<value_t>>(gr, "points");
      m        = point_mesh(pts);
    }

    private:
    std::vector<value_t> pts_;
    uint64_t mesh_hash_ = 0;
  };

  /**
   * @brief Linear interpolation of a function \f$ f \f$ defined on a triqs::mesh::point_mesh at a value \f$ x \f$.
   *
   * @details We calculate
   * \f[
   *   f(x) \approx f_{i_l} * w_l + f_{i_r} * w_r \; ,
   * \f]
   * where \f$ i_l \f$ and \f$ i_r \f$ are the indices of the mesh points bracketing \f$ x \f$, and \f$ w_l \f$ and
   * \f$ w_r \f$ are the interpolation weights based on the distances to the neighboring points.
   *
   * @tparam Value Value type of the mesh.
   * @param m Point mesh.
   * @param f Callable object \f$ f \f$ containing the function values at the mesh points.
   * @param x Value \f$ x \f$ at which to interpolate the function.
   * @return Linear interpolation of \f$ f(x) \f$.
   */
  template <typename Value>
  auto evaluate(point_mesh<Value> const &m, auto const &f, Value x) {
    EXPECTS(m.is_value_valid(x));

    auto itr_r = std::lower_bound(m.points().begin(), m.points().end(), x);
    long i_r   = itr_r - m.points().begin();

    if (i_r == 0) { return f(0); }
    if (i_r == m.size()) { return f(m.size() - 1); }

    long i_l = i_r - 1;

    Value x_l = m.points()[i_l];
    Value x_r = m.points()[i_r];
    auto del  = x_r - x_l;

    // The interpolation weights
    double w_r = (x - x_l) / del;
    double w_l = (x_r - x) / del;

    ASSERT(x_l <= x && x <= x_r);
    ASSERT(w_l + w_r - 1 < 1e-15);

    return f(i_l) * w_l + f(i_r) * w_r;
  }

  // Check mesh concepts.
  static_assert(Mesh<point_mesh<double>>);
  static_assert(MeshWithValues<point_mesh<double>>);

} // namespace triqs::mesh
