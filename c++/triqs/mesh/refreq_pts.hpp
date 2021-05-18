#pragma once

#include "./point_mesh.hpp"

namespace triqs::mesh {

  struct refreq_pts : point_mesh<R_domain> {

    // The point mesh base class
    using point_mesh_t = point_mesh<R_domain>;
    using point_mesh_t::point_mesh_t;

    using point_mesh<R_domain>::domain_t;

    using var_t = refreq_pts;

    static std::string hdf5_format() { return "MeshReFreqPts"; }

    friend void h5_write(h5::group fg, std::string const &subgroup_name, refreq_pts const &m) {
      h5_write_impl(fg, subgroup_name, m, hdf5_format().c_str());
    }

    friend void h5_read(h5::group fg, std::string const &subgroup_name, refreq_pts &m) { h5_read_impl(fg, subgroup_name, m, hdf5_format().c_str()); }
  };

} // namespace triqs::mesh
