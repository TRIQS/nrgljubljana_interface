#pragma once

#include "./point_mesh.hpp"

namespace triqs::mesh {

  using refreq_pts = point_mesh<double>;

} // namespace triqs::mesh

namespace triqs::gfs {

  using mesh::refreq_pts;

} // namespace triqs::gfs
