#include <gtest/gtest.h>

#include <triqs/gfs/hilbert_transform.hpp>
#include <triqs/gfs/meshes/refreq_pts.hpp>

#include <triqs/test_tools/gfs.hpp>

using namespace triqs::gfs;
using namespace triqs::arrays;
using namespace triqs::clef;

using namespace std::literals::complex_literals;

// Common globals
const auto w_     = placeholder<0>{};
const auto bl_    = placeholder<1>{};
const auto w_mesh = gf_mesh<refreq_pts>{-1.0, -0.5, 0.0, 0.5, 1.0};

TEST(hilbert_transform, gf_point) {

  // Construction
  auto rho = gf<refreq_pts, scalar_valued>{w_mesh, {}};

  // Compare against expected values
  EXPECT_EQ(hilbert_transform(rho, 1.0), 1.0);
  EXPECT_EQ(hilbert_transform(rho, 1i), 1.0);
}

TEST(hilbert_transform, gf_mesh) {

  // Construction
  auto rho = gf<refreq_pts, scalar_valued>{w_mesh, {}};

  // Initialize
  rho[w_] << 0.;

  // Full hilbert transform
  auto g = hilbert_transform(rho, w_mesh);

  // Compare against expected values
  EXPECT_EQ(g[0], 1.0);
}

TEST(hilbert_transform, block_gf_point) {

  // Construction
  auto rho = block_gf<refreq_pts, scalar_valued>(w_mesh, {{"up", {}}, {"dn", {}}});

  // Initialize
  rho[bl_][w_] << 0.;

  // Full hilbert transform
  auto vals = hilbert_transform(rho, 1.0);

  // Compare against expected values
  EXPECT_EQ(vals[0], 1.0);
  EXPECT_EQ(vals[1], 1.0);
}

TEST(hilbert_transform, block_gf_mesh) {

  // Construction
  auto rho = block_gf<refreq_pts, scalar_valued>(w_mesh, {{"up", {}}, {"dn", {}}});

  // Initialize
  rho[bl_][w_] << 0.;

  // Full hilbert transform
  auto bg = hilbert_transform(rho, w_mesh);

  // Compare against expected values
  EXPECT_EQ(bg[0][0], 1.0);
}
