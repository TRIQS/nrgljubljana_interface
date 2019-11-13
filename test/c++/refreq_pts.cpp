#include <gtest/gtest.h>

#include <triqs/h5.hpp>
#include <triqs/gfs.hpp>
#include <triqs/test_tools/gfs.hpp>

#include <triqs/gfs/meshes/refreq_pts.hpp>

using namespace triqs::gfs;
using namespace triqs::arrays;

TEST(refreq_pts, Base) {

  // Construction
  auto m = gf_mesh<refreq_pts>{-1.0, 0.0, 2.0};
  EXPECT_ANY_THROW((gf_mesh<refreq_pts>{-1.0, 2.0, 0.0}));
  auto G = gf<refreq_pts, scalar_valued>{m, {}};

  // Mesh Loop Initialization
  for (auto mp : m) G[mp] = double(mp);
  EXPECT_EQ(G.data(), (array<double, 1>{-1.0, 0.0, 2.0}));

  // Placeholder Initialization
  triqs::clef::placeholder<0> om_;
  G[om_] << 2.0 * om_;
  EXPECT_EQ(G.data(), (array<double, 1>{-2.0, 0.0, 4.0}));

  // Manual Initialization
  G[0] = 4.0;
  G[1] = 0.0;
  G[2] = 2.0;
  EXPECT_EQ(G.data(), (array<double, 1>{4.0, 0.0, 2.0}));

  // Test Linear Interpolation
  EXPECT_EQ(G(0.5), 0.5);
  EXPECT_EQ(G(-0.5), 2.0);
}

TEST(refreq_pts, h5) {

  // Construction
  auto m = gf_mesh<refreq_pts>{-1.0, 0.0, 2.0};
  EXPECT_ANY_THROW((gf_mesh<refreq_pts>{-1.0, 2.0, 0.0}));
  auto G = gf<refreq_pts, scalar_valued>{m, {}};

  // Mesh Loop Initialization
  for (auto mp : m) G[mp] = double(mp);
  EXPECT_EQ(G.data(), (array<double, 1>{-1.0, 0.0, 2.0}));

  // Store to file
  {
    auto archive = triqs::h5::file("anderson.out.h5", 'w');
    h5_write(archive, "G", G);
  }

  // Load from file
  {
    auto archive = triqs::h5::file("anderson.out.h5", 'r');
    auto G_h5    = gf<refreq_pts, scalar_valued>{m, {}};
    h5_read(archive, "G", G_h5);
    test_gfs_are_close(G, G_h5);
  }
}
