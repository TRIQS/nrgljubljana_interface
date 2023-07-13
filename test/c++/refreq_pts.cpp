#include <gtest/gtest.h>

#include <h5/h5.hpp>
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <triqs/test_tools/gfs.hpp>

#include <triqs/mesh/refreq_pts.hpp>

using namespace triqs::gfs;
using namespace triqs::mesh;
using namespace nda;

TEST(refreq_pts, Base) {

  // Construction
  auto m = refreq_pts{-1.0, 0.0, 3.0};
  EXPECT_ANY_THROW((refreq_pts{-1.0, 2.0, 0.0}));
  auto G = gf{m};

  // Mesh Loop Initialization
  for (auto mp : m) G[mp] = double(mp);
  EXPECT_EQ(G.data(), (array<double, 1>{-1.0, 0.0, 3.0}));

  // Placeholder Initialization
  nda::clef::placeholder<0> om_;
  G[om_] << 2.0 * om_;
  EXPECT_EQ(G.data(), (array<double, 1>{-2.0, 0.0, 6.0}));

  // Manual Initialization
  G[0] = 4.0;
  G[1] = 0.0;
  G[2] = 3.0;
  EXPECT_EQ(G.data(), (array<double, 1>{4.0, 0.0, 3.0}));

  // Test Linear Interpolation
  EXPECT_EQ(G(0.5), 0.5);
  EXPECT_EQ(G(-0.5), 2.0);

  // Test outside the support
#ifdef NRGLJUBLJANA_INTERFACE_DEBUG
  ASSERT_DEATH(G(4.0), "is_value_valid");
#endif

  // Test on mesh points
  EXPECT_EQ(G(-1.0), 4.0);
  EXPECT_EQ(G(0.0), 0.0);
  EXPECT_EQ(G(3.0), 3.0);
}

TEST(refreq_pts, h5) {

  // Construction
  auto m = refreq_pts{-1.0, 0.0, 2.0};
  EXPECT_ANY_THROW((refreq_pts{-1.0, 2.0, 0.0}));
  auto G = gf{m};

  // Mesh Loop Initialization
  for (auto mp : m) G[mp] = double(mp);
  EXPECT_EQ(G.data(), (array<double, 1>{-1.0, 0.0, 2.0}));

  // Store to file
  {
    auto archive = h5::file("refreq_pts.out.h5", 'w');
    h5_write(archive, "G", G);
  }

  // Load from file
  {
    auto archive = h5::file("refreq_pts.out.h5", 'r');
    auto G_h5    = gf{m};
    h5_read(archive, "G", G_h5);
    test_gfs_are_close(G, G_h5);
  }
}

TEST(refreq_pts, block_gf) {

  // Construction
  auto m   = refreq_pts{-1.0, 0.0, 2.0};
  auto Gbl = block_gf<refreq_pts>{m, {{"bl1", 2}, {"bl2", 2}}};

  // Mesh Loop Initialization
  for (auto mp : m) {
    Gbl[0][mp] = double(mp);
    Gbl[1][mp] = 2 * double(mp);
  }

  auto Gprod = gf{Gbl[0] * Gbl[1]};

  EXPECT_EQ(Gprod.data()(range(), 0, 0), (array<double, 1>{2.0, 0.0, 8.0}));
}

TEST(refreqs_pts, block_gf_scalar) {
  // Construction
//  auto m = gf_mesh<refreq_pts>{-1.0, 0.0, 2.0};
//  auto Gbl = block_gf<refreq_pts>{m, {{"bl1", {}}, {"bl2", {}}}};

//  for (auto mp : m) {
//    Gbl[0][mp] = double(mp);
//    Gbl[1][mp] = 2 * double(mp);
//  }
}
