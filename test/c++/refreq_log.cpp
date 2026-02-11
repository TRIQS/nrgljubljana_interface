#include <gtest/gtest.h>

#include <h5/h5.hpp>
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace triqs::gfs;
using namespace triqs::mesh;
using namespace nda;

TEST(refreq_log, Base) {

  // Construction: eps=0.1, w_max=10, ratio=2
  auto m = refreq_log{0.1, 10.0, 2.0};
  auto G = gf{m};

  // Mesh Loop Initialization
  for (auto mp : m) G[mp] = double(mp);

  // Placeholder Initialization
  nda::clef::placeholder<0> om_;
  G[om_] << 2.0 * om_;

  // Manual Initialization
  G[0] = 4.0;
  G[1] = 0.0;
}

TEST(refreq_log, h5) {

  // Construction
  auto m = refreq_log{0.1, 10.0, 2.0};
  auto G = gf{m};

  // Mesh Loop Initialization
  for (auto mp : m) G[mp] = double(mp);

  // Store to file
  {
    auto archive = h5::file("refreq_log.out.h5", 'w');
    h5_write(archive, "G", G);
  }

  // Load from file
  {
    auto archive = h5::file("refreq_log.out.h5", 'r');
    auto G_h5    = gf{m};
    h5_read(archive, "G", G_h5);
    test_gfs_are_close(G, G_h5);
  }
}

TEST(refreq_log, block_gf) {

  // Construction
  auto m   = refreq_log{0.1, 10.0, 2.0};
  auto Gbl = block_gf<refreq_log>{m, {{"bl1", 2}, {"bl2", 2}}};

  // Mesh Loop Initialization
  for (auto mp : m) {
    Gbl[0][mp] = double(mp);
    Gbl[1][mp] = 2 * double(mp);
  }

  auto Gprod = gf{Gbl[0] * Gbl[1]};
}
