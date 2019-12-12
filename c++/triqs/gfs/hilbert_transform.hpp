#include <triqs/gfs.hpp>
#include "./meshes/refreq_pts.hpp"

namespace triqs::gfs {

  /**
   * Calculate the hilbert transform of a given spectral function at fixed complex value z.
   *
   * @param gin The sectral function
   * @param z The complex value for which to evaluate the Hilbert transform
   * @tparam G The Green function type
   */
  template <typename G> typename G::target_t::value_t hilbert_transform(G const &gin, dcomplex z) REQUIRES(is_gf_v<G>) {
    static_assert(std::is_same_v<typename G::target_t, scalar_valued> or std::is_same_v<typename G::target_t, scalar_real_valued>,
                  "Hilbert transform only implemented for scalar valued Green functions");
    static_assert(std::is_same_v<typename G::variable_t, refreq> or std::is_same_v<typename G::variable_t, refreq_pts>,
                  "Hilbert transform only implemented for refreq and refreq_pts meshes");
    // === TODO ===
    // Calculate hilbert transform of gin at value z and return
    return {};
  }

  /**
   * Calculate the hilbert transform of a given spectral function on a given mesh
   *
   * @param gin The sectral function
   * @param mesh The mesh
   * @tparam G The Green function type
   * @tparam M The mesh type
   */
  template <typename G, typename M>
  typename G::regular_type hilbert_transform(G const &gin, M const &mesh) REQUIRES(is_gf_v<G> and is_instantiation_of_v<gf_mesh, M>) {
    auto gout = gf<typename M::var_t, typename G::target_t>{mesh, gin.target_shape()};
    for (auto mp : mesh) gout[mp] = hilbert_transform(gin, dcomplex{mp});
    return gout;
  }

  /**
   * Calculate the hilbert transform of a block spectral function at fixed complex value z.
   *
   * @param bgin The block sectral function
   * @param z The complex value for which to evaluate the Hilbert transform
   * @tparam BG The block Green function type
   */
  template <typename BG> auto hilbert_transform(BG const &bgin, dcomplex z) REQUIRES(is_block_gf_v<BG>) {
    using G = typename BG::g_t;
    auto l  = [z](G const &gin) { return hilbert_transform<G>(gin, z); };
    return map_block_gf(l, bgin);
  }

  /**
   * Calculate the hilbert transform of a block spectral function on a given mesh
   *
   * @param gin The sectral function
   * @param mesh The mesh
   * @tparam G The block Green function type
   * @tparam M The mesh type
   */
  template <typename BG, typename M>
  auto hilbert_transform(BG const &bgin, M const &mesh) REQUIRES(is_block_gf_v<BG> and is_instantiation_of_v<gf_mesh, M>) {
    using G = typename BG::g_t;
    auto l  = [&mesh](G const &gin) { return hilbert_transform<G, M>(gin, mesh); };
    return map_block_gf(l, bgin);
  }

} // namespace triqs::gfs
