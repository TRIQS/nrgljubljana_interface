// Rok Zitko, Nils Wentzell, Dec 2019

#include <triqs/gfs.hpp>
#include <itertools/itertools.hpp>
#include "./meshes/refreq_pts.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

#include <cmath>

namespace triqs::gfs {

  // Unwrap a lambda expression and evaluate it at x
  // https://martin-ueding.de/articles/cpp-lambda-into-gsl/index.html
  inline double unwrap(double x, void *p) {
    auto fp = static_cast<std::function<double(double)> *>(p);
    return (*fp)(x);
  }

  class integrator {
  private:
    size_t limit;                    // size of workspace
    bool exit_on_error;              // if true, integration error will trigger a hard error
    gsl_integration_workspace *work; // work space
    gsl_function F;                  // for evaluation of integrand
  public:
    integrator(size_t _limit = 1000, bool _exit_on_error = false) : limit(_limit), exit_on_error(_exit_on_error) {
      work = gsl_integration_workspace_alloc(limit);
      F.function = &unwrap;
    }
    ~integrator() { gsl_integration_workspace_free(work); }
      
    /**
     * Integrate function f on [a:b].
     *
     * @param f Function to be integrated
     * @param a Lower integration range boundary
     * @param b Upper integration range boundary
     * @param epsabs numeric integration epsilon (absolute)
     * @param epsrel numeric integration epsilon (relative)
     */
    double integrate(std::function<double(double)> f, double a, double b, double epsabs = 1e-14, double epsrel = 1e-10) {
      F.params = &f;
      double result, error;
      int status = gsl_integration_qag(&F, a, b, epsabs, epsrel, limit, GSL_INTEG_GAUSS15, work, &result, &error);
      if (status && abs(result) > epsabs && exit_on_error) TRIQS_RUNTIME_ERROR << "qag error: " << status << " -- " << gsl_strerror(status);
      return result;
    }
  };

  // Square of x
  inline double sqr(double x) { return x*x; }

  // Integrate[(-y/(y^2 + (x - omega)^2)), {omega, -B, B}]
  inline double atg(double x, double y, double B) { return atan((-B + x) / y) - atan((B + x) / y); }

  // Integrate[((x - omega)/(y^2 + (x - omega)^2)), {omega, -B, B}]
  inline double logs(double x, double y, double B) { return (-log(sqr(B - x) + sqr(y)) + log(sqr(B + x) + sqr(y))) / 2.0; }

  /**
   * Calculate the Hilbert transform of a given spectral function at fixed complex value z. This is the low-level routine, called from
   * other interfaces. The general strategy is to perform a direct integration of the defining integral for cases where z has a sufficiently
   * large imaginary part, otherwise the singularity is subtracted out and the calculation of the transformed integrand is performed using
   * an integration-variable substitution to better handle small values.
   *
   * @param Ain The spectral function
   * @param z The complex value for which to evaluate the Hilbert transform
   * @param lim_direct value of y=Im(z) above which rho(E)/(x+Iy-E) is directly integrated, and below which the singularity is removed
   * @tparam G The Green function type of Ain
   */
  template <typename G> dcomplex hilbert_transform(G const &Ain, dcomplex z, double lim_direct = 1e-3) REQUIRES(is_gf_v<G>) {
    static_assert(std::is_same_v<typename G::target_t, scalar_valued>,
                  "Hilbert transform only implemented for (complex) scalar-valued spectral functions");
    static_assert(std::is_same_v<typename G::variable_t, refreq> or std::is_same_v<typename G::variable_t, refreq_pts>,
                  "Hilbert transform only implemented for refreq and refreq_pts meshes");
    // Copy the input data for GSL interpolation routines.
    using DVEC = std::vector<double>;
    DVEC Xpts, Rpts, Ipts;
//    auto const & w_mesh = Ain.mesh();
//    auto Xpts = std::vector<double>(w_mesh.size());
//    std::copy(w_mesh.begin(), w_mesh.end(), Xpts.begin());
//    auto Rpts = make_regular(real(Ain.data()));
//    auto Ipts = make_regular(imag(Ain.data()));
    for (const auto &w : Ain.mesh()) {
      double x = w;
      double r = real(Ain[w]);
      double i = imag(Ain[w]);
      Xpts.push_back(x);
      Rpts.push_back(r);
      Ipts.push_back(i);
    }
    size_t len = Xpts.size();
    EXPECTS(std::is_sorted(Xpts.begin(), Xpts.end()));
    double Xmin = Xpts[0];
    double Xmax = Xpts[len-1];
    double B = std::max(abs(Xmin), abs(Xmax)); // half-bandwidth, i.e. support is [-B:B]
    // Initialize GSL and set up the cubic-spline interpolation
    gsl_set_error_handler_off();
    gsl_interp_accel *accr = gsl_interp_accel_alloc();
    gsl_interp_accel *acci = gsl_interp_accel_alloc();
    gsl_spline *spliner = gsl_spline_alloc(gsl_interp_cspline, len);
    gsl_spline *splinei = gsl_spline_alloc(gsl_interp_cspline, len);
    gsl_spline_init(spliner, &Xpts[0], &Rpts[0], len);
    gsl_spline_init(splinei, &Xpts[0], &Ipts[0], len);
    auto rhor = [Xmin, Xmax, spliner, accr](double x) -> double { return (Xmin <= x && x <= Xmax ? gsl_spline_eval(spliner, x, accr) : 0.0); };
    auto rhoi = [Xmin, Xmax, splinei, acci](double x) -> double { return (Xmin <= x && x <= Xmax ? gsl_spline_eval(splinei, x, acci) : 0.0); };
    double x = real(z);
    double y = imag(z);
    // Determine the integration limits for the case with removed singularity.
    auto limits = [x,y,B]() -> auto {
      const double W1 = (x - B) / abs(y); // Rescaled integration limits. Only the absolute value of y matters here.
      const double W2 = (B + x) / abs(y);
      assert(W2 >= W1);
      double lim1down = 1;
      double lim1up   = -1;
      double lim2down = 1;
      double lim2up   = -1;
      bool inside;
      if (W1 < 0 && W2 > 0) {       // x within the band
        const double ln1016 = -36.8; // \approx log(10^-16)
        lim1down = ln1016;
        lim1up   = log(-W1);
        lim2down = ln1016;
        lim2up   = log(W2);
        inside = true;
      } else
      if (W1 > 0 && W2 > 0) {       // x above the band
        lim2down = log(W1);
        lim2up   = log(W2);
        inside = false;
      } else
      if (W1 < 0 && W2 < 0) {       // x below the band
        lim1down = log(-W2);
        lim1up   = log(-W1);
        inside = false;
      } else {                      // special case: boundary points
        inside = true;
      }
      return std::make_tuple(lim1down, lim1up, lim2down, lim2up, inside);
    };
    // Low-level Hilbert-transform routines. calcA routine handles the case with removed singularity and 
    // perform the integration after a change of variables. calcB routine directly evaluates the defining 
    // integral of the Hilbert transform. Real and imaginary parts are determined in separate steps.
    integrator integr;
    auto calcA = [&integr,limits](auto f3p, auto f3m, auto d) -> double {
      auto [lim1down, lim1up, lim2down, lim2up, inside] = limits();
      auto result1 = (lim1down < lim1up ? integr.integrate(f3p, lim1down, lim1up) : 0.0);
      auto result2 = (lim2down < lim2up ? integr.integrate(f3m, lim2down, lim2up) : 0.0);
      auto result3 = (inside ? d : 0.0);
      return result1 + result2 + result3;
    };
    auto calcB = [&integr,B](auto f0) -> double {
      return integr.integrate(f0, -B, B); // direct integration
    };
    // Re part of rho(omega)/(z-omega)
    auto ref0 = [x,y,rhor,rhoi](double omega) -> double { return ( rhor(omega)*(x-omega) + rhoi(omega)*y )/(sqr(y)+sqr(x-omega)); };
    // Im part of rho(omega)/(z-omega)
    auto imf0 = [x,y,rhor,rhoi](double omega) -> double { return (rhor(omega)*(-y) + rhoi(omega)*(x-omega) )/(sqr(y)+sqr(x-omega)); };
    // Re part of rho(omega)/(z-omega) with the singularity subtracted out.
    auto ref1 = [x,y,rhor,rhoi](double omega) -> double { return ( (rhor(omega)-rhor(x))*(x-omega) + (rhoi(omega)-rhoi(x))*(y) )/(sqr(y)+sqr(x-omega)); };
    auto ref2 = [x,y,ref1](double W) -> double { return abs(y)*ref1(abs(y)*W+x); };
    auto ref3p = [x,y,ref2](double r) -> double { return ref2(exp(r)) * exp(r); };
    auto ref3m = [x,y,ref2](double r) -> double { return ref2(-exp(r)) * exp(r); };
    auto red = rhor(x) * logs(x,y,B) - rhoi(x) * atg(x,y,B);
    // Im part of rho(omega)/(z-omega) with the singularity subtracted out.
    auto imf1 = [x,y,rhor,rhoi](double omega) -> double { return ( (rhor(omega)-rhor(x))*(-y) + (rhoi(omega)-rhoi(x))*(x-omega) )/(sqr(y)+sqr(x-omega)); };
    auto imf2 = [x,y,imf1](double W) -> double { return abs(y)*imf1(abs(y)*W+x); };
    auto imf3p = [x,y,imf2](double r) -> double { return imf2(exp(r)) * exp(r); };
    auto imf3m = [x,y,imf2](double r) -> double { return imf2(-exp(r)) * exp(r); };
    auto imd = rhor(x) * atg(x,y,B) + rhoi(x) * logs(x,y,B);
    double re = (abs(y) < lim_direct ? calcA(ref3p, ref3m, red) : calcB(ref0));
    double im = (abs(y) < lim_direct ? calcA(imf3p, imf3m, imd) : calcB(imf0));
    // Deallocate GSL workspaces and return the result
    gsl_spline_free(spliner);
    gsl_spline_free(splinei);
    gsl_interp_accel_free(accr);
    gsl_interp_accel_free(acci);
    return dcomplex{re,im};
  }

  /**
   * Calculate the Hilbert transform of a given spectral function on a given mesh on real axis, i.e., the Green's function.
   *
   * @param Ain The spectral function
   * @param mesh The mesh
   * @param eps The (small) displacement from the real axis. For eps>0, the calculation produces a retarded Green's function.
   * @tparam G The Green function type of Ain
   * @tparam M The mesh type
   */
  template <typename G, typename M>
  typename G::regular_type hilbert_transform(G const &Ain, M const &mesh, double eps = 1e-16) REQUIRES(is_gf_v<G> and is_instantiation_of_v<gf_mesh, M>) {
    auto gout = gf<typename M::var_t, typename G::target_t>{mesh, Ain.target_shape()};
    for (const auto mp : mesh) gout[mp] = hilbert_transform(Ain, dcomplex{mp,eps});
    return gout;
  }

  /**
   * Calculate the Hilbert transform of a given spectral function on complex-valued grid stored as the values of an input Green's function object. The output is
   * an output Green's function on the same frequency domain as the input Green's function. This version of Hilbert transform is the one that is commonly
   * used in the context of DMFT.
   *
   * @param Ain The spectral function
   * @param gin Input Green's function.
   * @tparam G The Green function type of Ain and gin.
   */
  template <typename G>
  typename G::regular_type hilbert_transform(G const &Ain, G const &gin) REQUIRES(is_gf_v<G>) {
    auto gout = gin;
    for (const auto &mp : gin.mesh()) gout[mp] = hilbert_transform(Ain, gin[mp]);
    return gout;
  }

  /**
   * Calculate the Hilbert transform of a given spectral function at fixed complex value z. This version performs an element-wise
   * Hilbert transformation.
   *
   * @param Ain The matrix-valued spectral function
   * @param z The complex value for which to evaluate the Hilbert transform
   * @tparam G The Green function type of Ain
   */
  template <typename G> matrix<dcomplex> hilbert_transform_elementwise(G const &Ain, dcomplex z) REQUIRES(is_gf_v<G>) {
    static_assert(std::is_same_v<typename G::target_t, matrix_valued>,
                  "Hilbert transform only implemented for matrix-valued spectral functions");
    long size = Ain.target_shape()[0];
    auto mat = matrix<dcomplex>(size, size);
    for (auto [i, j] : itertools::product_range(size, size)) {
      auto gtemp = gf<refreq_pts, scalar_valued>{Ain.mesh(), {}};       // XXX: refreqs_pts -> typename G::mesh_t
      for (const auto &mp : Ain.mesh()) gtemp[mp] = Ain[mp](i,j);
      mat(i,j) = hilbert_transform(gtemp, z);
    }
    return mat;
  }

  /**
   * Calculate the Hilbert transform of a block spectral function at fixed complex value z.
   *
   * @param bAin The block spectral function
   * @param z The complex value for which to evaluate the Hilbert transform
   * @tparam BG The block Green function type
   */
  template <typename BG> auto hilbert_transform(BG const &bAin, dcomplex z) REQUIRES(is_block_gf_v<BG>) {
    using G = typename BG::g_t;
    auto l  = [z](G const &Ain) { return hilbert_transform<G>(Ain, z); };
    return map_block_gf(l, bAin);
  }

  /**
   * Calculate the Hilbert transform of a block spectral function on a given mesh
   *
   * @param Ain The spectral function
   * @param mesh The mesh
   * @tparam G The block Green function type
   * @tparam M The mesh type
   */
  template <typename BG, typename M>
  auto hilbert_transform(BG const &bAin, M const &mesh, double eps = 1e-16) REQUIRES(is_block_gf_v<BG> and is_instantiation_of_v<gf_mesh, M>) {
    using G = typename BG::g_t;
    auto l  = [&mesh, eps](G const &Ain) { return hilbert_transform<G, M>(Ain, mesh, eps); };
    return map_block_gf(l, bAin);
  }

} // namespace triqs::gfs
