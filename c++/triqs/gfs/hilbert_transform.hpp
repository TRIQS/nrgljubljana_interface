#include <triqs/gfs.hpp>
#include "./meshes/refreq_pts.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

#include <cmath>

namespace triqs::gfs {

  const double epsabs   = 1e-14;  // numeric integration epsilon (absolute)
  const double epsrel   = 1e-10;  // numeric integration epsilon (relative)
  const double ln1016   = -36.8;  // \approx log(10^-16)
  const size_t limit    = 1000;   // size of integration workspace

  // Unwrap a lambda expression and evaluate it at x.
  // https://martin-ueding.de/articles/cpp-lambda-into-gsl/index.html
  double unwrap(double x, void *p) {
    auto fp = static_cast<std::function<double(double)> *>(p);
    return (*fp)(x);
  }

  // Integrate function f on [a:b].
  double integrate(std::function<double(double)> f, double a, double b, gsl_integration_workspace *work) {
    gsl_function F;
    F.function = &unwrap;
    F.params = &f;
    double result, error;
    int status = gsl_integration_qag(&F, a, b, epsabs, epsrel, limit, GSL_INTEG_GAUSS15, work, &result, &error);
    bool exit_on_error = true;
    if (status && abs(result) > epsabs && exit_on_error) TRIQS_RUNTIME_ERROR << "qag error: " << status << " -- " << gsl_strerror(status);
    return result;
  }

  // Square of x
  double sqr(double x) { return x*x; }

  // Integrate[(-y/(y^2 + (x - omega)^2)), {omega, -B, B}]
  inline double atg(double x, double y, double B) { return atan((-B + x) / y) - atan((B + x) / y); }

  // Integrate[((x - omega)/(y^2 + (x - omega)^2)), {omega, -B, B}]
  inline double logs(double x, double y, double B) { return (-log(sqr(B - x) + sqr(y)) + log(sqr(B + x) + sqr(y))) / 2.0; }

  /**
   * Calculate the Hilbert transform of a given spectral function at fixed complex value z.
   *
   * @param gin The spectral function
   * @param z The complex value for which to evaluate the Hilbert transform
   * @tparam G The Green function type
   */
  template <typename G> typename G::target_t::value_t hilbert_transform(G const &gin, dcomplex z) REQUIRES(is_gf_v<G>) {
    static_assert(std::is_same_v<typename G::target_t, scalar_valued> or std::is_same_v<typename G::target_t, scalar_real_valued>,
                  "Hilbert transform only implemented for scalar valued Green functions");
    static_assert(std::is_same_v<typename G::variable_t, refreq> or std::is_same_v<typename G::variable_t, refreq_pts>,
                  "Hilbert transform only implemented for refreq and refreq_pts meshes");
//    std::cout << std::setprecision(16) << std::endl;
    gsl_set_error_handler_off();
    gsl_integration_workspace *work = gsl_integration_workspace_alloc(limit);
    using DVEC = std::vector<double>;
    DVEC Xpts, Rpts, Ipts;
    for (auto &w : gin.mesh()) {
      double x = w;
      double r = real(gin[w]);
      double i = imag(gin[w]);
      Xpts.push_back(x);
      Rpts.push_back(r);
      Ipts.push_back(i);
//      std::cout << x << " " << r << " " << i << std::endl;
    }
    size_t len = Xpts.size();
    double Xmin = Xpts[0]; // XXX: guaranteed sorted?
    double Xmax = Xpts[len-1];
    assert(Xmin < Xmax);
    double B = std::max(abs(Xmin), abs(Xmax)); // support is [-B:B]
    gsl_interp_accel *accr = gsl_interp_accel_alloc();
    gsl_interp_accel *acci = gsl_interp_accel_alloc();
    gsl_spline *spliner = gsl_spline_alloc(gsl_interp_cspline, len);
    gsl_spline *splinei = gsl_spline_alloc(gsl_interp_cspline, len);
    gsl_spline_init(spliner, &Xpts[0], &Rpts[0], len);
    gsl_spline_init(splinei, &Xpts[0], &Ipts[0], len);
//    std::cout << "Xmin=" << Xmin << " Xmax=" << Xmax << " B=" << B << std::endl;
//    double sum = gsl_spline_eval_integ(spliner, Xmin, Xmax, accr);
//    std::cout << "sum=" << sum << std::endl;
    auto rhor = [Xmin, Xmax, spliner, accr](double x) -> double { return (Xmin <= x && x <= Xmax ? gsl_spline_eval(spliner, x, accr) : 0.0); };
    auto rhoi = [Xmin, Xmax, splinei, acci](double x) -> double { return (Xmin <= x && x <= Xmax ? gsl_spline_eval(splinei, x, acci) : 0.0); };
//    for (auto &w : gin.mesh()) {
//      double x = w;
//      double y = rho(x);
//      std::cout << x << " " << y << std::endl;
//    }
    double x = real(z);
    double y = imag(z);
//    std::cout << x << " " << y << std::endl;
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
//      std::cout << lim1down << " " << lim1up << " " << lim2down << " " << lim2up << " " << W1 << " " << W2 << std::endl;
      return std::make_tuple(lim1down, lim1up, lim2down, lim2up, W1, W2, inside);
    };
    auto calcimA = [x,y,B,rhor,rhoi,limits,work]() -> double {
      // Im part of rho(omega)/(z-omega) with the singularity subtracted out.
      auto imf1 = [x,y,rhor,rhoi](double omega) -> double { return ( (rhor(omega)-rhor(x))*(-y) + (rhoi(omega)-rhoi(x))*(x-omega) )/(sqr(y)+sqr(x-omega)); };
      auto imf2 = [x,y,imf1](double W) -> double { return abs(y)*imf1(abs(y)*W+x); };
      auto imf3p = [x,y,imf2](double r) -> double { return imf2(exp(r)) * exp(r); };
      auto imf3m = [x,y,imf2](double r) -> double { return imf2(-exp(r)) * exp(r); };
      auto [lim1down, lim1up, lim2down, lim2up, W1, W2, inside] = limits();
      auto result1 = (lim1down < lim1up ? integrate(imf3p, lim1down, lim1up, work) : 0.0);
      auto result2 = (lim2down < lim2up ? integrate(imf3m, lim2down, lim2up, work) : 0.0);
      auto result3 = (inside ? rhor(x) * atg(x,y,B) : 0.0);
//      std::cout << "i: " << result1 << " " << result2 << " " << result3 << std::endl;
      return result1 + result2 + result3;
    };
    auto calcimB = [x,y,B,rhor,rhoi,work]() -> double {
      // Im part of rho(omega)/(z-omega)
      auto imf0 = [x,y,rhor,rhoi](double omega) -> double { return (rhor(omega)*(-y) + rhoi(omega)*(x-omega) )/(sqr(y)+sqr(x-omega)); };
      return integrate(imf0, -B, B, work);
    };
    auto calcreA = [x,y,B,rhor,rhoi,limits,work]() -> double {
      // Re part of rho(omega)/(z-omega) with the singularity subtracted out.
      auto ref1 = [x,y,rhor,rhoi](double omega) -> double { return ( (rhor(omega)-rhor(x))*(x-omega) + (rhoi(omega)-rhoi(x))*(y) )/(sqr(y)+sqr(x-omega)); };
      auto ref2 = [x,y,ref1](double W) -> double { return abs(y)*ref1(abs(y)*W+x); };
      auto ref3p = [x,y,ref2](double r) -> double { return ref2(exp(r)) * exp(r); };
      auto ref3m = [x,y,ref2](double r) -> double { return ref2(-exp(r)) * exp(r); };
      auto [lim1down, lim1up, lim2down, lim2up, W1, W2, inside] = limits();
      auto result1 = (lim1down < lim1up ? integrate(ref3p, lim1down, lim1up, work) : 0.0);
      auto result2 = (lim2down < lim2up ? integrate(ref3m, lim2down, lim2up, work) : 0.0);
      auto result3 = (inside ? rhor(x) * logs(x,y,B) : 0.0);
//      std::cout << "r: " << result1 << " " << result2 << " " << result3 << std::endl;
      return result1 + result2 + result3;
    };
    auto calcreB = [x,y,B,rhor,rhoi,work]() -> double {
      // Re part of rho(omega)/(z-omega)
      auto ref0 = [x,y,rhor,rhoi](double omega) -> double { return ( rhor(omega)*(x-omega) + rhoi(omega)*y )/(sqr(y)+sqr(x-omega)); };
      return integrate(ref0, -B, B, work);
    };
    const double LIM_DIRECT = 1e-3; // value y where we switch over to direct integration of rho(E)/(x+Iy-E)
    double re = (abs(y) < LIM_DIRECT ? calcreA() : calcreB());
    double im = (abs(y) < LIM_DIRECT ? calcimA() : calcimB());
    gsl_spline_free(spliner);
    gsl_spline_free(splinei);
    gsl_interp_accel_free(accr);
    gsl_interp_accel_free(acci);
    gsl_integration_workspace_free(work);
//    std::cout << "-> " << re << " " << im << std::endl;
    return dcomplex{re,im};
  }

  /**
   * Calculate the Hilbert transform of a given spectral function on a given mesh
   *
   * @param gin The spectral function
   * @param mesh The mesh
   * @tparam G The Green function type
   * @tparam M The mesh type
   */
  template <typename G, typename M>
  typename G::regular_type hilbert_transform(G const &gin, M const &mesh) REQUIRES(is_gf_v<G> and is_instantiation_of_v<gf_mesh, M>) {
    auto gout = gf<typename M::var_t, typename G::target_t>{mesh, gin.target_shape()};
//    for (auto mp : mesh) gout[mp] = hilbert_transform(gin, dcomplex{mp});
    return gout;
  }

  /**
   * Calculate the Hilbert transform of a block spectral function at fixed complex value z.
   *
   * @param bgin The block spectral function
   * @param z The complex value for which to evaluate the Hilbert transform
   * @tparam BG The block Green function type
   */
  template <typename BG> auto hilbert_transform(BG const &bgin, dcomplex z) REQUIRES(is_block_gf_v<BG>) {
    using G = typename BG::g_t;
    auto l  = [z](G const &gin) { return hilbert_transform<G>(gin, z); };
    return map_block_gf(l, bgin);
  }

  /**
   * Calculate the Hilbert transform of a block spectral function on a given mesh
   *
   * @param gin The spectral function
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
