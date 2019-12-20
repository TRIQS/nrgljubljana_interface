// Rok Zitko, Dec 2019

#include <gtest/gtest.h>

#include <triqs/hilbert_space/fundamental_operator_set.hpp>

#include <triqs/gfs/hilbert_transform.hpp>
#include <triqs/gfs/meshes/refreq_pts.hpp>

#include <triqs/test_tools/gfs.hpp>

#include <cmath>

using namespace triqs::gfs;
using namespace triqs::arrays;
using namespace triqs::clef;

using namespace std::literals::complex_literals;

// Common globals
const auto w_     = placeholder<0>{};
const auto bl_    = placeholder<1>{};
const auto w_mesh = gf_mesh<refreq_pts>{-1.0, -0.5, 0, 0.5, 1.0};

#define EXPECT_CPLX_EQ(a,b) { dcomplex A = a; dcomplex B = b; EXPECT_DOUBLE_EQ(real(A), real(B)); EXPECT_DOUBLE_EQ(imag(A), imag(B)); }

TEST(hilbert_transform, gf_point) {

  // Construction
  auto rho = gf<refreq_pts, scalar_valued>{w_mesh, {}};
  
  if (true) { // flat
    rho[w_] << 1.0;

    // Check all four quadrants
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.5+1.0i), dcomplex(0.4777557225137182,-1.446441332248135));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.5+0.1i), dcomplex(1.081219230625402,-2.877628929964088));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.5+0.01i), dcomplex(1.098434550385856,-3.114928751712775));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.5+0.001i), dcomplex(1.098610510894283,-3.138925989688554));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.5+0.0001i), dcomplex(1.098612270890332,-3.141325986925892));

    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.5-1.0i), dcomplex(0.4777557225137182,1.446441332248135));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.5-0.1i), dcomplex(1.081219230625402,2.877628929964088));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.5-0.01i), dcomplex(1.098434550385856,3.114928751712775));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.5-0.001i), dcomplex(1.098610510894283,3.138925989688554));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.5-0.0001i), dcomplex(1.098612270890332,3.141325986925892));

    EXPECT_CPLX_EQ(hilbert_transform(rho, -0.5+1.0i), dcomplex(-0.4777557225137182,-1.446441332248135));
    EXPECT_CPLX_EQ(hilbert_transform(rho, -0.5+0.1i), dcomplex(-1.081219230625402,-2.877628929964088));
    EXPECT_CPLX_EQ(hilbert_transform(rho, -0.5+0.01i), dcomplex(-1.098434550385856,-3.114928751712775));
    EXPECT_CPLX_EQ(hilbert_transform(rho, -0.5+0.001i), dcomplex(-1.098610510894283,-3.138925989688554));
    EXPECT_CPLX_EQ(hilbert_transform(rho, -0.5+0.0001i), dcomplex(-1.098612270890332,-3.141325986925892));

    EXPECT_CPLX_EQ(hilbert_transform(rho, -0.5-1.0i), dcomplex(-0.4777557225137182,1.446441332248135));
    EXPECT_CPLX_EQ(hilbert_transform(rho, -0.5-0.1i), dcomplex(-1.081219230625402,2.877628929964088));
    EXPECT_CPLX_EQ(hilbert_transform(rho, -0.5-0.01i), dcomplex(-1.098434550385856,3.114928751712775));
    EXPECT_CPLX_EQ(hilbert_transform(rho, -0.5-0.001i), dcomplex(-1.098610510894283,3.138925989688554));
    EXPECT_CPLX_EQ(hilbert_transform(rho, -0.5-0.0001i), dcomplex(-1.098612270890332,3.141325986925892));

    // Decreasing im part
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.0+1.0i), dcomplex(0,-1.570796326794897));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.0+0.1i), dcomplex(0,-2.942255348607469));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.0+0.01i), dcomplex(0,-3.121593320216463));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.0+0.001i), dcomplex(0,-3.13959265425646));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.0+0.0001i), dcomplex(0,-3.14139265359046));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.0+0.00001i), dcomplex(0,-3.141572653589794));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.0+0.000001i), dcomplex(0,-3.141590653589793));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.0+0.0000001i), dcomplex(0,-3.141592453589793));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.0+0.00000001i), dcomplex(0,-3.141592633589793));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.0+0.000000001i), dcomplex(0,-3.141592651589793));
    
    // Increasing im part
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.0+10i), dcomplex(0,-0.1993373049823241));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.0+100i), dcomplex(0,-0.01999933337333048));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.0+1000i), dcomplex(0,-0.001999999333333734));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.0+10000i), dcomplex(0,-0.0001999999993333333));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.0+100000i), dcomplex(0,-1.999999999933334e-05));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.0+1000000i), dcomplex(0,-1.999999999999333e-06));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.0+10000000i), dcomplex(0,-1.999999999999993e-07));
    
    // Close to boundary point
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.99+1.0i), dcomplex(0.8006629534117947,-1.115140355051159));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.99+0.1i), dcomplex(2.987005569216831,-1.620255957002947));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.99+0.01i), dcomplex(4.946743860228897,-2.351169406861533));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.99+0.001i), dcomplex(5.28832978555733,-3.041421488578111));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.99+0.0001i), dcomplex(5.293254828486919,-3.131542735646889));

    EXPECT_CPLX_EQ(hilbert_transform(rho, 1.01+1.0i), dcomplex(0.808662964077975,-1.099141080345021));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 1.01+0.1i), dcomplex(2.996980713653284,-1.421417417534515));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 1.01+0.01i), dcomplex(4.956743693557231,-0.7804230800665934));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 1.01+0.001i), dcomplex(5.298329866391787,-0.09917114009439884));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 1.01+0.0001i), dcomplex(5.303254911796501,-0.009949915442925172));

    // At the boundary point
    EXPECT_CPLX_EQ(hilbert_transform(rho, 1.0+1.0i), dcomplex(0.8047189562170501,-1.107148717794091));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 1.0+0.1i), dcomplex(2.996980713653285,-1.520837931072954));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 1.0+0.01i), dcomplex(5.298329866391789,-1.565796368460938));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 1.0+0.001i), dcomplex(7.600902584542076,-1.570296326836565));
  }
  
  if (true) { // linear
    rho[w_] << 2.0+w_;
  
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.0+1000.0i), dcomplex(-6.666662666669524e-07,-0.003999998666667467));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.0+100.0i), dcomplex(-6.666266695235873e-05,-0.03999866674666095));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.0+10.0i), dcomplex(-0.006626950176759454,-0.3986746099646482));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.0+1.0i), dcomplex(-0.4292036732051034,-3.141592653589793));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.0+0.1i), dcomplex(-1.705774465139253,-5.884510697214939));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.0+0.01i), dcomplex(-1.968784066797835,-6.243186640432926));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.0+0.001i), dcomplex(-1.996860407345744,-6.279185308512919));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.0+0.0001i), dcomplex(-1.999685860734642,-6.282785307180919));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.0+0.00001i), dcomplex(-1.999968584273802,-6.283145307179588));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.0+0.000001i), dcomplex(-1.999996858409345,-6.283181307179587));

    EXPECT_CPLX_EQ(hilbert_transform(rho, 1.5+1000.0i), dcomplex(5.333318733388548e-06,-0.003999991666694316));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 1.5+100.0i), dcomplex(0.0005331873885232592,-0.03999166943052282));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 1.5+10.0i), dcomplex(0.0519262218804473,-0.3919322278279866));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 1.5+1.0i), dcomplex(1.802893696398379,-1.664319233609853));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 1.5+0.1i), dcomplex(3.582935895160238,-0.3918963628891367));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 1.5+0.01i), dcomplex(3.632520806811094,-0.03989828404838237));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 1.5+0.001i), dcomplex(3.633027573530684,-0.003990554748917632));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 1.5+0.0001i), dcomplex(3.63303264231935,-0.000399056201417923));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 1.5+0.00001i), dcomplex(3.63303269300735,-3.990562086832032e-05));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 1.5+0.000001i), dcomplex(3.633032693514235,-3.990562087558569e-06));
  }
  
  if (true) { // boundary points
    rho[w_] << 1.0;
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.999999+0.0001i), dcomplex(9.903437056285835,-1.580745993456891));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 1.0+0.0001i), dcomplex(9.903487553786128,-1.570746326794938));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 1.000001+0.0001i), dcomplex(9.903438056285786,-1.560746660134039));
  }

  if (true) { // quadratic
    rho[w_] << w_*w_;
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.1+0.0001i), dcomplex(-0.1949353971060723,-0.02801896963898113));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.3+0.0001i), dcomplex(-0.5544130613827618,-0.2667433437990276));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.5+0.0001i), dcomplex(-0.7686654163590991,-0.7854361003849807));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.9+0.0001i), dcomplex(0.6316789446424058,-2.604504227768861));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.999+0.0001i), dcomplex(5.638584527557935,-3.035478473543426));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 1.0+0.0001i), dcomplex(9.903487553786128,-1.570746326794938));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 2.0+0.0001i), dcomplex(0.4017492503722853,-2.774977580331286e-05));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 10.0+0.0001i), dcomplex(0.06826901230949829,-6.909867907756972e-07));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 100.0+0.0001i), dcomplex(0.006786123243917645,-6.786941217796011e-09));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 1000.0+0.0001i), dcomplex(0.0006785718375024526,-6.785726553607633e-11));

    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.1+0.01i), dcomplex(-0.1893669676898915,-0.04687736359259637));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.3+0.01i), dcomplex(-0.5361892421811076,-0.2819199341519583));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.5+0.01i), dcomplex(-0.7358192647131685,-0.7888027308223274));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.9+0.01i), dcomplex(0.6781071531276157,-2.483757944444507));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.999+0.01i), dcomplex(3.367850434585618,-1.583988926839967));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 1.0+0.01i), dcomplex(3.379481167807754,-1.487056643611624));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 2.0+0.01i), dcomplex(0.4017278608919441,-0.002774801257511048));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 10.0+0.01i), dcomplex(0.06826894196135165,-6.90986070152268e-05));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 100.0+0.01i), dcomplex(0.006786123175586783,-6.786941149448809e-07));
    EXPECT_CPLX_EQ(hilbert_transform(rho, 1000.0+0.01i), dcomplex(0.0006785718374324345,-6.785726552907443e-09));
  }
  
  if (true) { // imaginary spectral function
    rho[w_] << 1.0i;
    EXPECT_CPLX_EQ(hilbert_transform(rho, 0.5+1.0i), dcomplex(1.446441332248135,0.4777557225137182));
  }
}

TEST(hilbert_transform, gf_mesh) {

  // Construction
  auto rho = gf<refreq_pts, scalar_valued>{w_mesh, {}};

  // Initialize
  rho[w_] << 1.0;

  // Full hilbert transform
  auto g = hilbert_transform(rho, w_mesh);

  // Compare against expected values
  EXPECT_CPLX_EQ(g(-1.0), dcomplex(-37.53450866846468,-1.570796326794897));
  EXPECT_CPLX_EQ(g(-0.5), dcomplex(-1.09861228866811,-3.141592653589793));
  EXPECT_CPLX_EQ(g(0), dcomplex(0,-3.141592653589793));
  EXPECT_CPLX_EQ(g(0.5), dcomplex(1.09861228866811,-3.141592653589793));
  EXPECT_CPLX_EQ(g(1.0), dcomplex(37.53450866846468,-1.570796326794897));
}

TEST(hilbert_transform, gf_gf) {
  auto rho = gf<refreq_pts, scalar_valued>{w_mesh, {}};
  rho[w_] << 1.0;
  
  auto gin = gf<refreq_pts, scalar_valued>{w_mesh, {}};
  gin[w_] << w_ + 1e-16i;
  
  auto g = hilbert_transform(rho, gin);
  
  // Same as above, since epsdefault=1e-16
  EXPECT_CPLX_EQ(g(-1.0), dcomplex(-37.53450866846468,-1.570796326794897));
  EXPECT_CPLX_EQ(g(-0.5), dcomplex(-1.09861228866811,-3.141592653589793));
  EXPECT_CPLX_EQ(g(0), dcomplex(0,-3.141592653589793));
  EXPECT_CPLX_EQ(g(0.5), dcomplex(1.09861228866811,-3.141592653589793));
  EXPECT_CPLX_EQ(g(1.0), dcomplex(37.53450866846468,-1.570796326794897));
}

TEST(hilbert_transform, matrix_gf_point) {
  auto rho = gf<refreq_pts, matrix_valued>{w_mesh, {2,2}};
  rho[w_] << matrix<double>{{1., 0.5}, {3., 2.}};

  auto mat = hilbert_transform_elementwise(rho, 0.5+1e-16i);
  
  EXPECT_CPLX_EQ(mat(0,0), dcomplex(1.09861228866811,-3.141592653589793));
  EXPECT_CPLX_EQ(mat(0,1), 0.5*dcomplex(1.09861228866811,-3.141592653589793));
  EXPECT_CPLX_EQ(mat(1,1), 2.*dcomplex(1.09861228866811,-3.141592653589793));
  EXPECT_CPLX_EQ(mat(1,0), 3.*dcomplex(1.09861228866811,-3.141592653589793));
}

TEST(hilbert_transform, block_gf_point) { // TO DO
  using namespace triqs::gfs;
  using namespace triqs::arrays;
  using namespace triqs::hilbert_space;
  using namespace triqs::utility;
  using namespace triqs::h5;

  // Construction
  //using g_w_t = block_gf<refreq_pts, matrix_valued>;
  //gf_struct_t gf_struct;
  //gf_struct.emplace_back("up", indices_t{1,1});
  //gf_struct.emplace_back("dn", indices_t{1,1});
  //auto rho = g_w_t{w_mesh, gf_struct};
  auto rho = block_gf<refreq_pts, matrix_valued>(w_mesh, {{"up", {1,1}}, {"dn", {1,1}}});

  // Initialize
//  rho[bl_][w_] << 0.;

  // Full hilbert transform
//  auto vals = hilbert_transform(rho, 1.0);

  // Compare against expected values
//  EXPECT_EQ(vals[0], 1.0);
//  EXPECT_EQ(vals[1], 1.0);
}

TEST(hilbert_transform, block_gf_mesh) { // TO DO

  // Construction
//  auto rho = block_gf<refreq_pts, scalar_valued>(w_mesh, {{"up", {}}, {"dn", {}}});
//  auto rho = block_gf<refreq_pts>(w_mesh, {{"up", {}}, {"dn", {}}});
//  auto rho = block_gf<refreq_pts,scalar_valued>(w_mesh, {{"up", {}}, {"dn", {}}});

  // Initialize
//  rho[bl_][w_] << 0.;

  // Full hilbert transform
//  auto bg = hilbert_transform(rho, w_mesh);

  // Compare against expected values
//  EXPECT_EQ(bg[0][0], 1.0);
}
