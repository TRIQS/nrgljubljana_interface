#include <nrgljubljana_interface/solver_core.hpp>

#include <triqs/gfs.hpp>
#include <triqs/h5.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace nrgljubljana_interface;

int main(int argc, char **argv) {

  auto env  = mpi::environment(argc, argv);
  auto comm = mpi::communicator();

  // System Parameters
  double U  = 1.0;
  double eps = -U / 2;

  // Construct Parameters
  constr_params_t cp;
  cp.model = "SIAM";
  cp.symtype = "QS";
  cp.mesh_max = 1.0;
  cp.mesh_min = 1e-3;
  cp.mesh_ratio = 1.1;

  // Set up the Solver
  solver_core S(cp);

  // Solve Parameters
  solve_params_t sp;
  sp.T = 1e-3;
  sp.Lambda = 4.0;
  sp.Nz = 2;
  sp.Tmin = 1e-4;
  sp.keep = 50;
  sp.keepenergy = 6.0;

  // Model parameters
  std::map<std::string, double> mp;
  mp["U1"] = 0.5;
  mp["eps1"] = -0.24;
  sp.model_parameters = mp;

  // Low-level NRG Parameters
  nrg_params_t np;
  np.bins = 50;
  S.set_nrg_params(np);

  // Initialize hybridization function
  for (auto const &w : S.log_mesh)
    S.Delta_w[0][w](0,0) = 0.05i;

  // Solve the impurity model
  S.solve(sp);

  if (comm.rank() == 0) h5_write(file("Solver.h5", 'w'), "S", S);

  // Rerun the Solver and Compare
  auto S_old = solver_core::h5_read_construct(file("Solver.h5", 'r'), "S");
  S_old.solve(S_old.last_solve_params.value());
//  assert_block_gfs_are_close(S_old.G_w, S.G_w, 1e-16);
}
