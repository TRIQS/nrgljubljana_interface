/*******************************************************************************
 *
 * nrgljubljana_interface: A TRIQS interface to the nrgljubliana impurity solver
 *
 * Copyright (c) 2019 The Simons foundation
 *   authors: Rok Zitko, Nils Wentzell
 *
 * nrgljubljana_interface is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * nrgljubljana_interface is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * nrgljubljana_interface. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#include "./solver_core.hpp"

//#include <boost/filesystem.hpp> // or C++17 filesystem ?

#include <triqs/utility/exceptions.hpp> // TRIQS_RUNTIME_ERROR

#include <algorithm>  // max
#include <cmath>      // pow, log
#include <cstdlib>    // system
#include <sys/stat.h> // mkdir
#include <cassert>    // assert
#include <iostream>
#include <fstream>
#include <utility>

#include <nrg-lib.h>

namespace nrgljubljana_interface {

  solver_core::solver_core(constr_params_t cp) : constr_params(cp) {

    // Initialize the hybridization function container TODO Generate Log Mesh
    Delta_w = g_w_t{{-1, -1e-99, 1e-99, 1}, {1, 1}};
  }

  // -------------------------------------------------------------------------------

  void solver_core::generate_param_file(double z) {
    const constr_params_t &cp = constr_params;
    const solve_params_t &sp  = solve_params;
    nrg_params_t &np          = nrg_params;
    np.z                      = z;
    // Automatically establish appropriate default values for the high-level interface
    set_params();
    // Generate the parameter file
    std::ofstream F("param");
    F << std::boolalpha; // important: we want to output true/false strings
    F << "[param]" << std::endl;
    F << "bandrescale=" << np.bandrescale << std::endl;
    F << "model=" << cp.problem << std::endl; // not required when using templates
    F << "symtype=QS" << std::endl;           // from template database TODO
    F << "Lambda=" << sp.Lambda << std::endl;
    F << "xmax=" << np.xmax << std::endl;
    F << "Nmax=" << np.Nmax << std::endl;
    F << "mMAX=" << np.mMAX << std::endl;
    //    F << "Tmin=" << sp.Tmin << std::endl;
    F << "keep=" << sp.keep << std::endl;
    F << "keepenergy=" << sp.keepenergy << std::endl;
    F << "keepmin=" << sp.keepmin << std::endl;
    F << "T=" << sp.T << std::endl;
    F << "ops=" << sp.ops << std::endl;
    F << "specs=" << sp.specs << std::endl;
    F << "specd=" << sp.specd << std::endl;
    F << "spect=" << sp.spect << std::endl;
    F << "specq=" << sp.specq << std::endl;
    F << "specot=" << sp.specot << std::endl;
    F << "specgt=" << sp.specgt << std::endl;
    F << "speci1t=" << sp.speci1t << std::endl;
    F << "speci2t=" << sp.speci2t << std::endl;
    F << "specchit=" << sp.specchit << std::endl;
    F << "specv3=" << sp.specv3 << std::endl;
    F << "v3mm=" << sp.v3mm << std::endl;
    F << "dmnrg=" << sp.dmnrg << std::endl;
    F << "cfs=" << sp.cfs << std::endl;
    F << "fdm=" << sp.fdm << std::endl;
    F << "fdmexpv=" << sp.fdmexpv << std::endl;
    F << "dmnrgmats=" << sp.dmnrgmats << std::endl;
    F << "fdmmats=" << sp.fdmmats << std::endl;
    F << "mats=" << sp.mats << std::endl;
    F << "alpha=" << sp.alpha << std::endl;
    F << "gamma=" << sp.gamma << std::endl; // ?
    F << "discretization=" << np.discretization << std::endl;
    F << "z=" << np.z << std::endl;
    F << "Nz=" << sp.Nz << std::endl; // !
    F << "polarized=" << np.polarized << std::endl;
    F << "pol2x2=" << np.pol2x2 << std::endl;
    F << "rungs=" << np.rungs << std::endl;
    F << "tri=" << np.tri << std::endl;
    F << "preccpp=" << np.preccpp << std::endl;
    F << "diag=" << np.diag << std::endl;
    F << "diagratio=" << np.diagratio << std::endl;
    F << "dsyevrlimit=" << np.dsyevrlimit << std::endl;
    F << "zheevrlimit=" << np.zheevrlimit << std::endl;
    F << "restart=" << np.restart << std::endl;
    F << "restartfactor=" << np.restartfactor << std::endl;
    F << "safeguard=" << np.safeguard << std::endl;
    F << "safeguardmax=" << np.safeguardmax << std::endl;
    F << "fixeps=" << np.fixeps << std::endl;
    F << "betabar=" << np.betabar << std::endl;
    F << "gtp=" << np.gtp << std::endl;
    F << "chitp=" << np.chitp << std::endl;
    F << "finite=" << np.finite << std::endl;
    F << "cfsgt=" << np.cfsgt << std::endl;
    F << "cfsls=" << np.cfsls << std::endl;
    F << "fdmgt=" << np.fdmgt << std::endl;
    F << "fdmls=" << np.fdmls << std::endl;
    F << "fdmexpvn=" << np.fdmexpvn << std::endl;
    F << "finitemats=" << np.finitemats << std::endl;
    F << "dm=" << np.dm << std::endl;
    F << "broaden_max=" << cp.mesh_max << std::endl;                // !
    F << "broaden_min=" << cp.mesh_min << std::endl;                // !
    F << "broaden_ratio=" << cp.mesh_ratio << std::endl;            // !
    F << "broaden_min_ratio=" << np.broaden_min_ratio << std::endl; // keep this one?
    F << "omega0=" << np.omega0 << std::endl;
    F << "omega0_ratio=" << np.omega0_ratio << std::endl;
    F << "diagth=" << np.diagth << std::endl;
    F << "substeps=" << np.substeps << std::endl;
    F << "strategy=" << np.strategy << std::endl;
    F << "Ninit=" << np.Ninit << std::endl;
    F << "reim=" << np.reim << std::endl;
    F << "dumpannotated=" << np.dumpannotated << std::endl;
    F << "dumpabs=" << np.dumpabs << std::endl;
    F << "dumpscaled=" << np.dumpscaled << std::endl;
    F << "dumpprecision=" << np.dumpprecision << std::endl;
    F << "dumpgroups=" << np.dumpgroups << std::endl;
    F << "grouptol=" << np.grouptol << std::endl;
    F << "dumpdiagonal=" << np.dumpdiagonal << std::endl;
    F << "savebins=" << np.savebins << std::endl;
    F << "broaden=" << np.broaden << std::endl;
    F << "emin=" << np.emin << std::endl;
    F << "emax=" << np.emax << std::endl;
    F << "bins=" << np.bins << std::endl;
    F << "accumulation=" << np.accumulation << std::endl;
    F << "linstep=" << np.linstep << std::endl;
    F << "discard_trim=" << np.discard_trim << std::endl;
    F << "discard_immediately=" << np.discard_immediately << std::endl;
    F << "goodE=" << np.goodE << std::endl;
    F << "NN1=" << np.NN1 << std::endl;
    F << "NN2even=" << np.NN2even << std::endl;
    F << "NN2avg=" << np.NN2avg << std::endl;
    F << "NNtanh=" << np.NNtanh << std::endl;
    F << "width_td=" << np.width_td << std::endl;
    F << "width_custom=" << np.width_custom << std::endl;
    F << "prec_td=" << np.prec_td << std::endl;
    F << "prec_custom=" << np.prec_custom << std::endl;
    F << "prec_xy=" << np.prec_xy << std::endl;
    F << "resume=" << np.resume << std::endl;
    F << "log=" << np.log << std::endl;
    F << "logall=" << np.logall << std::endl;
    F << "done=" << np.done << std::endl;
    F << "calc0=" << np.calc0 << std::endl;
    F << "lastall=" << np.lastall << std::endl;
    F << "lastalloverride=" << np.lastalloverride << std::endl;
    F << "dumpsubspaces=" << np.dumpsubspaces << std::endl;
    F << "dump_f=" << np.dump_f << std::endl;
    F << "dumpenergies=" << np.dumpenergies << std::endl;
    F << "logenumber=" << np.logenumber << std::endl;
    F << "stopafter=" << np.stopafter << std::endl;
    F << "forcestop=" << np.forcestop << std::endl;
    F << "removefiles=" << np.removefiles << std::endl;
    F << "noimag=" << np.noimag << std::endl;
    F << "checksumrules=" << np.checksumrules << std::endl;
    F << "checkdiag=" << np.checkdiag << std::endl;
    F << "checkrho=" << np.checkrho << std::endl;
    F << "dos=Delta.dat" << std::endl; // hard-coded
    F << "[extra]" << std::endl;
    for (const auto &i : sp.model_parameters) F << i.first << "=" << i.second << std::endl;
  }

  void solver_core::set_params() {
    const constr_params_t &cp = constr_params;
    const solve_params_t &sp  = solve_params;
    nrg_params_t &np          = nrg_params; // only nrg_params allowed to be changed!

    // Test if the low-level paramerers are sensible for use with the
    // high-level interface.
    if (np.discretization != "Z") TRIQS_RUNTIME_ERROR << "Must use discretization=Z in the high-level solver interface.";

    // Automatically set (override) some low-level parameters
    if (sp.Tmin > 0) { // If Tmin is set, determine the required length of the Wilson chain.
      np.Nmax    = 0;
      auto scale = [=](int n) {
        return (1. - 1. / sp.Lambda) / std::log(sp.Lambda) * std::pow(sp.Lambda, -(np.z - 1)) * std::pow(sp.Lambda, -(n - 1) / 2.);
      };
      while (scale(np.Nmax + 1) >= sp.Tmin) np.Nmax++;
    }
    if (np.mMAX < 0) // Determine the number of sites in the star representation
      np.mMAX = std::max(80, 2 * np.Nmax);
    if (np.xmax < 0) // Length of the x-interval in the discretization=Z (ODE) approach
      np.xmax = np.Nmax / 2. + 2.;
    if (np.bandrescale < 0) // Make the NRG energy window correspond to the extend of the frequency mesh
      np.bandrescale = cp.mesh_max;
  }

  void solver_core::set_nrg_params(nrg_params_t const &nrg_params_) { nrg_params = nrg_params_; }

  // Solve the problem for a given value of the twist parameter z
  void solver_core::solve_one_z(double z, std::string taskdir) {
    if (world.rank() == 0) {
      generate_param_file(z);
      if (mkdir(taskdir.c_str(), 0755) != 0) TRIQS_RUNTIME_ERROR << "failed to mkdir taskdir " << taskdir;
      std::string cmd = "./instantiate " + taskdir;
      if (system(cmd.c_str()) != 0) TRIQS_RUNTIME_ERROR << "Running " << cmd << " failed";
      if (chdir(taskdir.c_str()) != 0) TRIQS_RUNTIME_ERROR << "failed to chdir to taskdir " << taskdir;
      // Solve the impurity model
      set_workdir("."); // may be overridden by NRG_WORKDIR in environment
      run_nrg_master();
      if (chdir("..") != 0) TRIQS_RUNTIME_ERROR << "failed to return from taskdir " << taskdir;
    } else {
      //	 run_nrg_slave();
    }
  }

  std::string solver_core::create_tempdir() {
    std::string tempdir_template = "nrg_tempdir_XXXXXX";
    char x[tempdir_template.length() + 1];
    strncpy(x, tempdir_template.c_str(), tempdir_template.length() + 1);
    if (char *w = mkdtemp(x)) // create a unique directory
      return w;
    else
      TRIQS_RUNTIME_ERROR << "Failed to create a directory for temporary files.";
  }

  void solver_core::generate_hyb_file() {
    auto m = gf_mesh<refreq_pts>{-1, -1e-99, 1e-99, 1};
    auto g = gf<refreq_pts>{m, {1, 1}};
    g[0]   = 0.1;
    g[1]   = 0.11;
    g[2]   = 0.3;
    g[3]   = 0.2;
    //save_to_file(g, "Delta.dat");
    std::ofstream F("Delta.dat");
    for (auto w : g.mesh()) F << double(w) << " " << g[w](0, 0).real() << std::endl;
  }

  void solver_core::solve(solve_params_t const &sp) {
    solve_params = sp;
    std::string tempdir;
    // Reset the results
    container_set::operator=(container_set{});
    if (world.rank() == 0) {
      std::cout << "\nNRGLJUBLJANA_INTERFACE Solver\n";
      // Create a temporary directory
      tempdir = create_tempdir();
      if (chdir(tempdir.c_str()) != 0) TRIQS_RUNTIME_ERROR << "chdir to tempdir failed.";
      // Generate the hybdridisation function
      generate_hyb_file();
      generate_param_file(1.0); // we need a mock param file for 'adapt' tool
      // Copy files from the template library & discretize
      std::string templatedir = constr_params.templatedir;
      if (templatedir.empty())
        if (const char *env_tdir = std::getenv("NRGIF_TEMPLATE_DIR")) {
          templatedir = env_tdir;
        } else {
          templatedir = NRGIF_TEMPLATE_DIR;
        }
      std::string script = templatedir + "/" + constr_params.problem + "/prepare";
      if (system(script.c_str()) != 0) TRIQS_RUNTIME_ERROR << "Running prepare script failed: " << script;
      if (system("./discretize") != 0) TRIQS_RUNTIME_ERROR << "Running discretize script failed";
    }
    // Perform the calculations (this must run in all MPI processes)
    const double dz  = 1.0 / sp.Nz;
    const double eps = 1e-8;
    int cnt          = 1;
    for (double z = dz; z <= 1.0 + eps; z += dz, cnt++) {
      std::string taskdir = std::to_string(cnt);
      solve_one_z(z, taskdir);
    }
    assert(cnt == sp.Nz + 1);
    if (world.rank() == 0) {
      if (system("./process") != 0) TRIQS_RUNTIME_ERROR << "Running post-processing script failed";
      // Cleanup
      if (chdir("..") != 0) TRIQS_RUNTIME_ERROR << "failed to return from the tempdir";
      remove(tempdir.c_str());
    }
  }

  //  void solver_core::run_single(all_solve_params_t const &all_solve_params) {
  //  }

  // -------------------------------------------------------------------------------

  // Function that writes a solver object to hdf5 file

  void h5_write(triqs::h5::group h5group, std::string subgroup_name, solver_core const &s) {
    auto grp = h5group.create_group(subgroup_name);
    h5_write_attribute(grp, "TRIQS_HDF5_data_scheme", solver_core::hdf5_scheme());
    h5_write_attribute(grp, "TRIQS_GIT_HASH", std::string(AS_STRING(TRIQS_GIT_HASH)));
    h5_write_attribute(grp, "NRGLJUBLJANA_INTERFACE_GIT_HASH", std::string(AS_STRING(NRGLJUBLJANA_INTERFACE_GIT_HASH)));
    h5_write(grp, "", s.result_set());
    h5_write(grp, "constr_params", s.constr_params);
    h5_write(grp, "solve_params", s.solve_params);
    h5_write(grp, "nrg_params", s.nrg_params);
  }

  // Function that constructs a solver object from an hdf5 file
  solver_core solver_core::h5_read_construct(triqs::h5::group h5group, std::string subgroup_name) {
    auto grp           = h5group.open_group(subgroup_name);
    auto constr_params = h5_read<constr_params_t>(grp, "constr_params");
    auto s             = solver_core{constr_params};
    h5_read(grp, "", s.result_set());
    h5_read(grp, "solve_params", s.solve_params);
    h5_read(grp, "nrg_params", s.nrg_params);
    return s;
  }
} // namespace nrgljubljana_interface
