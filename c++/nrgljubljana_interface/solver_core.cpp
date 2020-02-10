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

#include <triqs/utility/exceptions.hpp> // TRIQS_RUNTIME_ERROR
#include <itertools/itertools.hpp>

#include <algorithm>  // max
#include <cmath>      // pow, log
#include <cstdlib>    // system
#include <sys/stat.h> // mkdir
#include <iostream>
#include <fstream>
#include <utility>
#include <deque>

#include <boost/lexical_cast.hpp>
#include <nrg-lib.h>

#include <mpi/mpi.hpp>
#include <mpi/string.hpp>

// std::filesystem is not production-ready as of 2020. As a workaround,
// we remove temporary files using "rm".
//
//#ifdef __has_include
//#  if __has_include(<filesystem>)
//#    include <filesystem>
//namespace fs = std::filesystem;
//#  elif __has_include(<experimental/filesystem>)
//#    include <experimental/filesystem>
//namespace fs = std::experimental::filesystem;
//#  else
//#    error "Missing <filesystem>"
//#  endif
//#endif

namespace nrgljubljana_interface {

  solver_core::solver_core(constr_params_t cp) : constr_params(cp) {
    gf_struct = read_structure("gf_struct", true); // true=mandatory
    chi_struct = read_structure("chi_struct", false); // false=optional

    // Read model-specific template info file
    std::ifstream F(constr_params.get_model_dir() + "/info");
    if (!F) TRIQS_RUNTIME_ERROR << "Failed to open template info file";
    auto getline = [&F](const std::string &keyword) -> std::string {
      std::string s;
      std::getline(F, s);
      auto pos = s.find(':');
      std::string keyword_found = s.substr(0, pos);
      TRIQS_ASSERT2(keyword == keyword_found, "template info file is corrupted");
      std::string value = s.substr(pos+1);
      return value;
    };
    constr_params.ops       = getline("ops");
    constr_params.specs     = getline("specs");
    constr_params.specd     = getline("specd");
    constr_params.spect     = getline("spect");
    constr_params.specq     = getline("specq");
    constr_params.specot    = getline("specot");
    constr_params.params    = getline("params");
    constr_params.polarized = getline("polarized") == "true";

    // Create the hybridization function on a logarithmic mesh
    std::vector<double> mesh_points;
    for (double w = cp.mesh_max; w > cp.mesh_min; w /= cp.mesh_ratio) {
      mesh_points.push_back(w);
      mesh_points.push_back(-w);
    }
    std::sort(begin(mesh_points), end(mesh_points));
    log_mesh = gf_mesh<refreq_pts>{mesh_points};
    Delta_w  = g_w_t{log_mesh, gf_struct};
    // We also construct G_w and Sigma_w here to enable their initialization in DMFT loops
    // prior to solver calls.
    G_w = g_w_t{log_mesh, gf_struct};
    Sigma_w = g_w_t{log_mesh, gf_struct};
  }

  // -------------------------------------------------------------------------------

  gf_struct_t solver_core::read_structure(const std::string &filename, bool mandatory = true) {
    std::ifstream F(constr_params.get_model_dir() + "/" + filename);
    if (mandatory && not F.is_open()) TRIQS_RUNTIME_ERROR << "Failed to open structure file " << filename;
    gf_struct_t _gf_struct;
    if (F) {
      std::string bl_name;
      int bl_size;
      while (F >> bl_name >> bl_size) {
        indices_t idx_lst;
        for (int i : range(bl_size)) idx_lst.push_back(i);
        _gf_struct.emplace_back(bl_name, idx_lst);
      }
    }
    return _gf_struct;
  }

  void solver_core::readGF(const std::string &name, std::optional<g_w_t> &G_w, gf_struct_t &_gf_struct) {
    G_w = g_w_t{log_mesh, _gf_struct};
    for (int bl_idx : range(_gf_struct.size())) {
      long bl_size = Delta_w[bl_idx].target_shape()[0];
      for (auto [i, j] : product_range(bl_size, bl_size)) {
        auto bl_name     = Delta_w.block_names()[bl_idx];
        auto file_ending = std::string{"_"} + bl_name + "_" + std::to_string(i) + std::to_string(j) + ".dat";
        std::string imGfilename = "im" + name + file_ending;
        std::ifstream imG(imGfilename);
        std::string reGfilename = "re" + name + file_ending;
        std::ifstream reG(reGfilename);
        if (imG && reG) {
          double w, re, im;
          for (auto const &mp : log_mesh) {
            imG >> w >> im;
            reG >> w >> re;
            TRIQS_ASSERT2(abs(w-double(mp)) < 1e-8*abs(w), "frequency mismatch");
            (*G_w)[bl_idx][mp](i, j) = re + 1i * im;
          }
        } else {
          // If the files cannot be read (do not exist), the corresponding GF matrix component
          // is zero-ed out.
          for (auto const &mp : log_mesh) (*G_w)[bl_idx][mp](i, j) = 0.0;
        }
      }
    }
  }
  
  void solver_core::readA(const std::string &name, std::optional<g_w_t> &A_w, gf_struct_t &_gf_struct) {
    A_w = g_w_t{log_mesh, _gf_struct};
    for (int bl_idx : range(_gf_struct.size())) {
      long bl_size = Delta_w[bl_idx].target_shape()[0];
      for (auto [i, j] : product_range(bl_size, bl_size)) {
        auto bl_name     = Delta_w.block_names()[bl_idx];
        auto file_ending = std::string{"_"} + bl_name + "_" + std::to_string(i) + std::to_string(j) + ".dat";
        std::string Afilename = name + file_ending;
        std::ifstream A(Afilename);
        if (A) {
          double w, re;
          for (auto const &mp : log_mesh) {
            A >> w >> re;
            TRIQS_ASSERT2(abs(w-double(mp)) < 1e-8*abs(w), "frequency mismatch");
            (*A_w)[bl_idx][mp](i, j) = re;
          }
        } else {
          for (auto const &mp : log_mesh) (*A_w)[bl_idx][mp](i, j) = 0.0;
        }
      }
    }
  }

  // Read expectation values and average over Nz runs
  void solver_core::readexpv(int Nz) {
    for (int cnt = 1; cnt <= Nz; cnt++) {
      const std::string expvfilename = std::to_string(cnt) + "/customfdm";
      std::ifstream F(expvfilename);
      if (!F) TRIQS_RUNTIME_ERROR << "Expectation values output file not found.";
      std::string snumber, skeyword, svalue;
      getline(F, snumber); // snumber not used
      getline(F, skeyword);
      getline(F, svalue);
      F.close();
      skeyword.erase(0,1); // drop #
      std::istringstream sskeyword(skeyword);
      std::deque<std::string> keywords{std::istream_iterator<std::string>{sskeyword},
        std::istream_iterator<std::string>{}};
      std::istringstream ssvalue(svalue);
      std::deque<double> values{std::istream_iterator<double>{ssvalue},
        std::istream_iterator<double>{}};
      TRIQS_ASSERT2(keywords.size() == values.size(), "customfdm corrupted");
      keywords.pop_front(); // ignore first column (temperature)
      values.pop_front();
      for (auto [k, v] : zip(keywords, values)) expv[k] += v; // zeroed by default constructor
    }
    for (auto &i : expv) { i.second /= Nz; } // calculate the average over discretization meshes
  }

  inline void call(std::string command, bool verbose = true, bool exit_on_failure = true) {
    std::string s = command + (verbose ? "" : " >/dev/null");
    std::cout << "running " << s << std::endl;
    if (system(s.c_str()) != 0) {
      if (exit_on_failure) 
        TRIQS_RUNTIME_ERROR << "Failure running NRGLjubljana_interface script: " << s;
      else
        std::cerr << "Warning: failure running NRGLjubljana_interface script: " << s << std::endl;
    }
  }
   
  // Write Gamma=-Im(Delta_w) to a file
  void solver_core::write_gamma() {
    for (int bl_idx : range(gf_struct.size())) {
      long bl_size = Delta_w[bl_idx].target_shape()[0];
      auto bl_name = Delta_w.block_names()[bl_idx];
      for (auto [i, j] : product_range(bl_size, bl_size)) {
        std::ofstream F("Gamma_" + bl_name + "_" + std::to_string(i) + std::to_string(j) + ".dat");
        const double cutoff = 1e-8; // ensure hybridisation function is positive
        for (auto const &w : Delta_w[bl_idx].mesh()) {
          double value = -Delta_w[bl_idx][w](i, j).imag();
          if (value < cutoff) { value = cutoff; }
          F << double(w) << " " << value << std::endl;
        }
      }
    }
  }
   
  void solver_core::solve(solve_params_t const &sp) {
    last_solve_params = sp;
    container_set::operator=(container_set{});  // Reset the results
    std::string tempdir{};
    if (world.rank() == 0) {
      std::cout << "NRGLJUBLJANA_INTERFACE Solver\n";
      // Create a temporary directory
      tempdir = create_tempdir();
    }
    world.barrier(); // Ensures temporary directory is created
    mpi::broadcast(tempdir, world);
    if (chdir(tempdir.c_str()) != 0) TRIQS_RUNTIME_ERROR << "chdir to tempdir failed.";
    if (world.rank() == 0) {
      write_gamma();
      // we need a mock param file for 'adapt' tool
      generate_param_file(1.0);
      call(constr_params.get_model_dir() + "/prepare", verbose);
      call("./discretize", verbose);
      for (int i: range(1,sp.Nz+1))
        instantiate(double(i)/sp.Nz, std::to_string(i));
    }
    world.barrier(); // Ensures the initialization scripts have completed

    // Perform the calculations in parallel
    auto chunk = mpi::chunk(range(1,sp.Nz+1), world);
    for (int i: chunk)
      solve_one(std::to_string(i));

    world.barrier(); // Ensures all subcalculations are completed
    // Post-Processing in perl script
    if (world.rank() == 0)
      call("./process", verbose);

    world.barrier(); // Ensures post-processing is completed
    // Read Results into TRIQS Containers
    readexpv(sp.Nz);
    readGF("G", G_w, gf_struct);
    readGF("F", F_w, gf_struct);
    readA("A", A_w, gf_struct);
    readA("B", B_w, gf_struct);
    readGF("SS", chi_SS_w, chi_struct);
    readGF("NN", chi_NN_w, chi_struct);

    // Post-Processing in C++ interface
    Sigma_w = (*F_w) / (*G_w);

    // Cleanup
    world.barrier(); // Ensures all processes have read the results before cleanup
    if (chdir("..") != 0) TRIQS_RUNTIME_ERROR << "failed to return from tempdir";
    world.barrier(); // Ensures all processes have chdired from temp dir before it is deleted
    if (world.rank() == 0) {
#ifdef NDEBUG
      // In production mode, we remove the temporary files.
      // In debug mode, we keep the temporary files for inspection.
//      std::error_code ec;
//      fs::remove_all(tempdir, ec);
//      if (ec) std::cout << "Warning: failed to remove the temporary directory." << std::endl;
//    Workaround:
      call("rm -rf " + tempdir, verbose, false);
#endif
    }
  }

  // -------------------------------------------------------------------------------

  void solver_core::set_params() {
    const constr_params_t &cp = constr_params;
    solve_params_t &sp        = *last_solve_params;
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
    if (sp.bandrescale < 0) // Make the NRG energy window correspond to the extent of the frequency mesh
      sp.bandrescale = cp.mesh_max;
    // Ensure the selected method is enabled. Other methods may be enabled as well, but only the output files for the selected method well be read-in by the nrglj-interface. 
    TRIQS_ASSERT2(sp.method == "fdm", "currently only method=fdm is supported"); // TODO
    if (sp.method == "fdm") { np.fdm = true; }
    if (sp.method == "dmnrg") { np.dmnrg = true; }
    if (sp.method == "cfs") { np.cfs = true; }
    if (sp.method == "finite") { np.finite = true; }
  }

  void solver_core::set_nrg_params(nrg_params_t const &nrg_params_) { nrg_params = nrg_params_; }

  // Creates input files for NRG iteration (param & data) by calling the "instantiate" script
  void solver_core::instantiate(double z, const std::string &taskdir) {
    std::cout << "Preparing z=" << z << " taskdir=" << taskdir << std::endl;
    TRIQS_ASSERT(world.rank() == 0); // should not run in parallel!
    generate_param_file(z);
    if (mkdir(taskdir.c_str(), 0755) != 0) TRIQS_RUNTIME_ERROR << "failed to mkdir taskdir " << taskdir;
    call("./instantiate " + taskdir, verbose);
  }
   
  // Solve the problem for a single value of the twist parameter z
  void solver_core::solve_one(const std::string &taskdir) {
    std::cout << "Solving taskdir=" << taskdir << " rank=" << world.rank() << std::endl;
    if (chdir(taskdir.c_str()) != 0) TRIQS_RUNTIME_ERROR << "failed to chdir to taskdir " << taskdir;
    // Solve the impurity model
    std::streambuf *old = std::cout.rdbuf();
    std::stringstream ss;
    if (!verbose) std::cout.rdbuf(ss.rdbuf());
    set_workdir("."); // may be overridden by NRG_WORKDIR in environment
    run_nrg_master();
    if (!verbose) std::cout.rdbuf(old);
    if (chdir("..") != 0) TRIQS_RUNTIME_ERROR << "failed to return from taskdir " << taskdir;
  }

  std::string solver_core::create_tempdir() {
    const std::string tempdir_template = "nrg_tempdir_XXXXXX";
    size_t len = tempdir_template.length()+1;
    auto x = std::make_unique<char[]>(len); // NOLINT
    strncpy(x.get(), tempdir_template.c_str(), len);
    if (auto w = mkdtemp(x.get())) // create a unique directory
      return w;
    else
      TRIQS_RUNTIME_ERROR << "Failed to create a directory for temporary files.";
  }
   
  void solver_core::check_model_params(const solve_params_t &sp) {
    std::istringstream ss(constr_params.params);
    do {
      std::string p;
      ss >> p;
      if (p != "" && sp.model_parameters.count(p) != 1)
        TRIQS_RUNTIME_ERROR << "Model parameter " << p << " has not been defined.";
    } while (ss);
  }

  void solver_core::generate_param_file(double z) {
    const constr_params_t &cp = constr_params;
    const solve_params_t &sp  = *last_solve_params;
    nrg_params_t &np          = nrg_params;
    np.z                      = z;
    // Automatically establish appropriate default values for the high-level interface
    set_params();
    // Check if all model parameters have been defined
    check_model_params(sp);
    // Generate the parameter file
    std::ofstream F("param");
    F << "[extra]" << std::endl;
    for (const auto &i : sp.model_parameters) F << i.first << "=" << i.second << std::endl;
    F << std::boolalpha; // important: we want to output true/false strings
    F << "[param]" << std::endl;
    F << "bandrescale=" << sp.bandrescale << std::endl; // !
    F << "model=" << cp.model << std::endl;
    F << "symtype=" << cp.symtype << std::endl;
    F << "Lambda=" << sp.Lambda << std::endl;
    F << "xmax=" << np.xmax << std::endl;
    F << "Nmax=" << np.Nmax << std::endl;
    F << "mMAX=" << np.mMAX << std::endl;
    F << "keep=" << sp.keep << std::endl;
    F << "keepenergy=" << sp.keepenergy << std::endl;
    F << "keepmin=" << sp.keepmin << std::endl;
    F << "T=" << sp.T << std::endl;
    F << "ops=" << cp.ops << std::endl;
    F << "specs=" << cp.specs << std::endl;
    F << "specd=" << cp.specd << std::endl;
    F << "spect=" << cp.spect << std::endl;
    F << "specq=" << cp.specq << std::endl;
    F << "specot=" << cp.specot << std::endl;
    F << "specgt=" << np.specgt << std::endl;
    F << "speci1t=" << np.speci1t << std::endl;
    F << "speci2t=" << np.speci2t << std::endl;
    F << "specchit=" << cp.specchit << std::endl;
    F << "specv3=" << cp.specv3 << std::endl;
    F << "v3mm=" << np.v3mm << std::endl;
    F << "dmnrg=" << np.dmnrg << std::endl;
    F << "cfs=" << np.cfs << std::endl;
    F << "fdm=" << np.fdm << std::endl;
    F << "fdmexpv=" << np.fdmexpv << std::endl;
    F << "dmnrgmats=" << np.dmnrgmats << std::endl;
    F << "fdmmats=" << np.fdmmats << std::endl;
    F << "mats=" << np.mats << std::endl;
    F << "alpha=" << sp.alpha << std::endl;
    F << "gamma=" << sp.gamma << std::endl; // ?
    F << "discretization=" << np.discretization << std::endl;
    F << "z=" << np.z << std::endl;
    F << "Nz=" << sp.Nz << std::endl; // !
    F << "polarized=" << cp.polarized << std::endl;
    F << "pol2x2=" << cp.pol2x2 << std::endl;
    F << "rungs=" << cp.rungs << std::endl;
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
    F << "data_has_rescaled_energies=false" << std::endl; // NRGLj-TRIQS interface uses unscaled energies
  }

  // -------------------------------------------------------------------------------

  // Function that writes a solver object to hdf5 file

  void h5_write(triqs::h5::group h5group, std::string subgroup_name, solver_core const &s) {
    auto grp = h5group.create_group(subgroup_name);
    h5_write_attribute(grp, "TRIQS_HDF5_data_scheme", solver_core::hdf5_scheme());
    h5_write_attribute(grp, "TRIQS_GIT_HASH", std::string(AS_STRING(TRIQS_GIT_HASH)));
    h5_write_attribute(grp, "NRGLJUBLJANA_INTERFACE_GIT_HASH", std::string(AS_STRING(NRGLJUBLJANA_INTERFACE_GIT_HASH)));
    h5_write(grp, "", s.result_set());
    h5_write(grp, "constr_params", s.constr_params);
    h5_write(grp, "nrg_params", s.nrg_params);
    h5_write(grp, "last_solve_params", s.last_solve_params);
    h5_write(grp, "Delta_w", s.Delta_w); // !!
  }

  // Function that constructs a solver object from an hdf5 file
  solver_core solver_core::h5_read_construct(triqs::h5::group h5group, std::string subgroup_name) {
    auto grp           = h5group.open_group(subgroup_name);
    auto constr_params = h5_read<constr_params_t>(grp, "constr_params");
    auto s             = solver_core{constr_params};
    h5_read(grp, "", s.result_set());
    h5_read(grp, "nrg_params", s.nrg_params);
    h5_read(grp, "last_solve_params", s.last_solve_params);
    h5_read(grp, "Delta_w", s.Delta_w);
    return s;
  }

  // Hilbert transform for refreq objects
  std::complex<double> hilbert_transform_refreq(const c_w_cvt &gf, std::complex<double> z){
    return hilbert_transform(gf, z);
  }

  matrix<std::complex<double>> hilbert_transform_elementwise(const m_w_cvt &gf, std::complex<double> z) {
    return hilbert_transform_elementwise(gf, z);
  }

} // namespace nrgljubljana_interface
