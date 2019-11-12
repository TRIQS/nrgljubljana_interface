/*******************************************************************************
 *
 * nrgljubljana_interface: A TRIQS based impurity solver
 *
 * Copyright (c) 2019 The Simons foundation
 *   authors: Nils Wentzell
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
#include "./post_process.hpp"

#include <nrg-lib.h>

namespace nrgljubljana_interface {

  solver_core::solver_core(constr_params_t const &p) : constr_params(p) {
  }

  // -------------------------------------------------------------------------------

  void solver_core::generate_param_file(constr_params_t const &cp,
                                        solve_params_t const &sp, 
					nrg_params_t const &np)
  {
    std::ofstream F("param");
    F << "[param]" << std::endl;
    F << "bandrescale=" << cp.bandrescale << std::endl; // related to mesh_max! autoconfig?
    F << "model=" << cp.problem << std::endl; // not required when using templates
    F << "symtype=QS" << std::endl; // from template database TODO
    F << "Lambda=" << sp.Lambda << std::endl;
//    F << "Nz=" << sp.Nz << std::endl;
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
    F << "broaden_max=" << cp.mesh_max << std::endl; // !
    F << "broaden_min=" << cp.mesh_min << std::endl; // !
    F << "broaden_ratio=" << cp.mesh_ratio << std::endl; // !
//    F << "broaden_max=" << np.broaden_max << std::endl;
//    F << "broaden_min=" << np.broaden_min << std::endl;
//    F << "broaden_ratio=" << np.broaden_ratio << std::endl;
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
    for (const auto &i : sp.model_parameters)
      F << i.first << "=" << i.second << std::endl;	
  }

  void solver_core::solve(solve_params_t const &solve_params) {

    last_solve_params = solve_params;

    if (world.rank() == 0)
      std::cout << "\n"
                   "NRGLJUBLJANA_INTERFACE Solver\n";

    // Reset the results
    container_set::operator=(container_set{});

    // Test if low-level paramerers are sensible as used in
    // the high-level interface
    if (np.discretization != "Z") {
    }
     
    // Automatically set (override) some low-level parameters
    if (sp.Tmin > 0) {
       np.Nmax = 0;
       auto scale = (1.0-1.0/sp.Lambda)/log(sp.Lambda)*pow(sp.Lambda
       while (scale(Nmax+1) >= sp.Tmin) np.Nmax++;
    }
    if (np.mMAX < 0)
       np.mMAX = max(80, 2*np.Nmax);
    if (np.xmax < 0)
       np.xmax = np.Nmax + 2.0;
     
    // Solve the impurity model
    generate_param_file(constr_params, solve_params, nrg_params);
    set_workdir(".");
     if (world.rank() == 0)
	run_nrg_master();
     else
	run_nrg_slave();
  }

  void solver_core::set_nrg_params(nrg_params_t const &nrg_params_) {
    nrg_params = nrg_params_;
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
      h5_write(grp, "last_solve_params", s.last_solve_params);
      h5_write(grp, "nrg_params", s.nrg_params);
    }

    // Function that constructs a solver object from an hdf5 file
    solver_core solver_core::h5_read_construct(triqs::h5::group h5group, std::string subgroup_name) {
      auto grp           = h5group.open_group(subgroup_name);
      auto constr_params = h5_read<constr_params_t>(grp, "constr_params");
      auto s             = solver_core{constr_params};
      h5_read(grp, "", s.result_set());
      h5_read(grp, "last_solve_params", s.last_solve_params);
      h5_read(grp, "nrg_params", s.nrg_params);
      return s;
    }
} // namespace nrgljubljana_interface
