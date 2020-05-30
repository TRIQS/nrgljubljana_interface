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
#include "./params.hpp"

namespace nrgljubljana_interface {

  void h5_write(h5::group h5group, std::string subgroup_name, constr_params_t const &cp) {
    auto grp = h5group.create_group(subgroup_name);
    h5_write(grp, "templatedir", cp.templatedir);
    h5_write(grp, "model", cp.model);
    h5_write(grp, "symtype", cp.symtype);
    h5_write(grp, "mesh_max", cp.mesh_max);
    h5_write(grp, "mesh_min", cp.mesh_min);
    h5_write(grp, "mesh_ratio", cp.mesh_ratio);
    h5_write(grp, "polarized", cp.polarized);
    h5_write(grp, "pol2x2", cp.pol2x2);
    h5_write(grp, "rungs", cp.rungs);
    h5_write(grp, "ops", cp.ops);
    h5_write(grp, "specs", cp.specs);
    h5_write(grp, "specd", cp.specd);
    h5_write(grp, "spect", cp.spect);
    h5_write(grp, "specq", cp.specq);
    h5_write(grp, "specot", cp.specot);
    h5_write(grp, "specchit", cp.specchit);
    h5_write(grp, "specv3", cp.specv3);
    h5_write(grp, "params", cp.params);
  }

  void h5_read(h5::group h5group, std::string subgroup_name, constr_params_t &cp) {
    auto grp = h5group.open_group(subgroup_name);
    h5_read(grp, "templatedir", cp.templatedir);
    h5_read(grp, "model", cp.model);
    h5_read(grp, "symtype", cp.symtype);
    h5_read(grp, "mesh_max", cp.mesh_max);
    h5_read(grp, "mesh_min", cp.mesh_min);
    h5_read(grp, "mesh_ratio", cp.mesh_ratio);
    h5_read(grp, "polarized", cp.polarized);
    h5_read(grp, "pol2x2", cp.pol2x2);
    h5_read(grp, "rungs", cp.rungs);
    h5_read(grp, "ops", cp.ops);
    h5_read(grp, "specs", cp.specs);
    h5_read(grp, "specd", cp.specd);
    h5_read(grp, "spect", cp.spect);
    h5_read(grp, "specq", cp.specq);
    h5_read(grp, "specot", cp.specot);
    h5_read(grp, "specchit", cp.specchit);
    h5_read(grp, "specv3", cp.specv3);
    h5_read(grp, "params", cp.params);
  }

  void h5_write(h5::group h5group, std::string subgroup_name, solve_params_t const &sp) {
    auto grp = h5group.create_group(subgroup_name);
    h5_write(grp, "Lambda", sp.Lambda);
    h5_write(grp, "Nz", sp.Nz);
    h5_write(grp, "Tmin", sp.Tmin);
    h5_write(grp, "keep", sp.keep);
    h5_write(grp, "keepenergy", sp.keepenergy);
    h5_write(grp, "keepmin", sp.keepmin);
    h5_write(grp, "T", sp.T);
    h5_write(grp, "alpha", sp.alpha);
    h5_write(grp, "gamma", sp.gamma);
    h5_write(grp, "method", sp.method);
    h5_write(grp, "bandrescale", sp.bandrescale);
    h5_write(grp, "model_parameters", sp.model_parameters);
  }

  void h5_read(h5::group h5group, std::string subgroup_name, solve_params_t &sp) {
    auto grp = h5group.open_group(subgroup_name);
    h5_read(grp, "Lambda", sp.Lambda);
    h5_read(grp, "Nz", sp.Nz);
    h5_read(grp, "Tmin", sp.Tmin);
    h5_read(grp, "keep", sp.keep);
    h5_read(grp, "keepenergy", sp.keepenergy);
    h5_read(grp, "keepmin", sp.keepmin);
    h5_read(grp, "T", sp.T);
    h5_read(grp, "alpha", sp.alpha);
    h5_read(grp, "gamma", sp.gamma);
    h5_read(grp, "method", sp.method);
    h5_read(grp, "bandrescale", sp.bandrescale);
    h5_read(grp, "model_parameters", sp.model_parameters);
  }

  void h5_write(h5::group h5group, std::string subgroup_name, nrg_params_t const &np) {
    auto grp = h5group.create_group(subgroup_name);
    h5_write(grp, "dmnrg", np.dmnrg);
    h5_write(grp, "cfs", np.cfs);
    h5_write(grp, "fdm", np.fdm);
    h5_write(grp, "fdmexpv", np.fdmexpv);
    h5_write(grp, "dmnrgmats", np.dmnrgmats);
    h5_write(grp, "fdmmats", np.fdmmats);
    h5_write(grp, "mats", np.mats);
    h5_write(grp, "specgt", np.specgt);
    h5_write(grp, "speci1t", np.speci1t);
    h5_write(grp, "speci2t", np.speci2t);
    h5_write(grp, "v3mm", np.v3mm);
    h5_write(grp, "mMAX", np.mMAX);
    h5_write(grp, "Nmax", np.Nmax);
    h5_write(grp, "xmax", np.xmax);
    h5_write(grp, "discretization", np.discretization);
    h5_write(grp, "z", np.z);
    h5_write(grp, "tri", np.tri);
    h5_write(grp, "preccpp", np.preccpp);
    h5_write(grp, "diag", np.diag);
    h5_write(grp, "diagratio", np.diagratio);
    h5_write(grp, "dsyevrlimit", np.dsyevrlimit);
    h5_write(grp, "zheevrlimit", np.zheevrlimit);
    h5_write(grp, "restart", np.restart);
    h5_write(grp, "restartfactor", np.restartfactor);
    h5_write(grp, "safeguard", np.safeguard);
    h5_write(grp, "safeguardmax", np.safeguardmax);
    h5_write(grp, "fixeps", np.fixeps);
    h5_write(grp, "betabar", np.betabar);
    h5_write(grp, "gtp", np.gtp);
    h5_write(grp, "chitp", np.chitp);
    h5_write(grp, "finite", np.finite);
    h5_write(grp, "cfsgt", np.cfsgt);
    h5_write(grp, "cfsls", np.cfsls);
    h5_write(grp, "fdmgt", np.fdmgt);
    h5_write(grp, "fdmls", np.fdmls);
    h5_write(grp, "fdmexpvn", np.fdmexpvn);
    h5_write(grp, "finitemats", np.finitemats);
    h5_write(grp, "dm", np.dm);
    h5_write(grp, "broaden_min_ratio", np.broaden_min_ratio);
    h5_write(grp, "omega0", np.omega0);
    h5_write(grp, "omega0_ratio", np.omega0_ratio);
    h5_write(grp, "diagth", np.diagth);
    h5_write(grp, "substeps", np.substeps);
    h5_write(grp, "strategy", np.strategy);
    h5_write(grp, "Ninit", np.Ninit);
    h5_write(grp, "reim", np.reim);
    h5_write(grp, "dumpannotated", np.dumpannotated);
    h5_write(grp, "dumpabs", np.dumpabs);
    h5_write(grp, "dumpscaled", np.dumpscaled);
    h5_write(grp, "dumpprecision", np.dumpprecision);
    h5_write(grp, "dumpgroups", np.dumpgroups);
    h5_write(grp, "grouptol", np.grouptol);
    h5_write(grp, "dumpdiagonal", np.dumpdiagonal);
    h5_write(grp, "savebins", np.savebins);
    h5_write(grp, "broaden", np.broaden);
    h5_write(grp, "emin", np.emin);
    h5_write(grp, "emax", np.emax);
    h5_write(grp, "bins", np.bins);
    h5_write(grp, "accumulation", np.accumulation);
    h5_write(grp, "linstep", np.linstep);
    h5_write(grp, "discard_trim", np.discard_trim);
    h5_write(grp, "discard_immediately", np.discard_immediately);
    h5_write(grp, "goodE", np.goodE);
    h5_write(grp, "NN1", np.NN1);
    h5_write(grp, "NN2even", np.NN2even);
    h5_write(grp, "NN2avg", np.NN2avg);
    h5_write(grp, "NNtanh", np.NNtanh);
    h5_write(grp, "width_td", np.width_td);
    h5_write(grp, "width_custom", np.width_custom);
    h5_write(grp, "prec_td", np.prec_td);
    h5_write(grp, "prec_custom", np.prec_custom);
    h5_write(grp, "prec_xy", np.prec_xy);
    h5_write(grp, "resume", np.resume);
    h5_write(grp, "log", np.log);
    h5_write(grp, "logall", np.logall);
    h5_write(grp, "done", np.done);
    h5_write(grp, "calc0", np.calc0);
    h5_write(grp, "lastall", np.lastall);
    h5_write(grp, "lastalloverride", np.lastalloverride);
    h5_write(grp, "dumpsubspaces", np.dumpsubspaces);
    h5_write(grp, "dump_f", np.dump_f);
    h5_write(grp, "dumpenergies", np.dumpenergies);
    h5_write(grp, "logenumber", np.logenumber);
    h5_write(grp, "stopafter", np.stopafter);
    h5_write(grp, "forcestop", np.forcestop);
    h5_write(grp, "removefiles", np.removefiles);
    h5_write(grp, "noimag", np.noimag);
    h5_write(grp, "checksumrules", np.checksumrules);
    h5_write(grp, "checkdiag", np.checkdiag);
    h5_write(grp, "checkrho", np.checkrho);
  }

  void h5_read(h5::group h5group, std::string subgroup_name, nrg_params_t &np) {
    auto grp = h5group.open_group(subgroup_name);
    h5_read(grp, "dmnrg", np.dmnrg);
    h5_read(grp, "cfs", np.cfs);
    h5_read(grp, "fdm", np.fdm);
    h5_read(grp, "fdmexpv", np.fdmexpv);
    h5_read(grp, "dmnrgmats", np.dmnrgmats);
    h5_read(grp, "fdmmats", np.fdmmats);
    h5_read(grp, "mats", np.mats);
    h5_read(grp, "specgt", np.specgt);
    h5_read(grp, "speci1t", np.speci1t);
    h5_read(grp, "speci2t", np.speci2t);
    h5_read(grp, "v3mm", np.v3mm);
    h5_read(grp, "mMAX", np.mMAX);
    h5_read(grp, "Nmax", np.Nmax);
    h5_read(grp, "xmax", np.xmax);
    h5_read(grp, "discretization", np.discretization);
    h5_read(grp, "z", np.z);
    h5_read(grp, "tri", np.tri);
    h5_read(grp, "preccpp", np.preccpp);
    h5_read(grp, "diag", np.diag);
    h5_read(grp, "diagratio", np.diagratio);
    h5_read(grp, "dsyevrlimit", np.dsyevrlimit);
    h5_read(grp, "zheevrlimit", np.zheevrlimit);
    h5_read(grp, "restart", np.restart);
    h5_read(grp, "restartfactor", np.restartfactor);
    h5_read(grp, "safeguard", np.safeguard);
    h5_read(grp, "safeguardmax", np.safeguardmax);
    h5_read(grp, "fixeps", np.fixeps);
    h5_read(grp, "betabar", np.betabar);
    h5_read(grp, "gtp", np.gtp);
    h5_read(grp, "chitp", np.chitp);
    h5_read(grp, "finite", np.finite);
    h5_read(grp, "cfsgt", np.cfsgt);
    h5_read(grp, "cfsls", np.cfsls);
    h5_read(grp, "fdmgt", np.fdmgt);
    h5_read(grp, "fdmls", np.fdmls);
    h5_read(grp, "fdmexpvn", np.fdmexpvn);
    h5_read(grp, "finitemats", np.finitemats);
    h5_read(grp, "dm", np.dm);
    h5_read(grp, "broaden_min_ratio", np.broaden_min_ratio);
    h5_read(grp, "omega0", np.omega0);
    h5_read(grp, "omega0_ratio", np.omega0_ratio);
    h5_read(grp, "diagth", np.diagth);
    h5_read(grp, "substeps", np.substeps);
    h5_read(grp, "strategy", np.strategy);
    h5_read(grp, "Ninit", np.Ninit);
    h5_read(grp, "reim", np.reim);
    h5_read(grp, "dumpannotated", np.dumpannotated);
    h5_read(grp, "dumpabs", np.dumpabs);
    h5_read(grp, "dumpscaled", np.dumpscaled);
    h5_read(grp, "dumpprecision", np.dumpprecision);
    h5_read(grp, "dumpgroups", np.dumpgroups);
    h5_read(grp, "grouptol", np.grouptol);
    h5_read(grp, "dumpdiagonal", np.dumpdiagonal);
    h5_read(grp, "savebins", np.savebins);
    h5_read(grp, "broaden", np.broaden);
    h5_read(grp, "emin", np.emin);
    h5_read(grp, "emax", np.emax);
    h5_read(grp, "bins", np.bins);
    h5_read(grp, "accumulation", np.accumulation);
    h5_read(grp, "linstep", np.linstep);
    h5_read(grp, "discard_trim", np.discard_trim);
    h5_read(grp, "discard_immediately", np.discard_immediately);
    h5_read(grp, "goodE", np.goodE);
    h5_read(grp, "NN1", np.NN1);
    h5_read(grp, "NN2even", np.NN2even);
    h5_read(grp, "NN2avg", np.NN2avg);
    h5_read(grp, "NNtanh", np.NNtanh);
    h5_read(grp, "width_td", np.width_td);
    h5_read(grp, "width_custom", np.width_custom);
    h5_read(grp, "prec_td", np.prec_td);
    h5_read(grp, "prec_custom", np.prec_custom);
    h5_read(grp, "prec_xy", np.prec_xy);
    h5_read(grp, "resume", np.resume);
    h5_read(grp, "log", np.log);
    h5_read(grp, "logall", np.logall);
    h5_read(grp, "done", np.done);
    h5_read(grp, "calc0", np.calc0);
    h5_read(grp, "lastall", np.lastall);
    h5_read(grp, "lastalloverride", np.lastalloverride);
    h5_read(grp, "dumpsubspaces", np.dumpsubspaces);
    h5_read(grp, "dump_f", np.dump_f);
    h5_read(grp, "dumpenergies", np.dumpenergies);
    h5_read(grp, "logenumber", np.logenumber);
    h5_read(grp, "stopafter", np.stopafter);
    h5_read(grp, "forcestop", np.forcestop);
    h5_read(grp, "removefiles", np.removefiles);
    h5_read(grp, "noimag", np.noimag);
    h5_read(grp, "checksumrules", np.checksumrules);
    h5_read(grp, "checkdiag", np.checkdiag);
    h5_read(grp, "checkrho", np.checkrho);
  }

} // namespace nrgljubljana_interface
