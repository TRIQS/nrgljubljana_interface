
// C.f. https://numpy.org/doc/1.21/reference/c-api/array.html#importing-the-api
#define PY_ARRAY_UNIQUE_SYMBOL _cpp2py_ARRAY_API
#ifndef CLAIR_C2PY_WRAP_GEN
#ifdef __clang__
// #pragma clang diagnostic ignored "-W#warnings"
#endif
#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#pragma GCC diagnostic ignored "-Wcast-function-type"
#pragma GCC diagnostic ignored "-Wcpp"
#endif

#define C2PY_VERSION_MAJOR 0
#define C2PY_VERSION_MINOR 1

#include <c2py/c2py.hpp>
#include <c2py/serialization/h5.hpp>

using c2py::operator""_a;

// ==================== enums =====================

// ==================== module classes =====================

// --------- class _c2py_cls_0 -----------
using _c2py_cls_0                                            = nrgljubljana_interface::constr_params_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_0>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_0> = "nrgljubljana_interface.solver_core.ConstrParamsT";

static int synth_constructor_0(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing nrgljubljana_interface::constr_params_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_0> *)self)->_c = new _c2py_cls_0{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing nrgljubljana_interface::constr_params_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_0> *)self)->_c);
  de("templatedir", self_c.templatedir, true);
  de("model", self_c.model, true);
  de("symtype", self_c.symtype, true);
  de("mesh_max", self_c.mesh_max, true);
  de("mesh_min", self_c.mesh_min, true);
  de("mesh_ratio", self_c.mesh_ratio, true);
  de("polarized", self_c.polarized, true);
  de("pol2x2", self_c.pol2x2, true);
  de("rungs", self_c.rungs, true);
  de("ops", self_c.ops, true);
  de("specs", self_c.specs, true);
  de("specd", self_c.specd, true);
  de("spect", self_c.spect, true);
  de("specq", self_c.specq, true);
  de("specot", self_c.specot, true);
  de("specchit", self_c.specchit, true);
  de("specv3", self_c.specv3, true);
  de("params", self_c.params, true);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_0> = synth_constructor_0;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_0> = c2py::replace_tags(
   R"DOC(Synthesized constructor with the following keyword arguments:

Parameters
----------
templatedir : {par_0}, default=NRGIF_TEMPLATE_DIR

model : {par_1}, default="SIAM"

symtype : {par_2}, default="QS"

mesh_max : {par_3}, default=10

mesh_min : {par_4}, default=1e-4

mesh_ratio : {par_5}, default=1.05

polarized : {par_6}, default=false

pol2x2 : {par_7}, default=false

rungs : {par_8}, default=false

ops : {par_9}, default=""

specs : {par_10}, default=""

specd : {par_11}, default=""

spect : {par_12}, default=""

specq : {par_13}, default=""

specot : {par_14}, default=""

specchit : {par_15}, default=""

specv3 : {par_16}, default=""

params : {par_17}, default=""

)DOC",
   "par",
   {c2py::python_typename<std::string>(), c2py::python_typename<std::string>(), c2py::python_typename<std::string>(), c2py::python_typename<double>(),
    c2py::python_typename<double>(), c2py::python_typename<double>(), c2py::python_typename<bool>(), c2py::python_typename<bool>(),
    c2py::python_typename<bool>(), c2py::python_typename<std::string>(), c2py::python_typename<std::string>(), c2py::python_typename<std::string>(),
    c2py::python_typename<std::string>(), c2py::python_typename<std::string>(), c2py::python_typename<std::string>(),
    c2py::python_typename<std::string>(), c2py::python_typename<std::string>(), c2py::python_typename<std::string>()});
// get_model_dir
static auto const _c2py_fun_0 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_0 const &self) -> decltype(auto) { return self.get_model_dir(); }, "self")};

static const auto _c2py_doc_0 = _c2py_fun_0.doc(R"DOC(
Resolve the template directory for the chosen model and symmetry.

Returns ``templatedir/model/symtype``, with the ``$NRGIF_TEMPLATE_DIR``
environment variable overriding ``templatedir`` when set.

Returns
-------
{ret_0}
   Absolute path to the model/symmetry template directory.
)DOC",
                                                {}, {c2py::python_typename<std::string>()});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_0>[] = {
   {"get_model_dir", (PyCFunction)c2py::pyfkw<_c2py_fun_0>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_0.c_str()},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_0  = R"DOC(Path to the template library (defaults to the bundled templates).)DOC";
constexpr auto _c2py_doc_member_1  = R"DOC(Impurity model to solve (selects a template directory).)DOC";
constexpr auto _c2py_doc_member_2  = R"DOC(Symmetry type (NRGLjubljana symmetry code, e.g. QS, QSZ, ISO).)DOC";
constexpr auto _c2py_doc_member_3  = R"DOC(Maximum frequency of the logarithmic mesh.)DOC";
constexpr auto _c2py_doc_member_4  = R"DOC(Minimum frequency of the logarithmic mesh.)DOC";
constexpr auto _c2py_doc_member_5  = R"DOC(Common ratio of the geometric (logarithmic) frequency mesh.)DOC";
constexpr auto _c2py_doc_member_6  = R"DOC(Use a spin-polarized Wilson chain.)DOC";
constexpr auto _c2py_doc_member_7  = R"DOC(Use a 2x2 spin structure in the Wilson chain.)DOC";
constexpr auto _c2py_doc_member_8  = R"DOC(Include channel-mixing terms in the Wilson chain.)DOC";
constexpr auto _c2py_doc_member_9  = R"DOC(Operators whose expectation values are to be calculated.)DOC";
constexpr auto _c2py_doc_member_10 = R"DOC(Spectral functions of singlet operators to compute.)DOC";
constexpr auto _c2py_doc_member_11 = R"DOC(Spectral functions of doublet operators to compute.)DOC";
constexpr auto _c2py_doc_member_12 = R"DOC(Spectral functions of triplet operators to compute.)DOC";
constexpr auto _c2py_doc_member_13 = R"DOC(Spectral functions of quadruplet operators to compute.)DOC";
constexpr auto _c2py_doc_member_14 = R"DOC(Spectral functions of orbital-triplet operators to compute.)DOC";
constexpr auto _c2py_doc_member_15 = R"DOC(Susceptibilities to compute.)DOC";
constexpr auto _c2py_doc_member_16 = R"DOC(3-leg vertex functions to compute.)DOC";
constexpr auto _c2py_doc_member_17 = R"DOC(List of model parameters that need to be specified.)DOC";
static PyObject *prop_get_dict_0(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_0> *)self)->_c);
  c2py::pydict dic;
  dic["templatedir"] = self_c.templatedir;
  dic["model"]       = self_c.model;
  dic["symtype"]     = self_c.symtype;
  dic["mesh_max"]    = self_c.mesh_max;
  dic["mesh_min"]    = self_c.mesh_min;
  dic["mesh_ratio"]  = self_c.mesh_ratio;
  dic["polarized"]   = self_c.polarized;
  dic["pol2x2"]      = self_c.pol2x2;
  dic["rungs"]       = self_c.rungs;
  dic["ops"]         = self_c.ops;
  dic["specs"]       = self_c.specs;
  dic["specd"]       = self_c.specd;
  dic["spect"]       = self_c.spect;
  dic["specq"]       = self_c.specq;
  dic["specot"]      = self_c.specot;
  dic["specchit"]    = self_c.specchit;
  dic["specv3"]      = self_c.specv3;
  dic["params"]      = self_c.params;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_0>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_0::templatedir, _c2py_cls_0>("templatedir", _c2py_doc_member_0),
   c2py::getsetdef_from_member<&_c2py_cls_0::model, _c2py_cls_0>("model", _c2py_doc_member_1),
   c2py::getsetdef_from_member<&_c2py_cls_0::symtype, _c2py_cls_0>("symtype", _c2py_doc_member_2),
   c2py::getsetdef_from_member<&_c2py_cls_0::mesh_max, _c2py_cls_0>("mesh_max", _c2py_doc_member_3),
   c2py::getsetdef_from_member<&_c2py_cls_0::mesh_min, _c2py_cls_0>("mesh_min", _c2py_doc_member_4),
   c2py::getsetdef_from_member<&_c2py_cls_0::mesh_ratio, _c2py_cls_0>("mesh_ratio", _c2py_doc_member_5),
   c2py::getsetdef_from_member<&_c2py_cls_0::polarized, _c2py_cls_0>("polarized", _c2py_doc_member_6),
   c2py::getsetdef_from_member<&_c2py_cls_0::pol2x2, _c2py_cls_0>("pol2x2", _c2py_doc_member_7),
   c2py::getsetdef_from_member<&_c2py_cls_0::rungs, _c2py_cls_0>("rungs", _c2py_doc_member_8),
   c2py::getsetdef_from_member<&_c2py_cls_0::ops, _c2py_cls_0>("ops", _c2py_doc_member_9),
   c2py::getsetdef_from_member<&_c2py_cls_0::specs, _c2py_cls_0>("specs", _c2py_doc_member_10),
   c2py::getsetdef_from_member<&_c2py_cls_0::specd, _c2py_cls_0>("specd", _c2py_doc_member_11),
   c2py::getsetdef_from_member<&_c2py_cls_0::spect, _c2py_cls_0>("spect", _c2py_doc_member_12),
   c2py::getsetdef_from_member<&_c2py_cls_0::specq, _c2py_cls_0>("specq", _c2py_doc_member_13),
   c2py::getsetdef_from_member<&_c2py_cls_0::specot, _c2py_cls_0>("specot", _c2py_doc_member_14),
   c2py::getsetdef_from_member<&_c2py_cls_0::specchit, _c2py_cls_0>("specchit", _c2py_doc_member_15),
   c2py::getsetdef_from_member<&_c2py_cls_0::specv3, _c2py_cls_0>("specv3", _c2py_doc_member_16),
   c2py::getsetdef_from_member<&_c2py_cls_0::params, _c2py_cls_0>("params", _c2py_doc_member_17),
   {"__dict__", (getter)prop_get_dict_0, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_0> =
   R"DOC(Construction parameters for the NRGLjubljana solver.)DOC" + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_0>;
// --------- class _c2py_cls_1 -----------
using _c2py_cls_1                                            = nrgljubljana_interface::solve_params_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_1>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_1> = "nrgljubljana_interface.solver_core.SolveParamsT";

static int synth_constructor_1(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing nrgljubljana_interface::solve_params_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_1> *)self)->_c = new _c2py_cls_1{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing nrgljubljana_interface::solve_params_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_1> *)self)->_c);
  de("Lambda", self_c.Lambda, true);
  de("Nz", self_c.Nz, true);
  de("Tmin", self_c.Tmin, true);
  de("keep", self_c.keep, true);
  de("keepenergy", self_c.keepenergy, true);
  de("keepmin", self_c.keepmin, true);
  de("T", self_c.T, true);
  de("alpha", self_c.alpha, true);
  de("gamma", self_c.gamma, true);
  de("method", self_c.method, true);
  de("bandrescale", self_c.bandrescale, true);
  de("model_parameters", self_c.model_parameters, false);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_1> = synth_constructor_1;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_1> =
   c2py::replace_tags(R"DOC(Synthesized constructor with the following keyword arguments:

Parameters
----------
model_parameters : {par_0}

Lambda : {par_1}, default=2.0

Nz : {par_2}, default=1

Tmin : {par_3}, default=1e-4

keep : {par_4}, default=100

keepenergy : {par_5}, default=-1.0

keepmin : {par_6}, default=0

T : {par_7}, default=0.001

alpha : {par_8}, default=0.3

gamma : {par_9}, default=0.2

method : {par_10}, default="fdm"

bandrescale : {par_11}, default=-1.0

)DOC",
                      "par",
                      {c2py::python_typename<std::map<std::string, double>>(), c2py::python_typename<double>(), c2py::python_typename<int>(),
                       c2py::python_typename<double>(), c2py::python_typename<unsigned long>(), c2py::python_typename<double>(),
                       c2py::python_typename<unsigned long>(), c2py::python_typename<double>(), c2py::python_typename<double>(),
                       c2py::python_typename<double>(), c2py::python_typename<std::string>(), c2py::python_typename<double>()});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_1>[] = {

   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_18 = R"DOC(Logarithmic discretization parameter.)DOC";
constexpr auto _c2py_doc_member_19 = R"DOC(Number of discretization meshes (interleaved twist parameters z).)DOC";
constexpr auto _c2py_doc_member_20 = R"DOC(Lowest energy scale on the Wilson chain.)DOC";
constexpr auto _c2py_doc_member_21 = R"DOC(Maximum number of states to keep at each NRG step.)DOC";
constexpr auto _c2py_doc_member_22 = R"DOC(Cut-off energy for truncation.)DOC";
constexpr auto _c2py_doc_member_23 = R"DOC(Minimum number of states to keep at each NRG step.)DOC";
constexpr auto _c2py_doc_member_24 = R"DOC(Temperature, :math:`k_B T / D`.)DOC";
constexpr auto _c2py_doc_member_25 = R"DOC(Width of the logarithmic gaussian used for broadening.)DOC";
constexpr auto _c2py_doc_member_26 = R"DOC(Parameter for the Gaussian convolution step.)DOC";
constexpr auto _c2py_doc_member_27 = R"DOC(Method for calculating the dynamical quantities.)DOC";
constexpr auto _c2py_doc_member_28 =
   R"DOC(Band rescaling factor (half-width of the support of the hybridisation function); set to mesh_max if negative.)DOC";
constexpr auto _c2py_doc_member_29 = R"DOC(Model parameters (name to value map, e.g. U1, eps1).)DOC";
static PyObject *prop_get_dict_1(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_1> *)self)->_c);
  c2py::pydict dic;
  dic["Lambda"]           = self_c.Lambda;
  dic["Nz"]               = self_c.Nz;
  dic["Tmin"]             = self_c.Tmin;
  dic["keep"]             = self_c.keep;
  dic["keepenergy"]       = self_c.keepenergy;
  dic["keepmin"]          = self_c.keepmin;
  dic["T"]                = self_c.T;
  dic["alpha"]            = self_c.alpha;
  dic["gamma"]            = self_c.gamma;
  dic["method"]           = self_c.method;
  dic["bandrescale"]      = self_c.bandrescale;
  dic["model_parameters"] = self_c.model_parameters;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_1>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_1::Lambda, _c2py_cls_1>("Lambda", _c2py_doc_member_18),
   c2py::getsetdef_from_member<&_c2py_cls_1::Nz, _c2py_cls_1>("Nz", _c2py_doc_member_19),
   c2py::getsetdef_from_member<&_c2py_cls_1::Tmin, _c2py_cls_1>("Tmin", _c2py_doc_member_20),
   c2py::getsetdef_from_member<&_c2py_cls_1::keep, _c2py_cls_1>("keep", _c2py_doc_member_21),
   c2py::getsetdef_from_member<&_c2py_cls_1::keepenergy, _c2py_cls_1>("keepenergy", _c2py_doc_member_22),
   c2py::getsetdef_from_member<&_c2py_cls_1::keepmin, _c2py_cls_1>("keepmin", _c2py_doc_member_23),
   c2py::getsetdef_from_member<&_c2py_cls_1::T, _c2py_cls_1>("T", _c2py_doc_member_24),
   c2py::getsetdef_from_member<&_c2py_cls_1::alpha, _c2py_cls_1>("alpha", _c2py_doc_member_25),
   c2py::getsetdef_from_member<&_c2py_cls_1::gamma, _c2py_cls_1>("gamma", _c2py_doc_member_26),
   c2py::getsetdef_from_member<&_c2py_cls_1::method, _c2py_cls_1>("method", _c2py_doc_member_27),
   c2py::getsetdef_from_member<&_c2py_cls_1::bandrescale, _c2py_cls_1>("bandrescale", _c2py_doc_member_28),
   c2py::getsetdef_from_member<&_c2py_cls_1::model_parameters, _c2py_cls_1>("model_parameters", _c2py_doc_member_29),
   {"__dict__", (getter)prop_get_dict_1, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_1> =
   R"DOC(Parameters for the solve() method.)DOC" + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_1>;
// --------- class _c2py_cls_2 -----------
using _c2py_cls_2                                            = nrgljubljana_interface::nrg_params_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_2>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_2> = "nrgljubljana_interface.solver_core.NrgParamsT";

static int synth_constructor_2(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing nrgljubljana_interface::nrg_params_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_2> *)self)->_c = new _c2py_cls_2{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError, ("Error in constructing nrgljubljana_interface::nrg_params_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_2> *)self)->_c);
  de("dmnrg", self_c.dmnrg, true);
  de("cfs", self_c.cfs, true);
  de("fdm", self_c.fdm, true);
  de("fdmexpv", self_c.fdmexpv, true);
  de("dmnrgmats", self_c.dmnrgmats, true);
  de("fdmmats", self_c.fdmmats, true);
  de("mats", self_c.mats, true);
  de("specgt", self_c.specgt, true);
  de("speci1t", self_c.speci1t, true);
  de("speci2t", self_c.speci2t, true);
  de("v3mm", self_c.v3mm, true);
  de("mMAX", self_c.mMAX, true);
  de("Nmax", self_c.Nmax, true);
  de("xmax", self_c.xmax, true);
  de("discretization", self_c.discretization, true);
  de("z", self_c.z, true);
  de("tri", self_c.tri, true);
  de("preccpp", self_c.preccpp, true);
  de("diag", self_c.diag, true);
  de("diagratio", self_c.diagratio, true);
  de("dsyevrlimit", self_c.dsyevrlimit, true);
  de("zheevrlimit", self_c.zheevrlimit, true);
  de("restart", self_c.restart, true);
  de("restartfactor", self_c.restartfactor, true);
  de("safeguard", self_c.safeguard, true);
  de("safeguardmax", self_c.safeguardmax, true);
  de("fixeps", self_c.fixeps, true);
  de("betabar", self_c.betabar, true);
  de("gtp", self_c.gtp, true);
  de("chitp", self_c.chitp, true);
  de("finite", self_c.finite, true);
  de("cfsgt", self_c.cfsgt, true);
  de("cfsls", self_c.cfsls, true);
  de("fdmgt", self_c.fdmgt, true);
  de("fdmls", self_c.fdmls, true);
  de("fdmexpvn", self_c.fdmexpvn, true);
  de("finitemats", self_c.finitemats, true);
  de("dm", self_c.dm, true);
  de("broaden_min_ratio", self_c.broaden_min_ratio, true);
  de("omega0", self_c.omega0, true);
  de("omega0_ratio", self_c.omega0_ratio, true);
  de("diagth", self_c.diagth, true);
  de("substeps", self_c.substeps, true);
  de("strategy", self_c.strategy, true);
  de("Ninit", self_c.Ninit, true);
  de("reim", self_c.reim, true);
  de("dumpannotated", self_c.dumpannotated, true);
  de("dumpabs", self_c.dumpabs, true);
  de("dumpscaled", self_c.dumpscaled, true);
  de("dumpprecision", self_c.dumpprecision, true);
  de("dumpgroups", self_c.dumpgroups, true);
  de("grouptol", self_c.grouptol, true);
  de("dumpdiagonal", self_c.dumpdiagonal, true);
  de("savebins", self_c.savebins, true);
  de("broaden", self_c.broaden, true);
  de("emin", self_c.emin, true);
  de("emax", self_c.emax, true);
  de("bins", self_c.bins, true);
  de("accumulation", self_c.accumulation, true);
  de("linstep", self_c.linstep, true);
  de("discard_trim", self_c.discard_trim, true);
  de("discard_immediately", self_c.discard_immediately, true);
  de("goodE", self_c.goodE, true);
  de("NN1", self_c.NN1, true);
  de("NN2even", self_c.NN2even, true);
  de("NN2avg", self_c.NN2avg, true);
  de("NNtanh", self_c.NNtanh, true);
  de("width_td", self_c.width_td, true);
  de("width_custom", self_c.width_custom, true);
  de("prec_td", self_c.prec_td, true);
  de("prec_custom", self_c.prec_custom, true);
  de("prec_xy", self_c.prec_xy, true);
  de("resume", self_c.resume, true);
  de("log", self_c.log, true);
  de("logall", self_c.logall, true);
  de("done", self_c.done, true);
  de("calc0", self_c.calc0, true);
  de("lastall", self_c.lastall, true);
  de("lastalloverride", self_c.lastalloverride, true);
  de("dumpsubspaces", self_c.dumpsubspaces, true);
  de("dump_f", self_c.dump_f, true);
  de("dumpenergies", self_c.dumpenergies, true);
  de("logenumber", self_c.logenumber, true);
  de("stopafter", self_c.stopafter, true);
  de("forcestop", self_c.forcestop, true);
  de("removefiles", self_c.removefiles, true);
  de("noimag", self_c.noimag, true);
  de("checksumrules", self_c.checksumrules, true);
  de("checkdiag", self_c.checkdiag, true);
  de("checkrho", self_c.checkrho, true);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_2> = synth_constructor_2;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_2> =
   c2py::replace_tags(R"DOC(Synthesized constructor with the following keyword arguments:

Parameters
----------
dmnrg : {par_0}, default=false

cfs : {par_1}, default=false

fdm : {par_2}, default=true

fdmexpv : {par_3}, default=true

dmnrgmats : {par_4}, default=false

fdmmats : {par_5}, default=false

mats : {par_6}, default=100

specgt : {par_7}, default=""

speci1t : {par_8}, default=""

speci2t : {par_9}, default=""

v3mm : {par_10}, default=false

mMAX : {par_11}, default=-1

Nmax : {par_12}, default=-1

xmax : {par_13}, default=-1.0

discretization : {par_14}, default="Z"

z : {par_15}, default=1.0

tri : {par_16}, default="old"

preccpp : {par_17}, default=2000

diag : {par_18}, default="default"

diagratio : {par_19}, default=1.0

dsyevrlimit : {par_20}, default=100

zheevrlimit : {par_21}, default=100

restart : {par_22}, default=true

restartfactor : {par_23}, default=2.0

safeguard : {par_24}, default=1e-5

safeguardmax : {par_25}, default=200

fixeps : {par_26}, default=1e-15

betabar : {par_27}, default=1.0

gtp : {par_28}, default=0.7

chitp : {par_29}, default=1.0

finite : {par_30}, default=false

cfsgt : {par_31}, default=false

cfsls : {par_32}, default=false

fdmgt : {par_33}, default=false

fdmls : {par_34}, default=false

fdmexpvn : {par_35}, default=0

finitemats : {par_36}, default=false

dm : {par_37}, default=false

broaden_min_ratio : {par_38}, default=3.0

omega0 : {par_39}, default=-1.0

omega0_ratio : {par_40}, default=1.0

diagth : {par_41}, default=1

substeps : {par_42}, default=false

strategy : {par_43}, default="kept"

Ninit : {par_44}, default=0

reim : {par_45}, default=false

dumpannotated : {par_46}, default=0

dumpabs : {par_47}, default=false

dumpscaled : {par_48}, default=true

dumpprecision : {par_49}, default=8

dumpgroups : {par_50}, default=true

grouptol : {par_51}, default=1e-6

dumpdiagonal : {par_52}, default=0

savebins : {par_53}, default=true

broaden : {par_54}, default=false

emin : {par_55}, default=-1.0

emax : {par_56}, default=-1.0

bins : {par_57}, default=1000

accumulation : {par_58}, default=0.0

linstep : {par_59}, default=0

discard_trim : {par_60}, default=1e-16

discard_immediately : {par_61}, default=1e-16

goodE : {par_62}, default=2.0

NN1 : {par_63}, default=false

NN2even : {par_64}, default=true

NN2avg : {par_65}, default=false

NNtanh : {par_66}, default=0.0

width_td : {par_67}, default=16

width_custom : {par_68}, default=16

prec_td : {par_69}, default=10

prec_custom : {par_70}, default=10

prec_xy : {par_71}, default=10

resume : {par_72}, default=false

log : {par_73}, default=""

logall : {par_74}, default=false

done : {par_75}, default=true

calc0 : {par_76}, default=true

lastall : {par_77}, default=false

lastalloverride : {par_78}, default=false

dumpsubspaces : {par_79}, default=false

dump_f : {par_80}, default=false

dumpenergies : {par_81}, default=false

logenumber : {par_82}, default=10

stopafter : {par_83}, default=""

forcestop : {par_84}, default=-1

removefiles : {par_85}, default=true

noimag : {par_86}, default=true

checksumrules : {par_87}, default=false

checkdiag : {par_88}, default=false

checkrho : {par_89}, default=false

)DOC",
                      "par", {c2py::python_typename<bool>(),          c2py::python_typename<bool>(),          c2py::python_typename<bool>(),
                              c2py::python_typename<bool>(),          c2py::python_typename<bool>(),          c2py::python_typename<bool>(),
                              c2py::python_typename<unsigned long>(), c2py::python_typename<std::string>(),   c2py::python_typename<std::string>(),
                              c2py::python_typename<std::string>(),   c2py::python_typename<bool>(),          c2py::python_typename<int>(),
                              c2py::python_typename<int>(),           c2py::python_typename<double>(),        c2py::python_typename<std::string>(),
                              c2py::python_typename<double>(),        c2py::python_typename<std::string>(),   c2py::python_typename<unsigned long>(),
                              c2py::python_typename<std::string>(),   c2py::python_typename<double>(),        c2py::python_typename<unsigned long>(),
                              c2py::python_typename<unsigned long>(), c2py::python_typename<bool>(),          c2py::python_typename<double>(),
                              c2py::python_typename<double>(),        c2py::python_typename<unsigned long>(), c2py::python_typename<double>(),
                              c2py::python_typename<double>(),        c2py::python_typename<double>(),        c2py::python_typename<double>(),
                              c2py::python_typename<bool>(),          c2py::python_typename<bool>(),          c2py::python_typename<bool>(),
                              c2py::python_typename<bool>(),          c2py::python_typename<bool>(),          c2py::python_typename<unsigned long>(),
                              c2py::python_typename<bool>(),          c2py::python_typename<bool>(),          c2py::python_typename<double>(),
                              c2py::python_typename<double>(),        c2py::python_typename<double>(),        c2py::python_typename<int>(),
                              c2py::python_typename<bool>(),          c2py::python_typename<std::string>(),   c2py::python_typename<unsigned long>(),
                              c2py::python_typename<bool>(),          c2py::python_typename<unsigned long>(), c2py::python_typename<bool>(),
                              c2py::python_typename<bool>(),          c2py::python_typename<unsigned long>(), c2py::python_typename<bool>(),
                              c2py::python_typename<double>(),        c2py::python_typename<unsigned long>(), c2py::python_typename<bool>(),
                              c2py::python_typename<bool>(),          c2py::python_typename<double>(),        c2py::python_typename<double>(),
                              c2py::python_typename<unsigned long>(), c2py::python_typename<double>(),        c2py::python_typename<double>(),
                              c2py::python_typename<double>(),        c2py::python_typename<double>(),        c2py::python_typename<double>(),
                              c2py::python_typename<bool>(),          c2py::python_typename<bool>(),          c2py::python_typename<bool>(),
                              c2py::python_typename<double>(),        c2py::python_typename<unsigned long>(), c2py::python_typename<unsigned long>(),
                              c2py::python_typename<unsigned long>(), c2py::python_typename<unsigned long>(), c2py::python_typename<unsigned long>(),
                              c2py::python_typename<bool>(),          c2py::python_typename<std::string>(),   c2py::python_typename<bool>(),
                              c2py::python_typename<bool>(),          c2py::python_typename<bool>(),          c2py::python_typename<bool>(),
                              c2py::python_typename<bool>(),          c2py::python_typename<bool>(),          c2py::python_typename<bool>(),
                              c2py::python_typename<bool>(),          c2py::python_typename<unsigned long>(), c2py::python_typename<std::string>(),
                              c2py::python_typename<int>(),           c2py::python_typename<bool>(),          c2py::python_typename<bool>(),
                              c2py::python_typename<bool>(),          c2py::python_typename<bool>(),          c2py::python_typename<bool>()});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_2>[] = {

   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_30  = R"DOC(Perform a DMNRG (density-matrix NRG) calculation.)DOC";
constexpr auto _c2py_doc_member_31  = R"DOC(Perform a CFS (complete Fock space) calculation.)DOC";
constexpr auto _c2py_doc_member_32  = R"DOC(Perform an FDM (full-density-matrix) calculation.)DOC";
constexpr auto _c2py_doc_member_33  = R"DOC(Calculate expectation values using the FDM algorithm.)DOC";
constexpr auto _c2py_doc_member_34  = R"DOC(Perform the DMNRG calculation on the Matsubara axis.)DOC";
constexpr auto _c2py_doc_member_35  = R"DOC(Perform the FDM calculation on the Matsubara axis.)DOC";
constexpr auto _c2py_doc_member_36  = R"DOC(Number of Matsubara points to collect.)DOC";
constexpr auto _c2py_doc_member_37  = R"DOC(Conductance curves to compute.)DOC";
constexpr auto _c2py_doc_member_38  = R"DOC(:math:`I_1` curves to compute.)DOC";
constexpr auto _c2py_doc_member_39  = R"DOC(:math:`I_2` curves to compute.)DOC";
constexpr auto _c2py_doc_member_40  = R"DOC(Compute the 3-leg vertex on the Matsubara/Matsubara axis.)DOC";
constexpr auto _c2py_doc_member_41  = R"DOC(Number of sites in the star representation (:math:`-1`: automatically determined).)DOC";
constexpr auto _c2py_doc_member_42  = R"DOC(Number of sites in the Wilson chain (:math:`-1`: automatically determined).)DOC";
constexpr auto _c2py_doc_member_43  = R"DOC(Largest :math:`x` in the discretization ODE solver (:math:`-1`: automatically determined).)DOC";
constexpr auto _c2py_doc_member_44  = R"DOC(Discretization scheme.)DOC";
constexpr auto _c2py_doc_member_45  = R"DOC(Parameter :math:`z` (twist) in the logarithmic discretization.)DOC";
constexpr auto _c2py_doc_member_46  = R"DOC(Tridiagonalisation approach.)DOC";
constexpr auto _c2py_doc_member_47  = R"DOC(Precision for tridiagonalisation.)DOC";
constexpr auto _c2py_doc_member_48  = R"DOC(Eigensolver routine (dsyev|dsyevr|zheev|zheevr|default).)DOC";
constexpr auto _c2py_doc_member_49  = R"DOC(Ratio of eigenstates computed in partial diagonalisation.)DOC";
constexpr auto _c2py_doc_member_50  = R"DOC(Minimal matrix size for dsyevr.)DOC";
constexpr auto _c2py_doc_member_51  = R"DOC(Minimal matrix size for zheevr.)DOC";
constexpr auto _c2py_doc_member_52  = R"DOC(Restart the calculation to achieve the truncation goal.)DOC";
constexpr auto _c2py_doc_member_53  = R"DOC(Rescale factor used when restart is true.)DOC";
constexpr auto _c2py_doc_member_54  = R"DOC(Additional states to keep in case of a near degeneracy.)DOC";
constexpr auto _c2py_doc_member_55  = R"DOC(Maximal number of additional states to keep.)DOC";
constexpr auto _c2py_doc_member_56  = R"DOC(Threshold value for eigenvalue splitting corrections.)DOC";
constexpr auto _c2py_doc_member_57  = R"DOC(Parameter :math:`\bar{\beta}` for thermodynamics.)DOC";
constexpr auto _c2py_doc_member_58  = R"DOC(Parameter :math:`p` for :math:`G(T)` calculations.)DOC";
constexpr auto _c2py_doc_member_59  = R"DOC(Parameter :math:`p` for :math:`\chi(T)` calculations.)DOC";
constexpr auto _c2py_doc_member_60  = R"DOC(Perform a Costi-Hewson-Zlatic finite-T calculation.)DOC";
constexpr auto _c2py_doc_member_61  = R"DOC(Compute the CFS greater correlation function.)DOC";
constexpr auto _c2py_doc_member_62  = R"DOC(Compute the CFS lesser correlation function.)DOC";
constexpr auto _c2py_doc_member_63  = R"DOC(Compute the FDM greater correlation function.)DOC";
constexpr auto _c2py_doc_member_64  = R"DOC(Compute the FDM lesser correlation function.)DOC";
constexpr auto _c2py_doc_member_65  = R"DOC(Iteration at which the expectation values are evaluated.)DOC";
constexpr auto _c2py_doc_member_66  = R"DOC(Perform a :math:`T > 0` calculation on the Matsubara axis.)DOC";
constexpr auto _c2py_doc_member_67  = R"DOC(Compute density matrices.)DOC";
constexpr auto _c2py_doc_member_68  = R"DOC(Auto-tune the `broaden_min` parameter.)DOC";
constexpr auto _c2py_doc_member_69  = R"DOC(Smallest energy scale in the problem, :math:`\omega_0`.)DOC";
constexpr auto _c2py_doc_member_70  = R"DOC(Sets :math:`\omega_0 = \mathtt{omega0\_ratio} \times T`.)DOC";
constexpr auto _c2py_doc_member_71  = R"DOC(Number of diagonalisation threads.)DOC";
constexpr auto _c2py_doc_member_72  = R"DOC(Use the interleaved diagonalization scheme.)DOC";
constexpr auto _c2py_doc_member_73  = R"DOC(Recalculation strategy.)DOC";
constexpr auto _c2py_doc_member_74  = R"DOC(Number of initial Wilson chain operators.)DOC";
constexpr auto _c2py_doc_member_75  = R"DOC(Output the imaginary parts of the correlators.)DOC";
constexpr auto _c2py_doc_member_76  = R"DOC(Number of eigenvalues to dump.)DOC";
constexpr auto _c2py_doc_member_77  = R"DOC(Dump in terms of absolute energies.)DOC";
constexpr auto _c2py_doc_member_78  = R"DOC(Dump using omega_N energy units.)DOC";
constexpr auto _c2py_doc_member_79  = R"DOC(Number of digits of precision used when dumping.)DOC";
constexpr auto _c2py_doc_member_80  = R"DOC(Dump by grouping degenerate states.)DOC";
constexpr auto _c2py_doc_member_81  = R"DOC(Energy tolerance for considering two states as degenerate.)DOC";
constexpr auto _c2py_doc_member_82  = R"DOC(Dump diagonal matrix elements.)DOC";
constexpr auto _c2py_doc_member_83  = R"DOC(Save binned (unbroadened) data.)DOC";
constexpr auto _c2py_doc_member_84  = R"DOC(Enable broadening of spectra.)DOC";
constexpr auto _c2py_doc_member_85  = R"DOC(Lower binning limit.)DOC";
constexpr auto _c2py_doc_member_86  = R"DOC(Upper binning limit.)DOC";
constexpr auto _c2py_doc_member_87  = R"DOC(Number of bins per decade for spectral data.)DOC";
constexpr auto _c2py_doc_member_88  = R"DOC(Shift of the accumulation points for binning.)DOC";
constexpr auto _c2py_doc_member_89  = R"DOC(Bin width for the linear mesh.)DOC";
constexpr auto _c2py_doc_member_90  = R"DOC(Peak clipping at the end of the run.)DOC";
constexpr auto _c2py_doc_member_91  = R"DOC(Peak clipping on the fly.)DOC";
constexpr auto _c2py_doc_member_92  = R"DOC(Energy window parameter for patching.)DOC";
constexpr auto _c2py_doc_member_93  = R"DOC(Perform N/N+1 patching.)DOC";
constexpr auto _c2py_doc_member_94  = R"DOC(Use even iterations in N/N+2 patching.)DOC";
constexpr auto _c2py_doc_member_95  = R"DOC(Average over even and odd N/N+2 spectra.)DOC";
constexpr auto _c2py_doc_member_96  = R"DOC(Parameter :math:`a` in the :math:`\tanh[a(x-0.5)]` window function.)DOC";
constexpr auto _c2py_doc_member_97  = R"DOC(Width of columns in the 'td' output file.)DOC";
constexpr auto _c2py_doc_member_98  = R"DOC(Width of columns in the 'custom' output file.)DOC";
constexpr auto _c2py_doc_member_99  = R"DOC(Precision of columns in the 'td' output file.)DOC";
constexpr auto _c2py_doc_member_100 = R"DOC(Precision of columns in the 'custom' output file.)DOC";
constexpr auto _c2py_doc_member_101 = R"DOC(Precision of the spectral function output.)DOC";
constexpr auto _c2py_doc_member_102 = R"DOC(Attempt to restart the calculation.)DOC";
constexpr auto _c2py_doc_member_103 = R"DOC(List of tokens defining what to log.)DOC";
constexpr auto _c2py_doc_member_104 = R"DOC(Log everything.)DOC";
constexpr auto _c2py_doc_member_105 = R"DOC(Create a DONE file.)DOC";
constexpr auto _c2py_doc_member_106 = R"DOC(Perform calculations at the 0-th iteration.)DOC";
constexpr auto _c2py_doc_member_107 = R"DOC(Keep all states in the last iteration for DMNRG.)DOC";
constexpr auto _c2py_doc_member_108 = R"DOC(Override the automatic lastall setting.)DOC";
constexpr auto _c2py_doc_member_109 = R"DOC(Save detailed subspace info.)DOC";
constexpr auto _c2py_doc_member_110 = R"DOC(Dump :math:`\langle f \rangle` matrix elements.)DOC";
constexpr auto _c2py_doc_member_111 = R"DOC(Dump all energies to a file.)DOC";
constexpr auto _c2py_doc_member_112 = R"DOC(Number of eigenvalues to show for log=e.)DOC";
constexpr auto _c2py_doc_member_113 = R"DOC(Stop the calculation at a given point.)DOC";
constexpr auto _c2py_doc_member_114 = R"DOC(Force stop at the given iteration (-1: disabled).)DOC";
constexpr auto _c2py_doc_member_115 = R"DOC(Remove temporary data files.)DOC";
constexpr auto _c2py_doc_member_116 = R"DOC(Do not output the imaginary parts of expectation values.)DOC";
constexpr auto _c2py_doc_member_117 = R"DOC(Check operator sum rules.)DOC";
constexpr auto _c2py_doc_member_118 = R"DOC(Test the diagonalisation results.)DOC";
constexpr auto _c2py_doc_member_119 = R"DOC(Test that :math:`\mathrm{tr}(\rho) = 1`.)DOC";
static PyObject *prop_get_dict_2(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_2> *)self)->_c);
  c2py::pydict dic;
  dic["dmnrg"]               = self_c.dmnrg;
  dic["cfs"]                 = self_c.cfs;
  dic["fdm"]                 = self_c.fdm;
  dic["fdmexpv"]             = self_c.fdmexpv;
  dic["dmnrgmats"]           = self_c.dmnrgmats;
  dic["fdmmats"]             = self_c.fdmmats;
  dic["mats"]                = self_c.mats;
  dic["specgt"]              = self_c.specgt;
  dic["speci1t"]             = self_c.speci1t;
  dic["speci2t"]             = self_c.speci2t;
  dic["v3mm"]                = self_c.v3mm;
  dic["mMAX"]                = self_c.mMAX;
  dic["Nmax"]                = self_c.Nmax;
  dic["xmax"]                = self_c.xmax;
  dic["discretization"]      = self_c.discretization;
  dic["z"]                   = self_c.z;
  dic["tri"]                 = self_c.tri;
  dic["preccpp"]             = self_c.preccpp;
  dic["diag"]                = self_c.diag;
  dic["diagratio"]           = self_c.diagratio;
  dic["dsyevrlimit"]         = self_c.dsyevrlimit;
  dic["zheevrlimit"]         = self_c.zheevrlimit;
  dic["restart"]             = self_c.restart;
  dic["restartfactor"]       = self_c.restartfactor;
  dic["safeguard"]           = self_c.safeguard;
  dic["safeguardmax"]        = self_c.safeguardmax;
  dic["fixeps"]              = self_c.fixeps;
  dic["betabar"]             = self_c.betabar;
  dic["gtp"]                 = self_c.gtp;
  dic["chitp"]               = self_c.chitp;
  dic["finite"]              = self_c.finite;
  dic["cfsgt"]               = self_c.cfsgt;
  dic["cfsls"]               = self_c.cfsls;
  dic["fdmgt"]               = self_c.fdmgt;
  dic["fdmls"]               = self_c.fdmls;
  dic["fdmexpvn"]            = self_c.fdmexpvn;
  dic["finitemats"]          = self_c.finitemats;
  dic["dm"]                  = self_c.dm;
  dic["broaden_min_ratio"]   = self_c.broaden_min_ratio;
  dic["omega0"]              = self_c.omega0;
  dic["omega0_ratio"]        = self_c.omega0_ratio;
  dic["diagth"]              = self_c.diagth;
  dic["substeps"]            = self_c.substeps;
  dic["strategy"]            = self_c.strategy;
  dic["Ninit"]               = self_c.Ninit;
  dic["reim"]                = self_c.reim;
  dic["dumpannotated"]       = self_c.dumpannotated;
  dic["dumpabs"]             = self_c.dumpabs;
  dic["dumpscaled"]          = self_c.dumpscaled;
  dic["dumpprecision"]       = self_c.dumpprecision;
  dic["dumpgroups"]          = self_c.dumpgroups;
  dic["grouptol"]            = self_c.grouptol;
  dic["dumpdiagonal"]        = self_c.dumpdiagonal;
  dic["savebins"]            = self_c.savebins;
  dic["broaden"]             = self_c.broaden;
  dic["emin"]                = self_c.emin;
  dic["emax"]                = self_c.emax;
  dic["bins"]                = self_c.bins;
  dic["accumulation"]        = self_c.accumulation;
  dic["linstep"]             = self_c.linstep;
  dic["discard_trim"]        = self_c.discard_trim;
  dic["discard_immediately"] = self_c.discard_immediately;
  dic["goodE"]               = self_c.goodE;
  dic["NN1"]                 = self_c.NN1;
  dic["NN2even"]             = self_c.NN2even;
  dic["NN2avg"]              = self_c.NN2avg;
  dic["NNtanh"]              = self_c.NNtanh;
  dic["width_td"]            = self_c.width_td;
  dic["width_custom"]        = self_c.width_custom;
  dic["prec_td"]             = self_c.prec_td;
  dic["prec_custom"]         = self_c.prec_custom;
  dic["prec_xy"]             = self_c.prec_xy;
  dic["resume"]              = self_c.resume;
  dic["log"]                 = self_c.log;
  dic["logall"]              = self_c.logall;
  dic["done"]                = self_c.done;
  dic["calc0"]               = self_c.calc0;
  dic["lastall"]             = self_c.lastall;
  dic["lastalloverride"]     = self_c.lastalloverride;
  dic["dumpsubspaces"]       = self_c.dumpsubspaces;
  dic["dump_f"]              = self_c.dump_f;
  dic["dumpenergies"]        = self_c.dumpenergies;
  dic["logenumber"]          = self_c.logenumber;
  dic["stopafter"]           = self_c.stopafter;
  dic["forcestop"]           = self_c.forcestop;
  dic["removefiles"]         = self_c.removefiles;
  dic["noimag"]              = self_c.noimag;
  dic["checksumrules"]       = self_c.checksumrules;
  dic["checkdiag"]           = self_c.checkdiag;
  dic["checkrho"]            = self_c.checkrho;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_2>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_2::dmnrg, _c2py_cls_2>("dmnrg", _c2py_doc_member_30),
   c2py::getsetdef_from_member<&_c2py_cls_2::cfs, _c2py_cls_2>("cfs", _c2py_doc_member_31),
   c2py::getsetdef_from_member<&_c2py_cls_2::fdm, _c2py_cls_2>("fdm", _c2py_doc_member_32),
   c2py::getsetdef_from_member<&_c2py_cls_2::fdmexpv, _c2py_cls_2>("fdmexpv", _c2py_doc_member_33),
   c2py::getsetdef_from_member<&_c2py_cls_2::dmnrgmats, _c2py_cls_2>("dmnrgmats", _c2py_doc_member_34),
   c2py::getsetdef_from_member<&_c2py_cls_2::fdmmats, _c2py_cls_2>("fdmmats", _c2py_doc_member_35),
   c2py::getsetdef_from_member<&_c2py_cls_2::mats, _c2py_cls_2>("mats", _c2py_doc_member_36),
   c2py::getsetdef_from_member<&_c2py_cls_2::specgt, _c2py_cls_2>("specgt", _c2py_doc_member_37),
   c2py::getsetdef_from_member<&_c2py_cls_2::speci1t, _c2py_cls_2>("speci1t", _c2py_doc_member_38),
   c2py::getsetdef_from_member<&_c2py_cls_2::speci2t, _c2py_cls_2>("speci2t", _c2py_doc_member_39),
   c2py::getsetdef_from_member<&_c2py_cls_2::v3mm, _c2py_cls_2>("v3mm", _c2py_doc_member_40),
   c2py::getsetdef_from_member<&_c2py_cls_2::mMAX, _c2py_cls_2>("mMAX", _c2py_doc_member_41),
   c2py::getsetdef_from_member<&_c2py_cls_2::Nmax, _c2py_cls_2>("Nmax", _c2py_doc_member_42),
   c2py::getsetdef_from_member<&_c2py_cls_2::xmax, _c2py_cls_2>("xmax", _c2py_doc_member_43),
   c2py::getsetdef_from_member<&_c2py_cls_2::discretization, _c2py_cls_2>("discretization", _c2py_doc_member_44),
   c2py::getsetdef_from_member<&_c2py_cls_2::z, _c2py_cls_2>("z", _c2py_doc_member_45),
   c2py::getsetdef_from_member<&_c2py_cls_2::tri, _c2py_cls_2>("tri", _c2py_doc_member_46),
   c2py::getsetdef_from_member<&_c2py_cls_2::preccpp, _c2py_cls_2>("preccpp", _c2py_doc_member_47),
   c2py::getsetdef_from_member<&_c2py_cls_2::diag, _c2py_cls_2>("diag", _c2py_doc_member_48),
   c2py::getsetdef_from_member<&_c2py_cls_2::diagratio, _c2py_cls_2>("diagratio", _c2py_doc_member_49),
   c2py::getsetdef_from_member<&_c2py_cls_2::dsyevrlimit, _c2py_cls_2>("dsyevrlimit", _c2py_doc_member_50),
   c2py::getsetdef_from_member<&_c2py_cls_2::zheevrlimit, _c2py_cls_2>("zheevrlimit", _c2py_doc_member_51),
   c2py::getsetdef_from_member<&_c2py_cls_2::restart, _c2py_cls_2>("restart", _c2py_doc_member_52),
   c2py::getsetdef_from_member<&_c2py_cls_2::restartfactor, _c2py_cls_2>("restartfactor", _c2py_doc_member_53),
   c2py::getsetdef_from_member<&_c2py_cls_2::safeguard, _c2py_cls_2>("safeguard", _c2py_doc_member_54),
   c2py::getsetdef_from_member<&_c2py_cls_2::safeguardmax, _c2py_cls_2>("safeguardmax", _c2py_doc_member_55),
   c2py::getsetdef_from_member<&_c2py_cls_2::fixeps, _c2py_cls_2>("fixeps", _c2py_doc_member_56),
   c2py::getsetdef_from_member<&_c2py_cls_2::betabar, _c2py_cls_2>("betabar", _c2py_doc_member_57),
   c2py::getsetdef_from_member<&_c2py_cls_2::gtp, _c2py_cls_2>("gtp", _c2py_doc_member_58),
   c2py::getsetdef_from_member<&_c2py_cls_2::chitp, _c2py_cls_2>("chitp", _c2py_doc_member_59),
   c2py::getsetdef_from_member<&_c2py_cls_2::finite, _c2py_cls_2>("finite", _c2py_doc_member_60),
   c2py::getsetdef_from_member<&_c2py_cls_2::cfsgt, _c2py_cls_2>("cfsgt", _c2py_doc_member_61),
   c2py::getsetdef_from_member<&_c2py_cls_2::cfsls, _c2py_cls_2>("cfsls", _c2py_doc_member_62),
   c2py::getsetdef_from_member<&_c2py_cls_2::fdmgt, _c2py_cls_2>("fdmgt", _c2py_doc_member_63),
   c2py::getsetdef_from_member<&_c2py_cls_2::fdmls, _c2py_cls_2>("fdmls", _c2py_doc_member_64),
   c2py::getsetdef_from_member<&_c2py_cls_2::fdmexpvn, _c2py_cls_2>("fdmexpvn", _c2py_doc_member_65),
   c2py::getsetdef_from_member<&_c2py_cls_2::finitemats, _c2py_cls_2>("finitemats", _c2py_doc_member_66),
   c2py::getsetdef_from_member<&_c2py_cls_2::dm, _c2py_cls_2>("dm", _c2py_doc_member_67),
   c2py::getsetdef_from_member<&_c2py_cls_2::broaden_min_ratio, _c2py_cls_2>("broaden_min_ratio", _c2py_doc_member_68),
   c2py::getsetdef_from_member<&_c2py_cls_2::omega0, _c2py_cls_2>("omega0", _c2py_doc_member_69),
   c2py::getsetdef_from_member<&_c2py_cls_2::omega0_ratio, _c2py_cls_2>("omega0_ratio", _c2py_doc_member_70),
   c2py::getsetdef_from_member<&_c2py_cls_2::diagth, _c2py_cls_2>("diagth", _c2py_doc_member_71),
   c2py::getsetdef_from_member<&_c2py_cls_2::substeps, _c2py_cls_2>("substeps", _c2py_doc_member_72),
   c2py::getsetdef_from_member<&_c2py_cls_2::strategy, _c2py_cls_2>("strategy", _c2py_doc_member_73),
   c2py::getsetdef_from_member<&_c2py_cls_2::Ninit, _c2py_cls_2>("Ninit", _c2py_doc_member_74),
   c2py::getsetdef_from_member<&_c2py_cls_2::reim, _c2py_cls_2>("reim", _c2py_doc_member_75),
   c2py::getsetdef_from_member<&_c2py_cls_2::dumpannotated, _c2py_cls_2>("dumpannotated", _c2py_doc_member_76),
   c2py::getsetdef_from_member<&_c2py_cls_2::dumpabs, _c2py_cls_2>("dumpabs", _c2py_doc_member_77),
   c2py::getsetdef_from_member<&_c2py_cls_2::dumpscaled, _c2py_cls_2>("dumpscaled", _c2py_doc_member_78),
   c2py::getsetdef_from_member<&_c2py_cls_2::dumpprecision, _c2py_cls_2>("dumpprecision", _c2py_doc_member_79),
   c2py::getsetdef_from_member<&_c2py_cls_2::dumpgroups, _c2py_cls_2>("dumpgroups", _c2py_doc_member_80),
   c2py::getsetdef_from_member<&_c2py_cls_2::grouptol, _c2py_cls_2>("grouptol", _c2py_doc_member_81),
   c2py::getsetdef_from_member<&_c2py_cls_2::dumpdiagonal, _c2py_cls_2>("dumpdiagonal", _c2py_doc_member_82),
   c2py::getsetdef_from_member<&_c2py_cls_2::savebins, _c2py_cls_2>("savebins", _c2py_doc_member_83),
   c2py::getsetdef_from_member<&_c2py_cls_2::broaden, _c2py_cls_2>("broaden", _c2py_doc_member_84),
   c2py::getsetdef_from_member<&_c2py_cls_2::emin, _c2py_cls_2>("emin", _c2py_doc_member_85),
   c2py::getsetdef_from_member<&_c2py_cls_2::emax, _c2py_cls_2>("emax", _c2py_doc_member_86),
   c2py::getsetdef_from_member<&_c2py_cls_2::bins, _c2py_cls_2>("bins", _c2py_doc_member_87),
   c2py::getsetdef_from_member<&_c2py_cls_2::accumulation, _c2py_cls_2>("accumulation", _c2py_doc_member_88),
   c2py::getsetdef_from_member<&_c2py_cls_2::linstep, _c2py_cls_2>("linstep", _c2py_doc_member_89),
   c2py::getsetdef_from_member<&_c2py_cls_2::discard_trim, _c2py_cls_2>("discard_trim", _c2py_doc_member_90),
   c2py::getsetdef_from_member<&_c2py_cls_2::discard_immediately, _c2py_cls_2>("discard_immediately", _c2py_doc_member_91),
   c2py::getsetdef_from_member<&_c2py_cls_2::goodE, _c2py_cls_2>("goodE", _c2py_doc_member_92),
   c2py::getsetdef_from_member<&_c2py_cls_2::NN1, _c2py_cls_2>("NN1", _c2py_doc_member_93),
   c2py::getsetdef_from_member<&_c2py_cls_2::NN2even, _c2py_cls_2>("NN2even", _c2py_doc_member_94),
   c2py::getsetdef_from_member<&_c2py_cls_2::NN2avg, _c2py_cls_2>("NN2avg", _c2py_doc_member_95),
   c2py::getsetdef_from_member<&_c2py_cls_2::NNtanh, _c2py_cls_2>("NNtanh", _c2py_doc_member_96),
   c2py::getsetdef_from_member<&_c2py_cls_2::width_td, _c2py_cls_2>("width_td", _c2py_doc_member_97),
   c2py::getsetdef_from_member<&_c2py_cls_2::width_custom, _c2py_cls_2>("width_custom", _c2py_doc_member_98),
   c2py::getsetdef_from_member<&_c2py_cls_2::prec_td, _c2py_cls_2>("prec_td", _c2py_doc_member_99),
   c2py::getsetdef_from_member<&_c2py_cls_2::prec_custom, _c2py_cls_2>("prec_custom", _c2py_doc_member_100),
   c2py::getsetdef_from_member<&_c2py_cls_2::prec_xy, _c2py_cls_2>("prec_xy", _c2py_doc_member_101),
   c2py::getsetdef_from_member<&_c2py_cls_2::resume, _c2py_cls_2>("resume", _c2py_doc_member_102),
   c2py::getsetdef_from_member<&_c2py_cls_2::log, _c2py_cls_2>("log", _c2py_doc_member_103),
   c2py::getsetdef_from_member<&_c2py_cls_2::logall, _c2py_cls_2>("logall", _c2py_doc_member_104),
   c2py::getsetdef_from_member<&_c2py_cls_2::done, _c2py_cls_2>("done", _c2py_doc_member_105),
   c2py::getsetdef_from_member<&_c2py_cls_2::calc0, _c2py_cls_2>("calc0", _c2py_doc_member_106),
   c2py::getsetdef_from_member<&_c2py_cls_2::lastall, _c2py_cls_2>("lastall", _c2py_doc_member_107),
   c2py::getsetdef_from_member<&_c2py_cls_2::lastalloverride, _c2py_cls_2>("lastalloverride", _c2py_doc_member_108),
   c2py::getsetdef_from_member<&_c2py_cls_2::dumpsubspaces, _c2py_cls_2>("dumpsubspaces", _c2py_doc_member_109),
   c2py::getsetdef_from_member<&_c2py_cls_2::dump_f, _c2py_cls_2>("dump_f", _c2py_doc_member_110),
   c2py::getsetdef_from_member<&_c2py_cls_2::dumpenergies, _c2py_cls_2>("dumpenergies", _c2py_doc_member_111),
   c2py::getsetdef_from_member<&_c2py_cls_2::logenumber, _c2py_cls_2>("logenumber", _c2py_doc_member_112),
   c2py::getsetdef_from_member<&_c2py_cls_2::stopafter, _c2py_cls_2>("stopafter", _c2py_doc_member_113),
   c2py::getsetdef_from_member<&_c2py_cls_2::forcestop, _c2py_cls_2>("forcestop", _c2py_doc_member_114),
   c2py::getsetdef_from_member<&_c2py_cls_2::removefiles, _c2py_cls_2>("removefiles", _c2py_doc_member_115),
   c2py::getsetdef_from_member<&_c2py_cls_2::noimag, _c2py_cls_2>("noimag", _c2py_doc_member_116),
   c2py::getsetdef_from_member<&_c2py_cls_2::checksumrules, _c2py_cls_2>("checksumrules", _c2py_doc_member_117),
   c2py::getsetdef_from_member<&_c2py_cls_2::checkdiag, _c2py_cls_2>("checkdiag", _c2py_doc_member_118),
   c2py::getsetdef_from_member<&_c2py_cls_2::checkrho, _c2py_cls_2>("checkrho", _c2py_doc_member_119),
   {"__dict__", (getter)prop_get_dict_2, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_2> =
   R"DOC(Low-level NRG parameters.)DOC" + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_2>;
// --------- class _c2py_cls_3 -----------
using _c2py_cls_3                                            = nrgljubljana_interface::solver_core;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_3>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_3> = "nrgljubljana_interface.solver_core.SolverCore";
static const auto _c2py_init_0 = c2py::dispatcher_c_kw_t{c2py::c_constructor<_c2py_cls_3, nrgljubljana_interface::constr_params_t>("cp")};
template <> constexpr initproc c2py::tp_init<_c2py_cls_3> = c2py::pyfkw_constructor<_c2py_init_0>;
template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_3> = _c2py_init_0.doc(R"DOC(
Construct an NRGLjubljana_interface solver.

Parameters
----------
cp : {par_0}
   Construction parameters.
)DOC",
                                                                    {{c2py::python_typename<nrgljubljana_interface::constr_params_t>()}});
// be_quiet
static auto const _c2py_fun_1 = c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_3 &self) -> decltype(auto) { return self.be_quiet(); }, "self")};

// check_model_params
static auto const _c2py_fun_2 = c2py::dispatcher_f_kw_t{c2py::cmethod(
   [](_c2py_cls_3 &self, const nrgljubljana_interface::solve_params_t &sp) -> decltype(auto) { return self.check_model_params(sp); }, "self", "sp")};

// create_tempdir
static auto const _c2py_fun_3 = c2py::dispatcher_f_kw_t{
   c2py::cmethod([](_c2py_cls_3 &self, const std::string &tempdir_) -> decltype(auto) { return self.create_tempdir(tempdir_); }, "self", "tempdir_")};

// generate_param_file
static auto const _c2py_fun_4 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_3 &self, double z) -> decltype(auto) { return self.generate_param_file(z); }, "self", "z")};

// instantiate
static auto const _c2py_fun_5 = c2py::dispatcher_f_kw_t{c2py::cmethod(
   [](_c2py_cls_3 &self, double z, const std::string &taskdir) -> decltype(auto) { return self.instantiate(z, taskdir); }, "self", "z", "taskdir")};

// read_structure
static auto const _c2py_fun_6 = c2py::dispatcher_f_kw_t{c2py::cmethod(
   [](_c2py_cls_3 &self, const std::string &filename, bool mandatory) -> decltype(auto) { return self.read_structure(filename, mandatory); }, "self",
   "filename", "mandatory")};

// readexpv
static auto const _c2py_fun_7 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_3 &self, int Nz) -> decltype(auto) { return self.readexpv(Nz); }, "self", "Nz")};

// readtdfdm
static auto const _c2py_fun_8 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_3 &self, int Nz) -> decltype(auto) { return self.readtdfdm(Nz); }, "self", "Nz")};

// set_nrg_params
static auto const _c2py_fun_9 = c2py::dispatcher_f_kw_t{c2py::cmethod(
   [](_c2py_cls_3 &self, const nrgljubljana_interface::nrg_params_t &nrg_params) -> decltype(auto) { return self.set_nrg_params(nrg_params); },
   "self", "nrg_params")};

// set_params
static auto const _c2py_fun_10 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_3 &self) -> decltype(auto) { return self.set_params(); }, "self")};

// set_verbosity
static auto const _c2py_fun_11 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_3 &self, bool v) -> decltype(auto) { return self.set_verbosity(v); }, "self", "v")};

// solve
static auto const _c2py_fun_12 = c2py::dispatcher_f_kw_t{c2py::cmethod(
   [](_c2py_cls_3 &self, const nrgljubljana_interface::solve_params_t &solve_params) -> decltype(auto) { return self.solve(solve_params); }, "self",
   "solve_params")};

// solve_one
static auto const _c2py_fun_13 = c2py::dispatcher_f_kw_t{
   c2py::cmethod([](_c2py_cls_3 &self, const std::string &taskdir) -> decltype(auto) { return self.solve_one(taskdir); }, "self", "taskdir")};

// write_gamma
static auto const _c2py_fun_14 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_3 &self) -> decltype(auto) { return self.write_gamma(); }, "self")};

static const auto _c2py_doc_1  = _c2py_fun_1.doc(R"DOC(
Suppress verbose output from the NRG solver.
)DOC");
static const auto _c2py_doc_2  = _c2py_fun_2.doc(R"DOC(
Check that all required model parameters have been defined.

Parameters
----------
sp : {par_0}
   Solve parameters to validate.
)DOC",
                                                 {{c2py::python_typename<const nrgljubljana_interface::solve_params_t &>()}});
static const auto _c2py_doc_3  = _c2py_fun_3.doc(R"DOC(
Create a temporary working directory for a series of NRG runs.

Parameters
----------
tempdir_ : {par_0}
   Base directory in which to create the temporary directory.

Returns
-------
{ret_0}
   Path to the created temporary directory.
)DOC",
                                                 {{c2py::python_typename<const std::string &>()}}, {c2py::python_typename<std::string>()});
static const auto _c2py_doc_4  = _c2py_fun_4.doc(R"DOC(
Produce the param file for a given value of the twist parameter :math:`z`.

Parameters
----------
z : {par_0}
   Discretization twist parameter.
)DOC",
                                                 {{c2py::python_typename<double>()}});
static const auto _c2py_doc_5  = _c2py_fun_5.doc(R"DOC(
Prepare the input files for an individual NRG calculation. Called from solve().

Parameters
----------
z : {par_0}
   Discretization twist parameter.
taskdir : {par_1}
   Directory in which to prepare the calculation.
)DOC",
                                                 {{c2py::python_typename<double>()}, {c2py::python_typename<const std::string &>()}});
static const auto _c2py_doc_6  = _c2py_fun_6.doc(R"DOC(
Read the block structure of Green's function objects from a file.

Parameters
----------
filename : {par_0}
   File describing the block (gf) structure.
mandatory : {par_1}
   If true, a missing file is an error.

Returns
-------
{ret_0}
   The Green's function block structure.
)DOC",
                                                 {{c2py::python_typename<const std::string &>()}, {c2py::python_typename<bool>()}},
                                                 {c2py::python_typename<triqs::gfs::gf_struct_t>()});
static const auto _c2py_doc_7  = _c2py_fun_7.doc(R"DOC(
Read expectation values from the NRG output files.

Parameters
----------
Nz : {par_0}
   Number of discretization twists.
)DOC",
                                                 {{c2py::python_typename<int>()}});
static const auto _c2py_doc_8  = _c2py_fun_8.doc(R"DOC(
Read thermodynamic variables (FDM algorithm) from the NRG output files.

Parameters
----------
Nz : {par_0}
   Number of discretization twists.
)DOC",
                                                 {{c2py::python_typename<int>()}});
static const auto _c2py_doc_9  = _c2py_fun_9.doc(R"DOC(
Adjust the advanced (low-level) NRG parameters.

Parameters
----------
nrg_params : {par_0}
   Low-level NRG parameters.
)DOC",
                                                 {{c2py::python_typename<const nrgljubljana_interface::nrg_params_t &>()}});
static const auto _c2py_doc_10 = _c2py_fun_10.doc(R"DOC(
Establish good defaults for the low-level NRG parameters.
)DOC");
static const auto _c2py_doc_11 = _c2py_fun_11.doc(R"DOC(
Set the verbosity (see also ``be_quiet()``).

Parameters
----------
v : {par_0}
   If true, enable verbose output.
)DOC",
                                                  {{c2py::python_typename<bool>()}});
static const auto _c2py_doc_12 = _c2py_fun_12.doc(R"DOC(
Solve the impurity problem.

Runs the full NRGLjubljana pipeline for the assigned hybridization function and
populates the result containers.

Parameters
----------
solve_params : {par_0}
   Parameters specific to the NRGLjubljana run.
)DOC",
                                                  {{c2py::python_typename<const nrgljubljana_interface::solve_params_t &>()}});
static const auto _c2py_doc_13 = _c2py_fun_13.doc(R"DOC(
Perform an individual NRG calculation. Called from ``solve()``.

Parameters
----------
taskdir : {par_0}
   Directory containing the prepared calculation.
)DOC",
                                                  {{c2py::python_typename<const std::string &>()}});
static const auto _c2py_doc_14 = _c2py_fun_14.doc(R"DOC(
Write :math:`\Gamma = -\mathrm{Im}\,\Delta(\omega)` to a file.
)DOC");

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_3>[] = {
   {"be_quiet", (PyCFunction)c2py::pyfkw<_c2py_fun_1>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_1.c_str()},
   {"check_model_params", (PyCFunction)c2py::pyfkw<_c2py_fun_2>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_2.c_str()},
   {"create_tempdir", (PyCFunction)c2py::pyfkw<_c2py_fun_3>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_3.c_str()},
   {"generate_param_file", (PyCFunction)c2py::pyfkw<_c2py_fun_4>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_4.c_str()},
   {"instantiate", (PyCFunction)c2py::pyfkw<_c2py_fun_5>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_5.c_str()},
   {"read_structure", (PyCFunction)c2py::pyfkw<_c2py_fun_6>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_6.c_str()},
   {"readexpv", (PyCFunction)c2py::pyfkw<_c2py_fun_7>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_7.c_str()},
   {"readtdfdm", (PyCFunction)c2py::pyfkw<_c2py_fun_8>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_8.c_str()},
   {"set_nrg_params", (PyCFunction)c2py::pyfkw<_c2py_fun_9>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_9.c_str()},
   {"set_params", (PyCFunction)c2py::pyfkw<_c2py_fun_10>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_10.c_str()},
   {"set_verbosity", (PyCFunction)c2py::pyfkw<_c2py_fun_11>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_11.c_str()},
   {"solve", (PyCFunction)c2py::pyfkw<_c2py_fun_12>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_12.c_str()},
   {"solve_one", (PyCFunction)c2py::pyfkw<_c2py_fun_13>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_13.c_str()},
   {"write_gamma", (PyCFunction)c2py::pyfkw<_c2py_fun_14>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_14.c_str()},
   {"__write_hdf5__", c2py::tpxx_write_h5<_c2py_cls_3>, METH_VARARGS, "  "},
   {"__getstate__", c2py::getstate_h5<_c2py_cls_3>, METH_NOARGS, ""},
   {"__setstate__", c2py::setstate_h5<_c2py_cls_3>, METH_O, ""},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_120 = R"DOC(Parameters used for the solver construction.)DOC";
constexpr auto _c2py_doc_member_121 = R"DOC(Low-level NRG parameters.)DOC";
constexpr auto _c2py_doc_member_122 = R"DOC(Parameters used for the most recent solve process.)DOC";
constexpr auto _c2py_doc_member_123 = R"DOC(If true, detailed output from NRGLjubljana and its tools is sent to stdout.)DOC";
constexpr auto _c2py_doc_member_124 = R"DOC(Keep the temporary directories after the calculation.)DOC";
constexpr auto _c2py_doc_member_125 = R"DOC(The Green's function structure object.)DOC";
constexpr auto _c2py_doc_member_126 = R"DOC(The hybridization function structure object.)DOC";
constexpr auto _c2py_doc_member_127 = R"DOC(The susceptibility structure object.)DOC";
constexpr auto _c2py_doc_member_128 = R"DOC(Logarithmic real-frequency mesh.)DOC";
constexpr auto _c2py_doc_member_129 = R"DOC(The hybridization function on the real-frequency axis.)DOC";
constexpr auto _c2py_doc_member_130 = R"DOC(The spectral function :math:`A(\omega)`.)DOC";
constexpr auto _c2py_doc_member_131 = R"DOC(The spectral function :math:`B_l(\omega)` of the auxiliary correlator :math:`F_l(\omega)`.)DOC";
constexpr auto _c2py_doc_member_132 = R"DOC(The spectral function :math:`B_r(\omega)` of the auxiliary correlator :math:`F_r(\omega)`.)DOC";
constexpr auto _c2py_doc_member_133 = R"DOC(The spectral function :math:`C(\omega)` of the auxiliary correlator :math:`I(\omega)`.)DOC";
constexpr auto _c2py_doc_member_134 = R"DOC(The retarded Green's function :math:`G(\omega)`.)DOC";
constexpr auto _c2py_doc_member_135 = R"DOC(The auxiliary Green's function :math:`F_l(\omega) = \Sigma(\omega)\, G(\omega)`.)DOC";
constexpr auto _c2py_doc_member_136 = R"DOC(The auxiliary Green's function :math:`F_r(\omega) = G(\omega)\, \Sigma(\omega)`.)DOC";
constexpr auto _c2py_doc_member_137 = R"DOC(The auxiliary Green's function :math:`I(\omega)`.)DOC";
constexpr auto _c2py_doc_member_138 = R"DOC(Constant Hartree shift to the self-energy, stored as a Green's function.)DOC";
constexpr auto _c2py_doc_member_139 =
   R"DOC(The retarded self-energy :math:`\Sigma(\omega)` (computed from :math:`F_l`, :math:`F_r`, :math:`G` and :math:`I`).)DOC";
constexpr auto _c2py_doc_member_140 = R"DOC(Expectation values of local impurity operators.)DOC";
constexpr auto _c2py_doc_member_141 = R"DOC(Thermodynamic variables (FDM algorithm).)DOC";
constexpr auto _c2py_doc_member_142 = R"DOC(Charge susceptibility :math:`\chi_{NN}(\omega)`.)DOC";
constexpr auto _c2py_doc_member_143 = R"DOC(Spin susceptibility :math:`\chi_{SS}(\omega)`.)DOC";

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_3>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_3::constr_params, _c2py_cls_3>("constr_params", _c2py_doc_member_120),
   c2py::getsetdef_from_member<&_c2py_cls_3::nrg_params, _c2py_cls_3>("nrg_params", _c2py_doc_member_121),
   c2py::getsetdef_from_member<&_c2py_cls_3::last_solve_params, _c2py_cls_3>("last_solve_params", _c2py_doc_member_122),
   c2py::getsetdef_from_member<&_c2py_cls_3::verbose, _c2py_cls_3>("verbose", _c2py_doc_member_123),
   c2py::getsetdef_from_member<&_c2py_cls_3::keep_temp_dir, _c2py_cls_3>("keep_temp_dir", _c2py_doc_member_124),
   c2py::getsetdef_from_member<&_c2py_cls_3::gf_struct, _c2py_cls_3>("gf_struct", _c2py_doc_member_125),
   c2py::getsetdef_from_member<&_c2py_cls_3::Delta_struct, _c2py_cls_3>("Delta_struct", _c2py_doc_member_126),
   c2py::getsetdef_from_member<&_c2py_cls_3::chi_struct, _c2py_cls_3>("chi_struct", _c2py_doc_member_127),
   c2py::getsetdef_from_member<&_c2py_cls_3::log_mesh, _c2py_cls_3>("log_mesh", _c2py_doc_member_128),
   c2py::getsetdef_from_member<&_c2py_cls_3::Delta_w, _c2py_cls_3>("Delta_w", _c2py_doc_member_129),
   c2py::getsetdef_from_member<&_c2py_cls_3::A_w, _c2py_cls_3>("A_w", _c2py_doc_member_130),
   c2py::getsetdef_from_member<&_c2py_cls_3::B_l_w, _c2py_cls_3>("B_l_w", _c2py_doc_member_131),
   c2py::getsetdef_from_member<&_c2py_cls_3::B_r_w, _c2py_cls_3>("B_r_w", _c2py_doc_member_132),
   c2py::getsetdef_from_member<&_c2py_cls_3::C_w, _c2py_cls_3>("C_w", _c2py_doc_member_133),
   c2py::getsetdef_from_member<&_c2py_cls_3::G_w, _c2py_cls_3>("G_w", _c2py_doc_member_134),
   c2py::getsetdef_from_member<&_c2py_cls_3::F_l_w, _c2py_cls_3>("F_l_w", _c2py_doc_member_135),
   c2py::getsetdef_from_member<&_c2py_cls_3::F_r_w, _c2py_cls_3>("F_r_w", _c2py_doc_member_136),
   c2py::getsetdef_from_member<&_c2py_cls_3::I_w, _c2py_cls_3>("I_w", _c2py_doc_member_137),
   c2py::getsetdef_from_member<&_c2py_cls_3::SigmaHartree_w, _c2py_cls_3>("SigmaHartree_w", _c2py_doc_member_138),
   c2py::getsetdef_from_member<&_c2py_cls_3::Sigma_w, _c2py_cls_3>("Sigma_w", _c2py_doc_member_139),
   c2py::getsetdef_from_member<&_c2py_cls_3::expv, _c2py_cls_3>("expv", _c2py_doc_member_140),
   c2py::getsetdef_from_member<&_c2py_cls_3::tdfdm, _c2py_cls_3>("tdfdm", _c2py_doc_member_141),
   c2py::getsetdef_from_member<&_c2py_cls_3::chi_NN_w, _c2py_cls_3>("chi_NN_w", _c2py_doc_member_142),
   c2py::getsetdef_from_member<&_c2py_cls_3::chi_SS_w, _c2py_cls_3>("chi_SS_w", _c2py_doc_member_143),

   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_3> = R"DOC(TRIQS interface to the NRGLjubljana numerical renormalization group impurity solver.

Given a hybridization function, drives the external NRGLjubljana solver to compute
spectral functions, Green's functions, the self-energy, susceptibilities and
thermodynamic/expectation values for quantum impurity models.)DOC"
   + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_3>;

// ==================== module functions ====================

// hilbert_transform_elementwise
static auto const _c2py_fun_15 = c2py::dispatcher_f_kw_t{c2py::cfun(
   [](const nrgljubljana_interface::m_w_cvt &gf, std::complex<double> z) { return nrgljubljana_interface::hilbert_transform_elementwise(gf, z); },
   "gf", "z")};

// hilbert_transform_refreq
static auto const _c2py_fun_16 = c2py::dispatcher_f_kw_t{c2py::cfun(
   [](const nrgljubljana_interface::c_w_cvt &gf, std::complex<double> z) { return nrgljubljana_interface::hilbert_transform_refreq(gf, z); }, "gf",
   "z")};

static const auto _c2py_doc_15 = _c2py_fun_15.doc(
   R"DOC(
Hilbert transform for refreq objects (matrix, elementwise).

Parameters
----------
gf : {par_0}
   Matrix-valued real-frequency Green's function.
z : {par_1}
   Complex frequency at which to evaluate the transform.

Returns
-------
{ret_0}
   Value of the elementwise Hilbert transform at :math:`z`.
)DOC",
   {{c2py::python_typename<const nrgljubljana_interface::m_w_cvt &>()}, {c2py::python_typename<std::complex<double>>()}},
   {c2py::python_typename<
      nda::basic_array<std::complex<double>, 2, nda::C_layout, 'M', nda::heap_basic<nda::mem::mallocator<nda::mem::AddressSpace::Host>>>>()});
static const auto _c2py_doc_16 =
   _c2py_fun_16.doc(R"DOC(
Hilbert transform for refreq objects (scalar).

Parameters
----------
gf : {par_0}
   Scalar real-frequency Green's function.
z : {par_1}
   Complex frequency at which to evaluate the transform.

Returns
-------
{ret_0}
   Value of the Hilbert transform at :math:`z`.
)DOC",
                    {{c2py::python_typename<const nrgljubljana_interface::c_w_cvt &>()}, {c2py::python_typename<std::complex<double>>()}},
                    {c2py::python_typename<std::complex<double>>()});
//--------------------- module function table  -----------------------------

static PyMethodDef module_methods[] = {
   {"hilbert_transform_elementwise", (PyCFunction)c2py::pyfkw<_c2py_fun_15>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_15.c_str()},
   {"hilbert_transform_refreq", (PyCFunction)c2py::pyfkw<_c2py_fun_16>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_16.c_str()},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

//--------------------- module struct & init error definition ------------

//// module doc directly in the code or "" if not present...
/// Or mandatory ?
static struct PyModuleDef module_def = {
   PyModuleDef_HEAD_INIT,
   "solver_core",                                                                                         /* name of module */
   R"RAWDOC(TRIQS interface to the NRGLjubljana numerical renormalization group impurity solver.)RAWDOC", /* module documentation, may be NULL */
   -1, /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
   module_methods,
   NULL,
   NULL,
   NULL,
   NULL};

//--------------------- module init function -----------------------------

extern "C" __attribute__((visibility("default"))) PyObject *PyInit_solver_core() {

  if (not c2py::check_python_version("solver_core")) return NULL;

  // import numpy iff 'numpy/arrayobject.h' included
#ifdef Py_ARRAYOBJECT_H
  import_array();
#endif

  PyObject *m;

  if (PyType_Ready(&c2py::wrap_pytype<c2py::py_range>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_0>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_1>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_2>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_3>) < 0) return NULL;

  m = PyModule_Create(&module_def);
  if (m == NULL) return NULL;

  auto &conv_table = *c2py::conv_table_sptr.get();

  conv_table[std::type_index(typeid(c2py::py_range)).name()] = &c2py::wrap_pytype<c2py::py_range>;
#define _add_type(T, N) c2py::add_type_object_to_main<T>(N, m, conv_table)
  _add_type(_c2py_cls_0, "ConstrParamsT");
  _add_type(_c2py_cls_1, "SolveParamsT");
  _add_type(_c2py_cls_2, "NrgParamsT");
  _add_type(_c2py_cls_3, "SolverCore");
#undef _add_type

  c2py::pyref module = c2py::pyref::module("h5.formats");
  if (not module) return nullptr;
  c2py::pyref register_class = module.attr("register_class");

  register_h5_type<_c2py_cls_3>(register_class);

  return m;
}
#endif
// CLAIR_WRAP_GEN
