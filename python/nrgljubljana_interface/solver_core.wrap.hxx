#include <c2py/c2py.hpp>

#ifndef C2PY_HXX_DECLARATION_solver_core_GUARDS
#define C2PY_HXX_DECLARATION_solver_core_GUARDS
template <> constexpr bool c2py::is_wrapped<nrgljubljana_interface::constr_params_t>     = true;
template <> inline constexpr auto c2py::tp_name<nrgljubljana_interface::constr_params_t> = "nrgljubljana_interface.solver_core.ConstrParamsT";
template <> constexpr bool c2py::is_wrapped<nrgljubljana_interface::solve_params_t>      = true;
template <> inline constexpr auto c2py::tp_name<nrgljubljana_interface::solve_params_t>  = "nrgljubljana_interface.solver_core.SolveParamsT";
template <> constexpr bool c2py::is_wrapped<nrgljubljana_interface::nrg_params_t>        = true;
template <> inline constexpr auto c2py::tp_name<nrgljubljana_interface::nrg_params_t>    = "nrgljubljana_interface.solver_core.NrgParamsT";
template <> constexpr bool c2py::is_wrapped<nrgljubljana_interface::solver_core>         = true;
template <> inline constexpr auto c2py::tp_name<nrgljubljana_interface::solver_core>     = "nrgljubljana_interface.solver_core.SolverCore";
#endif