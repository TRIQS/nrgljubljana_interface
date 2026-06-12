###############################################################################
#
# nrgljubljana_interface: A TRIQS based impurity solver
#
# Copyright (c) 2019 The Simons foundation
#   authors: Nils Wentzell
#
# nrgljubljana_interface is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# nrgljubljana_interface is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# nrgljubljana_interface. If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################
"""
High-level interface to the NRGLjubljana impurity solver.

Defines the :class:`Solver` convenience subclass over the compiled
:class:`~nrgljubljana_interface.solver_core.SolverCore`, accepting the construction,
solve, and low-level NRG parameters as keyword arguments and importing
:mod:`triqs.utility.mpi` so that MPI is initialized before the solver runs.
"""
from .solver_core import SolverCore, ConstrParamsT, SolveParamsT, NrgParamsT

from triqs.gf import *
from triqs.utility import mpi


# === The SolverCore Wrapper

class Solver(SolverCore):
    r"""
    NRGLjubljana impurity solver.

    Thin Python wrapper over the compiled
    :class:`~nrgljubljana_interface.solver_core.SolverCore`. Construct it with the model,
    symmetry type and frequency mesh, assign the hybridization function
    ``S.Delta_w['imp'] << ...``, then call :meth:`solve` and read the results
    (``S.A_w``, ``S.G_w``, ``S.Sigma_w``, ...).

    The :meth:`__init__`, :meth:`solve` and :meth:`set_nrg_params` methods each accept
    either keyword arguments (which are used to build the corresponding parameter object)
    or a single parameter object positionally, so a stored ``S.constr_params`` /
    ``S.last_solve_params`` / ``S.nrg_params`` can be passed straight back in.

    Parameters
    ----------
    args : ConstrParamsT, optional
        A construction-parameter object (e.g. a stored ``constr_params``) passed positionally,
        instead of the keyword arguments below.
    model : str, optional
        Impurity model to solve (selects a template directory), e.g. ``'SIAM'``. Default ``'SIAM'``.
    symtype : str, optional
        NRGLjubljana symmetry code, e.g. ``'QS'``, ``'QSZ'``, ``'ISO'``. Default ``'QS'``.
    mesh_max : float, optional
        Maximum frequency of the logarithmic mesh. Default ``10``.
    mesh_min : float, optional
        Minimum frequency of the logarithmic mesh. Default ``1e-4``.
    mesh_ratio : float, optional
        Common ratio of the geometric (logarithmic) frequency mesh. Default ``1.05``.

    Notes
    -----
    See :class:`~nrgljubljana_interface.solver_core.ConstrParamsT` for the full list of
    construction parameters.
    """

    def __init__(self, *args, **params_kw):
        cp = args[0] if args else ConstrParamsT(**params_kw)
        SolverCore.__init__(self, cp)

    def solve(self, *args, **params_kw):
        """
        Solve the impurity problem.

        Parameters
        ----------
        args : SolveParamsT, optional
            A solve-parameter object (e.g. a stored ``last_solve_params``).
        params_kw : dict {'param':value}
            Solve parameters (Lambda, Nz, T, model_parameters, ...); see SolveParamsT.
        """
        sp = args[0] if args else SolveParamsT(**params_kw)
        return SolverCore.solve(self, sp)

    def set_nrg_params(self, *args, **params_kw):
        """
        Set the advanced (low-level) NRG parameters.

        Parameters
        ----------
        args : NrgParamsT, optional
            A low-level NRG-parameter object (e.g. a stored ``nrg_params``).
        params_kw : dict {'param':value}
            Low-level NRG parameters; see NrgParamsT.
        """
        np = args[0] if args else NrgParamsT(**params_kw)
        return SolverCore.set_nrg_params(self, np)
