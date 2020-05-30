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
from solver_core import SolverCore

from triqs.gf import *
from triqs.utility import mpi


# === The SolverCore Wrapper

class Solver(SolverCore):
    def __init__(self, **params_kw):
        """
        Initialise the solver.

        Parameters
        ----------
        """
        # Initialise the core solver
        SolverCore.__init__(self, **params_kw)

    def solve(self, **params_kw):
        """
        Solve the impurity problem.

        Parameters
        ----------
        params_kw : dict {'param':value} that is passed to the core solver.
        """

        # Call the core solver's solve routine
        return SolverCore.solve(self, **params_kw)
