.. _documentation:

Documentation
*************

Examples
--------

.. toctree::
   :maxdepth: 1

   guide/aim
   guide/holstein
   guide/hubbard
   guide/2orb-UJ

Model templates
---------------

.. toctree::
   :maxdepth: 1

   templates/SIAM_QS
   templates/SIAM_QSZ
   templates/Holstein_QS
   templates/2orb-UJ_QS
   guide/defining_models

Miscelaneous
------------

.. toctree::
   :maxdepth: 1

   parallelization
   temporary

C++ reference manual
--------------------

The C++ reference manual can be found `here <./doxygen/index.html>`_.

Python reference manual
-----------------------

.. autosummary::
   :toctree: _ref
   :template: autosummary_module_template.rst
   :recursive:

   nrgljubljana_interface.solver
   nrgljubljana_interface.solver_core

Links to all relevant solver parameters:

* :doc:`Construction parameters <_ref/nrgljubljana_interface.solver_core.ConstrParamsT>`
* :doc:`Solve parameters <_ref/nrgljubljana_interface.solver_core.SolveParamsT>`
* :doc:`Low-level NRG parameters <_ref/nrgljubljana_interface.solver_core.NrgParamsT>`
