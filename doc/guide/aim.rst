.. _aim:

An example: the Anderson impurity model
=======================================

To illustrate how the NRG solver works in practice, we show the example of
one-orbital Anderson impurity embedded in a flat (Wilson) conduction bath. The
interacting part of the local Hamiltonian of this problem is simply:

.. math::

  \mathcal{H}_\mathrm{int} = U n_\uparrow n_\downarrow,

and the non-interacting Green's function is:

.. math::

    G^{-1}_{0,\sigma} (\omega) = \omega - \epsilon_f - V^2 g(\omega).

In this example, there is a Coulomb interaction :math:`U` on the impurity level,
which is at an energy :math:`\epsilon_f`. The bath Green's function is
:math:`g(\omega)`, and it has a flat density of states over the interval
:math:`[-1,1]`.  Finally, :math:`V` is the hybridization amplitude between the
impurity and the bath. Let us solve this problem with the NRG solver. Here is
the python :download:`script <aim.py>`:

.. literalinclude:: aim.py

Running this script takes a few minutes and generates an
HDF5 archie file called :file:`aim_solution.h5`. This file contains
the impurity spectral function, Green's function, auxiliary Green's function,
self-energy, and expectation values of local variables 
(:math:`\langle n \rangle` and :math:`\langle n^2 \rangle`). Let us plot
the spectral function:

.. plot:: guide/aim_plot.py
   :include-source:
   :scale: 100

As expected, the result shows a particle-hole symmetric impurity
spectral function.

Let us now go through the script in some more detail:

.. literalinclude:: aim.py
  :lines: 1-2

These lines import the classes to run the NRG solver and to store the results
in a HDF5 archive file.

.. literalinclude:: aim.py
  :lines: 5-7

These lines define the impurity parameters, the (half)bandwidth :math:`D=1`,
and the physical temperature :math:`T`.

.. literalinclude:: aim.py
  :lines: 10

Here we construct the Solver object. The most important parameters here
are the model name and its symmetry type, as well as the mesh parameters which
define the logarithmic discretization mesh on which all quantities (input and
output) are defined. The models are described in template files. The program
will automatically determine the correct block structure of hybridisation
and Green's functions from the templates. This ensures that all quantities
within Solver are correctly initialized, and in particular that the block structure
of all Green's function objects is consistent.

.. literalinclude:: aim.py
  :lines: 13

This line defines, in addition to the physical temperature,
the key solve parameters which affect the quality of the results: the
discretization parameter :math:`\Lambda`, the number of the interleaved
discretization meshes :math:`N_z`, the scale ``Tmin`` which defines the length
of the Wilson chain (and which should be lower than the lowest physical temperature
scale in the problem, in this case lower than the Kondo temperature), and
the trucation settings which are determined by two parameters, the maximum number
of states kept and the truncation energy cutoff. There are many other parameters,
but the default values are suitable for our purposes.

.. literalinclude:: aim.py
  :lines: 16-17

Finally, we have set the model parameters. These are contained in a separate Python
dictionary.

.. literalinclude:: aim.py
  :lines: 20-21

The low-level NRG parameters are contained in yet another Python dictionary.
Here we set the NRG energy window to the interval :math:`[-1:1]`, since
the default is that the energy window is the same as the extent of the
logarithmic mesh, which in this example is :math:`[-2:2]`, i.e., larger
than the support of the hybridisation function which we set in this line:

We set the hybridisation function:

.. literalinclude:: aim.py
  :lines: 24

We can now call the solver:

.. literalinclude:: aim.py
  :lines: 27

We store the results in a HDF5 archive:

.. literalinclude:: aim.py
  :lines: 30-35

The expectation values can be accessed through the Python dictionary
``S.expv``:

.. literalinclude:: aim.py
  :lines: 37

A version of this documentation example can also be found in {\tt tutorial/1_AIM}
in the form of a Python script and a Jupyter (IPython) notebook.
