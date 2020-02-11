.. _2orbital_aim:

Two-orbital Anderson impurity model
===================================

Here we illustrate calculations for multi-orbital models on the example
of the two-orbital Anderson impurity model with inter-orbital charge
and exchange (Hund's) coupling with the Hamiltonian

.. math::

  H_\mathrm{imp} = \sum_{i} \epsilon_i n_i + \sum_i U_i n_{i,\uparrow} n_{i,\downarrow} + U_{12} n_1 n_2 + J_{12} \mathbf{S}_1 \cdot \mathbf{S}_2

  where :math:\mathbf{S}_i are the spin operators:

  .. math::

    \mathbf{S}_i = \frac{1}{2} d_{i\alpha}^\dagger \boldsymbol{\sigma}_{\alpha\beta} d_{i\beta}

Let us solve this problem with the NRG solver. Here is
the python :download:`script <2orb-UJ.py>`:

.. literalinclude:: 2orb-UJ.py

Running this script takes quite a bit of time and generates an HDF5 archive file called :file:`2orb-UJ_solution.h5`.
This file contains the impurity spectral function, Green's function, auxiliary Green's function,
self-energy, and some expectation values of local variables. Let us plot the spectral function:

.. plot:: guide/2orb-UJ_plot.py
   :include-source:
   :scale: 100

