.. _2ord-UJ_QS:

2orb-UJ (QS)
============

This is the two-orbital Anderson impurity model with conserved
total charge (Q) and total spin (S) quantum numbers.

Hamiltonian
-----------

.. math::

  H_\mathrm{imp} = \sum_{i} \epsilon_i n_i + \sum_i U_i n_{i,\uparrow} n_{i,\downarrow} + U_{12} n_1 n_2 + J_{12} \mathbf{S}_1 \cdot \mathbf{S}_2

where :math:`\mathbf{S}_i` are the spin operators:

.. math::

  \mathbf{S}_i = \frac{1}{2} d_{i\alpha}^\dagger \boldsymbol{\sigma}_{\alpha\beta} d_{i\beta}


Parameters
----------

- :math:`\epsilon_i`, ``eps1`` and ``eps2``, energy levels
- :math:`U_i`, ``U1`` and ``U2``, electron-electron interaction
- :math:`U_{12}`, ``U12``, inter-level charge repulsion
- :math:`J_{12}`, ``J12``, inter-level exchange (Hund's) coupling

Expectation values
------------------

- :math:`\langle n_1 \rangle`, ``n_d1``, orbital 1 impurity occupancy
- :math:`\langle n_1^2 \rangle`, ``n_d1^2``, orbital 1 impurity occupancy squared
- :math:`\langle n_2 \rangle`, ``n_d2``, orbital 2 impurity occupancy
- :math:`\langle n_2^2 \rangle`, ``n_d2^2``, orbital 2 impurity occupancy squared
- :math:`\langle n_1 n_2 \rangle`, ``n_d1n_d2``, charge correlation
- :math:`\langle \mathbf{S}_1 \cdot \mathbf{S}_2 \rangle`, ``S_d1S_d2``, spin correlation

Structure of Green's functions
------------------------------

Single block ``imp``, 2x2 matrix-valued

Dynamic susceptibilities
------------------------

Dynamic spin and charge susceptibilities are calculated.
