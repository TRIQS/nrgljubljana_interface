.. _SIAM_QSZ:

SIAM (QSZ)
==========

This is the single-impurity single-orbital Anderson model with conserved
total charge (Q) and z-component of total spin (Sz) quantum numbers.

Hamiltonian
-----------

.. math::

  H_\mathrm{imp} = \epsilon_1 n + U_1 n_\uparrow n_\downarrow + B_1 S_z

with

.. math::

  S_z = \frac{1}{2} \left( n_\uparrow - n_\downarrow \right)

Parameters
----------

- :math:`\epsilon_1`, ``eps1``, energy level
- :math:`U_1`, ``U1``, electron-electron interaction
- :math:`B_1`, ``B1``, magnetic field (Zeeman splitting)

Expectation values
------------------

- :math:`\langle n \rangle`, ``n_d``, impurity occupancy
- :math:`\langle n^2 \rangle`, ``n_d^2``, impurity occupancy squared
- :math:`\langle \sum_\sigma d^\dagger_\sigma f_{0\sigma} + \text{h.c.} \rangle`, ``hop0``, hopping between the impurity and the zero-th site of the Wilson chain
- :math:`\langle S_z \rangle`, ``SZd``, spin polarization (magnetization)

Structure of Green's functions
------------------------------

Two blocks, ``up`` and ``dn``, scalar-valued (1x1 matrix)

Dynamic susceptibilities
------------------------

Dynamic spin and charge susceptibilities are calculated.
