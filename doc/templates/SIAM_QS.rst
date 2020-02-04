.. _SIAM_QS:

SIAM (QS)
=========

This is the simplest single-impurity (and single-orbital) Anderson model with conserved
total charge (Q) and total spin (S) quantum numbers.

Hamiltonian
-----------

.. math::

  H_\mathrm{imp} = \epsilon_1 n + U_1 n_\uparrow n_\downarrow

Parameters
----------

- :math:`\epsilon_1`, ``eps1``, energy level
- :math:`U_1`, ``U1``, electron-electron interaction

Expectation values
------------------

- :math:`\langle n \rangle`, ``n_d``, impurity occupancy
- :math:`\langle n^2 \rangle`, ``n_d^2``, impurity occupancy squared
- :math:`\langle \sum_\sigma d^\dagger_\sigma f_{0\sigma} + \text{h.c.} \rangle`,
``hop0``, hopping between the impurity and the zero-th site of the Wilson chain

Structure of Green's functions
------------------------------

Single block ``imp``, scalar-valued (1x1 matrix)

Dynamic susceptibilities
------------------------

Dynamic spin and charge susceptibilities are calculated.
