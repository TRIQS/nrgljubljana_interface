.. _SIAM_QS:

SIAM (QS)
=========

This is the simplest single-impurity single-orbital Anderson model with conserved
total charge (Q) and total spin (S) quantum numbers.

Hamiltonian
-----------

.. math::

  H_\mathrm{imp} = \epsilon_1 n + U_1 n_\uparrow n_\downarrow

where :math:`n_\sigma = d^\dagger_{\sigma} d_{\sigma}` and :math:`n=\sum_\sigma n_\sigma` with :math:`\sigma=\uparrow,\downarrow`.

Parameters
----------

- :math:`\epsilon_1`, ``eps1``, energy level
- :math:`U_1`, ``U1``, electron-electron interaction

The number 1 is an orbital index. The choice of including an index (which is not strictly needed here) is to ease writing
generic codes that support an arbitrary orbital degeneracy.

Expectation values
------------------

- :math:`\langle n \rangle`, ``n_d``, impurity occupancy
- :math:`\langle n^2 \rangle`, ``n_d^2``, impurity occupancy squared
- :math:`\langle \sum_\sigma d^\dagger_\sigma f_{0\sigma} + \text{h.c.} \rangle`, ``hop0``, hopping between the impurity and the zero-th site of the Wilson chain, :math:`f_0` (i.e., the local orbital of
the bath at the position of the impurity)

Structure of Green's functions
------------------------------

Single block ``imp``, scalar-valued (1x1 matrix)

Dynamic susceptibilities
------------------------

Dynamic spin and charge susceptibilities are calculated.
