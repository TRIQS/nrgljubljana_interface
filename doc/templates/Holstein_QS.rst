.. _Holstein_QS:

Holstein (QS)
=============

This is the single-orbital Anderson-Holstein model with a single phonon mode.
The phonon cutoff is part of the Hamiltonian specification and it is hard-coded to 10;
the model name is then e.g. ``Holstein/Nph=10``.
Conserved total charge (Q) and total spin (S) quantum numbers.

Hamiltonian
-----------

.. math::

  H_\mathrm{imp} = \epsilon_1 n + U_1 n_\uparrow n_\downarrow + g_1 (a+a^\dagger) (n - n_1) + \omega a^\dagger a

Parameters
----------

- :math:`\epsilon_1`, ``eps1``, energy level
- :math:`U_1`, ``U1``, electron-electron interaction
- :math:`g_1`, ``g1``, electron-phonon coupling
- :math:`n_1`, ``n1``, offset (charge reference point for e-ph coupling)
- :math:`\omega`, ``omega``, phonon frequency

Expectation values
------------------

- :math:`\langle n \rangle`, ``n_d``, impurity occupancy
- :math:`\langle n^2 \rangle`, ``n_d^2``, impurity occupancy squared
- :math:`\langle \sum_\sigma d^\dagger_\sigma f_{0\sigma} + \text{h.c.} \rangle`, ``hop0``, hopping between the impurity and the zero-th site of the Wilson chain
- :math:`\langle a^\dagger a \rangle`, ``nph``, phonon number expectation value
- :math:`\langle a+a^\dagger \rangle`, ``displ``, phonon displacement
- :math:`\langle (a+a^\dagger)^2 \rangle`, ``displ^2``, phonon displacement squared

Structure of Green's functions
------------------------------

Single block ``imp``, scalar-valued (1x1 matrix)

Dynamic susceptibilities
------------------------

Dynamic spin and charge susceptibilities are calculated.

Additional quantities
---------------------

Phonon spectral function.
