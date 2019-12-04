.. _aim:

An example: the Anderson impurity model
=======================================

To illustrate how the NRG solver works in practice, we show the example of a
one-orbital Anderson impurity embedded in a flat (Wilson) conduction bath. The
interacting part of the local Hamiltonian of this problem is simply:

.. math::

  \mathcal{H}_\mathrm{int} = U n_\uparrow n_\downarrow,
  
  and the non-interacting Green's function is:
  
  .. math::
  
    G^{-1}_{0,\sigma} (\omega) = \omega - \epsilon_f - V^2 \Delta_\sigma(\omega).
    
In this example, there is a Coulomb interaction :math:U on the impurity level,
which is at an energy :math:\epsilon_f. The bath Green's function is :math:\Delta(
\omega), and it has a flat density of states over the interval
:math:[-1,1].  Finally, :math:V is the hybridization amplitude between the
impurity and the bath. Let us solve this problem with the NRG solver. Here is
the python :download:script <aim.py>:

.. literalinclude:: aim.py

