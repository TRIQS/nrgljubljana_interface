.. _holstein:

Anderson-Holstein impurity model
================================

Here we consider an impurity model with electron-phonon coupling, to wit the Anderson-Holstein impurity model. 
The interacting part of the local Hamiltonian of this problem is:

.. math::

  \mathcal{H}_\mathrm{int} = U n_\uparrow n_\downarrow + g (a+a^\dagger)(n-n_1) + \omega a^\dagger a,

and the non-interacting Green's function is:

.. math::

    G^{-1}_{0,\sigma} (\omega) = \omega - \epsilon_f - V^2 g(\omega).

The new element (compared to the Anderson impurity model) is the presence of
a phonon mode with frequency :math:`\omega`, which couples through e-ph
coupling strength :math:`g`. Finally, :math:`n_1` is a shift, typically either
0, 1 or :math:`\langle n \rangle`.

Here is the python :download:`script <holstein.py>`:

.. literalinclude:: holstein.py

Running this script takes a few minutes and generates an HDF5 archive file called :file:`holstein_solution.h5`.
Let us plot the spectral function:

.. plot:: guide/holstein_plot.py
   :include-source:
   :scale: 100

As expected, the result shows a particle-hole symmetric impurity spectral function with a (charge) Kondo
resonance and side peaks at multiples of the phonon frequency :math:`\omega`;

Let us now go through the script in some more detail, focusing on the differences with respect to the
standard Anderson impurity model.

.. literalinclude:: holstein.py
  :lines: 5-10

Here we set the parameters. Note that we have set the Hubbard parameter to zero, thus this is
a pure Holstein impurity problem, with the correlation effects stemming solely from the electron-phonon
coupling. We set :math:`n_1` to 1, to ensure particle-hole symmetry.

.. literalinclude:: holstein.py
  :lines: 13

Here we construct the Solver object. Note the model name which contains the phonon number cutoff (10).
In NRGLjubljana_interface, the model names are actually paths to template files (either those
bundled with the interface in the directory ``templates/``, or custom templates created by the user). In the case of the Anderson-Holstein
model, the template files need to be generated for a specific value of the phonon cutoff.

.. literalinclude:: holstein.py
  :lines: 36-38

The Anderson-Holstein model defines several additional expectation values pertaining to the
phonon degree of freedom.
