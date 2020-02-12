.. _defining_models:

Generating new templates for custom models
==========================================

Here we discuss how to add new models to the template library.
As an example, we consider the case of two-orbital impurity model
(model ``2orb-UJ``) that comes predefined in the library of default templates.
The files are located in the directory ``templates/2orb-UJ/QS``.

The model is defined in ``2orb-UJ.m``:

.. literalinclude:: 2orb-UJ.m

.. _sneg: http://nrgljubljana.ijs.si/sneg

The first line declares that this is a two-channel problem with two impurity orbitals.
The second line defines how we label the impurity orbitals. Here we choose consecutively
numbered ``d[1]`` and ``d[2]``. The main part of the file then defines various parts
of the Hamiltonian, making use of sneg_ function
calls such as ``number``, ``Hubbard``, ``spinspin``, etc. Finally, in the last two lines
we define the auxiliary operators which are used in the "self-energy trick" to obtain
a high-quality estimate for the self-energy. They are obtained by taking the commutator
between the particle creation/annihilation operators and the interacting part of the
Hamiltonian.

The input file to NRG template generator is ``param``:

.. literalinclude:: param

Here we declare all parameters (although the values need not be given), instruct
``nrginit`` to generate a template file ``data.in`` rather than a specific instantiation
``data``, define the symmetry type (here ``QS``), and list all operators of interest.
The last line is mandatory for generating NRGLjubljana_interface template; it instructs
the code to generate Hamiltonians which are not scaled in terms of the initial energy
scale of the problem, as is customarily done in NRG calculations.

The operators of interest are defined in ``modeloperators.m``:

.. literalinclude:: modeloperators.m

Here we again make use of sneg functions, as well as of high-level functions defined
in ``initial.m`` which is part of the ``nrginit`` Mathematica code.

The template files are then obtained by running ``nrginit``. This will generate
``data.in``, ``ham_*`` and ``op.*`` template files. The log file with some details about
the generation process is saved as ``mmalog``.

The main template ``data.in`` is lengthy and we show only a small chunk:

.. literalinclude:: data.in_excerpt

The ``ham_*`` contain the symbolic expressions for the Hamiltonian matrices in various invariant subspaces.
For instance, in the half-filled :math:`Q=0` spin triplet :math:`(2S+1)=3` subspace,
we have ``ham_0.3``:

.. literalinclude:: ham_0.3

In a similar fashion, ``op.*`` files contain the symbolic forms of the operators. For instance
one part of the auxiliary operator ``self_d2`` takes the following form:

.. literalinclude:: op.self_d2_0.3_-1.4

We then need to define the structure of Green's functions, separately for
Green's functions and for dynamical susceptibility functions, in ``gf_stuct``
and ``chi_struct``, respectively. Here they are the same:

.. literalinclude:: gf_struct
.. literalinclude:: chi_struct

The quantities that need to be computed are listed in file ``info``:

.. literalinclude:: info

Finally, NRGLjubljana_interface requires a number of initialization, wrapper and processing scripts.
The first, ``prepare``, simply copies the required files to the temporary working directory:

.. literalinclude:: prepare

The second performs the discretization:

.. literalinclude:: discretize

Then comes the most complex one which instantiates the input file for NRG calculation, numerical-valued ``data``,
from the symbolic-expression template ``data.in``:

.. literalinclude:: instantiate

This script, in turn, calls another one that calculates the Wilson chain coefficients:

.. literalinclude:: wilson

The last script, ``process``, calculates averages over the twist parameter and performs the broadening
of spectral functions:

.. literalinclude:: process
