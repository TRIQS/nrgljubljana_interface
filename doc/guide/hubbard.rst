.. _hubbard:

DMFT for Hubbard and Hubbard-like models
========================================

We now present a Python script for performing DMFT calculations with the NRG
as the impurity solver. The script is fairly complete and robust, and as such it is suitable for production runs.
It supports a variety of models with different symmetry types (including cases
with block structure and with matrix-valued Green's functions), controls
the chemical potential to achieve the selected band occupancy, performs linear or Broyden mixing 
to accelerate the convergence, stores intermediate results (checkpoints) to enable
restarts, and supports arbitrary density of states in the band (e.g. tabulated in a file).
Here is the :download:`listing <hubbard.py>`:

.. literalinclude:: hubbard.py

Running this script takes between an hour or several hours, depending on the
degree of parallelization that is enabled and the available computing resources. 
During the run, some basic information about the DMFT iteration is stored in the file ``stats.dat``:

.. literalinclude:: stats.dat

The first column is the iteration number, the second the current estimate of the chemical potential,
the following three columns provide information about the convergence (difference between local lattice
and impurity Green's functions, difference between two consecutive local lattice Green's function, and
the difference between target and current band occupancy). This is followed by the occupancy calculed from
the spectral function and by the expectation values of a range of local operators (the set depends on the
model considered and on the symmetry type).

Note: NRG impurity solver produces big temporary files during the calculation. By default, these are
stored in the working directory that is created in the current directory. To override this behavior,
set the environment variable NRG_WORKDIR. A suitable choice is a local scratch hard drive or SSD,
or a global scratch mount on a high-performance filesystem.

When convergence is reached, the final results are stored in an HDF5 archive file called :file:`solution.h5`. 
This file contains all quantities of interest. Let us plot the spectral function:

.. plot:: guide/hubbard_plot.py
   :include-source:
   :scale: 100

Let us now go through the script in some more detail. We first import a number of useful modules:

.. literalinclude:: hubbard.py
  :lines: 4-12

We make use of scipy interpolation and integration routines to obtain highly accurate results based on
the real-frequency spectra, and we make use of the Broyden solver from ``scipy.optimize``.

The model name is a full path to the model definition in the template library, except the symmetry part
(e.g. "SIAM", "Holstein/Nph=10", etc.).

.. literalinclude:: hubbard.py
  :lines: 15

The code comes with a built-in support for the Bethe lattice (with half-bandwidth set by the parameter ``Bethe_unit``),
or one can use tabulated DOS (example :download:`script <dos.py>` for generating the DOS file).

.. literalinclude:: hubbard.py
  :lines: 16-18

We then set key parameters, such as Hubbard :math:`U`, the occupancy, the temperature, and optionally the magnetic field.
The parameters that are not used should be set to ``None``.

.. literalinclude:: hubbard.py
  :lines: 19-23

One can control the level of verbosity and the amount of intermediate results that are stored to disk:

.. literalinclude:: hubbard.py
  :lines: 27-29

The DMFT iteration parameters are set next. Note that there are several convergence criteria: they must all
be satisfied.

.. literalinclude:: hubbard.py
  :lines: 30-34

To accelerate the convergence to the converged result, one can make use of two different mixing schemes: linear mixing
(with parameter :math:`\alpha`) or Broyden mixing (with initial Jacobian :math:`J=-1/\alpha`). The direct
DMFT iteration corresponds to linear mixing with :math:`\alpha=1` which actually corresponds to no mixing at all.

.. literalinclude:: hubbard.py
  :lines: 35-37
  
For adjusting the chemical potential, it is useful to go through several DMFT self-consistency cycles (using the
same self-energy obtained in the last NRG run) to obtain an improved estimate. In practice order 10 cycles are quite sufficient in all cases.

.. literalinclude:: hubbard.py
  :lines: 38
  
The logarithmic discretization mesh should span the full frequency range where the spectral function has non-negligible
value. The common ratio should be sufficiently small to produce a fine discretization mesh.

.. literalinclude:: hubbard.py
  :lines: 39-41
  
The NRG has issues with regions of very low density of states in the band. To control those we define a minimal
value of the hybridisation strenth :math:`\Gamma=-\mathrm{Im}(\Delta)`:

.. literalinclude:: hubbard.py
  :lines: 42-43

The Python dictionary ``observables`` contains the list of all operators whose expectation values we want to evaluate.
Some standard choices are the impurity occupancy and impurity occupancy squared, in the presence of magnetic field
also the impurity magnetization, and in the presence of phonon modes we might be interested in the phonon number
and oscillator displacement operators.

.. literalinclude:: hubbard.py
  :lines: 49-53
  
We then set some further global variables that are used in the DMFT iteration:

.. literalinclude:: hubbard.py
  :lines: 56-59
  
One may notice that in the presence of the magnetic field, the code automatically switches to the appropriate
symmetry type (``QSZ``, meaning conserved total charge and the z-component of total spin). 

We are now ready to construct the impurity solver object:

.. literalinclude:: hubbard.py
  :lines: 62-63
  
We usually want to suppress the detailed output from the NRG by setting ``verbose_nrg`` to False. (In fact,
what this does is to redirect the output to a file named ``log`` in the temporary directory, so that 
the calculation can still be monitored.)

We define some useful functions for working with BlockGf objects:

.. literalinclude:: hubbard.py
  :lines: 65-68
  
The key NRG solver parameters that control the discretization and the truncation are set here:

.. literalinclude:: hubbard.py
  :lines: 70-71

The default value of the discretization parameter :math:`\Lambda=2` is an excellent choice for most cases, as long as the calculation remains
feasible in reasonable time with the given computing resources; if not, one should increase :math:`\Lambda`.
We average over ``Nz=4`` interleaves discretization grids; this is a good number for :math:`\Lambda=2`, but should be
increased for coarses discretization to reduce the artifacts. Furthermore, it needs to be increased for high-resolution calculations
where the broadening parameters are lower.
The parameter ``Tmin`` controls the length of the Wilson chain and it has to be smaller than the physical temperature 
:math:`T`. The truncation at rescaled energy 10 is sufficient to achieve good convergence of the spectral functions,
while beyond 12 there are hardly any improvements. ``keep`` (the maximum number of multiplets kept) should
be some large number and it is basically constrained by the available memory.

Model parameters are contained as a separate Python dictionary. We set it here:

.. literalinclude:: hubbard.py
  :lines: 74-77
  
We remove parameters that are not relevant (those set to ``None``). Furthermore, the interface checks if
all required parameters for a given template are in fact defined, to prevent user errors that are
difficult to detect (e.g. typos).

We provide a function for updating the chemical potential. It sets the global variable ``mu``, as well
as updates the impurity level to :math:`\epsilon_d=-\mu`.

.. literalinclude:: hubbard.py
  :lines: 80-83
  
If some details of the NRG calculation need to be tweaked, this can be achieve by setting the appropriate
parameters in ``nrg_params``.

.. literalinclude:: hubbard.py
  :lines: 87-89
  
Now we define the function ``ht0`` which evaluates the Hilbert transform of the density of states for
a given complex value :math:`z`. The analytical expression for the Bethe lattice is hard-coded, while the
generic case is handled by reading the DOS from a file and doing the Hilbert transform numerically.

.. literalinclude:: hubbard.py
  :lines: 93-103
  
We wrap this in another function ``ht`` which ensures that in cases of causality violation (imaginary
part of the self-energy being positive) things do no go terribly wrong. Tiny causality violations in NRG
commonly occur at very low temperatures close to the Fermi level, where :math:`\mathrm{Im}\Sigma` 
overshoots zero. This problem is well understood and, while annoying, it is usually not an issue, unless
one is explicitly interested in the low-temperature asymptotics.

.. literalinclude:: hubbard.py
  :lines: 106-107
  
In addition, another potential problem in NRG calculation is the artifacts arising from regions of very low DOS.
We alleviate this problem by providing a minimal cutoff value ``Delta_min`` for the hybridisation function:

.. literalinclude:: hubbard.py
  :lines: 281-290
  
In the following we will make use of a function which calculates a Green's function given hybridisation, self-energy and
chemical potential:

.. literalinclude:: hubbard.py
  :lines: 110-115

For accurate integration over spectral functions, we make use of cubic interpolation. Note that for
frequencies outside the mesh range, the spectral function is assummed to be strictly zero.

.. literalinclude:: hubbard.py
  :lines: 121-128

This interpolation object is then used e.g. to compute the occupancy given hybridisation, self-energy and
chemical potential:

.. literalinclude:: hubbard.py
  :lines: 131-137

And this, in turn, is used for occupancy control. The following function determines the optimal value of the
chemical potential for given hybridisation and self-energy. It makes use of ``root_scalar``, which implements
the secant method with initial values shifted by 0.1. This seems to be robust and works quite well in practice.

.. literalinclude:: hubbard.py
  :lines: 140-142
  
Now comes the key function which calculates the lattice Green's function, and from this a new
hybridisation function for the effective impurity model. The input is the self-energy (and implicitly the
chemical potential ``mu``):

.. literalinclude:: hubbard.py
  :lines: 146-162
  
Note that this implementation is limited to problems with diagonal self-energy (and we test for this
explicitly, as a precaution). This means that this code is only suitable for channel/orbital conserving models.
  
Now we are in position to provide an improved method for adjusting the chemical potential, which
repeatedly recalculates the hybridisation function and readjust the chemical potential. We do this
for ``max_mu_adjust`` cycles, with 10 being a good choice.

.. literalinclude:: hubbard.py
  :lines: 166-174
  
We now define some miscelaneous functions. First we need a way of estimating the difference between two
Green's functions. This will be used in the following for convergence testing. Here we again make use
of interpolation, and we define the difference as the integral over all frequencies (within the discretization
mesh window) of the difference squared. We ignore possible warnings about loss of accuracy in integration.

.. literalinclude:: hubbard.py
  :lines: 178-184

We now define some functions for saving Green's function object to tabulated ASCII files. The filename is just
the prefix. To this we append the block name and the suffix ``dat``.

.. literalinclude:: hubbard.py
  :lines: 186-208
  
In addition, the full set of results is stores as an HDF5 file. We store the relevant items from the solver
object, as well as some additional quantities (self-energy-improved impurity Green's function, local lattice Green's function,
chemical potential).

.. literalinclude:: hubbard.py
  :lines: 211-217
  
When restarting from stored results, we only need to read in the self-energy and the chemical potential, since this
is sufficient to initialize the DMFT calculation.

.. literalinclude:: hubbard.py
  :lines: 226-230

.. literalinclude:: hubbard.py
  :lines: 234-236
  
If not starting from a file, the initial self-energy includes the Hartree term, and the chemical potential is set to
:math:`U/2`:

.. literalinclude:: hubbard.py
  :lines: 239-246
  
The DMFT loop is terminated by raising exceptions. We define them here:

.. literalinclude:: hubbard.py
  :lines: 257-265

Now comes the main part, the function which performs a DMFT step. Let us disect it. First we update the iteration counter and
initialize the hybridisation function in the impurity solver object after fixing it to a minimal value.

.. literalinclude:: hubbard.py
  :lines: 294-300
  
We now solve the impurity model, calculate an improved estimate of the impurity GF using the "self-energy trick", and
apply the self-consistency equation:

.. literalinclude:: hubbard.py
  :lines: 302-305

We calculate some quantities for estimating the degree of convergence:

.. literalinclude:: hubbard.py
  :lines: 307-311

.. and we report some statistics to console and to file ``stats.dat``:

.. literalinclude:: hubbard.py
  :lines: 313-322

Furthermore, if ``store_steps=True``, we save the information about each DMFT step in a separate subdirectory:

.. literalinclude:: hubbard.py
  :lines: 324-330

At each step, we also store the HDF5 to permit the checkpoint/restart functionality:

.. literalinclude:: hubbard.py
  :lines: 332-333
  
We now test for convergence (and for failure to reach the convergence within the maximal number of iterations):

.. literalinclude:: hubbard.py
  :lines: 336-342
  
We are nearly done now. We adjust the chemical potential and return the optimized hybridisation function:

.. literalinclude:: hubbard.py
  :lines: 344-347

There are two types of driving routines for the full DMFT iteration. The first is based on linear mixing:

.. literalinclude:: hubbard.py
  :lines: 350-355
  
The second one is based on Broyden mixing. The idea is that the DMFT iteration can be thought of as a
procedure for finding the fixed point of the mutli-dimensional function :math:`F(\Delta_\mathrm{in}) = \Delta_\mathrm{out}[\Delta_\mathrm{in}]  - \Delta_\mathrm{in}`,
where :math:`\Delta_\mathrm{in}` is the input while :math:`\Delta_\mathrm{out}` the output of one DMFT step.
For ths purpose we make use of ``optimize.broyden1`` from scipy, but in order to do so, we
need to wrap/unwrap the Green's function objects to/from numpy arrays, which adds a degree of complexity:

.. literalinclude:: hubbard.py
  :lines: 357-383

Here we pick the method:

.. literalinclude:: hubbard.py
  :lines: 385-389

The problem is initialized here:

.. literalinclude:: hubbard.py
  :lines: 392-396

And the solver is started here:

.. literalinclude:: hubbard.py
  :lines: 402-413

This completes the exposition. Have fun running DMFT(NRG) calculations!


Note: this script has not yet been tested for the case of matrix-valued Green's functions.
