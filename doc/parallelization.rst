.. _parallelization:

Parallelization
===============

The code spends most of the time in LAPACK diagonalization routines (dsyev, dsyevr, zheev, etc.)
and in level-3 BLAS matrix-matrix multiplication routines (dgemm, zgemm). The execution time can
be significantly reduced by using a high-quality multithreaded BLAS/LAPACK libraries, such as
Intel MKL. Make sure that the relevant environment variables (e.g. OMP_NUM_THREADS, MKL_NUM_THREADS, MKL_DYNAMIC)
are set appropriately. Good choices are e.g. 4 or 8 threads.

The calculations over the different values of the twist parameter :math:`z` (different interleaved discretization grids)
are done in parallel using MPI parallelization. Good choice for the number of MPI processes is the value of ``Nz``, i.e.,
the number of different grids.

Finally, the diagonalisation of matrices can be OpenMP parallelized. This is, however, only beneficial for
large-scale multi-orbital calculations and is quite seldom used. This is configured using the parameter ``diagth``
in nrg_params_t ("low-level NRG parameters"), the default value being 1 (diagonalisations performed in series).
