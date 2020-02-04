.. _parallelization:

Parallelization
===============

The code spends most of the time in LAPACK diagonalization routines (dsyev, dsyevr, zheev, etc.)
and in level-3 BLAS matrix-matrix multiplication routines (dgemm, zgemm). The execution time can
be significantly reduced by using a high-quality multithreaded BLAS/LAPACK libraries, such as
Intel MKL. Make sure that the environment variables (OMP_NUM_THREADS, MKL_NUM_THREADS, MKL_DYNAMIC)
are set appropriately.

In addition, the diagonalisation of matrices can be MPI parallelized. This is, however, only beneficial for
large-scale multi-orbital calculations.

TODO: The z-averaging is planned to be implemented through nested MPI parallelization.
