.. _temporary:

Temporary files
===============

The code uses file-based inter-process communication to interface with the NRG Ljubljana code ("temporary files").
In addition, the NRG Ljubljana code itself stores temporary files during the calculation: in the first NRG sweep it stores
all eigenvalues and eigenvectors, then it calculates all the density matrices, which are then being read during
the second NRG sweep ("working directory"). Both default to the current directory, but can be configured using
environment variables NRG_TEMPDIR and NRG_WORKDIR. It is recommended that the first points to a fast storage
system (/dev/shm shared memory filesystem is perfectly adapted for this purpose, since in TEMPDIR we store small
files but we perform quite intensive I/O during the "instantiation" step of the calculation). The second should
point to a file system with sufficient storage space and ample bandwidth, i.e. either local or global
scratch storage. The storage space requirement for TEMPDIR can be estimated from the memory use at the NRG
steps with the largest number of retained states: this number should be multiplied by a factor of 20-30. (This
is only a rough estimate, YMMV.)


