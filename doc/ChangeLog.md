(changelog)=

# Changelog

## Version 3.1.0

nrgljubljana_interface version 3.1.0 is a compatibility
release for TRIQS version 3.1.0.


## Version 3.0.0

nrgljubljana_interface version 3.0.0 is a compatibility
release for TRIQS version 3.0.0 that introduces compatibility
with Python 3 (Python 2 no longer supported)

We provide a more detailed description of the changes below.

General
-------
* Update nrgljubljana to release 2020.09
* Link cpp2py_module against triqs_py library
* Minor cleaning in doc/guide/hubbard.py
* Rename h5 output file for refreq_pts test
* Do not use TRIQS_MAKE_NVP macro in serialization function of point_mesh
* improved handling of chdir to temp dir (and return); libmkl_rt
* report OpenMP settings
* add blas/lapack dependency
* Note on temporary storage, NRG_WORKDIR and NRG_TEMPDIR environment variables
* temporary directory configurable via env variable NRG_TEMPDIR.. it's best to set it to a RAM disk, e.g. /dev/shm/ on Linux systems
* update DMFT code
* comment out unused section of code
* remove unused variables
* pass by ref arg to call()
* copy&move constructors & operator= for Interpolation and Integration classes
* project description added
* DMFT tutorial improvements
* DMFT tutorial cleanup + fixes
* object-oriented DMFT(NRG) implementation
* remove old files

doc
---
* Update version number in landing page
* Small simplification in A_to_nparrays function

jenkins
-------
* Be sure to add custom boost to LD_LIBRARY_PATH
* Install perl-Math-Complex for centos 8

cmake
-----
* Synchronize c++/nrgljubljana_interface/CMakeLists.txt with app4triqs skeleton

py3
---
* In point_mesh, use new signature of write/assert_hdf5_format_as_string
* Port all python files and notebooks to Python3
* Apply ports_to_triqs3 porting script

Contributors: Philipp Dumitrescu, Dylan Simon, Nils Wentzell, Rok Zitko


## Version 2.2.0

nrgljubljana_interface provides a C++ and Python
interfaces between "NRG Ljubljana" numerical
renormalization group code and TRIQS. It is
based on problem template files and does not
require Mathematica during computation.
Mathematica is only required for generating new
template files for custom Hamiltonians and 
symmetry types.

This is the initial release for this project.
