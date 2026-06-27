(changelog)=

# Changelog

## Version 4.0.0

nrgljubljana_interface version 4.0.0 is a compatibility release for TRIQS version 4.0.0.
It is based on the latest [app4triqs](https://github.com/TRIQS/app4triqs) skeleton and
ports the Python bindings to the new clair + c2py binding generator. The most notable
changes are:

* Run the `port_to_triqs4` script for the TRIQS 4.0 module renames (`triqs.gf` -> `triqs.gfs`, mesh modules moved to `triqs.mesh`)
* Use the latest app4triqs skeleton, including the new Kubernetes-based Jenkins setup
* Port the Python bindings from cpp2py to clair + c2py
* Replace the custom `point_mesh` with the TRIQS `refreq_log` mesh and align with the latest mesh conventions
* Add improved self-energy estimators and new problem template files

We thank all contributors: Jennifer Coulter, Thomas Hahn, Alexander Hampel, Henri Menke, Don Rolih, Dylan Simon, Nils Wentzell, Rok Zitko

Find below an itemized list of changes in this release.

### General
* Port to TRIQS 4.0 via the `port_to_triqs4` script
* Port the Python bindings to clair + c2py and regenerate the `solver_core` module
* Replace the custom `point_mesh` with the TRIQS `refreq_log` mesh and align with the latest mesh conventions
* Add improved self-energy estimators and the definitions they require
* Add new problem template files (SIAM/ISO particle-hole-symmetric case, updates to SIAM/QS, SIAM/QSZ, Holstein/Nph=10/QS, 2orb-UJ/QS)
* Add a note and instructions on implementing new models in nrgljubljana_interface
* Update the bundled NRG Ljubljana dependency to b2dc884
* Update the doxygen documentation in `container_set.hpp`, `params.hpp` and `solver_core.hpp`
* Clean up references to the old estimator
* Fix unrecognized Python tests

### cmake
* Fetch the NRG Ljubljana dependency through the CPM package manager, pinned to a fixed `rokzitko/nrgljubljana` commit
* Allow both static and dynamic builds of the nrgljubljana library
* Disable the Mathematica search on OSX
* Fix the test environments

### cpp
* Replace the deprecated OpenMP function `omp_get_nested`
* Fix the calls to NRG Ljubljana functions

### python
* Update test scripts and reference `.h5` files for the new self-energy estimator
* Replace deprecated `numpy.complex_` with `numpy.complex128`

### jenkins
* Migrate to the new Kubernetes-based Jenkins system
* Add the GSL and Boost dependencies to `jenkins/Dockerfile`, bump Boost to 1.88
* Disable IntelLLVM builds

### ghactions
* Modernize and simplify `build.yml`, building in parallel with Ninja
* Add the GSL, Boost-MPI and Boost-system dependencies
* Fix the macOS + clang build (libc++/libomp/unwind paths, OpenMP detection)

### doc
* Add a Flatiron Institute support notice to `README.md`


## Version 3.3.0

nrgljubljana_interface version 3.3.0 is a compatibility release for TRIQS version 3.3.0.

We thank all contributors: Thomas Hahn, Alexander Hampel, Henri Menke, Nils Wentzell

Find below an itemized list of changes in this release.

### General
* Add intel-linux toolset to boost build in the Dockerfile

### ghactions
* Remove clang build and add libboost-serialization-dev as a dependency


## Version 3.2.0

nrgljubljana_interface version 3.2.0 is a compatibility release for TRIQS version 3.2.0.

We thank all contributors: Thomas Hahn, Alexander Hampel, Dylan Simon, Nils Wentzell

Find below an itemized list of changes in this release.

### General
* Remove use of std::unary_function
* Use std::lower_bound over std::ranges::lower_bound
* Raise hilbert test tolerance
* Port to TRIQS 3.2


## Version 3.1.0

nrgljubljana_interface version 3.1.0 is a compatibility
release for TRIQS version 3.1.0.
* uses the latest [app4triqs/3.1.x](https://github.com/TRIQS/app4triqs) skeleton
* switches documentation to read-the-docs theme
* improves ghactions and jenkins configuration

We provide a more detailed description of the changes below.

### General
* C++20 fix: REQUIRES -> requires
* Disable osx gcc builds due to boost standard library incompatibility
* Fix get_interpolation_data in point_mesh.hpp for boundary points
* Porting to triqs 3.1 + nda + meshes
* Partial implementation of the template 2orb-UJt/QS
* Support from different structure of GF and Delta
* CLEAN problem : single-channel no-impurity problem
* read-in tdfdm (thermodynamics) results ; has_struct() function ; support for no-impurity problems ; keep_temp_dir variable

### jenkins
* Disable sanitization build due to compilation timeouts
* Fix link to boost tarball

We thank all contributors: Alexander Hampel, Dylan Simon, Nils Wentzell, Rok Zitko


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
