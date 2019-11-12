# Generated automatically using the command :
# c++2py ../../c++/nrgljubljana_interface/solver_core.hpp -p --members_read_only -N nrgljubljana_interface -a nrgljubljana_interface -m solver_core -o solver_core --moduledoc="The nrgljubljana_interface solve_core module" -C pytriqs --cxxflags="-std=c++17" --target_file_only
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "solver_core", doc = r"The nrgljubljana_interface solve_core module", app_name = "nrgljubljana_interface")

# Imports

# Add here all includes
module.add_include("nrgljubljana_interface/solver_core.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/optional.hpp>
#include <cpp2py/converters/string.hpp>
#include <cpp2py/converters/vector.hpp>
#include <triqs/cpp2py_converters/h5.hpp>

using namespace nrgljubljana_interface;
""")


# The class solver_core
c = class_(
        py_type = "SolverCore",  # name of the python class
        c_type = "nrgljubljana_interface::solver_core",   # name of the C++ class
        doc = r"""The Solver class""",   # doc of the C++ class
        hdf5 = True,
)

c.add_member(c_name = "constr_params",
             c_type = "nrgljubljana_interface::constr_params_t",
             read_only= True,
             doc = r"""""")

c.add_member(c_name = "nrg_params",
             c_type = "nrgljubljana_interface::nrg_params_t",
             read_only= True,
             doc = r"""""")

c.add_member(c_name = "last_solve_params",
             c_type = "std::optional<solve_params_t>",
             read_only= True,
             doc = r"""""")

c.add_constructor("""(**nrgljubljana_interface::constr_params_t)""", doc = r"""Construct a NRGLJUBLJANA_INTERFACE solver



+----------------+-------------+---------+-----------------------------------------+
| Parameter Name | Type        | Default | Documentation                           |
+================+=============+=========+=========================================+
| bandrescale    | double      | 1.0     | Band rescaling factor                   |
+----------------+-------------+---------+-----------------------------------------+
| problem        | std::string | "SIAM"  | Model considered (templated)            |
+----------------+-------------+---------+-----------------------------------------+
| mesh_max       | double      | 10      | Mesh maximum frequency                  |
+----------------+-------------+---------+-----------------------------------------+
| mesh_min       | double      | 1e-4    | Mesh minimum frequency                  |
+----------------+-------------+---------+-----------------------------------------+
| mesh_ratio     | double      | 1.05    | Common ratio of the geometric sequence  |
+----------------+-------------+---------+-----------------------------------------+
""")

c.add_method("""void solve (**nrgljubljana_interface::solve_params_t)""",
             doc = r"""Solve method that performs NRGLJUBLJANA_INTERFACE calculation



+----------------+-------------+---------+------------------------------------------------------+
| Parameter Name | Type        | Default | Documentation                                        |
+================+=============+=========+======================================================+
| Lambda         | double      | 2.0     | Logarithmic discretization parameter                 |
+----------------+-------------+---------+------------------------------------------------------+
| Nz             | int         | 1       | Number of discretization meshes                      |
+----------------+-------------+---------+------------------------------------------------------+
| xmax           | double      | 20      | Largest x in discretization ODE solver               |
+----------------+-------------+---------+------------------------------------------------------+
| Tmin           | double      | 1e-4    | Lowest scale on the Wilson chain                     |
+----------------+-------------+---------+------------------------------------------------------+
| keep           | size_t      | 100     | Maximum number of states to keep at each step        |
+----------------+-------------+---------+------------------------------------------------------+
| keepenergy     | double      | -1.0    | Cut-off energy for truncation                        |
+----------------+-------------+---------+------------------------------------------------------+
| keepmin        | size_t      | 0       | Minimum number of states to keep at each step        |
+----------------+-------------+---------+------------------------------------------------------+
| T              | double      | 0.001   | Temperature, k_B T/D,                                |
+----------------+-------------+---------+------------------------------------------------------+
| ops            | std::string | ""      | Operators to be calculated                           |
+----------------+-------------+---------+------------------------------------------------------+
| specs          | std::string | ""      | Spectral functions (singlet ops) to compute          |
+----------------+-------------+---------+------------------------------------------------------+
| specd          | std::string | ""      | Spectral functions (doublet ops) to compute          |
+----------------+-------------+---------+------------------------------------------------------+
| spect          | std::string | ""      | Spectral functions (triplet ops) to compute          |
+----------------+-------------+---------+------------------------------------------------------+
| specq          | std::string | ""      | Spectral functions (quadruplet ops) to compute       |
+----------------+-------------+---------+------------------------------------------------------+
| specot         | std::string | ""      | Spectral functions (orbital triplet ops) to compute  |
+----------------+-------------+---------+------------------------------------------------------+
| specgt         | std::string | ""      | Conductance curves to compute                        |
+----------------+-------------+---------+------------------------------------------------------+
| speci1t        | std::string | ""      | I_1 curves to compute                                |
+----------------+-------------+---------+------------------------------------------------------+
| speci2t        | std::string | ""      | I_2 curves to compute                                |
+----------------+-------------+---------+------------------------------------------------------+
| specchit       | std::string | ""      | Susceptibilities to compute                          |
+----------------+-------------+---------+------------------------------------------------------+
| specv3         | std::string | ""      | 3-leg vertex functions to compute?                   |
+----------------+-------------+---------+------------------------------------------------------+
| v3mm           | bool        | false   | Compute 3-leg vertex on matsubara/matsubara axis?    |
+----------------+-------------+---------+------------------------------------------------------+
| dmnrg          | bool        | false   | Perform DMNRG (density-matrix NRG) calculation       |
+----------------+-------------+---------+------------------------------------------------------+
| cfs            | bool        | false   | Perform CFS (complete Fock space) calculation        |
+----------------+-------------+---------+------------------------------------------------------+
| fdm            | bool        | false   | Perform FDM (full-density-matrix) calculation        |
+----------------+-------------+---------+------------------------------------------------------+
| fdmexpv        | bool        | false   | Calculate expectation values using FDM               |
+----------------+-------------+---------+------------------------------------------------------+
| dmnrgmats      | bool        | false   | DMNRG calculation on Matsubara axis                  |
+----------------+-------------+---------+------------------------------------------------------+
| fdmmats        | bool        | false   | FDM calculation on Matsubara axis                    |
+----------------+-------------+---------+------------------------------------------------------+
| mats           | size_t      | 100     | Number of Matsubara points to collect                |
+----------------+-------------+---------+------------------------------------------------------+
| alpha          | double      | 0.3     | Width of logarithmic gaussian                        |
+----------------+-------------+---------+------------------------------------------------------+
| gamma          | double      | 0.2     | Parameter for Gaussian convolution step              |
+----------------+-------------+---------+------------------------------------------------------+
| post_process   | bool        | true    | Perform post processing                              |
+----------------+-------------+---------+------------------------------------------------------+
""")

c.add_method("""void set_nrg_params (**nrgljubljana_interface::nrg_params_t)""",
             doc = r"""



+---------------------+-------------+-----------+------------------------------------------------------------+
| Parameter Name      | Type        | Default   | Documentation                                              |
+=====================+=============+===========+============================================================+
| discretization      | std::string | "Z"       | Discretization scheme                                      |
+---------------------+-------------+-----------+------------------------------------------------------------+
| z                   | double      | 1.0       | Parameter z in the logarithmic discretization              |
+---------------------+-------------+-----------+------------------------------------------------------------+
| polarized           | bool        | false     | Spin-polarized Wilson chain                                |
+---------------------+-------------+-----------+------------------------------------------------------------+
| pol2x2              | bool        | false     | 2x2 spin structure in Wilson chain                         |
+---------------------+-------------+-----------+------------------------------------------------------------+
| rungs               | bool        | false     | Channel-mixing terms in Wilson chain                       |
+---------------------+-------------+-----------+------------------------------------------------------------+
| tri                 | std::string | "old"     | Tridiagonalisation approach                                |
+---------------------+-------------+-----------+------------------------------------------------------------+
| preccpp             | size_t      | 2000      | Precision for tridiagonalisation                           |
+---------------------+-------------+-----------+------------------------------------------------------------+
| diag                | std::string | "default" | Eigensolver routine (dsyev|dsyevr|zheev|zheevr|default)    |
+---------------------+-------------+-----------+------------------------------------------------------------+
| diagratio           | double      | 1.0       | Ratio of eigenstates computed in partial diagonalisation   |
+---------------------+-------------+-----------+------------------------------------------------------------+
| dsyevrlimit         | size_t      | 100       | Minimal matrix size for dsyevr                             |
+---------------------+-------------+-----------+------------------------------------------------------------+
| zheevrlimit         | size_t      | 100       | Minimal matrix size for zheevr                             |
+---------------------+-------------+-----------+------------------------------------------------------------+
| restart             | bool        | true      | Restart calculation to achieve truncation goal?            |
+---------------------+-------------+-----------+------------------------------------------------------------+
| restartfactor       | double      | 2.0       | Rescale factor for restart=true                            |
+---------------------+-------------+-----------+------------------------------------------------------------+
| safeguard           | double      | 1e-5      | Additional states to keep in case of a near degeneracy     |
+---------------------+-------------+-----------+------------------------------------------------------------+
| safeguardmax        | size_t      | 200       | Maximal number of additional states                        |
+---------------------+-------------+-----------+------------------------------------------------------------+
| fixeps              | double      | 1e-15     | Threshold value for eigenvalue splitting corrections       |
+---------------------+-------------+-----------+------------------------------------------------------------+
| betabar             | double      | 1.0       | Parameter \bar{\beta} for thermodynamics                   |
+---------------------+-------------+-----------+------------------------------------------------------------+
| gtp                 | double      | 0.7       | Parameter p for G(T) calculations                          |
+---------------------+-------------+-----------+------------------------------------------------------------+
| chitp               | double      | 1.0       | Parameter p for chi(T) calculations                        |
+---------------------+-------------+-----------+------------------------------------------------------------+
| finite              | bool        | false     | Perform Costi-Hewson-Zlatic finite-T calculation           |
+---------------------+-------------+-----------+------------------------------------------------------------+
| cfsgt               | bool        | false     | CFS greater correlation function                           |
+---------------------+-------------+-----------+------------------------------------------------------------+
| cfsls               | bool        | false     | CFS lesser correlation function                            |
+---------------------+-------------+-----------+------------------------------------------------------------+
| fdmgt               | bool        | false     | FDM greater correlation function?                          |
+---------------------+-------------+-----------+------------------------------------------------------------+
| fdmls               | bool        | false     | FDM lesser correlation function?                           |
+---------------------+-------------+-----------+------------------------------------------------------------+
| fdmexpvn            | size_t      | 0         | Iteration where we evaluate expv                           |
+---------------------+-------------+-----------+------------------------------------------------------------+
| finitemats          | bool        | false     | T>0 calculation on Matsubara axis                          |
+---------------------+-------------+-----------+------------------------------------------------------------+
| dm                  | bool        | false     | Compute density matrixes?                                  |
+---------------------+-------------+-----------+------------------------------------------------------------+
| broaden_max         | double      | 10        | Broadening mesh maximum frequency                          |
+---------------------+-------------+-----------+------------------------------------------------------------+
| broaden_min         | double      | -99.      | Broadening mesh minimum frequency                          |
+---------------------+-------------+-----------+------------------------------------------------------------+
| broaden_min_ratio   | double      | 3.0       | Auto-tune broaden_min parameter                            |
+---------------------+-------------+-----------+------------------------------------------------------------+
| broaden_ratio       | double      | 1.05      | Common ration of the geometric sequence                    |
+---------------------+-------------+-----------+------------------------------------------------------------+
| omega0              | double      | -1.0      | Smallest energy scale in the problem                       |
+---------------------+-------------+-----------+------------------------------------------------------------+
| omega0_ratio        | double      | 1.0       | omega0 = omega0_ratio x T                                  |
+---------------------+-------------+-----------+------------------------------------------------------------+
| diagth              | int         | 1         | Diagonalisation threads                                    |
+---------------------+-------------+-----------+------------------------------------------------------------+
| substeps            | bool        | false     | Interleaved diagonalization                                |
+---------------------+-------------+-----------+------------------------------------------------------------+
| strategy            | std::string | "kept"    | Recalculation strategy                                     |
+---------------------+-------------+-----------+------------------------------------------------------------+
| Ninit               | size_t      | 0         | Initial Wilson chain ops                                   |
+---------------------+-------------+-----------+------------------------------------------------------------+
| reim                | bool        | false     | Output imaginary parts of correlators?                     |
+---------------------+-------------+-----------+------------------------------------------------------------+
| dumpannotated       | size_t      | 0         | Number of eigenvalues to dump                              |
+---------------------+-------------+-----------+------------------------------------------------------------+
| dumpabs             | bool        | false     | Dump in terms of absolute energies                         |
+---------------------+-------------+-----------+------------------------------------------------------------+
| dumpscaled          | bool        | true      | Dump using omega_N energy units                            |
+---------------------+-------------+-----------+------------------------------------------------------------+
| dumpprecision       | size_t      | 8         | Dump with # digits of precision                            |
+---------------------+-------------+-----------+------------------------------------------------------------+
| dumpgroups          | bool        | true      | Dump by grouping degenerate states                         |
+---------------------+-------------+-----------+------------------------------------------------------------+
| grouptol            | double      | 1e-6      | Energy tolerance for considering two states as degenerate  |
+---------------------+-------------+-----------+------------------------------------------------------------+
| dumpdiagonal        | size_t      | 0         | Dump diagonal matrix elements                              |
+---------------------+-------------+-----------+------------------------------------------------------------+
| savebins            | bool        | false     | Save binned (unbroadened) data                             |
+---------------------+-------------+-----------+------------------------------------------------------------+
| broaden             | bool        | true      | Enable broadening of spectra                               |
+---------------------+-------------+-----------+------------------------------------------------------------+
| emin                | double      | -1.0      | Lower binning limit                                        |
+---------------------+-------------+-----------+------------------------------------------------------------+
| emax                | double      | -1.0      | Upper binning limit                                        |
+---------------------+-------------+-----------+------------------------------------------------------------+
| bins                | size_t      | 1000      | bins/decade for spectral data                              |
+---------------------+-------------+-----------+------------------------------------------------------------+
| accumulation        | double      | 0.0       | Shift of the accumulation points for binning               |
+---------------------+-------------+-----------+------------------------------------------------------------+
| linstep             | double      | 0         | Bin width for linear mesh                                  |
+---------------------+-------------+-----------+------------------------------------------------------------+
| DISCARD_TRIM        | double      | 1e-16     | Peak clipping at the end of the run                        |
+---------------------+-------------+-----------+------------------------------------------------------------+
| DISCARD_IMMEDIATELY | double      | 1e-16     | Peak clipping on the fly                                   |
+---------------------+-------------+-----------+------------------------------------------------------------+
| goodE               | double      | 2.0       | Energy window parameter for patching                       |
+---------------------+-------------+-----------+------------------------------------------------------------+
| NN1                 | bool        | false     | Do N/N+1 patching?                                         |
+---------------------+-------------+-----------+------------------------------------------------------------+
| NN2even             | bool        | true      | Use even iterations in N/N+2 patching                      |
+---------------------+-------------+-----------+------------------------------------------------------------+
| NN2avg              | bool        | false     | Average over even and odd N/N+2 spectra                    |
+---------------------+-------------+-----------+------------------------------------------------------------+
| NNtanh              | double      | 0.0       | a in tanh[a(x-0.5)] window function                        |
+---------------------+-------------+-----------+------------------------------------------------------------+
| width_td            | size_t      | 16        | Widht of columns in 'td'                                   |
+---------------------+-------------+-----------+------------------------------------------------------------+
| width_custom        | size_t      | 16        | Width of columns in 'custom'                               |
+---------------------+-------------+-----------+------------------------------------------------------------+
| prec_td             | size_t      | 10        | Precision of columns in 'td'                               |
+---------------------+-------------+-----------+------------------------------------------------------------+
| prec_custom         | size_t      | 10        | Precision of columns in 'custom'                           |
+---------------------+-------------+-----------+------------------------------------------------------------+
| prec_xy             | size_t      | 10        | Precision of spectral function output                      |
+---------------------+-------------+-----------+------------------------------------------------------------+
| resume              | bool        | false     | Attempt restart?                                           |
+---------------------+-------------+-----------+------------------------------------------------------------+
| log                 | std::string | ""        | list of tokens to define what to log                       |
+---------------------+-------------+-----------+------------------------------------------------------------+
| logall              | bool        | false     | Log everything                                             |
+---------------------+-------------+-----------+------------------------------------------------------------+
| done                | bool        | true      | Create DONE file?                                          |
+---------------------+-------------+-----------+------------------------------------------------------------+
| calc0               | bool        | true      | Perform calculations at 0-th iteration?                    |
+---------------------+-------------+-----------+------------------------------------------------------------+
| lastall             | bool        | false     | Keep all states in the last iteratio for DMNRG             |
+---------------------+-------------+-----------+------------------------------------------------------------+
| lastalloverride     | bool        | false     | Override automatic lastall setting                         |
+---------------------+-------------+-----------+------------------------------------------------------------+
| dumpsubspaces       | bool        | false     | Save detailed subspace info                                |
+---------------------+-------------+-----------+------------------------------------------------------------+
| dump_f              | bool        | false     | Dump <f> matrix elements                                   |
+---------------------+-------------+-----------+------------------------------------------------------------+
| dumpenergies        | bool        | false     | Dump (all) energies to file?                               |
+---------------------+-------------+-----------+------------------------------------------------------------+
| logenumber          | size_t      | 10        | # of eigenvalues to show for log=e                         |
+---------------------+-------------+-----------+------------------------------------------------------------+
| stopafter           | std::string | ""        | Stop calculation at some point?                            |
+---------------------+-------------+-----------+------------------------------------------------------------+
| forcestop           | int         | -1        | Stop iteration?                                            |
+---------------------+-------------+-----------+------------------------------------------------------------+
| removefiles         | bool        | true      | Remove temporary data files?                               |
+---------------------+-------------+-----------+------------------------------------------------------------+
| noimag              | bool        | true      | Do not output imaginary parts of expvs                     |
+---------------------+-------------+-----------+------------------------------------------------------------+
| checksumrules       | bool        | false     | Check operator sumrules                                    |
+---------------------+-------------+-----------+------------------------------------------------------------+
| checkdiag           | bool        | false     | Test diag results                                          |
+---------------------+-------------+-----------+------------------------------------------------------------+
| checkrho            | bool        | false     | Test tr(rho)=1                                             |
+---------------------+-------------+-----------+------------------------------------------------------------+
""")

c.add_method("""std::string hdf5_scheme ()""",
             is_static = True,
             doc = r"""""")

module.add_class(c)


# Converter for solve_params_t
c = converter_(
        c_type = "nrgljubljana_interface::solve_params_t",
        doc = r"""The parameters for the solve function""",
)
c.add_member(c_name = "Lambda",
             c_type = "double",
             initializer = """ 2.0 """,
             doc = r"""Logarithmic discretization parameter""")

c.add_member(c_name = "Nz",
             c_type = "int",
             initializer = """ 1 """,
             doc = r"""Number of discretization meshes""")

c.add_member(c_name = "xmax",
             c_type = "double",
             initializer = """ 20 """,
             doc = r"""Largest x in discretization ODE solver""")

c.add_member(c_name = "Tmin",
             c_type = "double",
             initializer = """ 1e-4 """,
             doc = r"""Lowest scale on the Wilson chain""")

c.add_member(c_name = "keep",
             c_type = "size_t",
             initializer = """ 100 """,
             doc = r"""Maximum number of states to keep at each step""")

c.add_member(c_name = "keepenergy",
             c_type = "double",
             initializer = """ -1.0 """,
             doc = r"""Cut-off energy for truncation""")

c.add_member(c_name = "keepmin",
             c_type = "size_t",
             initializer = """ 0 """,
             doc = r"""Minimum number of states to keep at each step""")

c.add_member(c_name = "T",
             c_type = "double",
             initializer = """ 0.001 """,
             doc = r"""Temperature, k_B T/D,""")

c.add_member(c_name = "ops",
             c_type = "std::string",
             initializer = """ "" """,
             doc = r"""Operators to be calculated""")

c.add_member(c_name = "specs",
             c_type = "std::string",
             initializer = """ "" """,
             doc = r"""Spectral functions (singlet ops) to compute""")

c.add_member(c_name = "specd",
             c_type = "std::string",
             initializer = """ "" """,
             doc = r"""Spectral functions (doublet ops) to compute""")

c.add_member(c_name = "spect",
             c_type = "std::string",
             initializer = """ "" """,
             doc = r"""Spectral functions (triplet ops) to compute""")

c.add_member(c_name = "specq",
             c_type = "std::string",
             initializer = """ "" """,
             doc = r"""Spectral functions (quadruplet ops) to compute""")

c.add_member(c_name = "specot",
             c_type = "std::string",
             initializer = """ "" """,
             doc = r"""Spectral functions (orbital triplet ops) to compute""")

c.add_member(c_name = "specgt",
             c_type = "std::string",
             initializer = """ "" """,
             doc = r"""Conductance curves to compute""")

c.add_member(c_name = "speci1t",
             c_type = "std::string",
             initializer = """ "" """,
             doc = r"""I_1 curves to compute""")

c.add_member(c_name = "speci2t",
             c_type = "std::string",
             initializer = """ "" """,
             doc = r"""I_2 curves to compute""")

c.add_member(c_name = "specchit",
             c_type = "std::string",
             initializer = """ "" """,
             doc = r"""Susceptibilities to compute""")

c.add_member(c_name = "specv3",
             c_type = "std::string",
             initializer = """ "" """,
             doc = r"""3-leg vertex functions to compute?""")

c.add_member(c_name = "v3mm",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Compute 3-leg vertex on matsubara/matsubara axis?""")

c.add_member(c_name = "dmnrg",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Perform DMNRG (density-matrix NRG) calculation""")

c.add_member(c_name = "cfs",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Perform CFS (complete Fock space) calculation""")

c.add_member(c_name = "fdm",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Perform FDM (full-density-matrix) calculation""")

c.add_member(c_name = "fdmexpv",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Calculate expectation values using FDM""")

c.add_member(c_name = "dmnrgmats",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""DMNRG calculation on Matsubara axis""")

c.add_member(c_name = "fdmmats",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""FDM calculation on Matsubara axis""")

c.add_member(c_name = "mats",
             c_type = "size_t",
             initializer = """ 100 """,
             doc = r"""Number of Matsubara points to collect""")

c.add_member(c_name = "alpha",
             c_type = "double",
             initializer = """ 0.3 """,
             doc = r"""Width of logarithmic gaussian""")

c.add_member(c_name = "gamma",
             c_type = "double",
             initializer = """ 0.2 """,
             doc = r"""Parameter for Gaussian convolution step""")

c.add_member(c_name = "post_process",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Perform post processing""")

module.add_converter(c)

# Converter for nrg_params_t
c = converter_(
        c_type = "nrgljubljana_interface::nrg_params_t",
        doc = r"""NRG low-level parameters""",
)
c.add_member(c_name = "discretization",
             c_type = "std::string",
             initializer = """ "Z" """,
             doc = r"""Discretization scheme""")

c.add_member(c_name = "z",
             c_type = "double",
             initializer = """ 1.0 """,
             doc = r"""Parameter z in the logarithmic discretization""")

c.add_member(c_name = "polarized",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Spin-polarized Wilson chain""")

c.add_member(c_name = "pol2x2",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""2x2 spin structure in Wilson chain""")

c.add_member(c_name = "rungs",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Channel-mixing terms in Wilson chain""")

c.add_member(c_name = "tri",
             c_type = "std::string",
             initializer = """ "old" """,
             doc = r"""Tridiagonalisation approach""")

c.add_member(c_name = "preccpp",
             c_type = "size_t",
             initializer = """ 2000 """,
             doc = r"""Precision for tridiagonalisation""")

c.add_member(c_name = "diag",
             c_type = "std::string",
             initializer = """ "default" """,
             doc = r"""Eigensolver routine (dsyev|dsyevr|zheev|zheevr|default)""")

c.add_member(c_name = "diagratio",
             c_type = "double",
             initializer = """ 1.0 """,
             doc = r"""Ratio of eigenstates computed in partial diagonalisation""")

c.add_member(c_name = "dsyevrlimit",
             c_type = "size_t",
             initializer = """ 100 """,
             doc = r"""Minimal matrix size for dsyevr""")

c.add_member(c_name = "zheevrlimit",
             c_type = "size_t",
             initializer = """ 100 """,
             doc = r"""Minimal matrix size for zheevr""")

c.add_member(c_name = "restart",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Restart calculation to achieve truncation goal?""")

c.add_member(c_name = "restartfactor",
             c_type = "double",
             initializer = """ 2.0 """,
             doc = r"""Rescale factor for restart=true""")

c.add_member(c_name = "safeguard",
             c_type = "double",
             initializer = """ 1e-5 """,
             doc = r"""Additional states to keep in case of a near degeneracy""")

c.add_member(c_name = "safeguardmax",
             c_type = "size_t",
             initializer = """ 200 """,
             doc = r"""Maximal number of additional states""")

c.add_member(c_name = "fixeps",
             c_type = "double",
             initializer = """ 1e-15 """,
             doc = r"""Threshold value for eigenvalue splitting corrections""")

c.add_member(c_name = "betabar",
             c_type = "double",
             initializer = """ 1.0 """,
             doc = r"""Parameter \bar{\beta} for thermodynamics""")

c.add_member(c_name = "gtp",
             c_type = "double",
             initializer = """ 0.7 """,
             doc = r"""Parameter p for G(T) calculations""")

c.add_member(c_name = "chitp",
             c_type = "double",
             initializer = """ 1.0 """,
             doc = r"""Parameter p for chi(T) calculations""")

c.add_member(c_name = "finite",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Perform Costi-Hewson-Zlatic finite-T calculation""")

c.add_member(c_name = "cfsgt",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""CFS greater correlation function""")

c.add_member(c_name = "cfsls",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""CFS lesser correlation function""")

c.add_member(c_name = "fdmgt",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""FDM greater correlation function?""")

c.add_member(c_name = "fdmls",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""FDM lesser correlation function?""")

c.add_member(c_name = "fdmexpvn",
             c_type = "size_t",
             initializer = """ 0 """,
             doc = r"""Iteration where we evaluate expv""")

c.add_member(c_name = "finitemats",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""T>0 calculation on Matsubara axis""")

c.add_member(c_name = "dm",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Compute density matrixes?""")

c.add_member(c_name = "broaden_max",
             c_type = "double",
             initializer = """ 10 """,
             doc = r"""Broadening mesh maximum frequency""")

c.add_member(c_name = "broaden_min",
             c_type = "double",
             initializer = """ -99. """,
             doc = r"""Broadening mesh minimum frequency""")

c.add_member(c_name = "broaden_min_ratio",
             c_type = "double",
             initializer = """ 3.0 """,
             doc = r"""Auto-tune broaden_min parameter""")

c.add_member(c_name = "broaden_ratio",
             c_type = "double",
             initializer = """ 1.05 """,
             doc = r"""Common ration of the geometric sequence""")

c.add_member(c_name = "omega0",
             c_type = "double",
             initializer = """ -1.0 """,
             doc = r"""Smallest energy scale in the problem""")

c.add_member(c_name = "omega0_ratio",
             c_type = "double",
             initializer = """ 1.0 """,
             doc = r"""omega0 = omega0_ratio x T""")

c.add_member(c_name = "diagth",
             c_type = "int",
             initializer = """ 1 """,
             doc = r"""Diagonalisation threads""")

c.add_member(c_name = "substeps",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Interleaved diagonalization""")

c.add_member(c_name = "strategy",
             c_type = "std::string",
             initializer = """ "kept" """,
             doc = r"""Recalculation strategy""")

c.add_member(c_name = "Ninit",
             c_type = "size_t",
             initializer = """ 0 """,
             doc = r"""Initial Wilson chain ops""")

c.add_member(c_name = "reim",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Output imaginary parts of correlators?""")

c.add_member(c_name = "dumpannotated",
             c_type = "size_t",
             initializer = """ 0 """,
             doc = r"""Number of eigenvalues to dump""")

c.add_member(c_name = "dumpabs",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Dump in terms of absolute energies""")

c.add_member(c_name = "dumpscaled",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Dump using omega_N energy units""")

c.add_member(c_name = "dumpprecision",
             c_type = "size_t",
             initializer = """ 8 """,
             doc = r"""Dump with # digits of precision""")

c.add_member(c_name = "dumpgroups",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Dump by grouping degenerate states""")

c.add_member(c_name = "grouptol",
             c_type = "double",
             initializer = """ 1e-6 """,
             doc = r"""Energy tolerance for considering two states as degenerate""")

c.add_member(c_name = "dumpdiagonal",
             c_type = "size_t",
             initializer = """ 0 """,
             doc = r"""Dump diagonal matrix elements""")

c.add_member(c_name = "savebins",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Save binned (unbroadened) data""")

c.add_member(c_name = "broaden",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Enable broadening of spectra""")

c.add_member(c_name = "emin",
             c_type = "double",
             initializer = """ -1.0 """,
             doc = r"""Lower binning limit""")

c.add_member(c_name = "emax",
             c_type = "double",
             initializer = """ -1.0 """,
             doc = r"""Upper binning limit""")

c.add_member(c_name = "bins",
             c_type = "size_t",
             initializer = """ 1000 """,
             doc = r"""bins/decade for spectral data""")

c.add_member(c_name = "accumulation",
             c_type = "double",
             initializer = """ 0.0 """,
             doc = r"""Shift of the accumulation points for binning""")

c.add_member(c_name = "linstep",
             c_type = "double",
             initializer = """ 0 """,
             doc = r"""Bin width for linear mesh""")

c.add_member(c_name = "DISCARD_TRIM",
             c_type = "double",
             initializer = """ 1e-16 """,
             doc = r"""Peak clipping at the end of the run""")

c.add_member(c_name = "DISCARD_IMMEDIATELY",
             c_type = "double",
             initializer = """ 1e-16 """,
             doc = r"""Peak clipping on the fly""")

c.add_member(c_name = "goodE",
             c_type = "double",
             initializer = """ 2.0 """,
             doc = r"""Energy window parameter for patching""")

c.add_member(c_name = "NN1",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Do N/N+1 patching?""")

c.add_member(c_name = "NN2even",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Use even iterations in N/N+2 patching""")

c.add_member(c_name = "NN2avg",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Average over even and odd N/N+2 spectra""")

c.add_member(c_name = "NNtanh",
             c_type = "double",
             initializer = """ 0.0 """,
             doc = r"""a in tanh[a(x-0.5)] window function""")

c.add_member(c_name = "width_td",
             c_type = "size_t",
             initializer = """ 16 """,
             doc = r"""Widht of columns in 'td'""")

c.add_member(c_name = "width_custom",
             c_type = "size_t",
             initializer = """ 16 """,
             doc = r"""Width of columns in 'custom'""")

c.add_member(c_name = "prec_td",
             c_type = "size_t",
             initializer = """ 10 """,
             doc = r"""Precision of columns in 'td'""")

c.add_member(c_name = "prec_custom",
             c_type = "size_t",
             initializer = """ 10 """,
             doc = r"""Precision of columns in 'custom'""")

c.add_member(c_name = "prec_xy",
             c_type = "size_t",
             initializer = """ 10 """,
             doc = r"""Precision of spectral function output""")

c.add_member(c_name = "resume",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Attempt restart?""")

c.add_member(c_name = "log",
             c_type = "std::string",
             initializer = """ "" """,
             doc = r"""list of tokens to define what to log""")

c.add_member(c_name = "logall",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Log everything""")

c.add_member(c_name = "done",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Create DONE file?""")

c.add_member(c_name = "calc0",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Perform calculations at 0-th iteration?""")

c.add_member(c_name = "lastall",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Keep all states in the last iteratio for DMNRG""")

c.add_member(c_name = "lastalloverride",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Override automatic lastall setting""")

c.add_member(c_name = "dumpsubspaces",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Save detailed subspace info""")

c.add_member(c_name = "dump_f",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Dump <f> matrix elements""")

c.add_member(c_name = "dumpenergies",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Dump (all) energies to file?""")

c.add_member(c_name = "logenumber",
             c_type = "size_t",
             initializer = """ 10 """,
             doc = r"""# of eigenvalues to show for log=e""")

c.add_member(c_name = "stopafter",
             c_type = "std::string",
             initializer = """ "" """,
             doc = r"""Stop calculation at some point?""")

c.add_member(c_name = "forcestop",
             c_type = "int",
             initializer = """ -1 """,
             doc = r"""Stop iteration?""")

c.add_member(c_name = "removefiles",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Remove temporary data files?""")

c.add_member(c_name = "noimag",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Do not output imaginary parts of expvs""")

c.add_member(c_name = "checksumrules",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Check operator sumrules""")

c.add_member(c_name = "checkdiag",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Test diag results""")

c.add_member(c_name = "checkrho",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Test tr(rho)=1""")

module.add_converter(c)

# Converter for constr_params_t
c = converter_(
        c_type = "nrgljubljana_interface::constr_params_t",
        doc = r"""The parameters for the solver construction""",
)
c.add_member(c_name = "bandrescale",
             c_type = "double",
             initializer = """ 1.0 """,
             doc = r"""Band rescaling factor""")

c.add_member(c_name = "problem",
             c_type = "std::string",
             initializer = """ "SIAM" """,
             doc = r"""Model considered (templated)""")

c.add_member(c_name = "mesh_max",
             c_type = "double",
             initializer = """ 10 """,
             doc = r"""Mesh maximum frequency""")

c.add_member(c_name = "mesh_min",
             c_type = "double",
             initializer = """ 1e-4 """,
             doc = r"""Mesh minimum frequency""")

c.add_member(c_name = "mesh_ratio",
             c_type = "double",
             initializer = """ 1.05 """,
             doc = r"""Common ratio of the geometric sequence""")

module.add_converter(c)


module.generate_code()