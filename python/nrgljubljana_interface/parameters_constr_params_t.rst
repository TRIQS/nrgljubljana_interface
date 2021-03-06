+----------------+-------------+--------------------+--------------------------------------------------------------+
| Parameter Name | Type        | Default            | Documentation                                                |
+================+=============+====================+==============================================================+
| templatedir    | std::string | NRGIF_TEMPLATE_DIR | Path to the template library (default to bundled templates)  |
+----------------+-------------+--------------------+--------------------------------------------------------------+
| model          | std::string | "SIAM"             | Model considered (templated)                                 |
+----------------+-------------+--------------------+--------------------------------------------------------------+
| symtype        | std::string | "QS"               | Symmetry                                                     |
+----------------+-------------+--------------------+--------------------------------------------------------------+
| mesh_max       | double      | 10                 | Mesh maximum frequency                                       |
+----------------+-------------+--------------------+--------------------------------------------------------------+
| mesh_min       | double      | 1e-4               | Mesh minimum frequency                                       |
+----------------+-------------+--------------------+--------------------------------------------------------------+
| mesh_ratio     | double      | 1.05               | Common ratio of the geometric sequence                       |
+----------------+-------------+--------------------+--------------------------------------------------------------+
| polarized      | bool        | false              | Spin-polarized Wilson chain                                  |
+----------------+-------------+--------------------+--------------------------------------------------------------+
| pol2x2         | bool        | false              | 2x2 spin structure in Wilson chain                           |
+----------------+-------------+--------------------+--------------------------------------------------------------+
| rungs          | bool        | false              | Channel-mixing terms in Wilson chain                         |
+----------------+-------------+--------------------+--------------------------------------------------------------+
| ops            | std::string | ""                 | Operators to be calculated                                   |
+----------------+-------------+--------------------+--------------------------------------------------------------+
| specs          | std::string | ""                 | Spectral functions (singlet ops) to compute                  |
+----------------+-------------+--------------------+--------------------------------------------------------------+
| specd          | std::string | ""                 | Spectral functions (doublet ops) to compute                  |
+----------------+-------------+--------------------+--------------------------------------------------------------+
| spect          | std::string | ""                 | Spectral functions (triplet ops) to compute                  |
+----------------+-------------+--------------------+--------------------------------------------------------------+
| specq          | std::string | ""                 | Spectral functions (quadruplet ops) to compute               |
+----------------+-------------+--------------------+--------------------------------------------------------------+
| specot         | std::string | ""                 | Spectral functions (orbital triplet ops) to compute          |
+----------------+-------------+--------------------+--------------------------------------------------------------+
| specchit       | std::string | ""                 | Susceptibilities to compute                                  |
+----------------+-------------+--------------------+--------------------------------------------------------------+
| specv3         | std::string | ""                 | 3-leg vertex functions to compute?                           |
+----------------+-------------+--------------------+--------------------------------------------------------------+
| params         | std::string | ""                 | List of model parameters that need to be specified           |
+----------------+-------------+--------------------+--------------------------------------------------------------+