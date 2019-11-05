# Generated automatically using the command :
# c++2py ../../c++/nrgljubljana_interface/post_process.hpp --members_read_only -N nrgljubljana_interface -a nrgljubljana_interface -m post_process -o post_process -C pytriqs --moduledoc="The nrgljubljana_interface postprocess functionality" --cxxflags="-std=c++17" --target_file_only
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "post_process", doc = r"The nrgljubljana_interface postprocess functionality", app_name = "nrgljubljana_interface")

# Imports

# Add here all includes
module.add_include("nrgljubljana_interface/post_process.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""

using namespace nrgljubljana_interface;
""")




module.generate_code()