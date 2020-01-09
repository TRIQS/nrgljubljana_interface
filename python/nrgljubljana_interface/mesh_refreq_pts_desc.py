from cpp2py.wrap_generator import *
import re

module = module_(full_name = "mesh_refreq_pts", doc = "The refreq_pts mesh", app_name="nrgljubljana_interface")

module.add_imports(*['pytriqs.gf']) 

module.add_include("<triqs/gfs.hpp>")
module.add_include("<triqs/gfs/meshes/refreq_pts.hpp>")

module.add_include("<cpp2py/converters/string.hpp>")
module.add_include("<cpp2py/converters/vector.hpp>")
module.add_include("<cpp2py/converters/function.hpp>")
module.add_include("<cpp2py/converters/optional.hpp>")

module.add_include("<triqs/cpp2py_converters.hpp>")

module.add_using("namespace triqs::arrays")
module.add_using("namespace triqs::gfs")
module.add_using("triqs::utility::mini_vector")
module.add_preamble("""
""")

########################
##   Mesh refreq_pts
########################

m = class_( py_type = "MeshReFreqPts",
        c_type = "gf_mesh<refreq_pts>",
        c_type_absolute = "triqs::gfs::gf_mesh<triqs::gfs::refreq_pts>",
        hdf5 = True,
        serializable= "tuple",
        is_printable= True,
        comparisons = "== !="
       )

m.add_constructor(signature = "(std::vector<gf_mesh<refreq_pts>::domain_pt_t> pts)")
m.add_method("long index_to_linear(long i)", doc = "index -> linear index")
m.add_len(calling_pattern = "int result = self_c.size()", doc = "Size of the mesh")
m.add_iterator()
m.add_method("PyObject * values()",
             calling_pattern = """
                static auto cls = pyref::get_class("pytriqs.gf", "MeshValueGenerator", /* raise_exception */ true);
                pyref args = PyTuple_Pack(1, self);
                auto result = PyObject_CallObject(cls, args);
             """, doc = "A numpy array of all the values of the mesh points")

m.add_method_copy()
m.add_method_copy_from()

# m.add_property(name = "omega_min",
               # getter = cfunction(calling_pattern="double result = self_c.x_min()",
               # signature = "double()",
               # doc = "Inverse temperature"))

# m.add_property(name = "omega_max",
               # getter = cfunction(calling_pattern="double result = self_c.x_max()",
               # signature = "double()",
               # doc = "Inverse temperature"))

# m.add_property(name = "delta",
               # getter = cfunction(calling_pattern="double result = self_c.delta()",
               # signature = "double()",
               # doc = "The mesh-spacing"))

module.add_class(m)

##   Code generation
module.generate_code()
