from cpp2py.wrap_generator import *
import re

module = module_(full_name = "mesh_refreq_pts", doc = "The refreq_pts mesh", app_name="nrgljubljana_interface")

module.add_imports(*['triqs.gf']) 

module.add_include("<triqs/gfs.hpp>")
module.add_include("<triqs/mesh/refreq_pts.hpp>")

module.add_include("<cpp2py/converters/string.hpp>")
module.add_include("<cpp2py/converters/vector.hpp>")
module.add_include("<cpp2py/converters/function.hpp>")
module.add_include("<cpp2py/converters/optional.hpp>")

module.add_include("<triqs/cpp2py_converters.hpp>")

module.add_using("namespace nda")
module.add_using("namespace triqs::gfs")
module.add_using("namespace triqs::mesh")
module.add_using("std::array")
module.add_preamble("""
""")

########################
##   Mesh refreq_pts
########################

m = class_( py_type = "MeshReFreqPts",
        c_type = "refreq_pts",
        c_type_absolute = "triqs::mesh::refreq_pts",
        hdf5 = True,
        serializable= "h5",
        is_printable= True,
        comparisons = "== !="
       )

m.add_method(f"refreq_pts::data_index_t to_data_index(refreq_pts::index_t index)", doc = "Function to convert an index to a data index")
m.add_method(f"refreq_pts::index_t to_index(refreq_pts::data_index_t data_index)", doc = "Function to convert a data index to an index")
m.add_getitem(signature = f"refreq_pts::mesh_point_t operator[](refreq_pts::data_index_t data_index)", doc = "Get a mesh-point given the data index")
m.add_call(signature = f"refreq_pts::mesh_point_t operator()(refreq_pts::index_t index)", calling_pattern = " auto result = self_c(index)", doc = "Get a mesh-point given the index")
m.add_len(calling_pattern = "int result = self_c.size()", doc = "Size of the mesh")
m.add_property(name = "mesh_hash",
           getter = cfunction(calling_pattern="uint64_t result = self_c.mesh_hash()",
           signature = "uint64_t()",
           doc = "The hash encoding the mesh configuration"))

m.add_iterator()

m.add_method("PyObject * values()",
             calling_pattern = """
                static auto cls = pyref::get_class("triqs.gf", "MeshValueGenerator", /* raise_exception */ true);
                pyref args = PyTuple_Pack(1, self);
                auto result = PyObject_CallObject(cls, args);
             """, doc = "A numpy array of all the values of the mesh points")
m.add_method(f"refreq_pts::value_t to_value(refreq_pts::index_t index)", doc = "index -> value")

m.add_method_copy()
m.add_method_copy_from()

m.add_constructor(signature = "(std::vector<refreq_pts::value_t> pts)")

module.add_class(m)

##   Code generation
module.generate_code()
