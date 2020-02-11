# Copy h5 files to binary dir
file(GLOB_RECURSE all_h5_ref_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.ref.h5)
foreach(file ${all_h5_ref_files})
  configure_file(${file} ${file} COPYONLY)
endforeach()

# List of all tests
set(all_tests 1_SIAM_QS 2_SIAM_QSZ 3_2orb-UJ_QS 4_mesh_max=10 5_bandrescale 6_solver_h5io
  7_solve_twice 8_solve_last A_refreq_pts)

foreach(test ${all_tests})
  get_filename_component(test_name ${test} NAME_WE)
  get_filename_component(test_dir ${test} DIRECTORY)
  add_test(NAME Py_${test_name} COMMAND ${TRIQS_PYTHON_INTERPRETER} ${CMAKE_CURRENT_SOURCE_DIR}/${test_dir}/${test_name}.py WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
  set_property(TEST Py_${test_name} APPEND PROPERTY ENVIRONMENT
    PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH}
    PATH=${CMAKE_BINARY_DIR}/nrgljubljana/tools:$ENV{PATH}
    NRGIF_TEMPLATE_DIR=${PROJECT_SOURCE_DIR}/templates
    ${SANITIZER_RT_PRELOAD}
  )
endforeach()

# MPI testing
macro(add_mpi_test TestName NumProc)
  get_filename_component(test_name ${TestName} NAME_WE)
  get_filename_component(test_dir ${TestName} DIRECTORY)
  add_test(NAME ${TestName}_np${NumProc}
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NumProc} ${MPIEXEC_PREFLAGS} ${TRIQS_PYTHON_INTERPRETER} ${CMAKE_CURRENT_SOURCE_DIR}/${test_dir}/${test_name}.py 
           ${MPIEXEC_POSTFLAGS} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
  set_property(TEST ${TestName}_np${NumProc} APPEND PROPERTY ENVIRONMENT
      PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH}
      PATH=${CMAKE_BINARY_DIR}/nrgljubljana/tools:$ENV{PATH}
      NRGIF_TEMPLATE_DIR=${PROJECT_SOURCE_DIR}/templates
      ${SANITIZER_RT_PRELOAD}
  )
endmacro()

add_mpi_test(MPI_SIAM_QS 1)
add_mpi_test(MPI_SIAM_QS 2)
add_mpi_test(MPI_SIAM_QS 3)
add_mpi_test(MPI_SIAM_QS 4)
