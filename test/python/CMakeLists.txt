# Copy h5 files to binary dir
file(GLOB_RECURSE all_h5_ref_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.ref.h5)
foreach(file ${all_h5_ref_files})
  configure_file(${file} ${file} COPYONLY)
endforeach()

# List of all tests
set(all_tests test1_siam)

foreach(test ${all_tests})
  set(test_name ${test})
  set(test_dir ${test})
  set(test_script test.py)
  configure_file(${test_dir}/${test_script} ${test_dir}/${test_script} COPYONLY)
  add_test(NAME Py_${test_name} COMMAND ${TRIQS_PYTHON_INTERPRETER} ${test_script} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
  set_property(TEST Py_${test_name} APPEND PROPERTY ENVIRONMENT 
     PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH} 
     NRGIF_TEMPLATE_DIR=${PROJECT_SOURCE_DIR}/templates ${SANITIZER_RT_PRELOAD}
     PATH=${CMAKE_BINARY_DIR}/nrgljubljana/tools:$ENV{PATH}
  )
endforeach()
