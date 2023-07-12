file(GLOB_RECURSE sources *.cpp)
add_library(${PROJECT_NAME}_c ${sources})
add_library(${PROJECT_NAME}::${PROJECT_NAME}_c ALIAS ${PROJECT_NAME}_c)

# Link against triqs and enable warnings
target_link_libraries(${PROJECT_NAME}_c PUBLIC triqs PRIVATE $<BUILD_INTERFACE:${PROJECT_NAME}_warnings>)

# Configure target and compilation
set_target_properties(${PROJECT_NAME}_c PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
)
target_include_directories(${PROJECT_NAME}_c PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++>)
target_include_directories(${PROJECT_NAME}_c SYSTEM INTERFACE $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>)
target_compile_definitions(${PROJECT_NAME}_c PUBLIC
				NRGIF_TEMPLATE_DIR="${CMAKE_INSTALL_PREFIX}/templates"
				NRGLJUBLJANA_INTERFACE_GIT_HASH=${PROJECT_GIT_HASH}
				TRIQS_GIT_HASH=${TRIQS_GIT_HASH}
				$<$<CONFIG:Debug>:NRGLJUBLJANA_INTERFACE_DEBUG>
				$<$<CONFIG:Debug>:TRIQS_DEBUG>
				$<$<CONFIG:Debug>:TRIQS_ARRAYS_ENFORCE_BOUNDCHECK>
			  )

# ############
# Dependencies

# Link against NRG
target_link_libraries(${PROJECT_NAME}_c PRIVATE nrgljubljana_c)

# Link against other libraries
target_link_libraries(${PROJECT_NAME}_c PUBLIC gsl blas_lapack)
# Note: we link against blas_lapack to prevent problem reported in https://software.intel.com/en-us/forums/intel-math-kernel-library/topic/759670

# Add tool dependencies
add_dependencies(${PROJECT_NAME}_c adapt nrgchain matrix diag unitary broaden kk)

# C++ filesystem library. Currently not used, keep commented out.
#if((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") AND (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0.0"))
#  target_link_libraries(${PROJECT_NAME}_c PRIVATE stdc++fs)
#  message("[filesystem library] Using stdc++fs")
#elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
#  target_link_libraries(${PROJECT_NAME}_c PRIVATE c++experimental c++fs)
#  message("[filesystem library] Using c++experimenal c++fs")
#else()
#  message(FATAL_ERROR "Don't know how to handle filesystem dependency.")
#endif()

# Install library and headers
install(TARGETS ${PROJECT_NAME}_c EXPORT ${PROJECT_NAME}-targets DESTINATION lib)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION include FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h")


# ========= Static Analyzer Checks ==========

option(ANALYZE_SOURCES OFF "Run static analyzer checks if found (clang-tidy, cppcheck)")
if(ANALYZE_SOURCES)

  # Locate static analyzer tools
  find_program(CPPCHECK_EXECUTABLE NAMES "cppcheck" PATHS ENV PATH)
  find_program(CLANG_TIDY_EXECUTABLE NAMES "clang-tidy" PATHS ENV PATH)

  # Run clang-tidy if found
  if(CLANG_TIDY_EXECUTABLE)
    message(STATUS "clang-tidy found: ${CLANG_TIDY_EXECUTABLE}")
    set_target_properties(${PROJECT_NAME}_c PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXECUTABLE}")
  else()
    message(STATUS "clang-tidy not found in $PATH. Please consider installing clang-tidy for additional checks!")
  endif()

  # Run cppcheck if found
  if(CPPCHECK_EXECUTABLE)
    message(STATUS "cppcheck found: ${CPPCHECK_EXECUTABLE}")
    add_custom_command(
	TARGET ${PROJECT_NAME}_c
	COMMAND ${CPPCHECK_EXECUTABLE}
	  --enable=warning,style,performance,portability
	  --std=c++23
	  --template=gcc
	  --verbose
	  --force
	  --quiet
          ${sources}
	WORKING_DIRECTORY
	  ${CMAKE_CURRENT_SOURCE_DIR}
    )
  else()
    message(STATUS "cppcheck not found in $PATH. Please consider installing cppcheck for additional checks!")
  endif()

endif()

# ========= Dynamic Analyzer Checks ==========

option(ASAN OFF "Compile library and executables with LLVM Address Sanitizer")
if(ASAN)
  find_package(sanitizer REQUIRED COMPONENTS asan)
  target_link_libraries(${PROJECT_NAME}_c PUBLIC $<BUILD_INTERFACE:asan>)
endif()

option(UBSAN OFF "Compile library and executables with LLVM Undefined Behavior Sanitizer")
if(UBSAN)
  find_package(sanitizer REQUIRED COMPONENTS ubsan)
  target_link_libraries(${PROJECT_NAME}_c PUBLIC $<BUILD_INTERFACE:ubsan>)
endif()

option(MSAN OFF "Compile library and executables with LLVM Memory Sanitizer")
if(MSAN)
  find_package(sanitizer REQUIRED COMPONENTS msan)
  target_link_libraries(${PROJECT_NAME}_c PUBLIC $<BUILD_INTERFACE:msan>)
endif()

option(TSAN OFF "Compile library and executables with LLVM Thread Sanitizer")
if(TSAN)
  find_package(sanitizer REQUIRED COMPONENTS tsan)
  target_link_libraries(${PROJECT_NAME}_c PUBLIC $<BUILD_INTERFACE:tsan>)
endif()
