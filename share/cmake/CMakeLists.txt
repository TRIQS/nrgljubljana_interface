configure_file(nrgljubljana_interface-config.cmake.in nrgljubljana_interface-config.cmake @ONLY)
configure_file(nrgljubljana_interface-config-version.cmake.in nrgljubljana_interface-config-version.cmake @ONLY)
install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/nrgljubljana_interface-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/nrgljubljana_interface-config-version.cmake
  DESTINATION lib/cmake/nrgljubljana_interface
)

install(EXPORT nrgljubljana_interface-targets NAMESPACE nrgljubljana_interface:: DESTINATION lib/cmake/nrgljubljana_interface)
