add_subdirectory(cmake)

if(NOT CMAKE_INSTALL_PREFIX STREQUAL TRIQS_ROOT)

  configure_file(nrgljubljana_interface.modulefile.in nrgljubljana_interface.modulefile @ONLY)
  configure_file(nrgljubljana_interfacevars.sh.in nrgljubljana_interfacevars.sh @ONLY)

  install(
    FILES
      ${CMAKE_CURRENT_BINARY_DIR}/nrgljubljana_interface.modulefile
      ${CMAKE_CURRENT_BINARY_DIR}/nrgljubljana_interfacevars.sh
    DESTINATION share
  )
 
  message(STATUS "***************************************************************")
  message(STATUS "* Custom install Location. Use:                                ")
  message(STATUS "*                                                              ")
  message(STATUS "*   source ${CMAKE_INSTALL_PREFIX}/share/nrgljubljana_interfacevars.sh      ")
  message(STATUS "*                                                              ")
  message(STATUS "* to set up the environment variables                          ")
  message(STATUS "***************************************************************")

endif()
