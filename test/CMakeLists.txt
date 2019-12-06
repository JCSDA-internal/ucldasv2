list( APPEND mymodel_test_input
  testinput/interface.yml
  )


# link the input files for the tests
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
foreach(FILENAME ${mymodel_test_input})
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
    ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
    ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach()


#================================================================================
# Tests of class interfaces
#================================================================================

ecbuild_add_test(
  TARGET  test_mymodel_geometry
  SOURCES executables/TestGeometry.cc
  ARGS    testinput/interface.yml
  LIBS    mymodel )

# ecbuild_add_test(
#   TARGET  test_mymodel_state
#   SOURCES executables/TestState.cc
#   ARGS    testinput/interface.yml
#   LIBS    mymodel )

# ecbuild_add_test(
#   TARGET  test_mymodel_modelauxcontrol
#   SOURCES executables/TestModelAuxControl.cc
#   ARGS    testinput/interface.yml
#   LIBS    mymodel )

# ecbuild_add_test(
#   TARGET  test_mymodel_model
#   SOURCES executables/TestModel.cc
#   ARGS    testinput/interface.yml
#   LIBS    mymodel )

# ecbuild_add_test(
#   TARGET  test_mymodel_increment
#   SOURCES executables/TestIncrement.cc
#   ARGS    testinput/interface.yml
#   LIBS    mymodel )

# ecbuild_add_test(
#   TARGET  test_mymodel_errorcovariance
#   SOURCES executables/TestErrorCovariance.cc
#   ARGS    testinput/interface.yml
#   LIBS    mymodel )

# ecbuild_add_test(
#   TARGET  test_mymodel_modelauxcovariance
#   SOURCES executables/TestModelAuxCovariance.cc
#   ARGS    testinput/interface.yml
#   LIBS    mymodel )

# ecbuild_add_test(
#   TARGET  test_mymodel_modelauxincrement
#   SOURCES executables/TestModelAuxIncrement.cc
#   ARGS    testinput/interface.yml
#   LIBS    mymodel )

# ecbuild_add_test(
#   TARGET  test_mymodel_linearmodel
#   SOURCES executables/TestLinearModel.cc
#   ARGS    testinput/interface.yml
#   LIBS    mymodel )  

# ecbuild_add_test(
#   TARGET  test_mymodel_geometryiterator
#   SOURCES executables/TestGeometryIterator.cc
#   ARGS    testinput/interface.yml
#   LIBS    mymodel )  
