list( APPEND ucldasv2_test_input
  testinput/geometry.yml
  testinput/errorcovariance.yml
  testinput/geometryiterator.yml
  testinput/getvalues.yml
  testinput/lineargetvalues.yml
  testinput/linearmodel.yml
  testinput/increment.yml
  testinput/model.yml
  testinput/modelaux.yml
  testinput/state.yml
  )

set( ucldasv2_data_files
  Data/fields_metadata.yml
)

set( ucldasv2_model_param
  Data/20x20x4/ucldasv2_gridspec.nc
  Data/20x20x4/LND_input
)

set( ucldasv2_model_namelist
  Data/20x20x4/input.nml
)

set( ucldasv2_model_restarts
  Data/20x20x4/INPUT/lnd.res.nc
)

# link the input files for the tests
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
foreach(FILENAME ${ucldasv2_test_input})
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
    ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
    ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

# Create Data directory for model input and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/INPUT)
foreach(FILENAME ${ucldasv2_model_restarts})
     get_filename_component(filename ${FILENAME} NAME )
     execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/INPUT/${filename} )
endforeach(FILENAME)
install(FILES ${ucldasv2_model_restarts}
        DESTINATION ${INSTALL_DATA_DIR}/testdata/20x20x4/INPUT/ )

# LND's resource files
foreach(FILENAME ${ucldasv2_model_param} ${ucldasv2_data_files})
     get_filename_component(filename ${FILENAME} NAME )
     execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${filename} )
endforeach(FILENAME)
install(FILES ${ucldasv2_data_files}
        DESTINATION ${INSTALL_DATA_DIR}/testdata )

# FMS input.nml that may be overwritten during the testing
foreach(FILENAME ${ucldasv2_model_namelist})
     get_filename_component(filename ${FILENAME} NAME )
     execute_process( COMMAND ${CMAKE_COMMAND} -E copy
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/inputnml/${filename} )
endforeach(FILENAME)

#================================================================================
# Tests of class interfaces
#================================================================================

ecbuild_add_test(
  TARGET  test_ucldasv2_geometry
  SOURCES executables/TestGeometry.cc
  ARGS    testinput/geometry.yml
  LIBS    ucldasv2 )

ecbuild_add_test(
  TARGET  test_ucldasv2_state
  SOURCES executables/TestState.cc
  ARGS    testinput/state.yml
  LIBS    ucldasv2 )

ecbuild_add_test(
  TARGET  test_ucldasv2_increment
  SOURCES executables/TestIncrement.cc
  ARGS    testinput/increment.yml
  LIBS    ucldasv2 )

ecbuild_add_test(
  TARGET  test_ucldasv2_geometryiterator
  SOURCES executables/TestGeometryIterator.cc
  ARGS    testinput/geometryiterator.yml
  LIBS    ucldasv2 )

# TODO(template_impl) ecbuild_add_test(
# TODO(template_impl)   TARGET  test_ucldasv2_model
# TODO(template_impl)   SOURCES executables/TestModel.cc
# TODO(template_impl)   ARGS    testinput/model.yml
# TODO(template_impl)   LIBS    ucldasv2 )

# TODO(template_impl) ecbuild_add_test(
# TODO(template_impl)   TARGET  test_ucldasv2_modelauxcontrol
# TODO(template_impl)   SOURCES executables/TestModelAuxControl.cc
# TODO(template_impl)   ARGS    testinput/modelaux.yml
# TODO(template_impl)   LIBS    ucldasv2 )

# TODO(template_impl) ecbuild_add_test(
# TODO(template_impl)   TARGET  test_ucldasv2_getvalues
# TODO(template_impl)   SOURCES executables/TestGetValues.cc
# TODO(template_impl)   ARGS    testinput/getvalues.yml
# TODO(template_impl)   LIBS    ucldasv2 )

# TODO(template_impl) ecbuild_add_test(
# TODO(template_impl)   TARGET  test_ucldasv2_lineargetvalues
# TODO(template_impl)   SOURCES executables/TestLinearGetValues.cc
# TODO(template_impl)   ARGS    testinput/lineargetvalues.yml
# TODO(template_impl)   LIBS    ucldasv2 )

# TODO(template_impl) ecbuild_add_test(
# TODO(template_impl)   TARGET  test_ucldasv2_errorcovariance
# TODO(template_impl)   SOURCES executables/TestErrorCovariance.cc
# TODO(template_impl)   ARGS    testinput/errorcovariance.yml
# TODO(template_impl)   LIBS    ucldasv2 )

# TODO(template_impl) ecbuild_add_test(
# TODO(template_impl)   TARGET  test_ucldasv2_modelauxcovariance
# TODO(template_impl)   SOURCES executables/TestModelAuxCovariance.cc
# TODO(template_impl)   ARGS    testinput/modelaux.yml
# TODO(template_impl)   LIBS    ucldasv2 )

# TODO(template_impl) ecbuild_add_test(
# TODO(template_impl)   TARGET  test_ucldasv2_modelauxincrement
# TODO(template_impl)   SOURCES executables/TestModelAuxIncrement.cc
# TODO(template_impl)   ARGS    testinput/modelaux.yml
# TODO(template_impl)   LIBS    ucldasv2 )

# TODO(template_impl) ecbuild_add_test(
# TODO(template_impl)   TARGET  test_ucldasv2_linearmodel
# TODO(template_impl)   SOURCES executables/TestLinearModel.cc
# TODO(template_impl)   ARGS    testinput/linearmodel.yml
# TODO(template_impl)   LIBS    ucldasv2 )
