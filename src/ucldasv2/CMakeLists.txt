#include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

# the source code in this directory,
# all other files are handled through the subdirectory's respective CMakeLists.txt
list( APPEND ucldasv2_src_files
  Traits.h
  Fortran.h
)

# a convenience frunction to make adding source files from subdirectories easier
# (To keep from having to use the full path every time)
function(ucldasv2_target_sources)
  foreach( src IN LISTS ARGN)
    target_sources(ucldasv2 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/${src})
  endforeach()
endfunction()

# the main ucldasv2 library for this interface
ecbuild_add_library( TARGET   ucldasv2
                     SOURCES  ${ucldasv2_src_files}
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${UCLDASV2_LINKER_LANGUAGE}
                    )

target_include_directories(ucldasv2 PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
  $<INSTALL_INTERFACE:include/ucldasv2>)

target_compile_features( ucldasv2 PUBLIC cxx_std_11 )

target_link_libraries( ucldasv2 PUBLIC NetCDF::NetCDF_Fortran )
target_link_libraries( ucldasv2 PUBLIC fckit )
target_link_libraries( ucldasv2 PUBLIC atlas )
target_link_libraries( ucldasv2 PUBLIC atlas_f )
target_link_libraries( ucldasv2 PUBLIC fms )
target_link_libraries( ucldasv2 PUBLIC ucland )
target_link_libraries( ucldasv2 PUBLIC oops )
target_link_libraries( ucldasv2 PUBLIC saber )
target_link_libraries( ucldasv2 PUBLIC ioda )
target_link_libraries( ucldasv2 PUBLIC ufo )

# add source code in the subdirectories
# TODO(template_impl) add_subdirectory(Covariance)
add_subdirectory(Fields)
add_subdirectory(Geometry)
add_subdirectory(GeometryIterator)
add_subdirectory(GetValues)
add_subdirectory(Increment)
add_subdirectory(Model)
# TODO(template_impl) add_subdirectory(ModelAux)
add_subdirectory(State)
add_subdirectory(Utils)

## Fortran modules
set(MODULE_DIR ${PROJECT_NAME}/module)
set_target_properties(${PROJECT_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR})
install(DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR}/ DESTINATION ${CMAKE_INSTALL_LIBDIR}/${MODULE_DIR})
target_include_directories(${PROJECT_NAME} INTERFACE
                            $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${MODULE_DIR}>
                            $<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/${MODULE_DIR}>)
