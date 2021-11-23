# (C) Copyright 2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# UCLAND
################################################################################

cmake_minimum_required( VERSION 3.12)
project( ucland VERSION 2020.4.0 LANGUAGES Fortran )

# ecbuild integration
find_package(ecbuild 3.3.2 REQUIRED)
include( ecbuild_system NO_POLICY_SCOPE )
ecbuild_declare_project()

list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include( ucland_compiler_flags )

# land_bgc option
option( ENABLE_LAND_BGC "Compile and use land bgc modules in UCLAND")


################################################################################
# Dependencies
################################################################################

find_package(OpenMP COMPONENTS C Fortran)
find_package(MPI REQUIRED COMPONENTS C Fortran)
find_package(NetCDF REQUIRED COMPONENTS Fortran )
find_package(fms 2020.4.0 REQUIRED )


################################################################################
# Sources
################################################################################

set( UCLAND_LINKER_LANGUAGE Fortran )

# ucland library
#-------------------------------------------------------------------------------
add_subdirectory( src )

list( APPEND ucland_src_files
    ${src_files}
)

if (${ENABLE_LAND_BGC})
    add_definitions( -D_USE_GENERIC_TRACER )
    add_subdirectory( pkg/land_BGC )
    list( APPEND ucland_src_files ${land_bgc_src_files} )
endif()

ecbuild_add_library( TARGET   ucland
                     SOURCES  ${ucland_src_files}
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${UCLAND_LINKER_LANGUAGE}
                    )

target_link_libraries(ucland PUBLIC NetCDF::NetCDF_Fortran)
target_link_libraries(ucland PUBLIC MPI::MPI_Fortran)
target_link_libraries(ucland PUBLIC OpenMP::OpenMP_C OpenMP::OpenMP_Fortran)
target_link_libraries(ucland PUBLIC fms)

# Fortran module output directory for build and install interfaces
set(MODULE_DIR module/${PROJECT_NAME}/${CMAKE_Fortran_COMPILER_ID}/${CMAKE_Fortran_COMPILER_VERSION})
set_target_properties(${PROJECT_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR})
install(DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR}/ DESTINATION ${MODULE_DIR})
target_include_directories(${PROJECT_NAME} INTERFACE
                            $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${MODULE_DIR}>
                            $<INSTALL_INTERFACE:${MODULE_DIR}>)

# Installed .h include locations
target_include_directories(ucland PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/framework>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)


################################################################################
# Finalise configuration
################################################################################

ecbuild_install_project( NAME ucland )
ecbuild_print_summary()
