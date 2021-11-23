# (C) Copyright 2017-2020 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# UCLDASV2
################################################################################

cmake_minimum_required( VERSION 3.12 FATAL_ERROR )

project( ucldasv2 VERSION 0.1 LANGUAGES C CXX Fortran)

find_package(ecbuild 3.3.2 REQUIRED)

set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH})

set( ENABLE_OS_TESTS           OFF CACHE BOOL "Disable OS tests" FORCE )
set( ENABLE_LARGE_FILE_SUPPORT OFF CACHE BOOL "Disable testing of large file support" FORCE )
set( ENABLE_MPI                ON CACHE BOOL "Compile with MPI" )

include( ecbuild_system NO_POLICY_SCOPE )

ecbuild_requires_macro_version( 2.5 )


################################################################################
# Project
################################################################################

ecbuild_declare_project()

ecbuild_enable_fortran( REQUIRED )

set( UCLDASV2_LINKER_LANGUAGE CXX )


################################################################################
# Dependencies
################################################################################
# Boost Headers
include_directories( ${Boost_INCLUDE_DIR} )
find_package( OpenMP COMPONENTS CXX Fortran )
find_package( MPI REQUIRED COMPONENTS C CXX Fortran )
find_package( NetCDF REQUIRED COMPONENTS Fortran )
find_package( eckit  1.11.6     REQUIRED )
find_package( fckit  0.7.0      REQUIRED )
find_package( oops   0.1.0      REQUIRED )
find_package( saber  0.0.1      REQUIRED )
find_package( ioda   0.1.0      REQUIRED )
find_package( ufo    0.1.0      REQUIRED )
find_package( fms    2020.4.0   REQUIRED)

if(OpenMP_FOUND)
    find_package( atlas 0.20.2 REQUIRED COMPONENTS OMP OMP_Fortran )
else()
    find_package( atlas 0.20.2 REQUIRED )
endif()

################################################################################
# Export package info
################################################################################
list( APPEND UCLDASV2_TPLS atlas atlas_f LAPACK MPI NetCDF )

set( UCLDASV2_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src
                       ${CMAKE_CURRENT_BINARY_DIR}/src
                       ${CMAKE_Fortran_MODULE_DIRECTORY} )
set( UCLDASV2_LIBRARIES ucldasv2 )

get_directory_property( UCLDASV2_DEFINITIONS COMPILE_DEFINITIONS )

foreach( _tpl ${UCLDASV2_TPLS} )
string( TOUPPER ${_tpl} TPL )
list( APPEND UCLDASV2_EXTRA_DEFINITIONS   ${${TPL}_DEFINITIONS}  ${${TPL}_TPL_DEFINITIONS}  )
list( APPEND UCLDASV2_EXTRA_INCLUDE_DIRS  ${${TPL}_INCLUDE_DIRS} ${${TPL}_TPL_INCLUDE_DIRS} )
list( APPEND UCLDASV2_EXTRA_LIBRARIES     ${${TPL}_LIBRARIES}    ${${TPL}_TPL_LIBRARIES}    )
endforeach()


################################################################################
# Sources
################################################################################
include( ucldasv2_compiler_flags )
include_directories( ${UCLDASV2_INCLUDE_DIRS} ${UCLDASV2_EXTRA_INCLUDE_DIRS} )

add_subdirectory( src )
add_subdirectory( test )

# Build Doxygen documentation
add_subdirectory( Documents )

if(ECBUILD_INSTALL_FORTRAN_MODULES)
  install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/${CMAKE_CFG_INTDIR} DESTINATION ${INSTALL_INCLUDE_DIR} )
endif()

################################################################################
# Finalise configuration
################################################################################

# prepares a tar.gz of the sources and/or binaries
ecbuild_install_project( NAME ucldasv2 )

# print the summary of the configuration
ecbuild_print_summary()
