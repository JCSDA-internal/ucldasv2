/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UCLDASV2_GEOMETRY_GEOMETRYFORTRAN_H_
#define UCLDASV2_GEOMETRY_GEOMETRYFORTRAN_H_

#include "ucldasv2/Fortran.h"

#include "oops/base/Variables.h"

// Forward declarations
namespace atlas {
  namespace field {
    class FieldSetImpl;
  }
  namespace functionspace {
    class FunctionSpaceImpl;
  }
}
namespace eckit {
  class Configuration;
}

namespace ucldasv2 {

  extern "C" {
    void ucldasv2_geo_setup_f90(F90geom &,
                            const eckit::Configuration * const &,
                            const eckit::mpi::Comm *);
    void ucldasv2_geo_set_atlas_lonlat_f90(const F90geom &,
                                       atlas::field::FieldSetImpl *);
    void ucldasv2_geo_set_atlas_functionspace_pointer_f90(const F90geom &,
                      atlas::functionspace::FunctionSpaceImpl *);
    void ucldasv2_geo_fill_atlas_fieldset_f90(const F90geom &,
                                          atlas::field::FieldSetImpl *);
    void ucldasv2_geo_clone_f90(F90geom &, const F90geom &);
    void ucldasv2_geo_delete_f90(F90geom &);
    void ucldasv2_geo_start_end_f90(const F90geom &, int &, int &, int &,
                                int &);
    void ucldasv2_geo_get_num_levels_f90(const F90geom &, const oops::Variables &,
                                    const size_t &, size_t[]);
  }
}  // namespace ucldasv2
#endif  // UCLDASV2_GEOMETRY_GEOMETRYFORTRAN_H_
