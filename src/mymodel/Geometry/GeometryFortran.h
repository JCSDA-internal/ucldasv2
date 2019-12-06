/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef MYMODEL_GEOMETRY_GEOMETRYFORTRAN_H_
#define MYMODEL_GEOMETRY_GEOMETRYFORTRAN_H_

#include "eckit/config/Configuration.h"

namespace mymodel {

  typedef int F90geom;

  extern "C" {
    void mymodel_geometry_setup_f90(F90geom &,
                                    const eckit::Configuration * const *);
    void mymodel_geometry_clone_f90(F90geom &, const F90geom &);
    void mymodel_geometry_delete_f90(F90geom &);
  }
}  // namespace mymodel

#endif  // MYMODEL_GEOMETRY_GEOMETRYFORTRAN_H_
