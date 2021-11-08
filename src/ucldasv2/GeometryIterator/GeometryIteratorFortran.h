/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UCLDASV2_GEOMETRYITERATOR_GEOMETRYITERATORFORTRAN_H_
#define UCLDASV2_GEOMETRYITERATOR_GEOMETRYITERATORFORTRAN_H_

#include "ucldasv2/Fortran.h"

namespace ucldasv2 {

  extern "C" {
    void ucldasv2_geom_iter_setup_f90(F90iter &, const F90geom &,
                                  const int &, const int &);
    void ucldasv2_geom_iter_clone_f90(F90iter &, const F90iter &);
    void ucldasv2_geom_iter_delete_f90(F90iter &);
  }
}  // namespace ucldasv2
#endif  // UCLDASV2_GEOMETRYITERATOR_GEOMETRYITERATORFORTRAN_H_
