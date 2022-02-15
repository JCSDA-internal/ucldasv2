/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#pragma once

#include "ucldasv2/Fortran.h"

namespace ucldasv2 {
  extern "C" {
    void ucldasv2_model2geovals_changevar_f90(const F90geom &,
                                          const F90flds &,
                                          F90flds &);
  }
}  // namespace ucldasv2
