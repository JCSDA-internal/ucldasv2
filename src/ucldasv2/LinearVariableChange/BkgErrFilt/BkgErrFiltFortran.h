/*
 * (C) Copyright 2017-2021  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "soca/Fortran.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace soca {

  extern "C" {
    void soca_bkgerrfilt_setup_f90(F90balopmat &,
                                   const eckit::Configuration * const *,
                                   const F90flds &,
                                   const F90geom &);
    void soca_bkgerrfilt_delete_f90(F90balopmat &);
    void soca_bkgerrfilt_mult_f90(const F90balopmat &,
                                  const F90flds &,
                                  F90flds &);
    void soca_bkgerrfilt_multad_f90(const F90balopmat,
                                    F90flds &,
                                    const F90flds &);
  }
}  // namespace soca
