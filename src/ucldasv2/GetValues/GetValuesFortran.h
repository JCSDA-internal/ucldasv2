/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UCLDASV2_GETVALUES_GETVALUESFORTRAN_H_
#define UCLDASV2_GETVALUES_GETVALUESFORTRAN_H_

#include "ucldasv2/Fortran.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace ufo {
  class Locations;
}

namespace util {
  class DateTime;
  class Duration;
}

namespace ucldasv2 {

extern "C" {
  void ucldasv2_getvalues_create_f90(F90getval &,
                                 const F90geom &,
                                 const ufo::Locations &);
  void ucldasv2_getvalues_delete_f90(F90getval &);
  void ucldasv2_getvalues_fill_geovals_f90(const F90getval &,
                                       const F90geom &,
                                       const F90flds &,
                                       const util::DateTime &,
                                       const util::DateTime &,
                                       const ufo::Locations &,
                                       const F90goms &);
  void ucldasv2_getvalues_fill_geovals_tl_f90(const F90getval &,
                                          const F90geom &,
                                          const F90flds &,
                                          const util::DateTime &,
                                          const util::DateTime &,
                                          const ufo::Locations &,
                                          const F90goms &);
  void ucldasv2_getvalues_fill_geovals_ad_f90(const F90getval &,
                                          const F90geom &,
                                          const F90flds &,
                                          const util::DateTime &,
                                          const util::DateTime &,
                                          const ufo::Locations &,
                                          const F90goms &);
};  // extern "C"

// -------------------------------------------------------------------------------------------------

}  // namespace ucldasv2
#endif  // UCLDASV2_GETVALUES_GETVALUESFORTRAN_H_
