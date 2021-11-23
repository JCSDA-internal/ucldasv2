/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UCLDASV2_FORTRAN_H_
#define UCLDASV2_FORTRAN_H_

namespace ucldasv2 {

  /// key type for ucldasv2_geom_mod::ucldasv2_geom
  typedef int F90geom;

  /// key type for ucldasv2_geom_iter_mod::ucldasv2_geom_iter
  typedef int F90iter;

  /// key type for ucldasv2_model_mod::ucldasv2_model
  typedef int F90model;

  /// key type for ufo_geovals_mod::ufo_geovals
  typedef int F90goms;

  /// key type for ucldasv2_fields_mod::ucldasv2_fields
  typedef int F90flds;

  /// key type for ucldasv2_getvalues_mod::ucldasv2_getvalues
  typedef int F90getval;

  /// key type for ucldasv2_covariance_mod::ucldasv2_cov
  typedef int F90bmat;

}  // namespace ucldasv2
#endif  // UCLDASV2_FORTRAN_H_
