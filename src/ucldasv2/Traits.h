/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UCLDASV2_TRAITS_H_
#define UCLDASV2_TRAITS_H_

#include <string>

// TODO(template_impl) #include "ucldasv2/Covariance/ErrorCovariance.h"
#include "ucldasv2/Geometry/Geometry.h"
#include "ucldasv2/GeometryIterator/GeometryIterator.h"
#include "ucldasv2/GetValues/GetValues.h"
#include "ucldasv2/GetValues/LinearGetValues.h"
#include "ucldasv2/Increment/Increment.h"
#include "ucldasv2/LinearVariableChange/LinearVariableChange.h"
// TODO(template_impl) #include "ucldasv2/ModelBias/ModelBias.h"
// TODO(template_impl) #include "ucldasv2/ModelBias/ModelBiasCovariance.h"
// TODO(template_impl) #include "ucldasv2/ModelBias/ModelBiasIncrement.h"
#include "ucldasv2/State/State.h"
#include "ucldasv2/VariableChange/VariableChange.h"

namespace ucldasv2 {

/**
 * \brief The main traits structure for UCLDASV2.
 *
 * This structure is responsible for supplying UCLDASV2 specific code to the JEDI
 * applications within \ref src/mains directory.
 */
struct Traits {
  static std::string name() {return "UCLDASV2";}
  static std::string nameCovar() {return "SocaError";}

  typedef ucldasv2::Geometry             Geometry;
  typedef ucldasv2::GeometryIterator     GeometryIterator;
  typedef ucldasv2::State                State;
  typedef ucldasv2::Increment            Increment;
// TODO(template_impl)   typedef ucldasv2::ErrorCovariance      Covariance;
  typedef ucldasv2::GetValues            GetValues;
  typedef ucldasv2::LinearGetValues      LinearGetValues;

// TODO(template_impl)   typedef ucldasv2::ModelBias            ModelAuxControl;
// TODO(template_impl)   typedef ucldasv2::ModelBiasIncrement   ModelAuxIncrement;
// TODO(template_impl)   typedef ucldasv2::ModelBiasCovariance  ModelAuxCovariance;

  typedef ucldasv2::LinearVariableChange LinearVariableChange;
  typedef ucldasv2::VariableChange       VariableChange;
};

}  // namespace ucldasv2

#endif  // UCLDASV2_TRAITS_H_
