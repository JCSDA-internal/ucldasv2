/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UCLDASV2_TRAITS_H_
#define UCLDASV2_TRAITS_H_

#include <string>

// TODO(template_impl) #include "ucldasv2/Covariance/Covariance.h"
#include "ucldasv2/Geometry/Geometry.h"
// TODO(template_impl) #include "ucldasv2/GeometryIterator/GeometryIterator.h"
// TODO(template_impl) #include "ucldasv2/GetValues/GetValues.h"
// TODO(template_impl) #include "ucldasv2/GetValues/LinearGetValues.h"
// TODO(template_impl) #include "ucldasv2/Increment/Increment.h"
// TODO(template_impl) #include "ucldasv2/ModelAux/ModelAuxControl.h"
// TODO(template_impl) #include "ucldasv2/ModelAux/ModelAuxCovariance.h"
// TODO(template_impl) #include "ucldasv2/ModelAux/ModelAuxIncrement.h"
// TODO(template_impl) #include "ucldasv2/State/State.h"

namespace ucldasv2 {

  struct Traits{
    static std::string name() {return "ucldasv2";}
    static std::string nameCovar() {return "ucldasv2Covar";}
    static std::string nameCovar4D() {return "ucldasv2Covar";}

    // Interfaces that ucldasv2 has to implement
    // ---------------------------------------------------
// TODO(template_impl) typedef ucldasv2::Covariance          Covariance;
    typedef ucldasv2::Geometry            Geometry;
// TODO(template_impl) typedef ucldasv2::GeometryIterator    GeometryIterator;
// TODO(template_impl) typedef ucldasv2::GetValues           GetValues;
// TODO(template_impl) typedef ucldasv2::Increment           Increment;
// TODO(template_impl) typedef ucldasv2::LinearGetValues     LinearGetValues;
// TODO(template_impl) typedef ucldasv2::ModelAuxControl     ModelAuxControl;
// TODO(template_impl) typedef ucldasv2::ModelAuxCovariance  ModelAuxCovariance;
// TODO(template_impl) typedef ucldasv2::ModelAuxIncrement   ModelAuxIncrement;
// TODO(template_impl) typedef ucldasv2::State               State;
  };
}  // namespace ucldasv2

#endif  // UCLDASV2_TRAITS_H_
