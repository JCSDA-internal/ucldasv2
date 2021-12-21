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
#include "ucldasv2/GeometryIterator/GeometryIterator.h"
#include "ucldasv2/GetValues/GetValues.h"
#include "ucldasv2/GetValues/LinearGetValues.h"
#include "ucldasv2/Increment/Increment.h"
// TODO(template_impl) #include "ucldasv2/ModelAux/ModelAuxControl.h"
// TODO(template_impl) #include "ucldasv2/ModelAux/ModelAuxCovariance.h"
// TODO(template_impl) #include "ucldasv2/ModelAux/ModelAuxIncrement.h"
#include "ucldasv2/State/State.h"

namespace ucldasv2 {

  struct Traits{
    static std::string name() {return "ucldasv2";}
    static std::string nameCovar() {return "ucldasv2Error";}

    // Interfaces that ucldasv2 has to implement
    // ---------------------------------------------------
// TODO(template_impl) typedef ucldasv2::Covariance          Covariance;
    typedef ucldasv2::Geometry            Geometry;
    typedef ucldasv2::GeometryIterator    GeometryIterator;
    typedef ucldasv2::GetValues           GetValues;
    typedef ucldasv2::Increment           Increment;
    typedef ucldasv2::LinearGetValues     LinearGetValues;
// TODO(template_impl) typedef ucldasv2::ModelAuxControl     ModelAuxControl;
// TODO(template_impl) typedef ucldasv2::ModelAuxCovariance  ModelAuxCovariance;
// TODO(template_impl) typedef ucldasv2::ModelAuxIncrement   ModelAuxIncrement;
    typedef ucldasv2::State               State;
  };
}  // namespace ucldasv2

#endif  // UCLDASV2_TRAITS_H_
