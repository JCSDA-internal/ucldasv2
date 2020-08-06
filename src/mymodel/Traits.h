/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef MYMODEL_TRAITS_H_
#define MYMODEL_TRAITS_H_

#include <string>

// TODO(template_impl) #include "mymodel/Covariance/Covariance.h"
#include "mymodel/Geometry/Geometry.h"
// TODO(template_impl) #include "mymodel/GeometryIterator/GeometryIterator.h"
// TODO(template_impl) #include "mymodel/GetValues/GetValues.h"
// TODO(template_impl) #include "mymodel/GetValues/LinearGetValues.h"
// TODO(template_impl) #include "mymodel/Increment/Increment.h"
// TODO(template_impl) #include "mymodel/ModelAux/ModelAuxControl.h"
// TODO(template_impl) #include "mymodel/ModelAux/ModelAuxCovariance.h"
// TODO(template_impl) #include "mymodel/ModelAux/ModelAuxIncrement.h"
// TODO(template_impl) #include "mymodel/State/State.h"

namespace mymodel {

  struct Traits{
    static std::string name() {return "mymodel";}
    static std::string nameCovar() {return "mymodelCovar";}
    static std::string nameCovar4D() {return "mymodelCovar";}

    // Interfaces that mymodel has to implement
    // ---------------------------------------------------
// TODO(template_impl) typedef mymodel::Covariance          Covariance;
    typedef mymodel::Geometry            Geometry;
// TODO(template_impl) typedef mymodel::GeometryIterator    GeometryIterator;
// TODO(template_impl) typedef mymodel::GetValues           GetValues;
// TODO(template_impl) typedef mymodel::Increment           Increment;
// TODO(template_impl) typedef mymodel::LinearGetValues     LinearGetValues;
// TODO(template_impl) typedef mymodel::ModelAuxControl     ModelAuxControl;
// TODO(template_impl) typedef mymodel::ModelAuxCovariance  ModelAuxCovariance;
// TODO(template_impl) typedef mymodel::ModelAuxIncrement   ModelAuxIncrement;
// TODO(template_impl) typedef mymodel::State               State;
  };
}  // namespace mymodel

#endif  // MYMODEL_TRAITS_H_
