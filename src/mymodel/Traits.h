/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef MYMODEL_TRAITS_H_
#define MYMODEL_TRAITS_H_

#include <string>

// #include "mymodel/Covariance/Covariance.h"
#include "mymodel/Geometry/Geometry.h"
// #include "mymodel/GeometryIterator/GeometryIterator.h"
// #include "mymodel/GetValuesTraj/GetValuesTraj.h"
// #include "mymodel/Increment/Increment.h"
// #include "mymodel/ModelAux/ModelAuxCovariance.h"
// #include "mymodel/ModelAux/ModelAuxControl.h"
// #include "mymodel/ModelAux/ModelAuxIncrement.h"
// #include "mymodel/State/State.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasCovariance.h"
#include "ufo/ObsBiasIncrement.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsOperator.h"

namespace mymodel {

  struct Traits{
    static std::string name() {return "mymodel";}
    static std::string nameCovar() {return "mymodelCovar";}
    static std::string nameCovar4D() {return "mymodelCovar";}

    // Interfaces that mymodel has to implement
    // ---------------------------------------------------
    // typedef mymodel::Covariance          Covariance;
    typedef mymodel::Geometry            Geometry;
    // typedef mymodel::GeometryIterator    GeometryIterator;
    // typedef mymodel::GetValuesTraj       InterpolatorTraj;
    // typedef mymodel::Increment           Increment;
    // typedef mymodel::State               State;
    // typedef mymodel::ModelAuxCovariance  ModelAuxCovariance;
    // typedef mymodel::ModelAuxControl     ModelAuxControl;
    // typedef mymodel::ModelAuxIncrement   ModelAuxIncrement;

    // Interfaces that are already provided by JEDI
    typedef ufo::GeoVaLs              GeoVaLs;
    typedef ufo::Locations            Locations;
    typedef ufo::ObsBias              ObsAuxControl;
    typedef ufo::ObsBiasCovariance    ObsAuxCovariance;
    typedef ufo::ObsBiasIncrement     ObsAuxIncrement;
    typedef ufo::ObsDiagnostics       ObsDiagnostics;
    typedef ufo::ObsOperator          ObsOperator;
    typedef ioda::ObsSpace            ObsSpace;
    typedef ioda::ObsVector           ObsVector;
    template <typename DATA> using ObsDataVector = ioda::ObsDataVector<DATA>;
  };
}  // namespace mymodel

#endif  // MYMODEL_TRAITS_H_
