/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UCLDASV2_MODELAUX_MODELAUXCONTROL_H_
#define UCLDASV2_MODELAUX_MODELAUXCONTROL_H_

#include <ostream>
#include <string>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// forward declarations
namespace eckit {
  class Configuration;
}
namespace ucldasv2 {
  class Geometry;
  class ModelAuxIncrement;
}

//-----------------------------------------------------------------------------

namespace ucldasv2 {

  // ModelAuxControl class
  class ModelAuxControl : public util::Printable,
                          private util::ObjectCounter<ModelAuxControl> {
   public:
    static const std::string classname() {return "ucldasv2::ModelAuxControl";}

    // constructors / destructor
    ModelAuxControl(const Geometry &, const eckit::Configuration &);
    ModelAuxControl(const Geometry &, const ModelAuxControl &);
    ModelAuxControl(const ModelAuxControl &, const bool);
    ~ModelAuxControl();

    // math operators
    ModelAuxControl & operator +=(const ModelAuxIncrement &);

   private:
    void print(std::ostream & os) const;
  };
}  // namespace ucldasv2
#endif  // UCLDASV2_MODELAUX_MODELAUXCONTROL_H_
