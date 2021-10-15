/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UCLDASV2_GETVALUES_GETVALUES_H_
#define UCLDASV2_GETVALUES_GETVALUES_H_

#include <memory>
#include <ostream>
#include <string>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "ufo/Locations.h"

// forward declarations
namespace ufo {
  class GeoVaLs;
  class Locations;
}
namespace ucldasv2 {
  class Geometry;
  class State;
}

// ----------------------------------------------------------------------------

namespace ucldasv2 {

  // GetValues class: interpolate state to observation locations
  class GetValues : public util::Printable,
                    private util::ObjectCounter<GetValues> {
   public:
    static const std::string classname() {return "ucldasv2::GetValues";}

    // constructors, destructors
    GetValues(const Geometry &, const ufo::Locations & locs);
    virtual ~GetValues();

    // fills in geovals for all observations in the timeframe (t1, t2],
    void fillGeoVaLs(const State &,
                     const util::DateTime & t1,
                     const util::DateTime & t2,
                     ufo::GeoVaLs &) const;

   private:
    void print(std::ostream &) const;

    std::shared_ptr<const Geometry> geom_;
    ufo::Locations locs_;
  };
}  // namespace ucldasv2

#endif  // UCLDASV2_GETVALUES_GETVALUES_H_
