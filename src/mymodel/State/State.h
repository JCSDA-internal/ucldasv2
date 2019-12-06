/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef MYMODEL_STATE_STATE_H_
#define MYMODEL_STATE_STATE_H_

#include <memory>
#include <ostream>
#include <string>

#include "eckit/mpi/Comm.h"

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// forward declarations
namespace eckit {
  class Configuration;
}
namespace oops {
  class Variables;
}
namespace ufo {
  class GeoVaLs;
  class Locations;
}
namespace mymodel {
  class Fields;
  class GetValuesTraj;
  class Geometry;
  class Increment;
}

// ----------------------------------------------------------------------------

namespace mymodel {

  // State class
  class State : public util::Printable,
                private util::ObjectCounter<State> {
   public:
    static const std::string classname() {return "mymodel::State";}

    // constructor, destructor
    State(const Geometry &, const oops::Variables &,
          const eckit::Configuration &);
    State(const Geometry &, const State &);
    State(const State &);
    State & operator=(const State &);
    ~State();

    // interpolate state to observations locations
    void getValues(const ufo::Locations &,
                   const oops::Variables &,
                   ufo::GeoVaLs &) const;
    void getValues(const ufo::Locations &,
                   const oops::Variables &,
                   ufo::GeoVaLs &,
                   GetValuesTraj &) const;

    // interactions with increment
    State & operator+=(const Increment &);

    const util::DateTime & validTime() const;
    util::DateTime & validTime();
    double norm() const;
    void write(const eckit::Configuration &) const;
    void zero();
    void accumul(const double &, const State &);

   private:
    void print(std::ostream &) const;

    std::unique_ptr<Fields> fields_;
  };
}  // namespace mymodel

#endif  // MYMODEL_STATE_STATE_H_
