/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef MYMODEL_FIELDS_FIELDS_H_
#define MYMODEL_FIELDS_FIELDS_H_

#include <ostream>
#include <string>

#include <boost/shared_ptr.hpp>

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
namespace mymodel {
  class Geometry;
}

// ----------------------------------------------------------------------------

namespace mymodel {

  // Fields class
  class Fields : public util::Printable,
                 private util::ObjectCounter<Fields> {
   public:
    static const std::string classname() {return "mymodel::Fields";}

    Fields(const Geometry &, const oops::Variables &,
           const eckit::Configuration &);
    ~Fields();

    const util::DateTime & time() const { return time_; }
    util::DateTime & time() { return time_; }

    double norm() const;

    boost::shared_ptr<const Geometry> geometry() const;

   private:
    void print(std::ostream &) const;

    boost::shared_ptr<const Geometry> geom_;
    util::DateTime time_;
  };

}  // namespace mymodel

#endif  // MYMODEL_FIELDS_FIELDS_H_
