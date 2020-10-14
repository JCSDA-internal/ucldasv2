/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef MYMODEL_GEOMETRY_GEOMETRY_H_
#define MYMODEL_GEOMETRY_GEOMETRY_H_

#include <ostream>
#include <string>
#include <vector>

#include "eckit/mpi/Comm.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// forward declarations
namespace eckit {
  class Configuration;
}
namespace mymodel {
  class GeometryIterator;
}

// ----------------------------------------------------------------------------

namespace mymodel {

  // Geometry class
  class Geometry : public util::Printable,
                   private util::ObjectCounter<Geometry> {
   public:
    static const std::string classname() {return "mymodel::Geometry";}

    // constructors and destructor
    explicit Geometry(const eckit::Configuration &, const eckit::mpi::Comm &);
    Geometry(const Geometry &);
    ~Geometry();

    // accessors
    const eckit::mpi::Comm & getComm() const {return comm_;}

    // These are needed for the GeometryIterator Interface
    // TODO(template_impl) GeometryIterator begin() const; 
    // TODO(template_impl) GeometryIterator end() const; 

    // vertical coordinate (only needed for GETKF?)
    std::vector<double> verticalCoord(std::string &) const;

   private:
    void print(std::ostream &) const;

    const eckit::mpi::Comm & comm_;
  };
}  // namespace mymodel

#endif  // MYMODEL_GEOMETRY_GEOMETRY_H_
