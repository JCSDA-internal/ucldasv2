/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "mymodel/Geometry/Geometry.h"
#include "mymodel/GeometryIterator/GeometryIterator.h"

#include "eckit/config/Configuration.h"

#include "oops/util/abor1_cpp.h"

namespace mymodel {

// ----------------------------------------------------------------------------

  Geometry::Geometry(const eckit::Configuration & conf,
                     const eckit::mpi::Comm & comm)
    : comm_(comm) {
    util::abor1_cpp("Geometry::Geometry() needs to be implemented.",
                    __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  Geometry::Geometry(const Geometry & other)
    : comm_(other.comm_) {
    util::abor1_cpp("Geometry::Geometry() needs to be implemented.",
                     __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  Geometry::~Geometry() {
    util::abor1_cpp("Geometry::~Geometry() needs to be implemented.",
                     __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

/* TODO(template_impl)
  GeometryIterator Geometry::begin() const {
    util::abor1_cpp("Geometry::begin() needs to be implemented.",
                    __FILE__, __LINE__);
    return GeometryIterator(*this, 0, 0);
  }
TODO(template_impl) */

// ----------------------------------------------------------------------------

/* TODO(template_impl)
  GeometryIterator Geometry::end() const {
    util::abor1_cpp("Geometry::end() needs to be implemented.",
                    __FILE__, __LINE__);
    return GeometryIterator(*this, 0, 0);
  }
TODO(template_impl) */

// ----------------------------------------------------------------------------
  void Geometry::print(std::ostream & os) const {
    util::abor1_cpp("Geometry::print() needs to be implemented.",
                    __FILE__, __LINE__);
    os << "Geometry: "
       << "(TODO, print diagnostic info about the geometry here)"
       << std::endl;
  }

// ----------------------------------------------------------------------------
  std::vector<double> Geometry::verticalCoord(std::string &) const {
    util::abor1_cpp("Geometry::verticalCoord() needs to be implemented.",
                    __FILE__, __LINE__);
    return {};
  }

// ----------------------------------------------------------------------------

}  // namespace mymodel
