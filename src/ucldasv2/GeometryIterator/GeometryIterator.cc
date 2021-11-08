/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include "ucldasv2/Geometry/Geometry.h"
#include "ucldasv2/GeometryIterator/GeometryIterator.h"
#include "ucldasv2/GeometryIterator/GeometryIteratorFortran.h"

#include "eckit/config/Configuration.h"
#include "eckit/geometry/Point2.h"

#include "oops/util/abor1_cpp.h"

namespace ucldasv2 {

// ----------------------------------------------------------------------------

  GeometryIterator::GeometryIterator(const Geometry& geom,
                                     const int & iindex, const int & jindex) {
    ucldasv2_geom_iter_setup_f90(keyIter_, geom.toFortran(), iindex, jindex);
  }

// ----------------------------------------------------------------------------

  GeometryIterator::GeometryIterator(const GeometryIterator& iter) {
    ucldasv2_geom_iter_clone_f90(keyIter_, iter.toFortran());
  }

// ----------------------------------------------------------------------------

  GeometryIterator::~GeometryIterator() {
    ucldasv2_geom_iter_delete_f90(keyIter_);
  }

// ----------------------------------------------------------------------------

  bool GeometryIterator::operator==(const GeometryIterator &) const {
    util::abor1_cpp(
      "GeometryIterator::operator==() needs to be implemented.",
      __FILE__, __LINE__);
    return false;
  }

// ----------------------------------------------------------------------------

  bool GeometryIterator::operator!=(const GeometryIterator &) const {
    util::abor1_cpp(
      "GeometryIterator::operator!=() needs to be implemented.",
      __FILE__, __LINE__);
    return false;
  }

// ----------------------------------------------------------------------------

  GeometryIterator& GeometryIterator::operator++() {
    util::abor1_cpp(
      "GeometryIterator::operator++() needs to be implemented.",
      __FILE__, __LINE__);
    return *this;
  }

// ----------------------------------------------------------------------------

  eckit::geometry::Point2 GeometryIterator::operator*() const {
    util::abor1_cpp("GeometryIterator::operator*() needs to be implemented.",
                     __FILE__, __LINE__);
    return eckit::geometry::Point2(0.0, 0.0);
  }

// ----------------------------------------------------------------------------

  void GeometryIterator::print(std::ostream  & os) const {
    util::abor1_cpp("GeometryIterator::print() needs to be implemented.",
                    __FILE__, __LINE__);
    os << "(TODO, print diagnostic info about the GeometryIterator here)"
       << std::endl;
  }

// ----------------------------------------------------------------------------

}  // namespace ucldasv2
