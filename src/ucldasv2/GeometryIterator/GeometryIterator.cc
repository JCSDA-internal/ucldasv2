/*
 * (C) Copyright 2019-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ucldasv2/Geometry/Geometry.h"
#include "ucldasv2/GeometryIterator/GeometryIterator.h"
#include "ucldasv2/GeometryIterator/GeometryIteratorFortran.h"

#include "eckit/config/Configuration.h"
#include "eckit/geometry/Point3.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------

namespace ucldasv2 {


// -----------------------------------------------------------------------------

GeometryIterator::GeometryIterator(const GeometryIterator& iter) {
  ucldasv2_geom_iter_clone_f90(keyIter_, iter.toFortran());
}

// -----------------------------------------------------------------------------

GeometryIterator::GeometryIterator(const Geometry& geom,
                                   const int & iindex, const int & jindex,
                                   const int & kindex) {
  ucldasv2_geom_iter_setup_f90(keyIter_, geom.toFortran(),
  iindex, jindex, kindex);
}


// -----------------------------------------------------------------------------

GeometryIterator::~GeometryIterator() {
  ucldasv2_geom_iter_delete_f90(keyIter_);
}

// -----------------------------------------------------------------------------

bool GeometryIterator::operator==(const GeometryIterator & other) const {
  int equals = 0;
  ucldasv2_geom_iter_equals_f90(keyIter_, other.toFortran(), equals);
  return (equals == 1);
}

// -----------------------------------------------------------------------------

bool GeometryIterator::operator!=(const GeometryIterator & other) const {
  int equals = 0;
  ucldasv2_geom_iter_equals_f90(keyIter_, other.toFortran(), equals);
  return (equals == 0);
}

// -----------------------------------------------------------------------------

eckit::geometry::Point3 GeometryIterator::operator*() const {
  double lat, lon, dep;
  ucldasv2_geom_iter_current_f90(keyIter_, lon, lat, dep);
  return eckit::geometry::Point3(lon, lat, dep);
}

// -----------------------------------------------------------------------------

GeometryIterator& GeometryIterator::operator++() {
  ucldasv2_geom_iter_next_f90(keyIter_);
  return *this;
}
// -----------------------------------------------------------------------------

int GeometryIterator::iteratorDimension() const {
  int dimension;
  ucldasv2_geom_iter_dimension_f90(keyIter_, dimension);
  return dimension;
}

// -----------------------------------------------------------------------------

void GeometryIterator::print(std::ostream & os) const {
  double lat, lon, dep;
  ucldasv2_geom_iter_current_f90(keyIter_, lon, lat, dep);
  os << "GeometryIterator, lat/lon/depth: " << lat << " / " << lon
     << " / " << dep << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ucldasv2
