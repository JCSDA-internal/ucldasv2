/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UCLDASV2_GEOMETRYITERATOR_GEOMETRYITERATOR_H_
#define UCLDASV2_GEOMETRYITERATOR_GEOMETRYITERATOR_H_

#include <iterator>
#include <ostream>
#include <string>
#include "ucldasv2/Fortran.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// forward declarations
namespace eckit {
  class Configuration;
  namespace geometry {
    class Point2;
  }
}
namespace ucldasv2 {
  class Geometry;
}

// ----------------------------------------------------------------------------

namespace ucldasv2 {

  // Geometry class
  class GeometryIterator : public util::Printable,
                           private util::ObjectCounter<GeometryIterator> {
   public:
    static const std::string classname() {return "ucldasv2::GeometryIterator";}

    // constructors / destructor
    GeometryIterator(const GeometryIterator &);
    explicit GeometryIterator(const Geometry & geom,
                              const int & iindex = 1, const int & jindex = 1);
    ~GeometryIterator();

    // other operators
    bool operator==(const GeometryIterator &) const;
    bool operator!=(const GeometryIterator &) const;
    GeometryIterator& operator++();
    eckit::geometry::Point2 operator*() const;

    F90iter & toFortran() {return keyIter_;}
    const F90iter & toFortran() const {return keyIter_;}

   private:
    void print(std::ostream &) const;
    F90iter keyIter_;
  };
}  // namespace ucldasv2

#endif  // UCLDASV2_GEOMETRYITERATOR_GEOMETRYITERATOR_H_
