/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UCLDASV2_GEOMETRY_GEOMETRY_H_
#define UCLDASV2_GEOMETRY_GEOMETRY_H_

#include <fstream>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"

#include "ucldasv2/Fortran.h"
#include "ucldasv2/Geometry/FmsInput.h"
#include "ucldasv2/Geometry/GeometryFortran.h"
#include "ucldasv2/GeometryIterator/GeometryIterator.h"
#include "ucldasv2/GeometryIterator/GeometryIteratorFortran.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// forward declarations
namespace atlas {
  class FieldSet;
  class FunctionSpace;
  namespace functionspace {
    class PointCloud;
  }
}
namespace oops {
  class Variables;
}

// ---------------------------------------------------------------------------

namespace ucldasv2 {

  /// Geometry handles geometry for UCLDASV2 model.
  class Geometry : public util::Printable,
                   private util::ObjectCounter<Geometry> {
   public:
    static const std::string classname() {return "ucldasv2::Geometry";}

    explicit Geometry(const eckit::Configuration &, const eckit::mpi::Comm &);
    Geometry(const Geometry &);
    ~Geometry();

    GeometryIterator begin() const;
    GeometryIterator end() const;
    std::vector<size_t> variableSizes(const oops::Variables & vars) const;
    std::vector<double> verticalCoord(std::string &) const {return {};}

    int& toFortran() {return keyGeom_;}
    const int& toFortran() const {return keyGeom_;}
    const eckit::mpi::Comm & getComm() const {return comm_;}

    atlas::FunctionSpace * atlasFunctionSpace() const;
    atlas::FieldSet * atlasFieldSet() const;

   private:
    Geometry & operator=(const Geometry &);
    void print(std::ostream &) const;
    int keyGeom_;
    const eckit::mpi::Comm & comm_;
    FmsInput fmsinput_;
    std::unique_ptr<atlas::functionspace::PointCloud> atlasFunctionSpace_;
    std::unique_ptr<atlas::FieldSet> atlasFieldSet_;
  };
}  // namespace ucldasv2

#endif  // UCLDASV2_GEOMETRY_GEOMETRY_H_
