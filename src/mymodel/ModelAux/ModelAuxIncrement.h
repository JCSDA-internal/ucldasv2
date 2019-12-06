/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef MYMODEL_MODELAUX_MODELAUXINCREMENT_H_
#define MYMODEL_MODELAUX_MODELAUXINCREMENT_H_

#include <ostream>
#include <string>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// forward declarations
namespace eckit {
  class Configuration;
}
namespace mymodel {
  class Geometry;
}

//-----------------------------------------------------------------------------

namespace mymodel {

  // ModelAuxControl class
  class ModelAuxIncrement : public util::Printable,
                            private util::ObjectCounter<ModelAuxIncrement> {
   public:
    static const std::string classname() {return "mymodel::ModelAuxIncrement";}

    ModelAuxIncrement(const ModelAuxIncrement &, const eckit::Configuration &);
    ModelAuxIncrement(const ModelAuxIncrement &, const bool);
    ModelAuxIncrement(const Geometry &, const eckit::Configuration &);
    ~ModelAuxIncrement();

    // Linear algebra operators
    void zero();
    ModelAuxIncrement & operator*=(const double);
    ModelAuxIncrement & operator+=(const ModelAuxIncrement &);
    ModelAuxIncrement & operator-=(const ModelAuxIncrement &);
    double norm() const;
    void axpy(const double, const ModelAuxIncrement &);
    double dot_product_with(const ModelAuxIncrement &) const;

   private:
    void print(std::ostream &) const;
  };
}  // namespace mymodel
#endif  // MYMODEL_MODELAUX_MODELAUXINCREMENT_H_
