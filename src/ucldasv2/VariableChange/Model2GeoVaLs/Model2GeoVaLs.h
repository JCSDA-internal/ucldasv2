/*
 * (C) Copyright 2021-2021  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>

#include "ucldasv2/Traits.h"

#include "ucldasv2/VariableChange/Base/VariableChangeBase.h"

namespace ucldasv2 {

class Model2GeoVaLs: public VariableChangeBase {
 public:
  const std::string classname() {return "ucldasv2::Model2GeoVaLs";}

  Model2GeoVaLs(const Geometry &, const eckit::Configuration &);
  ~Model2GeoVaLs();

  void changeVar(const State &, State &) const override;
  void changeVarInverse(const State &, State &) const override;

 private:
  std::unique_ptr<Geometry> geom_;
  void print(std::ostream &) const override {}
};

}  // namespace ucldasv2
