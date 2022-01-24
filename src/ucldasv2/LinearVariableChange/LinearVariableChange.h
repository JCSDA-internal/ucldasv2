/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include <boost/ptr_container/ptr_vector.hpp>

#include "oops/base/LinearVariableChangeParametersBase.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/Printable.h"

#include "ucldasv2/LinearVariableChange/Base/LinearVariableChangeBase.h"

namespace ucldasv2 {

// -----------------------------------------------------------------------------

class LinearVariableChangeParameters :
  public oops::LinearVariableChangeParametersBase {
  OOPS_CONCRETE_PARAMETERS(LinearVariableChangeParameters,
                           oops::LinearVariableChangeParametersBase)
 public:
  oops::OptionalParameter<std::vector<LinearVariableChangeParametersWrapper>>
         linearVariableChangesWrapper{"linear variable changes", this};
};

// -----------------------------------------------------------------------------

class LinearVariableChange : public util::Printable {
 public:
  static const std::string classname() {return "ucldasv2::LinearVariableChange";}

  typedef LinearVariableChangeParameters Parameters_;

  // Vector of variable changes typedefs
  typedef typename boost::ptr_vector<LinearVariableChangeBase> LinVarChaVec_;
  typedef typename LinVarChaVec_::const_iterator icst_;
  typedef typename LinVarChaVec_::const_reverse_iterator ircst_;

  explicit LinearVariableChange(const Geometry &, const Parameters_ &);
  ~LinearVariableChange();

  void setTrajectory(const State &, const State &);

  void multiply(Increment &, const oops::Variables &) const;
  void multiplyInverse(Increment &, const oops::Variables &) const;
  void multiplyAD(Increment &, const oops::Variables &) const;
  void multiplyInverseAD(Increment &, const oops::Variables &) const;

 private:
  void print(std::ostream &) const override;
  Parameters_ params_;
  std::shared_ptr<const Geometry> geom_;
  LinVarChaVec_ linVarChas_;
};

// -----------------------------------------------------------------------------

}  // namespace ucldasv2