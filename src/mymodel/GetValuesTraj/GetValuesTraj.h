/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef MYMODEL_GETVALUESTRAJ_GETVALUESTRAJ_H_
#define MYMODEL_GETVALUESTRAJ_GETVALUESTRAJ_H_

#include <ostream>

#include "oops/util/Printable.h"

// ----------------------------------------------------------------------------

namespace mymodel {

  // GetValuesTraj class
  class GetValuesTraj : public util::Printable {
   public:
    GetValuesTraj();
    ~GetValuesTraj();

   private:
    void print(std::ostream & os) const;
  };

}  // namespace mymodel

#endif  // MYMODEL_GETVALUESTRAJ_GETVALUESTRAJ_H_
