/*
 * (C) Copyright 2019-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ucldasv2/Traits.h"
#include "oops/runs/Run.h"
#include "oops/runs/Variational.h"
#include "ufo/instantiateObsFilterFactory.h"
#include "ufo/ObsTraits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  ufo::instantiateObsFilterFactory<ufo::ObsTraits>();
  oops::Variational<ucldasv2::Traits, ufo::ObsTraits> var;
  return run.execute(var);
}
