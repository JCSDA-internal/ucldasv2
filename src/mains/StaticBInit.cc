/*
 * (C) Copyright 2019-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ucldasv2/Traits.h"
#include "oops/runs/Run.h"
#include "oops/runs/StaticBInit.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  oops::StaticBInit<ucldasv2::Traits> bmat;
  return run.execute(bmat);
}
