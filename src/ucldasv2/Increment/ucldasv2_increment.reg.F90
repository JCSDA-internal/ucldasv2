! (C) Copyright 2020-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> registry for ucldasv2_increment_mod::ucldasv2_increment instances for use in
!! Fortran/C++ interface of ucldasv2_increment_mod_c
module ucldasv2_increment_reg

use ucldasv2_increment_mod

implicit none
private

#define LISTED_TYPE ucldasv2_increment

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry for ucldasv2_increment instances
type(registry_t), public :: ucldasv2_increment_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

end module ucldasv2_increment_reg
