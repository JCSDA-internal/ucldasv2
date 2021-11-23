! (C) Copyright 2020-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> registry for ucldasv2_state_mod::ucldasv2_state instances for use in
!! Fortran/C++ interfaces of ucldasv2_state_mod_c
module ucldasv2_state_reg

use ucldasv2_state_mod

implicit none
private

!> Linked list interface - defines registry_t type
#define LISTED_TYPE ucldasv2_state
#include "oops/util/linkedList_i.f"

!> Global registry for ucldasv2_state instances
type(registry_t), public:: ucldasv2_state_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

end module ucldasv2_state_reg
