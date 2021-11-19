!
! (C) Copyright 2019-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence
! Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> C++ interfaces for ucldasv2_geom_iter_mod::ucldasv2_geom_iter
module ucldasv2_geom_iter_mod_c

use iso_c_binding
use kinds
use ucldasv2_geom_iter_mod
use ucldasv2_geom_mod_c,  only : ucldasv2_geom_registry
use ucldasv2_geom_mod, only: ucldasv2_geom

implicit none
private


#define LISTED_TYPE ucldasv2_geom_iter

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry for ucldasv2_geom_iter instances
type(registry_t), public :: ucldasv2_geom_iter_registry


contains


!> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_geom_iter_mod::ucldasv2_geom_iter::setup()
subroutine ucldasv2_geom_iter_setup_c(c_key_self, c_key_geom, c_iindex, c_jindex) bind(c, name='ucldasv2_geom_iter_setup_f90')
  integer(c_int), intent(inout) :: c_key_self !< Geometry iterator
  integer(c_int), intent(   in) :: c_key_geom !< Geometry
  integer(c_int), intent(   in) :: c_iindex    !< Index
  integer(c_int), intent(   in) :: c_jindex    !< Index

  ! Local variables
  type(ucldasv2_geom_iter),     pointer :: self
  type(ucldasv2_geom),          pointer :: geom

  ! Interface
  call ucldasv2_geom_iter_registry%init()
  call ucldasv2_geom_iter_registry%add(c_key_self)
  call ucldasv2_geom_iter_registry%get(c_key_self, self)
  call ucldasv2_geom_registry%get(c_key_geom, geom)

  ! Call Fortran
  call self%setup(geom, c_iindex, c_jindex)

end subroutine ucldasv2_geom_iter_setup_c

! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_geom_iter_mod::ucldasv2_geom_iter::clone()
  subroutine ucldasv2_geom_iter_clone_c(c_key_self, c_key_other) bind(c, name='ucldasv2_geom_iter_clone_f90')
  integer(c_int), intent(inout) :: c_key_self  !< Geometry iterator
  integer(c_int), intent(   in) :: c_key_other !< Other geometry iterator

  ! Local variables
  type(ucldasv2_geom_iter), pointer :: self, other

  ! Interface
  call ucldasv2_geom_iter_registry%get(c_key_other, other)
  call ucldasv2_geom_iter_registry%init()
  call ucldasv2_geom_iter_registry%add(c_key_self)
  call ucldasv2_geom_iter_registry%get(c_key_self, self)

  ! Call Fortran
  call self%clone(other)

end subroutine ucldasv2_geom_iter_clone_c


! ------------------------------------------------------------------------------
!> !> C++ interface for deleting ucldasv2_geom_iter_mod::ucldasv2_geom_iter
subroutine ucldasv2_geom_iter_delete_c(c_key_self) bind(c, name='ucldasv2_geom_iter_delete_f90')
    integer(c_int), intent(inout) :: c_key_self !< Geometry iterator

    ! Clear interface
    call ucldasv2_geom_iter_registry%remove(c_key_self)

end subroutine ucldasv2_geom_iter_delete_c

! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_geom_iter_mod::ucldasv2_geom_iter::equals()
subroutine ucldasv2_geom_iter_equals_c(c_key_self, c_key_other, c_equals) bind(c, name='ucldasv2_geom_iter_equals_f90')
  integer(c_int), intent(inout) :: c_key_self  !< Geometry iterator
  integer(c_int), intent(   in) :: c_key_other !< Other geometry iterator
  integer(c_int), intent(inout) :: c_equals    !< Equality flag

  ! Local variables
  type(ucldasv2_geom_iter),pointer :: self,other

  ! Interface
  call ucldasv2_geom_iter_registry%get(c_key_self, self)
  call ucldasv2_geom_iter_registry%get(c_key_other, other)

  ! Call Fortran
  call self%equals(other, c_equals)

end subroutine ucldasv2_geom_iter_equals_c

! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_geom_iter_mod::ucldasv2_geom_iter::current()
subroutine ucldasv2_geom_iter_current_c(c_key_self, c_lon, c_lat) bind(c, name='ucldasv2_geom_iter_current_f90')
  integer(c_int), intent(   in) :: c_key_self !< Geometry iterator
  real(c_double), intent(inout) :: c_lat      !< Latitude
  real(c_double), intent(inout) :: c_lon      !< Longitude

  ! Local variables
  type(ucldasv2_geom_iter), pointer :: self

  ! Interface
  call ucldasv2_geom_iter_registry%get(c_key_self, self)

  ! Call Fortran
  call self%current(c_lon, c_lat)

end subroutine ucldasv2_geom_iter_current_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_geom_iter_mod::ucldasv2_geom_iter::next()
subroutine ucldasv2_geom_iter_next_c(c_key_self) bind(c, name='ucldasv2_geom_iter_next_f90')
  integer(c_int), intent(in) :: c_key_self !< Geometry iterator

  ! Local variables
  type(ucldasv2_geom_iter), pointer :: self

  ! Interface
  call ucldasv2_geom_iter_registry%get(c_key_self, self)

  ! Call Fortran
  call self%next()
end subroutine ucldasv2_geom_iter_next_c


end module ucldasv2_geom_iter_mod_c
