!
! (C) Copyright 2019-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module ucldasv2_geom_iter_mod

  use iso_c_binding
  use kinds
  use ucldasv2_geom_mod, only: ucldasv2_geom

  implicit none

  private
  public :: ucldasv2_geom_iter
  public :: ucldasv2_geom_iter_registry
  public :: ucldasv2_geom_iter_setup, ucldasv2_geom_iter_clone

  type :: ucldasv2_geom_iter
    type(ucldasv2_geom), pointer :: geom => null() !< Geometry
    integer :: iind = 1  !< index e.g. lat(iind,jind)
    integer :: jind = 1  !< 
  end type ucldasv2_geom_iter

#define LISTED_TYPE ucldasv2_geom_iter

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: ucldasv2_geom_iter_registry

contains

  ! ------------------------------------------------------------------------------
  ! Public
  ! ------------------------------------------------------------------------------

  !> Linked list implementation
#include "oops/util/linkedList_c.f"

  ! ------------------------------------------------------------------------------
  !> Setup for the geometry iterator
  subroutine ucldasv2_geom_iter_setup(self, geom, iind, jind)

    ! Passed variables
    type(ucldasv2_geom_iter),     intent(inout) :: self !< Geometry iterator
    type(ucldasv2_geom), pointer, intent(   in) :: geom !< Geometry
    integer,                  intent(   in) :: iind, jind  !< Index

    ! Associate geometry
    self%geom => geom

    ! Define iind/jind for local tile
    self%iind = iind
    self%jind = jind

  end subroutine ucldasv2_geom_iter_setup

  ! ------------------------------------------------------------------------------
  !> Clone for the geometry iterator
  subroutine ucldasv2_geom_iter_clone(self, other)

    ! Passed variables
    type(ucldasv2_geom_iter), intent(inout) :: self  !< Geometry iterator
    type(ucldasv2_geom_iter), intent(   in) :: other !< Other geometry iterator

    ! Associate geometry
    self%geom => other%geom

    ! Copy iind/jind
    self%iind = other%iind
    self%jind = other%jind

  end subroutine ucldasv2_geom_iter_clone
end module ucldasv2_geom_iter_mod
