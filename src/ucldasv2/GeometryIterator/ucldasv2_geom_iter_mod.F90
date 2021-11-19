! (C) Copyright 2019-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Geometry iterator
module ucldasv2_geom_iter_mod

use kinds, only : kind_real
use ucldasv2_geom_mod, only: ucldasv2_geom

implicit none
private


! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> Geometry iterator
!!
!! When initialized, the iterator points to the first valid local grid cell.
!! Calls to ucldasv2_geom_iter::next() moves the iterator forward, and calls to
!! ucldasv2_geom_iter::current() retrieves the lat/lon of the current grid cell.
!! The iterator is mainly used by ucldasv2_increment_mod::ucldasv2_increment::getpoint()
!! and ucldasv2_increment_mod::ucldasv2_increment::setpoint()
type, public :: ucldasv2_geom_iter
  type(ucldasv2_geom), pointer :: geom => null() !< Geometry

  integer :: iind = 1  !< i index of current grid point
  integer :: jind = 1  !< j index of current grid point

contains

  !> \copybrief ucldasv2_geom_iter_setup \see ucldasv2_geom_iter_setup
  procedure :: setup => ucldasv2_geom_iter_setup

  !> \copybrief ucldasv2_geom_iter_clone \see ucldasv2_geom_iter_clone
  procedure :: clone => ucldasv2_geom_iter_clone

  !> \copybrief ucldasv2_geom_iter_equals \see ucldasv2_geom_iter_equals
  procedure :: equals => ucldasv2_geom_iter_equals

  !> \copybrief ucldasv2_geom_iter_current \see ucldasv2_geom_iter_current
  procedure :: current => ucldasv2_geom_iter_current

  !> \copybrief ucldasv2_geom_iter_next \see ucldasv2_geom_iter_next
  procedure :: next => ucldasv2_geom_iter_next

end type ucldasv2_geom_iter


! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
!> Setup for the geometry iterator
!!
!! \relates ucldasv2_geom_iter_mod::ucldasv2_geom_iter
subroutine ucldasv2_geom_iter_setup(self, geom, iind, jind)
  class(ucldasv2_geom_iter),    intent(inout) :: self
  type(ucldasv2_geom), pointer, intent(   in) :: geom !< Pointer to geometry
  integer,                  intent(   in) :: iind, jind  !< starting index

  ! Associate geometry
  self%geom => geom

  ! Define iind/jind for local tile
  self%iind = iind
  self%jind = jind

end subroutine ucldasv2_geom_iter_setup


! ------------------------------------------------------------------------------
!> Clone for the geometry iterator from \p other to \p self
!!
!! \relates ucldasv2_geom_iter_mod::ucldasv2_geom_iter
subroutine ucldasv2_geom_iter_clone(self, other)
  class(ucldasv2_geom_iter), intent(inout) :: self
  type(ucldasv2_geom_iter),  intent(   in) :: other !< Other geometry iterator to clone from

  ! Associate geometry
  self%geom => other%geom

  ! Copy iind/jind
  self%iind = other%iind
  self%jind = other%jind

end subroutine ucldasv2_geom_iter_clone


! ------------------------------------------------------------------------------
!> Check for the geometry iterator equality (pointing to same i/j location)
!!
!! \relates ucldasv2_geom_iter_mod::ucldasv2_geom_iter
subroutine ucldasv2_geom_iter_equals(self, other, equals)
  class(ucldasv2_geom_iter), intent( in) :: self
  type(ucldasv2_geom_iter),  intent( in) :: other  !< Other geometry iterator
  integer,               intent(out) :: equals !< Equality flag

  ! Initialization
  equals = 0

  ! Check equality
  if (associated(self%geom, other%geom) .and. (self%iind==other%iind) &
      .and. (self%jind==other%jind)) equals = 1

end subroutine ucldasv2_geom_iter_equals


! ------------------------------------------------------------------------------
!> Get geometry iterator current lat/lon
!!
!! \throws abor1_ftn aborts if iterator is out of bounds
!! \relates ucldasv2_geom_iter_mod::ucldasv2_geom_iter
subroutine ucldasv2_geom_iter_current(self, lon, lat)
  class(ucldasv2_geom_iter), intent( in) :: self
  real(kind_real),    intent(out) :: lat  !< Latitude
  real(kind_real),    intent(out) :: lon  !< Longitude

  ! Check iind/jind
  if (self%iind == -1 .AND. self%jind == -1) then
    ! special case of {-1,-1} means end of the grid
    lat = self%geom%lat(self%geom%iec,self%geom%jec)
    lon = self%geom%lon(self%geom%iec,self%geom%jec)
  elseif (self%iind < self%geom%isc .OR. self%iind > self%geom%iec .OR. &
          self%jind < self%geom%jsc .OR. self%jind > self%geom%jec) then
    ! outside of the grid
    call abor1_ftn('ucldasv2_geom_iter_current: iterator out of bounds')
  else
    ! inside of the grid
    lat = self%geom%lat(self%iind,self%jind)
    lon = self%geom%lon(self%iind,self%jind)
  endif

end subroutine ucldasv2_geom_iter_current


! ------------------------------------------------------------------------------
!> Update geometry iterator to next point
!!
!! \todo skip over masked points
!! \relates ucldasv2_geom_iter_mod::ucldasv2_geom_iter
subroutine ucldasv2_geom_iter_next(self)
  class(ucldasv2_geom_iter), intent(inout) :: self
  integer :: iind, jind

  iind = self%iind
  jind = self%jind

  ! do while ((iind.lt.self%geom%iec).and.(jind.lt.self%geom%jec))

    ! increment by 1
    if (iind.lt.self%geom%iec) then
      iind = iind + 1
    elseif (iind.eq.self%geom%iec) then
      iind = self%geom%isc
      jind = jind + 1
    end if

    ! ! skip this point if it is on land
    ! if (self%geom%mask2d(iind,jind).lt.1) then
    !   cycle
    ! else
    !   exit
    ! endif

  ! end do

  if (jind > self%geom%jec) then
      iind=-1
      jind=-1
  end if

  self%iind = iind
  self%jind = jind

end subroutine ucldasv2_geom_iter_next
! ------------------------------------------------------------------------------

end module ucldasv2_geom_iter_mod
