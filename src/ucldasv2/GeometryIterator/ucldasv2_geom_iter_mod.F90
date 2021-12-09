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

  integer :: iindex = 1  !< i index of current grid point
  integer :: jindex = 1  !< j index of current grid point
  integer :: kindex = 1  !< k index of current grid point 

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
subroutine ucldasv2_geom_iter_setup(self, geom, iindex, jindex, kindex)
  class(ucldasv2_geom_iter),    intent(inout) :: self
  type(ucldasv2_geom), pointer, intent(   in) :: geom !< Pointer to geometry
  integer,         intent(   in) :: iindex, jindex, kindex  !< starting index

  ! Associate geometry
  self%geom => geom

  ! Define iindex/jindex/kindex for local tile
  self%iindex = iindex
  self%jindex = jindex
  self%kindex = kindex

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

  ! Copy iindex/jindex/kindex
  self%iindex = other%iindex
  self%jindex = other%jindex
  self%kindex = other%kindex

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
  if (associated(self%geom, other%geom)) then
    select case(self%geom%iterator_dimension)
    case (2) ! 2-d iterator
      if ((self%iindex==other%iindex) .and. (self%jindex==other%jindex)) equals = 1
    case (3) ! 3-d iterator
      if ((self%iindex==other%iindex) .and. (self%jindex==other%jindex) .and. &
          (self%kindex==other%kindex) ) equals = 1
    case default
      call abor1_ftn('ucldasv2_geom_iter_equals: unknown geom%iterator_dimension')
    end select
  endif

end subroutine ucldasv2_geom_iter_equals


! ------------------------------------------------------------------------------
!> Get geometry iterator current lat/lon
!!
!! \throws abor1_ftn aborts if iterator is out of bounds
!! \relates ucldasv2_geom_iter_mod::ucldasv2_geom_iter
subroutine ucldasv2_geom_iter_current(self, lon, lat, lev)

  ! Passed variables
  class(ucldasv2_geom_iter), intent( in) :: self !< Geometry iterator
  real(kind_real),    intent(out) :: lat  !< Latitude
  real(kind_real),    intent(out) :: lon  !< Longitude
  real(kind_real),    intent(out) :: lev  !< Level

  ! Check iindex/jindex
  if (self%iindex == -1 .AND. self%jindex == -1) then
    ! special case of {-1,-1} means end of the grid
    lat = self%geom%lat(self%geom%iec,self%geom%jec)
    lon = self%geom%lon(self%geom%iec,self%geom%jec)
  elseif (self%iindex < self%geom%isc .OR. self%iindex > self%geom%iec .OR. &
          self%jindex < self%geom%jsc .OR. self%jindex > self%geom%jec) then
    ! outside of the grid
    call abor1_ftn('ucldasv2_geom_iter_current: lat/lon iterator out of bounds')
  else
    ! inside of the grid
    lat = self%geom%lat(self%iindex,self%jindex)
    lon = self%geom%lon(self%iindex,self%jindex)
  endif

  ! check kindex
  select case(self%geom%iterator_dimension)
  case (2) ! 2-d iterator
    lev = -99999
  case (3) ! 3-d iterator
    lev = self%kindex
  case default
    call abor1_ftn('ucldasv2_geom_iter_current: unknown geom%iterator_dimension')
  end select

end subroutine ucldasv2_geom_iter_current


! ------------------------------------------------------------------------------
!> Update geometry iterator to next point
!!
!! \todo skip over masked points
!! \relates ucldasv2_geom_iter_mod::ucldasv2_geom_iter
subroutine ucldasv2_geom_iter_next(self)
  class(ucldasv2_geom_iter), intent(inout) :: self
  integer :: iindex, jindex, kindex

  iindex = self%iindex
  jindex = self%jindex
  kindex = self%kindex

  ! increment by 1
  select case(self%geom%iterator_dimension)
  case (2) ! 2-d iterator
    if (iindex.lt.self%geom%iec) then
      iindex = iindex + 1
    elseif (iindex.eq.self%geom%iec) then
      iindex = self%geom%isc
      jindex = jindex + 1
    end if

    if (jindex > self%geom%jec) then
      iindex=-1
      jindex=-1
    end if
  case (3) ! 3-d iterator
    if (iindex.lt.self%geom%iec) then
      iindex = iindex + 1
    elseif (iindex.eq.self%geom%iec) then
      iindex = self%geom%isc
      if (jindex.lt.self%geom%jec) then
        jindex = jindex + 1
      elseif (jindex.eq.self%geom%jec) then
        jindex = self%geom%jsc
        kindex = kindex + 1
      end if !j loop
    end if !iloop

    if (kindex > self%geom%nsnow) then
      iindex=-1
      jindex=-1
      kindex=-1
    end if !kloop
  case default
    call abor1_ftn('ucldasv2_geom_iter_next: unknown geom%iterator_dimension')
  end select

  self%iindex = iindex
  self%jindex = jindex
  self%kindex = kindex

end subroutine ucldasv2_geom_iter_next
! ------------------------------------------------------------------------------

end module ucldasv2_geom_iter_mod
