! (C) Copyright 2020-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! ------------------------------------------------------------------------------
!> C++ interfaces for ucldasv2_getvalues_mod::ucldasv2_getvalues
module ucldasv2_getvalues_mod_c

use datetime_mod, only: datetime, c_f_datetime
use iso_c_binding
use ufo_geovals_mod_c, only: ufo_geovals_registry
use ufo_geovals_mod, only: ufo_geovals
use ufo_locations_mod, only: ufo_locations

! ucldasv2 modules
use ucldasv2_geom_mod, only: ucldasv2_geom
use ucldasv2_geom_mod_c, only: ucldasv2_geom_registry
use ucldasv2_getvalues_mod, only: ucldasv2_getvalues
use ucldasv2_getvalues_reg, only: ucldasv2_getvalues_registry
use ucldasv2_state_mod, only: ucldasv2_state
use ucldasv2_state_reg, only: ucldasv2_state_registry
use ucldasv2_increment_mod, only: ucldasv2_increment
use ucldasv2_increment_reg, only: ucldasv2_increment_registry

implicit none
private


! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_getvalues_mod::ucldasv2_getvalues::create()
subroutine ucldasv2_getvalues_create_c(c_key_self, c_key_geom, c_locs) &
           bind (c, name='ucldasv2_getvalues_create_f90')
integer(c_int),     intent(inout) :: c_key_self      !< Key to self
integer(c_int),     intent(in)    :: c_key_geom      !< Key to geometry
type(c_ptr), value, intent(in)    :: c_locs

type(ucldasv2_getvalues), pointer :: self
type(ucldasv2_geom),      pointer :: geom
type(ufo_locations)           :: locs

! Create object
call ucldasv2_getvalues_registry%init()
call ucldasv2_getvalues_registry%add(c_key_self)
call ucldasv2_getvalues_registry%get(c_key_self, self)

! Others
call ucldasv2_geom_registry%get(c_key_geom, geom)
locs = ufo_locations(c_locs)

! Call method
call self%create(geom, locs)

end subroutine ucldasv2_getvalues_create_c


! --------------------------------------------------------------------------------------------------
!> C++ interface for ucldasv2_getvalues_mod::ucldasv2_getvalues::delete()
subroutine ucldasv2_getvalues_delete_c(c_key_self) bind (c, name='ucldasv2_getvalues_delete_f90')
integer(c_int), intent(inout) :: c_key_self !< Key to self

type(ucldasv2_getvalues), pointer :: self

! Get object
call ucldasv2_getvalues_registry%get(c_key_self, self)

! Call method
call self%delete()

! Remove object
call ucldasv2_getvalues_registry%remove(c_key_self)

end subroutine ucldasv2_getvalues_delete_c


! --------------------------------------------------------------------------------------------------
!> C++ interface for ucldasv2_getvalues_mod::ucldasv2_getvalues::fill_geovals()
subroutine ucldasv2_getvalues_fill_geovals_c(c_key_self, c_key_geom, c_key_state, c_t1, c_t2, &
                                         c_locs, c_key_geovals) &
           bind (c, name='ucldasv2_getvalues_fill_geovals_f90')

integer(c_int),     intent(in) :: c_key_self
integer(c_int),     intent(in) :: c_key_geom
integer(c_int),     intent(in) :: c_key_state
type(c_ptr), value, intent(in) :: c_t1
type(c_ptr), value, intent(in) :: c_t2
type(c_ptr), value, intent(in) :: c_locs
integer(c_int),     intent(in) :: c_key_geovals

type(ucldasv2_getvalues), pointer :: self
type(ucldasv2_geom),      pointer :: geom
type(ucldasv2_state),     pointer :: state
type(datetime)                :: t1
type(datetime)                :: t2
type(ufo_locations)           :: locs
type(ufo_geovals),    pointer :: geovals

! Get objects
call ucldasv2_getvalues_registry%get(c_key_self, self)
call ucldasv2_geom_registry%get(c_key_geom, geom)
call ucldasv2_state_registry%get(c_key_state, state)
call c_f_datetime(c_t1, t1)
call c_f_datetime(c_t2, t2)
locs = ufo_locations(c_locs)
call ufo_geovals_registry%get(c_key_geovals, geovals)

! Call method
call self%fill_geovals(geom, state, t1, t2, locs, geovals)

end subroutine ucldasv2_getvalues_fill_geovals_c


! --------------------------------------------------------------------------------------------------
!> C++ interface for ucldasv2_getvalues_mod::ucldasv2_getvalues::fill_geovals()
subroutine ucldasv2_getvalues_fill_geovals_tl_c(c_key_self, c_key_geom, c_key_incr, c_t1, c_t2, &
                                         c_locs, c_key_geovals) &
           bind (c, name='ucldasv2_getvalues_fill_geovals_tl_f90')

integer(c_int),     intent(in) :: c_key_self
integer(c_int),     intent(in) :: c_key_geom
integer(c_int),     intent(in) :: c_key_incr
type(c_ptr), value, intent(in) :: c_t1
type(c_ptr), value, intent(in) :: c_t2
type(c_ptr), value, intent(in) :: c_locs
integer(c_int),     intent(in) :: c_key_geovals

type(ucldasv2_getvalues), pointer :: self
type(ucldasv2_geom),      pointer :: geom
type(ucldasv2_increment), pointer :: incr
type(datetime)                :: t1
type(datetime)                :: t2
type(ufo_locations)           :: locs
type(ufo_geovals),    pointer :: geovals

! Get objects
call ucldasv2_getvalues_registry%get(c_key_self, self)
call ucldasv2_geom_registry%get(c_key_geom, geom)
call ucldasv2_increment_registry%get(c_key_incr, incr)
call c_f_datetime(c_t1, t1)
call c_f_datetime(c_t2, t2)
locs = ufo_locations(c_locs)
call ufo_geovals_registry%get(c_key_geovals, geovals)

! Call method
call self%fill_geovals(geom, incr, t1, t2, locs, geovals)

end subroutine ucldasv2_getvalues_fill_geovals_tl_c


! --------------------------------------------------------------------------------------------------
!> C++ interface for ucldasv2_getvalues_mod::ucldasv2_getvalues::fill_geovals_ad()
subroutine ucldasv2_getvalues_fill_geovals_ad_c(c_key_self, c_key_geom, c_key_incr, c_t1, c_t2, &
                                            c_locs, c_key_geovals) &
           bind (c, name='ucldasv2_getvalues_fill_geovals_ad_f90')

integer(c_int),     intent(in) :: c_key_self
integer(c_int),     intent(in) :: c_key_geom
integer(c_int),     intent(in) :: c_key_incr
type(c_ptr), value, intent(in) :: c_t1
type(c_ptr), value, intent(in) :: c_t2
type(c_ptr), value, intent(in) :: c_locs
integer(c_int),     intent(in) :: c_key_geovals

type(ucldasv2_getvalues), pointer :: self
type(ucldasv2_geom),      pointer :: geom
type(ucldasv2_increment), pointer :: incr
type(datetime)                :: t1
type(datetime)                :: t2
type(ufo_locations)           :: locs
type(ufo_geovals),    pointer :: geovals

! Get objects
call ucldasv2_getvalues_registry%get(c_key_self, self)
call ucldasv2_geom_registry%get(c_key_geom, geom)
call ucldasv2_increment_registry%get(c_key_incr, incr)
call c_f_datetime(c_t1, t1)
call c_f_datetime(c_t2, t2)
locs = ufo_locations(c_locs)
call ufo_geovals_registry%get(c_key_geovals, geovals)

! Call method
call self%fill_geovals_ad(geom, incr, t1, t2, locs, geovals)

end subroutine ucldasv2_getvalues_fill_geovals_ad_c

end module ucldasv2_getvalues_mod_c
