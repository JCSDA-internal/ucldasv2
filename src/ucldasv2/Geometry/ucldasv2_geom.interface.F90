! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence
! Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module ucldasv2_geom_mod_c

use atlas_module, only: atlas_fieldset, atlas_functionspace_pointcloud
use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module,           only: fckit_mpi_comm
use iso_c_binding
use oops_variables_mod, only: oops_variables

! ucldasv2 modules
use ucldasv2_geom_mod, only: ucldasv2_geom
use ucldasv2_fields_metadata_mod, only : ucldasv2_field_metadata


implicit none
private


#define LISTED_TYPE ucldasv2_geom

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry for ucldasv2_geom instances
type(registry_t), public :: ucldasv2_geom_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_geom_mod::ucldasv2_geom::init()
subroutine ucldasv2_geo_setup_c(c_key_self, c_conf, c_comm) bind(c,name='ucldasv2_geo_setup_f90')
  integer(c_int),  intent(inout) :: c_key_self
  type(c_ptr),        intent(in) :: c_conf
  type(c_ptr), value, intent(in) :: c_comm

  type(ucldasv2_geom), pointer :: self

  call ucldasv2_geom_registry%init()
  call ucldasv2_geom_registry%add(c_key_self)
  call ucldasv2_geom_registry%get(c_key_self,self)

  call self%init(fckit_configuration(c_conf), fckit_mpi_comm(c_comm) )
end subroutine ucldasv2_geo_setup_c


! --------------------------------------------------------------------------------------------------
!> C++ interface for ucldasv2_geom_mod::ucldasv2_geom::set_atlas_lonlat()
subroutine ucldasv2_geo_set_atlas_lonlat_c(c_key_self, c_afieldset) bind(c,name='ucldasv2_geo_set_atlas_lonlat_f90')
  integer(c_int), intent(in) :: c_key_self
  type(c_ptr), intent(in), value :: c_afieldset

  type(ucldasv2_geom), pointer :: self
  type(atlas_fieldset) :: afieldset

  call ucldasv2_geom_registry%get(c_key_self,self)
  afieldset = atlas_fieldset(c_afieldset)

  call self%set_atlas_lonlat(afieldset)
end subroutine ucldasv2_geo_set_atlas_lonlat_c


! --------------------------------------------------------------------------------------------------
!> C++ interface to get atlas functionspace pointr from ucldasv2_geom_mod::ucldasv2_geom
subroutine ucldasv2_geo_set_atlas_functionspace_pointer_c(c_key_self,c_afunctionspace) &
  bind(c,name='ucldasv2_geo_set_atlas_functionspace_pointer_f90')
  integer(c_int), intent(in)     :: c_key_self
  type(c_ptr), intent(in), value :: c_afunctionspace

  type(ucldasv2_geom),pointer :: self

  call ucldasv2_geom_registry%get(c_key_self,self)

  self%afunctionspace = atlas_functionspace_pointcloud(c_afunctionspace)
end subroutine ucldasv2_geo_set_atlas_functionspace_pointer_c


! --------------------------------------------------------------------------------------------------
!> C++ interface for ucldasv2_geom_mod::ucldasv2_geom::fill_atlas_fieldset()
subroutine ucldasv2_geo_fill_atlas_fieldset_c(c_key_self, c_afieldset) &
  bind(c,name='ucldasv2_geo_fill_atlas_fieldset_f90')

  integer(c_int),     intent(in) :: c_key_self
  type(c_ptr), value, intent(in) :: c_afieldset

  type(ucldasv2_geom), pointer :: self
  type(atlas_fieldset) :: afieldset

  call ucldasv2_geom_registry%get(c_key_self,self)
  afieldset = atlas_fieldset(c_afieldset)

  call self%fill_atlas_fieldset(afieldset)
end subroutine ucldasv2_geo_fill_atlas_fieldset_c

! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_geom_mod::ucldasv2_geom::clone()
subroutine ucldasv2_geo_clone_c(c_key_self, c_key_other) bind(c,name='ucldasv2_geo_clone_f90')
  integer(c_int), intent(inout) :: c_key_self
  integer(c_int), intent(in)    :: c_key_other

  type(ucldasv2_geom), pointer :: self, other

  call ucldasv2_geom_registry%add(c_key_self)
  call ucldasv2_geom_registry%get(c_key_self, self)
  call ucldasv2_geom_registry%get(c_key_other, other )

  call self%clone(other)
end subroutine ucldasv2_geo_clone_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_geom_mod::ucldasv2_geom::end()
subroutine ucldasv2_geo_delete_c(c_key_self) bind(c,name='ucldasv2_geo_delete_f90')
  integer(c_int), intent(inout) :: c_key_self

  type(ucldasv2_geom), pointer :: self

  call ucldasv2_geom_registry%get(c_key_self, self)
  call self%end()
  call ucldasv2_geom_registry%remove(c_key_self)
end subroutine ucldasv2_geo_delete_c


! ------------------------------------------------------------------------------
!> C++ interface to return begin and end of local geometry in ucldasv2_geom
subroutine ucldasv2_geo_start_end_c(c_key_self, ist, iend, jst, jend) bind(c, name='ucldasv2_geo_start_end_f90')
  integer(c_int), intent( in) :: c_key_self
  integer(c_int), intent(out) :: ist, iend, jst, jend

  type(ucldasv2_geom), pointer :: self
  call ucldasv2_geom_registry%get(c_key_self, self)

  ist  = self%isc
  iend = self%iec
  jst  = self%jsc
  jend = self%jec
end subroutine ucldasv2_geo_start_end_c


! ------------------------------------------------------------------------------
!> C++ interface to get number of levels for ucldasv2_geom
subroutine ucldasv2_geo_get_num_levels_c(c_key_self, c_vars, c_levels_size, c_levels) &
           bind(c, name='ucldasv2_geo_get_num_levels_f90')
  integer(c_int),     intent(in)  :: c_key_self
  type(c_ptr), value, intent(in)  :: c_vars
  integer(c_size_t),  intent(in)  :: c_levels_size
  integer(c_size_t),  intent(out) :: c_levels(c_levels_size)

  type(ucldasv2_geom), pointer :: self
  type(oops_variables)     :: vars
  integer :: i
  character(len=:), allocatable :: field_name
  type(ucldasv2_field_metadata) :: field

  call ucldasv2_geom_registry%get(c_key_self, self)
  vars = oops_variables(c_vars)

  do i = 1,vars%nvars()
    field_name = vars%variable(i)
    field = self%fields_metadata%get(field_name)
    select case(field%levels)
    case ("1")
      c_levels(i) = 1
    case ("full_lnd")
      if (field_name == field%getval_name_surface) then
        c_levels(i) = 1
      else
        c_levels(i) = self%nsoil
      end if
    case ("full_snw")
      if (field_name == field%getval_name_surface) then
        c_levels(i) = 1
      else
        c_levels(i) = self%nsnow
      end if
    case default
      call abor1_ftn('ERROR in ucldasv2_geo_get_num_levels_c, unknown "levels" '//field%levels)
    end select
  end do
end subroutine ucldasv2_geo_get_num_levels_c

! ------------------------------------------------------------------------------

end module ucldasv2_geom_mod_c
