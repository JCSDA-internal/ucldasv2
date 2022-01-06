! (C) Copyright 2020-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> C++ interfaces for ucldasv2_state_mod::ucldasv2_state
module ucldasv2_state_mod_c

use datetime_mod, only: datetime, c_f_datetime
use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use kinds, only: kind_real
use oops_variables_mod, only: oops_variables

! ucldasv2 modules
use ucldasv2_geom_mod_c, only: ucldasv2_geom_registry
use ucldasv2_geom_mod, only: ucldasv2_geom
use ucldasv2_increment_mod, only: ucldasv2_increment
use ucldasv2_increment_reg, only: ucldasv2_increment_registry
use ucldasv2_state_mod, only: ucldasv2_state
use ucldasv2_state_reg, only: ucldasv2_state_registry

implicit none
private


! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_state_mod::ucldasv2_state version of
!! ucldasv2_fields_mod::ucldasv2_fields::create()
subroutine ucldasv2_state_create_c(c_key_self, c_key_geom, c_vars) bind(c,name='ucldasv2_state_create_f90')
    integer(c_int), intent(inout) :: c_key_self !< Handle to field
    integer(c_int),    intent(in) :: c_key_geom !< Geometry
    type(c_ptr),value, intent(in) :: c_vars     !< List of variables

    type(ucldasv2_state), pointer :: self
    type(ucldasv2_geom),  pointer :: geom
    type(oops_variables)          :: vars

    call ucldasv2_geom_registry%get(c_key_geom, geom)
    call ucldasv2_state_registry%init()
    call ucldasv2_state_registry%add(c_key_self)
    call ucldasv2_state_registry%get(c_key_self,self)

    vars = oops_variables(c_vars)
    call self%create(geom, vars)

end subroutine ucldasv2_state_create_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_state_mod::ucldasv2_state version of
!! ucldasv2_fields_mod::ucldasv2_fields::delete()
subroutine ucldasv2_state_delete_c(c_key_self) bind(c,name='ucldasv2_state_delete_f90')
    integer(c_int), intent(inout) :: c_key_self

    type(ucldasv2_state),    pointer :: self

    call ucldasv2_state_registry%get(c_key_self,self)
    call self%delete( )
    call ucldasv2_state_registry%remove(c_key_self)

end subroutine ucldasv2_state_delete_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_state_mod::ucldasv2_state version of
!! ucldasv2_fields_mod::ucldasv2_fields::zeros()
subroutine ucldasv2_state_zero_c(c_key_self) bind(c,name='ucldasv2_state_zero_f90')
    integer(c_int), intent(in) :: c_key_self

    type(ucldasv2_state), pointer :: self

    call ucldasv2_state_registry%get(c_key_self,self)
    call self%zeros()

end subroutine ucldasv2_state_zero_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_state_mod::ucldasv2_state version of
!! ucldasv2_fields_mod::ucldasv2_fields::copy()
subroutine ucldasv2_state_copy_c(c_key_self,c_key_rhs) bind(c,name='ucldasv2_state_copy_f90')
    integer(c_int), intent(in) :: c_key_self
    integer(c_int), intent(in) :: c_key_rhs

    type(ucldasv2_state), pointer :: self
    type(ucldasv2_state), pointer :: rhs

    call ucldasv2_state_registry%get(c_key_self,self)
    call ucldasv2_state_registry%get(c_key_rhs,rhs)

    call self%copy(rhs)

end subroutine ucldasv2_state_copy_c

! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_state_mod::ucldasv2_state version of
!! ucldasv2_fields_mod::ucldasv2_fields::axpy()
subroutine ucldasv2_state_axpy_c(c_key_self,c_zz,c_key_rhs) bind(c,name='ucldasv2_state_axpy_f90')
    integer(c_int), intent(in) :: c_key_self
    real(c_double), intent(in) :: c_zz
    integer(c_int), intent(in) :: c_key_rhs

    type(ucldasv2_state), pointer :: self
    real(kind=kind_real)          :: zz
    type(ucldasv2_state), pointer :: rhs

    call ucldasv2_state_registry%get(c_key_self,self)
    call ucldasv2_state_registry%get(c_key_rhs,rhs)
    zz = c_zz

    call self%axpy(zz,rhs)

end subroutine ucldasv2_state_axpy_c

! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_state_mod::ucldasv2_state::add_incr()
subroutine ucldasv2_state_add_incr_c(c_key_self,c_key_rhs) bind(c,name='ucldasv2_state_add_incr_f90')
    integer(c_int), intent(in) :: c_key_self
    integer(c_int), intent(in) :: c_key_rhs

    type(ucldasv2_state),     pointer :: self
    type(ucldasv2_increment), pointer :: rhs

    call ucldasv2_state_registry%get(c_key_self,self)
    call ucldasv2_increment_registry%get(c_key_rhs,rhs)

    call self%add_incr(rhs)

end subroutine ucldasv2_state_add_incr_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_state_mod::ucldasv2_state version of
!! ucldasv2_fields_mod::ucldasv2_fields::read()
subroutine ucldasv2_state_read_file_c(c_key_fld, c_conf, c_dt) bind(c,name='ucldasv2_state_read_file_f90')
    integer(c_int), intent(in) :: c_key_fld  !< Fields
    type(c_ptr),    intent(in) :: c_conf     !< Configuration
    type(c_ptr), intent(inout) :: c_dt       !< DateTime

    type(ucldasv2_state), pointer :: fld
    type(datetime)                :: fdate

    call ucldasv2_state_registry%get(c_key_fld,fld)
    call c_f_datetime(c_dt, fdate)
    call fld%read(fckit_configuration(c_conf), fdate)

end subroutine ucldasv2_state_read_file_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_state_mod::ucldasv2_state version of
!! ucldasv2_fields_mod::ucldasv2_fields::write_rst()
subroutine ucldasv2_state_write_file_c(c_key_fld, c_conf, c_dt) bind(c,name='ucldasv2_state_write_file_f90')
    integer(c_int), intent(in) :: c_key_fld  !< Fields
    type(c_ptr),    intent(in) :: c_conf     !< Configuration
    type(c_ptr),    intent(in) :: c_dt       !< DateTime

    type(ucldasv2_state), pointer :: fld
    type(datetime)                :: fdate

    call ucldasv2_state_registry%get(c_key_fld,fld)
    call c_f_datetime(c_dt, fdate)
    call fld%write_rst(fckit_configuration(c_conf), fdate)

end subroutine ucldasv2_state_write_file_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_state_mod::ucldasv2_state version of
!! ucldasv2_fields_mod::ucldasv2_fields::gpnorm()
subroutine ucldasv2_state_gpnorm_c(c_key_fld, kf, pstat) bind(c,name='ucldasv2_state_gpnorm_f90')
    integer(c_int),    intent(in) :: c_key_fld
    integer(c_int),    intent(in) :: kf
    real(c_double), intent(inout) :: pstat(3*kf)

    type(ucldasv2_state), pointer :: fld
    real(kind=kind_real)          :: zstat(3, kf)
    integer :: jj, js, jf

    call ucldasv2_state_registry%get(c_key_fld,fld)

    call fld%gpnorm(kf, zstat)
    jj=0
    do jf = 1, kf
        do js = 1, 3
        jj=jj+1
        pstat(jj) = zstat(js,jf)
        enddo
    enddo

end subroutine ucldasv2_state_gpnorm_c

! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_state_mod::ucldasv2_state RMS
subroutine ucldasv2_state_rms_c(c_key_fld, prms) bind(c,name='ucldasv2_state_rms_f90')
    integer(c_int),    intent(in) :: c_key_fld
    real(c_double), intent(inout) :: prms

    type(ucldasv2_state), pointer :: fld
    real(kind=kind_real)          :: zz

    call ucldasv2_state_registry%get(c_key_fld,fld)

    call fld%dot_prod(fld, zz)
    prms = sqrt(zz)

end subroutine ucldasv2_state_rms_c


! ------------------------------------------------------------------------------
!> C++ interface to get ucldasv2_state_mod::ucldasv2_state dimensions sizes
subroutine ucldasv2_state_sizes_c(c_key_fld, nx, ny, nz, nf) bind(c,name='ucldasv2_state_sizes_f90')
    integer(c_int),         intent(in) :: c_key_fld
    integer(kind=c_int), intent(inout) :: nx, ny, nz, nf

    type(ucldasv2_state), pointer :: fld

    call ucldasv2_state_registry%get(c_key_fld,fld)

    nx = size(fld%geom%lon,1)
    ny = size(fld%geom%lon,2)
    nz = fld%geom%nsnow
    nf = size(fld%fields)

end subroutine ucldasv2_state_sizes_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_state_mod::ucldasv2_state::convert()
subroutine ucldasv2_state_change_resol_c(c_key_fld,c_key_rhs) bind(c,name='ucldasv2_state_change_resol_f90')
    integer(c_int), intent(in) :: c_key_fld
    integer(c_int), intent(in) :: c_key_rhs

    type(ucldasv2_state), pointer :: fld, rhs

    call ucldasv2_state_registry%get(c_key_fld,fld)
    call ucldasv2_state_registry%get(c_key_rhs,rhs)

    ! TODO (Guillaume or Travis) implement == in geometry or something to that
    ! effect.
     if (size(fld%geom%lon,1)==size(rhs%geom%lon,1) .and. size(fld%geom%lat,2)==size(rhs%geom%lat,2) .and. &
       fld%geom%nsnow==rhs%geom%nsnow ) then
       call fld%copy(rhs)
     else
      call fld%convert(rhs)
     endif

end subroutine ucldasv2_state_change_resol_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_state_mod::ucldasv2_state version of
!! ucldasv2_fields_mod::ucldasv2_fields::serial_size()
subroutine ucldasv2_state_serial_size_c(c_key_self,c_key_geom,c_vec_size) bind (c,name='ucldasv2_state_serial_size_f90')

  implicit none
  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(in) :: c_key_geom
  integer(c_size_t), intent(out) :: c_vec_size

  type(ucldasv2_state), pointer :: self
  type(ucldasv2_geom),  pointer :: geom
  integer :: vec_size

  call ucldasv2_state_registry%get(c_key_self,self)
  call ucldasv2_geom_registry%get(c_key_geom,geom)

  call self%serial_size(geom, vec_size)
  c_vec_size = vec_size

end subroutine ucldasv2_state_serial_size_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_state_mod::ucldasv2_state version of
!! ucldasv2_fields_mod::ucldasv2_fields::serialize()
subroutine ucldasv2_state_serialize_c(c_key_self,c_key_geom,c_vec_size,c_vec) bind (c,name='ucldasv2_state_serialize_f90')

  implicit none
  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(in) :: c_key_geom
  integer(c_size_t), intent(in) :: c_vec_size
  real(c_double), intent(out) :: c_vec(c_vec_size)

  type(ucldasv2_state), pointer :: self
  type(ucldasv2_geom),  pointer :: geom

  integer :: vec_size

  vec_size = c_vec_size
  call ucldasv2_state_registry%get(c_key_self,self)
  call ucldasv2_geom_registry%get(c_key_geom,geom)

  call self%serialize(geom, vec_size, c_vec)

end subroutine ucldasv2_state_serialize_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_state_mod::ucldasv2_state version of
!! ucldasv2_fields_mod::ucldasv2_fields::deserialize()
subroutine ucldasv2_state_deserialize_c(c_key_self,c_key_geom,c_vec_size,c_vec,c_index) bind (c,name='ucldasv2_state_deserialize_f90')

  implicit none
  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(in) :: c_key_geom
  integer(c_size_t), intent(in) :: c_vec_size
  real(c_double), intent(in) :: c_vec(c_vec_size)
  integer(c_size_t), intent(inout) :: c_index

  type(ucldasv2_state), pointer :: self
  type(ucldasv2_geom),  pointer :: geom
  integer :: vec_size, idx

  vec_size = c_vec_size
  idx = c_index
  call ucldasv2_state_registry%get(c_key_self,self)
  call ucldasv2_geom_registry%get(c_key_geom,geom)

  call self%deserialize(geom,vec_size,c_vec, idx)
  c_index=idx

end subroutine ucldasv2_state_deserialize_c

! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_state_mod::ucldasv2_state version of
!! ucldasv2_fields_mod::ucldasv2_fields::update_fields()
subroutine ucldasv2_state_update_fields_c(c_key_self, c_vars) &
           bind (c,name='ucldasv2_state_update_fields_f90')

integer(c_int),     intent(in) :: c_key_self
type(c_ptr), value, intent(in) :: c_vars

type(ucldasv2_state), pointer :: f_self
type(oops_variables)      :: f_vars

! LinkedList
! ----------
call ucldasv2_state_registry%get(c_key_self, f_self)

! Fortrain APIs
! -------------
f_vars = oops_variables(c_vars)

! Call implementation
! -------------------
call f_self%update_fields(f_vars)

end subroutine ucldasv2_state_update_fields_c

end module ucldasv2_state_mod_c
