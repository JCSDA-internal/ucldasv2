! (C) Copyright 2020-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> C++ interfaces for ucldasv2_increment_mod::ucldasv2_increment
module ucldasv2_increment_mod_c

use atlas_module, only: atlas_fieldset
use datetime_mod, only: datetime, c_f_datetime
use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use kinds, only: kind_real
use oops_variables_mod, only : oops_variables

! ucldasv2 modules
use ucldasv2_geom_iter_mod_c, only: ucldasv2_geom_iter_registry
use ucldasv2_geom_iter_mod, only : ucldasv2_geom_iter
use ucldasv2_geom_mod_c, only: ucldasv2_geom_registry
use ucldasv2_geom_mod, only: ucldasv2_geom
use ucldasv2_increment_mod, only : ucldasv2_increment
use ucldasv2_increment_reg, only : ucldasv2_increment_registry
use ucldasv2_state_mod, only : ucldasv2_state
use ucldasv2_state_reg, only : ucldasv2_state_registry

implicit none
private

contains


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_increment_mod::ucldasv2_increment version of
!! ucldasv2_fields_mod::ucldasv2_fields::create()
subroutine ucldasv2_increment_create_c(c_key_self, c_key_geom, c_vars) bind(c,name='ucldasv2_increment_create_f90')
  integer(c_int), intent(inout) :: c_key_self !< Handle to field
  integer(c_int),    intent(in) :: c_key_geom !< Geometry
  type(c_ptr),value, intent(in) :: c_vars     !< List of variables

  type(ucldasv2_increment),pointer :: self
  type(ucldasv2_geom),  pointer :: geom
  type(oops_variables)      :: vars

  call ucldasv2_geom_registry%get(c_key_geom, geom)
  call ucldasv2_increment_registry%init()
  call ucldasv2_increment_registry%add(c_key_self)
  call ucldasv2_increment_registry%get(c_key_self,self)

  vars = oops_variables(c_vars)
  call self%create(geom, vars)

end subroutine ucldasv2_increment_create_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_increment_mod::ucldasv2_increment version of
!! ucldasv2_fields_mod::ucldasv2_fields::delete()
subroutine ucldasv2_increment_delete_c(c_key_self) bind(c,name='ucldasv2_increment_delete_f90')
  integer(c_int), intent(inout) :: c_key_self

  type(ucldasv2_increment),    pointer :: self

  call ucldasv2_increment_registry%get(c_key_self,self)
  call self%delete( )
  call ucldasv2_increment_registry%remove(c_key_self)

end subroutine ucldasv2_increment_delete_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_increment_mod::ucldasv2_increment version of
!! ucldasv2_fields_mod::ucldasv2_fields::ones()
subroutine ucldasv2_increment_ones_c(c_key_self) bind(c,name='ucldasv2_increment_ones_f90')
  integer(c_int), intent(in) :: c_key_self

  type(ucldasv2_increment), pointer :: self

  call ucldasv2_increment_registry%get(c_key_self,self)
  call self%ones()

end subroutine ucldasv2_increment_ones_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_increment_mod::ucldasv2_increment version of
!! ucldasv2_fields_mod::ucldasv2_fields::zeros()
subroutine ucldasv2_increment_zero_c(c_key_self) bind(c,name='ucldasv2_increment_zero_f90')
  integer(c_int), intent(in) :: c_key_self

  type(ucldasv2_increment), pointer :: self

  call ucldasv2_increment_registry%get(c_key_self,self)
  call self%zeros()

end subroutine ucldasv2_increment_zero_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_increment_mod::ucldasv2_increment::dirac()
subroutine ucldasv2_increment_dirac_c(c_key_self,c_conf) bind(c,name='ucldasv2_increment_dirac_f90')
  integer(c_int), intent(in) :: c_key_self
  type(c_ptr),    intent(in) :: c_conf !< Configuration

  type(ucldasv2_increment), pointer :: self

  call ucldasv2_increment_registry%get(c_key_self,self)
  call self%dirac(fckit_configuration(c_conf))

end subroutine ucldasv2_increment_dirac_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_increment_mod::ucldasv2_increment::random()
subroutine ucldasv2_increment_random_c(c_key_self) bind(c,name='ucldasv2_increment_random_f90')
  integer(c_int), intent(in) :: c_key_self

  type(ucldasv2_increment), pointer :: self

  call ucldasv2_increment_registry%get(c_key_self,self)
  call self%random()

end subroutine ucldasv2_increment_random_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_increment_mod::ucldasv2_increment version of
!! ucldasv2_fields_mod::ucldasv2_fields::copy()
subroutine ucldasv2_increment_copy_c(c_key_self,c_key_rhs) bind(c,name='ucldasv2_increment_copy_f90')
  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(in) :: c_key_rhs

  type(ucldasv2_increment), pointer :: self
  type(ucldasv2_increment), pointer :: rhs

  call ucldasv2_increment_registry%get(c_key_self,self)
  call ucldasv2_increment_registry%get(c_key_rhs,rhs)

  call self%copy(rhs)

end subroutine ucldasv2_increment_copy_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_increment_mod::ucldasv2_increment version of
!! ucldasv2_fields_mod::ucldasv2_fields::add()
subroutine ucldasv2_increment_self_add_c(c_key_self,c_key_rhs) bind(c,name='ucldasv2_increment_self_add_f90')
  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(in) :: c_key_rhs

  type(ucldasv2_increment), pointer :: self
  type(ucldasv2_increment), pointer :: rhs

  call ucldasv2_increment_registry%get(c_key_self,self)
  call ucldasv2_increment_registry%get(c_key_rhs,rhs)

  call self%add(rhs)

end subroutine ucldasv2_increment_self_add_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_increment_mod::ucldasv2_increment::schur()
subroutine ucldasv2_increment_self_schur_c(c_key_self,c_key_rhs) bind(c,name='ucldasv2_increment_self_schur_f90')
  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(in) :: c_key_rhs

  type(ucldasv2_increment), pointer :: self
  type(ucldasv2_increment), pointer :: rhs

  call ucldasv2_increment_registry%get(c_key_self,self)
  call ucldasv2_increment_registry%get(c_key_rhs,rhs)

  call self%schur(rhs)

end subroutine ucldasv2_increment_self_schur_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_increment_mod::ucldasv2_increment version of
!! ucldasv2_fields_mod::ucldasv2_fields::sub()
subroutine ucldasv2_increment_self_sub_c(c_key_self,c_key_rhs) bind(c,name='ucldasv2_increment_self_sub_f90')
  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(in) :: c_key_rhs

  type(ucldasv2_increment), pointer :: self
  type(ucldasv2_increment), pointer :: rhs

  call ucldasv2_increment_registry%get(c_key_self,self)
  call ucldasv2_increment_registry%get(c_key_rhs,rhs)

  call self%sub(rhs)

end subroutine ucldasv2_increment_self_sub_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_increment_mod::ucldasv2_increment version of
!! ucldasv2_fields_mod::ucldasv2_fields::mul()
subroutine ucldasv2_increment_self_mul_c(c_key_self,c_zz) bind(c,name='ucldasv2_increment_self_mul_f90')
  integer(c_int), intent(in) :: c_key_self
  real(c_double), intent(in) :: c_zz

  type(ucldasv2_increment), pointer :: self
  real(kind=kind_real) :: zz

  call ucldasv2_increment_registry%get(c_key_self,self)
  zz = c_zz

  call self%mul(zz)

end subroutine ucldasv2_increment_self_mul_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_increment_mod::ucldasv2_increment version of
!! ucldasv2_fields_mod::ucldasv2_fields::axpy()
subroutine ucldasv2_increment_accumul_c(c_key_self,c_zz,c_key_rhs) bind(c,name='ucldasv2_increment_accumul_f90')
  integer(c_int), intent(in) :: c_key_self
  real(c_double), intent(in) :: c_zz
  integer(c_int), intent(in) :: c_key_rhs

  type(ucldasv2_increment), pointer :: self
  real(kind=kind_real)          :: zz
  type(ucldasv2_state),     pointer :: rhs

  call ucldasv2_increment_registry%get(c_key_self,self)
  call ucldasv2_state_registry%get(c_key_rhs,rhs)
  zz = c_zz

  call self%axpy(zz,rhs)

end subroutine ucldasv2_increment_accumul_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_increment_mod::ucldasv2_increment version of
!! ucldasv2_fields_mod::ucldasv2_fields::axpy()
subroutine ucldasv2_increment_axpy_c(c_key_self,c_zz,c_key_rhs) bind(c,name='ucldasv2_increment_axpy_f90')
  integer(c_int), intent(in) :: c_key_self
  real(c_double), intent(in) :: c_zz
  integer(c_int), intent(in) :: c_key_rhs

  type(ucldasv2_increment), pointer :: self
  real(kind=kind_real)      :: zz
  type(ucldasv2_increment), pointer :: rhs

  call ucldasv2_increment_registry%get(c_key_self,self)
  call ucldasv2_increment_registry%get(c_key_rhs,rhs)
  zz = c_zz

  call self%axpy(zz,rhs)

end subroutine ucldasv2_increment_axpy_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_increment_mod::ucldasv2_increment version of
!! ucldasv2_fields_mod::ucldasv2_fields::dot_prod()
subroutine ucldasv2_increment_dot_prod_c(c_key_fld1,c_key_fld2,c_prod) bind(c,name='ucldasv2_increment_dot_prod_f90')
  integer(c_int),    intent(in) :: c_key_fld1, c_key_fld2
  real(c_double), intent(inout) :: c_prod

  type(ucldasv2_increment), pointer :: fld1, fld2
  real(kind=kind_real) :: zz

  call ucldasv2_increment_registry%get(c_key_fld1,fld1)
  call ucldasv2_increment_registry%get(c_key_fld2,fld2)

  call fld1%dot_prod(fld2,zz)

  c_prod = zz

end subroutine ucldasv2_increment_dot_prod_c


! ------------------------------------------------------------------------------
!> C++ interface for subtracting two ucldasv2_increment_mod::ucldasv2_increment
!! using  ucldasv2_state_mod::ucldasv2_state::diff_incr()
subroutine ucldasv2_increment_diff_incr_c(c_key_lhs,c_key_x1,c_key_x2) bind(c,name='ucldasv2_increment_diff_incr_f90')
  integer(c_int), intent(in) :: c_key_lhs
  integer(c_int), intent(in) :: c_key_x1
  integer(c_int), intent(in) :: c_key_x2

  type(ucldasv2_increment), pointer :: lhs
  type(ucldasv2_state),     pointer :: x1
  type(ucldasv2_state),     pointer :: x2

  call ucldasv2_increment_registry%get(c_key_lhs,lhs)
  call ucldasv2_state_registry%get(c_key_x1,x1)
  call ucldasv2_state_registry%get(c_key_x2,x2)
  call x1%diff_incr(x2, lhs)

end subroutine ucldasv2_increment_diff_incr_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_increment_mod::ucldasv2_increment::convert()
subroutine ucldasv2_increment_change_resol_c(c_key_fld,c_key_rhs) bind(c,name='ucldasv2_increment_change_resol_f90')
  integer(c_int), intent(in) :: c_key_fld
  integer(c_int), intent(in) :: c_key_rhs

  type(ucldasv2_increment), pointer :: fld, rhs

  call ucldasv2_increment_registry%get(c_key_fld,fld)
  call ucldasv2_increment_registry%get(c_key_rhs,rhs)

  ! TODO (Guillaume or Travis) implement == in geometry or something to that effect.
  if ( size(fld%geom%lon,1)==size(rhs%geom%lon,1) .and. &
        size(fld%geom%lat,2)==size(rhs%geom%lat,2) .and. &
        fld%geom%nsnow==rhs%geom%nsnow ) then
      call fld%copy(rhs)
  else
      call fld%convert(rhs)
  endif

end subroutine ucldasv2_increment_change_resol_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_increment_mod::ucldasv2_increment::set_atlas()
subroutine ucldasv2_increment_set_atlas_c(c_key_self,c_key_geom,c_vars,c_afieldset) &
  & bind (c,name='ucldasv2_increment_set_atlas_f90')

  implicit none
  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(in) :: c_key_geom
  type(c_ptr), value, intent(in) :: c_vars
  type(c_ptr), intent(in), value :: c_afieldset

  type(ucldasv2_increment), pointer :: self
  type(ucldasv2_geom),  pointer :: geom
  type(oops_variables) :: vars
  type(atlas_fieldset) :: afieldset

  call ucldasv2_increment_registry%get(c_key_self,self)
  call ucldasv2_geom_registry%get(c_key_geom, geom)
  vars = oops_variables(c_vars)
  afieldset = atlas_fieldset(c_afieldset)

  call self%set_atlas(geom, vars, afieldset)

end subroutine ucldasv2_increment_set_atlas_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_increment_mod::ucldasv2_increment::to_atlas()
subroutine ucldasv2_increment_to_atlas_c(c_key_self,c_key_geom,c_vars,c_afieldset) &
  & bind (c,name='ucldasv2_increment_to_atlas_f90')

  implicit none
  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(in) :: c_key_geom
  type(c_ptr), value, intent(in) :: c_vars
  type(c_ptr), intent(in), value :: c_afieldset

  type(ucldasv2_increment), pointer :: self
  type(ucldasv2_geom),  pointer :: geom
  type(oops_variables) :: vars
  type(atlas_fieldset) :: afieldset

  call ucldasv2_increment_registry%get(c_key_self,self)
  call ucldasv2_geom_registry%get(c_key_geom, geom)
  vars = oops_variables(c_vars)
  afieldset = atlas_fieldset(c_afieldset)

  call self%to_atlas(geom, vars, afieldset)

end subroutine ucldasv2_increment_to_atlas_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_increment_mod::ucldasv2_increment::getpoint()
subroutine ucldasv2_increment_from_atlas_c(c_key_self,c_key_geom,c_vars,c_afieldset) &
  & bind (c,name='ucldasv2_increment_from_atlas_f90')

  implicit none
  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(in) :: c_key_geom
  type(c_ptr), value, intent(in) :: c_vars
  type(c_ptr), intent(in), value :: c_afieldset

  type(ucldasv2_increment), pointer :: self
  type(ucldasv2_geom),  pointer :: geom
  type(oops_variables) :: vars
  type(atlas_fieldset) :: afieldset

  call ucldasv2_increment_registry%get(c_key_self, self)
  call ucldasv2_geom_registry%get(c_key_geom, geom)
  vars = oops_variables(c_vars)
  afieldset = atlas_fieldset(c_afieldset)

  call self%from_atlas(geom, vars, afieldset)

end subroutine ucldasv2_increment_from_atlas_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_increment_mod::ucldasv2_increment version of
!! ucldasv2_fields_mod::ucldasv2_fields::read()
subroutine ucldasv2_increment_read_file_c(c_key_fld, c_conf, c_dt) bind(c,name='ucldasv2_increment_read_file_f90')
  integer(c_int), intent(in) :: c_key_fld  !< Fields
  type(c_ptr),    intent(in) :: c_conf     !< Configuration
  type(c_ptr), intent(inout) :: c_dt       !< DateTime

  type(ucldasv2_increment), pointer :: fld
  type(datetime)            :: fdate

  call ucldasv2_increment_registry%get(c_key_fld,fld)
  call c_f_datetime(c_dt, fdate)
  call fld%read(fckit_configuration(c_conf), fdate)

end subroutine ucldasv2_increment_read_file_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_increment_mod::ucldasv2_increment version of
!! ucldasv2_fields_mod::ucldasv2_fields::write()
subroutine ucldasv2_increment_write_file_c(c_key_fld, c_conf, c_dt) bind(c,name='ucldasv2_increment_write_file_f90')
  integer(c_int), intent(in) :: c_key_fld  !< Fields
  type(c_ptr),    intent(in) :: c_conf     !< Configuration
  type(c_ptr),    intent(in) :: c_dt       !< DateTime

  type(ucldasv2_increment), pointer :: fld
  type(datetime)            :: fdate

  call ucldasv2_increment_registry%get(c_key_fld,fld)
  call c_f_datetime(c_dt, fdate)
  call fld%write_rst(fckit_configuration(c_conf), fdate)

end subroutine ucldasv2_increment_write_file_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_increment_mod::ucldasv2_increment version of
!! ucldasv2_fields_mod::ucldasv2_fields::gpnorm()
subroutine ucldasv2_increment_gpnorm_c(c_key_fld, kf, pstat) bind(c,name='ucldasv2_increment_gpnorm_f90')
  integer(c_int),    intent(in) :: c_key_fld
  integer(c_int),    intent(in) :: kf
  real(c_double), intent(inout) :: pstat(3*kf)

  type(ucldasv2_increment), pointer :: fld
  real(kind=kind_real)      :: zstat(3, kf)
  integer :: jj, js, jf

  call ucldasv2_increment_registry%get(c_key_fld,fld)

  call fld%gpnorm(kf, zstat)
  jj=0
  do jf = 1, kf
      do js = 1, 3
        jj=jj+1
        pstat(jj) = zstat(js,jf)
      enddo
  enddo

end subroutine ucldasv2_increment_gpnorm_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_increment_mod::ucldasv2_increment RMS
subroutine ucldasv2_increment_rms_c(c_key_fld, prms) bind(c,name='ucldasv2_increment_rms_f90')
  integer(c_int),    intent(in) :: c_key_fld
  real(c_double), intent(inout) :: prms

  type(ucldasv2_increment), pointer :: fld
  real(kind=kind_real)      :: zz

  call ucldasv2_increment_registry%get(c_key_fld,fld)

  call fld%dot_prod(fld, zz)
  prms = sqrt(zz)

end subroutine ucldasv2_increment_rms_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_increment_mod::ucldasv2_increment::getpoint()
subroutine ucldasv2_increment_getpoint_c(c_key_fld,c_key_iter,values, values_len) bind(c,name='ucldasv2_increment_getpoint_f90')
  integer(c_int), intent(in) :: c_key_fld
  integer(c_int), intent(in) :: c_key_iter
  integer(c_int), intent(in) :: values_len
  real(c_double), intent(inout) :: values(values_len)

  type(ucldasv2_increment),      pointer :: fld
  type(ucldasv2_geom_iter), pointer :: iter

  call ucldasv2_increment_registry%get(c_key_fld,fld)
  call ucldasv2_geom_iter_registry%get(c_key_iter,iter)

  call fld%getpoint(iter, values)

end subroutine ucldasv2_increment_getpoint_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_increment_mod::ucldasv2_increment::setpoint()
subroutine ucldasv2_increment_setpoint_c(c_key_fld,c_key_iter,values, values_len) bind(c,name='ucldasv2_increment_setpoint_f90')
  integer(c_int), intent(inout) :: c_key_fld
  integer(c_int), intent(in) :: c_key_iter
  integer(c_int), intent(in) :: values_len
  real(c_double), intent(in) :: values(values_len)

  type(ucldasv2_increment),      pointer :: fld
  type(ucldasv2_geom_iter), pointer :: iter

  call ucldasv2_increment_registry%get(c_key_fld,fld)
  call ucldasv2_geom_iter_registry%get(c_key_iter,iter)

  call fld%setpoint(iter, values)

end subroutine ucldasv2_increment_setpoint_c


! ------------------------------------------------------------------------------
!> C++ interface to get ucldasv2_increment_mod::ucldasv2_increment dimension sizes
subroutine ucldasv2_incrementnum_c(c_key_fld, nx, ny, nz, nf) bind(c,name='ucldasv2_increment_sizes_f90')
  integer(c_int),         intent(in) :: c_key_fld
  integer(kind=c_int), intent(inout) :: nx, ny, nz, nf

  type(ucldasv2_increment), pointer :: fld

  call ucldasv2_increment_registry%get(c_key_fld,fld)

  nx = size(fld%geom%lon,1)
  ny = size(fld%geom%lon,2)
  nz = fld%geom%nsnow
  nf = size(fld%fields)

end subroutine ucldasv2_incrementnum_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_increment_mod::ucldasv2_increment version of
!! ucldasv2_fields_mod::ucldasv2_fields::serial_size()
subroutine ucldasv2_increment_serial_size_c(c_key_self,c_key_geom,c_vec_size) bind (c,name='ucldasv2_increment_serial_size_f90')

  implicit none
  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(in) :: c_key_geom
  integer(c_size_t), intent(out) :: c_vec_size

  type(ucldasv2_increment), pointer :: self
  type(ucldasv2_geom),  pointer :: geom

  integer :: vec_size

  call ucldasv2_increment_registry%get(c_key_self,self)
  call ucldasv2_geom_registry%get(c_key_geom,geom)

  call self%serial_size(geom,vec_size)
  c_vec_size = vec_size

end subroutine ucldasv2_increment_serial_size_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_increment_mod::ucldasv2_increment version of
!! ucldasv2_fields_mod::ucldasv2_fields::serialize()
subroutine ucldasv2_increment_serialize_c(c_key_self,c_key_geom,c_vec_size,c_vec) bind (c,name='ucldasv2_increment_serialize_f90')

  implicit none
  integer(c_int),    intent(in) :: c_key_self
  integer(c_int),    intent(in) :: c_key_geom
  integer(c_size_t), intent(in) :: c_vec_size
  real(c_double),   intent(out) :: c_vec(c_vec_size)

  integer :: vec_size
  type(ucldasv2_increment), pointer :: self
  type(ucldasv2_geom),  pointer :: geom

  vec_size = c_vec_size
  call ucldasv2_increment_registry%get(c_key_self,self)
  call ucldasv2_geom_registry%get(c_key_geom,geom)

  call self%serialize(geom, vec_size, c_vec)

end subroutine ucldasv2_increment_serialize_c


! ------------------------------------------------------------------------------
!> C++ interface for ucldasv2_increment_mod::ucldasv2_increment version of
!! ucldasv2_fields_mod::ucldasv2_fields::deserialize()
subroutine ucldasv2_increment_deserialize_c(c_key_self,c_key_geom,c_vec_size,c_vec,c_index) bind (c,name='ucldasv2_increment_deserialize_f90')

  implicit none
  integer(c_int),    intent(in) :: c_key_self
  integer(c_int),    intent(in) :: c_key_geom
  integer(c_size_t), intent(in) :: c_vec_size
  real(c_double),    intent(in) :: c_vec(c_vec_size)
  integer(c_size_t), intent(inout) :: c_index

  integer :: vec_size, idx
  type(ucldasv2_increment), pointer :: self
  type(ucldasv2_geom),  pointer :: geom

  vec_size = c_vec_size
  idx = c_index
  call ucldasv2_increment_registry%get(c_key_self,self)
  call ucldasv2_geom_registry%get(c_key_geom,geom)

  call self%deserialize(geom, vec_size, c_vec, idx)
  c_index=idx

end subroutine ucldasv2_increment_deserialize_c

end module
