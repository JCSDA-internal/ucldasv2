! (C) Copyright 2020-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> State fields
module ucldasv2_state_mod

use fckit_log_module, only: fckit_log
use kinds, only: kind_real
use oops_variables_mod

! ucldasv2 modules
use ucldasv2_geom_mod
use ucldasv2_fields_mod
!!use ucldasv2_increment_mod
!!use ucldasv2_convert_state_mod

implicit none
private


!-------------------------------------------------------------------------------
!> State fields.
!!
!! Any procedures that are shared with ucldasv2_increment are implemented
!! in the ucldasv2_fields base class
type, public, extends(ucldasv2_fields) :: ucldasv2_state

contains

  !> \name interactions with increment
  !! \{

  !> \copybrief ucldasv2_state_diff_incr \see ucldasv2_state_diff_incr
!  procedure :: diff_incr=> ucldasv2_state_diff_incr

  !> \copybrief ucldasv2_state_add_incr \see ucldasv2_state_add_incr
!  procedure :: add_incr => ucldasv2_state_add_incr

  !> \}


  !> \name misc
  !! \{

  !> \copybrief ucldasv2_state_convert \see ucldasv2_state_convert
  procedure :: convert => ucldasv2_state_convert

  !> \}

end type


!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
!> add a set of increments to the set of fields
!!
!! \throws abor1_ftn aborts if \p rhs is not a subset of \p self
!! \relates ucldasv2_state_mod::ucldasv2_state
!subroutine ucldasv2_state_add_incr(self, rhs)
!  class(ucldasv2_state),  intent(inout) :: self
! class(ucldasv2_increment), intent(in) :: rhs !< increment to add to \p self

! type(ucldasv2_field), pointer :: fld, fld_r
! integer :: i, k

!  real(kind=kind_real) :: min_ice = 1e-6_kind_real
!  real(kind=kind_real) :: amin = 1e-6_kind_real
!  real(kind=kind_real) :: amax = 10.0_kind_real
!  real(kind=kind_real), allocatable :: alpha(:,:), aice_bkg(:,:), aice_ana(:,:)
!  type(ucldasv2_fields) :: incr

  ! make sure rhs is a subset of self
! call rhs%check_subset(self)

  ! Make a copy of the increment
! call incr%copy(rhs)


  ! for each field that exists in incr, add to self
! do i=1,size(incr%fields)
!   fld_r => incr%fields(i)
!   call self%get(fld_r%name, fld)
!   fld%val = fld%val + fld_r%val
! end do

!end subroutine ucldasv2_state_add_incr


! ------------------------------------------------------------------------------
!> subtract two sets of fields, saving the results in \p inc
!!
!! \f$ inc = x1 - x2 \f$
!! \throws abor1_ftn aborts if \p inc and \p x2 are not subsets of \p x1
!! \relates ucldasv2_state_mod::ucldasv2_state
!subroutine ucldasv2_state_diff_incr(x1, x2, inc)
!  class(ucldasv2_state),      intent(in)    :: x1
!  class(ucldasv2_state),      intent(in)    :: x2
! class(ucldasv2_increment), intent(inout)  :: inc

! integer :: i
! type(ucldasv2_field), pointer :: f1, f2

  ! make sure fields correct shapes
! call inc%check_subset(x2)
! call x2%check_subset(x1)

  ! subtract
! do i=1,size(inc%fields)
!   call x1%get(inc%fields(i)%name, f1)
!   call x2%get(inc%fields(i)%name, f2)
!   inc%fields(i)%val = f1%val - f2%val
! end do
!end subroutine ucldasv2_state_diff_incr


! ------------------------------------------------------------------------------
!> Change resolution of \p rhs to \p self
!!
!! \p self must have valid "layer_depth" and "hocn" fields. The other fields
!! are interpolated from \p rhs to \p self. Any variables that are marked as
!! "positive definite" in the metadata configuration file are forced to be >= 0.0
!! after interpolation.
!!
!! \relates ucldasv2_state_mod::ucldasv2_state
subroutine ucldasv2_state_convert(self, rhs)
  class(ucldasv2_state), intent(inout) :: self
  class(ucldasv2_state), intent(in)    :: rhs   !< source

! integer :: n
! type(ucldasv2_convertstate_type) :: convert_state
! type(ucldasv2_field), pointer :: field1, field2, hocn1, hocn2, layer_depth

! call rhs%get("hocn", hocn1)
! call self%get("hocn", hocn2)
! call convert_state%setup(rhs%geom, self%geom, hocn1, hocn2)
! do n = 1, size(rhs%fields)
!   if (rhs%fields(n)%name=='layer_depth') cycle ! skip layer_depth interpolation
!   field1 => rhs%fields(n)
!   call self%get(trim(field1%name),field2)
!   if (field1%metadata%io_file=="ocn" .or. field1%metadata%io_file=="sfc" .or. field1%metadata%io_file=="ice")  &
!   call convert_state%change_resol(field1, field2, rhs%geom, self%geom)
!   ! Insure that positive definite variables are still >0
!   if (rhs%fields(n)%metadata%property=='positive_definite') then
!      where (field2%val<0.0)
!         field2%val=0.0
!      end where
!   end if
! end do !n

! ! Set layer depth for new grid
! call self%get("layer_depth", layer_depth)
! call self%geom%thickness2depth(hocn2%val, layer_depth%val)
! call convert_state%clean()
end subroutine ucldasv2_state_convert


! ------------------------------------------------------------------------------


end module
