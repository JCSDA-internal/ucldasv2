!
! (C) Copyright 2019-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence
! Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module ucldasv2_geom_iter_interface

  use iso_c_binding
  use kinds
  use ucldasv2_geom_iter_mod
  use ucldasv2_geom_mod_c,  only : ucldasv2_geom_registry
  use ucldasv2_geom_mod, only: ucldasv2_geom

  implicit none

  private

contains

  ! ------------------------------------------------------------------------------
  !> Setup geometry iterator
  subroutine ucldasv2_geom_iter_setup_c(c_key_self, c_key_geom, c_iindex, c_jindex) bind(c, name='ucldasv2_geom_iter_setup_f90')

    ! Passed variables
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
    call ucldasv2_geom_iter_setup(self, geom, c_iindex, c_jindex)

  end subroutine ucldasv2_geom_iter_setup_c

  ! ------------------------------------------------------------------------------
  !> Clone geometry iterator
  subroutine ucldasv2_geom_iter_clone_c(c_key_self, c_key_other) bind(c, name='ucldasv2_geom_iter_clone_f90')

    ! Passed variables
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
    call ucldasv2_geom_iter_clone(self, other)

  end subroutine ucldasv2_geom_iter_clone_c

  ! ------------------------------------------------------------------------------
  !> Delete geometry iterator
  subroutine ucldasv2_geom_iter_delete_c(c_key_self) bind(c, name='ucldasv2_geom_iter_delete_f90')

      ! Passed variables
      integer(c_int), intent(inout) :: c_key_self !< Geometry iterator

      ! Clear interface
      call ucldasv2_geom_iter_registry%remove(c_key_self)

  end subroutine ucldasv2_geom_iter_delete_c

end module ucldasv2_geom_iter_interface
