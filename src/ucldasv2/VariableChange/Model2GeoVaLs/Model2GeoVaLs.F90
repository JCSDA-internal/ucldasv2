! (C) Copyright 2020-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> C++ interface for converting model variables to geovals (mostly identity function)
module ucldasv2_model2geovals_mod_c

use iso_c_binding
use kinds, only: kind_real

! ucldasv2 modules
use ucldasv2_fields_mod, only: ucldasv2_field
use ucldasv2_geom_mod_c, only: ucldasv2_geom_registry
use ucldasv2_geom_mod, only: ucldasv2_geom
use ucldasv2_state_mod, only: ucldasv2_state
use ucldasv2_state_reg, only: ucldasv2_state_registry

implicit none
private


!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!> C++ interface for the non-linear change of variables from model to geovals
!!
!! This is *mostly* an identity operator, except for a small number of derived variables
!! that are to be calculated here ("distance_from_coast", "sea_area_fraction", etc.)
!! \throws abor1_ftn aborts if field name is not handled.
subroutine ucldasv2_model2geovals_changevar_f90(c_key_geom, c_key_xin, c_key_xout) &
  bind(c,name='ucldasv2_model2geovals_changevar_f90')
  integer(c_int), intent(in) :: c_key_geom, c_key_xin, c_key_xout

  type(ucldasv2_geom),  pointer :: geom
  type(ucldasv2_state), pointer :: xin, xout
  type(ucldasv2_field), pointer :: field
  integer :: i

  call ucldasv2_geom_registry%get(c_key_geom, geom)
  call ucldasv2_state_registry%get(c_key_xin, xin)
  call ucldasv2_state_registry%get(c_key_xout, xout)
!
  do i=1, size(xout%fields)

    ! identity operators
   
    write(*,*) (xout%fields(i)%name)

    ! special cases
    select case (xout%fields(i)%name)
       
      ! identity operators
      case ('snowd')
        call xin%get(xout%fields(i)%metadata%name, field)
        if (xout%fields(i)%name == field%metadata%name .or. &
            xout%fields(i)%name == field%metadata%getval_name ) then
          xout%fields(i)%val(:,:,:) =  field%val(:,:,:) !< full field
        elseif (field%metadata%getval_name_surface == xout%fields(i)%name) then
          xout%fields(i)%val(:,:,1) = field%val(:,:,1) !< surface only of a 3D field
        else
          call abor1_ftn( 'error in ucldasv2_model2geovals_changevar_f90 processing ' &
                          // xout%fields(i)%name )
        endif

      end select
  
  end do
end subroutine

!-------------------------------------------------------------------------------

end module
