! (C) Copyright 2019-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module mymodel_geometry_mod_c

use iso_c_binding

use fckit_configuration_module, only: fckit_configuration
use mymodel_geometry_mod, only: mymodel_geometry

implicit none
private

! Setup the C/Fortran interface registry
#define LISTED_TYPE mymodel_geometry
#include "oops/util/linkedList_i.f"
type(registry_t) :: mymodel_geometry_registry

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

! Setup the C/Fortran interface registry
#include "oops/util/linkedList_c.f"

! -----------------------------------------------------------------------------

subroutine c_mymodel_geometry_setup(c_key_self, c_conf) bind(c,name='mymodel_geometry_setup_f90')
  integer(c_int), intent(inout) :: c_key_self
  type(c_ptr),       intent(in) :: c_conf
  
  type(mymodel_geometry), pointer :: self
  
  call mymodel_geometry_registry%init()
  call mymodel_geometry_registry%add(c_key_self)
  call mymodel_geometry_registry%get(c_key_self, self)  
  call self%init(fckit_configuration(c_conf))
  
end subroutine c_mymodel_geometry_setup

! -----------------------------------------------------------------------------

subroutine c_mymodel_geometry_clone(c_key_self, c_key_other) bind(c,name='mymodel_geometry_clone_f90')

  integer(c_int), intent(inout) :: c_key_self
  integer(c_int), intent(in)    :: c_key_other

  type(mymodel_geometry), pointer :: self, other
  
  call mymodel_geometry_registry%add(c_key_self)
  call mymodel_geometry_registry%get(c_key_self , self )
  call mymodel_geometry_registry%get(c_key_other, other)
  call self%clone(other)

end subroutine c_mymodel_geometry_clone

! -----------------------------------------------------------------------------

subroutine c_mymodel_geometry_delete(c_key_self) bind(c,name='mymodel_geometry_delete_f90')
  integer(c_int), intent(inout) :: c_key_self

  type(mymodel_geometry), pointer :: self
  
  call mymodel_geometry_registry%get(c_key_self , self )
  call self%delete()
  call mymodel_geometry_registry%remove(c_key_self)

end subroutine c_mymodel_geometry_delete

!------------------------------------------------------------------------------

end module