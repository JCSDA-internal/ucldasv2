module mymodel_geometry_mod

use fckit_configuration_module, only: fckit_configuration

implicit none
private

public :: mymodel_geometry

!------------------------------------------------------------------------------

type :: mymodel_geometry
contains
  procedure :: init   => mymodel_geometry_init
  procedure :: clone  => mymodel_geometry_clone
  procedure :: delete => mymodel_geometry_delete
end type

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine mymodel_geometry_init(self, f_conf)
  class(mymodel_geometry),   intent(out) :: self
  type(fckit_configuration),  intent(in) :: f_conf

  call abor1_ftn("ERROR: mymodel_geometry_init() needs to be implemented.")
end subroutine

!------------------------------------------------------------------------------

subroutine mymodel_geometry_clone(self, other)
  class(mymodel_geometry),  intent(out) :: self
  class(mymodel_geometry),   intent(in) :: other
      
  call abor1_ftn("ERROR: mymodel_geometry_clone() needs to be implemented.")
end subroutine
  
!------------------------------------------------------------------------------

subroutine mymodel_geometry_delete(self)
  class(mymodel_geometry),  intent(out) :: self
       
  call abor1_ftn("ERROR: mymodel_geometry_delete() needs to be implemented.")
end subroutine

!------------------------------------------------------------------------------

end module