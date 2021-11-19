! (C) Copyright 2017-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.


!> Handle fields for the model.
!!
!! ucldasv2_fields represents a state or increment, and contains one or more
!! ucldasv2_field instances for each of the fields. The metadata associated
!! with a given field is stored in ucldasv2_fields_metadata_mod::ucldasv2_fields_metadata
module ucldasv2_fields_mod

! JEDI modules
use datetime_mod, only: datetime, datetime_set, datetime_to_string, &
                        datetime_create, datetime_diff
use duration_mod, only: duration, duration_to_string
use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module, only: fckit_mpi_min, fckit_mpi_max, fckit_mpi_sum
use kinds, only: kind_real
use oops_variables_mod, only: oops_variables
use tools_const, only: deg2rad

! FMS modules
use fms_io_mod, only: fms_io_init, fms_io_exit, register_restart_field, &
                      restart_file_type, restore_state, free_restart_type, save_restart
use fms_mod,    only: write_data, set_domain
use horiz_interp_mod, only : horiz_interp_type
use horiz_interp_spherical_mod, only : horiz_interp_spherical, horiz_interp_spherical_del, &
                                       horiz_interp_spherical_new
use mpp_domains_mod, only : mpp_update_domains

! UCLDASV2 modules
use ucldasv2_fields_metadata_mod, only : ucldasv2_field_metadata
use ucldasv2_geom_mod, only : ucldasv2_geom

implicit none
private


! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
!> Holds all data and metadata related to a single field variable.
!!
!! Instances of these types are to be held by ucldasv2_fields.
!! The members ucldasv2_field::mask can remain \c null, in which it is assumed that
!! no mask is used.
type, public :: ucldasv2_field

  !> The internally used name of the field.
  character(len=:),     allocatable :: name

  !> The number of vertical levels.
  integer                           :: nz

  !> The actual field data.
  real(kind=kind_real), allocatable :: val(:,:,:)

  !> Pointer to the relevant mask in ucldasv2_geom_mod::ucldasv2_geom
  !!
  !! If \c null, it is assumed that no mask is present
  real(kind=kind_real),     pointer :: mask(:,:) => null()!!

  !> Pointer to the relevant longitudes in ucldasv2_geom_mod::ucldasv2_geom
  !!
  !! \note This should never remain \c null() after initialization of the class.
  real(kind=kind_real),     pointer :: lon(:,:) => null()

  !> Pointer to the relevant latitudes in ucldasv2_geom_mod::ucldasv2_geom
  !!
  !! \note This should never remain \c null() after initialization of the class.
  real(kind=kind_real),     pointer :: lat(:,:) => null()

  !> Parameters for the field as determined by the configuration yaml.
  !!
  !! see ucldasv2_fields_metadata_mod::ucldasv2_field_metadata
  type(ucldasv2_field_metadata)         :: metadata

contains

  !>\copybrief ucldasv2_field_copy \see ucldasv2_field_copy
  procedure :: copy            => ucldasv2_field_copy

  !>\copybrief ucldasv2_field_delete \see ucldasv2_field_delete
  procedure :: delete          => ucldasv2_field_delete

  !>\copybrief ucldasv2_field_check_congruent \see ucldasv2_field_check_congruent
  procedure :: check_congruent => ucldasv2_field_check_congruent

  !>\copybrief ucldasv2_field_update_halo \see ucldasv2_field_update_halo
  procedure :: update_halo     => ucldasv2_field_update_halo

  !>\copybrief ucldasv2_field_stencil_interp \see ucldasv2_field_stencil_interp
  procedure :: stencil_interp  => ucldasv2_field_stencil_interp

end type ucldasv2_field


! ------------------------------------------------------------------------------
!> A collection of ucldasv2_field types representing a collective state or increment.
!!
!! The base class for ucldasv2_increment_mod::ucldasv2_increment and ucldasv2_state_mod::ucldasv2_state
type, public :: ucldasv2_fields

  !> Pointer to the relevant ucldasv2_geom_mod::ucldasv2_geom
  !!
  !! \note This should never remain \c null() after initialization of the class.
  type(ucldasv2_geom),  pointer :: geom => null()

  !> The ucldasv2_field instances that make up the fields
  type(ucldasv2_field), pointer :: fields(:) => null()

contains
  !> \name constructors / destructors
  !! \{

  !> \copybrief ucldasv2_fields_create \see ucldasv2_fields_create
  procedure :: create => ucldasv2_fields_create

  !> \copybrief ucldasv2_fields_copy \see ucldasv2_fields_copy
  procedure :: copy   => ucldasv2_fields_copy

  !> \copybrief ucldasv2_fields_delete \see ucldasv2_fields_delete
  procedure :: delete => ucldasv2_fields_delete

  !> \}

  !> \name field getters/checkers
  !! \{

  !> \copybrief ucldasv2_fields_get \see ucldasv2_fields_get
  procedure :: get    => ucldasv2_fields_get

  !> \copybrief ucldasv2_fields_has \see ucldasv2_fields_has
  procedure :: has    => ucldasv2_fields_has

  !> \copybrief ucldasv2_fields_check_congruent \see ucldasv2_fields_check_congruent
  procedure :: check_congruent => ucldasv2_fields_check_congruent

  !> \copybrief ucldasv2_fields_check_subset \see ucldasv2_fields_check_subset
  procedure :: check_subset    => ucldasv2_fields_check_subset

  !> \}

  !> \name math operators
  !! \{

  !> \copybrief ucldasv2_fields_add \see ucldasv2_fields_add
  procedure :: add      => ucldasv2_fields_add

  !> \copybrief ucldasv2_fields_axpy \see ucldasv2_fields_axpy
  procedure :: axpy     => ucldasv2_fields_axpy

  !> \copybrief ucldasv2_fields_dotprod \see ucldasv2_fields_dotprod
  procedure :: dot_prod => ucldasv2_fields_dotprod

  !> \copybrief ucldasv2_fields_gpnorm \see ucldasv2_fields_gpnorm
  procedure :: gpnorm   => ucldasv2_fields_gpnorm

 !> \copybrief ucldasv2_fields_mul \see ucldasv2_fields_mul
  procedure :: mul      => ucldasv2_fields_mul

  !> \copybrief ucldasv2_fields_sub \see ucldasv2_fields_sub
  procedure :: sub      => ucldasv2_fields_sub

  !> \copybrief ucldasv2_fields_ones \see ucldasv2_fields_ones
  procedure :: ones     => ucldasv2_fields_ones

  !> \copybrief ucldasv2_fields_zeros \see ucldasv2_fields_zeros
  procedure :: zeros    => ucldasv2_fields_zeros

  !> \}

  !> \name I/O
  !! \{

  !> \copybrief ucldasv2_fields_read \see ucldasv2_fields_read
  procedure :: read      => ucldasv2_fields_read

  !> \copybrief ucldasv2_fields_write_file \see ucldasv2_fields_write_file
  procedure :: write_file=> ucldasv2_fields_write_file

  !> \copybrief ucldasv2_fields_write_rst \see ucldasv2_fields_write_rst
  procedure :: write_rst => ucldasv2_fields_write_rst

  !> \}

  !> \name misc
  !! \{

  !> \copybrief ucldasv2_fields_update_halos \see ucldasv2_fields_update_halos
  procedure :: update_halos => ucldasv2_fields_update_halos

  !> \copybrief ucldasv2_fields_colocate \see ucldasv2_fields_colocate
  procedure :: colocate  => ucldasv2_fields_colocate
  !> \}

  !> \name serialization
  !! \{

  !> \copybrief ucldasv2_fields_serial_size \see ucldasv2_fields_serial_size
  procedure :: serial_size => ucldasv2_fields_serial_size

  !> \copybrief ucldasv2_fields_serialize \see ucldasv2_fields_serialize
  procedure :: serialize   => ucldasv2_fields_serialize

  !> \copybrief ucldasv2_fields_deserialize \see ucldasv2_fields_deserialize
  procedure :: deserialize => ucldasv2_fields_deserialize

  !> \}

end type ucldasv2_fields


! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
contains


! ------------------------------------------------------------------------------
! ucldasv2_field subroutines
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
!> Copy a field from \p rhs to \p self.
!!
!! If the fields are not congruent, this subroutine will throw an error.
!! \p self must be allocated first.
!! \relates ucldasv2_fields_mod::ucldasv2_field
subroutine ucldasv2_field_copy(self, rhs)
  class(ucldasv2_field), intent(inout) :: self !< The field to copy \b to
  type(ucldasv2_field),  intent(in)    :: rhs !< The field to copy \b from

  call self%check_congruent(rhs)

  ! the only variable that should be different is %val
  self%val = rhs%val

  ! NOTE: the pointers (mask, lat, lon) will be different, but should NOT
  ! be changed to point to rhs pointers. Bad things happen
end subroutine ucldasv2_field_copy


! ------------------------------------------------------------------------------
!> Update the data in the halo region of the field.
!!
!! \relates ucldasv2_fields_mod::ucldasv2_field
!! \todo have field keep a pointer to its relevant sections of ucldasv2_geom?
subroutine ucldasv2_field_update_halo(self, geom)
  class(ucldasv2_field),     intent(inout) :: self
  type(ucldasv2_geom), pointer, intent(in) :: geom !< ucldasv2_geom from ucldasv2_fields

  call mpp_update_domains(self%val, geom%Domain%mpp_domain)
end subroutine ucldasv2_field_update_halo


! ------------------------------------------------------------------------------
!> Perform spatial interpolation between two grids.
!!
!! Interpolation used is inverse distance weidghted, taking into
!! consideration the mask.
!! \param[in] geom: The geometry to interpolate to
!! \param[in] interp2d: interpolation object created by calling
!!     \c horiz_interp_spherical_new() in FMS
!!
!! \relates ucldasv2_fields_mod::ucldasv2_field
subroutine ucldasv2_field_stencil_interp(self, geom, interp2d)
  class(ucldasv2_field),     intent(inout) :: self
  type(ucldasv2_geom), pointer, intent(in) :: geom
  type(horiz_interp_type),  intent(in) :: interp2d

  integer :: k
  real(kind=kind_real), allocatable :: val(:,:,:)

  allocate(val, mold=self%val)
  val = self%val
  do k = 1, self%nz
     call horiz_interp_spherical(interp2d, &
          & val(geom%isd:geom%ied, geom%jsd:geom%jed,k), &
          & self%val(geom%isc:geom%iec, geom%jsc:geom%jec,k))
  end do
  call self%update_halo(geom)
end subroutine ucldasv2_field_stencil_interp


! ------------------------------------------------------------------------------
!> Make sure the two fields are the same in terms of name, size, shape.
!!
!! \throws abor1_ftn Halts program if fields are not congruent
!! \relates ucldasv2_fields_mod::ucldasv2_field
subroutine ucldasv2_field_check_congruent(self, rhs)
  class(ucldasv2_field), intent(in) :: self
  type(ucldasv2_field),  intent(in) :: rhs !< other field to check for congruency
  integer :: i

  if ( self%nz /= rhs%nz ) call abor1_ftn("ucldasv2_field:  self%nz /= rhs%nz")
  if ( self%name /= rhs%name ) call abor1_ftn("ucldasv2_field:  self%name /= rhs%name")
  if ( size(shape(self%val)) /= size(shape(rhs%val)) ) &
    call abor1_ftn("ucldasv2_field: shape of self%val /= rhs%val")
  do i =1, size(shape(self%val))
    if (size(self%val, dim=i) /= size(rhs%val, dim=i)) &
      call abor1_ftn("ucldasv2_field: shape of self%val /= rhs%val")
  end do
end subroutine ucldasv2_field_check_congruent


! ------------------------------------------------------------------------------
!> Delete the ucldasv2_field object.
!!
!! \relates ucldasv2_fields_mod::ucldasv2_field
subroutine ucldasv2_field_delete(self)
  class(ucldasv2_field), intent(inout) :: self

  deallocate(self%val)
end subroutine


! ------------------------------------------------------------------------------
! ucldasv2_fields subroutines
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
!> For a given list of field names, initialize the properties of those fields
!!
!! \param[in] vars: List of variables to initialize. They must be present in the
!!   configuration file used to create ucldasv2_fields_metadata_mod::ucldasv2_fields_metadata
!!
!! \throws abor1_ftn aborts if illegal grid or levels specified
!! \relates ucldasv2_fields_mod::ucldasv2_fields
subroutine ucldasv2_fields_init_vars(self, vars)
  class(ucldasv2_fields),         intent(inout) :: self
  character(len=:), allocatable, intent(in) :: vars(:)

  integer :: i, nz

  allocate(self%fields(size(vars)))
  do i=1,size(vars)
    self%fields(i)%name = trim(vars(i))

    ! get the field metadata parameters that are read in from a config file
    self%fields(i)%metadata = self%geom%fields_metadata%get(self%fields(i)%name)

    ! Set grid location and masks
    select case(self%fields(i)%metadata%grid)
    case ('h')
      self%fields(i)%lon => self%geom%lon
      self%fields(i)%lat => self%geom%lat
      if (self%fields(i)%metadata%masked) &
          self%fields(i)%mask => self%geom%mask2d
    case default
      call abor1_ftn('ucldasv2_fields::create(): Illegal grid '// &
                     self%fields(i)%metadata%grid // &
                     ' given for ' // self%fields(i)%name)
    end select

    ! determine number of levels
    if (self%fields(i)%name == self%fields(i)%metadata%getval_name_surface) then
      ! if this field is a surface getval, override the number of levels with 1
      nz = 1
    else
      select case(self%fields(i)%metadata%levels)
      case ('full_lnd')
        nz = self%geom%nsoil
      case ('full_snw')
        nz = self%geom%nsnow
      case ('1') ! TODO, generalize to work with any number?
        nz = 1
      case default
        call abor1_ftn('ucldasv2_fields::create(): Illegal levels '//self%fields(i)%metadata%levels// &
                       ' given for ' // self%fields(i)%name)
      end select
    endif

    ! allocate space
    self%fields(i)%nz = nz
    allocate(self%fields(i)%val(&
      self%geom%isd:self%geom%ied, &
      self%geom%jsd:self%geom%jed, &
      nz ))

  end do
end subroutine


! ------------------------------------------------------------------------------
!> Create a new set of fields, allocate space for them, and initialize to zero
!!
!! \see ucldasv2_fields_init_vars
!! \relates ucldasv2_fields_mod::ucldasv2_fields
subroutine ucldasv2_fields_create(self, geom, vars)
  class(ucldasv2_fields),        intent(inout) :: self
  type(ucldasv2_geom),  pointer, intent(inout) :: geom !< geometry to associate with the fields
  type(oops_variables),      intent(inout) :: vars !< list of field names to create

  character(len=:), allocatable :: vars_str(:)
  integer :: i

  ! make sure current object has not already been allocated
  if (associated(self%fields)) &
    call abor1_ftn("ucldasv2_fields::create(): object already allocated")

  ! associate geometry
  self%geom => geom

  ! initialize the variable parameters
  allocate(character(len=1024) :: vars_str(vars%nvars()))
  do i=1,vars%nvars()
    vars_str(i) = trim(vars%variable(i))
  end do
  call ucldasv2_fields_init_vars(self, vars_str)

  ! set everything to zero
  call self%zeros()
end subroutine ucldasv2_fields_create


! ------------------------------------------------------------------------------
!> delete all the fields
!!
!! \relates ucldasv2_fields_mod::ucldasv2_fields
subroutine ucldasv2_fields_delete(self)
  class(ucldasv2_fields), intent(inout) :: self
  integer :: i

  ! clear the fields and nullify pointers
  nullify(self%geom)
  do i = 1, size(self%fields)
    call self%fields(i)%delete()
  end do
  deallocate(self%fields)
  nullify(self%fields)

end subroutine


! ------------------------------------------------------------------------------
!> Copy the contents of \p rhs to \p self.
!!
!! \p self will be initialized with the variable names in \p rhs if
!! not already initialized.
!!
!! \see ucldasv2_fields_init_vars
!! \relates ucldasv2_fields_mod::ucldasv2_fields
subroutine ucldasv2_fields_copy(self, rhs)
  class(ucldasv2_fields), intent(inout) :: self
  class(ucldasv2_fields),  intent(in)    :: rhs !< fields to copy from

  character(len=:), allocatable :: vars_str(:)
  integer :: i
  type(ucldasv2_field), pointer :: rhs_fld

  ! initialize the variables based on the names in rhs
  if (.not. associated(self%fields)) then
    self%geom => rhs%geom
    allocate(character(len=1024) :: vars_str(size(rhs%fields)))
    do i=1, size(vars_str)
      vars_str(i) = rhs%fields(i)%name
    end do
    call ucldasv2_fields_init_vars(self, vars_str)
  end if

  ! copy values from rhs to self, only if the variable exists
  !  in self
  do i=1,size(self%fields)
    call rhs%get(self%fields(i)%name, rhs_fld)
    call self%fields(i)%copy(rhs_fld)
  end do

end subroutine


! ------------------------------------------------------------------------------
!> Get a pointer to the ucldasv2_field with the given name.
!!
!! \note use ucldasv2_fields::has() if you need to check for optional fields
!! \throws abor1_ftn If no field exists with that name, the prorgam aborts
!! \relates ucldasv2_fields_mod::ucldasv2_fields
subroutine ucldasv2_fields_get(self, name, field)
  class(ucldasv2_fields),         intent(in) :: self
  character(len=*),           intent(in) :: name !< name of field to find
  type(ucldasv2_field), pointer, intent(out) :: field  !< a pointer to the resulting field

  integer :: i

  ! find the field with the given name
  do i=1,size(self%fields)
    if (trim(name) == self%fields(i)%name) then
      field => self%fields(i)
      return
    end if
  end do

  ! oops, the field was not found
  call abor1_ftn("ucldasv2_fields::get():  cannot find field "//trim(name))
end subroutine


! ------------------------------------------------------------------------------
!> Returns whether a field with the given name exists
!!
!! \relates ucldasv2_fields_mod::ucldasv2_fields
function ucldasv2_fields_has(self, name) result(res)
  class(ucldasv2_fields), intent(in) :: self
  character(len=*),   intent(in) :: name !< name of field to find

  logical :: res
  integer :: i

  res = .false.
  do i=1,size(self%fields)
    if (trim(name) == self%fields(i)%name) then
      res = .true.
      return
    end if
  end do
end function


! ------------------------------------------------------------------------------
!> Update the halo region of all fields.
!!
!! \relates ucldasv2_fields_mod::ucldasv2_fields
subroutine ucldasv2_fields_update_halos(self)
  class(ucldasv2_fields), intent(inout) :: self
  integer :: i

  do i=1,size(self%fields)
    call self%fields(i)%update_halo(self%geom)
  end do
end subroutine ucldasv2_fields_update_halos


! ------------------------------------------------------------------------------
!> Set the value of all fields to one.
!!
!! \relates ucldasv2_fields_mod::ucldasv2_fields
subroutine ucldasv2_fields_ones(self)
  class(ucldasv2_fields), intent(inout) :: self
  integer :: i

  do i = 1, size(self%fields)
    self%fields(i)%val = 1.0_kind_real
  end do

end subroutine ucldasv2_fields_ones


! ------------------------------------------------------------------------------
!> Reset the value of all fields to zero.
!!
!! \relates ucldasv2_fields_mod::ucldasv2_fields
subroutine ucldasv2_fields_zeros(self)
  class(ucldasv2_fields), intent(inout) :: self
  integer :: i

  do i = 1, size(self%fields)
    self%fields(i)%val = 0.0_kind_real
  end do

end subroutine ucldasv2_fields_zeros


! ------------------------------------------------------------------------------
!> Add two sets of fields together
!!
!! \f$ self = self + rhs \f$
!!
!! \throws abor1_ftn aborts if two fields are not congruent
!! \relates ucldasv2_fields_mod::ucldasv2_fields
subroutine ucldasv2_fields_add(self, rhs)
  class(ucldasv2_fields), intent(inout) :: self
  class(ucldasv2_fields),     intent(in) :: rhs !< other field to add
  integer :: i

  ! make sure fields are same shape
  call self%check_congruent(rhs)

  ! add
  do i=1,size(self%fields)
    self%fields(i)%val = self%fields(i)%val + rhs%fields(i)%val
  end do
end subroutine ucldasv2_fields_add


! ------------------------------------------------------------------------------
!> subtract two sets of fields
!!
!! \f$ self = self - rhs \f$
!!
!! \throws abor1_ftn aborts if two fields are not congruent
!! \relates ucldasv2_fields_mod::ucldasv2_fields
subroutine ucldasv2_fields_sub(self, rhs)
  class(ucldasv2_fields), intent(inout) :: self
  class(ucldasv2_fields),     intent(in) :: rhs !< other field to subtract
  integer :: i

  ! make sure fields are same shape
  call self%check_congruent(rhs)

  ! subtract
  do i=1,size(self%fields)
    self%fields(i)%val = self%fields(i)%val - rhs%fields(i)%val
  end do
end subroutine ucldasv2_fields_sub


! ------------------------------------------------------------------------------
!> Multiply a set of fields by a constant.
!!
!! \f$ self = zz * self \f$
!! \relates ucldasv2_fields_mod::ucldasv2_fields
subroutine ucldasv2_fields_mul(self, zz)
  class(ucldasv2_fields), intent(inout) :: self
  real(kind=kind_real),  intent(in) :: zz !< the constant by which to multipy the field
  integer :: i

  do i=1,size(self%fields)
    self%fields(i)%val = zz * self%fields(i)%val
  end do
end subroutine ucldasv2_fields_mul


! ------------------------------------------------------------------------------
!> Add two fields (multiplying the rhs first)
!!
!! \f$self = self + zz * rhs\f$
!!
!! \throws abor1_ftn aborts if \p is not a subset of \rhs
!! \relates ucldasv2_fields_mod::ucldasv2_fields
subroutine ucldasv2_fields_axpy(self, zz, rhs)
  class(ucldasv2_fields), intent(inout) :: self
  real(kind=kind_real),  intent(in) :: zz !< constant by which to multiply other rhs
  class(ucldasv2_fields),    intent(in) :: rhs !< other field to add

  type(ucldasv2_field), pointer :: f_rhs, f_lhs
  integer :: i

  ! make sure fields are correct shape
  call self%check_subset(rhs)

  do i=1,size(self%fields)
    f_lhs => self%fields(i)
    if (.not. rhs%has(f_lhs%name)) cycle
    call rhs%get(f_lhs%name, f_rhs)
    f_lhs%val = f_lhs%val + zz *f_rhs%val
  end do
end subroutine ucldasv2_fields_axpy

! ------------------------------------------------------------------------------
!> Calculate the global dot product of two sets of fields.
!!
!! \throws abor1_ftn aborts if two fields are not congruent
!! \relates ucldasv2_fields_mod::ucldasv2_fields
subroutine ucldasv2_fields_dotprod(self, rhs, zprod)
  class(ucldasv2_fields),     intent(in) :: self
  class(ucldasv2_fields),      intent(in) :: rhs !< field 2 of dot product
  real(kind=kind_real),  intent(out) :: zprod !< The resulting dot product

  real(kind=kind_real) :: local_zprod
  integer :: ii, jj, kk, n
  type(ucldasv2_field), pointer :: field1, field2

  ! make sure fields are same shape
  call self%check_congruent(rhs)

  ! loop over (almost) all fields
  local_zprod = 0.0_kind_real
  do n=1,size(self%fields)
    field1 => self%fields(n)
    field2 => rhs%fields(n)

    ! add the given field to the dot product (only using the compute domain)
    do ii = self%geom%isc, self%geom%iec
      do jj = self%geom%jsc, self%geom%jec
        ! masking
        if (associated(field1%mask)) then
          if (field1%mask(ii,jj) < 1) cycle
        endif

        ! add to dot product
        do kk=1,field1%nz
          local_zprod = local_zprod + field1%val(ii,jj,kk) * field2%val(ii,jj,kk)
        end do
      end do
    end do
  end do

  ! Get global dot product
  call self%geom%f_comm%allreduce(local_zprod, zprod, fckit_mpi_sum())
end subroutine ucldasv2_fields_dotprod

! ------------------------------------------------------------------------------
!> read a set of fields from a file
!!
!! \param[in] f_conf : Configuration with the following parameters
!!    - "read_from_file" :
!!      - 0 = Invent the state
!!      - 1 = read state
!!      - 2 = (nothing??)
!!      - 3 = read increment
!!    - "date" : (required if read_from_file == 0)
!!    - "basename" : The common part of the path prepended to the following
!!       \c *_filename parameters
!!    - "lnd_filename" : land filename
!! \param[inout] vdate : If fields are being invented (read_from_file == 0),
!!    the \p vdate is used as the valid date of the fields. If the fields are
!!    being read in as a state (read_from_file == 1), \p vdate is set the the
!!    date from the files
!! \relates ucldasv2_fields_mod::ucldasv2_fields
subroutine ucldasv2_fields_read(self, f_conf, vdate)
  class(ucldasv2_fields),        intent(inout) :: self
  type(fckit_configuration), intent(in)    :: f_conf
  type(datetime),            intent(inout) :: vdate

  integer, parameter :: max_string_length=800
  character(len=max_string_length) :: lnd_filename, filename
  character(len=:), allocatable :: basename, incr_filename
  integer :: iread = 0
  integer :: ii
  type(restart_file_type), target :: land_restart
  type(restart_file_type), pointer :: restart
  integer :: idr
  integer :: isd, ied, jsd, jed
  integer :: isc, iec, jsc, jec
  integer :: i, j, nz, n
  character(len=:), allocatable :: str
  type(ucldasv2_field), pointer :: field, field2, snowd, mld, layer_depth

  if ( f_conf%has("read_from_file") ) &
      call f_conf%get_or_die("read_from_file", iread)

  call self%get("snowd", snowd)

  ! Get Indices for data domain and allocate common layer depth array
  isd = self%geom%isd ; ied = self%geom%ied
  jsd = self%geom%jsd ; jed = self%geom%jed

  nz = snowd%nz

  ! iread = 0: Invent state
  if (iread==0) then
     call self%zeros()
     call f_conf%get_or_die("date", str)
     call datetime_set(str, vdate)
  end if

  ! iread = 0: Invent state
  if (iread==0) then
     call self%zeros()
     call f_conf%get_or_die("date", str)
     call datetime_set(str, vdate)
  end if

  ! TODO redo this to be generic

  ! iread = 1 (state) or 3 (increment): Read restart file
  if ((iread==1).or.(iread==3)) then

    ! filename for land
    call f_conf%get_or_die("basename", str)
    basename = str
    call f_conf%get_or_die("lnd_filename", str)
    lnd_filename = trim(basename) // trim(str)

    call fms_io_init()

    ! built-in variables
    do i=1,size(self%fields)
      if(self%fields(i)%metadata%io_name /= "") then
      ! which file are we reading from?
        select case(self%fields(i)%metadata%io_file)
        case ('lnd')
          filename = lnd_filename
          restart => land_restart
        case default
          call abor1_ftn('read_file(): illegal io_file: '//self%fields(i)%metadata%io_file)
        end select

      ! setup to read
        if (self%fields(i)%nz == 1) then
          idr = register_restart_field(restart, filename, self%fields(i)%metadata%io_name, &
              self%fields(i)%val(:,:,1), domain=self%geom%Domain%mpp_domain)
        else
          idr = register_restart_field(restart, filename, self%fields(i)%metadata%io_name, &
              self%fields(i)%val(:,:,:), domain=self%geom%Domain%mpp_domain)
        end if
      end if
    end do

    call restore_state(land_restart, directory='')
    call free_restart_type(land_restart)

    call fms_io_exit()

    ! Indices for compute domain
    isc = self%geom%isc ; iec = self%geom%iec
    jsc = self%geom%jsc ; jec = self%geom%jec

    ! Update halo
    do n=1,size(self%fields)
      field => self%fields(n)
      call mpp_update_domains(field%val, self%geom%Domain%mpp_domain)
    end do

    ! Set vdate if reading state
    if (iread==1) then
      call f_conf%get_or_die("date", str)
      call datetime_set(str, vdate)
    end if

    return
  end if

end subroutine ucldasv2_fields_read


! ------------------------------------------------------------------------------
!> calculate global statistics for each field (min, max, average)
!!
!! \param[in] nf: The number of fields, should be equal to the size of
!!     ucldasv2_fields::fields
!! \param[out] pstat: a 2D array with shape (i,j). For each field index
!!     i is set as 0 = min, 1 = max, 2 = average, for j number of fields.
!! \relates ucldasv2_fields_mod::ucldasv2_fields
subroutine ucldasv2_fields_gpnorm(self, nf, pstat)
  class(ucldasv2_fields),      intent(in) :: self
  integer,                 intent(in) :: nf
  real(kind=kind_real),   intent(out) :: pstat(3, nf)

  logical :: mask(self%geom%isc:self%geom%iec, self%geom%jsc:self%geom%jec)
  real(kind=kind_real) :: lnd_count, local_lnd_count, tmp(3)
  integer :: jj, isc, iec, jsc, jec
  type(ucldasv2_field), pointer :: field

  ! Indices for compute domain
  isc = self%geom%isc ; iec = self%geom%iec
  jsc = self%geom%jsc ; jec = self%geom%jec

  ! calculate global min, max, mean for each field
  do jj=1, size(self%fields)
    call self%get(self%fields(jj)%name, field)

    ! get the mask and the total number of grid cells
    if (.not. associated(field%mask)) then
       mask = .true.
     else
       mask = field%mask(isc:iec, jsc:jec) > 0.0
     end if
    local_lnd_count = count(mask)
    call self%geom%f_comm%allreduce(local_lnd_count, lnd_count, fckit_mpi_sum())

    ! calculate global min/max/mean
    call fldinfo(field%val(isc:iec,jsc:jec,:), mask, tmp)
    call self%geom%f_comm%allreduce(tmp(1), pstat(1,jj), fckit_mpi_min())
    call self%geom%f_comm%allreduce(tmp(2), pstat(2,jj), fckit_mpi_max())
    call self%geom%f_comm%allreduce(tmp(3), pstat(3,jj), fckit_mpi_sum())
    pstat(3,jj) = pstat(3,jj)/lnd_count
  end do
end subroutine ucldasv2_fields_gpnorm


! ------------------------------------------------------------------------------
!> Make sure two sets of fields are the same shape (same variables, same resolution)
!!
!! \throws abor1_ftn aborts if two fields are not congruent.
!!
!! \relates ucldasv2_fields_mod::ucldasv2_fields
subroutine ucldasv2_fields_check_congruent(self, rhs)
  class(ucldasv2_fields), intent(in) :: self
  class(ucldasv2_fields), intent(in) :: rhs !< other fields to check for congruency

  integer :: i, j

  ! number of fields should be the same
  if (size(self%fields) /= size(rhs%fields)) &
    call abor1_ftn("ucldasv2_fields: contains different number of fields")

  ! each field should match (name, size, shape)
  do i=1,size(self%fields)
    if (self%fields(i)%name /= rhs%fields(i)%name) &
      call abor1_ftn("ucldasv2_fields: field have different names")
    do j = 1, size(shape(self%fields(i)%val))
      if (size(self%fields(i)%val, dim=j) /= size(rhs%fields(i)%val, dim=j) ) then
        call abor1_ftn("ucldasv2_fields: field '"//self%fields(i)%name//"' has different dimensions")
      end if
    end do
  end do
end subroutine ucldasv2_fields_check_congruent


! ------------------------------------------------------------------------------
!> make sure two sets of fields are the same shape for fields they have in common
!!
!! \throws abor1_ftn aborts if \p self is not a subset of \p rhs
!! \relates ucldasv2_fields_mod::ucldasv2_fields
subroutine ucldasv2_fields_check_subset(self, rhs)
  class(ucldasv2_fields), intent(in) :: self
  class(ucldasv2_fields), intent(in) :: rhs !< other field that \p self should be subset of

  type(ucldasv2_field), pointer :: fld
  integer :: i, j

  ! each field should match (name, size, shape)
  do i=1,size(self%fields)
    if (.not. rhs%has(self%fields(i)%name)) &
      call abor1_ftn("ucldasv2_fields: self is not a subset of rhs")
    call rhs%get(self%fields(i)%name, fld)
    do j = 1, size(shape(fld%val))
      if (size(self%fields(i)%val, dim=j) /= size(fld%val, dim=j) ) then
        call abor1_ftn("ucldasv2_fields: field '"//self%fields(i)%name//"' has different dimensions")
      end if
    end do
  end do
end subroutine ucldasv2_fields_check_subset


! ------------------------------------------------------------------------------
!> Save ucldasv2 fields to file using fms write_data
!!
!! \param[in] filename : The name of the file to save to
!!
!! \relates ucldasv2_fields_mod::ucldasv2_fields
subroutine ucldasv2_fields_write_file(self, filename)
  class(ucldasv2_fields),  intent(in) :: self
  character(len=*),   intent(in) :: filename

  integer :: ii

  call fms_io_init()
  call set_domain( self%geom%Domain%mpp_domain )

  ! write out all fields
  do ii = 1, size(self%fields)
    call write_data( filename, self%fields(ii)%name, self%fields(ii)%val(:,:,:), self%geom%Domain%mpp_domain)
  end do

  call fms_io_exit()
end subroutine ucldasv2_fields_write_file


! ------------------------------------------------------------------------------
!> Save ucldasv2 fields in a restart format
!!
!! TODO this can be generalized even more
!! \relates ucldasv2_fields_mod::ucldasv2_fields
subroutine ucldasv2_fields_write_rst(self, f_conf, vdate)
  class(ucldasv2_fields),    intent(inout) :: self      !< Fields
  type(fckit_configuration), intent(in)    :: f_conf   !< Configuration
  type(datetime),            intent(inout) :: vdate    !< DateTime

  integer, parameter :: max_string_length=800
  character(len=max_string_length) :: lnd_filename, filename
  type(restart_file_type), target :: land_restart
  type(restart_file_type), pointer :: restart
  integer :: idr, i
  type(ucldasv2_field), pointer :: field

  call fms_io_init()

  ! filenames
  lnd_filename = ucldasv2_genfilename(f_conf,max_string_length,vdate,"lnd")

  ! built in variables
  do i=1,size(self%fields)
    field => self%fields(i)
    if (len_trim(field%metadata%io_file) /= 0) then
      ! which file are we writing to
      select case(field%metadata%io_file)
      case ('lnd')
        filename = lnd_filename
        restart => land_restart
      case default
        call abor1_ftn('ucldasv2_write_restart(): illegal io_file: '//field%metadata%io_file)
      end select

      ! write
      if (field%nz == 1) then
        idr = register_restart_field( restart, filename, field%metadata%io_name, &
          field%val(:,:,1), domain=self%geom%Domain%mpp_domain)
      else
        idr = register_restart_field( restart, filename, field%metadata%io_name, &
        field%val(:,:,:), domain=self%geom%Domain%mpp_domain)
      end if
    end if
  end do

  ! write out and cleanup
  call save_restart(land_restart, directory='')
  call free_restart_type(land_restart)
  call fms_io_exit()

end subroutine ucldasv2_fields_write_rst

! ------------------------------------------------------------------------------
!> Colocate by interpolating from one c-grid location to another.
!!
!! \warning only works on the "h" grid currently (not the "u" or "v" grid)
!! \relates ucldasv2_fields_mod::ucldasv2_fields
subroutine ucldasv2_fields_colocate(self, cgridlocout)
  class(ucldasv2_fields),    intent(inout) :: self !< self
  character(len=1),         intent(in) :: cgridlocout !< colocate to cgridloc (u, v or h)

  integer :: i, k
  real(kind=kind_real), allocatable :: val(:,:,:)
  real(kind=kind_real), pointer :: lon_out(:,:) => null()
  real(kind=kind_real), pointer :: lat_out(:,:) => null()
  type(ucldasv2_geom),  pointer :: g => null()
  type(horiz_interp_type) :: interp2d

  ! Associate lon_out and lat_out according to cgridlocout
  select case(cgridlocout)
  case ('h')
    lon_out => self%geom%lon
    lat_out => self%geom%lat
  case default
    call abor1_ftn('ucldasv2_fields::colocate(): unknown c-grid location '// cgridlocout)
  end select

  ! Apply interpolation to all fields, when necessary
  do i=1,size(self%fields)

    ! Check if already colocated
    !if (self%fields(i)%metadata%grid == cgridlocout) cycle

    ! Initialize fms spherical idw interpolation
     g => self%geom
     call horiz_interp_spherical_new(interp2d, &
       & real(deg2rad*self%fields(i)%lon(g%isd:g%ied,g%jsd:g%jed), 8), &
       & real(deg2rad*self%fields(i)%lat(g%isd:g%ied,g%jsd:g%jed), 8), &
       & real(deg2rad*lon_out(g%isc:g%iec,g%jsc:g%jec), 8), &
       & real(deg2rad*lat_out(g%isc:g%iec,g%jsc:g%jec), 8))

    ! Make a temporary copy of field
    if (allocated(val)) deallocate(val)
    allocate(val, mold=self%fields(i)%val)
    val = self%fields(i)%val

    ! Interpolate all levels
    do k = 1, self%fields(i)%nz
      call self%fields(i)%stencil_interp(self%geom, interp2d)
    end do

    ! Update c-grid location
    self%fields(i)%metadata%grid = cgridlocout
    select case(cgridlocout)
    case ('h')
      self%fields(i)%lon => self%geom%lon
      self%fields(i)%lat => self%geom%lat
    end select

 end do
 call horiz_interp_spherical_del(interp2d)

end subroutine ucldasv2_fields_colocate


! ------------------------------------------------------------------------------
!> Number of elements to return in the serialized array
!!
!! \see ucldasv2_fields_serialize
!! \relates ucldasv2_fields_mod::ucldasv2_fields
subroutine ucldasv2_fields_serial_size(self, geom, vec_size)
  class(ucldasv2_fields),    intent(in)  :: self
  type(ucldasv2_geom),       intent(in)  :: geom !< todo remove, not needed?
  integer,               intent(out) :: vec_size !< resulting size of vector

  integer :: i

  ! Loop over fields
  vec_size = 0
  do i=1,size(self%fields)
    vec_size = vec_size + size(self%fields(i)%val)
  end do

end subroutine ucldasv2_fields_serial_size


! ------------------------------------------------------------------------------
!> Return the fields as a serialized array
!!
!! \see ucldasv2_fields_serial_size
!! \relates ucldasv2_fields_mod::ucldasv2_fields
subroutine ucldasv2_fields_serialize(self, geom, vec_size, vec)
  class(ucldasv2_fields),    intent(in)  :: self
  type(ucldasv2_geom),       intent(in)  :: geom  !< todo remove this, not needed?
  integer,               intent(in)  :: vec_size !< size of vector to return
  real(kind=kind_real),  intent(out) :: vec(vec_size) !< fields as a serialized vector

  integer :: index, i, nn

  ! Loop over fields, levels and horizontal points
  index = 1
  do i=1,size(self%fields)
    nn = size(self%fields(i)%val)
    vec(index:index+nn-1) = reshape(self%fields(i)%val, (/ nn /) )
    index = index + nn
  end do

end subroutine ucldasv2_fields_serialize

! ------------------------------------------------------------------------------
!> Deserialize, creating fields from a single serialized array
!!
!! \see ucldasv2_fields_serialize
!! \relates ucldasv2_fields_mod::ucldasv2_fields
subroutine ucldasv2_fields_deserialize(self, geom, vec_size, vec, index)
  class(ucldasv2_fields), intent(inout) :: self
  type(ucldasv2_geom),       intent(in)    :: geom !< todo remove this, not needed?
  integer,               intent(in)    :: vec_size !< size of \p vec
  real(kind=kind_real),  intent(in)    :: vec(vec_size) !< vector to deserialize
  integer,               intent(inout) :: index !< index in \p vec at which to start deserializing

  integer :: i, nn

  ! Loop over fields, levels and horizontal points
  do i=1,size(self%fields)
    nn = size(self%fields(i)%val)
    self%fields(i)%val = reshape(vec(index+1:index+1+nn), shape(self%fields(i)%val))
    index = index + nn
  end do

end subroutine ucldasv2_fields_deserialize


! ------------------------------------------------------------------------------
! Internal module functions/subroutines
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> Calculate min/max/mean statistics for a given field, using a mask.
!!
!! \param[in] fld : the field to calculate the statistics on
!! \param[in] mask : statistics are only calculated where \p mask is \c .true.
!! \param[out] info : [0] = min, [1] = max, [2] = average
subroutine fldinfo(fld, mask, info)
  real(kind=kind_real),  intent(in) :: fld(:,:,:)
  logical,               intent(in) :: mask(:,:)
  real(kind=kind_real), intent(out) :: info(3)

  integer :: z
  real(kind=kind_real) :: tmp(3,size(fld, dim=3))

  ! calculate the min/max/sum separately for each masked level
  do z = 1, size(tmp, dim=2)
     tmp(1,z) = minval(fld(:,:,z), mask=mask)
     tmp(2,z) = maxval(fld(:,:,z), mask=mask)
     tmp(3,z) = sum(   fld(:,:,z), mask=mask) / size(fld, dim=3)
  end do

  ! then combine the min/max/sum over all levels
  info(1) = minval(tmp(1,:))
  info(2) = maxval(tmp(2,:))
  info(3) = sum(   tmp(3,:))
end subroutine fldinfo

! ------------------------------------------------------------------------------
!> Generate filename (based on oops/qg)
!!
!! The configuration \p f_conf is expected to provide the following
!! - "datadir" : the directory the filenames should be prefixed with
!! - "exp" : experiment name
!! - "type" : one of "fc", "an", "incr", "ens"
!! - "member" : required only if "type == ens"
function ucldasv2_genfilename (f_conf,length,vdate,domain_type)
  type(fckit_configuration),  intent(in) :: f_conf
  integer,                    intent(in) :: length
  type(datetime),             intent(in) :: vdate
  character(len=3), optional, intent(in) :: domain_type

  character(len=length)                  :: ucldasv2_genfilename
  character(len=length) :: fdbdir, expver, typ, validitydate, referencedate, sstep, &
       & prefix, mmb
  type(datetime) :: rdate
  type(duration) :: step
  integer lenfn
  character(len=:), allocatable :: str

  call f_conf%get_or_die("datadir", str)
  fdbdir = str
  call f_conf%get_or_die("exp", str)
  expver = str
  call f_conf%get_or_die("type", str)
  typ = str

  if (present(domain_type)) then
     expver = trim(domain_type)//"."//expver
  else
     expver = "lnd.snw."//expver
  end if
  if (typ=="ens") then
     call f_conf%get_or_die("member", str)
     mmb = str
     lenfn = LEN_TRIM(fdbdir) + 1 + LEN_TRIM(expver) + 1 + LEN_TRIM(typ) + 1 + LEN_TRIM(mmb)
     prefix = TRIM(fdbdir) // "/" // TRIM(expver) // "." // TRIM(typ) // "." // TRIM(mmb)
  else
     lenfn = LEN_TRIM(fdbdir) + 1 + LEN_TRIM(expver) + 1 + LEN_TRIM(typ)
     prefix = TRIM(fdbdir) // "/" // TRIM(expver) // "." // TRIM(typ)
  endif

  if (typ=="fc" .or. typ=="ens") then
     call f_conf%get_or_die("date", str)
     referencedate = str
     call datetime_to_string(vdate, validitydate)
     call datetime_create(TRIM(referencedate), rdate)
     call datetime_diff(vdate, rdate, step)
     call duration_to_string(step, sstep)
     lenfn = lenfn + 1 + LEN_TRIM(referencedate) + 1 + LEN_TRIM(sstep)
     ucldasv2_genfilename = TRIM(prefix) // "." // TRIM(referencedate) // "." // TRIM(sstep)
  endif

  if (typ=="an" .or. typ=="incr") then
     call datetime_to_string(vdate, validitydate)
     lenfn = lenfn + 1 + LEN_TRIM(validitydate)
     ucldasv2_genfilename = TRIM(prefix) // "." // TRIM(validitydate)
  endif

  if (lenfn>length) &
       & call abor1_ftn("fields:genfilename: filename too long")

   if ( allocated(str) ) deallocate(str)

end function ucldasv2_genfilename


end module ucldasv2_fields_mod
