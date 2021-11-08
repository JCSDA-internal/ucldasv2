!
! (C) Copyright    2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence
! Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module ucldasv2_geom_mod

use atlas_module, only: atlas_functionspace_pointcloud, atlas_fieldset, &
    atlas_field, atlas_real, atlas_integer, atlas_geometry, atlas_indexkdtree

! jedi modules
use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module, only: fckit_mpi_comm, fckit_mpi_sum
use kinds, only: kind_real
use fms_io_mod, only : fms_io_init, fms_io_exit, &
                       register_restart_field, restart_file_type, &
                       restore_state, free_restart_type
use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain, &
                            mpp_get_global_domain, mpp_update_domains

! ucland / fms modules
use LND_domains, only : LND_domain_type

! ucldasv2 modules
use ucldasv2_ucland, only: ucldasv2_geomdomain_init

implicit none

private
public :: ucldasv2_geom, &
          geom_get_domain_indices

!> Geometry data structure
type :: ucldasv2_geom
    type(LND_domain_type), pointer :: Domain !< Land model domain
    integer :: nsoil, nsnow
    integer :: isc, iec, jsc, jec  !< indices of compute domain
    integer :: isd, ied, jsd, jed  !< indices of data domain
    integer :: isg, ieg, jsg, jeg  !< indices of global domain
    integer :: iscl, iecl, jscl, jecl  !< indices of local compute domain
    integer :: isdl, iedl, jsdl, jedl  !< indices of local data domain
    real(kind=kind_real), allocatable, dimension(:,:) :: lon, lat !< Tracer point grid
    real(kind=kind_real), allocatable, dimension(:,:) :: mask2d
    real(kind=kind_real), allocatable, dimension(:,:) :: cell_elev
    real(kind=kind_real), allocatable, dimension(:,:,:) :: h
    logical :: save_local_domain = .false. ! If true, save the local geometry for each pe.
    character(len=:), allocatable :: geom_grid_file
    type(fckit_mpi_comm) :: f_comm !< MPI communicator
    type(atlas_functionspace_pointcloud) :: afunctionspace

    contains
    procedure :: init => ucldasv2_geom_init
    procedure :: end => ucldasv2_geom_end
    procedure :: clone => ucldasv2_geom_clone
    procedure :: set_atlas_lonlat => geom_set_atlas_lonlat
    procedure :: fill_atlas_fieldset => geom_fill_atlas_fieldset
end type ucldasv2_geom

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> Setup geometry object
subroutine ucldasv2_geom_init(self, f_conf, f_comm)
  class(ucldasv2_geom),      intent(out) :: self
  type(fckit_configuration), intent(in)  :: f_conf
  type(fckit_mpi_comm),      intent(in)  :: f_comm !< MPI communicator for this geometry

  character(len=:), allocatable :: str
  logical :: full_init = .false.

  ! MPI communicator
  self%f_comm = f_comm

  ! Domain decomposition
  call ucldasv2_geomdomain_init(self%Domain, self%nsnow, f_comm)

  ! User-defined grid filename
  if ( .not. f_conf%get("geom_grid_file", self%geom_grid_file) ) &
     self%geom_grid_file = "ucldasv2_gridspec.nc" ! default if not found

  ! Allocate geometry arrays
  call ucldasv2_geom_allocate(self)

  ! Check if a full initialization is required, default to false
  if ( .not. f_conf%get("full_init", full_init) ) full_init = .false.

  ! Read the geometry from file by default,
  ! skip this step if a full init is required
  if ( .not. full_init) call ucldasv2_geom_read(self)

end subroutine ucldasv2_geom_init

! ------------------------------------------------------------------------------
!> Geometry destructor
subroutine ucldasv2_geom_end(self)
  class(ucldasv2_geom), intent(out)  :: self

  if (allocated(self%lon))           deallocate(self%lon)
  if (allocated(self%lat))           deallocate(self%lat)
  if (allocated(self%mask2d))        deallocate(self%mask2d)
  if (allocated(self%cell_elev))     deallocate(self%cell_elev)
  if (allocated(self%h))             deallocate(self%h)
  nullify(self%Domain)
  call self%afunctionspace%final()

end subroutine ucldasv2_geom_end

! --------------------------------------------------------------------------------------------------
!> Set ATLAS lonlat fieldset
subroutine geom_set_atlas_lonlat(self, afieldset)
  class(ucldasv2_geom),  intent(inout) :: self
  type(atlas_fieldset), intent(inout) :: afieldset

  real(kind_real), pointer :: real_ptr(:,:)
  type(atlas_field) :: afield

  ! Create lon/lat field
  afield = atlas_field(name="lonlat", kind=atlas_real(kind_real), shape=(/2,(self%iec-self%isc+1)*(self%jec-self%jsc+1)/))
  call afield%data(real_ptr)
  real_ptr(1,:) = reshape(self%lon(self%isc:self%iec,self%jsc:self%jec),(/(self%iec-self%isc+1)*(self%jec-self%jsc+1)/))
  real_ptr(2,:) = reshape(self%lat(self%isc:self%iec,self%jsc:self%jec),(/(self%iec-self%isc+1)*(self%jec-self%jsc+1)/))
  call afieldset%add(afield)

end subroutine geom_set_atlas_lonlat

! --------------------------------------------------------------------------------------------------
!> Fill ATLAS fieldset
subroutine geom_fill_atlas_fieldset(self, afieldset)
  class(ucldasv2_geom),  intent(inout) :: self
  type(atlas_fieldset), intent(inout) :: afieldset

  integer :: i, jz, n
  integer, pointer :: int_ptr_2(:,:)
  real(kind=kind_real), pointer :: real_ptr_1(:), real_ptr_2(:,:)
  type(atlas_field) :: afield

  ! Add elev
  afield = self%afunctionspace%create_field(name='elev', kind=atlas_real(kind_real), levels=0)
  call afield%data(real_ptr_1)
  real_ptr_1 = pack(self%cell_elev(self%isc:self%iec,self%jsc:self%jec),.true.)
  call afieldset%add(afield)
  call afield%final()

 ! Add geographical mask
  afield = self%afunctionspace%create_field(name='gmask', kind=atlas_integer(kind(0)), levels=self%nsnow)
  call afield%data(int_ptr_2)
  do jz=1,self%nsnow
    int_ptr_2(jz,:) = int(reshape(self%mask2d(self%isc:self%iec,self%jsc:self%jec), &
  & (/(self%iec-self%isc+1)*(self%jec-self%jsc+1)/)))
  end do
  call afieldset%add(afield)
  call afield%final()

end subroutine geom_fill_atlas_fieldset

! ------------------------------------------------------------------------------
!> Clone, self = other
subroutine ucldasv2_geom_clone(self, other)
  class(ucldasv2_geom), intent(inout) :: self
  class(ucldasv2_geom), intent(in) :: other

  ! Clone communicator
  self%f_comm = other%f_comm

  ! Clone fms domain and vertical levels
  self%Domain => other%Domain
  self%nsnow = other%nsnow

  !
  self%geom_grid_file = other%geom_grid_file

  ! Allocate and clone geometry
  call ucldasv2_geom_allocate(self)
  self%lon = other%lon
  self%lat = other%lat
  self%mask2d = other%mask2d
  self%cell_elev = other%cell_elev
  self%h = other%h
! call other%fields_metadata%clone(self%fields_metadata)
end subroutine ucldasv2_geom_clone

! ------------------------------------------------------------------------------
!> Allocate memory and point to ucland data structure
subroutine ucldasv2_geom_allocate(self)
  class(ucldasv2_geom), intent(inout) :: self

  integer :: nsnow
  integer :: isd, ied, jsd, jed

  ! Get domain shape (number of levels, indices of data and compute domain)
  call geom_get_domain_indices(self, "compute", self%isc, self%iec, self%jsc, self%jec)
  call geom_get_domain_indices(self, "data", isd, ied, jsd, jed)
  self%isd = isd ;  self%ied = ied ; self%jsd = jsd; self%jed = jed
  call geom_get_domain_indices(self, "global", self%isg, self%ieg, self%jsg, self%jeg)
  call geom_get_domain_indices(self, "compute", self%iscl, self%iecl, self%jscl, self%jecl, local=.true.)
  call geom_get_domain_indices(self, "data", self%isdl, self%iedl, self%jsdl, self%jedl, local=.true.)
  nsnow = self%nsnow

  ! Allocate arrays on compute domain
  allocate(self%lon(isd:ied,jsd:jed));           self%lon = 0.0_kind_real
  allocate(self%lat(isd:ied,jsd:jed));           self%lat = 0.0_kind_real
  allocate(self%mask2d(isd:ied,jsd:jed));        self%mask2d = 0.0_kind_real
  allocate(self%cell_elev(isd:ied,jsd:jed));     self%cell_elev = 0.0_kind_real
  allocate(self%h(isd:ied,jsd:jed,1:nsnow));     self%h = 0.0_kind_real

end subroutine ucldasv2_geom_allocate

! ------------------------------------------------------------------------------
!> Read geometry from file
subroutine ucldasv2_geom_read(self)
  class(ucldasv2_geom), intent(inout) :: self

  integer :: idr_geom
  type(restart_file_type) :: geom_restart

  call fms_io_init()
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'lon', &
                                   &self%lon(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'lat', &
                                   &self%lat(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'elevation', &
                                   &self%cell_elev(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'mask2d', &
                                   &self%mask2d(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'h', &
                                   &self%h(:,:,:), &
                                   domain=self%Domain%mpp_domain)
  call restore_state(geom_restart, directory='')
  call free_restart_type(geom_restart)
  call fms_io_exit()

end subroutine ucldasv2_geom_read

! ------------------------------------------------------------------------------
!> Get indices for compute or data domain
subroutine geom_get_domain_indices(self, domain_type, is, ie, js, je, local)
  class(ucldasv2_geom), intent(in) :: self
  character(len=*),       intent(in) :: domain_type
  integer,               intent(out) :: is, ie, js, je
  logical,      optional, intent(in) :: local

  integer :: isc, iec, jsc, jec
  integer :: isd, ied, jsd, jed
  integer :: isg, ieg, jsg, jeg

  call mpp_get_compute_domain(self%Domain%mpp_domain,isc,iec,jsc,jec)
  call mpp_get_data_domain(self%Domain%mpp_domain,isd,ied,jsd,jed)
  call mpp_get_global_domain(self%Domain%mpp_domain, isg, ieg, jsg, jeg)
  if (present(local)) then
     isc = isc - (isd-1) ; iec = iec - (isd-1) ; ied = ied - (isd-1) ; isd = 1
     jsc = jsc - (jsd-1) ; jec = jec - (jsd-1) ; jed = jed - (jsd-1) ; jsd = 1
  end if

  select case (trim(domain_type))
  case ("compute")
     is = isc; ie = iec; js = jsc; je = jec;
  case ("data")
     is = isd; ie = ied; js = jsd; je = jed;
  case ("global")
     is = isg; ie = ieg; js = jsg; je = jeg;
  end select

end subroutine geom_get_domain_indices

end module ucldasv2_geom_mod
