!> Describes the decomposed LND domain and has routines for communications across PEs
module LND_domains

! This file is part of UCLAND. See LICENSE.md for the license.
use LND_coms, only : PE_here, root_PE, num_PEs
use LND_coms, only : broadcast, sum_across_PEs, min_across_PEs, max_across_PEs
use LND_error_handler, only : LND_error, LND_mesg, NOTE, WARNING, FATAL, is_root_pe
use LND_file_parser, only : get_param, log_param, log_version
use LND_file_parser, only : param_file_type
use LND_string_functions, only : slasher

use mpp_domains_mod, only : mpp_define_layout, mpp_get_boundary
use mpp_domains_mod, only : LND_define_io_domain => mpp_define_io_domain
use mpp_domains_mod, only : LND_define_domain => mpp_define_domains
use mpp_domains_mod, only : domain2D, domain1D, mpp_get_data_domain
use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_global_domain
use mpp_domains_mod, only : mpp_update_domains, CYCLIC_GLOBAL_DOMAIN, FOLD_NORTH_EDGE
use fms_io_mod,      only : file_exist, parse_mask_table
use fms_affinity_mod, only : fms_affinity_init, fms_affinity_set, fms_affinity_get

implicit none ; private

public :: LND_domains_init
public :: LND_define_domain
public :: PE_here, root_PE, num_PEs, sum_across_pes
public :: domain2D

!> The LND_domain_type contains information about the domain decompositoin.
type, public :: LND_domain_type
  type(domain2D), pointer :: mpp_domain => NULL() !< The FMS domain with halos
                                !! on this processor, centered at h points.
  type(domain2D), pointer :: mpp_domain_d2 => NULL() !< A coarse FMS domain with halos
                                !! on this processor, centered at h points.
  integer :: niglobal           !< The total horizontal i-domain size.
  integer :: njglobal           !< The total horizontal j-domain size.
  integer :: nihalo             !< The i-halo size in memory.
  integer :: njhalo             !< The j-halo size in memory.
  logical :: symmetric          !< True if symmetric memory is used with
                                !! this domain.
  logical :: nonblocking_updates  !< If true, non-blocking halo updates are
                                !! allowed.  The default is .false. (for now).
  logical :: thin_halo_updates  !< If true, optional arguments may be used to
                                !! specify the width of the halos that are
                                !! updated with each call.
  integer :: layout(2)          !< This domain's processor layout.  This is
                                !! saved to enable the construction of related
                                !! new domains with different resolutions or
                                !! other properties.
  integer :: io_layout(2)       !< The IO-layout used with this domain.
  integer :: X_FLAGS            !< Flag that specifies the properties of the
                                !! domain in the i-direction in a define_domain
                                !call.
  integer :: Y_FLAGS            !< Flag that specifies the properties of the
                                !! domain in the j-direction in a define_domain
                                !call.
  logical, pointer :: maskmap(:,:) => NULL() !< A pointer to an array indicating
                                !! which logical processors are actually used for
                                !! the ocean code. The other logical processors
                                !! would be contain only land points and are not
                                !! assigned to actual processors. This need not be
                                !! assigned if all logical processors are used.
end type LND_domain_type

contains

!> LND_domains_init initalizes a LND_domain_type variable, based on the information
!! read in from a param_file_type, and optionally returns data describing various'
!! properties of the domain type.
subroutine LND_domains_init(LND_dom, param_file, symmetric, static_memory, &
                            NIHALO, NJHALO, NIGLOBAL, NJGLOBAL, NIPROC, NJPROC, &
                            min_halo, domain_name, include_name, param_suffix)
  type(LND_domain_type),           pointer       :: LND_dom      !< A pointer to the LND_domain_type
                                                                 !! being defined here.
  type(param_file_type),           intent(in)    :: param_file   !< A structure to parse for
                                                                 !! run-time parameters
  logical, optional,               intent(in)    :: symmetric    !< If present, this specifies
                                            !! whether this domain is symmetric, regardless of
                                            !! whether the macro SYMMETRIC_MEMORY_ is defined.
  logical, optional,               intent(in)    :: static_memory !< If present and true, this
                         !! domain type is set up for static memory and error checking of
                         !! various input values is performed against those in the input file.
  integer, optional,               intent(in)    :: NIHALO       !< Default halo sizes, required
                                                                 !! with static memory.
  integer, optional,               intent(in)    :: NJHALO       !< Default halo sizes, required
                                                                 !! with static memory.
  integer, optional,               intent(in)    :: NIGLOBAL     !< Total domain sizes, required
                                                                 !! with static memory.
  integer, optional,               intent(in)    :: NJGLOBAL     !< Total domain sizes, required
                                                                 !! with static memory.
  integer, optional,               intent(in)    :: NIPROC       !< Processor counts, required with
                                                                 !! static memory.
  integer, optional,               intent(in)    :: NJPROC       !< Processor counts, required with
                                                                 !! static memory.
  integer, dimension(2), optional, intent(inout) :: min_halo     !< If present, this sets the
                                        !! minimum halo size for this domain in the i- and j-
                                        !! directions, and returns the actual halo size used.
  character(len=*),      optional, intent(in)    :: domain_name  !< A name for this domain, "LND"
                                                                 !! if missing.
  character(len=*),      optional, intent(in)    :: include_name !< A name for model's include file,
                                                                 !! "LND_memory.h" if missing.
  character(len=*),      optional, intent(in)    :: param_suffix !< A suffix to apply to
                                                                 !! layout-specific parameters.

  ! Local variables
  integer, dimension(2) :: layout = (/ 1, 1 /)
  integer, dimension(2) :: io_layout = (/ 0, 0 /)
  integer, dimension(4) :: global_indices
!$ integer :: land_nthreads       ! Number of Openmp threads
!$ integer :: get_cpu_affinity, omp_get_thread_num, omp_get_num_threads
!$ logical :: land_omp_hyper_thread
  integer :: nihalo_dflt, njhalo_dflt
  integer :: pe, proc_used
  integer :: X_FLAGS, Y_FLAGS
  logical :: reentrant_x, reentrant_y, tripolar_N, is_static
  logical            :: mask_table_exists
  character(len=128) :: mask_table, inputdir
  character(len=64)  :: dom_name, inc_nm
  character(len=200) :: mesg

  integer :: xsiz, ysiz, nip_parsed, njp_parsed
  integer :: isc,iec,jsc,jec ! The bounding indices of the computational domain.
  character(len=8) :: char_xsiz, char_ysiz, char_niglobal, char_njglobal
  character(len=40) :: nihalo_nm, njhalo_nm, layout_nm, io_layout_nm, masktable_nm
  character(len=40) :: niproc_nm, njproc_nm
  integer :: xhalo_d2,yhalo_d2
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl ! This module's name.

  if (.not.associated(LND_dom)) then
    allocate(LND_dom)
    allocate(LND_dom%mpp_domain)
    allocate(LND_dom%mpp_domain_d2)
  endif

  pe = PE_here()
  proc_used = num_PEs()

  mdl = "LND_domains"

  LND_dom%symmetric = .true.
  if (present(symmetric)) then ; LND_dom%symmetric = symmetric ; endif
  if (present(min_halo)) mdl = trim(mdl)//" min_halo"

  dom_name = "LND" ; inc_nm = "LND_memory.h"
  if (present(domain_name)) dom_name = trim(domain_name)
  if (present(include_name)) inc_nm = trim(include_name)

  nihalo_nm = "NIHALO" ; njhalo_nm = "NJHALO"
  layout_nm = "LAYOUT" ; io_layout_nm = "IO_LAYOUT" ; masktable_nm = "MASKTABLE"
  niproc_nm = "NIPROC" ; njproc_nm = "NJPROC"
  if (present(param_suffix)) then ; if (len(trim(adjustl(param_suffix))) > 0) then
    nihalo_nm = "NIHALO"//(trim(adjustl(param_suffix)))
    njhalo_nm = "NJHALO"//(trim(adjustl(param_suffix)))
    layout_nm = "LAYOUT"//(trim(adjustl(param_suffix)))
    io_layout_nm = "IO_LAYOUT"//(trim(adjustl(param_suffix)))
    masktable_nm = "MASKTABLE"//(trim(adjustl(param_suffix)))
    niproc_nm = "NIPROC"//(trim(adjustl(param_suffix)))
    njproc_nm = "NJPROC"//(trim(adjustl(param_suffix)))
  endif ; endif

  is_static = .false. ; if (present(static_memory)) is_static = static_memory
  if (is_static) then
    if (.not.present(NIHALO)) call LND_error(FATAL, "NIHALO must be "// &
      "present in the call to LND_domains_init with static memory.")
    if (.not.present(NJHALO)) call LND_error(FATAL, "NJHALO must be "// &
      "present in the call to LND_domains_init with static memory.")
    if (.not.present(NIGLOBAL)) call LND_error(FATAL, "NIGLOBAL must be "// &
      "present in the call to LND_domains_init with static memory.")
    if (.not.present(NJGLOBAL)) call LND_error(FATAL, "NJGLOBAL must be "// &
      "present in the call to LND_domains_init with static memory.")
    if (.not.present(NIPROC)) call LND_error(FATAL, "NIPROC must be "// &
      "present in the call to LND_domains_init with static memory.")
    if (.not.present(NJPROC)) call LND_error(FATAL, "NJPROC must be "// &
      "present in the call to LND_domains_init with static memory.")
  endif

  ! Read all relevant parameters and write them to the model log.
  call get_param(param_file, mdl, "REENTRANT_X", reentrant_x, &
                 "If true, the domain is zonally reentrant.", default=.true.)
  call get_param(param_file, mdl, "REENTRANT_Y", reentrant_y, &
                 "If true, the domain is meridionally reentrant.", &
                 default=.false.)
  call get_param(param_file, mdl, "TRIPOLAR_N", tripolar_N, &
                 "Use tripolar connectivity at the northern edge of the "//&
                 "domain.  With TRIPOLAR_N, NIGLOBAL must be even.", &
                 default=.false.)

#ifndef NOT_SET_AFFINITY
!$  call fms_affinity_init
!$OMP PARALLEL
!$OMP master
!$ land_nthreads = omp_get_num_threads()
!$OMP END MASTER
!$OMP END PARALLEL
!$ if(land_nthreads < 2 ) then
!$   call get_param(param_file, mdl, "LAND_OMP_THREADS", land_nthreads, &
!$              "The number of OpenMP threads that UCLAND will use.", &
!$              default = 1, layoutParam=.true.)
!$   call get_param(param_file, mdl, "LAND_OMP_HYPER_THREAD", land_omp_hyper_thread, &
!$              "If True, use hyper-threading.", default = .false., layoutParam=.true.)
!$   call fms_affinity_set('LAND', land_omp_hyper_thread, land_nthreads)
!$   call omp_set_num_threads(land_nthreads)
!$   write(6,*) "LND_domains_mod OMPthreading ", fms_affinity_get(), omp_get_thread_num(), omp_get_num_threads()
!$   call flush(6)
!$ endif
#endif
  call log_param(param_file, mdl, "!SYMMETRIC_MEMORY_", LND_dom%symmetric, &
                 "If defined, the velocity point data domain includes "//&
                 "every face of the thickness points. In other words, "//&
                 "some arrays are larger than others, depending on where "//&
                 "they are on the staggered grid.  Also, the starting "//&
                 "index of the velocity-point arrays is usually 0, not 1. "//&
                 "This can only be set at compile time.",&
                 layoutParam=.true.)
  call get_param(param_file, mdl, "NONBLOCKING_UPDATES", LND_dom%nonblocking_updates, &
                 "If true, non-blocking halo updates may be used.", &
                 default=.false., layoutParam=.true.)
  call get_param(param_file, mdl, "THIN_HALO_UPDATES", LND_dom%thin_halo_updates, &
                 "If true, optional arguments may be used to specify the "//&
                 "the width of the halos that are updated with each call.", &
                 default=.true., layoutParam=.true.)

  nihalo_dflt = 4 ; njhalo_dflt = 4
  if (present(NIHALO)) nihalo_dflt = NIHALO
  if (present(NJHALO)) njhalo_dflt = NJHALO

  call log_param(param_file, mdl, "!STATIC_MEMORY_", is_static, &
                 "If STATIC_MEMORY_ is defined, the principle variables "//&
                 "will have sizes that are statically determined at "//&
                 "compile time.  Otherwise the sizes are not determined "//&
                 "until run time. The STATIC option is substantially "//&
                 "faster, but does not allow the PE count to be changed "//&
                 "at run time.  This can only be set at compile time.",&
                 layoutParam=.true.)

  if (is_static) then
    call get_param(param_file, mdl, "NIGLOBAL", LND_dom%niglobal, &
                 "The total number of thickness grid points in the "//&
                 "x-direction in the physical domain. With STATIC_MEMORY_ "//&
                 "this is set in "//trim(inc_nm)//" at compile time.", &
                 static_value=NIGLOBAL)
    call get_param(param_file, mdl, "NJGLOBAL", LND_dom%njglobal, &
                 "The total number of thickness grid points in the "//&
                 "y-direction in the physical domain. With STATIC_MEMORY_ "//&
                 "this is set in "//trim(inc_nm)//" at compile time.", &
                 static_value=NJGLOBAL)
    if (LND_dom%niglobal /= NIGLOBAL) call LND_error(FATAL,"LND_domains_init: " // &
     "static mismatch for NIGLOBAL_ domain size. Header file does not match input namelist")
    if (LND_dom%njglobal /= NJGLOBAL) call LND_error(FATAL,"LND_domains_init: " // &
     "static mismatch for NJGLOBAL_ domain size. Header file does not match input namelist")

  else
    call get_param(param_file, mdl, "NIGLOBAL", LND_dom%niglobal, &
                 "The total number of thickness grid points in the "//&
                 "x-direction in the physical domain. With STATIC_MEMORY_ "//&
                 "this is set in "//trim(inc_nm)//" at compile time.", &
                 fail_if_missing=.true.)
    call get_param(param_file, mdl, "NJGLOBAL", LND_dom%njglobal, &
                 "The total number of thickness grid points in the "//&
                 "y-direction in the physical domain. With STATIC_MEMORY_ "//&
                 "this is set in "//trim(inc_nm)//" at compile time.", &
                 fail_if_missing=.true.)
  endif

  call get_param(param_file, mdl, trim(nihalo_nm), LND_dom%nihalo, &
                 "The number of halo points on each side in the x-direction. How this is set "//&
                 "varies with the calling component and static or dynamic memory configuration.", &
                 default=nihalo_dflt, static_value=nihalo_dflt)
  call get_param(param_file, mdl, trim(njhalo_nm), LND_dom%njhalo, &
                 "The number of halo points on each side in the y-direction. How this is set "//&
                 "varies with the calling component and static or dynamic memory configuration.", &
                 default=njhalo_dflt, static_value=njhalo_dflt)
  if (present(min_halo)) then
    LND_dom%nihalo = max(LND_dom%nihalo, min_halo(1))
    min_halo(1) = LND_dom%nihalo
    LND_dom%njhalo = max(LND_dom%njhalo, min_halo(2))
    min_halo(2) = LND_dom%njhalo
    ! These are generally used only with static memory, so they are considerd
    ! layout params.
    call log_param(param_file, mdl, "!NIHALO min_halo", LND_dom%nihalo, layoutParam=.true.)
    call log_param(param_file, mdl, "!NJHALO min_halo", LND_dom%nihalo, layoutParam=.true.)
  endif
  if (is_static .and. .not.present(min_halo)) then
    if (LND_dom%nihalo /= NIHALO) call LND_error(FATAL,"LND_domains_init: " // &
           "static mismatch for "//trim(nihalo_nm)//" domain size")
    if (LND_dom%njhalo /= NJHALO) call LND_error(FATAL,"LND_domains_init: " // &
           "static mismatch for "//trim(njhalo_nm)//" domain size")
  endif

  global_indices(1) = 1 ; global_indices(2) = LND_dom%niglobal
  global_indices(3) = 1 ; global_indices(4) = LND_dom%njglobal

  call get_param(param_file, mdl, "INPUTDIR", inputdir, do_not_log=.true., default=".")
  inputdir = slasher(inputdir)

  call get_param(param_file, mdl, trim(masktable_nm), mask_table, &
                 "A text file to specify n_mask, layout and mask_list. "//&
                 "This feature masks out processors that contain only land points. "//&
                 "The first line of mask_table is the number of regions to be masked out. "//&
                 "The second line is the layout of the model and must be "//&
                 "consistent with the actual model layout. "//&
                 "The following (n_mask) lines give the logical positions "//&
                 "of the processors that are masked out. The mask_table "//&
                 "can be created by tools like check_mask. The "//&
                 "following example of mask_table masks out 2 processors, "//&
                 "(1,2) and (3,6), out of the 24 in a 4x6 layout: \n"//&
                 " 2\n 4,6\n 1,2\n 3,6\n", default="LND_mask_table", &
                 layoutParam=.true.)
  mask_table = trim(inputdir)//trim(mask_table)
  mask_table_exists = file_exist(mask_table)

  if (is_static) then
    layout(1) = NIPROC ; layout(2) = NJPROC
  else
    call get_param(param_file, mdl, trim(layout_nm), layout, &
                 "The processor layout to be used, or 0, 0 to automatically "//&
                 "set the layout based on the number of processors.", default=0, &
                 do_not_log=.true.)
    call get_param(param_file, mdl, trim(niproc_nm), nip_parsed, &
                 "The number of processors in the x-direction.", default=-1, &
                 do_not_log=.true.)
    call get_param(param_file, mdl, trim(njproc_nm), njp_parsed, &
                 "The number of processors in the y-direction.", default=-1, &
                 do_not_log=.true.)
    if (nip_parsed > -1) then
      if ((layout(1) > 0) .and. (layout(1) /= nip_parsed)) &
        call LND_error(FATAL, trim(layout_nm)//" and "//trim(niproc_nm)//" set inconsistently. "//&
                              "Only LAYOUT should be used.")
      layout(1) = nip_parsed
      call LND_mesg(trim(niproc_nm)//" used to set "//trim(layout_nm)//" in dynamic mode.  "//&
                    "Shift to using "//trim(layout_nm)//" instead.")
    endif
    if (njp_parsed > -1) then
      if ((layout(2) > 0) .and. (layout(2) /= njp_parsed)) &
        call LND_error(FATAL, trim(layout_nm)//" and "//trim(njproc_nm)//" set inconsistently. "//&
                              "Only "//trim(layout_nm)//" should be used.")
      layout(2) = njp_parsed
      call LND_mesg(trim(njproc_nm)//" used to set "//trim(layout_nm)//" in dynamic mode.  "//&
                    "Shift to using "//trim(layout_nm)//" instead.")
    endif

    if ( layout(1)==0 .and. layout(2)==0 ) &
      call mpp_define_layout(global_indices, proc_used, layout)
    if ( layout(1)/=0 .and. layout(2)==0 ) layout(2) = proc_used/layout(1)
    if ( layout(1)==0 .and. layout(2)/=0 ) layout(1) = proc_used/layout(2)

    if (layout(1)*layout(2) /= proc_used .and. (.not. mask_table_exists) ) then
      write(mesg,'("LND_domains_init: The product of the two components of layout, ", &
            &      2i4,", is not the number of PEs used, ",i5,".")') &
            layout(1),layout(2),proc_used
      call LND_error(FATAL, mesg)
    endif
  endif
  call log_param(param_file, mdl, trim(niproc_nm), layout(1), &
                 "The number of processors in the x-direction. With "//&
                 "STATIC_MEMORY_ this is set in "//trim(inc_nm)//" at compile time.",&
                 layoutParam=.true.)
  call log_param(param_file, mdl, trim(njproc_nm), layout(2), &
                 "The number of processors in the y-direction. With "//&
                 "STATIC_MEMORY_ this is set in "//trim(inc_nm)//" at compile time.",&
                 layoutParam=.true.)
  call log_param(param_file, mdl, trim(layout_nm), layout, &
                 "The processor layout that was actually used.",&
                 layoutParam=.true.)

  ! Idiot check that fewer PEs than columns have been requested
  if (layout(1)*layout(2)>LND_dom%niglobal*LND_dom%njglobal)  then
    write(mesg,'(a,2(i5,x,a))') 'You requested to use',layout(1)*layout(2), &
      'PEs but there are only',LND_dom%niglobal*LND_dom%njglobal,'columns in the model'
    call LND_error(FATAL, mesg)
  endif

  if (mask_table_exists) then
    call LND_error(NOTE, 'LND_domains_init: reading maskmap information from '//&
                         trim(mask_table))
    allocate(LND_dom%maskmap(layout(1), layout(2)))
    call parse_mask_table(mask_table, LND_dom%maskmap, dom_name)
  endif

  !   Set up the I/O layout, and check that it uses an even multiple of the
  ! number of PEs in each direction.
  io_layout(:) = (/ 1, 1 /)
  call get_param(param_file, mdl, trim(io_layout_nm), io_layout, &
                 "The processor layout to be used, or 0,0 to automatically "//&
                 "set the io_layout to be the same as the layout.", default=1, &
                 layoutParam=.true.)

  if (io_layout(1) < 0) then
    write(mesg,'("LND_domains_init: IO_LAYOUT(1) = ",i4,".  Negative values "//&
         &"are not allowed in ")') io_layout(1)
    call LND_error(FATAL, mesg//trim(IO_layout_nm))
  elseif (io_layout(1) > 0) then ; if (modulo(layout(1), io_layout(1)) /= 0) then
    write(mesg,'("LND_domains_init: The i-direction I/O-layout, IO_LAYOUT(1)=",i4, &
         &", does not evenly divide the i-direction layout, NIPROC=,",i4,".")') &
          io_layout(1),layout(1)
    call LND_error(FATAL, mesg)
  endif ; endif

  if (io_layout(2) < 0) then
    write(mesg,'("LND_domains_init: IO_LAYOUT(2) = ",i4,".  Negative values "//&
         &"are not allowed in ")') io_layout(2)
    call LND_error(FATAL, mesg//trim(IO_layout_nm))
  elseif (io_layout(2) /= 0) then ; if (modulo(layout(2), io_layout(2)) /= 0) then
    write(mesg,'("LND_domains_init: The j-direction I/O-layout, IO_LAYOUT(2)=",i4, &
         &", does not evenly divide the j-direction layout, NJPROC=,",i4,".")') &
          io_layout(2),layout(2)
    call LND_error(FATAL, mesg)
  endif ; endif

  if (io_layout(2) == 0) io_layout(2) = layout(2)
  if (io_layout(1) == 0) io_layout(1) = layout(1)

  X_FLAGS = 0 ; Y_FLAGS = 0
  if (reentrant_x) X_FLAGS = CYCLIC_GLOBAL_DOMAIN
  if (reentrant_y) Y_FLAGS = CYCLIC_GLOBAL_DOMAIN
  if (tripolar_N) then
    Y_FLAGS = FOLD_NORTH_EDGE
    if (reentrant_y) call LND_error(FATAL,"LND_domains: "// &
      "TRIPOLAR_N and REENTRANT_Y may not be defined together.")
  endif

  global_indices(1) = 1 ; global_indices(2) = LND_dom%niglobal
  global_indices(3) = 1 ; global_indices(4) = LND_dom%njglobal

  if (mask_table_exists) then
    call LND_define_domain( global_indices, layout, LND_dom%mpp_domain, &
                xflags=X_FLAGS, yflags=Y_FLAGS, &
                xhalo=LND_dom%nihalo, yhalo=LND_dom%njhalo, &
                symmetry = LND_dom%symmetric, name=dom_name, &
                maskmap=LND_dom%maskmap )
  else
    call LND_define_domain( global_indices, layout, LND_dom%mpp_domain, &
                xflags=X_FLAGS, yflags=Y_FLAGS, &
                xhalo=LND_dom%nihalo, yhalo=LND_dom%njhalo, &
                symmetry = LND_dom%symmetric, name=dom_name)
  endif

  if ((io_layout(1) > 0) .and. (io_layout(2) > 0) .and. &
      (layout(1)*layout(2) > 1)) then
    call LND_define_io_domain(LND_dom%mpp_domain, io_layout)
  endif

! Save the extra data for creating other domains of different resolution that
! overlay this domain
  LND_dom%X_FLAGS = X_FLAGS
  LND_dom%Y_FLAGS = Y_FLAGS
  LND_dom%layout = layout
  LND_dom%io_layout = io_layout

  if (is_static) then
  !   A requirement of equal sized compute domains is necessary when
  !   STATIC_MEMORY_
  ! is used.
    call mpp_get_compute_domain(LND_dom%mpp_domain,isc,iec,jsc,jec)
    xsiz = iec - isc + 1
    ysiz = jec - jsc + 1
    if (xsiz*NIPROC /= LND_dom%niglobal .OR. ysiz*NJPROC /= LND_dom%njglobal) then
       write( char_xsiz,'(i4)' ) NIPROC
       write( char_ysiz,'(i4)' ) NJPROC
       write( char_niglobal,'(i4)' ) LND_dom%niglobal
       write( char_njglobal,'(i4)' ) LND_dom%njglobal
       call LND_error(WARNING,'LND_domains: Processor decomposition (NIPROC_,NJPROC_) = (' &
           //trim(char_xsiz)//','//trim(char_ysiz)// &
           ') does not evenly divide size set by preprocessor macro ('&
           //trim(char_niglobal)//','//trim(char_njglobal)// '). ')
       call LND_error(FATAL,'LND_domains:  #undef STATIC_MEMORY_ in "//trim(inc_nm)//" to use &
           &dynamic allocation, or change processor decomposition to evenly divide the domain.')
    endif
  endif

  global_indices(1) = 1 ; global_indices(2) = int(LND_dom%niglobal/2)
  global_indices(3) = 1 ; global_indices(4) = int(LND_dom%njglobal/2)
  !For downsampled domain, recommend a halo of 1 (or 0?) since we're not doing
  !wide-stencil computations.
  !But that does not work because the downsampled field would not have the
  !correct size to pass the checks, e.g., we get
  !error: downsample_diag_indices_get: peculiar size 28 in i-direction\ndoes not
  !match one of 24 25 26 27
  xhalo_d2 = int(LND_dom%nihalo/2)
  yhalo_d2 = int(LND_dom%njhalo/2)
  if (mask_table_exists) then
    call LND_define_domain( global_indices, layout, LND_dom%mpp_domain_d2, &
                xflags=X_FLAGS, yflags=Y_FLAGS, &
                xhalo=xhalo_d2, yhalo=yhalo_d2, &
                symmetry = LND_dom%symmetric, name=trim("LNDc"), &
                maskmap=LND_dom%maskmap )
  else
    call LND_define_domain( global_indices, layout, LND_dom%mpp_domain_d2, &
                xflags=X_FLAGS, yflags=Y_FLAGS, &
                xhalo=xhalo_d2, yhalo=yhalo_d2, &
                symmetry = LND_dom%symmetric, name=trim("LNDc"))
  endif

  if ((io_layout(1) > 0) .and. (io_layout(2) > 0) .and. &
      (layout(1)*layout(2) > 1)) then
    call LND_define_io_domain(LND_dom%mpp_domain_d2, io_layout)
  endif

end subroutine LND_domains_init

end module LND_domains
