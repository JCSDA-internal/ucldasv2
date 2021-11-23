! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module ucldasv2_ucland

use fckit_mpi_module, only: fckit_mpi_comm
use mpp_mod,    only : mpp_init
use fms_io_mod, only : fms_io_init, fms_io_exit
use fms_mod,    only : read_data, write_data, fms_init, fms_end
use time_manager_mod,         only: time_type

use kinds, only: kind_real

use LND_domains,         only : LND_domain_type, LND_domains_init
use LND_error_handler,   only : LND_error, LND_mesg, WARNING, FATAL, is_root_pe
use LND_file_parser,     only : get_param, param_file_type, close_param_file
use LND_get_input,       only : directories, Get_LND_Input, directories
use LND_string_functions,only : uppercase
use LND_time_manager,    only : time_type, set_date, get_date, &
                                real_to_time, time_type_to_real, &
                                operator(+), operator(-), operator(*), operator(/), &
                                operator(>), operator(<), operator(>=), &
                                increment_date, set_calendar_type, month_name, &
                                JULIAN, GREGORIAN, NOLEAP, THIRTY_DAY_MONTHS, &
                                NO_CALENDAR

implicit none

private
public :: ucldasv2_geomdomain_init

contains

! ------------------------------------------------------------------------------
!> Initialize ucland's domain
subroutine ucldasv2_geomdomain_init(Domain, nk, f_comm)
  type(LND_domain_type), pointer, intent(in) :: Domain !< Land model domain
  integer, intent(out)                       :: nk
  type(fckit_mpi_comm),           intent(in) :: f_comm

  type(param_file_type) :: param_file                !< Structure to parse for run-time parameters
  type(directories)     :: dirs                      !< Structure containing several relevant directory paths
  character(len=40)  :: mod_name = "ucldasv2_ucland" ! This module's name.

  call mpp_init(localcomm=f_comm%communicator())

  ! Initialize fms
  call fms_init()

  ! Initialize fms io
  call fms_io_init()

  ! Parse grid inputs
  call Get_LND_Input(param_file, dirs)

  ! Domain decomposition/Inintialize mpp domains
  call LND_domains_init(Domain, param_file)

  ! Get number of levels
  call get_param(param_file, mod_name, "NSOIL", nk, fail_if_missing=.true.)

  call close_param_file(param_file)
  call fms_io_exit()

end subroutine ucldasv2_geomdomain_init

end module ucldasv2_ucland
