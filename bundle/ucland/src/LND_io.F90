!> This module contains I/O framework code
module LND_io

! This file is part of UCLAND. See LICENSE.md for the license.

use LND_error_handler,    only : LND_error, NOTE, FATAL, WARNING
use LND_domains,          only : LND_domain_type
use LND_string_functions, only : lowercase, slasher

use ensemble_manager_mod, only : get_ensemble_id
use fms_io_mod,           only : file_exist
use mpp_domains_mod,      only : domain2d
use mpp_io_mod,           only : open_file => mpp_open, close_file => mpp_close
use fms_mod,              only : open_namelist_file, check_nml_error

use netcdf

implicit none ; private

public :: close_file, file_exists, slasher, ensembler
public :: open_namelist_file, check_nml_error

!> Indicate whether a file exists, perhaps with domain decomposition
interface file_exists
  module procedure FMS_file_exists
  module procedure LND_file_exists
end interface

contains

!> Returns a name with "%#E" or "%E" replaced with the ensemble member number.
function ensembler(name, ens_no_in) result(en_nm)
  character(len=*),  intent(in) :: name       !< The name to be modified
  integer, optional, intent(in) :: ens_no_in  !< The number of the current ensemble member
  character(len=len(name)) :: en_nm  !< The name encoded with the ensemble number

  ! This function replaces "%#E" or "%E" with the ensemble number anywhere it
  ! occurs in name, with %E using 4 or 6 digits (depending on the ensemble size)
  ! and %#E using # digits, where # is a number from 1 to 9.

  character(len=len(name)) :: tmp
  character(10) :: ens_num_char
  character(3)  :: code_str
  integer :: ens_no
  integer :: n, is, ie

  en_nm = trim(name)
  if (index(name,"%") == 0) return

  if (present(ens_no_in)) then
    ens_no = ens_no_in
  else
    ens_no = get_ensemble_id()
  endif

  write(ens_num_char, '(I10)') ens_no ; ens_num_char = adjustl(ens_num_char)
  do
    is = index(en_nm,"%E")
    if (is == 0) exit
    if (len(en_nm) < len(trim(en_nm)) + len(trim(ens_num_char)) - 2) &
      call LND_error(FATAL, "LND_io ensembler: name "//trim(name)// &
      " is not long enough for %E expansion for ens_no "//trim(ens_num_char))
    tmp = en_nm(1:is-1)//trim(ens_num_char)//trim(en_nm(is+2:))
    en_nm = tmp
  enddo

  if (index(name,"%") == 0) return

  write(ens_num_char, '(I10.10)') ens_no
  do n=1,9 ; do
    write(code_str, '("%",I1,"E")') n

    is = index(en_nm,code_str)
    if (is == 0) exit
    if (ens_no < 10**n) then
      if (len(en_nm) < len(trim(en_nm)) + n-3) call LND_error(FATAL, &
        "LND_io ensembler: name "//trim(name)//" is not long enough for %E expansion.")
      tmp = en_nm(1:is-1)//trim(ens_num_char(11-n:10))//trim(en_nm(is+3:))
    else
      call LND_error(FATAL, "LND_io ensembler: Ensemble number is too large "//&
          "to be encoded with "//code_str//" in "//trim(name))
    endif
    en_nm = tmp
  enddo ; enddo

end function ensembler

!> Returns true if the named file or its domain-decomposed variant exists.
function LND_file_exists(filename, LND_Domain)
  character(len=*),       intent(in) :: filename   !< The name of the file being inquired about
  type(LND_domain_type),  intent(in) :: LND_Domain !< The LND_Domain that describes the decomposition

! This function uses the fms_io function file_exist to determine whether
! a named file (or its decomposed variant) exists.

  logical :: LND_file_exists

  LND_file_exists = file_exist(filename, LND_Domain%mpp_domain)

end function LND_file_exists

!> Returns true if the named file or its domain-decomposed variant exists.
function FMS_file_exists(filename, domain, no_domain)
  character(len=*), intent(in)         :: filename  !< The name of the file being inquired about
  type(domain2d), optional, intent(in) :: domain    !< The mpp domain2d that describes the decomposition
  logical,        optional, intent(in) :: no_domain !< This file does not use domain decomposition
! This function uses the fms_io function file_exist to determine whether
! a named file (or its decomposed variant) exists.

  logical :: FMS_file_exists

  FMS_file_exists = file_exist(filename, domain, no_domain)

end function FMS_file_exists

end module LND_io
