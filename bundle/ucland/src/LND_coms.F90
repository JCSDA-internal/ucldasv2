!> Interfaces to non-domain-oriented communication subroutines, including the
!! UCLAND reproducing sums facility
module LND_coms

! This file is part of UCLAND. See LICENSE.md for the license.

use mpp_mod, only : PE_here => mpp_pe, root_PE => mpp_root_pe, num_PEs => mpp_npes
use mpp_mod, only : broadcast => mpp_broadcast
use mpp_mod, only : sum_across_PEs => mpp_sum, max_across_PEs => mpp_max, min_across_PEs => mpp_min

implicit none ; private

public :: PE_here, root_PE, num_PEs
public :: broadcast, sum_across_PEs, min_across_PEs, max_across_PEs

end module LND_coms
