!>
! Module: mo_gas_optics

! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2018,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!
! Description: Permute arrays
!

module mo_util_reorder
  use mo_rte_kind, only: wp
  use mo_reorder_kernels, only: reorder_123x312_kernel, reorder_123x321_kernel
  implicit none
  private
  public :: reorder123x312, reorder123x321
contains

  function reorder123x312(array)
    real(wp), dimension(:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=3),size(array,dim=1),size(array,dim=2)) :: reorder123x312

    call reorder_123x312_kernel(size(array,dim=1), size(array,dim=2), size(array,dim=3), array, reorder123x312)
  end function reorder123x312

  function reorder123x321(array)
    real(wp), dimension(:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=3),size(array,dim=2),size(array,dim=1)) :: reorder123x321

    call reorder_123x321_kernel(size(array,dim=1), size(array,dim=2), size(array,dim=3), array, reorder123x321)
  end function reorder123x321

end module
