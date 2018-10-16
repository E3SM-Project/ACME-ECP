
#include "directives.inc"

module zero_mod
  implicit none

contains

  subroutine zero(ncrms)
    use vars
    use microphysics, only : total_water
    implicit none
    integer, intent(in) :: ncrms
    integer k,icrm, j, i
    !dudt(nxp1, ny  , nzm, 3, ncrms)
    !dvdt(nx  , nyp1, nzm, 3, ncrms)
    !dwdt(nx  , ny  , nz , 3, ncrms)
    !misc(nx  , ny  , nz ,    ncrms)
    
    !_dir _par _loop _gang _vector collapse(4) _async(1)
    do icrm = 1 , ncrms
      do k = 1 , nz
        do j = 1 , nyp1
          do i = 1 , nxp1
            if (i <= nxp1 .and. j <= ny   .and. k <= nzm) dudt(i,j,k,na,icrm) = 0.
            if (i <= nx   .and. j <= nyp1 .and. k <= nzm) dvdt(i,j,k,na,icrm) = 0.
            if (i <= nx   .and. j <= ny   .and. k <= nz ) dwdt(i,j,k,na,icrm) = 0.
            if (i <= nx   .and. j <= ny   .and. k <= nz ) misc(i,j,k,   icrm) = 0.
          enddo
        enddo
      enddo
    enddo
  end

end module zero_mod
