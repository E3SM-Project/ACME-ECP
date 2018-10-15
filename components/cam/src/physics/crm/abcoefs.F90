
#include "directives.inc"

module abcoefs_mod
  implicit none

contains

  subroutine abcoefs(ncrms)
    !      coefficients for the Adams-Bashforth scheme
    use grid
    use params, only: crm_rknd
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) alpha, beta
    integer :: icrm

    !_dir _par _loop _gang _vector _kin(dt3) _async(1)
    do icrm = 1 , ncrms
      if(nstep.ge.3.and.nadams.eq.3.or.nrestart.eq.2) then
        alpha = dt3(nb) / dt3(na)
        beta = dt3(nc) / dt3(na)
        ct = (2.+3.* alpha) / (6.* (alpha + beta) * beta)
        bt = -(1.+2.*(alpha + beta) * ct)/(2. * alpha)
        at = 1. - bt - ct
      else if(nstep.ge.2) then
        at = 3./2.
        bt = -1./2.
        ct = 0.
      else
        at = 1.
        bt = 0.
        ct = 0.
      end if
    enddo

  end subroutine abcoefs

end module abcoefs_mod
