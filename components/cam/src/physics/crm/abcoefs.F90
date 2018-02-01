
module abcoefs_mod
	implicit none

contains

	subroutine abcoefs(ncrms,icrm)
		!      coefficients for the Adams-Bashforth scheme
		use grid
		use params, only: crm_rknd
		implicit none
    integer, intent(in) :: ncrms,icrm
		real(crm_rknd) alpha, beta

		if(nstep(icrm).ge.3.and.nadams.eq.3.or.nrestart.eq.2) then
			alpha = dt3(icrm,nb(icrm)) / dt3(icrm,na(icrm))
			beta = dt3(icrm,nc(icrm)) / dt3(icrm,na(icrm))
			ct(icrm) = (2.+3.* alpha) / (6.* (alpha + beta) * beta)
			bt(icrm) = -(1.+2.*(alpha + beta) * ct(icrm))/(2. * alpha)
			at(icrm) = 1. - bt(icrm) - ct(icrm)
		else if(nstep(icrm).ge.2) then
			at(icrm) = 3./2.
			bt(icrm) = -1./2.
			ct(icrm) = 0.
		else
			at(icrm) = 1.
			bt(icrm) = 0.
			ct(icrm) = 0.
		end if

	end subroutine abcoefs

end module abcoefs_mod
