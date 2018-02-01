module adams_mod
	implicit none

contains

  subroutine adams(ncrms,icrm)
    !       Adams-Bashforth scheme
    use vars
    use params, only: crm_rknd
    implicit none
    integer, intent(in) :: ncrms,icrm

    real(crm_rknd) dtdx, dtdy, dtdz, rhox, rhoy, rhoz
    integer i,j,k

    dtdx = dtn(icrm)/dx
    dtdy = dtn(icrm)/dy
    dtdz = dtn(icrm)/dz(icrm)

    do k=1,nzm
      rhox = rho(icrm,k)*dtdx
      rhoy = rho(icrm,k)*dtdy
      rhoz = rhow(icrm,k)*dtdz
      do j=1,ny
        do i=1,nx

          dudt(icrm,i,j,k,nc(icrm)) = u(icrm,i,j,k) + dt3(icrm,na(icrm)) &
          *(at(icrm)*dudt(icrm,i,j,k,na(icrm))+bt(icrm)*dudt(icrm,i,j,k,nb(icrm))+ct(icrm)*dudt(icrm,i,j,k,nc(icrm)))

          dvdt(icrm,i,j,k,nc(icrm)) = v(icrm,i,j,k) + dt3(icrm,na(icrm)) &
          *(at(icrm)*dvdt(icrm,i,j,k,na(icrm))+bt(icrm)*dvdt(icrm,i,j,k,nb(icrm))+ct(icrm)*dvdt(icrm,i,j,k,nc(icrm)))

          dwdt(icrm,i,j,k,nc(icrm)) = w(icrm,i,j,k) + dt3(icrm,na(icrm)) &
          *(at(icrm)*dwdt(icrm,i,j,k,na(icrm))+bt(icrm)*dwdt(icrm,i,j,k,nb(icrm))+ct(icrm)*dwdt(icrm,i,j,k,nc(icrm)))

          u(icrm,i,j,k) = 0.5*(u(icrm,i,j,k)+dudt(icrm,i,j,k,nc(icrm))) * rhox
          v(icrm,i,j,k) = 0.5*(v(icrm,i,j,k)+dvdt(icrm,i,j,k,nc(icrm))) * rhoy
          misc(icrm,i,j,k) = 0.5*(w(icrm,i,j,k)+dwdt(icrm,i,j,k,nc(icrm)))
          w(icrm,i,j,k) = 0.5*(w(icrm,i,j,k)+dwdt(icrm,i,j,k,nc(icrm))) * rhoz


        end do
      end do
    end do

  end subroutine adams

end module adams_mod
