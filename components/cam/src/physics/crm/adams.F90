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

    dtdx = dtn/dx
    dtdy = dtn/dy
    dtdz = dtn/dz(icrm)

    do k=1,nzm
      rhox = rho(k,icrm)*dtdx
      rhoy = rho(k,icrm)*dtdy
      rhoz = rhow(k,icrm)*dtdz
      do j=1,ny
        do i=1,nx

          dudt(i,j,k,nc(icrm),icrm) = u(i,j,k,icrm) + dt3(na(icrm),icrm) &
          *(at*dudt(i,j,k,na(icrm),icrm)+bt*dudt(i,j,k,nb(icrm),icrm)+ct*dudt(i,j,k,nc(icrm),icrm))

          dvdt(i,j,k,nc(icrm),icrm) = v(i,j,k,icrm) + dt3(na(icrm),icrm) &
          *(at*dvdt(i,j,k,na(icrm),icrm)+bt*dvdt(i,j,k,nb(icrm),icrm)+ct*dvdt(i,j,k,nc(icrm),icrm))

          dwdt(i,j,k,nc(icrm),icrm) = w(i,j,k,icrm) + dt3(na(icrm),icrm) &
          *(at*dwdt(i,j,k,na(icrm),icrm)+bt*dwdt(i,j,k,nb(icrm),icrm)+ct*dwdt(i,j,k,nc(icrm),icrm))

          u(i,j,k,icrm) = 0.5*(u(i,j,k,icrm)+dudt(i,j,k,nc(icrm),icrm)) * rhox
          v(i,j,k,icrm) = 0.5*(v(i,j,k,icrm)+dvdt(i,j,k,nc(icrm),icrm)) * rhoy
          misc(i,j,k,icrm) = 0.5*(w(i,j,k,icrm)+dwdt(i,j,k,nc(icrm),icrm))
          w(i,j,k,icrm) = 0.5*(w(i,j,k,icrm)+dwdt(i,j,k,nc(icrm),icrm)) * rhoz


        end do
      end do
    end do

  end subroutine adams

end module adams_mod
