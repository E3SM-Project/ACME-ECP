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
    dtdz = dtn/dz

    do k=1,nzm
      rhox = rho(k)*dtdx
      rhoy = rho(k)*dtdy
      rhoz = rhow(k)*dtdz
      do j=1,ny
        do i=1,nx

          dudt(icrm,i,j,k,nc) = u(icrm,i,j,k) + dt3(na) &
          *(at*dudt(icrm,i,j,k,na)+bt*dudt(icrm,i,j,k,nb)+ct*dudt(icrm,i,j,k,nc))

          dvdt(icrm,i,j,k,nc) = v(icrm,i,j,k) + dt3(na) &
          *(at*dvdt(icrm,i,j,k,na)+bt*dvdt(icrm,i,j,k,nb)+ct*dvdt(icrm,i,j,k,nc))

          dwdt(icrm,i,j,k,nc) = w(icrm,i,j,k) + dt3(na) &
          *(at*dwdt(icrm,i,j,k,na)+bt*dwdt(icrm,i,j,k,nb)+ct*dwdt(icrm,i,j,k,nc))

          u(icrm,i,j,k) = 0.5*(u(icrm,i,j,k)+dudt(icrm,i,j,k,nc)) * rhox
          v(icrm,i,j,k) = 0.5*(v(icrm,i,j,k)+dvdt(icrm,i,j,k,nc)) * rhoy
          misc(icrm,i,j,k) = 0.5*(w(icrm,i,j,k)+dwdt(icrm,i,j,k,nc))
          w(icrm,i,j,k) = 0.5*(w(icrm,i,j,k)+dwdt(icrm,i,j,k,nc)) * rhoz


        end do
      end do
    end do

  end subroutine adams

end module adams_mod
