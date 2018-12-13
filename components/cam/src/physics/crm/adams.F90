module adams_mod
  implicit none

contains

  subroutine adams(ncrms)
    !       Adams-Bashforth scheme
    use vars
    use params, only: crm_rknd
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) dtdx, dtdy, dtdz, rhox, rhoy, rhoz
    integer i,j,k,icrm

    dtdx = dtn/dx
    dtdy = dtn/dy

    !$acc parallel loop collapse(4) copyin(dz,rho,rhow,dt3) copy(dudt,dvdt,dwdt,u,v,w,misc) async(1)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            dtdz = dtn/dz(icrm)
            rhox = rho(k,icrm)*dtdx
            rhoy = rho(k,icrm)*dtdy
            rhoz = rhow(k,icrm)*dtdz
            dudt(i,j,k,nc,icrm) = u(i,j,k,icrm) + dt3(na) *(at*dudt(i,j,k,na,icrm)+bt*dudt(i,j,k,nb,icrm)+ct*dudt(i,j,k,nc,icrm))
            dvdt(i,j,k,nc,icrm) = v(i,j,k,icrm) + dt3(na) *(at*dvdt(i,j,k,na,icrm)+bt*dvdt(i,j,k,nb,icrm)+ct*dvdt(i,j,k,nc,icrm))
            dwdt(i,j,k,nc,icrm) = w(i,j,k,icrm) + dt3(na) *(at*dwdt(i,j,k,na,icrm)+bt*dwdt(i,j,k,nb,icrm)+ct*dwdt(i,j,k,nc,icrm))
            u(i,j,k,icrm) = 0.5*(u(i,j,k,icrm)+dudt(i,j,k,nc,icrm)) * rhox
            v(i,j,k,icrm) = 0.5*(v(i,j,k,icrm)+dvdt(i,j,k,nc,icrm)) * rhoy
            misc(i,j,k,icrm) = 0.5*(w(i,j,k,icrm)+dwdt(i,j,k,nc,icrm))
            w(i,j,k,icrm) = 0.5*(w(i,j,k,icrm)+dwdt(i,j,k,nc,icrm)) * rhoz
          end do
        end do
      end do
    end do

  end subroutine adams

end module adams_mod
