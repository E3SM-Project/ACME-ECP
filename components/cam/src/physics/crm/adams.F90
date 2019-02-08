module adams_mod
  use params, only: asyncid
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

    !$acc parallel loop collapse(4) copyin(dz,rho,rhow,dt3) copy(dudt,dvdt,dwdt,u,v,w,misc) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            dtdz = dtn/dz(icrm)
            rhox = rho(k,icrm)*dtdx
            rhoy = rho(k,icrm)*dtdy
            rhoz = rhow(k,icrm)*dtdz
            dudt(i,j,k,nc,icrm) = u(icrm,i,j,k) + dt3(na) *(at*dudt(i,j,k,na,icrm)+bt*dudt(i,j,k,nb,icrm)+ct*dudt(i,j,k,nc,icrm))
            dvdt(i,j,k,nc,icrm) = v(icrm,i,j,k) + dt3(na) *(at*dvdt(i,j,k,na,icrm)+bt*dvdt(i,j,k,nb,icrm)+ct*dvdt(i,j,k,nc,icrm))
            dwdt(i,j,k,nc,icrm) = w(icrm,i,j,k) + dt3(na) *(at*dwdt(i,j,k,na,icrm)+bt*dwdt(i,j,k,nb,icrm)+ct*dwdt(i,j,k,nc,icrm))
            u(icrm,i,j,k) = 0.5*(u(icrm,i,j,k)+dudt(i,j,k,nc,icrm)) * rhox
            v(icrm,i,j,k) = 0.5*(v(icrm,i,j,k)+dvdt(i,j,k,nc,icrm)) * rhoy
            misc(i,j,k,icrm) = 0.5*(w(icrm,i,j,k)+dwdt(i,j,k,nc,icrm))
            w(icrm,i,j,k) = 0.5*(w(icrm,i,j,k)+dwdt(i,j,k,nc,icrm)) * rhoz
          end do
        end do
      end do
    end do

  end subroutine adams

end module adams_mod
