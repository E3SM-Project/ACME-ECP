module adams_mod
  implicit none

contains

  subroutine adams(ncrms)
    !       Adams-Bashforth scheme
    use vars
    use params, only: crm_rknd
    implicit none
    integer, intent(in) :: ncrms

    real(crm_rknd) dtdx, dtdy, dtdz(ncrms), rhox(ncrms), rhoy(ncrms), rhoz(ncrms)
    integer i,j,k,icrm

    dtdx = dtn/dx
    dtdy = dtn/dy
    dtdz(:) = dtn/dz(:)

    do k=1,nzm
      rhox(:) = rho(:,k)*dtdx
      rhoy(:) = rho(:,k)*dtdy
      rhoz(:) = rhow(:,k)*dtdz(:)
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms

            dudt(icrm,i,j,k,nc) = u(icrm,i,j,k) + dt3(na) &
            *(at*dudt(icrm,i,j,k,na)+bt*dudt(icrm,i,j,k,nb)+ct*dudt(icrm,i,j,k,nc))

            dvdt(icrm,i,j,k,nc) = v(icrm,i,j,k) + dt3(na) &
            *(at*dvdt(icrm,i,j,k,na)+bt*dvdt(icrm,i,j,k,nb)+ct*dvdt(icrm,i,j,k,nc))

            dwdt(icrm,i,j,k,nc) = w(icrm,i,j,k) + dt3(na) &
            *(at*dwdt(icrm,i,j,k,na)+bt*dwdt(icrm,i,j,k,nb)+ct*dwdt(icrm,i,j,k,nc))

            u(icrm,i,j,k) = 0.5*(u(icrm,i,j,k)+dudt(icrm,i,j,k,nc)) * rhox(icrm)
            v(icrm,i,j,k) = 0.5*(v(icrm,i,j,k)+dvdt(icrm,i,j,k,nc)) * rhoy(icrm)
            misc(icrm,i,j,k) = 0.5*(w(icrm,i,j,k)+dwdt(icrm,i,j,k,nc))
            w(icrm,i,j,k) = 0.5*(w(icrm,i,j,k)+dwdt(icrm,i,j,k,nc)) * rhoz(icrm)


          end do
        end do
      end do
    end do

  end subroutine adams

end module adams_mod
