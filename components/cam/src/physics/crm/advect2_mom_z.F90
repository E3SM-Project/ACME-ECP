module advect2_mom_z_mod
  implicit none

contains

  subroutine advect2_mom_z(ncrms)
    !       momentum tendency due to the 2nd-order-central vertical advection
    use vars
    use params, only: crm_rknd
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) fuz(nx,ny,nz),fvz(nx,ny,nz),fwz(nx,ny,nzm)
    integer i, j, k, kc, kb,icrm
    real(crm_rknd) dz2, dz25, www, rhoi

    do icrm = 1 , ncrms

      dz25=1./(4.*dz(icrm))
      dz2=dz25*2.

      do j=1,ny
        do i=1,nx
          fuz(i,j,1) = 0.
          fvz(i,j,1) = 0.
          fuz(i,j,nz) = 0.
          fvz(i,j,nz) = 0.
          fwz(i,j,1) = 0.
          fwz(i,j,nzm) = 0.
        end do
      end do

      uwle(1,icrm) = 0.
      vwle(1,icrm) = 0.

      if(RUN3D) then

        do k=2,nzm
          kb = k-1
          rhoi = dz25 * rhow(k,icrm)
          uwle(k,icrm) = 0.
          vwle(k,icrm) = 0.
          do j=1,ny
            do i=1,nx
              fuz(i,j,k) = rhoi*(w(i,j,k,icrm)+w(i-1,j,k,icrm))*(u(i,j,k,icrm)+u(i,j,kb,icrm))
              fvz(i,j,k) = rhoi*(w(i,j,k,icrm)+w(i,j-1,k,icrm))*(v(i,j,k,icrm)+v(i,j,kb,icrm))
              uwle(k,icrm) = uwle(k,icrm)+fuz(i,j,k)
              vwle(k,icrm) = vwle(k,icrm)+fvz(i,j,k)
            end do
          end do
        end do

      else

        do k=2,nzm
          kb = k-1
          rhoi = dz25 * rhow(k,icrm)
          uwle(k,icrm) = 0.
          vwle(k,icrm) = 0.
          do j=1,ny
            do i=1,nx
              www = rhoi*(w(i,j,k,icrm)+w(i-1,j,k,icrm))
              fuz(i,j,k) = www*(u(i,j,k,icrm)+u(i,j,kb,icrm))
              fvz(i,j,k) = www*(v(i,j,k,icrm)+v(i,j,kb,icrm))
              uwle(k,icrm) = uwle(k,icrm)+fuz(i,j,k)
              vwle(k,icrm) = vwle(k,icrm)+fvz(i,j,k)
            end do
          end do
        end do


      endif

      do k=1,nzm
        kc = k+1
        rhoi = 1./(rho(k,icrm)*adz(k,icrm))
        do j=1,ny
          do i=1,nx
            dudt(i,j,k,na(icrm),icrm)=dudt(i,j,k,na(icrm),icrm)-(fuz(i,j,kc)-fuz(i,j,k))*rhoi
            dvdt(i,j,k,na(icrm),icrm)=dvdt(i,j,k,na(icrm),icrm)-(fvz(i,j,kc)-fvz(i,j,k))*rhoi
            fwz(i,j,k)=dz25*(w(i,j,kc,icrm)*rhow(kc,icrm)+w(i,j,k,icrm)*rhow(k,icrm))*(w(i,j,kc,icrm)+w(i,j,k,icrm))
          end do
        end do
      end do

      do k=2,nzm
        kb=k-1
        rhoi = 1./(rhow(k,icrm)*adzw(k,icrm))
        do j=1,ny
          do i=1,nx
            dwdt(i,j,k,na(icrm),icrm)=dwdt(i,j,k,na(icrm),icrm)-(fwz(i,j,k)-fwz(i,j,kb))*rhoi
          end do
        end do
      end do ! k
    enddo

  end subroutine advect2_mom_z

end module advect2_mom_z_mod
