module advect2_mom_z_mod
  implicit none

contains

  subroutine advect2_mom_z(ncrms)
    !       momentum tendency due to the 2nd-order-central vertical advection
    use vars
    use params, only: crm_rknd
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) :: fuz(nx,ny,nz ,ncrms)
    real(crm_rknd) :: fvz(nx,ny,nz ,ncrms)
    real(crm_rknd) :: fwz(nx,ny,nzm,ncrms)
    integer i, j, k, kc, kb,icrm
    real(crm_rknd) dz25, www, rhoi

    !$acc enter data create(fuz,fvz,fwz) async(1)

    !$acc parallel loop collapse(2) async(1)
    do icrm = 1 , ncrms
      do k = 1 , nz
        uwle(k,icrm) = 0.
        vwle(k,icrm) = 0.
      enddo
    enddo

    !$acc parallel loop collapse(3) async(1)
    do icrm = 1 , ncrms
      do j=1,ny
        do i=1,nx
          dz25=1./(4.*dz(icrm))
          fuz(i,j,1  ,icrm) = 0.
          fuz(i,j,nz ,icrm) = 0.
          fvz(i,j,1  ,icrm) = 0.
          fvz(i,j,nz ,icrm) = 0.
          fwz(i,j,1  ,icrm) = 0.
          fwz(i,j,nzm,icrm) = 0.
        end do
      end do
    enddo

    if(RUN3D) then

      !$acc parallel loop collapse(4) async(1)
      do icrm = 1 , ncrms
        do k=2,nzm
          do j=1,ny
            do i=1,nx
              dz25=1./(4.*dz(icrm))
              kb = k-1
              rhoi = dz25 * rhow(k,icrm)
              fuz(i,j,k,icrm) = rhoi*(w(i,j,k,icrm)+w(i-1,j  ,k,icrm))*(u(i,j,k,icrm)+u(i,j,kb,icrm))
              fvz(i,j,k,icrm) = rhoi*(w(i,j,k,icrm)+w(i  ,j-1,k,icrm))*(v(i,j,k,icrm)+v(i,j,kb,icrm))
              !$acc atomic update
              uwle(k,icrm) = uwle(k,icrm)+fuz(i,j,k,icrm)
              !$acc atomic update
              vwle(k,icrm) = vwle(k,icrm)+fvz(i,j,k,icrm)
            end do
          end do
        end do
      end do

    else

      !$acc parallel loop collapse(4) async(1)
      do icrm = 1 , ncrms
        do k=2,nzm
          do j=1,ny
            do i=1,nx
              dz25=1./(4.*dz(icrm))
              kb = k-1
              rhoi = dz25 * rhow(k,icrm)
              www = rhoi*(w(i,j,k,icrm)+w(i-1,j,k,icrm))
              fuz(i,j,k,icrm) = www*(u(i,j,k,icrm)+u(i,j,kb,icrm))
              fvz(i,j,k,icrm) = www*(v(i,j,k,icrm)+v(i,j,kb,icrm))
              !$acc atomic update
              uwle(k,icrm) = uwle(k,icrm)+fuz(i,j,k,icrm)
              !$acc atomic update
              vwle(k,icrm) = vwle(k,icrm)+fvz(i,j,k,icrm)
            end do
          end do
        end do
      end do

    endif

    !$acc parallel loop collapse(4) async(1)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            dz25=1./(4.*dz(icrm))
            kc = k+1
            rhoi = 1./(rho(k,icrm)*adz(k,icrm))
            dudt(i,j,k,na,icrm)=dudt(i,j,k,na,icrm)-(fuz(i,j,kc,icrm)-fuz(i,j,k,icrm))*rhoi
            dvdt(i,j,k,na,icrm)=dvdt(i,j,k,na,icrm)-(fvz(i,j,kc,icrm)-fvz(i,j,k,icrm))*rhoi
            fwz(i,j,k,icrm)=dz25*(w(i,j,kc,icrm)*rhow(kc,icrm)+w(i,j,k,icrm)*rhow(k,icrm))*(w(i,j,kc,icrm)+w(i,j,k,icrm))
          end do
        end do
      end do
    end do

    !$acc parallel loop collapse(4) async(1)
    do icrm = 1 , ncrms
      do k=2,nzm
        do j=1,ny
          do i=1,nx
            kb=k-1
            rhoi = 1./(rhow(k,icrm)*adzw(k,icrm))
            dwdt(i,j,k,na,icrm)=dwdt(i,j,k,na,icrm)-(fwz(i,j,k,icrm)-fwz(i,j,kb,icrm))*rhoi
          end do
        end do
      end do ! k
    end do

    !$acc exit data delete(fuz,fvz,fwz) async(1)

  end subroutine advect2_mom_z

end module advect2_mom_z_mod
