module advect2_mom_z_mod
  use params, only: asyncid
  implicit none

contains

  subroutine advect2_mom_z(ncrms)
    !       momentum tendency due to the 2nd-order-central vertical advection
    use vars
    use params, only: crm_rknd
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd), allocatable :: fuz(:,:,:,:)
    real(crm_rknd), allocatable :: fvz(:,:,:,:)
    real(crm_rknd), allocatable :: fwz(:,:,:,:)
    integer i, j, k, kc, kb,icrm
    real(crm_rknd) dz25, www, rhoi

    allocate( fuz(ncrms,nx,ny,nz ) )
    allocate( fvz(ncrms,nx,ny,nz ) )
    allocate( fwz(ncrms,nx,ny,nzm) )

    !$acc enter data create(fuz,fvz,fwz) async(asyncid)

    !$acc parallel loop collapse(2) copy(uwle,vwle) async(asyncid)
    do k = 1 , nz
      do icrm = 1 , ncrms
        uwle(icrm,k) = 0.
        vwle(icrm,k) = 0.
      enddo
    enddo

    !$acc parallel loop collapse(3) copyin(dz) copy(fuz,fvz,fwz) async(asyncid)
    do j=1,ny
      do i=1,nx
        do icrm = 1 , ncrms
          dz25=1./(4.*dz(icrm))
          fuz(icrm,i,j,1  ) = 0.
          fuz(icrm,i,j,nz ) = 0.
          fvz(icrm,i,j,1  ) = 0.
          fvz(icrm,i,j,nz ) = 0.
          fwz(icrm,i,j,1  ) = 0.
          fwz(icrm,i,j,nzm) = 0.
        end do
      end do
    enddo

    if(RUN3D) then

      !$acc parallel loop collapse(4) copyin(rhow,dz,w,u,v) copy(fuz,fvz,uwle,vwle) async(asyncid)
      do k=2,nzm
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              dz25=1./(4.*dz(icrm))
              kb = k-1
              rhoi = dz25 * rhow(icrm,k)
              fuz(icrm,i,j,k) = rhoi*(w(icrm,i,j,k)+w(icrm,i-1,j  ,k))*(u(icrm,i,j,k)+u(icrm,i,j,kb))
              fvz(icrm,i,j,k) = rhoi*(w(icrm,i,j,k)+w(icrm,i  ,j-1,k))*(v(icrm,i,j,k)+v(icrm,i,j,kb))
              !$acc atomic update
              uwle(icrm,k) = uwle(icrm,k)+fuz(icrm,i,j,k)
              !$acc atomic update
              vwle(icrm,k) = vwle(icrm,k)+fvz(icrm,i,j,k)
            end do
          end do
        end do
      end do

    else

      !$acc parallel loop collapse(4) copyin(dz,rhow,w,u,v) copy(fuz,fvz,uwle,vwle) async(asyncid)
      do k=2,nzm
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              dz25=1./(4.*dz(icrm))
              kb = k-1
              rhoi = dz25 * rhow(icrm,k)
              www = rhoi*(w(icrm,i,j,k)+w(icrm,i-1,j,k))
              fuz(icrm,i,j,k) = www*(u(icrm,i,j,k)+u(icrm,i,j,kb))
              fvz(icrm,i,j,k) = www*(v(icrm,i,j,k)+v(icrm,i,j,kb))
              !$acc atomic update
              uwle(icrm,k) = uwle(icrm,k)+fuz(icrm,i,j,k)
              !$acc atomic update
              vwle(icrm,k) = vwle(icrm,k)+fvz(icrm,i,j,k)
            end do
          end do
        end do
      end do

    endif

    !$acc parallel loop collapse(4) copyin(dz,rho,adz,fuz,fvz,rhow,w) copy(dudt,dvdt,fwz) async(asyncid)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            dz25=1./(4.*dz(icrm))
            kc = k+1
            rhoi = 1./(rho(icrm,k)*adz(icrm,k))
            dudt(icrm,i,j,k,na)=dudt(icrm,i,j,k,na)-(fuz(icrm,i,j,kc)-fuz(icrm,i,j,k))*rhoi
            dvdt(icrm,i,j,k,na)=dvdt(icrm,i,j,k,na)-(fvz(icrm,i,j,kc)-fvz(icrm,i,j,k))*rhoi
            fwz(icrm,i,j,k)=dz25*(w(icrm,i,j,kc)*rhow(icrm,kc)+w(icrm,i,j,k)*rhow(icrm,k))*(w(icrm,i,j,kc)+w(icrm,i,j,k))
          end do
        end do
      end do
    end do

    !$acc parallel loop collapse(4) copyin(fwz,adzw,rhow) copy(dwdt) async(asyncid)
    do k=2,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            kb=k-1
            rhoi = 1./(rhow(icrm,k)*adzw(icrm,k))
            dwdt(icrm,i,j,k,na)=dwdt(icrm,i,j,k,na)-(fwz(icrm,i,j,k)-fwz(icrm,i,j,kb))*rhoi
          end do
        end do
      end do ! k
    end do

    !$acc exit data delete(fuz,fvz,fwz) async(asyncid)

    deallocate( fuz )
    deallocate( fvz )
    deallocate( fwz )

  end subroutine advect2_mom_z

end module advect2_mom_z_mod
