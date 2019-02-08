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
    real(crm_rknd) :: fuz(nx,ny,nz ,ncrms)
    real(crm_rknd) :: fvz(nx,ny,nz ,ncrms)
    real(crm_rknd) :: fwz(nx,ny,nzm,ncrms)
    integer i, j, k, kc, kb,icrm
    real(crm_rknd) dz25, www, rhoi

    !$acc enter data create(fuz,fvz,fwz) async(asyncid)

    !$acc parallel loop collapse(2) copyout(vwle,uwle) async(asyncid)
    do icrm = 1 , ncrms
      do k = 1 , nz
        uwle(k,icrm) = 0.
        vwle(k,icrm) = 0.
      enddo
    enddo

    !$acc parallel loop collapse(3) copy(fuz,fwz,fvz) async(asyncid)
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

      !$acc parallel loop collapse(4) copyin(rhow,dz,w,u,v) copy(fuz,fvz,uwle,vwle) async(asyncid)
      do icrm = 1 , ncrms
        do k=2,nzm
          do j=1,ny
            do i=1,nx
              dz25=1./(4.*dz(icrm))
              kb = k-1
              rhoi = dz25 * rhow(k,icrm)
              fuz(i,j,k,icrm) = rhoi*(w(icrm,i,j,k)+w(icrm,i-1,j  ,k))*(u(icrm,i,j,k)+u(icrm,i,j,kb))
              fvz(i,j,k,icrm) = rhoi*(w(icrm,i,j,k)+w(icrm,i  ,j-1,k))*(v(icrm,i,j,k)+v(icrm,i,j,kb))
              !$acc atomic update
              uwle(k,icrm) = uwle(k,icrm)+fuz(i,j,k,icrm)
              !$acc atomic update
              vwle(k,icrm) = vwle(k,icrm)+fvz(i,j,k,icrm)
            end do
          end do
        end do
      end do

    else

      !$acc parallel loop collapse(4) copyin(u,v,w,rhow,dz) copy(fvz,vwle,uwle,fuz) async(asyncid)
      do icrm = 1 , ncrms
        do k=2,nzm
          do j=1,ny
            do i=1,nx
              dz25=1./(4.*dz(icrm))
              kb = k-1
              rhoi = dz25 * rhow(k,icrm)
              www = rhoi*(w(icrm,i,j,k)+w(icrm,i-1,j,k))
              fuz(i,j,k,icrm) = www*(u(icrm,i,j,k)+u(icrm,i,j,kb))
              fvz(i,j,k,icrm) = www*(v(icrm,i,j,k)+v(icrm,i,j,kb))
              !$acc atomic update
              uwle(k,icrm) = uwle(k,icrm)+fuz(i,j,k,icrm)
              !$acc atomic update
              vwle(k,icrm) = vwle(k,icrm)+fvz(i,j,k,icrm)
            end do
          end do
        end do
      end do

    endif

    !$acc parallel loop collapse(4) copyin(fvz,rho,w,rhow,fuz,dz,adz) copy(dudt,dvdt,fwz) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            dz25=1./(4.*dz(icrm))
            kc = k+1
            rhoi = 1./(rho(k,icrm)*adz(k,icrm))
            dudt(icrm,i,j,k,na)=dudt(icrm,i,j,k,na)-(fuz(i,j,kc,icrm)-fuz(i,j,k,icrm))*rhoi
            dvdt(icrm,i,j,k,na)=dvdt(icrm,i,j,k,na)-(fvz(i,j,kc,icrm)-fvz(i,j,k,icrm))*rhoi
            fwz(i,j,k,icrm)=dz25*(w(icrm,i,j,kc)*rhow(kc,icrm)+w(icrm,i,j,k)*rhow(k,icrm))*(w(icrm,i,j,kc)+w(icrm,i,j,k))
          end do
        end do
      end do
    end do

    !$acc parallel loop collapse(4) copyin(rhow,fwz,adzw) copy(dwdt) async(asyncid)
    do icrm = 1 , ncrms
      do k=2,nzm
        do j=1,ny
          do i=1,nx
            kb=k-1
            rhoi = 1./(rhow(k,icrm)*adzw(icrm,k))
            dwdt(icrm,i,j,k,na)=dwdt(icrm,i,j,k,na)-(fwz(i,j,k,icrm)-fwz(i,j,kb,icrm))*rhoi
          end do
        end do
      end do ! k
    end do

    !$acc exit data delete(fuz,fvz,fwz) async(asyncid)

  end subroutine advect2_mom_z

end module advect2_mom_z_mod
