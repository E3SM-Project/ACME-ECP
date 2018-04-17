module advect2_mom_z_mod
  implicit none

contains

  subroutine advect2_mom_z(ncrms)
    !       momentum tendency due to the 2nd-order-central vertical advection
    use vars
    use params, only: crm_rknd
    implicit none
    integer, intent(in) :: ncrms

    real(crm_rknd) fuz(ncrms,nx,ny,nz)
    real(crm_rknd) fvz(ncrms,nx,ny,nz)
    real(crm_rknd) fwz(ncrms,nx,ny,nzm)
    integer i, j, k, kc, kb,icrm
    real(crm_rknd) :: www, rhoi, tmp1, tmp2

    if(RUN3D) then

      do k=1,nz
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              fuz(icrm,i,j,k) = 0.
              fvz(icrm,i,j,k) = 0.
              fuz(icrm,i,j,k) = 0.
              fvz(icrm,i,j,k) = 0.
              fwz(icrm,i,j,k) = 0.
              fwz(icrm,i,j,k) = 0.
              if (i == 1 .and. j == 1) then
                uwle(icrm,k) = 0.
                vwle(icrm,k) = 0.
              endif
            end do
          end do
        end do
      end do

      do k=2,nzm
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              kb = k-1
              rhoi = rhow(icrm,k)/(4.*dz(icrm))
              fuz(icrm,i,j,k) = rhoi*(w(icrm,i,j,k)+w(icrm,i-1,j,k))*(u(icrm,i,j,k)+u(icrm,i,j,kb))
              fvz(icrm,i,j,k) = rhoi*(w(icrm,i,j,k)+w(icrm,i,j-1,k))*(v(icrm,i,j,k)+v(icrm,i,j,kb))
              uwle(icrm,k) = uwle(icrm,k)+fuz(icrm,i,j,k)
              vwle(icrm,k) = vwle(icrm,k)+fvz(icrm,i,j,k)
            end do
          end do
        end do
      end do

      do k=1,nzm
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              kc = k+1
              rhoi = 1./(rho(icrm,k)*adz(icrm,k))
              dudt(icrm,i,j,k,na)=dudt(icrm,i,j,k,na)-(fuz(icrm,i,j,kc)-fuz(icrm,i,j,k))*rhoi
              dvdt(icrm,i,j,k,na)=dvdt(icrm,i,j,k,na)-(fvz(icrm,i,j,kc)-fvz(icrm,i,j,k))*rhoi
              fwz(icrm,i,j,k)=(w(icrm,i,j,kc)*rhow(icrm,kc)+w(icrm,i,j,k)*rhow(icrm,k))*(w(icrm,i,j,kc)+w(icrm,i,j,k))/(4.*dz(icrm))
            end do
          end do
        end do
      end do

      do k=2,nzm
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              kb=k-1
              rhoi = 1./(rhow(icrm,k)*adzw(icrm,k))
              dwdt(icrm,i,j,k,na)=dwdt(icrm,i,j,k,na)-(fwz(icrm,i,j,k)-fwz(icrm,i,j,kb))*rhoi
            end do
          end do
        end do
      end do ! k

    else

      !$acc parallel loop gang vector collapse(2) default(present) async(1)
      do k=1,nz
        do icrm = 1 , ncrms
          uwle(icrm,k) = 0.
          vwle(icrm,k) = 0.
        end do
      end do

      !$acc parallel loop gang vector collapse(4) default(present) async(1)
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              kc = k+1
              kb = k-1
              rhoi = 1./(rho(icrm,k)*adz(icrm,k))

              tmp1 = rhow(icrm,kc)/(4.*dz(icrm))*(w(icrm,i,j,kc)+w(icrm,i-1,j,kc))*(u(icrm,i,j,kc)+u(icrm,i,j,k ))
              if (k == 1) then
                tmp2 = 0
              else
                tmp2 = rhow(icrm,k)/(4.*dz(icrm))*(w(icrm,i,j,k)+w(icrm,i-1,j,k))*(u(icrm,i,j,k)+u(icrm,i,j,kb))
              endif
              dudt(icrm,i,j,k,na)=dudt(icrm,i,j,k,na)-(tmp1-tmp2)*rhoi
              !$acc atomic update
              uwle(icrm,k) = uwle(icrm,k)+tmp2

              tmp1 = rhow(icrm,kc)/(4.*dz(icrm))*(w(icrm,i,j,kc)+w(icrm,i-1,j,kc))*(v(icrm,i,j,kc)+v(icrm,i,j,k ))
              if (k == 1) then
                tmp2 = 0
              else
                tmp2 = rhow(icrm,k)/(4.*dz(icrm))*(w(icrm,i,j,k)+w(icrm,i-1,j,k))*(v(icrm,i,j,k)+v(icrm,i,j,kb))
              endif
              dvdt(icrm,i,j,k,na)=dvdt(icrm,i,j,k,na)-(tmp1-tmp2)*rhoi
              !$acc atomic update
              vwle(icrm,k) = vwle(icrm,k)+tmp2

              if (k >= 2) then
                rhoi = 1./(rhow(icrm,k)*adzw(icrm,k))
                tmp1=(w(icrm,i,j,kc)*rhow(icrm,kc)+w(icrm,i,j,k )*rhow(icrm,k ))*(w(icrm,i,j,kc)+w(icrm,i,j,k ))/(4.*dz(icrm))
                tmp2=(w(icrm,i,j,k )*rhow(icrm,k )+w(icrm,i,j,kb)*rhow(icrm,kb))*(w(icrm,i,j,k )+w(icrm,i,j,kb))/(4.*dz(icrm))
                dwdt(icrm,i,j,k,na)=dwdt(icrm,i,j,k,na)-(tmp1-tmp2)*rhoi
              endif

            end do
          end do
        end do
      end do

    endif

  end subroutine advect2_mom_z

end module advect2_mom_z_mod
