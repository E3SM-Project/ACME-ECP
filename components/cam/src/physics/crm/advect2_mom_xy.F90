module advect2_mom_xy_mod
  implicit none

contains

  subroutine advect2_mom_xy(ncrms,icrm)

    !        momentum tendency due to 2nd-order-central horizontal advection

    use vars
    use params, only: crm_rknd

    implicit none
    integer, intent(in) :: ncrms,icrm

    real(crm_rknd) fu(0:nx,1-YES3D:ny,nzm)
    real(crm_rknd) fv(0:nx,1-YES3D:ny,nzm)
    real(crm_rknd) fw(0:nx,1-YES3D:ny,nzm)
    real(crm_rknd) dx25, dy25, irho

    integer i, j, k, kc, kcu, ic, jb, ib, jc

    dx25 = 0.25 / dx
    dy25 = 0.25 / dy


    if(RUN3D) then

      do k = 1,nzm
        kc= k+1
        kcu =min(kc, nzm)
        irho = 1./(rhow(icrm,kc)*adzw(icrm,kc))

        do j = 1, ny
          jb = j-1
          do i = 0, nx
            ic = i+1
            fu(i,j,k)=dx25*(u(icrm,ic,j,k)+u(icrm,i,j,k))*(u(icrm,i,j,k)+u(icrm,ic,j,k))
            fv(i,j,k)=dx25*(u(icrm,ic,j,k)+u(icrm,ic,jb,k))*(v(icrm,i,j,k)+v(icrm,ic,j,k))
            fw(i,j,k)=dx25*(u(icrm,ic,j,k)*rho(icrm,k)*adz(icrm,k)+ &
            u(icrm,ic,j,kcu)*rho(icrm,kcu)*adz(icrm,kcu))*(w(icrm,i,j,kc)+w(icrm,ic,j,kc))
          end do
          do i = 1, nx
            ib = i-1
            dudt(icrm,i,j,k,na(icrm))  = dudt(icrm,i,j,k,na(icrm))  - (fu(i,j,k)-fu(ib,j,k))
            dvdt(icrm,i,j,k,na(icrm))  = dvdt(icrm,i,j,k,na(icrm))  - (fv(i,j,k)-fv(ib,j,k))
            dwdt(icrm,i,j,kc,na(icrm)) = dwdt(icrm,i,j,kc,na(icrm))-irho*(fw(i,j,k)-fw(ib,j,k))
          end do
        end do

        do j = 0, ny
          jc = j+1
          do i = 1, nx
            ib = i-1
            fu(i,j,k)=dy25*(v(icrm,i,jc,k)+v(icrm,ib,jc,k))*(u(icrm,i,j,k)+u(icrm,i,jc,k))
            fv(i,j,k)=dy25*(v(icrm,i,jc,k)+v(icrm,i,j,k))*(v(icrm,i,j,k)+v(icrm,i,jc,k))
            fw(i,j,k)=dy25*(v(icrm,i,jc,k)*rho(icrm,k)*adz(icrm,k)+ &
            v(icrm,i,jc,kcu)*rho(icrm,kcu)*adz(icrm,kcu))*(w(icrm,i,j,kc)+w(icrm,i,jc,kc))
          end do
        end do
        do j = 1,ny
          jb = j-1
          do i = 1, nx
            dudt(icrm,i,j,k,na(icrm)) = dudt(icrm,i,j,k,na(icrm)) - (fu(i,j,k) - fu(i,jb,k))
            dvdt(icrm,i,j,k,na(icrm)) = dvdt(icrm,i,j,k,na(icrm)) - (fv(i,j,k) - fv(i,jb,k))
            dwdt(icrm,i,j,kc,na(icrm))= dwdt(icrm,i,j,kc,na(icrm))-irho*(fw(i,j,k)-fw(i,jb,k))
          end do
        end do

      end do ! k


    else

      j=1

      do k = 1,nzm
        kc= k+1
        kcu =min(kc, nzm)
        irho = 1./(rhow(icrm,kc)*adzw(icrm,kc))

        do i = 0, nx
          ic = i+1
          fu(i,j,k)=dx25*(u(icrm,ic,j,k)+u(icrm,i,j,k))*(u(icrm,i,j,k)+u(icrm,ic,j,k))
          fv(i,j,k)=dx25*(u(icrm,ic,j,k)+u(icrm,i,j,k))*(v(icrm,i,j,k)+v(icrm,ic,j,k))
          fw(i,j,k)=dx25*(u(icrm,ic,j,k)*rho(icrm,k)*adz(icrm,k)+ &
          u(icrm,ic,j,kcu)*rho(icrm,kcu)*adz(icrm,kcu))*(w(icrm,i,j,kc)+w(icrm,ic,j,kc))
        end do
        do i = 1, nx
          ib = i-1
          dudt(icrm,i,j,k,na(icrm))  = dudt(icrm,i,j,k,na(icrm))  - (fu(i,j,k)-fu(ib,j,k))
          dvdt(icrm,i,j,k,na(icrm))  = dvdt(icrm,i,j,k,na(icrm))  - (fv(i,j,k)-fv(ib,j,k))
          dwdt(icrm,i,j,kc,na(icrm)) = dwdt(icrm,i,j,kc,na(icrm))-irho*(fw(i,j,k)-fw(ib,j,k))
        end do

      end do ! k

    endif

  end subroutine advect2_mom_xy

end module advect2_mom_xy_mod
