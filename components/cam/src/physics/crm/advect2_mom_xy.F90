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
        irho = 1./(rhow(kc,icrm)*adzw(kc,icrm))

        do j = 1, ny
          jb = j-1
          do i = 0, nx
            ic = i+1
            fu(i,j,k)=dx25*(u(ic,j,k,icrm)+u(i,j,k,icrm))*(u(i,j,k,icrm)+u(ic,j,k,icrm))
            fv(i,j,k)=dx25*(u(ic,j,k,icrm)+u(ic,jb,k,icrm))*(v(i,j,k,icrm)+v(ic,j,k,icrm))
            fw(i,j,k)=dx25*(u(ic,j,k,icrm)*rho(k,icrm)*adz(k,icrm)+ &
            u(ic,j,kcu,icrm)*rho(kcu,icrm)*adz(kcu,icrm))*(w(i,j,kc,icrm)+w(ic,j,kc,icrm))
          end do
          do i = 1, nx
            ib = i-1
            dudt(i,j,k,na,icrm)  = dudt(i,j,k,na,icrm)  - (fu(i,j,k)-fu(ib,j,k))
            dvdt(i,j,k,na,icrm)  = dvdt(i,j,k,na,icrm)  - (fv(i,j,k)-fv(ib,j,k))
            dwdt(i,j,kc,na,icrm) = dwdt(i,j,kc,na,icrm)-irho*(fw(i,j,k)-fw(ib,j,k))
          end do
        end do

        do j = 0, ny
          jc = j+1
          do i = 1, nx
            ib = i-1
            fu(i,j,k)=dy25*(v(i,jc,k,icrm)+v(ib,jc,k,icrm))*(u(i,j,k,icrm)+u(i,jc,k,icrm))
            fv(i,j,k)=dy25*(v(i,jc,k,icrm)+v(i,j,k,icrm))*(v(i,j,k,icrm)+v(i,jc,k,icrm))
            fw(i,j,k)=dy25*(v(i,jc,k,icrm)*rho(k,icrm)*adz(k,icrm)+ &
            v(i,jc,kcu,icrm)*rho(kcu,icrm)*adz(kcu,icrm))*(w(i,j,kc,icrm)+w(i,jc,kc,icrm))
          end do
        end do
        do j = 1,ny
          jb = j-1
          do i = 1, nx
            dudt(i,j,k,na,icrm) = dudt(i,j,k,na,icrm) - (fu(i,j,k) - fu(i,jb,k))
            dvdt(i,j,k,na,icrm) = dvdt(i,j,k,na,icrm) - (fv(i,j,k) - fv(i,jb,k))
            dwdt(i,j,kc,na,icrm)= dwdt(i,j,kc,na,icrm)-irho*(fw(i,j,k)-fw(i,jb,k))
          end do
        end do

      end do ! k


    else

      j=1

      do k = 1,nzm
        kc= k+1
        kcu =min(kc, nzm)
        irho = 1./(rhow(kc,icrm)*adzw(kc,icrm))

        do i = 0, nx
          ic = i+1
          fu(i,j,k)=dx25*(u(ic,j,k,icrm)+u(i,j,k,icrm))*(u(i,j,k,icrm)+u(ic,j,k,icrm))
          fv(i,j,k)=dx25*(u(ic,j,k,icrm)+u(i,j,k,icrm))*(v(i,j,k,icrm)+v(ic,j,k,icrm))
          fw(i,j,k)=dx25*(u(ic,j,k,icrm)*rho(k,icrm)*adz(k,icrm)+ &
          u(ic,j,kcu,icrm)*rho(kcu,icrm)*adz(kcu,icrm))*(w(i,j,kc,icrm)+w(ic,j,kc,icrm))
        end do
        do i = 1, nx
          ib = i-1
          dudt(i,j,k,na,icrm)  = dudt(i,j,k,na,icrm)  - (fu(i,j,k)-fu(ib,j,k))
          dvdt(i,j,k,na,icrm)  = dvdt(i,j,k,na,icrm)  - (fv(i,j,k)-fv(ib,j,k))
          dwdt(i,j,kc,na,icrm) = dwdt(i,j,kc,na,icrm)-irho*(fw(i,j,k)-fw(ib,j,k))
        end do

      end do ! k

    endif

  end subroutine advect2_mom_xy

end module advect2_mom_xy_mod
