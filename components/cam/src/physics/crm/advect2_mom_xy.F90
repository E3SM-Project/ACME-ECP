module advect2_mom_xy_mod
  implicit none

contains

  subroutine advect2_mom_xy(ncrms)

    !        momentum tendency due to 2nd-order-central horizontal advection

    use vars
    use params, only: crm_rknd

    implicit none
    integer, intent(in) :: ncrms

    real(crm_rknd) fu(ncrms,0:nx,1-YES3D:ny,nzm)
    real(crm_rknd) fv(ncrms,0:nx,1-YES3D:ny,nzm)
    real(crm_rknd) fw(ncrms,0:nx,1-YES3D:ny,nzm)
    real(crm_rknd) dx25, dy25, irho(ncrms), tmp1, tmp2

    integer i, j, k, kc, kcu, ic, jb, ib, jc, icrm

    dx25 = 0.25 / dx
    dy25 = 0.25 / dy


    if(RUN3D) then


      do k = 1,nzm
        kc= k+1
        kcu =min(kc, nzm)
        irho(:) = 1./(rhow(:,kc)*adzw(:,kc))

        do j = 1, ny
          jb = j-1
          do i = 0, nx
            ic = i+1
            fu(:,i,j,k)=dx25*(u(:,ic,j,k)+u(:,i,j,k))*(u(:,i,j,k)+u(:,ic,j,k))
            fv(:,i,j,k)=dx25*(u(:,ic,j,k)+u(:,ic,jb,k))*(v(:,i,j,k)+v(:,ic,j,k))
            fw(:,i,j,k)=dx25*(u(:,ic,j,k)*rho(:,k)*adz(:,k)+ &
            u(:,ic,j,kcu)*rho(:,kcu)*adz(:,kcu))*(w(:,i,j,kc)+w(:,ic,j,kc))
          end do
          do i = 1, nx
            ib = i-1
            dudt(:,i,j,k,na)  = dudt(:,i,j,k,na)  - (fu(:,i,j,k)-fu(:,ib,j,k))
            dvdt(:,i,j,k,na)  = dvdt(:,i,j,k,na)  - (fv(:,i,j,k)-fv(:,ib,j,k))
            dwdt(:,i,j,kc,na) = dwdt(:,i,j,kc,na)-irho(:)*(fw(:,i,j,k)-fw(:,ib,j,k))
          end do
        end do

        do j = 0, ny
          jc = j+1
          do i = 1, nx
            ib = i-1
            fu(:,i,j,k)=dy25*(v(:,i,jc,k)+v(:,ib,jc,k))*(u(:,i,j,k)+u(:,i,jc,k))
            fv(:,i,j,k)=dy25*(v(:,i,jc,k)+v(:,i,j,k))*(v(:,i,j,k)+v(:,i,jc,k))
            fw(:,i,j,k)=dy25*(v(:,i,jc,k)*rho(:,k)*adz(:,k)+ &
            v(:,i,jc,kcu)*rho(:,kcu)*adz(:,kcu))*(w(:,i,j,kc)+w(:,i,jc,kc))
          end do
        end do
        do j = 1,ny
          jb = j-1
          do i = 1, nx
            dudt(:,i,j,k,na) = dudt(:,i,j,k,na) - (fu(:,i,j,k) - fu(:,i,jb,k))
            dvdt(:,i,j,k,na) = dvdt(:,i,j,k,na) - (fv(:,i,j,k) - fv(:,i,jb,k))
            dwdt(:,i,j,kc,na)= dwdt(:,i,j,kc,na)-irho(:)*(fw(:,i,j,k)-fw(:,i,jb,k))
          end do
        end do

      end do ! k


    else

      j=1

      !$acc parallel loop gang vector collapse(3) default(present) async(1)
      do k = 1,nzm
        do i = 1, nx
          do icrm = 1 , ncrms
            kc= k+1
            kcu =min(kc, nzm)
            ib = i-1
            ic = i+1
            tmp1 = dx25*(u(icrm,ic,j,k)+u(icrm,i ,j,k))*(u(icrm,i ,j,k)+u(icrm,ic,j,k))
            tmp2 = dx25*(u(icrm,i ,j,k)+u(icrm,ib,j,k))*(u(icrm,ib,j,k)+u(icrm,i ,j,k))
            dudt(icrm,i,j,k,na)  = dudt(icrm,i,j,k,na)  - (tmp1-tmp2)

            tmp1 = dx25*(u(icrm,ic,j,k)+u(icrm,i ,j,k))*(v(icrm,i ,j,k)+v(icrm,ic,j,k))
            tmp2 = dx25*(u(icrm,i ,j,k)+u(icrm,ib,j,k))*(v(icrm,ib,j,k)+v(icrm,i ,j,k))
            dvdt(icrm,i,j,k,na)  = dvdt(icrm,i,j,k,na)  - (tmp1-tmp2)

            tmp1 = dx25*(u(icrm,ic,j,k)*rho(icrm,k)*adz(icrm,k)+u(icrm,ic,j,kcu)*rho(icrm,kcu)*adz(icrm,kcu))*(w(icrm,i ,j,kc)+w(icrm,ic,j,kc))
            tmp2 = dx25*(u(icrm,i ,j,k)*rho(icrm,k)*adz(icrm,k)+u(icrm,i ,j,kcu)*rho(icrm,kcu)*adz(icrm,kcu))*(w(icrm,ib,j,kc)+w(icrm,i ,j,kc))
            dwdt(icrm,i,j,kc,na) = dwdt(icrm,i,j,kc,na)-(tmp1-tmp2)/(rhow(icrm,kc)*adzw(icrm,kc))
          end do
        end do
      end do ! k

    endif

  end subroutine advect2_mom_xy

end module advect2_mom_xy_mod
