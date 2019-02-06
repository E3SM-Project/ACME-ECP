module advect2_mom_xy_mod
  use params, only: asyncid
  implicit none

contains

  subroutine advect2_mom_xy(ncrms)

    !        momentum tendency due to 2nd-order-central horizontal advection

    use vars
    use params, only: crm_rknd

    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) dx25, dy25, irho, fu1, fu2, fv1, fv2, fw1, fw2
    integer i, j, k, kc, kcu, ic, jb, ib, jc,icrm

    dx25 = 0.25 / dx
    dy25 = 0.25 / dy

    if(RUN3D) then

      !$acc parallel loop collapse(4) copyin(rhow,adzw,u,v,w,adz,rho) copy(dudt,dvdt,dwdt) async(asyncid)
      do icrm = 1 , ncrms
        do k = 1,nzm
          do j = 1, ny
            do i = 1, nx
              kc= k+1
              kcu =min(kc, nzm)
              irho = 1./(rhow(kc,icrm)*adzw(kc,icrm))
              jb = j-1
              ic = i+1
              fu1 = dx25*(u(ic-1,j,k,icrm)+u(i-1,j,k,icrm))*(u(i-1,j,k,icrm)+u(ic-1,j,k,icrm))
              fu2 = dx25*(u(ic  ,j,k,icrm)+u(i  ,j,k,icrm))*(u(i  ,j,k,icrm)+u(ic  ,j,k,icrm))
              dudt(i,j,k,na,icrm)  = dudt(i,j,k,na,icrm)  - (fu2-fu1)
              fv1 = dx25*(u(ic-1,j,k,icrm)+u(ic-1,jb,k,icrm))*(v(i-1,j,k,icrm)+v(ic-1,j,k,icrm))
              fv2 = dx25*(u(ic  ,j,k,icrm)+u(ic  ,jb,k,icrm))*(v(i  ,j,k,icrm)+v(ic  ,j,k,icrm))
              dvdt(i,j,k,na,icrm)  = dvdt(i,j,k,na,icrm)  - (fv2-fv1)
              fw1 = dx25*(u(ic-1,j,k,icrm)*rho(k,icrm)*adz(k,icrm)+u(ic-1,j,kcu,icrm)*rho(kcu,icrm)*adz(kcu,icrm))*(w(i-1,j,kc,icrm)+w(ic-1,j,kc,icrm))
              fw2 = dx25*(u(ic  ,j,k,icrm)*rho(k,icrm)*adz(k,icrm)+u(ic  ,j,kcu,icrm)*rho(kcu,icrm)*adz(kcu,icrm))*(w(i  ,j,kc,icrm)+w(ic  ,j,kc,icrm))
              dwdt(i,j,kc,na,icrm) = dwdt(i,j,kc,na,icrm)-irho*(fw2-fw1)

              jc = j+1
              ib = i-1
              fu1 = dy25*(v(i,jc-1,k,icrm)+v(ib,jc-1,k,icrm))*(u(i,j-1,k,icrm)+u(i,jc-1,k,icrm))
              fu2 = dy25*(v(i,jc  ,k,icrm)+v(ib,jc  ,k,icrm))*(u(i,j  ,k,icrm)+u(i,jc  ,k,icrm))
              dudt(i,j,k,na,icrm) = dudt(i,j,k,na,icrm) - (fu2-fu1)
              fv1 = dy25*(v(i,jc-1,k,icrm)+v(i,j-1,k,icrm))*(v(i,j-1,k,icrm)+v(i,jc-1,k,icrm))
              fv2 = dy25*(v(i,jc  ,k,icrm)+v(i,j  ,k,icrm))*(v(i,j  ,k,icrm)+v(i,jc  ,k,icrm))
              dvdt(i,j,k,na,icrm) = dvdt(i,j,k,na,icrm) - (fv2-fv1)
              fw1 = dy25*(v(i,jc-1,k,icrm)*rho(k,icrm)*adz(k,icrm)+v(i,jc-1,kcu,icrm)*rho(kcu,icrm)*adz(kcu,icrm))*(w(i,j-1,kc,icrm)+w(i,jc-1,kc,icrm))
              fw2 = dy25*(v(i,jc  ,k,icrm)*rho(k,icrm)*adz(k,icrm)+v(i,jc  ,kcu,icrm)*rho(kcu,icrm)*adz(kcu,icrm))*(w(i,j  ,kc,icrm)+w(i,jc  ,kc,icrm))
              dwdt(i,j,kc,na,icrm)= dwdt(i,j,kc,na,icrm)-irho*(fw2-fw1)
            end do
          end do
        end do ! k
      end do

    else

      j=1
      !$acc parallel loop collapse(3) copy(dudt,dvdt,dwdt) copyin(w,v,u,adzw,rhow,rho,adz) async(asyncid)
      do icrm = 1 , ncrms
        do k = 1,nzm
          do i = 1, nx
            kc= k+1
            kcu =min(kc, nzm)
            irho = 1./(rhow(kc,icrm)*adzw(kc,icrm))
            ic = i+1
            fu1 = dx25*(u(ic-1,j,k,icrm)+u(i-1,j,k,icrm))*(u(i-1,j,k,icrm)+u(ic-1,j,k,icrm))
            fu2 = dx25*(u(ic  ,j,k,icrm)+u(i  ,j,k,icrm))*(u(i  ,j,k,icrm)+u(ic  ,j,k,icrm))
            dudt(i,j,k,na,icrm)  = dudt(i,j,k,na,icrm)  - (fu2-fu1)
            fv1 = dx25*(u(ic-1,j,k,icrm)+u(i-1,j,k,icrm))*(v(i-1,j,k,icrm)+v(ic-1,j,k,icrm))
            fv2 = dx25*(u(ic  ,j,k,icrm)+u(i  ,j,k,icrm))*(v(i  ,j,k,icrm)+v(ic  ,j,k,icrm))
            dvdt(i,j,k,na,icrm)  = dvdt(i,j,k,na,icrm)  - (fv2-fv1)
            fw1 = dx25*(u(ic-1,j,k,icrm)*rho(k,icrm)*adz(k,icrm)+u(ic-1,j,kcu,icrm)*rho(kcu,icrm)*adz(kcu,icrm))*(w(i-1,j,kc,icrm)+w(ic-1,j,kc,icrm))
            fw2 = dx25*(u(ic  ,j,k,icrm)*rho(k,icrm)*adz(k,icrm)+u(ic  ,j,kcu,icrm)*rho(kcu,icrm)*adz(kcu,icrm))*(w(i  ,j,kc,icrm)+w(ic  ,j,kc,icrm))
            dwdt(i,j,kc,na,icrm) = dwdt(i,j,kc,na,icrm)-irho*(fw2-fw1)
          end do
        end do ! k
      end do

    endif

  end subroutine advect2_mom_xy

end module advect2_mom_xy_mod
