module advect_scalar2D_mod
  implicit none

contains

  subroutine advect_scalar2D (f, u, w, rho, rhow, flux, ncrms)

    !     positively definite monotonic advection with non-oscillatory option

    use grid
    use params, only: dowallx, crm_rknd
    implicit none
    integer, intent(in) :: ncrms


    real(crm_rknd) f   (ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real(crm_rknd) u   (ncrms,dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm)
    real(crm_rknd) w   (ncrms,dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz )
    real(crm_rknd) rho (ncrms,nzm)
    real(crm_rknd) rhow(ncrms,nz)
    real(crm_rknd) flux(ncrms,nz)

    real(crm_rknd) mx (0:nxp1,1,nzm)
    real(crm_rknd) mn (0:nxp1,1,nzm)
    real(crm_rknd) uuu(-1:nxp3,1,nzm)
    real(crm_rknd) www(-1:nxp2,1,nz)

    real(crm_rknd) eps, dd
    integer i,j,k,ic,ib,kc,kb,icrm
    logical nonos
    real(crm_rknd) iadz(nzm),irho(nzm),irhow(nzm)

    real(crm_rknd) x1, x2, a, b, a1, a2, y
    real(crm_rknd) andiff,across,pp,pn
    andiff(x1,x2,a,b)=(abs(a)-a*a*b)*0.5*(x2-x1)
    across(x1,a1,a2)=0.03125*a1*a2*x1
    pp(y)= max(real(0.,crm_rknd),y)
    pn(y)=-min(real(0.,crm_rknd),y)

    do icrm = 1 , ncrms

      nonos = .true.
      eps = 1.e-10

      j=1

      www(:,:,nz)=0.

      if(dowallx) then

        if(mod(rank,nsubdomains_x).eq.0) then
          do k=1,nzm
            do i=dimx1_u,1
              u(icrm,i,j,k) = 0.
            end do
          end do
        end if
        if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
          do k=1,nzm
            do i=nx+1,dimx2_u
              u(icrm,i,j,k) = 0.
            end do
          end do
        end if

      end if

      !-----------------------------------------

      if(nonos) then

        do k=1,nzm
          kc=min(nzm,k+1)
          kb=max(1,k-1)
          do i=0,nxp1
            ib=i-1
            ic=i+1
            mx(i,j,k)=max(f(icrm,ib,j,k),f(icrm,ic,j,k),f(icrm,i,j,kb),f(icrm,i,j,kc),f(icrm,i,j,k))
            mn(i,j,k)=min(f(icrm,ib,j,k),f(icrm,ic,j,k),f(icrm,i,j,kb),f(icrm,i,j,kc),f(icrm,i,j,k))
          end do
        end do

      end if  ! nonos

      do k=1,nzm
        kb=max(1,k-1)
        do i=-1,nxp3
          uuu(i,j,k)=max(real(0.,crm_rknd),u(icrm,i,j,k))*f(icrm,i-1,j,k)+min(real(0.,crm_rknd),u(icrm,i,j,k))*f(icrm,i,j,k)
        end do
        do i=-1,nxp2
          www(i,j,k)=max(real(0.,crm_rknd),w(icrm,i,j,k))*f(icrm,i,j,kb)+min(real(0.,crm_rknd),w(icrm,i,j,k))*f(icrm,i,j,k)
        end do
        flux(icrm,k) = 0.
        do i=1,nx
          flux(icrm,k) = flux(icrm,k) + www(i,j,k)
        end do
      end do

      do k=1,nzm
        irho(k) = 1./rho(icrm,k)
        iadz(k) = 1./adz(icrm,k)
        do i=-1,nxp2
          f(icrm,i,j,k) = f(icrm,i,j,k) - (uuu(i+1,j,k)-uuu(i,j,k) &
          + (www(i,j,k+1)-www(i,j,k))*iadz(k))*irho(k)
        end do
      end do


      do k=1,nzm
        kc=min(nzm,k+1)
        kb=max(1,k-1)
        dd=2./(kc-kb)/adz(icrm,k)
        irhow(k)=1./(rhow(icrm,k)*adz(icrm,k))
        do i=0,nxp2
          ib=i-1
          uuu(i,j,k)=andiff(f(icrm,ib,j,k),f(icrm,i,j,k),u(icrm,i,j,k),irho(k)) &
          - across(dd*(f(icrm,ib,j,kc)+f(icrm,i,j,kc)-f(icrm,ib,j,kb)-f(icrm,i,j,kb)), &
          u(icrm,i,j,k), w(icrm,ib,j,k)+w(icrm,ib,j,kc)+w(icrm,i,j,k)+w(icrm,i,j,kc)) *irho(k)
        end do


        do i=0,nxp1
          ib=i-1
          ic=i+1
          www(i,j,k)=andiff(f(icrm,i,j,kb),f(icrm,i,j,k),w(icrm,i,j,k),irhow(k)) &
          -across(f(icrm,ic,j,kb)+f(icrm,ic,j,k)-f(icrm,ib,j,kb)-f(icrm,ib,j,k), &
          w(icrm,i,j,k), u(icrm,i,j,kb)+u(icrm,i,j,k)+u(icrm,ic,j,k)+u(icrm,ic,j,kb)) *irho(k)
        end do
      end do
      www(:,:,1) = 0.
      !---------- non-osscilatory option ---------------

      if(nonos) then

        do k=1,nzm
          kc=min(nzm,k+1)
          kb=max(1,k-1)
          do i=0,nxp1
            ib=i-1
            ic=i+1
            mx(i,j,k)=max(f(icrm,ib,j,k),f(icrm,ic,j,k),f(icrm,i,j,kb),f(icrm,i,j,kc),f(icrm,i,j,k),mx(i,j,k))
            mn(i,j,k)=min(f(icrm,ib,j,k),f(icrm,ic,j,k),f(icrm,i,j,kb),f(icrm,i,j,kc),f(icrm,i,j,k),mn(i,j,k))
          end do
        end do

        do k=1,nzm
          kc=min(nzm,k+1)
          do i=0,nxp1
            ic=i+1
            mx(i,j,k)=rho(icrm,k)*(mx(i,j,k)-f(icrm,i,j,k))/(pn(uuu(ic,j,k)) + pp(uuu(i,j,k))+&
            iadz(k)*(pn(www(i,j,kc)) + pp(www(i,j,k)))+eps)
            mn(i,j,k)=rho(icrm,k)*(f(icrm,i,j,k)-mn(i,j,k))/(pp(uuu(ic,j,k)) + pn(uuu(i,j,k))+&
            iadz(k)*(pp(www(i,j,kc)) + pn(www(i,j,k)))+eps)
          end do
        end do

        do k=1,nzm
          kb=max(1,k-1)
          do i=1,nxp1
            ib=i-1
            uuu(i,j,k)= pp(uuu(i,j,k))*min(real(1.,crm_rknd),mx(i,j,k), mn(ib,j,k)) &
            - pn(uuu(i,j,k))*min(real(1.,crm_rknd),mx(ib,j,k),mn(i,j,k))
          end do
          do i=1,nx
            www(i,j,k)= pp(www(i,j,k))*min(real(1.,crm_rknd),mx(i,j,k), mn(i,j,kb)) &
            - pn(www(i,j,k))*min(real(1.,crm_rknd),mx(i,j,kb),mn(i,j,k))
            flux(icrm,k) = flux(icrm,k) + www(i,j,k)
          end do
        end do


      endif ! nonos


      do k=1,nzm
        kc=k+1
        do i=1,nx
          ! MK: added fix for very small negative values (relative to positive values)
          !     especially  when such large numbers as
          !     hydrometeor concentrations are advected. The reason for negative values is
          !     most likely truncation error.
          f(icrm,i,j,k)= max(real(0.,crm_rknd), f(icrm,i,j,k) - (uuu(i+1,j,k)-uuu(i,j,k) &
          +(www(i,j,k+1)-www(i,j,k))*iadz(k))*irho(k))
        end do
      end do
    end do

  end subroutine advect_scalar2D

end module advect_scalar2D_mod
