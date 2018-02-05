module advect_scalar3D_mod
  implicit none

contains

  subroutine advect_scalar3D (f, u, v, w, rho, rhow, flux, ncrms)

    !     positively definite monotonic advection with non-oscillatory option

    use grid
    use params, only: dowallx, dowally, crm_rknd
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) f   (ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real(crm_rknd) u   (ncrms,dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm)
    real(crm_rknd) v   (ncrms,dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm)
    real(crm_rknd) w   (ncrms,dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz )
    real(crm_rknd) rho (ncrms,nzm)
    real(crm_rknd) rhow(ncrms,nz)
    real(crm_rknd) flux(ncrms,nz)

    real(crm_rknd) mx (ncrms,0:nxp1,0:nyp1,nzm)
    real(crm_rknd) mn (ncrms,0:nxp1,0:nyp1,nzm)
    real(crm_rknd) uuu(ncrms,-1:nxp3,-1:nyp2,nzm)
    real(crm_rknd) vvv(ncrms,-1:nxp2,-1:nyp3,nzm)
    real(crm_rknd) www(ncrms,-1:nxp2,-1:nyp2,nz)

    real(crm_rknd) eps, dd(ncrms)
    real(crm_rknd) iadz(ncrms,nzm),irho(ncrms,nzm),irhow(ncrms,nzm)
    integer i,j,k,ic,ib,jc,jb,kc,kb,icrm
    logical nonos

    real(crm_rknd) x1, x2, a, b, a1, a2, y
    real(crm_rknd) andiff,across,pp,pn
    andiff(x1,x2,a,b)=(abs(a)-a*a*b)*0.5*(x2-x1)
    across(x1,a1,a2)=0.03125*a1*a2*x1
    pp(y)= max(real(0.,crm_rknd),y)
    pn(y)=-min(real(0.,crm_rknd),y)


    nonos = .true.
    eps = 1.e-10

    www(:,:,:,nz)=0.

    if(dowallx) then
      if(mod(rank,nsubdomains_x).eq.0) then
        do k=1,nzm
          do j=dimy1_u,dimy2_u
            do i=dimx1_u,1
              do icrm = 1 , ncrms
                u(icrm,i,j,k) = 0.
              end do
            end do
          end do
        end do
      end if
      if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
        do k=1,nzm
          do j=dimy1_u,dimy2_u
            do i=nx+1,dimx2_u
              do icrm = 1 , ncrms
                u(icrm,i,j,k) = 0.
              end do
            end do
          end do
        end do
      end if

    end if

    if(dowally) then

      if(rank.lt.nsubdomains_x) then
        do k=1,nzm
          do j=dimy1_v,1
            do i=dimx1_v,dimx2_v
              do icrm = 1 , ncrms
                v(icrm,i,j,k) = 0.
              end do
            end do
          end do
        end do
      end if
      if(rank.gt.nsubdomains-nsubdomains_x-1) then
        do k=1,nzm
          do j=ny+1,dimy2_v
            do i=dimx1_v,dimx2_v
              do icrm = 1 , ncrms
                v(icrm,i,j,k) = 0.
              end do
            end do
          end do
        end do
      end if

    end if

    !-----------------------------------------

    if(nonos) then

      do k=1,nzm
        kc=min(nzm,k+1)
        kb=max(1,k-1)
        do j=0,nyp1
          jb=j-1
          jc=j+1
          do i=0,nxp1
            do icrm = 1 , ncrms
              ib=i-1
              ic=i+1
              mx(icrm,i,j,k)=max(f(icrm,ib,j,k),f(icrm,ic,j,k),f(icrm,i,jb,k), &
              f(icrm,i,jc,k),f(icrm,i,j,kb),f(icrm,i,j,kc),f(icrm,i,j,k))
              mn(icrm,i,j,k)=min(f(icrm,ib,j,k),f(icrm,ic,j,k),f(icrm,i,jb,k), &
              f(icrm,i,jc,k),f(icrm,i,j,kb),f(icrm,i,j,kc),f(icrm,i,j,k))
            end do
          end do
        end do
      end do

    end if  ! nonos

    do k=1,nzm
      do j=-1,nyp2
        do i=-1,nxp3
          do icrm = 1 , ncrms
            uuu(icrm,i,j,k)=max(real(0.,crm_rknd),u(icrm,i,j,k))*f(icrm,i-1,j,k)+min(real(0.,crm_rknd),u(icrm,i,j,k))*f(icrm,i,j,k)
          end do
        end do
      end do
    end do

    do k=1,nzm
      do j=-1,nyp3
        do i=-1,nxp2
          do icrm = 1 , ncrms
            vvv(icrm,i,j,k)=max(real(0.,crm_rknd),v(icrm,i,j,k))*f(icrm,i,j-1,k)+min(real(0.,crm_rknd),v(icrm,i,j,k))*f(icrm,i,j,k)
          end do
        end do
      end do
    end do

    do k=1,nzm
      kb=max(1,k-1)
      do j=-1,nyp2
        do i=-1,nxp2
          do icrm = 1 , ncrms
            www(icrm,i,j,k)=max(real(0.,crm_rknd),w(icrm,i,j,k))*f(icrm,i,j,kb)+min(real(0.,crm_rknd),w(icrm,i,j,k))*f(icrm,i,j,k)
          end do
        end do
      end do
      flux(:,k) = 0.
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            flux(icrm,k) = flux(icrm,k) + www(icrm,i,j,k)
          end do
        end do
      end do
    end do


    do k=1,nzm
      irho(:,k) = 1./rho(:,k)
      iadz(:,k) = 1./adz(:,k)
      do j=-1,nyp2
        do i=-1,nxp2
          do icrm = 1 , ncrms
            f(icrm,i,j,k)=f(icrm,i,j,k) -(uuu(icrm,i+1,j,k)-uuu(icrm,i,j,k)+vvv(icrm,i,j+1,k)-vvv(icrm,i,j,k) &
            +(www(icrm,i,j,k+1)-www(icrm,i,j,k))*iadz(icrm,k))*irho(icrm,k)
          end do
        end do
      end do
    end do


    do k=1,nzm
      kc=min(nzm,k+1)
      kb=max(1,k-1)
      dd(:)=2./(kc-kb)/adz(:,k)
      do j=0,nyp1
        jb=j-1
        jc=j+1
        do i=0,nxp2
          do icrm = 1 , ncrms
            ib=i-1
            uuu(icrm,i,j,k)=andiff(f(icrm,ib,j,k),f(icrm,i,j,k),u(icrm,i,j,k),irho(icrm,k)) &
            -(across(f(icrm,ib,jc,k)+f(icrm,i,jc,k)-f(icrm,ib,jb,k)-f(icrm,i,jb,k), &
            u(icrm,i,j,k), v(icrm,ib,j,k)+v(icrm,ib,jc,k)+v(icrm,i,jc,k)+v(icrm,i,j,k)) &
            +across(dd(icrm)*(f(icrm,ib,j,kc)+f(icrm,i,j,kc)-f(icrm,ib,j,kb)-f(icrm,i,j,kb)), &
            u(icrm,i,j,k), w(icrm,ib,j,k)+w(icrm,ib,j,kc)+w(icrm,i,j,k)+w(icrm,i,j,kc))) *irho(icrm,k)
          end do
        end do
      end do
    end do

    do k=1,nzm
      kc=min(nzm,k+1)
      kb=max(1,k-1)
      dd(:)=2./(kc-kb)/adz(:,k)
      do j=0,nyp2
        jb=j-1
        do i=0,nxp1
          do icrm = 1 , ncrms
            ib=i-1
            ic=i+1
            vvv(icrm,i,j,k)=andiff(f(icrm,i,jb,k),f(icrm,i,j,k),v(icrm,i,j,k),irho(icrm,k)) &
            -(across(f(icrm,ic,jb,k)+f(icrm,ic,j,k)-f(icrm,ib,jb,k)-f(icrm,ib,j,k), &
            v(icrm,i,j,k), u(icrm,i,jb,k)+u(icrm,i,j,k)+u(icrm,ic,j,k)+u(icrm,ic,jb,k)) &
            +across(dd(icrm)*(f(icrm,i,jb,kc)+f(icrm,i,j,kc)-f(icrm,i,jb,kb)-f(icrm,i,j,kb)), &
            v(icrm,i,j,k), w(icrm,i,jb,k)+w(icrm,i,j,k)+w(icrm,i,j,kc)+w(icrm,i,jb,kc))) *irho(icrm,k)
          end do
        end do
      end do
    end do

    do k=1,nzm
      kb=max(1,k-1)
      irhow(:,k)=1./(rhow(:,k)*adz(:,k))
      do j=0,nyp1
        jb=j-1
        jc=j+1
        do i=0,nxp1
          do icrm = 1 , ncrms
            ib=i-1
            ic=i+1
            www(icrm,i,j,k)=andiff(f(icrm,i,j,kb),f(icrm,i,j,k),w(icrm,i,j,k),irhow(icrm,k)) &
            -(across(f(icrm,ic,j,kb)+f(icrm,ic,j,k)-f(icrm,ib,j,kb)-f(icrm,ib,j,k), &
            w(icrm,i,j,k), u(icrm,i,j,kb)+u(icrm,i,j,k)+u(icrm,ic,j,k)+u(icrm,ic,j,kb)) &
            +across(f(icrm,i,jc,k)+f(icrm,i,jc,kb)-f(icrm,i,jb,k)-f(icrm,i,jb,kb), &
            w(icrm,i,j,k), v(icrm,i,j,kb)+v(icrm,i,jc,kb)+v(icrm,i,jc,k)+v(icrm,i,j,k))) *irho(icrm,k)
          end do
        end do
      end do
    end do

    www(:,:,:,1) = 0.

    !---------- non-osscilatory option ---------------

    if(nonos) then

      do k=1,nzm
        kc=min(nzm,k+1)
        kb=max(1,k-1)
        do j=0,nyp1
          jb=j-1
          jc=j+1
          do i=0,nxp1
            do icrm = 1 , ncrms
              ib=i-1
              ic=i+1
              mx(icrm,i,j,k)=max(f(icrm,ib,j,k),f(icrm,ic,j,k),f(icrm,i,jb,k), &
              f(icrm,i,jc,k),f(icrm,i,j,kb),f(icrm,i,j,kc),f(icrm,i,j,k),mx(icrm,i,j,k))
              mn(icrm,i,j,k)=min(f(icrm,ib,j,k),f(icrm,ic,j,k),f(icrm,i,jb,k), &
              f(icrm,i,jc,k),f(icrm,i,j,kb),f(icrm,i,j,kc),f(icrm,i,j,k),mn(icrm,i,j,k))
            end do
          end do
        end do
      end do

      do k=1,nzm
        kc=min(nzm,k+1)
        do j=0,nyp1
          jc=j+1
          do i=0,nxp1
            do icrm = 1 , ncrms
              ic=i+1
              mx(icrm,i,j,k)=rho(icrm,k)*(mx(icrm,i,j,k)-f(icrm,i,j,k))/ &
              (pn(uuu(icrm,ic,j,k)) + pp(uuu(icrm,i,j,k))+ &
              pn(vvv(icrm,i,jc,k)) + pp(vvv(icrm,i,j,k))+ &
              iadz(icrm,k)*(pn(www(icrm,i,j,kc)) + pp(www(icrm,i,j,k)))+eps)
              mn(icrm,i,j,k)=rho(icrm,k)*(f(icrm,i,j,k)-mn(icrm,i,j,k))/ &
              (pp(uuu(icrm,ic,j,k)) + pn(uuu(icrm,i,j,k))+ &
              pp(vvv(icrm,i,jc,k)) + pn(vvv(icrm,i,j,k))+ &
              iadz(icrm,k)*(pp(www(icrm,i,j,kc)) + pn(www(icrm,i,j,k)))+eps)
            end do
          end do
        end do
      end do

      do k=1,nzm
        do j=1,ny
          do i=1,nxp1
            do icrm = 1 , ncrms
              ib=i-1
              uuu(icrm,i,j,k)=pp(uuu(icrm,i,j,k))*min(real(1.,crm_rknd),mx(icrm,i,j,k), mn(icrm,ib,j,k)) &
              - pn(uuu(icrm,i,j,k))*min(real(1.,crm_rknd),mx(icrm,ib,j,k),mn(icrm,i,j,k))
            end do
          end do
        end do
      end do

      do k=1,nzm
        do j=1,nyp1
          jb=j-1
          do i=1,nx
            do icrm = 1 , ncrms
              vvv(icrm,i,j,k)=pp(vvv(icrm,i,j,k))*min(real(1.,crm_rknd),mx(icrm,i,j,k), mn(icrm,i,jb,k)) &
              - pn(vvv(icrm,i,j,k))*min(real(1.,crm_rknd),mx(icrm,i,jb,k),mn(icrm,i,j,k))
            end do
          end do
        end do
      end do

      do k=1,nzm
        kb=max(1,k-1)
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              www(icrm,i,j,k)=pp(www(icrm,i,j,k))*min(real(1.,crm_rknd),mx(icrm,i,j,k), mn(icrm,i,j,kb)) &
              - pn(www(icrm,i,j,k))*min(real(1.,crm_rknd),mx(icrm,i,j,kb),mn(icrm,i,j,k))
              flux(icrm,k) = flux(icrm,k) + www(icrm,i,j,k)
            end do
          end do
        end do
      end do


    endif ! nonos


    do k=1,nzm
      kc=k+1
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            ! MK: added fix for very small negative values (relative to positive values)
            !     especially  when such large numbers as
            !     hydrometeor concentrations are advected. The reason for negative values is
            !     most likely truncation error.

            f(icrm,i,j,k)=max(real(0.,crm_rknd),f(icrm,i,j,k) -(uuu(icrm,i+1,j,k)-uuu(icrm,i,j,k)+vvv(icrm,i,j+1,k)-vvv(icrm,i,j,k) &
            +(www(icrm,i,j,k+1)-www(icrm,i,j,k))*iadz(icrm,k))*irho(icrm,k))
          end do
        end do
      end do
    end do

  end subroutine advect_scalar3D

end module advect_scalar3D_mod
