module advect_scalar2D_mod
  implicit none

contains

  subroutine advect_scalar2D (ncrms, f, u, w, rho, rhow, flux)
    !     positively definite monotonic advection with non-oscillatory option
    use grid
    use params, only: dowallx, crm_rknd
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) f(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm,ncrms)
    real(crm_rknd) u(dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm,ncrms)
    real(crm_rknd) w(dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz ,ncrms)
    real(crm_rknd) rho(nzm,ncrms)
    real(crm_rknd) rhow(nz,ncrms)
    real(crm_rknd) flux(nz,ncrms)

    real(crm_rknd) mx (0:nxp1 ,1,nzm,ncrms)
    real(crm_rknd) mn (0:nxp1 ,1,nzm,ncrms)
    real(crm_rknd) uuu(-1:nxp3,1,nzm,ncrms)
    real(crm_rknd) www(-1:nxp2,1,nz ,ncrms)
    real(crm_rknd) eps, dd
    integer i,j,k,ic,ib,kc,kb,icrm
    logical nonos
    real(crm_rknd) iadz(nzm,ncrms),irho(nzm,ncrms),irhow(nzm,ncrms)
    real(crm_rknd) x1, x2, a, b, a1, a2, y
    real(crm_rknd) andiff,across,pp,pn

    !Statement functions
    andiff(x1,x2,a,b)=(abs(a)-a*a*b)*0.5*(x2-x1)
    across(x1,a1,a2)=0.03125*a1*a2*x1
    pp(y)= max(real(0.,crm_rknd),y)
    pn(y)=-min(real(0.,crm_rknd),y)

    nonos = .true.
    eps = 1.e-10

    j=1

    !$acc enter data copyin(f,u,w,rho,rhow,flux) async(1)

    !$acc enter data create(mx,mn,uuu,www,iadz,irho,irhow) async(1)

    !$acc parallel loop collapse(2) copy(www) async(1)
    do icrm = 1 , ncrms
      do i = -1 , nxp2
        www(i,j,nz,icrm)=0.
      enddo
    enddo

    if (dowallx) then
      if (mod(rank,nsubdomains_x).eq.0) then
        !$acc parallel loop collapse(3) copy(u) async(1)
        do icrm = 1 , ncrms
          do k=1,nzm
            do i=dimx1_u,1
              u(i,j,k,icrm) = 0.
            enddo
          enddo
        enddo
      endif
      if (mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
        !$acc parallel loop collapse(3) copy(u) async(1)
        do icrm = 1 , ncrms
          do k=1,nzm
            do i=nx+1,dimx2_u
              u(i,j,k,icrm) = 0.
            enddo
          enddo
        enddo
      endif
    endif

    !-----------------------------------------

    if (nonos) then
      !$acc parallel loop collapse(3) copyin(f) copy(mx,mn) async(1)
      do icrm = 1 , ncrms
        do k=1,nzm
          do i=0,nxp1
            kc=min(nzm,k+1)
            kb=max(1,k-1)
            ib=i-1
            ic=i+1
            mx(i,j,k,icrm)=max(f(ib,j,k,icrm),f(ic,j,k,icrm),f(i,j,kb,icrm),f(i,j,kc,icrm),f(i,j,k,icrm))
            mn(i,j,k,icrm)=min(f(ib,j,k,icrm),f(ic,j,k,icrm),f(i,j,kb,icrm),f(i,j,kc,icrm),f(i,j,k,icrm))
          enddo
        enddo
      enddo
    endif  ! nonos

    !$acc parallel loop collapse(3) copyin(u,f,w) copy(uuu,www,flux) async(1)
    do icrm = 1 , ncrms
      do k=1,nzm
        do i=-1,nxp3
          kb=max(1,k-1)
          uuu(i,j,k,icrm)=max(real(0.,crm_rknd),u(i,j,k,icrm))*f(i-1,j,k,icrm)+min(real(0.,crm_rknd),u(i,j,k,icrm))*f(i,j,k,icrm)
          if (i <= nxp2) www(i,j,k,icrm)=max(real(0.,crm_rknd),w(i,j,k,icrm))*f(i,j,kb,icrm)+min(real(0.,crm_rknd),w(i,j,k,icrm))*f(i,j,k,icrm)
          if (i == 0) flux(k,icrm) = 0.
        enddo
      enddo
    enddo

    !$acc parallel loop collapse(3) copyin(rho,adz,rhow,uuu,www) copy(irho,iadz,irhow,f,flux) async(1)
    do icrm = 1 , ncrms
      do k=1,nzm
        do i=-1,nxp2
          if (i >= 1 .and. i <= nx) then
            !$acc atomic update
            flux(k,icrm) = flux(k,icrm) + www(i,j,k,icrm)
          endif
          if (i == -1) then
            irho(k,icrm) = 1./rho(k,icrm)
            iadz(k,icrm) = 1./adz(k,icrm)
            irhow(k,icrm)=1./(rhow(k,icrm)*adz(k,icrm))
          endif
          f(i,j,k,icrm) = f(i,j,k,icrm) - (uuu(i+1,j,k,icrm)-uuu(i,j,k,icrm)  + (www(i,j,k+1,icrm)-www(i,j,k,icrm))*iadz(k,icrm))*irho(k,icrm)
        enddo
      enddo
    enddo

    !$acc exit data copyout(mx,mn,uuu,www,iadz,irho,irhow) async(1)

    !$acc exit data copyout(f,u,w,rho,rhow,flux) async(1)
    !$acc wait(1)

    !!$acc parallel loop collapse(3) copyin(adz,f,u,irho,w,irhow) copy(uuu,www) async(1)
    do icrm = 1 , ncrms
      do k=1,nzm
        do i=0,nxp2
          kc=min(nzm,k+1)
          kb=max(1,k-1)
          dd=2./(kc-kb)/adz(k,icrm)
          ib=i-1
          uuu(i,j,k,icrm)=andiff(f(ib,j,k,icrm),f(i,j,k,icrm),u(i,j,k,icrm),irho(k,icrm)) &
          - across(dd*(f(ib,j,kc,icrm)+f(i,j,kc,icrm)-f(ib,j,kb,icrm)-f(i,j,kb,icrm)), &
          u(i,j,k,icrm), w(ib,j,k,icrm)+w(ib,j,kc,icrm)+w(i,j,k,icrm)+w(i,j,kc,icrm)) *irho(k,icrm)
          if (i <= nxp1) then
            ic=i+1
            www(i,j,k,icrm)=andiff(f(i,j,kb,icrm),f(i,j,k,icrm),w(i,j,k,icrm),irhow(k,icrm)) &
            -across(f(ic,j,kb,icrm)+f(ic,j,k,icrm)-f(ib,j,kb,icrm)-f(ib,j,k,icrm), &
            w(i,j,k,icrm), u(i,j,kb,icrm)+u(i,j,k,icrm)+u(ic,j,k,icrm)+u(ic,j,kb,icrm)) *irho(k,icrm)
          endif
        enddo
      enddo
    enddo
    !!$acc parallel loop collapse(2) copy(www) async(1)
    do icrm = 1 , ncrms
      do i = -1 , nxp2
        www(i,j,1,icrm) = 0.
      enddo
    enddo
    !---------- non-osscilatory option ---------------

    if (nonos) then
      !!$acc parallel loop collapse(3) copyin(f) copy(mx,mn) async(1)
      do icrm = 1 , ncrms
        do k=1,nzm
          do i=0,nxp1
            kc=min(nzm,k+1)
            kb=max(1,k-1)
            ib=i-1
            ic=i+1
            mx(i,j,k,icrm)=max(f(ib,j,k,icrm),f(ic,j,k,icrm),f(i,j,kb,icrm),f(i,j,kc,icrm),f(i,j,k,icrm),mx(i,j,k,icrm))
            mn(i,j,k,icrm)=min(f(ib,j,k,icrm),f(ic,j,k,icrm),f(i,j,kb,icrm),f(i,j,kc,icrm),f(i,j,k,icrm),mn(i,j,k,icrm))
          enddo
        enddo
      enddo

      !!$acc parallel loop collapse(3) copyin(f,rho,uuu,www,iadz) copy(mx,mn) async(1)
      do icrm = 1 , ncrms
        do k=1,nzm
          do i=0,nxp1
            kc=min(nzm,k+1)
            ic=i+1
            mx(i,j,k,icrm)=rho(k,icrm)*(mx(i,j,k,icrm)-f(i,j,k,icrm))/(pn(uuu(ic,j,k,icrm)) + pp(uuu(i,j,k,icrm))+iadz(k,icrm)*(pn(www(i,j,kc,icrm)) + pp(www(i,j,k,icrm)))+eps)
            mn(i,j,k,icrm)=rho(k,icrm)*(f(i,j,k,icrm)-mn(i,j,k,icrm))/(pp(uuu(ic,j,k,icrm)) + pn(uuu(i,j,k,icrm))+iadz(k,icrm)*(pp(www(i,j,kc,icrm)) + pn(www(i,j,k,icrm)))+eps)
          enddo
        enddo
      enddo

      !!$acc parallel loop collapse(3) copyin(mx,mn) copy(uuu,www,flux) async(1)
      do icrm = 1 , ncrms
        do k=1,nzm
          do i=1,nxp1
            ib=i-1
            uuu(i,j,k,icrm)= pp(uuu(i,j,k,icrm))*min(real(1.,crm_rknd),mx(i,j,k,icrm), mn(ib,j,k,icrm)) - pn(uuu(i,j,k,icrm))*min(real(1.,crm_rknd),mx(ib,j,k,icrm),mn(i,j,k,icrm))
            if (i <= nx) then
              kb=max(1,k-1)
              www(i,j,k,icrm)= pp(www(i,j,k,icrm))*min(real(1.,crm_rknd),mx(i,j,k,icrm), mn(i,j,kb,icrm)) - pn(www(i,j,k,icrm))*min(real(1.,crm_rknd),mx(i,j,kb,icrm),mn(i,j,k,icrm))
              !$acc atomic update
              flux(k,icrm) = flux(k,icrm) + www(i,j,k,icrm)
            endif
          enddo
        enddo
      enddo
    endif ! nonos

    !!$acc parallel loop collapse(3) copyin(uuu,www,iadz,irho) copy(f) async(1)
    do icrm = 1 , ncrms
      do k=1,nzm
        do i=1,nx
          kc=k+1
          ! MK: added fix for very small negative values (relative to positive values)
          !     especially  when such large numbers as
          !     hydrometeor concentrations are advected. The reason for negative values is
          !     most likely truncation error.
          f(i,j,k,icrm)= max(real(0.,crm_rknd), f(i,j,k,icrm) - (uuu(i+1,j,k,icrm)-uuu(i,j,k,icrm) + (www(i,j,k+1,icrm)-www(i,j,k,icrm))*iadz(k,icrm))*irho(k,icrm))
        enddo
      enddo
    enddo

  end subroutine advect_scalar2D

end module advect_scalar2D_mod
