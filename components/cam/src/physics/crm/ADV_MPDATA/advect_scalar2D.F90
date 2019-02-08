module advect_scalar2D_mod
  use params, only: asyncid
  implicit none

contains

  subroutine advect_scalar2D (ncrms, f, u, w, rho, rhow, flux)
    !     positively definite monotonic advection with non-oscillatory option
    use grid
    use params, only: dowallx, crm_rknd
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) f(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real(crm_rknd) u(ncrms,dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm)
    real(crm_rknd) w(ncrms,dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz )
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

    !$acc enter data create(mx,mn,uuu,www,iadz,irho,irhow) async(asyncid)

    !$acc parallel loop collapse(2) copy(www) async(asyncid)
    do icrm = 1 , ncrms
      do i = -1 , nxp2
        www(i,j,nz,icrm)=0.
      enddo
    enddo

    if (dowallx) then
      if (mod(rank,nsubdomains_x).eq.0) then
        !$acc parallel loop collapse(3) copy(u) async(asyncid)
        do icrm = 1 , ncrms
          do k=1,nzm
            do i=dimx1_u,1
              u(icrm,i,j,k) = 0.
            enddo
          enddo
        enddo
      endif
      if (mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
        !$acc parallel loop collapse(3) copy(u) async(asyncid)
        do icrm = 1 , ncrms
          do k=1,nzm
            do i=nx+1,dimx2_u
              u(icrm,i,j,k) = 0.
            enddo
          enddo
        enddo
      endif
    endif

    !-----------------------------------------

    if (nonos) then
      !$acc parallel loop collapse(3) copyin(f) copy(mx,mn) async(asyncid)
      do icrm = 1 , ncrms
        do k=1,nzm
          do i=0,nxp1
            kc=min(nzm,k+1)
            kb=max(1,k-1)
            ib=i-1
            ic=i+1
            mx(i,j,k,icrm)=max(f(icrm,ib,j,k),f(icrm,ic,j,k),f(icrm,i,j,kb),f(icrm,i,j,kc),f(icrm,i,j,k))
            mn(i,j,k,icrm)=min(f(icrm,ib,j,k),f(icrm,ic,j,k),f(icrm,i,j,kb),f(icrm,i,j,kc),f(icrm,i,j,k))
          enddo
        enddo
      enddo
    endif  ! nonos

    !$acc parallel loop collapse(3) copyin(u,f,w) copy(uuu,www,flux) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        do i=-1,nxp3
          kb=max(1,k-1)
          uuu(i,j,k,icrm)=max(real(0.,crm_rknd),u(icrm,i,j,k))*f(icrm,i-1,j,k)+&
                          min(real(0.,crm_rknd),u(icrm,i,j,k))*f(icrm,i,j,k)
          if (i <= nxp2) www(i,j,k,icrm)=max(real(0.,crm_rknd),w(icrm,i,j,k))*&
                                          f(icrm,i,j,kb)+min(real(0.,crm_rknd),w(icrm,i,j,k))*f(icrm,i,j,k)
          if (i == 0) flux(k,icrm) = 0.
        enddo
      enddo
    enddo
    !$acc parallel loop collapse(2) copyin(rho,adz,rhow) copy(irho,iadz,irhow) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        irho(k,icrm) = 1./rho(k,icrm)
        iadz(k,icrm) = 1./adz(icrm,k)
        irhow(k,icrm)=1./(rhow(k,icrm)*adz(icrm,k))
      enddo
    enddo
    !$acc parallel loop collapse(3) copyin(uuu,www,iadz,irho) copy(f,flux) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        do i=-1,nxp2
          if (i >= 1 .and. i <= nx) then
            !$acc atomic update
            flux(k,icrm) = flux(k,icrm) + www(i,j,k,icrm)
          endif
          f(icrm,i,j,k) = f(icrm,i,j,k) - (uuu(i+1,j,k,icrm)-uuu(i,j,k,icrm)  + &
                                          (www(i,j,k+1,icrm)-www(i,j,k,icrm))*iadz(k,icrm))*irho(k,icrm)
        enddo
      enddo
    enddo

    !$acc parallel loop collapse(3) copyin(adz,f,u,irho,w,irhow) copy(uuu,www) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        do i=0,nxp2
          kc=min(nzm,k+1)
          kb=max(1,k-1)
          dd=2./(kc-kb)/adz(icrm,k)
          ib=i-1
          uuu(i,j,k,icrm)=andiff(f(icrm,ib,j,k),f(icrm,i,j,k),u(icrm,i,j,k),irho(k,icrm)) &
          - across(dd*(f(icrm,ib,j,kc)+f(icrm,i,j,kc)-f(icrm,ib,j,kb)-f(icrm,i,j,kb)), &
          u(icrm,i,j,k), w(icrm,ib,j,k)+w(icrm,ib,j,kc)+w(icrm,i,j,k)+w(icrm,i,j,kc)) *irho(k,icrm)
          if (i <= nxp1) then
            ic=i+1
            www(i,j,k,icrm)=andiff(f(icrm,i,j,kb),f(icrm,i,j,k),w(icrm,i,j,k),irhow(k,icrm)) &
            -across(f(icrm,ic,j,kb)+f(icrm,ic,j,k)-f(icrm,ib,j,kb)-f(icrm,ib,j,k), &
            w(icrm,i,j,k), u(icrm,i,j,kb)+u(icrm,i,j,k)+u(icrm,ic,j,k)+u(icrm,ic,j,kb)) *irho(k,icrm)
          endif
        enddo
      enddo
    enddo

    !$acc parallel loop collapse(2) copy(www) async(asyncid)
    do icrm = 1 , ncrms
      do i = -1 , nxp2
        www(i,j,1,icrm) = 0.
      enddo
    enddo
    !---------- non-osscilatory option ---------------

    if (nonos) then
      !$acc parallel loop collapse(3) copyin(f) copy(mx,mn) async(asyncid)
      do icrm = 1 , ncrms
        do k=1,nzm
          do i=0,nxp1
            kc=min(nzm,k+1)
            kb=max(1,k-1)
            ib=i-1
            ic=i+1
            mx(i,j,k,icrm)=max(f(icrm,ib,j,k),f(icrm,ic,j,k),f(icrm,i,j,kb),f(icrm,i,j,kc),f(icrm,i,j,k),mx(i,j,k,icrm))
            mn(i,j,k,icrm)=min(f(icrm,ib,j,k),f(icrm,ic,j,k),f(icrm,i,j,kb),f(icrm,i,j,kc),f(icrm,i,j,k),mn(i,j,k,icrm))
          enddo
        enddo
      enddo

      !$acc parallel loop collapse(3) copyin(f,rho,uuu,www,iadz) copy(mx,mn) async(asyncid)
      do icrm = 1 , ncrms
        do k=1,nzm
          do i=0,nxp1
            kc=min(nzm,k+1)
            ic=i+1
            mx(i,j,k,icrm)=rho(k,icrm)*(mx(i,j,k,icrm)-f(icrm,i,j,k))/(pn(uuu(ic,j,k,icrm)) + &
                           pp(uuu(i,j,k,icrm))+iadz(k,icrm)*(pn(www(i,j,kc,icrm)) + pp(www(i,j,k,icrm)))+eps)
            mn(i,j,k,icrm)=rho(k,icrm)*(f(icrm,i,j,k)-mn(i,j,k,icrm))/(pp(uuu(ic,j,k,icrm)) + &
                           pn(uuu(i,j,k,icrm))+iadz(k,icrm)*(pp(www(i,j,kc,icrm)) + pn(www(i,j,k,icrm)))+eps)
          enddo
        enddo
      enddo

      !$acc parallel loop collapse(3) copyin(mx,mn) copy(uuu,www,flux) async(asyncid)
      do icrm = 1 , ncrms
        do k=1,nzm
          do i=1,nxp1
            ib=i-1
            uuu(i,j,k,icrm)= pp(uuu(i,j,k,icrm))*min(real(1.,crm_rknd),mx(i,j,k,icrm), mn(ib,j,k,icrm)) - &
                             pn(uuu(i,j,k,icrm))*min(real(1.,crm_rknd),mx(ib,j,k,icrm),mn(i,j,k,icrm))
            if (i <= nx) then
              kb=max(1,k-1)
              www(i,j,k,icrm)= pp(www(i,j,k,icrm))*min(real(1.,crm_rknd),mx(i,j,k,icrm), mn(i,j,kb,icrm)) - &
                               pn(www(i,j,k,icrm))*min(real(1.,crm_rknd),mx(i,j,kb,icrm),mn(i,j,k,icrm))
              !$acc atomic update
              flux(k,icrm) = flux(k,icrm) + www(i,j,k,icrm)
            endif
          enddo
        enddo
      enddo
    endif ! nonos

    !$acc parallel loop collapse(3) copyin(uuu,www,iadz,irho) copy(f) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        do i=1,nx
          kc=k+1
          ! MK: added fix for very small negative values (relative to positive values)
          !     especially  when such large numbers as
          !     hydrometeor concentrations are advected. The reason for negative values is
          !     most likely truncation error.
          f(icrm,i,j,k)= max(real(0.,crm_rknd), f(icrm,i,j,k) - (uuu(i+1,j,k,icrm)-uuu(i,j,k,icrm) + &
                         (www(i,j,k+1,icrm)-www(i,j,k,icrm))*iadz(k,icrm))*irho(k,icrm))
        enddo
      enddo
    enddo

    !$acc exit data delete(mx,mn,uuu,www,iadz,irho,irhow) async(asyncid)

  end subroutine advect_scalar2D

end module advect_scalar2D_mod
