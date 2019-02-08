module advect_scalar3D_mod
  use params, only: asyncid
  implicit none

contains

  subroutine advect_scalar3D (ncrms, f, u, v, w, rho, rhow, flux)
    !     positively definite monotonic advection with non-oscillatory option
    use grid
    use params, only: dowallx, dowally, crm_rknd
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) f(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real(crm_rknd) u(ncrms,dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm)
    real(crm_rknd) v(ncrms,dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm)
    real(crm_rknd) w(ncrms,dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz )
    real(crm_rknd) rho(ncrms,nzm)
    real(crm_rknd) rhow(nz,ncrms)
    real(crm_rknd) flux(nz,ncrms)
    real(crm_rknd) mx (0:nxp1 ,0:nyp1 ,nzm,ncrms)
    real(crm_rknd) mn (0:nxp1 ,0:nyp1 ,nzm,ncrms)
    real(crm_rknd) uuu(-1:nxp3,-1:nyp2,nzm,ncrms)
    real(crm_rknd) vvv(-1:nxp2,-1:nyp3,nzm,ncrms)
    real(crm_rknd) www(-1:nxp2,-1:nyp2,nz ,ncrms)
    real(crm_rknd) eps, dd
    real(crm_rknd) iadz(nzm,ncrms),irho(nzm,ncrms),irhow(nzm,ncrms)
    integer i,j,k,ic,ib,jc,jb,kc,kb, icrm
    logical nonos
    real(crm_rknd) x1, x2, a, b, a1, a2, y
    real(crm_rknd) andiff,across,pp,pn

    !Statement functions
    andiff(x1,x2,a,b)=(abs(a)-a*a*b)*0.5*(x2-x1)
    across(x1,a1,a2)=0.03125*a1*a2*x1
    pp(y)= max(real(0.,crm_rknd),y)
    pn(y)=-min(real(0.,crm_rknd),y)

    nonos = .true.
    eps = 1.e-10

    !$acc enter data create(mx,mn,uuu,vvv,www,iadz,irho,irhow) async(asyncid)

    !$acc parallel loop collapse(3) copy(www) async(asyncid)
    do icrm = 1 , ncrms
      do j = -1 , nyp2
        do i = -1 , nxp2
          www(i,j,nz,icrm)=0.
        enddo
      enddo
    enddo

    if (dowallx) then
      if (mod(rank,nsubdomains_x).eq.0) then
        !$acc parallel loop collapse(4) copy(u) async(asyncid)
        do icrm = 1 , ncrms
          do k=1,nzm
            do j=dimy1_u,dimy2_u
              do i=dimx1_u,1
                u(icrm,i,j,k) = 0.
              enddo
            enddo
          enddo
        enddo
      endif
      if (mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
        !$acc parallel loop collapse(4) copy(u) async(asyncid)
        do icrm = 1 , ncrms
          do k=1,nzm
            do j=dimy1_u,dimy2_u
              do i=nx+1,dimx2_u
                u(icrm,i,j,k) = 0.
              enddo
            enddo
          enddo
        enddo
      endif
    endif

    if (dowally) then
      if (rank.lt.nsubdomains_x) then
        !$acc parallel loop collapse(4) copy(v) async(asyncid)
        do icrm = 1 , ncrms
          do k=1,nzm
            do j=dimy1_v,1
              do i=dimx1_v,dimx2_v
                v(icrm,i,j,k) = 0.
              enddo
            enddo
          enddo
        enddo
      endif
      if (rank.gt.nsubdomains-nsubdomains_x-1) then
        !$acc parallel loop collapse(4) copy(v) async(asyncid)
        do icrm = 1 , ncrms
          do k=1,nzm
            do j=ny+1,dimy2_v
              do i=dimx1_v,dimx2_v
                v(icrm,i,j,k) = 0.
              enddo
            enddo
          enddo
        enddo
      endif
    endif

    !-----------------------------------------

    if (nonos) then
      !$acc parallel loop collapse(4) copyin(f) copy(mn,mx) async(asyncid)
      do icrm = 1 , ncrms
        do k=1,nzm
          do j=0,nyp1
            do i=0,nxp1
              kc=min(nzm,k+1)
              kb=max(1,k-1)
              jb=j-1
              jc=j+1
              ib=i-1
              ic=i+1
              mx(i,j,k,icrm)=max(f(icrm,ib,j,k),f(icrm,ic,j,k),f(icrm,i,jb,k),f(icrm,i,jc,k),&
                                 f(icrm,i,j,kb),f(icrm,i,j,kc),f(icrm,i,j,k))
              mn(i,j,k,icrm)=min(f(icrm,ib,j,k),f(icrm,ic,j,k),f(icrm,i,jb,k),f(icrm,i,jc,k),&
                                 f(icrm,i,j,kb),f(icrm,i,j,kc),f(icrm,i,j,k))
            enddo
          enddo
        enddo
      enddo
    endif  ! nonos

    !$acc parallel loop collapse(4) copyin(u,v,w,f) copy(uuu,vvv,www,flux) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=-1,nyp3
          do i=-1,nxp3
            kb=max(1,k-1)
            if (j <= nyp2                ) uuu(i,j,k,icrm)=max(real(0.,crm_rknd),u(icrm,i,j,k))*f(icrm,i-1,j,k)+&
                                                           min(real(0.,crm_rknd),u(icrm,i,j,k))*f(icrm,i,j,k)
            if (i <= nxp2                ) vvv(i,j,k,icrm)=max(real(0.,crm_rknd),v(icrm,i,j,k))*f(icrm,i,j-1,k)+&
                                                           min(real(0.,crm_rknd),v(icrm,i,j,k))*f(icrm,i,j,k)
            if (i <= nxp2 .and. j <= nyp2) www(i,j,k,icrm)=max(real(0.,crm_rknd),w(icrm,i,j,k))*f(icrm,i,j,kb )+&
                                                           min(real(0.,crm_rknd),w(icrm,i,j,k))*f(icrm,i,j,k)
            if (i == -1 .and. j == -1) then
              flux(k,icrm) = 0.
            endif
          enddo
        enddo
      enddo
    enddo

    !$acc parallel loop collapse(2) copyin(rho,adz,rhow) copy(irho,iadz,irhow) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        irho(k,icrm) = 1./rho(icrm,k)
        iadz(k,icrm) = 1./adz(icrm,k)
        irhow(k,icrm)=1./(rhow(k,icrm)*adz(icrm,k))
      enddo
    enddo

    !$acc parallel loop collapse(4) copyin(irho,iadz,www,uuu,vvv) copy(flux,f) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=-1,nyp2
          do i=-1,nxp2
            if (i == -1 .and. j == -1) then
            endif
            if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny) then
              !$acc atomic update
              flux(k,icrm) = flux(k,icrm) + www(i,j,k,icrm)
            endif
            f(icrm,i,j,k)=f(icrm,i,j,k)-( uuu(i+1,j,k,icrm)-uuu(i,j,k,icrm)  & 
                              + vvv(i,j+1,k,icrm)-vvv(i,j,k,icrm)  &
                              +(www(i,j,k+1,icrm)-www(i,j,k,icrm) )*iadz(k,icrm))*irho(k,icrm)
          enddo
        enddo
      enddo
    enddo

    !$acc parallel loop collapse(4) copyin(f,u,irho,v,w,irhow,adz) copy(uuu,vvv,www) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=0,nyp2
          do i=0,nxp2
            if (j <= nyp1) then
              kc=min(nzm,k+1)
              kb=max(1,k-1)
              dd=2./(kc-kb)/adz(icrm,k)
              jb=j-1
              jc=j+1
              ib=i-1
              uuu(i,j,k,icrm)=andiff(f(icrm,ib,j,k),f(icrm,i,j,k),u(icrm,i,j,k),irho(k,icrm)) &
              -(across(f(icrm,ib,jc,k)+f(icrm,i,jc,k)-f(icrm,ib,jb,k)-f(icrm,i,jb,k), &
              u(icrm,i,j,k), v(icrm,ib,j,k)+v(icrm,ib,jc,k)+v(icrm,i,jc,k)+v(icrm,i,j,k)) &
              +across(dd*(f(icrm,ib,j,kc)+f(icrm,i,j,kc)-f(icrm,ib,j,kb)-f(icrm,i,j,kb)), &
              u(icrm,i,j,k), w(icrm,ib,j,k)+w(icrm,ib,j,kc)+w(icrm,i,j,k)+w(icrm,i,j,kc))) *irho(k,icrm)
            endif
            if (i <= nxp1) then
              kc=min(nzm,k+1)
              kb=max(1,k-1)
              dd=2./(kc-kb)/adz(icrm,k)
              jb=j-1
              ib=i-1
              ic=i+1
              vvv(i,j,k,icrm)=andiff(f(icrm,i,jb,k),f(icrm,i,j,k),v(icrm,i,j,k),irho(k,icrm)) &
              -(across(f(icrm,ic,jb,k)+f(icrm,ic,j,k)-f(icrm,ib,jb,k)-f(icrm,ib,j,k), &
              v(icrm,i,j,k), u(icrm,i,jb,k)+u(icrm,i,j,k)+u(icrm,ic,j,k)+u(icrm,ic,jb,k)) &
              +across(dd*(f(icrm,i,jb,kc)+f(icrm,i,j,kc)-f(icrm,i,jb,kb)-f(icrm,i,j,kb)), &
              v(icrm,i,j,k), w(icrm,i,jb,k)+w(icrm,i,j,k)+w(icrm,i,j,kc)+w(icrm,i,jb,kc))) *irho(k,icrm)
            endif
            if (i <= nxp1 .and. j <= nyp1) then
              kb=max(1,k-1)
              jb=j-1
              jc=j+1
              ib=i-1
              ic=i+1
              www(i,j,k,icrm)=andiff(f(icrm,i,j,kb),f(icrm,i,j,k),w(icrm,i,j,k),irhow(k,icrm)) &
              -(across(f(icrm,ic,j,kb)+f(icrm,ic,j,k)-f(icrm,ib,j,kb)-f(icrm,ib,j,k), &
              w(icrm,i,j,k), u(icrm,i,j,kb)+u(icrm,i,j,k)+u(icrm,ic,j,k)+u(icrm,ic,j,kb)) &
              +across(f(icrm,i,jc,k)+f(icrm,i,jc,kb)-f(icrm,i,jb,k)-f(icrm,i,jb,kb), &
              w(icrm,i,j,k), v(icrm,i,j,kb)+v(icrm,i,jc,kb)+v(icrm,i,jc,k)+v(icrm,i,j,k))) *irho(k,icrm)
            endif
          enddo
        enddo
      enddo
    enddo

    !$acc parallel loop collapse(3) copy(www) async(asyncid)
    do icrm = 1 , ncrms
      do j = -1 , nyp2
        do i = -1 , nxp2
          www(i,j,1,icrm) = 0.
        enddo
      enddo
    enddo

    !---------- non-osscilatory option ---------------
    if (nonos) then
      !$acc parallel loop collapse(4) copyin(f) copy(mx,mn) async(asyncid)
      do icrm = 1 , ncrms
        do k=1,nzm
          do j=0,nyp1
            do i=0,nxp1
              kc=min(nzm,k+1)
              kb=max(1,k-1)
              jb=j-1
              jc=j+1
              ib=i-1
              ic=i+1
              mx(i,j,k,icrm)=max(f(icrm,ib,j,k),f(icrm,ic,j,k),f(icrm,i,jb,k),f(icrm,i,jc,k),&
                                 f(icrm,i,j,kb),f(icrm,i,j,kc),f(icrm,i,j,k),mx(i,j,k,icrm))
              mn(i,j,k,icrm)=min(f(icrm,ib,j,k),f(icrm,ic,j,k),f(icrm,i,jb,k),f(icrm,i,jc,k),&
                                 f(icrm,i,j,kb),f(icrm,i,j,kc),f(icrm,i,j,k),mn(i,j,k,icrm))
            enddo
          enddo
        enddo
      enddo

      !$acc parallel loop collapse(4) copyin(rho,f,uuu,vvv,www,iadz) copy(mx,mn) async(asyncid)
      do icrm = 1 , ncrms
        do k=1,nzm
          do j=0,nyp1
            do i=0,nxp1
              kc=min(nzm,k+1)
              jc=j+1
              ic=i+1
              mx(i,j,k,icrm)=rho(icrm,k)*(mx(i,j,k,icrm)-f(icrm,i,j,k))/ &
                        ( pn(uuu(ic,j,k,icrm)) + pp(uuu(i,j,k,icrm))+ &
                          pn(vvv(i,jc,k,icrm)) + pp(vvv(i,j,k,icrm))+ &
                         (pn(www(i,j,kc,icrm)) + pp(www(i,j,k,icrm)))*iadz(k,icrm)+eps)
              mn(i,j,k,icrm)=rho(icrm,k)*(f(icrm,i,j,k)-mn(i,j,k,icrm))/ &
                        ( pp(uuu(ic,j,k,icrm)) + pn(uuu(i,j,k,icrm))+ &
                          pp(vvv(i,jc,k,icrm)) + pn(vvv(i,j,k,icrm))+ &
                         (pp(www(i,j,kc,icrm)) + pn(www(i,j,k,icrm)))*iadz(k,icrm)+eps)
            enddo
          enddo
        enddo
      enddo

      !$acc parallel loop collapse(4) copyin(mx,mn) copy(uuu,vvv,www,flux) async(asyncid)
      do icrm = 1 , ncrms
        do k=1,nzm
          do j=1,nyp1
            do i=1,nxp1
              if (j <= ny) then
                ib=i-1
                uuu(i,j,k,icrm) = pp(uuu(i,j,k,icrm))*min(real(1.,crm_rknd),mx(i,j,k,icrm), mn(ib,j,k,icrm)) &
                                 -pn(uuu(i,j,k,icrm))*min(real(1.,crm_rknd),mx(ib,j,k,icrm),mn(i,j,k,icrm))
              endif
              if (i <= nx) then
                jb=j-1
                vvv(i,j,k,icrm) = pp(vvv(i,j,k,icrm))*min(real(1.,crm_rknd),mx(i,j,k,icrm), mn(i,jb,k,icrm)) &
                                 -pn(vvv(i,j,k,icrm))*min(real(1.,crm_rknd),mx(i,jb,k,icrm),mn(i,j,k,icrm))
              endif
              if (i <= nx .and. j <= ny) then
                kb=max(1,k-1)
                www(i,j,k,icrm) = pp(www(i,j,k,icrm))*min(real(1.,crm_rknd),mx(i,j,k,icrm), mn(i,j,kb,icrm)) &
                                 -pn(www(i,j,k,icrm))*min(real(1.,crm_rknd),mx(i,j,kb,icrm),mn(i,j,k,icrm))
                !$acc atomic update
                flux(k,icrm) = flux(k,icrm) + www(i,j,k,icrm)
              endif
            enddo
          enddo
        enddo
      enddo
    endif ! nonos

    !$acc parallel loop collapse(4) copyin(uuu,vvv,www,w,irho,iadz) copy(f) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            ! MK: added fix for very small negative values (relative to positive values)
            !     especially  when such large numbers as
            !     hydrometeor concentrations are advected. The reason for negative values is
            !     most likely truncation error.
            kc=k+1
            f(icrm,i,j,k)=max(real(0.,crm_rknd),f(icrm,i,j,k) -(uuu(i+1,j,k,icrm)-uuu(i,j,k,icrm)+&
                             vvv(i,j+1,k,icrm)-vvv(i,j,k,icrm)+(www(i,j,k+1,icrm)-www(i,j,k,icrm))*iadz(k,icrm))*irho(k,icrm))
          enddo
        enddo
      enddo
    enddo

    !$acc exit data delete(mx,mn,uuu,vvv,www,iadz,irho,irhow) async(asyncid)

  end subroutine advect_scalar3D

end module advect_scalar3D_mod
