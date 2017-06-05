
subroutine advect_scalar2D (f, u, w, rho, rhow, flux)
  !     positively definite monotonic advection with non-oscillatory option
  use grid
  use params, only: dowallx
  implicit none
  real    , intent(inout) :: f(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
  real    , intent(inout) :: u(dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm)
  real    , intent(in   ) :: w(dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz )
  real    , intent(in   ) :: rho (nzm)
  real    , intent(in   ) :: rhow(nz )
  real    , intent(  out) :: flux(nz )
  real    :: mx ( 0:nxp1,1,nzm)
  real    :: mn ( 0:nxp1,1,nzm)
  real    :: uuu(-1:nxp3,1,nzm)
  real    :: www(-1:nxp2,1,nz )
  real    :: eps, dd
  integer :: i,j,k,ic,ib,kc,kb
  logical :: nonos
  real    :: iadz(nzm),irho(nzm),irhow(nzm)
  real    :: x1, x2, a, b, a1, a2, y
  real    :: andiff,across,pp,pn

  !Statement functions
  andiff(x1,x2,a,b) = (abs(a)-a*a*b)*0.5*(x2-x1)
  across(x1,a1,a2)  = 0.03125*a1*a2*x1
  pp    (y)         =  max(0.,y)
  pn    (y)         = -min(0.,y)
  
  !Because we're in 2-D
  j=1

  !Initialization
  nonos = .true.
  eps = 1.e-10

  !$acc enter data pcreate(uuu,www,u,mn,mx,f,w,irho,rho,iadz,adz,flux,rhow,irhow) async(1)
  !$acc update device(f,u,w,rho,adz,rhow) async(1)

  !$acc parallel loop gang vector present(www) async(1)
  do i = -1 , nxp2
    www(i,j,nz)=0.
  enddo
  
  if (dowallx) then
    if ( mod(rank,nsubdomains_x) == 0 ) then
      !$acc parallel loop gang vector collapse(2) present(u) async(1)
      do k = 1 , nzm
        do i = dimx1_u , 1
          u(i,j,k) = 0.
        enddo
      enddo
    endif
    if ( mod(rank,nsubdomains_x) == nsubdomains_x-1 ) then
      !$acc parallel loop gang vector collapse(2) present(u) async(1)
      do k = 1 , nzm
        do i = nx+1 , dimx2_u
          u(i,j,k) = 0.
        enddo
      enddo
    endif
    !$acc update host(u) async(1)
  endif
  
  if (nonos) then
    !$acc parallel loop gang vector collapse(2) present(mx,mn,f) private(kc,kb,ib,ic) async(1)
    do k = 1 , nzm
      do i = 0 , nxp1
        kc = min(nzm,k+1)
        kb = max(1  ,k-1)
        ib = i-1
        ic = i+1
        mx(i,j,k) = max( f(ib,j,k) , f(ic,j,k) , f(i,j,kb) , f(i,j,kc) , f(i,j,k) )
        mn(i,j,k) = min( f(ib,j,k) , f(ic,j,k) , f(i,j,kb) , f(i,j,kc) , f(i,j,k) )
      enddo
    enddo
  endif  ! nonos
  
  !$acc parallel loop gang vector collapse(2) present(uuu,u,f,w,www) private(kb) async(1)
  do k = 1 , nzm
    do i = -1 , nxp3
      kb = max(1,k-1)
      uuu(i,j,k) = max(0.,u(i,j,k)) * f(i-1,j,k ) + min(0.,u(i,j,k)) * f(i,j,k)
      if (i <= nxp2) www(i,j,k) = max(0.,w(i,j,k)) * f(i  ,j,kb) + min(0.,w(i,j,k)) * f(i,j,k)
    enddo
  enddo

  !$acc parallel loop gang vector present(www,flux) async(1)
  do k = 1 , nzm
    flux(k) = 0.
    !$acc loop seq
    do i = 1 , nx
      flux(k) = flux(k) + www(i,j,k)
    enddo
  enddo
  
  !$acc parallel loop gang vector collapse(2) present(irho,rho,iadz,adz,f,uuu,www) async(1)
  do k=1,nzm
    do i=-1,nxp2
      irho(k) = 1./rho(k)
      iadz(k) = 1./adz(k)
      f(i,j,k) = f(i,j,k) - (uuu(i+1,j,k)-uuu(i,j,k) + (www(i,j,k+1)-www(i,j,k))*iadz(k))*irho(k)            
    enddo
  enddo 
  
  !$acc parallel loop gang vector collapse(2) present(irhow,uuu,f,u,irho,w,www,adz,rhow) private(kc,kb,dd,ib) async(1)
  do k=1,nzm
    do i=0,nxp2
      kc=min(nzm,k+1)
      kb=max(1,k-1)
      dd=2./(kc-kb)/adz(k)
      irhow(k)=1./(rhow(k)*adz(k))
      ib=i-1
      uuu(i,j,k) = andiff(f(ib,j,k),f(i,j,k),u(i,j,k),irho(k)) - &
                   across(dd*(f(ib,j,kc)+f(i,j,kc)-f(ib,j,kb)-f(i,j,kb))  ,  &
                          u(i,j,k), w(ib,j,k)+w(ib,j,kc)+w(i,j,k)+w(i,j,kc)) * irho(k)
      if (i <= nxp1) then
        ib=i-1
        ic=i+1
        www(i,j,k)=andiff(f(i,j,kb),f(i,j,k),w(i,j,k),irhow(k)) - &
                   across(f(ic,j,kb)+f(ic,j,k)-f(ib,j,kb)-f(ib,j,k), &
                          w(i,j,k), u(i,j,kb)+u(i,j,k)+u(ic,j,k)+u(ic,j,kb)) * irho(k)
      endif
    enddo
  enddo

  !$acc parallel loop gang vector present(www) async(1)
  do i = -1 , nxp2
    www(i,j,1)=0.
  enddo
  
  if(nonos) then
    !$acc parallel loop gang vector collapse(2) present(mx,f,mn) private(kc,kb,ib,ic) async(1)
    do k=1,nzm
      do i=0,nxp1
        kc=min(nzm,k+1)
        kb=max(1  ,k-1)
        ib=i-1
        ic=i+1
        mx(i,j,k) = max( f(ib,j,k) , f(ic,j,k) , f(i,j,kb) , f(i,j,kc) , f(i,j,k) , mx(i,j,k) )
        mn(i,j,k) = min( f(ib,j,k) , f(ic,j,k) , f(i,j,kb) , f(i,j,kc) , f(i,j,k) , mn(i,j,k) )
      enddo
    enddo
    !$acc parallel loop gang vector collapse(2) present(mx,mn,rho,f,uuu,www,iadz) private(kc,ic) async(1)
    do k=1,nzm
      do i=0,nxp1
        kc=min(nzm,k+1)
        ic=i+1
        mx(i,j,k)=rho(k)*(mx(i,j,k)-f(i,j,k))/(pn(uuu(ic,j,k)) + pp(uuu(i,j,k))+&
                  iadz(k)*(pn(www(i,j,kc)) + pp(www(i,j,k)))+eps)
        mn(i,j,k)=rho(k)*(f(i,j,k)-mn(i,j,k))/(pp(uuu(ic,j,k)) + pn(uuu(i,j,k))+&
                  iadz(k)*(pp(www(i,j,kc)) + pn(www(i,j,k)))+eps)
      enddo
    enddo
    !$acc parallel loop gang vector collapse(2) present(uuu,mx,mn,www) private(kb,ib) async(1)
    do k=1,nzm
      do i=1,nxp1
        kb=max(1,k-1)
        ib=i-1
        uuu(i,j,k)= pp(uuu(i,j,k))*min(1.,mx(i,j,k), mn(ib,j,k)) - pn(uuu(i,j,k))*min(1.,mx(ib,j,k),mn(i,j,k))
        if (i <= nx) then
          www(i,j,k)= pp(www(i,j,k))*min(1.,mx(i,j,k), mn(i,j,kb)) - pn(www(i,j,k))*min(1.,mx(i,j,kb),mn(i,j,k))
        endif
      enddo
    enddo
    !$acc parallel loop gang vector present(flux,www) async(1)
    do k=1,nzm
      do i=1,nx
        flux(k) = flux(k) + www(i,j,k)
      enddo
    enddo
  endif ! nonos

  !$acc parallel loop gang vector collapse(2) present(f,uuu,www,iadz,irho) private(kc) async(1)
  do k=1,nzm
    do i=1,nx
      kc=k+1
      ! MK: added fix for very small negative values (relative to positive values) 
      !     especially  when such large numbers as
      !     hydrometeor concentrations are advected. The reason for negative values is
      !     most likely truncation error.
      f(i,j,k)= max(0., f(i,j,k) - (uuu(i+1,j,k)-uuu(i,j,k) + (www(i,j,k+1)-www(i,j,k))*iadz(k))*irho(k))
    enddo
  enddo 

  !$acc update host(f,flux) async(1)
  !$acc wait(1)

end subroutine advect_scalar2D


