module pressure_mod
  use task_util_mod
  implicit none

contains

  ! Non-blocking receives before blocking sends
  subroutine pressure(ncrms)
    !       Original pressure solver based on horizontal slabs
    !       (C) 1998, 2002 Marat Khairoutdinov
    !       Works only when the number of slabs is equal to the number of processors.
    !       Therefore, the number of processors shouldn't exceed the number of levels nzm
    !       Also, used for a 2D version
    !       For more processors for the given number of levels and 3D, use pressure_big
    use vars
    use params, only: dowallx, dowally, docolumn, crm_rknd
    use press_rhs_mod
    use press_grad_mod
    use fft_mod
    implicit none
    integer, intent(in) :: ncrms
    integer, parameter :: npressureslabs = nsubdomains
    integer, parameter :: nzslab = max(1,nzm / npressureslabs)
    integer, parameter :: nx2=nx_gl+2, ny2=ny_gl+2*YES3D
    integer, parameter :: n3i=3*nx_gl/2+1,n3j=3*ny_gl/2+1
    real(crm_rknd) f(nx2,ny2,nzslab,ncrms) ! global rhs and array for FTP coefficeients
    real(crm_rknd) ff(nx+1,ny+2*YES3D,nzm,ncrms)  ! local (subdomain's) version of f
    real(crm_rknd) work(nx2,ny2),trigxi(n3i),trigxj(n3j) ! FFT stuff
    integer ifaxj(100),ifaxi(100)
    real(8) a(nzm,ncrms),b,c(nzm,ncrms),e
    real(8) xi,xj,xnx,xny,ddx2,ddy2,pii,factx,facty
    real(8) alfa(nzm-1),beta(nzm-1)
    integer i, j, k, id, jd, m, n, it, jt, ii, jj, icrm
    integer nyp22
    integer iii(0:nx_gl),jjj(0:ny_gl)
    integer iwall,jwall
    integer :: numgangs  !For working aroung PGI OpenACC bug where it didn't create enough gangs
    real(8), allocatable :: eign(:,:)

    !$acc enter data create(iii,jjj,f,ff,trigxi,trigxj,ifaxi,ifaxj,a,c) async(asyncid)

    it = 0
    jt = 0

    !-----------------------------------------------------------------

    if(docolumn) return

    if(dowallx) then
      iwall=1
    else
      iwall=0
    endif
    if(RUN2D) then
      nyp22=1
      jwall=0
    else
      nyp22=nyp2
      if(dowally) then
        jwall=2
      else
        jwall=0
      endif
    endif

    allocate(eign(nxp1-iwall,nyp22-jwall))
    !$acc enter data create(eign) async(asyncid)

    !-----------------------------------------------------------------
    !  Compute the r.h.s. of the Poisson equation for pressure
    call press_rhs(ncrms)

    !-----------------------------------------------------------------
    !   Form the horizontal slabs of right-hand-sides of Poisson equation
    n = 0
    !$acc parallel loop collapse(4) copyin(p) copyout(f) async(asyncid)
    do icrm = 1 , ncrms
      do k = 1,nzslab
        do j = 1,ny
          do i = 1,nx
            f(i,j,k,icrm) = p(i,j,k,icrm)
          enddo
        enddo
      enddo
    enddo

    !-------------------------------------------------
    ! Perform Fourier transformation for a slab:
    !$acc parallel loop copyout(ifaxi,trigxi,ifaxj,trigxj) async(asyncid)
    do icrm = 1 , 1
      call fftfax_crm(nx_gl,ifaxi,trigxi)
      if(RUN3D) call fftfax_crm(ny_gl,ifaxj,trigxj)
    enddo
    !$acc parallel loop collapse(2) copyin(trigxi,ifaxi,trigxj,ifaxj) copy(f) private(work) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzslab
        call fft991_crm(f(1,1,k,icrm),work,trigxi,ifaxi,1,nx2,nx_gl,ny_gl,-1)
        if(RUN3D) then
          call fft991_crm(f(1,1,k,icrm),work,trigxj,ifaxj,nx2,1,ny_gl,nx_gl+1,-1)
        endif
      enddo
    enddo

    !-------------------------------------------------
    !   Send Fourier coeffiecients back to subdomains:
    !$acc parallel loop collapse(4) copyin(f) copyout(ff) async(asyncid)
    do icrm = 1 , ncrms
      do k = 1,nzslab
        do j = 1,nyp22-jwall
          do i = 1,nxp1-iwall
            ff(i,j,k,icrm) = f(i,j,k,icrm)
          enddo
        enddo
      enddo
    enddo

    !-------------------------------------------------
    !   Solve the tri-diagonal system for Fourier coeffiecients
    !   in the vertical for each subdomain:
    !$acc parallel loop collapse(2) copyin(adz,adzw,dz,rhow) copyout(a,c) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        a(k,icrm)=rhow(k,icrm)/(adz(icrm,k)*adzw(icrm,k)*dz(icrm)*dz(icrm))
        c(k,icrm)=rhow(k+1,icrm)/(adz(icrm,k)*adzw(icrm,k+1)*dz(icrm)*dz(icrm))
      enddo
    enddo

      do j=1,nyp22-jwall
        do i=1,nxp1-iwall
      enddo
    enddo

    !$acc parallel loop collapse(2) copy(eign) async(asyncid)
    do j=1,nyp22-jwall
      do i=1,nxp1-iwall
        ddx2=1._8/(dx*dx)
        ddy2=1._8/(dy*dy)
        pii = 3.14159265358979323846D0
        xnx=pii/nx_gl
        xny=pii/ny_gl
        if(dowally) then
          jd=j+jt-1
          facty = 1.d0
        else
          jd=(j+jt-0.1)/2.
          facty = 2.d0
        endif
        xj=jd
        if(dowallx) then
          id=i+it-1
          factx = 1.d0
        else
          id=(i+it-0.1)/2.
          factx = 2.d0
        endif
        xi=id
        eign(i,j)=(2._8*cos(factx*xnx*xi)-2._8)*ddx2+(2._8*cos(facty*xny*xj)-2._8)*ddy2
      enddo
    enddo

    !For working aroung PGI OpenACC bug where it didn't create enough gangs
    numgangs = ceiling(ncrms*(nyp22-jwall)*(nxp2-iwall)/128.)
    !$acc parallel loop collapse(3) vector_length(128) num_gangs(numgangs) private(alfa,beta) copyin(a,c,rho,eign) copy(ff) async(asyncid)
    do icrm = 1 , ncrms
      do j=1,nyp22-jwall
        do i=1,nxp1-iwall
          if(dowally) then
            jd=j+jt-1
          else
            jd=(j+jt-0.1)/2.
          endif
          if(dowallx) then
            id=i+it-1
          else
            id=(i+it-0.1)/2.
          endif
          if(id+jd.eq.0) then
            b=1._8/(eign(i,j)*rho(icrm,1)-a(1,icrm)-c(1,icrm))
            alfa(1)=-c(1,icrm)*b
            beta(1)=ff(i,j,1,icrm)*b
          else
            b=1._8/(eign(i,j)*rho(icrm,1)-c(1,icrm))
            alfa(1)=-c(1,icrm)*b
            beta(1)=ff(i,j,1,icrm)*b
          endif
          do k=2,nzm-1
            e=1._8/(eign(i,j)*rho(icrm,k)-a(k,icrm)-c(k,icrm)+a(k,icrm)*alfa(k-1))
            alfa(k)=-c(k,icrm)*e
            beta(k)=(ff(i,j,k,icrm)-a(k,icrm)*beta(k-1))*e
          enddo
          ff(i,j,nzm,icrm)=(ff(i,j,nzm,icrm)-a(nzm,icrm)*beta(nzm-1))/(eign(i,j)*rho(icrm,nzm)-a(nzm,icrm)+a(nzm,icrm)*alfa(nzm-1))
          do k=nzm-1,1,-1
            ff(i,j,k,icrm)=alfa(k)*ff(i,j,k+1,icrm)+beta(k)
          enddo
        enddo
      enddo
    enddo

    !-----------------------------------------------------------------
    n = 0
    !$acc parallel loop collapse(4) copyin(ff) copyout(f) async(asyncid)
    do icrm = 1 , ncrms
      do k = 1,nzslab
        do j = 1,nyp22-jwall
          do i = 1,nxp1-iwall
            f(i,j,k,icrm) = ff(i,j,k,icrm)
          enddo
        enddo
      enddo
    enddo

    !-------------------------------------------------
    !   Perform inverse Fourier transformation:
    !$acc parallel loop collapse(2) copyin(trigxi,ifaxi,trigxj,ifaxj) copy(f) private(work) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzslab
        if(RUN3D) then
          call fft991_crm(f(1,1,k,icrm),work,trigxj,ifaxj,nx2,1,ny_gl,nx_gl+1,+1)
        endif
        call fft991_crm(f(1,1,k,icrm),work,trigxi,ifaxi,1,nx2,nx_gl,ny_gl,+1)
      enddo
    enddo

    !-----------------------------------------------------------------
    !   Fill the pressure field for each subdomain:
    !$acc parallel loop copyout(iii,jjj) async(asyncid)
    do icrm = 1,1
      do i=1,nx_gl
        iii(i)=i
      enddo
      iii(0)=nx_gl
      do j=1,ny_gl
        jjj(j)=j
      enddo
      jjj(0)=ny_gl
    enddo

    n = 0
    !$acc parallel loop collapse(4) copyin(iii,jjj,f) copyout(p) async(asyncid)
    do icrm = 1 , ncrms
      do k = 1,nzslab
        do j = 1-YES3D,ny
          do i = 0,nx
            jj=jjj(j)
            ii=iii(i)
            p(i,j,k,icrm) = f(ii,jj,k,icrm)
          enddo
        enddo
      enddo
    enddo

    !  Add pressure gradient term to the rhs of the momentum equation:
    call press_grad(ncrms)

    !$acc exit data delete(eign) async(asyncid)
    deallocate(eign)

    !$acc exit data delete(iii,jjj,f,ff,trigxi,trigxj,ifaxi,ifaxj,a,c) async(asyncid)

  end subroutine pressure

end module pressure_mod
