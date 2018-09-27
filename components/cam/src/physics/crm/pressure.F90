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
    real(crm_rknd) f(nx2,ny2,nzslab) ! global rhs and array for FTP coefficeients
    real(crm_rknd) ff(nx+1,ny+2*YES3D,nzm)  ! local (subdomain's) version of f
    real(crm_rknd) work(nx2,ny2),trigxi(n3i),trigxj(n3j) ! FFT stuff
    integer ifaxj(100),ifaxi(100)
    real(8) a(nzm),b,c(nzm),e,fff(nzm)
    real(8) xi,xj,xnx,xny,ddx2,ddy2,pii,factx,facty,eign
    real(8) alfa(nzm-1),beta(nzm-1)
    integer i, j, k, id, jd, m, n, it, jt, ii, jj, icrm
    integer nyp22
    integer iii(0:nx_gl),jjj(0:ny_gl)
    integer iwall,jwall

    ! check if the grid size allows the computation:
    if(nsubdomains.gt.nzm) then
      if(masterproc) print*,'pressure_orig: nzm < nsubdomains. STOP'
      call task_abort
    endif

    if(mod(nzm,npressureslabs).ne.0) then
      if(masterproc) print*,'pressure_orig: nzm/npressureslabs is not round number. STOP'
      call task_abort
    endif

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

    !-----------------------------------------------------------------
    !  Compute the r.h.s. of the Poisson equation for pressure
    call press_rhs(ncrms)

    do icrm = 1 , ncrms
      !-----------------------------------------------------------------
      !   Form the horizontal slabs of right-hand-sides of Poisson equation
      !   for the global domain. Request sending and receiving tasks.
      n = rank*nzslab
      do k = 1,nzslab
        do j = 1,ny
          do i = 1,nx
            f(i+it,j+jt,k) = p(i,j,k+n,icrm)
          enddo
        enddo
      enddo

      !-------------------------------------------------
      ! Perform Fourier transformation for a slab:
      if(rank.lt.npressureslabs) then
        call fftfax_crm(nx_gl,ifaxi,trigxi)
        if(RUN3D) call fftfax_crm(ny_gl,ifaxj,trigxj)
        do k=1,nzslab
          call fft991_crm(f(1,1,k),work,trigxi,ifaxi,1,nx2,nx_gl,ny_gl,-1)
          if(RUN3D) then
            call fft991_crm(f(1,1,k),work,trigxj,ifaxj,nx2,1,ny_gl,nx_gl+1,-1)
          endif
        enddo
      endif

      !-------------------------------------------------
      !   Send Fourier coeffiecients back to subdomains:
      n = rank*nzslab
      do k = 1,nzslab
        do j = 1,nyp22-jwall
          do i = 1,nxp1-iwall
            ff(i,j,k+n) = f(i+it,j+jt,k)
          enddo
        enddo
      enddo

      !-------------------------------------------------
      !   Solve the tri-diagonal system for Fourier coeffiecients
      !   in the vertical for each subdomain:
      do k=1,nzm
        a(k)=rhow(k,icrm)/(adz(k,icrm)*adzw(k,icrm)*dz(icrm)*dz(icrm))
        c(k)=rhow(k+1,icrm)/(adz(k,icrm)*adzw(k+1,icrm)*dz(icrm)*dz(icrm))
      enddo

      ddx2=1._8/(dx*dx)
      ddy2=1._8/(dy*dy)
      pii = acos(-1._8)
      xnx=pii/nx_gl
      xny=pii/ny_gl
      do j=1,nyp22-jwall
        if(dowally) then
          jd=j+jt-1
          facty = 1.d0
        else
          jd=(j+jt-0.1)/2.
          facty = 2.d0
        endif
        xj=jd
        do i=1,nxp1-iwall
          if(dowallx) then
            id=i+it-1
            factx = 1.d0
          else
            id=(i+it-0.1)/2.
            factx = 2.d0
          endif
          fff(1:nzm) = ff(i,j,1:nzm)
          xi=id
          eign=(2._8*cos(factx*xnx*xi)-2._8)*ddx2+(2._8*cos(facty*xny*xj)-2._8)*ddy2
          if(id+jd.eq.0) then
            b=1._8/(eign*rho(1,icrm)-a(1)-c(1))
            alfa(1)=-c(1)*b
            beta(1)=fff(1)*b
          else
            b=1._8/(eign*rho(1,icrm)-c(1))
            alfa(1)=-c(1)*b
            beta(1)=fff(1)*b
          endif
          do k=2,nzm-1
            e=1._8/(eign*rho(k,icrm)-a(k)-c(k)+a(k)*alfa(k-1))
            alfa(k)=-c(k)*e
            beta(k)=(fff(k)-a(k)*beta(k-1))*e
          enddo
          fff(nzm)=(fff(nzm)-a(nzm)*beta(nzm-1))/(eign*rho(nzm,icrm)-a(nzm)+a(nzm)*alfa(nzm-1))
          do k=nzm-1,1,-1
            fff(k)=alfa(k)*fff(k+1)+beta(k)
          enddo
          ff(i,j,1:nzm) = fff(1:nzm)
        enddo
      enddo

      !-----------------------------------------------------------------
      !   Send the Fourier coefficient to the tasks performing
      !   the inverse Fourier transformation:
      n = rank*nzslab
      do k = 1,nzslab
        do j = 1,nyp22-jwall
          do i = 1,nxp1-iwall
            f(i+it,j+jt,k) = ff(i,j,k+n)
          enddo
        enddo
      enddo

      !-------------------------------------------------
      !   Perform inverse Fourier transformation:
      if(rank.lt.npressureslabs) then
        do k=1,nzslab
          if(RUN3D) then
            call fft991_crm(f(1,1,k),work,trigxj,ifaxj,nx2,1,ny_gl,nx_gl+1,+1)
          endif
          call fft991_crm(f(1,1,k),work,trigxi,ifaxi,1,nx2,nx_gl,ny_gl,+1)
        enddo
      endif

      !-----------------------------------------------------------------
      !   Fill the pressure field for each subdomain:
      do i=1,nx_gl
        iii(i)=i
      enddo
      iii(0)=nx_gl
      do j=1,ny_gl
        jjj(j)=j
      enddo
      jjj(0)=ny_gl

      n = rank*nzslab
      do k = 1,nzslab
        do j = 1-YES3D,ny
          jj=jjj(j+jt)
          do i = 0,nx
            ii=iii(i+it)
            p(i,j,k+n,icrm) = f(ii,jj,k)
          enddo
        enddo
      enddo

    enddo

    !  Add pressure gradient term to the rhs of the momentum equation:
    call press_grad(ncrms)

  end subroutine pressure

end module pressure_mod
