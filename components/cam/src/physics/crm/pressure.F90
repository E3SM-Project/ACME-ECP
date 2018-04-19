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

    use vars, only: nzm, nx_gl, ny_gl, yes3d, nsubdomains, nx, ny, masterproc, run2d, nyp2, rank, run3d, nxp1, dx, dy, &
                    p, adz, adzw, dz, rho, rhow
                    !p: inout
                    !adz: in
                    !adzw: in
                    !dz: in
                    !rho: in
                    !rhow: in
    use params, only: dowallx, dowally, docolumn, crm_rknd
    use press_rhs_mod
    use press_grad_mod
    use vars
    implicit none
    integer, intent(in) :: ncrms


    integer, parameter :: npressureslabs = nsubdomains
    integer, parameter :: nzslab = max(1,nzm / npressureslabs)
    integer, parameter :: nx2=nx_gl+2, ny2=ny_gl+2*YES3D
    integer, parameter :: n3i=3*nx_gl/2+1,n3j=3*ny_gl/2+1

    real(crm_rknd), allocatable :: f      (:,:,:,:) ! global rhs and array for FTP coefficeients
    real(crm_rknd), allocatable :: ff     (:,:,:,:)  ! local (subdomain's) version of f
    real(crm_rknd), allocatable :: work   (:,:)
    real(crm_rknd), allocatable :: trigxi (:)
    real(crm_rknd), allocatable :: trigxj (:)
    integer       , allocatable :: ifaxj  (:)
    integer       , allocatable :: ifaxi  (:)
    real(8)       , allocatable :: a      (:,:)
    real(8)       , allocatable :: c      (:,:)
    real(8)       , allocatable :: fff    (:)
    real(8)       , allocatable :: alfa   (:)
    real(8)       , allocatable :: beta   (:)
    integer       , allocatable :: reqs_in(:)
    integer       , allocatable :: iii    (:)
    integer       , allocatable :: jjj    (:)
    logical       , allocatable :: flag   (:)

    real(8) xi,xj,xnx,xny,ddx2,ddy2,pii,factx,facty,eign,e,b
    integer i, j, k, id, jd, m, n, it, jt, ii, jj, tag, rf, icrm
    integer nyp22, n_in, count
    integer iwall,jwall

    allocate( f      (ncrms,nx2,ny2,nzslab)      )
    allocate( ff     (ncrms,nx+1,ny+2*YES3D,nzm) )
    allocate( work   (nx2,ny2)                   )
    allocate( trigxi (n3i)                       )
    allocate( trigxj (n3j)                       )
    allocate( ifaxj  (100)                       )
    allocate( ifaxi  (100)                       )
    allocate( a      (ncrms,nzm)                 )
    allocate( c      (ncrms,nzm)                 )
    allocate( fff    (nzm)                       )
    allocate( alfa   (nzm-1)                     )
    allocate( beta   (nzm-1)                     )
    allocate( reqs_in(nsubdomains)               )
    allocate( iii    (0:nx_gl)                   )
    allocate( jjj    (0:ny_gl)                   )
    allocate( flag   (nsubdomains)               )

    !$acc enter data create(f,ff,work,trigxi,trigxj,ifaxj,ifaxi,a,c,fff,alfa,beta,reqs_in,iii,jjj,flag) async(1)

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
    end if
    if(RUN2D) then
      nyp22=1
      jwall=0
    else
      nyp22=nyp2
      if(dowally) then
        jwall=2
      else
        jwall=0
      end if
    endif

    !-----------------------------------------------------------------
    !  Compute the r.h.s. of the Poisson equation for pressure

    call press_rhs(ncrms)

    !-----------------------------------------------------------------
    !   Form the horizontal slabs of right-hand-sides of Poisson equation
    !   for the global domain. Request sending and receiving tasks.

    n = rank*nzslab
    !$acc parallel loop gang vector collapse(4) default(present) async(1)
    do k = 1,nzslab
      do j = 1,ny
        do i = 1,nx
          do icrm = 1 , ncrms
            f(icrm,i,j,k) = p(icrm,i,j,k+n)
          enddo
        end do
      end do
    end do

    !$acc update host(f) async(1)
    !$acc wait(1)

    !-------------------------------------------------
    ! Perform Fourier transformation for a slab:
    call fftfax_crm(nx_gl,ifaxi,trigxi)
    if(RUN3D) call fftfax_crm(ny_gl,ifaxj,trigxj)
    do k=1,nzslab
      do icrm = 1 , ncrms
        call fft991_crm(f(icrm,:,:,k),work,trigxi,ifaxi,1,nx2,nx_gl,ny_gl,-1)
        if(RUN3D) call fft991_crm(f(icrm,:,:,k),work,trigxj,ifaxj,nx2,1,ny_gl,nx_gl+1,-1)
      end do
    end do

    !$acc update device(f) async(1)

    !-------------------------------------------------
    !   Send Fourier coeffiecients back to subdomains:
    n = rank*nzslab
    !$acc parallel loop gang vector collapse(4) default(present) async(1)
    do k = 1,nzslab
      do j = 1,nyp22-jwall
        do i = 1,nxp1-iwall
          do icrm = 1 , ncrms
            ff(icrm,i,j,k+n) = f(icrm,i+it,j+jt,k)
          end do
        end do
      end do
    end do

    !-------------------------------------------------
    !   Solve the tri-diagonal system for Fourier coeffiecients
    !   in the vertical for each subdomain:

    !$acc parallel loop gang vector collapse(2) default(present) async(1)
    do k=1,nzm
      do icrm = 1 , ncrms
        a(icrm,k)=rhow(icrm,k)/(adz(icrm,k)*adzw(icrm,k)*dz(icrm)*dz(icrm))
        c(icrm,k)=rhow(icrm,k+1)/(adz(icrm,k)*adzw(icrm,k+1)*dz(icrm)*dz(icrm))
      end do
    end do

    ddx2=1._8/(dx*dx)
    ddy2=1._8/(dy*dy)
    pii = acos(-1._8)
    xnx=pii/nx_gl
    xny=pii/ny_gl
    !$acc parallel loop gang vector collapse(3) private(fff,alfa,beta) default(present) async(1)
    do j=1,nyp22-jwall
      do i=1,nxp1-iwall
        do icrm = 1 , ncrms
          if(dowally) then
            jd=j+jt-1
            facty = 1.d0
          else
            jd=(j+jt-0.1)/2.
            facty = 2.d0
          end if
          xj=jd
          if(dowallx) then
            id=i+it-1
            factx = 1.d0
          else
            id=(i+it-0.1)/2.
            factx = 2.d0
          end if
          fff(1:nzm) = ff(icrm,i,j,1:nzm)
          xi=id
          eign=(2._8*cos(factx*xnx*xi)-2._8)*ddx2+ &
          (2._8*cos(facty*xny*xj)-2._8)*ddy2
          if(id+jd.eq.0) then
            b=1._8/(eign*rho(icrm,1)-a(icrm,1)-c(icrm,1))
            alfa(1)=-c(icrm,1)*b
            beta(1)=fff(1)*b
          else
            b=1._8/(eign*rho(icrm,1)-c(icrm,1))
            alfa(1)=-c(icrm,1)*b
            beta(1)=fff(1)*b
          end if
          do k=2,nzm-1
            e=1._8/(eign*rho(icrm,k)-a(icrm,k)-c(icrm,k)+a(icrm,k)*alfa(k-1))
            alfa(k)=-c(icrm,k)*e
            beta(k)=(fff(k)-a(icrm,k)*beta(k-1))*e
          end do

          fff(nzm)=(fff(nzm)-a(icrm,nzm)*beta(nzm-1))/ &
          (eign*rho(icrm,nzm)-a(icrm,nzm)+a(icrm,nzm)*alfa(nzm-1))

          do k=nzm-1,1,-1
            fff(k)=alfa(k)*fff(k+1)+beta(k)
          end do
          ff(icrm,i,j,1:nzm) = fff(1:nzm)

        enddo
      end do
    end do

    !-----------------------------------------------------------------
    !   Send the Fourier coefficient to the tasks performing
    !   the inverse Fourier transformation:
    n = rank*nzslab
    !$acc parallel loop gang vector collapse(4) default(present) async(1)
    do k = 1,nzslab
      do j = 1,nyp22-jwall
        do i = 1,nxp1-iwall
          do icrm = 1 , ncrms
            f(icrm,i+it,j+jt,k) = ff(icrm,i,j,k+n)
          end do
        end do
      end do
    end do

    !$acc update host(f) async(1)
    !$acc wait(1)

    !-------------------------------------------------
    !   Perform inverse Fourier transformation:

    if(rank.lt.npressureslabs) then
      do k=1,nzslab
        do icrm = 1 , ncrms
          if(RUN3D) call fft991_crm(f(icrm,:,:,k),work,trigxj,ifaxj,nx2,1,ny_gl,nx_gl+1,+1)
          call fft991_crm(f(icrm,:,:,k),work,trigxi,ifaxi,1,nx2,nx_gl,ny_gl,+1)
        end do
      end do
    endif

    !$acc update device(f) async(1)

    !-----------------------------------------------------------------
    !   Fill the pressure field for each subdomain:

    !$acc parallel loop gang vector default(present) async(1)
    do i=0,nx_gl
      iii(i)=i
      if (i == 0) iii(0) = nx_gl
    end do
    !$acc parallel loop gang vector default(present) async(1)
    do j=0,ny_gl
      jjj(j)=j
      if (j == 0) jjj(0) = ny_gl
    end do

    n = rank*nzslab
    !$acc parallel loop gang vector collapse(4) default(present) async(1)
    do k = 1,nzslab
      do j = 1-YES3D,ny
        do i = 0,nx
          do icrm = 1 , ncrms
            jj=jjj(j+jt)
            ii=iii(i+it)
            p(icrm,i,j,k+n) = f(icrm,ii,jj,k)
          end do
        end do
      end do
    end do

    !  Add pressure gradient term to the rhs of the momentum equation:
    call press_grad(ncrms)

    !$acc exit data delete(f,ff,work,trigxi,trigxj,ifaxj,ifaxi,a,c,fff,alfa,beta,reqs_in,iii,jjj,flag) async(1)

    deallocate( f       )
    deallocate( ff      )
    deallocate( work    )
    deallocate( trigxi  )
    deallocate( trigxj  )
    deallocate( ifaxj   )
    deallocate( ifaxi   )
    deallocate( a       )
    deallocate( c       )
    deallocate( fff     )
    deallocate( alfa    )
    deallocate( beta    )
    deallocate( reqs_in )
    deallocate( iii     )
    deallocate( jjj     )
    deallocate( flag    )

  end subroutine pressure

end module pressure_mod
