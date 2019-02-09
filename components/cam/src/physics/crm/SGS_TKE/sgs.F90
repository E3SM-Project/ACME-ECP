
module sgs

  ! module for original SAM subgrid-scale SGS closure (Smagorinsky or 1st-order TKE)
  ! Marat Khairoutdinov, 2012

  use grid, only: nx,nxp1,ny,nyp1,YES3D,nzm,nz,dimx1_s,dimx2_s,dimy1_s,dimy2_s
  use params, only: dosgs, crm_rknd, asyncid
  use vars, only: tke2, tk2
  implicit none

  !----------------------------------------------------------------------
  ! Required definitions:

  !!! prognostic scalar (need to be advected arround the grid):

  integer, parameter :: nsgs_fields = 1   ! total number of prognostic sgs vars


  !!! sgs diagnostic variables that need to exchange boundary information (via MPI):

  integer, parameter :: nsgs_fields_diag = 2   ! total number of diagnostic sgs vars

  ! diagnostic fields' boundaries:
  integer, parameter :: dimx1_d=0, dimx2_d=nxp1, dimy1_d=1-YES3D, dimy2_d=nyp1

  integer, parameter :: flag_sgs3Dout(nsgs_fields) = (/0/)
  integer, parameter :: flag_sgsdiag3Dout(nsgs_fields_diag) = (/0,0/)


  logical:: advect_sgs = .false. ! advect prognostics or not, default - not (Smagorinsky)
  logical, parameter:: do_sgsdiag_bound = .true.  ! exchange boundaries for diagnostics fields

  ! SGS fields that output by default (if =1).

  !!! these arrays may be needed for output statistics:


  !------------------------------------------------------------------
  ! internal (optional) definitions:

  ! make aliases for prognostic variables:


  ! make aliases for diagnostic variables:


  logical:: dosmagor   ! if true, then use Smagorinsky closure

  ! whannah
  ! logical:: doscalar   ! if true, transport a passive scalar in the place of prognostic SGS TKE only if dosmagor=.true.

  ! Local diagnostics:

  real(crm_rknd), allocatable, target :: sgs_field     (:,:,:,:,:)
  real(crm_rknd), allocatable, target :: sgs_field_diag(:,:,:,:,:)
  real(crm_rknd), allocatable :: fluxbsgs (:,:,:,:) ! surface fluxes
  real(crm_rknd), allocatable :: fluxtsgs (:,:,:,:) ! top boundary fluxes
  real(crm_rknd), allocatable :: sgswle   (:,:,:)  ! resolved vertical flux
  real(crm_rknd), allocatable :: sgswsb   (:,:,:)  ! SGS vertical flux
  real(crm_rknd), allocatable :: sgsadv   (:,:,:)  ! tendency due to vertical advection
  real(crm_rknd), allocatable :: sgslsadv (:,:,:)  ! tendency due to large-scale vertical advection
  real(crm_rknd), allocatable :: sgsdiff  (:,:,:)  ! tendency due to vertical diffusion
  real(crm_rknd), allocatable :: grdf_x(:,:)! grid factor for eddy diffusion in x
  real(crm_rknd), allocatable :: grdf_y(:,:)! grid factor for eddy diffusion in y
  real(crm_rknd), allocatable :: grdf_z(:,:)! grid factor for eddy diffusion in z
  real(crm_rknd), allocatable :: tkesbbuoy (:,:)
  real(crm_rknd), allocatable :: tkesbshear(:,:)
  real(crm_rknd), allocatable :: tkesbdiss (:,:)
  real(crm_rknd), allocatable :: tkesbdiff (:,:)
  real(crm_rknd), pointer :: tke (:,:,:,:)   ! SGS TKE
  real(crm_rknd), pointer :: tk  (:,:,:,:) ! SGS eddy viscosity
  real(crm_rknd), pointer :: tkh (:,:,:,:) ! SGS eddy conductivity

CONTAINS


  subroutine allocate_sgs(ncrms)
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) :: zero
    allocate( sgs_field(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, nsgs_fields) )
    allocate( sgs_field_diag(ncrms,dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm, nsgs_fields_diag) )
    allocate( fluxbsgs (nx,ny,1:nsgs_fields,ncrms)  )
    allocate( fluxtsgs (nx,ny,1:nsgs_fields,ncrms)  )
    allocate( sgswle(ncrms,nz,1:nsgs_fields)   )
    allocate( sgswsb(nz,1:nsgs_fields,ncrms)   )
    allocate( sgsadv(ncrms,nz,1:nsgs_fields)   )
    allocate( sgslsadv(nz,1:nsgs_fields,ncrms)   )
    allocate( sgsdiff(nz,1:nsgs_fields,ncrms)   )
    allocate( grdf_x(ncrms,nzm) )
    allocate( grdf_y(ncrms,nzm) )
    allocate( grdf_z(ncrms,nzm) )
    allocate( tkesbbuoy(ncrms,nz) )
    allocate( tkesbshear(ncrms,nz) )
    allocate( tkesbdiss(ncrms,nz) )
    allocate( tkesbdiff(nz,ncrms) )

    tke(1:,dimx1_s:,dimy1_s:,1:) => sgs_field     (:,:,:,:,1)
    tk (1:,dimx1_d:,dimy1_d:,1:) => sgs_field_diag(:,:,:,:,1)
    tkh(1:,dimx1_d:,dimy1_d:,1:) => sgs_field_diag(:,:,:,:,2)

    zero = 0

    sgs_field = zero
    sgs_field_diag = zero
    fluxbsgs  = zero
    fluxtsgs  = zero
    sgswle = zero
    sgswsb = zero
    sgsadv = zero
    sgslsadv = zero
    sgsdiff = zero
    grdf_x = zero
    grdf_y = zero
    grdf_z = zero
    tkesbbuoy = zero
    tkesbshear = zero
    tkesbdiss = zero
    tkesbdiff = zero
  end subroutine allocate_sgs


  subroutine deallocate_sgs()
    implicit none
    deallocate( sgs_field  )
    deallocate( sgs_field_diag  )
    deallocate( fluxbsgs   )
    deallocate( fluxtsgs   )
    deallocate( sgswle  )
    deallocate( sgswsb  )
    deallocate( sgsadv  )
    deallocate( sgslsadv  )
    deallocate( sgsdiff  )
    deallocate( grdf_x  )
    deallocate( grdf_y  )
    deallocate( grdf_z  )
    deallocate( tkesbbuoy  )
    deallocate( tkesbshear  )
    deallocate( tkesbdiss  )
    deallocate( tkesbdiff  )
    nullify(tke)
    nullify(tk )
    nullify(tkh)
  end subroutine deallocate_sgs


  ! required microphysics subroutines and function:
  !----------------------------------------------------------------------
  !!! Read microphysics options from prm (namelist) file

  subroutine sgs_setparm()

    use grid, only: case
    implicit none

    integer ierr, ios, ios_missing_namelist, place_holder

    !======================================================================
    NAMELIST /SGS_TKE/ &
    dosmagor ! Diagnostic Smagorinsky closure

    NAMELIST /BNCUIODSBJCB/ place_holder

    dosmagor = .true.  ! default
    ! doscalar = .false. ! default ! whannah

    !----------------------------------
    !  Read namelist for microphysics options from prm file:
    !------------
    !open(55,file='./'//trim(case)//'/prm', status='old',form='formatted')

    !read (UNIT=55,NML=BNCUIODSBJCB,IOSTAT=ios_missing_namelist)
    !rewind(55) !note that one must rewind before searching for new namelists

    !read (55,SGS_TKE,IOSTAT=ios)

    advect_sgs = .not.dosmagor

    !if (ios.ne.0) then
    !   !namelist error checking
    !   if(ios.ne.ios_missing_namelist) then
    !      write(*,*) '****** ERROR: bad specification in SGS_TKE namelist'
    !      call task_abort()
    !   end if
    !end if
    !close(55)

    ! END UW ADDITION
    !======================================================================

  end subroutine sgs_setparm

  !----------------------------------------------------------------------
  !!! Initialize sgs:


  subroutine sgs_init(ncrms,icrm)
    use grid, only: nrestart, dx, dy, dz, adz, masterproc
    use params, only: LES
    implicit none
    integer, intent(in) :: ncrms,icrm
    integer k

    if(nrestart.eq.0) then

      sgs_field(icrm,:,:,:,:) = 0.
      sgs_field_diag(icrm,:,:,:,:) = 0.

      fluxbsgs(:,:,:,icrm) = 0.
      fluxtsgs(:,:,:,icrm) = 0.

    end if

    !  if(masterproc) then
    !     if(dosmagor) then
    !        write(*,*) 'Smagorinsky SGS Closure'
    !     else
    !        write(*,*) 'Prognostic TKE 1.5-order SGS Closure'
    !     end if
    !  end if

    if(LES) then
      do k=1,nzm
        grdf_x(icrm,k) = dx**2/(adz(icrm,k)*dz(icrm))**2
        grdf_y(icrm,k) = dy**2/(adz(icrm,k)*dz(icrm))**2
        grdf_z(icrm,k) = 1.
      end do
    else
      do k=1,nzm
        grdf_x(icrm,k) = min( real(16.,crm_rknd), dx**2/(adz(icrm,k)*dz(icrm))**2)
        grdf_y(icrm,k) = min( real(16.,crm_rknd), dy**2/(adz(icrm,k)*dz(icrm))**2)
        grdf_z(icrm,k) = 1.
      end do
    end if

    sgswle(icrm,:,:) = 0.
    sgswsb  (:,:,icrm) = 0.
    sgsadv(icrm,:,:) = 0.
    sgsdiff (:,:,icrm) = 0.
    sgslsadv(:,:,icrm) = 0.


  end subroutine sgs_init

  !----------------------------------------------------------------------
  !!! make some initial noise in sgs:
  !
  subroutine setperturb_sgs(ncrms,icrm,ptype)

    use vars, only: q0, z
    integer, intent(in) :: ncrms,icrm
    integer, intent(in) :: ptype
    integer i,j,k

    select case (ptype)

    case(0)

      do k=1,nzm
        do j=1,ny
          do i=1,nx
            if(k.le.4.and..not.dosmagor) then
              tke(icrm,i,j,k)=0.04*(5-k)
            endif
          end do
        end do
      end do

    case(1)

      do k=1,nzm
        do j=1,ny
          do i=1,nx
            if(q0(k,icrm).gt.6.e-3.and..not.dosmagor) then
              tke(icrm,i,j,k)=1.
            endif
          end do
        end do
      end do

    case(2)

    case(3)   ! gcss wg1 smoke-cloud case

      do k=1,nzm
        do j=1,ny
          do i=1,nx
            if(q0(k,icrm).gt.0.5e-3.and..not.dosmagor) then
              tke(icrm,i,j,k)=1.
            endif
          end do
        end do
      end do


    case(4)  ! gcss wg1 arm case

      do k=1,nzm
        do j=1,ny
          do i=1,nx
            if(z(icrm,k).le.150..and..not.dosmagor) then
              tke(icrm,i,j,k)=0.15*(1.-z(icrm,k)/150.)
            endif
          end do
        end do
      end do


    case(5)  ! gcss wg1 BOMEX case

      do k=1,nzm
        do j=1,ny
          do i=1,nx
            if(z(icrm,k).le.3000..and..not.dosmagor) then
              tke(icrm,i,j,k)=1.-z(icrm,k)/3000.
            endif
          end do
        end do
      end do

    case(6)  ! GCSS Lagragngian ASTEX


      do k=1,nzm
        do j=1,ny
          do i=1,nx
            if(q0(k,icrm).gt.6.e-3.and..not.dosmagor) then
              tke(icrm,i,j,k)=1.
            endif
          end do
        end do
      end do


    case default

    end select

  end subroutine setperturb_sgs

  !----------------------------------------------------------------------
  !!! Estimate Courant number limit for SGS
  !

  subroutine kurant_sgs(ncrms,cfl)
    use grid, only: dt, dx, dy, dz, adz, adzw
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd), intent(inout) :: cfl
    integer k,icrm, j, i
    real(crm_rknd) tkhmax(ncrms,nz), tmp

    !$acc enter data create(tkhmax) async(asyncid)

    !$acc parallel loop collapse(2) copyout(tkhmax) async(asyncid)
    do k = 1,nzm
      do icrm = 1 , ncrms
        tkhmax(icrm,k) = 0.
      enddo
    enddo

    !$acc parallel loop collapse(4) copy(tkhmax) copyin(tkh) async(asyncid)
    do k = 1,nzm
      do j = 1 , ny
        do i = 1 , nx
          do icrm = 1 , ncrms
            !$acc atomic update
            tkhmax(icrm,k) = max(tkhmax(icrm,k),tkh(icrm,i,j,k))
          enddo
        enddo
      end do
    end do

    !$acc parallel loop collapse(2) private(tmp) copy(cfl) copyin(tkhmax,grdf_x,grdf_y,grdf_z,dz,adzw) async(asyncid)
    do k=1,nzm
      do icrm = 1 , ncrms
        tmp = max( 0.5*tkhmax(icrm,k)*grdf_z(icrm,k)*dt/(dz(icrm)*adzw(icrm,k))**2  , &
                   0.5*tkhmax(icrm,k)*grdf_x(icrm,k)*dt/dx**2  , &
                   YES3D*0.5*tkhmax(icrm,k)*grdf_y(icrm,k)*dt/dy**2  )
        !$acc atomic update
        cfl = max( cfl , tmp )
      end do
    end do

    !$acc exit data delete(tkhmax) async(asyncid)

  end subroutine kurant_sgs


  !----------------------------------------------------------------------
  !!! compute sgs diffusion of momentum:
  !
  subroutine sgs_mom(ncrms)
    use diffuse_mom_mod, only: diffuse_mom
    implicit none
    integer, intent(in) :: ncrms

#ifdef __PGI
    !Passing tk via first element to avoid PGI pointer bug
    call diffuse_mom(ncrms,grdf_x, grdf_y, grdf_z, dimx1_d, dimx2_d, dimy1_d, dimy2_d, tk(1,dimx1_d,dimy1_d,1))
#else
    call diffuse_mom(ncrms,grdf_x, grdf_y, grdf_z, dimx1_d, dimx2_d, dimy1_d, dimy2_d, tk)
#endif
  end subroutine sgs_mom

  !----------------------------------------------------------------------
  !!! compute sgs diffusion of scalars:
  !
  subroutine sgs_scalars(ncrms)
    use diffuse_scalar_mod, only: diffuse_scalar
    use vars
    use microphysics
    use crmtracers
    use scalar_momentum_mod
    use params, only: dotracers
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) dummy(nz,ncrms)
    real(crm_rknd) fluxbtmp(nx,ny,ncrms), fluxttmp(nx,ny,ncrms), difftmp(nz,ncrms), wsbtmp(nz,ncrms)
    integer i,j,kk,k,icrm

    !$acc enter data create(dummy,fluxbtmp,fluxttmp,difftmp,wsbtmp) async(asyncid)
    
#ifdef __PGI
    !Passing tkh via first element to avoid PGI pointer bug
    call diffuse_scalar(ncrms,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_y,grdf_z,tkh(1,dimx1_d,dimy1_d,1),t,fluxbt,fluxtt,tdiff,twsb)
#else
    call diffuse_scalar(ncrms,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_y,grdf_z,tkh,t,fluxbt,fluxtt,tdiff,twsb)
#endif

    if(advect_sgs) then
      !$acc parallel loop collapse(2) copyin(sgswsb) copy(wsbtmp) async(asyncid)
      do icrm = 1, ncrms
        do k = 1 , nz
          wsbtmp(k,icrm) = sgswsb(k,1,icrm)
        enddo
      enddo
#ifdef __PGI
      !Passing tkh via first element to avoid PGI pointer bug
      call diffuse_scalar(ncrms,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_y,grdf_z,tkh(1,dimx1_d,dimy1_d,1),tke,fzero,fzero,dummy,wsbtmp)
#else
      call diffuse_scalar(ncrms,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_y,grdf_z,tkh,tke,fzero,fzero,dummy,wsbtmp)
#endif
      !$acc parallel loop collapse(2) copyin(wsbtmp) copy(sgswsb) async(asyncid)
      do icrm = 1, ncrms
        do k = 1 , nz
          sgswsb(k,1,icrm) = wsbtmp(k,icrm)
        enddo
      enddo
    end if

    !    diffusion of microphysics prognostics:
    call micro_flux(ncrms)
    !do icrm = 1, ncrms
    !  total_water_evap(icrm) = total_water_evap(icrm) - total_water(ncrms,icrm)
    !enddo

    do k = 1,nmicro_fields
      if(   k.eq.index_water_vapor             &! transport water-vapor variable no metter what
      .or. docloud.and.flag_precip(k).ne.1    & ! transport non-precipitation vars
      .or. doprecip.and.flag_precip(k).eq.1 ) then
        !$acc parallel loop collapse(2) copyin(fluxbmk,fluxtmk) copy(fluxbtmp,fluxttmp) async(asyncid)
        do icrm = 1 , ncrms
          do j = 1 , ny
            do i = 1 , nx
              fluxbtmp(i,j,icrm) = fluxbmk(i,j,k,icrm)
              fluxttmp(i,j,icrm) = fluxtmk(i,j,k,icrm)
            enddo
          enddo
        enddo
        !$acc parallel loop collapse(2) copyin(mkdiff,mkwsb) copy(difftmp,wsbtmp) async(asyncid)
        do icrm = 1 , ncrms
          do kk = 1 , nz
            difftmp(kk,icrm) = mkdiff(kk,k,icrm)
            wsbtmp (kk,icrm) = mkwsb (kk,k,icrm)
          enddo
        enddo
#ifdef __PGI
        !Passing tkh via first element to avoid PGI pointer bug
        call diffuse_scalar(ncrms,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_y,grdf_z,tkh(1,dimx1_d,dimy1_d,1),micro_field(:,:,:,:,k),fluxbtmp,fluxttmp,difftmp,wsbtmp)
#else
        call diffuse_scalar(ncrms,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_y,grdf_z,tkh,micro_field(:,:,:,:,k),fluxbtmp,fluxttmp,difftmp,wsbtmp)
#endif
        !$acc parallel loop collapse(2) copyin(difftmp,wsbtmp) copy(mkdiff,mkwsb) async(asyncid)
        do icrm = 1 , ncrms
          do kk = 1 , nz
            mkdiff(kk,k,icrm) = difftmp(kk,icrm)
            mkwsb (kk,k,icrm) = wsbtmp (kk,icrm)
          enddo
        enddo
      end if
    end do

    !$acc exit data delete(dummy,fluxbtmp,fluxttmp,difftmp,wsbtmp) async(asyncid)

    !if(dotracers) then
    !  call tracers_flux()
    !  do k = 1,ntracers
    !    fluxbtmp(1:nx,1:ny,icrm) = fluxbtr(:,:,k,icrm)
    !    fluxttmp(1:nx,1:ny,icrm) = fluxttr(:,:,k,icrm)
    !    call diffuse_scalar(ncrms,icrm,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_y,grdf_z,tkh,tracer(:,:,:,k,icrm),fluxbtmp(:,:,icrm),fluxttmp(:,:,icrm), &
    !    trdiff(:,k,icrm),trwsb(:,k,icrm), &
    !    dummy,dummy,dummy,.false.)
    !    !!$          call diffuse_scalar(ncrms,icrm,tracer(:,:,:,k,icrm),fluxbtr(:,:,k,icrm),fluxttr(:,:,k,icrm),trdiff(:,k,icrm),trwsb(:,k,icrm), &
    !    !!$                           dummy,dummy,dummy,.false.)
    !  end do
    !end if

    !do icrm = 1 , ncrms
    !  total_water_evap(icrm) = total_water_evap(icrm) + total_water(ncrms,icrm)
    !enddo

#if defined(SP_ESMT)
    ! diffusion of scalar momentum tracers
    !Passing tkh via first element to avoid PGI pointer bug
    call diffuse_scalar(ncrms,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_y,grdf_z,tkh(1,dimx1_d,dimy1_d,1),u_esmt,fluxb_u_esmt,fluxt_u_esmt,u_esmt_diff,u_esmt_sgs)
    !Passing tkh via first element to avoid PGI pointer bug
    call diffuse_scalar(ncrms,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_y,grdf_z,tkh(1,dimx1_d,dimy1_d,1),v_esmt,fluxb_v_esmt,fluxt_v_esmt,v_esmt_diff,v_esmt_sgs)
#endif
  end subroutine sgs_scalars

!----------------------------------------------------------------------
!!! compute sgs processes (beyond advection):
!
subroutine sgs_proc(ncrms)
  use tke_full_mod, only: tke_full
  use grid, only: dt,icycle
  use params, only: dosmoke
  implicit none
  integer, intent(in) :: ncrms
  integer :: icrm, k, j, i
  !    SGS TKE equation:

#ifdef __PGI
  !Passing tke, tk, and tkh via first element to avoid PGI pointer bug
  if(dosgs) call tke_full(ncrms,dimx1_d, dimx2_d, dimy1_d, dimy2_d, &
                          grdf_x, grdf_y, grdf_z, dosmagor,   &
                          tkesbdiss, tkesbshear, tkesbbuoy,   &
                          tke(1,dimx1_s,dimy1_s,1), tk(1,dimx1_d,dimy1_d,1), tkh(1,dimx1_d,dimy1_d,1))
#else
  if(dosgs) call tke_full(ncrms,dimx1_d, dimx2_d, dimy1_d, dimy2_d, &
                          grdf_x, grdf_y, grdf_z, dosmagor,   &
                          tkesbdiss, tkesbshear, tkesbbuoy,   &
                          tke, tk, tkh)
#endif
  !$acc parallel loop collapse(4) copyin(tke) copy(tke2) async(asyncid)
  do k = 1 , nzm
    do j = dimy1_s,dimy2_s
      do i = dimx1_s,dimx2_s
        do icrm = 1 , ncrms
          tke2(icrm,i,j,k) = tke(icrm,i,j,k)
        enddo
      enddo
    enddo
  enddo
  !$acc parallel loop collapse(4) copyin(tk) copy(tk2) async(asyncid)
  do k = 1 , nzm
    do j = dimy1_d,dimy2_d
      do i = dimx1_d,dimx2_d
        do icrm = 1 , ncrms
          tk2(icrm,i,j,k) = tk(icrm,i,j,k)
        enddo
      enddo
    enddo
  enddo

end subroutine sgs_proc

!----------------------------------------------------------------------
!!! Diagnose arrays nessesary for dynamical core and statistics:
!
subroutine sgs_diagnose()
  ! None

end subroutine sgs_diagnose

!----------------------------------------------------------------------
!!! Initialize the list of sgs statistics
!
subroutine sgs_hbuf_init(namelist,deflist,unitlist,status,average_type,count,sgscount)
  character(*) namelist(*), deflist(*), unitlist(*)
  integer status(*),average_type(*),count,sgscount

end subroutine sgs_hbuf_init


end module sgs
