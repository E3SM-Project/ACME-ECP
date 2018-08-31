module sgs

  ! module for original SAM subgrid-scale SGS closure (Smagorinsky or 1st-order TKE)
  ! Marat Khairoutdinov, 2012

  use grid, only: nx,nxp1,ny,nyp1,YES3D,nzm,nz,dimx1_s,dimx2_s,dimy1_s,dimy2_s
  use params, only: dosgs, crm_rknd
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

  real(crm_rknd), allocatable, target :: sgs_field     (:,:,:,:)
  real(crm_rknd), allocatable, target :: sgs_field_diag(:,:,:,:)
  real(crm_rknd), allocatable :: fluxbsgs (:,:,:) ! surface fluxes
  real(crm_rknd), allocatable :: fluxtsgs (:,:,:) ! top boundary fluxes
  real(crm_rknd), allocatable :: sgswle   (:,:)  ! resolved vertical flux
  real(crm_rknd), allocatable :: sgswsb   (:,:)  ! SGS vertical flux
  real(crm_rknd), allocatable :: sgsadv   (:,:)  ! tendency due to vertical advection
  real(crm_rknd), allocatable :: sgslsadv (:,:)  ! tendency due to large-scale vertical advection
  real(crm_rknd), allocatable :: sgsdiff  (:,:)  ! tendency due to vertical diffusion
  real(crm_rknd), allocatable :: grdf_x(:)! grid factor for eddy diffusion in x
  real(crm_rknd), allocatable :: grdf_y(:)! grid factor for eddy diffusion in y
  real(crm_rknd), allocatable :: grdf_z(:)! grid factor for eddy diffusion in z
  real(crm_rknd), allocatable :: tkesbbuoy (:)
  real(crm_rknd), allocatable :: tkesbshear(:)
  real(crm_rknd), allocatable :: tkesbdiss (:)
  real(crm_rknd), allocatable :: tkesbdiff (:)
  real(crm_rknd), pointer :: tke (:,:,:)   ! SGS TKE
  real(crm_rknd), pointer :: tk  (:,:,:) ! SGS eddy viscosity
  real(crm_rknd), pointer :: tkh (:,:,:) ! SGS eddy conductivity

CONTAINS


  subroutine allocate_sgs()
    implicit none
    real(crm_rknd) :: zero
    allocate( sgs_field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, nsgs_fields) )
    allocate( sgs_field_diag(dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm, nsgs_fields_diag) )
    allocate( fluxbsgs (nx,ny,1:nsgs_fields)  )
    allocate( fluxtsgs (nx,ny,1:nsgs_fields)  )
    allocate( sgswle(nz,1:nsgs_fields)   )
    allocate( sgswsb(nz,1:nsgs_fields)   )
    allocate( sgsadv(nz,1:nsgs_fields)   )
    allocate( sgslsadv(nz,1:nsgs_fields)   )
    allocate( sgsdiff(nz,1:nsgs_fields)   )
    allocate( grdf_x(nzm) )
    allocate( grdf_y(nzm) )
    allocate( grdf_z(nzm) )
    allocate( tkesbbuoy(nz) )
    allocate( tkesbshear(nz) )
    allocate( tkesbdiss(nz) )
    allocate( tkesbdiff(nz) )

    tke(dimx1_s:,dimy1_s:,1:) => sgs_field     (:,:,:,1)
    tk (dimx1_d:,dimy1_d:,1:) => sgs_field_diag(:,:,:,1)
    tkh(dimx1_d:,dimy1_d:,1:) => sgs_field_diag(:,:,:,2)

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
    tke = zero
    tk   = zero
    tkh  = zero
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

      sgs_field = 0.
      sgs_field_diag = 0.

      fluxbsgs = 0.
      fluxtsgs = 0.

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
        grdf_x(k) = dx**2/(adz(k)*dz(icrm))**2
        grdf_y(k) = dy**2/(adz(k)*dz(icrm))**2
        grdf_z(k) = 1.
      end do
    else
      do k=1,nzm
        grdf_x(k) = min( real(16.,crm_rknd), dx**2/(adz(k)*dz(icrm))**2)
        grdf_y(k) = min( real(16.,crm_rknd), dy**2/(adz(k)*dz(icrm))**2)
        grdf_z(k) = 1.
      end do
    end if

    sgswle = 0.
    sgswsb = 0.
    sgsadv = 0.
    sgsdiff = 0.
    sgslsadv = 0.


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
              tke(i,j,k)=0.04*(5-k)
            endif
          end do
        end do
      end do

    case(1)

      do k=1,nzm
        do j=1,ny
          do i=1,nx
            if(q0(k,icrm).gt.6.e-3.and..not.dosmagor) then
              tke(i,j,k)=1.
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
              tke(i,j,k)=1.
            endif
          end do
        end do
      end do


    case(4)  ! gcss wg1 arm case

      do k=1,nzm
        do j=1,ny
          do i=1,nx
            if(z(k).le.150..and..not.dosmagor) then
              tke(i,j,k)=0.15*(1.-z(k)/150.)
            endif
          end do
        end do
      end do


    case(5)  ! gcss wg1 BOMEX case

      do k=1,nzm
        do j=1,ny
          do i=1,nx
            if(z(k).le.3000..and..not.dosmagor) then
              tke(i,j,k)=1.-z(k)/3000.
            endif
          end do
        end do
      end do

    case(6)  ! GCSS Lagragngian ASTEX


      do k=1,nzm
        do j=1,ny
          do i=1,nx
            if(q0(k,icrm).gt.6.e-3.and..not.dosmagor) then
              tke(i,j,k)=1.
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

  subroutine kurant_sgs(ncrms,icrm,cfl)

    use grid, only: dt, dx, dy, dz, adz, adzw
    implicit none
    integer, intent(in) :: ncrms,icrm
    real(crm_rknd), intent(out) :: cfl

    integer k
    real(crm_rknd) tkhmax(nz)

    do k = 1,nzm
      tkhmax(k) = maxval(tkh(1:nx,1:ny,k))
    end do

    cfl = 0.
    do k=1,nzm
      cfl = max(cfl,        &
      0.5*tkhmax(k)*grdf_z(k)*dt/(dz(icrm)*adzw(k))**2, &
      0.5*tkhmax(k)*grdf_x(k)*dt/dx**2, &
      YES3D*0.5*tkhmax(k)*grdf_y(k)*dt/dy**2)
    end do

  end subroutine kurant_sgs


  !----------------------------------------------------------------------
  !!! compute sgs diffusion of momentum:
  !
  subroutine sgs_mom(ncrms,icrm)
    use diffuse_mom_mod, only: diffuse_mom
    implicit none
    integer, intent(in) :: ncrms,icrm

    call diffuse_mom(ncrms,icrm,grdf_x, grdf_y, grdf_z, dimx1_d, dimx2_d, dimy1_d, dimy2_d, tk)

  end subroutine sgs_mom

  !----------------------------------------------------------------------
  !!! compute sgs diffusion of scalars:
  !
  subroutine sgs_scalars(ncrms,icrm)
    use diffuse_scalar_mod, only: diffuse_scalar
    use vars
    use microphysics
    use crmtracers
    use scalar_momentum_mod
    use params, only: dotracers
    implicit none
    integer, intent(in) :: ncrms,icrm
    real(crm_rknd) dummy(nz)
    real(crm_rknd) fluxbtmp(nx,ny), fluxttmp(nx,ny) !bloss
    integer k


    call diffuse_scalar(ncrms,icrm,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_y,grdf_z,tkh,t,fluxbt,fluxtt,tdiff(:,icrm),twsb(:,icrm), &
    t2lediff(:,icrm),t2lediss(:,icrm),twlediff(:,icrm),.true.)

    if(advect_sgs) then
      call diffuse_scalar(ncrms,icrm,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_y,grdf_z,tkh,tke,fzero,fzero,dummy,sgswsb, &
      dummy,dummy,dummy,.false.)
    end if


    !
    !    diffusion of microphysics prognostics:
    !
    call micro_flux(ncrms,icrm)

    total_water_evap(icrm) = total_water_evap(icrm) - total_water(ncrms,icrm)

    do k = 1,nmicro_fields
      if(   k.eq.index_water_vapor             &! transport water-vapor variable no metter what
      .or. docloud.and.flag_precip(k).ne.1    & ! transport non-precipitation vars
      .or. doprecip.and.flag_precip(k).eq.1 ) then
      fluxbtmp(1:nx,1:ny) = fluxbmk(1:nx,1:ny,k)
      fluxttmp(1:nx,1:ny) = fluxtmk(1:nx,1:ny,k)
      call diffuse_scalar(ncrms,icrm,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_y,grdf_z,tkh,micro_field(:,:,:,k),fluxbtmp,fluxttmp, &
      mkdiff(:,k),mkwsb(:,k), dummy,dummy,dummy,.false.)
    end if
  end do

  total_water_evap(icrm) = total_water_evap(icrm) + total_water(ncrms,icrm)

  ! diffusion of tracers:

  if(dotracers) then

    call tracers_flux()

    do k = 1,ntracers

      fluxbtmp = fluxbtr(:,:,k)
      fluxttmp = fluxttr(:,:,k)
      call diffuse_scalar(ncrms,icrm,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_y,grdf_z,tkh,tracer(:,:,:,k),fluxbtmp,fluxttmp, &
      trdiff(:,k),trwsb(:,k), &
      dummy,dummy,dummy,.false.)
      !!$          call diffuse_scalar(ncrms,icrm,tracer(:,:,:,k),fluxbtr(:,:,k),fluxttr(:,:,k),trdiff(:,k),trwsb(:,k), &
      !!$                           dummy,dummy,dummy,.false.)

    end do

  end if


#if defined(SP_ESMT)

    ! diffusion of scalar momentum tracers

    call diffuse_scalar(ncrms,icrm,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_y,grdf_z,tkh,   &
                        u_esmt,fluxb_u_esmt,fluxt_u_esmt,u_esmt_diff,u_esmt_sgs,    &
                        dummy,dummy,dummy,.false.)

    call diffuse_scalar(ncrms,icrm,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_y,grdf_z,tkh,   &
                        v_esmt,fluxb_v_esmt,fluxt_v_esmt,v_esmt_diff,v_esmt_sgs,    &
                        dummy,dummy,dummy,.false.)

#endif



end subroutine sgs_scalars

!----------------------------------------------------------------------
!!! compute sgs processes (beyond advection):
!
subroutine sgs_proc(ncrms,icrm)
  use tke_full_mod, only: tke_full
  use grid, only: nstep,dt,icycle
  use params, only: dosmoke
  implicit none
  integer, intent(in) :: ncrms,icrm
  !    SGS TKE equation:

  if(dosgs) call tke_full(ncrms,icrm,dimx1_d, dimx2_d, dimy1_d, dimy2_d, &
                          grdf_x, grdf_y, grdf_z, dosmagor,   &
                          tkesbdiss, tkesbshear, tkesbbuoy,   &
                          tke, tk, tkh)

  tke2(:,:,:,icrm) = tke
  tk2(:,:,:,icrm) = tk

end subroutine sgs_proc

!----------------------------------------------------------------------
!!! Diagnose arrays nessesary for dynamical core and statistics:
!
subroutine sgs_diagnose()
  ! None

end subroutine sgs_diagnose

!----------------------------------------------------------------------
! called when stepout() called

subroutine sgs_print()
  use utils, only: fminmax_print

  call fminmax_print('tke:',tke,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
  call fminmax_print('tk:',tk,0,nxp1,1-YES3D,nyp1,nzm)
  call fminmax_print('tkh:',tkh,0,nxp1,1-YES3D,nyp1,nzm)

end subroutine sgs_print

!----------------------------------------------------------------------
!!! Initialize the list of sgs statistics
!
subroutine sgs_hbuf_init(namelist,deflist,unitlist,status,average_type,count,sgscount)
  character(*) namelist(*), deflist(*), unitlist(*)
  integer status(*),average_type(*),count,sgscount

end subroutine sgs_hbuf_init


end module sgs
