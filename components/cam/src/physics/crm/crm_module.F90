
module crm_module
  use task_init_mod, only: task_init
  use abcoefs_mod, only: abcoefs
  use kurant_mod, only: kurant
  use setperturb_mod, only: setperturb
  use boundaries_mod, only: boundaries
  use forcing_mod, only: forcing
  use advect_mom_mod, only: advect_mom
  use adams_mod, only: adams
  use advect_all_scalars_mod, only: advect_all_scalars
  use sat_mod
  use crmsurface_mod
#ifdef sam1mom
  use precip_init_mod
#endif
  use zero_mod
  use buoyancy_mod
  use pressure_mod
  use uvw_mod
  use diagnose_mod
  use damping_mod
  use ice_fall_mod
  use coriolis_mod

  use crm_state_module,       only: crm_state_type
  use crm_rad_module,         only: crm_rad_type
  use crm_input_module,       only: crm_input_type
  use crm_output_module,      only: crm_output_type
  use crm_ecpp_output_module, only: crm_ecpp_output_type

!---------------------------------------------------------------
!  Super-parameterization's main driver
!  Marat Khairoutdinov, 2001-2009
!---------------------------------------------------------------
use setparm_mod, only : setparm

contains

subroutine crm(lchnk, icol, ncrms, phys_stage, dt_gl, plev, &
                crm_input, crm_state, crm_rad,  &
#ifdef CLUBB_CRM
                clubb_buffer,           &
                crm_cld, clubb_tk,      &
                clubb_tkh, relvar,      &
                accre_enhan, qclvar,    &
#endif
                crm_ecpp_output, crm_output )
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    use shr_kind_mod          , only: r8 => shr_kind_r8
    use phys_grid             , only: get_rlon_p, get_rlat_p, get_gcol_p  !, get_gcol_all_p
    use ppgrid                , only: pcols
    use vars
    use params
    use microphysics
    use sgs
    use crmtracers
    use scalar_momentum_mod
#ifdef MODAL_AERO
    use modal_aero_data       , only: ntot_amode
#endif
    use crmdims               , only: nclubbvars, crm_nx_rad, crm_ny_rad
#ifdef CLUBB_CRM
    use clubb_sgs             , only: advance_clubb_sgs, clubb_sgs_setup, clubb_sgs_cleanup, apply_clubb_sgs_tndcy, apply_clubb_sgs_tndcy_scalars, &
                                      apply_clubb_sgs_tndcy_mom, t2thetal, total_energy
    use clubb_precision       , only: time_precision, core_rknd
    use clubbvars             , only: up2, vp2, wprtp, wpthlp, wp2, wp3, rtp2, thlp2, rtpthlp, upwp, vpwp, cloud_frac, t_tndcy, qc_tndcy, qv_tndcy, &
                                      u_tndcy, v_tndcy, lrestart_clubb, rho_ds_zt, rho_ds_zm, thv_ds_zt, thv_ds_zm, invrs_rho_ds_zt, invrs_rho_ds_zm, &
                                      tracer_tndcy, sclrp2, sclrprtp, sclrpthlp, wpsclrp, relvarg, accre_enhang, qclvarg, edsclr_dim, sclr_dim, rho_ds_zt, &
                                      rho_ds_zm, rtm_spurious_source, thlm_spurious_source
    use fill_holes            , only: vertical_integral
    use numerical_check       , only: calculate_spurious_source
    use grid_class            , only: gr
#endif
#ifdef ECPP
    use ecppvars              , only: qlsink, precr, precsolid, &
                                      area_bnd_final, area_bnd_sum, area_cen_final, area_cen_sum, &
                                      mass_bnd_final, mass_bnd_sum, rh_cen_sum, qcloud_cen_sum, qice_cen_sum, &
                                      qlsink_cen_sum, precr_cen_sum, precsolid_cen_sum, xkhvsum, wup_thresh, wdown_thresh, &
                                      wwqui_cen_sum, wwqui_bnd_sum, wwqui_cloudy_cen_sum, wwqui_cloudy_bnd_sum, &
                                      qlsink_bf_cen_sum, qlsink_avg_cen_sum, prain_cen_sum, qlsink_bf, prain
    use module_ecpp_crm_driver, only: ecpp_crm_stat, ecpp_crm_init, ecpp_crm_cleanup, ntavg1_ss, ntavg2_ss
    use ecppvars              , only: NCLASS_CL, ncls_ecpp_in, NCLASS_PR
#endif /* ECPP */
    use cam_abortutils        , only: endrun
    use time_manager          , only: get_nstep

    implicit none

    !-----------------------------------------------------------------------------------------------
    ! Interface variable declarations
    !-----------------------------------------------------------------------------------------------

    integer , intent(in   ) :: lchnk                            ! chunk identifier (only for lat/lon and random seed)
    integer , intent(in   ) :: ncrms                            ! Number of "vector" GCM columns to push down into CRM for SIMD vectorization / more threading
    integer , intent(in   ) :: phys_stage                       ! physics run stage indicator (1 or 2 = bc or ac)
    integer , intent(in   ) :: plev                             ! number of levels in parent model
    real(r8), intent(in   ) :: dt_gl                            ! global model's time step
    integer , intent(in   ) :: icol                (ncrms)      ! column identifier (only for lat/lon and random seed)
    type(crm_input_type),      intent(in   ) :: crm_input
    type(crm_state_type),      intent(inout) :: crm_state
    type(crm_rad_type),        intent(inout) :: crm_rad
#ifdef CLUBB_CRM
    real(r8), intent(inout), target :: clubb_buffer(ncrms,crm_nx, crm_ny, crm_nz+1,1:nclubbvars)
    real(r8), intent(  out) :: crm_cld             (ncrms,crm_nx, crm_ny, crm_nz+1)
    real(r8), intent(  out) :: clubb_tk            (ncrms,crm_nx, crm_ny, crm_nz)
    real(r8), intent(  out) :: clubb_tkh           (ncrms,crm_nx, crm_ny, crm_nz)
    real(r8), intent(  out) :: relvar              (ncrms,crm_nx, crm_ny, crm_nz)
    real(r8), intent(  out) :: accre_enhan         (ncrms,crm_nx, crm_ny, crm_nz)
    real(r8), intent(  out) :: qclvar              (ncrms,crm_nx, crm_ny, crm_nz)
#endif /* CLUBB_CRM */
    type(crm_ecpp_output_type),intent(inout) :: crm_ecpp_output
    type(crm_output_type),     intent(inout) :: crm_output

    !-----------------------------------------------------------------------------------------------
    ! Local variable declarations
    !-----------------------------------------------------------------------------------------------

    real(r8),       parameter :: umax = 0.5*crm_dx/crm_dt       ! maxumum ampitude of the l.s. wind
    real(r8),       parameter :: wmin = 2.                      ! minimum up/downdraft velocity for stat
    real(crm_rknd), parameter :: cwp_threshold = 0.001          ! threshold for cloud condensate for shaded fraction calculation
    integer,        parameter :: perturb_seed_scale = 1000      ! scaling value for setperturb() seed value (seed = gcol * perturb_seed_scale)
    real(r8)        :: crm_run_time                             ! length of CRM integration (=dt_gl*0.5 if SP_CRM_SPLIT is defined)
    real(r8)        :: icrm_run_time                            ! = 1 / crm_run_time
    real(r8)        :: factor_xy, factor_xyt, idt_gl
    real(crm_rknd)  :: tmp1, tmp2, tmp
    real(crm_rknd)  :: u2z,v2z,w2z
    integer         :: i,j,k,l,ptop,nn,icyc,icrm
    integer         :: kx
    real(crm_rknd)  :: ustar(ncrms), bflx(ncrms), wnd(ncrms), qsat, omg
    real(crm_rknd)  :: colprec,colprecs
    real(r8)        :: qtot(ncrms,20)    ! Total water for water conservation check

    !!! These should all be inputs
    integer         :: igstep            ! GCM time steps
    integer         :: iseed             ! seed for random perturbation
    !!! variables for radiation grouping method
    real(crm_rknd) :: crm_nx_rad_fac
    real(crm_rknd) :: crm_ny_rad_fac
    integer        :: i_rad
    integer        :: j_rad

    !!! Arrays
    real(crm_rknd), allocatable :: t00(:,:)
    real(crm_rknd), allocatable :: fluxbtmp(:,:,:)
    real(crm_rknd), allocatable :: fluxttmp(:,:,:)    !bloss
    real(crm_rknd), allocatable :: tln  (:,:)
    real(crm_rknd), allocatable :: qln  (:,:)
    real(crm_rknd), allocatable :: qccln(:,:)
    real(crm_rknd), allocatable :: qiiln(:,:)
    real(crm_rknd), allocatable :: uln  (:,:)
    real(crm_rknd), allocatable :: vln  (:,:)
#if defined(SP_ESMT)
    real(crm_rknd), allocatable  :: uln_esmt(:,:)
    real(crm_rknd), allocatable  :: vln_esmt(:,:)     ! tempoerary variables for expliciit scalar momentum transport
#endif
    real(crm_rknd), allocatable  :: cwp     (:,:,:)
    real(crm_rknd), allocatable  :: cwph    (:,:,:)
    real(crm_rknd), allocatable  :: cwpm    (:,:,:)
    real(crm_rknd), allocatable  :: cwpl    (:,:,:)
    logical       , allocatable  :: flag_top(:,:,:)
    real(crm_rknd), allocatable  :: cltemp  (:,:,:)
    real(crm_rknd), allocatable  :: cmtemp  (:,:,:)
    real(crm_rknd), allocatable  :: chtemp  (:,:,:)
    real(crm_rknd), allocatable  :: cttemp  (:,:,:)
    integer       , allocatable  :: gcolindex(:,:)  ! array of global latitude indices
#ifdef CLUBB_CRM
    !Array indicies for spurious RTM check
    real(kind=core_rknd), allocatable :: thlm_flux_top, thlm_flux_sfc, rtm_flux_top, rtm_flux_sfc
    real(kind=core_rknd), allocatable :: rtm_integral_before (:,:)
    real(kind=core_rknd), allocatable :: rtm_integral_after  (:,:)
    real(kind=core_rknd), allocatable :: thlm_integral_before(:,:)
    real(kind=core_rknd), allocatable :: thlm_integral_after (:,:)
    real(kind=core_rknd), allocatable :: thlm_before(:)
    real(kind=core_rknd), allocatable :: thlm_after (:)
    real(kind=core_rknd), allocatable :: rtm_column (:) ! Total water (vapor + liquid)     [kg/kg]
#endif /* CLUBB_CRM */
    real(crm_rknd) :: zeroval

    real(r8), allocatable :: dd_crm (:,:)     ! mass entraiment from downdraft
    real(r8), allocatable :: mui_crm(:,:)     ! mass flux up at the interface
    real(r8), allocatable :: mdi_crm(:,:)     ! mass flux down at the interface

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------

    allocate( t00(nz,ncrms) )
    allocate( fluxbtmp(nx,ny,ncrms) )
    allocate( fluxttmp(nx,ny,ncrms) )
    allocate( tln(plev,ncrms) )
    allocate( qln(plev,ncrms) )
    allocate( qccln(plev,ncrms) )
    allocate( qiiln(plev,ncrms) )
    allocate( uln(plev,ncrms) )
    allocate( vln(plev,ncrms) )
#if defined(SP_ESMT)
    allocate( uln_esmt(plev,ncrms) )
    allocate( vln_esmt(plev,ncrms) )
#endif
    allocate( cwp(nx,ny,ncrms) )
    allocate( cwph(nx,ny,ncrms) )
    allocate( cwpm(nx,ny,ncrms) )
    allocate( cwpl(nx,ny,ncrms) )
    allocate( flag_top(nx,ny,ncrms) )
    allocate( gcolindex(pcols,ncrms) )
    allocate( cltemp(nx,ny,ncrms) )
    allocate( cmtemp(nx,ny,ncrms) )
    allocate( chtemp(nx,ny,ncrms) )
    allocate( cttemp(nx,ny,ncrms) )
#ifdef CLUBB_CRM
    allocate( rtm_integral_before (nx,ny) )
    allocate( rtm_integral_after (nx,ny) )
    allocate( thlm_integral_before(nx,ny) )
    allocate( thlm_integral_after(nx,ny) )
    allocate( thlm_before(nzm) )
    allocate( thlm_after(nzm) )
    allocate( rtm_column(nzm) )
#endif /* CLUBB_CRM */
    allocate( dd_crm (ncrms,plev)   )
    allocate( mui_crm(ncrms,plev+1) )
    allocate( mdi_crm(ncrms,plev+1) )

    zeroval = 0

    t00  = zeroval
    fluxbtmp  = zeroval
    fluxttmp  = zeroval
    tln  = zeroval
    qln  = zeroval
    qccln  = zeroval
    qiiln  = zeroval
    uln  = zeroval
    vln  = zeroval
#if defined(SP_ESMT)
    uln_esmt = zeroval
    vln_esmt = zeroval
#endif
    cwp = zeroval
    cwph = zeroval
    cwpm = zeroval
    cwpl = zeroval
    flag_top = .false.
    gcolindex = zeroval
    cltemp = zeroval
    cmtemp = zeroval
    chtemp = zeroval
    cttemp = zeroval
#ifdef CLUBB_CRM
    rtm_integral_before  = zeroval
    rtm_integral_after  = zeroval
    thlm_integral_before = zeroval
    thlm_integral_after = zeroval
    thlm_before = zeroval
    thlm_after = zeroval
    rtm_column = zeroval
#endif /* CLUBB_CRM */

  call allocate_params(ncrms)
  call allocate_vars(ncrms)
  call allocate_grid(ncrms)
  call allocate_tracers(ncrms)
  call allocate_sgs(ncrms)
  call allocate_micro(ncrms)
#ifdef sam1mom
  call allocate_micro_params(ncrms)
#endif
#if defined(SP_ESMT)
  call allocate_scalar_momentum(ncrms)
#endif

  !Loop over "vector columns"
  do icrm = 1 , ncrms
    latitude0 (icrm) = get_rlat_p(lchnk, icol(icrm)) * 57.296_r8
    longitude0(icrm) = get_rlon_p(lchnk, icol(icrm)) * 57.296_r8

    igstep = get_nstep()

!-----------------------------------------------

    dostatis  = .false.    ! no statistics are collected.
    idt_gl    = 1._r8/dt_gl
    ptop      = plev-nzm+1
    factor_xy = 1._r8/dble(nx*ny)
    crm_rad%temperature  (icrm,:,:,:) = 0.
    crm_rad%qv (icrm,:,:,:) = 0.
    crm_rad%qc (icrm,:,:,:) = 0.
    crm_rad%qi (icrm,:,:,:) = 0.
    crm_rad%cld(icrm,:,:,:) = 0.
#ifdef m2005
    crm_rad%nc(icrm,:,:,:) = 0.0
    crm_rad%ni(icrm,:,:,:) = 0.0
    crm_rad%qs(icrm,:,:,:) = 0.0
    crm_rad%ns(icrm,:,:,:) = 0.0
#endif /* m2005 */
    bflx(icrm) = crm_input%bflxls(icrm)
    wnd(icrm) = crm_input%wndls(icrm)

!-----------------------------------------

#ifdef CLUBB_CRM
    if(igstep == 1) then
      lrestart_clubb = .false.
    else
     lrestart_clubb = .true.
    endif
#endif /* CLUBB_CRM */

    call task_init ()
    call setparm()

    fcor(icrm)= 4*pi/86400.*sin(latitude0(icrm)*pi/180.)
    fcorz(icrm) = sqrt(4.*(2*pi/(3600.*24.))**2-fcor(icrm)**2)
    fcory(:,icrm) = fcor(icrm)
    fcorzy(:,icrm) = fcorz(icrm)
    do j=1,ny
      do i=1,nx
        latitude (i,j,icrm) = latitude0(icrm)
        longitude(i,j,icrm) = longitude0(icrm)
      end do
    end do

    if(crm_input%ocnfrac(icrm).gt.0.5) then
       OCEAN(icrm) = .true.
    else
       LAND(icrm) = .true.
    end if

    ! Create CRM vertical grid and initialize some vertical reference arrays:
    do k = 1, nzm
      z(k,icrm) = crm_input%zmid(icrm,plev-k+1) - crm_input%zint(icrm,plev+1)
      zi(k,icrm) = crm_input%zint(icrm,plev-k+2)- crm_input%zint(icrm,plev+1)
      pres(k,icrm) = crm_input%pmid(icrm,plev-k+1)/100.
      presi(k,icrm) = crm_input%pint(icrm,plev-k+2)/100.
      prespot(k,icrm)=(1000./pres(k,icrm))**(rgas/cp)
      bet(k,icrm) = ggr/crm_input%tl(icrm,plev-k+1)
      gamaz(k,icrm)=ggr/cp*z(k,icrm)
    end do ! k
   ! zi(nz,icrm) =  crm_input%zint(plev-nz+2)
    zi(nz,icrm) = crm_input%zint(icrm,plev-nz+2)-crm_input%zint(icrm,plev+1) !+++mhwang, 2012-02-04
    presi(nz,icrm) = crm_input%pint(icrm, plev-nz+2)/100.

    dz(icrm) = 0.5*(z(1,icrm)+z(2,icrm))
    do k=2,nzm
      adzw(k,icrm) = (z(k,icrm)-z(k-1,icrm))/dz(icrm)
    end do
    adzw(1,icrm)  = 1.
    adzw(nz,icrm) = adzw(nzm,icrm)
    !+++mhwang fix the adz bug. (adz needs to be consistent with zi)
    !2012-02-04 Minghuai Wang (minghuai.wang@pnnl.gov)
    do k=1, nzm
      adz(k,icrm)=(zi(k+1,icrm)-zi(k,icrm))/dz(icrm)
    end do

    do k = 1,nzm
      rho(k,icrm) = crm_input%pdel(icrm,plev-k+1)/ggr/(adz(k,icrm)*dz(icrm))
    end do
    do k=2,nzm
    ! rhow(k,icrm) = 0.5*(rho(k,icrm)+rho(k-1,icrm))
    !+++mhwang fix the rhow bug (rhow needes to be consistent with crm_input%pmid)
    !2012-02-04 Minghuai Wang (minghuai.wang@pnnl.gov)
      rhow(k,icrm) = (crm_input%pmid(icrm,plev-k+2)-crm_input%pmid(icrm,plev-k+1))/ggr/(adzw(k,icrm)*dz(icrm))
    end do
    rhow(1,icrm) = 2.*rhow(2,icrm) - rhow(3,icrm)
#ifdef CLUBB_CRM /* Fix extrapolation for 30 point grid */
    if (  2.*rhow(nzm,icrm) - rhow(nzm-1,icrm) > 0. ) then
       rhow(nz,icrm)= 2.*rhow(nzm,icrm) - rhow(nzm-1,icrm)
    else
       rhow(nz,icrm)= sqrt( rhow(nzm,icrm) )
    endif
#else
    rhow(nz,icrm)= 2.*rhow(nzm,icrm) - rhow(nzm-1,icrm)
#endif /*CLUBB_CRM*/

    !  Initialize CRM fields:
    u   (1:nx,1:ny,1:nzm,icrm) = crm_state%u_wind(icrm,1:nx,1:ny,1:nzm)
    v   (1:nx,1:ny,1:nzm,icrm) = crm_state%v_wind(icrm,1:nx,1:ny,1:nzm)*YES3D
    w   (1:nx,1:ny,1:nzm,icrm) = crm_state%w_wind(icrm,1:nx,1:ny,1:nzm)
    tabs(1:nx,1:ny,1:nzm,icrm) = crm_state%temperature(icrm,1:nx,1:ny,1:nzm)

    ! limit the velocity at the very first step:
    if(u(1,1,1,icrm).eq.u(2,1,1,icrm).and.u(3,1,2,icrm).eq.u(4,1,2,icrm)) then
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            u(i,j,k,icrm) = min( umax, max(-umax,u(i,j,k,icrm)) )
            v(i,j,k,icrm) = min( umax, max(-umax,v(i,j,k,icrm)) )*YES3D
          enddo
        enddo
      enddo
    endif


#if defined(SP_ESMT)
    do k=1,nzm
      u_esmt(:,:,k,icrm) = crm_input%ul_esmt(icrm,plev-k+1)
      v_esmt(:,:,k,icrm) = crm_input%vl_esmt(icrm,plev-k+1)
    end do
#endif

      ! Populate microphysics array from crm_state
#ifdef m2005
      micro_field(1:nx,1:ny,1:nzm,1,icrm)  = crm_state%qt(icrm,1:nx,1:ny,1:nzm)
      micro_field(1:nx,1:ny,1:nzm,2,icrm)  = crm_state%nc(icrm,1:nx,1:ny,1:nzm)
      micro_field(1:nx,1:ny,1:nzm,3,icrm)  = crm_state%qr(icrm,1:nx,1:ny,1:nzm)
      micro_field(1:nx,1:ny,1:nzm,4,icrm)  = crm_state%nr(icrm,1:nx,1:ny,1:nzm)
      micro_field(1:nx,1:ny,1:nzm,5,icrm)  = crm_state%qi(icrm,1:nx,1:ny,1:nzm)
      micro_field(1:nx,1:ny,1:nzm,6,icrm)  = crm_state%ni(icrm,1:nx,1:ny,1:nzm)
      micro_field(1:nx,1:ny,1:nzm,7,icrm)  = crm_state%qs(icrm,1:nx,1:ny,1:nzm)
      micro_field(1:nx,1:ny,1:nzm,8,icrm)  = crm_state%ns(icrm,1:nx,1:ny,1:nzm)
      micro_field(1:nx,1:ny,1:nzm,9,icrm)  = crm_state%qg(icrm,1:nx,1:ny,1:nzm)
      micro_field(1:nx,1:ny,1:nzm,10,icrm) = crm_state%ng(icrm,1:nx,1:ny,1:nzm)
      cloudliq(1:nx,1:ny,1:nzm,icrm) = crm_state%qc(icrm,1:nx,1:ny,1:nzm)
#else
      micro_field(1:nx,1:ny,1:nzm,1,icrm) = crm_state%qt(icrm,1:nx,1:ny,1:nzm)
      micro_field(1:nx,1:ny,1:nzm,2,icrm) = crm_state%qp(icrm,1:nx,1:ny,1:nzm)
      qn(1:nx,1:ny,1:nzm,icrm) = crm_state%qn(icrm,1:nx,1:ny,1:nzm)
#endif

#ifdef m2005
    do k=1, nzm
#ifdef MODAL_AERO
      ! set aerosol data
      l=plev-k+1
      naer (k, 1:ntot_amode,icrm) = crm_input%naermod (icrm,l, 1:ntot_amode)
      vaer (k, 1:ntot_amode,icrm) = crm_input%vaerosol(icrm,l, 1:ntot_amode)
      hgaer(k, 1:ntot_amode,icrm) = crm_input%hygro   (icrm,l, 1:ntot_amode)
#endif /* MODAL_AERO */
      do j=1, ny
        do i=1, nx
          if(cloudliq(i,j,k,icrm).gt.0) then
            if(dopredictNc) then
              if( micro_field(i,j,k,incl,icrm).eq.0) micro_field(i,j,k,incl,icrm) = 1.0e6*Nc0/rho(k,icrm)
            endif
          endif
        enddo
      enddo
    enddo
#endif /* m2005 */

    w(:,:,nz,icrm)=0.
    wsub (:,icrm) = 0.      !used in clubb, +++mhwang
    dudt(1:nx,1:ny,1:nzm,1:3,icrm) = 0.
    dvdt(1:nx,1:ny,1:nzm,1:3,icrm) = 0.
    dwdt(1:nx,1:ny,1:nz,1:3,icrm) = 0.
    tke (1:nx,1:ny,1:nzm,icrm) = 0.
    tk  (1:nx,1:ny,1:nzm,icrm) = 0.
    tkh (1:nx,1:ny,1:nzm,icrm) = 0.
    p   (1:nx,1:ny,1:nzm,icrm) = 0.

    CF3D(1:nx,1:ny,1:nzm,icrm) = 1.

    call micro_init(ncrms,icrm)

    ! initialize sgs fields
    call sgs_init(ncrms,icrm)

    colprec=0
    colprecs=0
    do k=1,nzm
      u0(k,icrm)=0.
      v0(k,icrm)=0.
      t0(k,icrm)=0.
      t00(k,icrm)=0.
      tabs0(k,icrm)=0.
      q0(k,icrm)=0.
      qv0(k,icrm)=0.
      !+++mhwang these are not initialized ??
      qn0(k,icrm) = 0.0
      qp0(k,icrm) = 0.0
      tke0(k,icrm) = 0.0
      !---mhwang
      do j=1,ny
        do i=1,nx
          t(i,j,k,icrm) = tabs(i,j,k,icrm)+gamaz(k,icrm) &
                    -fac_cond*qcl(i,j,k,icrm)-fac_sub*qci(i,j,k,icrm) &
                    -fac_cond*qpl(i,j,k,icrm)-fac_sub*qpi(i,j,k,icrm)
          colprec=colprec+(qpl(i,j,k,icrm)+qpi(i,j,k,icrm))*crm_input%pdel(icrm,plev-k+1)
          colprecs=colprecs+qpi(i,j,k,icrm)*crm_input%pdel(icrm,plev-k+1)
          u0(k,icrm)=u0(k,icrm)+u(i,j,k,icrm)
          v0(k,icrm)=v0(k,icrm)+v(i,j,k,icrm)
          t0(k,icrm)=t0(k,icrm)+t(i,j,k,icrm)
          t00(k,icrm)=t00(k,icrm)+t(i,j,k,icrm)+fac_cond*qpl(i,j,k,icrm)+fac_sub*qpi(i,j,k,icrm)
          tabs0(k,icrm)=tabs0(k,icrm)+tabs(i,j,k,icrm)
          q0(k,icrm)=q0(k,icrm)+qv(i,j,k,icrm)+qcl(i,j,k,icrm)+qci(i,j,k,icrm)
          qv0(k,icrm) = qv0(k,icrm) + qv(i,j,k,icrm)
          qn0(k,icrm) = qn0(k,icrm) + qcl(i,j,k,icrm) + qci(i,j,k,icrm)
          qp0(k,icrm) = qp0(k,icrm) + qpl(i,j,k,icrm) + qpi(i,j,k,icrm)
          tke0(k,icrm)=tke0(k,icrm)+tke(i,j,k,icrm)
        enddo
      enddo

      u0(k,icrm) = u0(k,icrm) * factor_xy
      v0(k,icrm) = v0(k,icrm) * factor_xy
      t0(k,icrm) = t0(k,icrm) * factor_xy
      t00(k,icrm) = t00(k,icrm) * factor_xy
      tabs0(k,icrm) = tabs0(k,icrm) * factor_xy
      q0(k,icrm) = q0(k,icrm) * factor_xy
      qv0(k,icrm) = qv0(k,icrm) * factor_xy
      qn0(k,icrm) = qn0(k,icrm) * factor_xy
      qp0(k,icrm) = qp0(k,icrm) * factor_xy
      tke0(k,icrm) = tke0(k,icrm) * factor_xy
#ifdef CLUBB_CRM
      ! Update thetav for CLUBB.  This is needed when we have a higher model top
      ! than is in the sounding, because we subsequently use tv0 to initialize
      ! thv_ds_zt/zm, which appear in CLUBB's anelastic buoyancy terms.
      ! -dschanen UWM 11 Feb 2010
      tv0(k,icrm) = tabs0(k,icrm)*prespot(k,icrm)*(1.+epsv*q0(k,icrm))
#endif /* CLUBB_CRM */

      l = plev-k+1
      uln(l,icrm) = min( umax, max(-umax,crm_input%ul(icrm,l)) )
      vln(l,icrm) = min( umax, max(-umax,crm_input%vl(icrm,l)) )*YES3D
      ttend(k,icrm) = (crm_input%tl(icrm,l)+gamaz(k,icrm)- fac_cond*(crm_input%qccl(icrm,l)+crm_input%qiil(icrm,l))-fac_fus*crm_input%qiil(icrm,l)-t00(k,icrm))*idt_gl
      qtend(k,icrm) = (crm_input%ql(icrm,l)+crm_input%qccl(icrm,l)+crm_input%qiil(icrm,l)-q0(k,icrm))*idt_gl
      utend(k,icrm) = (uln(l,icrm)-u0(k,icrm))*idt_gl
      vtend(k,icrm) = (vln(l,icrm)-v0(k,icrm))*idt_gl
      ug0(k,icrm) = uln(l,icrm)
      vg0(k,icrm) = vln(l,icrm)
      tg0(k,icrm) = crm_input%tl(icrm,l)+gamaz(k,icrm)-fac_cond*crm_input%qccl(icrm,l)-fac_sub*crm_input%qiil(icrm,l)
      qg0(k,icrm) = crm_input%ql(icrm,l)+crm_input%qccl(icrm,l)+crm_input%qiil(icrm,l)

    end do ! k

    uhl(icrm) = u0(1,icrm)
    vhl(icrm) = v0(1,icrm)

! estimate roughness length assuming logarithmic profile of velocity near the surface:

    ustar(icrm) = sqrt(crm_input%tau00(icrm)/rho(1,icrm))
    z0(icrm) = z0_est(z(1,icrm),bflx(icrm),wnd(icrm),ustar(icrm))
    z0(icrm) = max(real(0.00001,crm_rknd),min(real(1.,crm_rknd),z0(icrm)))

    crm_output%timing_factor = 0.

    crm_output%prectend(icrm)=colprec
    crm_output%precstend(icrm)=colprecs


#ifdef CLUBB_CRM
    if(doclubb) then
      fluxbu(:, :,icrm) = crm_input%fluxu00(icrm)/rhow(1,icrm)
      fluxbv(:, :,icrm) = crm_input%fluxv00(icrm)/rhow(1,icrm)
      fluxbt(:, :,icrm) = crm_input%fluxt00(icrm)/rhow(1,icrm)
      fluxbq(:, :,icrm) = crm_input%fluxq00(icrm)/rhow(1,icrm)
    else
      fluxbu(:, :,icrm) = 0.
      fluxbv(:, :,icrm) = 0.
      fluxbt(:, :,icrm) = 0.
      fluxbq(:, :,icrm) = 0.
    endif
#else
    fluxbu(:,:,icrm)=0.
    fluxbv(:,:,icrm)=0.
    fluxbt(:,:,icrm)=0.
    fluxbq(:,:,icrm)=0.
#endif /* CLUBB_CRM */
    fluxtu  (:,:,icrm)=0.
    fluxtv  (:,:,icrm)=0.
    fluxtt  (:,:,icrm)=0.
    fluxtq  (:,:,icrm)=0.
    fzero   (:,:,icrm) =0.
    precsfc (:,:,icrm)=0.
    precssfc(:,:,icrm)=0.

!---------------------------------------------------
    crm_output%cld   (icrm,:) = 0.
    crm_output%cldtop(icrm,:) = 0.
    crm_output%gicewp(icrm,:) = 0
    crm_output%gliqwp(icrm,:) = 0
    crm_output%mctot (icrm,:) = 0.
    crm_output%mcup  (icrm,:) = 0.
    crm_output%mcdn  (icrm,:) = 0.
    crm_output%mcuup (icrm,:) = 0.
    crm_output%mcudn (icrm,:) = 0.
    crm_output%qc_mean(icrm,:) = 0.
    crm_output%qi_mean(icrm,:) = 0.
    crm_output%qs_mean(icrm,:) = 0.
    crm_output%qg_mean(icrm,:) = 0.
    crm_output%qr_mean(icrm,:) = 0.
#ifdef m2005
    crm_output%nc_mean(icrm,:) = 0.
    crm_output%ni_mean(icrm,:) = 0.
    crm_output%ns_mean(icrm,:) = 0.
    crm_output%ng_mean(icrm,:) = 0.
    crm_output%nr_mean(icrm,:) = 0.
    ! hm 8/31/11 add new variables
    crm_output%aut_a (icrm,:) = 0.
    crm_output%acc_a (icrm,:) = 0.
    crm_output%evpc_a(icrm,:) = 0.
    crm_output%evpr_a(icrm,:) = 0.
    crm_output%mlt_a (icrm,:) = 0.
    crm_output%sub_a (icrm,:) = 0.
    crm_output%dep_a (icrm,:) = 0.
    crm_output%con_a (icrm,:) = 0.

    ! hm 8/31/11 add new output
    ! these are increments added to calculate gcm-grid and time-step avg
    ! note - these values are also averaged over the icycle loop following
    ! the approach for precsfc
    aut1a (:,:,:,icrm) = 0.
    acc1a (:,:,:,icrm) = 0.
    evpc1a(:,:,:,icrm) = 0.
    evpr1a(:,:,:,icrm) = 0.
    mlt1a (:,:,:,icrm) = 0.
    sub1a (:,:,:,icrm) = 0.
    dep1a (:,:,:,icrm) = 0.
    con1a (:,:,:,icrm) = 0.
#endif /* m2005 */

    crm_output%mu_crm (icrm,:) = 0.
    crm_output%md_crm (icrm,:) = 0.
    crm_output%eu_crm (icrm,:) = 0.
    crm_output%du_crm (icrm,:) = 0.
    crm_output%ed_crm (icrm,:) = 0.
    dd_crm (icrm,:) = 0.
    crm_output%jt_crm (icrm)   = 0.
    crm_output%mx_crm (icrm)   = 0.

    mui_crm(icrm,:) = 0.
    mdi_crm(icrm,:) = 0.

    crm_output%flux_qt   (icrm,:) = 0.
    crm_output%flux_u    (icrm,:) = 0.
    crm_output%flux_v    (icrm,:) = 0.
    crm_output%fluxsgs_qt(icrm,:) = 0.
    crm_output%tkez      (icrm,:) = 0.
    crm_output%tkesgsz   (icrm,:) = 0.
    crm_output%tkz       (icrm,:) = 0.
    crm_output%flux_qp   (icrm,:) = 0.
    crm_output%precflux  (icrm,:) = 0.
    crm_output%qt_trans  (icrm,:) = 0.
    crm_output%qp_trans  (icrm,:) = 0.
    crm_output%qp_fall   (icrm,:) = 0.
    crm_output%qp_evp    (icrm,:) = 0.
    crm_output%qp_src    (icrm,:) = 0.
    crm_output%qt_ls     (icrm,:) = 0.
    crm_output%t_ls      (icrm,:) = 0.

    uwle(:,icrm)     = 0.
    uwsb(:,icrm)     = 0.
    vwle(:,icrm)     = 0.
    vwsb(:,icrm)     = 0.
    qpsrc(:,icrm)    = 0.
    qpevp(:,icrm)    = 0.
    qpfall(:,icrm)   = 0.
    precflux(:,icrm) = 0.

!--------------------------------------------------
#ifdef sam1mom
    if(doprecip) call precip_init(ncrms,icrm)
#endif

    if ( igstep <= 1 ) then
        iseed = get_gcol_p(lchnk,icol(icrm)) * perturb_seed_scale
        call setperturb(ncrms,icrm,iseed)
    end if

    !--------------------------
    ! whannah - sanity check for new method to calculate radiation
    ! over averaged groups of columns instead of each individually
    if ( mod(nx,crm_nx_rad)==0 .or. mod(nx,crm_nx_rad)==0  ) then
      crm_nx_rad_fac = real(crm_nx_rad,crm_rknd)/real(nx,crm_rknd)
      crm_ny_rad_fac = real(crm_ny_rad,crm_rknd)/real(ny,crm_rknd)
    else
      write(0,*) "crm_nx_rad and crm_ny_rad need to be divisible by nx and ny"
      call endrun('crm main')
    end if

#ifndef CLUBB_CRM
    !--------------------------
    ! do a CLUBB sanity check
    if ( doclubb .or. doclubbnoninter ) then
      write(0,*) "Cannot call CLUBB if -DCLUBB is not in FFLAGS"
      call endrun('crm main')
    endif
#endif
#ifdef CLUBB_CRM
    !------------------------------------------------------------------
    ! Do initialization for UWM CLUBB
    !------------------------------------------------------------------
    up2       (1:nx, 1:ny, 1:nz ) = clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  1)
    vp2       (1:nx, 1:ny, 1:nz ) = clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  2)
    wprtp     (1:nx, 1:ny, 1:nz ) = clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  3)
    wpthlp    (1:nx, 1:ny, 1:nz ) = clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  4)
    wp2       (1:nx, 1:ny, 1:nz ) = clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  5)
    wp3       (1:nx, 1:ny, 1:nz ) = clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  6)
    rtp2      (1:nx, 1:ny, 1:nz ) = clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  7)
    thlp2     (1:nx, 1:ny, 1:nz ) = clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  8)
    rtpthlp   (1:nx, 1:ny, 1:nz ) = clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  9)
    upwp      (1:nx, 1:ny, 1:nz ) = clubb_buffer(icrm,1:nx, 1:ny, 1:nz , 10)
    vpwp      (1:nx, 1:ny, 1:nz ) = clubb_buffer(icrm,1:nx, 1:ny, 1:nz , 11)
    cloud_frac(1:nx, 1:ny, 1:nz ) = clubb_buffer(icrm,1:nx, 1:ny, 1:nz , 12)
    t_tndcy   (1:nx, 1:ny, 1:nzm) = clubb_buffer(icrm,1:nx, 1:ny, 1:nzm, 13)
    qc_tndcy  (1:nx, 1:ny, 1:nzm) = clubb_buffer(icrm,1:nx, 1:ny, 1:nzm, 14)
    qv_tndcy  (1:nx, 1:ny, 1:nzm) = clubb_buffer(icrm,1:nx, 1:ny, 1:nzm, 15)
    u_tndcy   (1:nx, 1:ny, 1:nzm) = clubb_buffer(icrm,1:nx, 1:ny, 1:nzm, 16)
    v_tndcy   (1:nx, 1:ny, 1:nzm) = clubb_buffer(icrm,1:nx, 1:ny, 1:nzm, 17)

    ! since no tracer is carried in the current version of MMF, these
    ! tracer-related restart varialbes are set to zero. +++mhwang, 2011-08
    tracer_tndcy = 0.0
    sclrp2       = 0.0
    sclrprtp     = 0.0
    sclrpthlp    = 0.0
    wpsclrp      = 0.0

    if((doclubb.and.docloud).or.(.not.doclubb .and. .not.docloud)) then
      write(0, *) 'doclubb and docloud can not both be true or be false'
      call endrun('crm_clubb2')
    endif
    if((doclubb_sfc_fluxes.and.docam_sfc_fluxes)) then
      write(0, *) 'doclubb_sfc_fluxes and dosam_sfc_fluxes can not both be true'
      call endrun('crm_clubb_fluxes')
    endif

    if ( doclubb .or. doclubbnoninter ) then
      call clubb_sgs_setup( real( dt*real( nclubb ), kind=time_precision), &
                            latitude(:,:,icrm), longitude(:,:,icrm), z(:,icrm), rho(:,icrm), zi(:,icrm), rhow(:,icrm), tv0(:,icrm), tke(:,:,:,icrm) )
    endif
#endif /* CLUBB_CRM */
  enddo

#ifdef ECPP
  ntavg1_ss = min(600._r8, dt_gl)   ! 10 minutes  or the GCM timestep, whichever smaller
  ntavg2_ss = dt_gl                 ! # of seconds to average between computing categories, must be a multiple of ntavgt1_ss.

  !!! ecpp_crm_init has to be called after ntavg1_ss and ntavg2_ss are set
  call ecpp_crm_init(ncrms)

  qlsink    = 0.0
  qlsink_bf = 0.0
  prain     = 0.0
  precr     = 0.0
  precsolid = 0.0
#endif /* ECPP */

  nstop = dt_gl/dt
  dt = dt_gl/nstop

#if defined( SP_CRM_SPLIT )
  nstop  = ceiling( nstop * 0.5 )
  crm_run_time  = dt_gl * 0.5
  icrm_run_time = 1._r8/crm_run_time
#else
  crm_run_time  = dt_gl
  icrm_run_time = 1._r8/crm_run_time
#endif
  factor_xyt = factor_xy / real(nstop,crm_rknd)


  !========================================================================================
  !----------------------------------------------------------------------------------------
  !   Main time loop
  !----------------------------------------------------------------------------------------
  !========================================================================================
  do nstep = 1 , nstop
    do icrm = 1 , ncrms
      crm_output%timing_factor(icrm) = crm_output%timing_factor(icrm)+1
    enddo

    !------------------------------------------------------------------
    !  Check if the dynamical time step should be decreased
    !  to handle the cases when the flow being locally linearly unstable
    !------------------------------------------------------------------
    call kurant(ncrms)

    do icyc=1,ncycle
      icycle = icyc
      dtn = dt/ncycle
      do icrm = 1 , ncrms
        dt3(na(icrm),icrm) = dtn
      enddo
      dtfactor = dtn/dt

      !---------------------------------------------
      !  	the Adams-Bashforth scheme in time
      call abcoefs(ncrms)

      !---------------------------------------------
      !  	initialize stuff:
      call zero(ncrms)

      !-----------------------------------------------------------
      !       Buoyancy term:
      call buoyancy(ncrms)

      !------------------------------------------------------------
      !       Large-scale and surface forcing:
      call forcing(ncrms)

      !!! Apply radiative tendency
      do icrm = 1 , ncrms
        do k=1,nzm
          do j=1,ny
            do i=1,nx
              i_rad = ceiling( real(i,crm_rknd) * crm_nx_rad_fac )
              j_rad = ceiling( real(j,crm_rknd) * crm_ny_rad_fac )
              tmp = crm_rad%qrad(icrm,i_rad,j_rad,k)*dtn
              !$acc atomic update
              t(i,j,k,icrm) = t(i,j,k,icrm) + tmp
            enddo
          enddo
        enddo
      enddo

      !----------------------------------------------------------
      !   	suppress turbulence near the upper boundary (spange):
      if (dodamping) call damping(ncrms)

      !---------------------------------------------------------
      !   Ice fall-out
#ifdef CLUBB_CRM
      if ( docloud .or. doclubb ) then
        call ice_fall(ncrms)
      endif
#else
      if(docloud) then
        call ice_fall(ncrms)
      endif
#endif

      !----------------------------------------------------------
      !     Update scalar boundaries after large-scale processes:
      call boundaries(ncrms,3)

      !---------------------------------------------------------
      !     Update boundaries for velocities:
      call boundaries(ncrms,0)

      !-----------------------------------------------
      !     surface fluxes:
      if (dosurface) call crmsurface(ncrms,bflx)

      !-----------------------------------------------------------
      !  SGS physics:
      if (dosgs) call sgs_proc(ncrms)

      !----------------------------------------------------------
      !     Fill boundaries for SGS diagnostic fields:
      call boundaries(ncrms,4)

      !-----------------------------------------------
      !       advection of momentum:
      call advect_mom(ncrms)

      do icrm = 1 , ncrms
        !----------------------------------------------------------
        !	SGS effects on momentum:
        if(dosgs) call sgs_mom(ncrms,icrm)

        !-----------------------------------------------------------
        !       Coriolis force:
        if (docoriolis) call coriolis(ncrms,icrm)

        !---------------------------------------------------------
        !       compute rhs of the Poisson equation and solve it for pressure.
        call pressure(ncrms,icrm)

        !---------------------------------------------------------
        !       find velocity field at n+1/2 timestep needed for advection of scalars:
        !  Note that at the end of the call, the velocities are in nondimensional form.
        call adams(ncrms,icrm)
      enddo

      !----------------------------------------------------------
      !     Update boundaries for all prognostic scalar fields for advection:
      call boundaries(ncrms,2)

      do icrm = 1 , ncrms
        !---------------------------------------------------------
        !      advection of scalars :
        call advect_all_scalars(ncrms,icrm)

        !-----------------------------------------------------------
        !    Convert velocity back from nondimensional form:
        call uvw(ncrms,icrm)
      enddo

      !----------------------------------------------------------
      !     Update boundaries for scalars to prepare for SGS effects:
      call boundaries(ncrms,3)

      do icrm = 1 , ncrms
        !---------------------------------------------------------
        !      SGS effects on scalars :
        if (dosgs) call sgs_scalars(ncrms,icrm)

        !-----------------------------------------------------------
        !       Calculate PGF for scalar momentum tendency
#if defined( SP_ESMT ) && defined( SP_ESMT_PGF )
            call scalar_momentum_tend()
#endif

        !-----------------------------------------------------------
        !       Cloud condensation/evaporation and precipitation processes:
#ifdef CLUBB_CRM
        if(docloud.or.dosmoke.or.doclubb) call micro_proc(ncrms,icrm)
#else
        if(docloud.or.dosmoke) call micro_proc(ncrms,icrm)
#endif /*CLUBB_CRM*/

        !-----------------------------------------------------------
        !    Compute diagnostics fields:
        call diagnose(ncrms,icrm)

        !----------------------------------------------------------
        ! Rotate the dynamic tendency arrays for Adams-bashforth scheme:
        nn=na(icrm)
        na(icrm)=nc(icrm)
        nc(icrm)=nb(icrm)
        nb(icrm)=nn
      enddo ! icrm
    enddo ! icycle

    do icrm = 1 , ncrms
      !----------------------------------------------------------
      !----------------------------------------------------------
#ifdef ECPP
      ! Here ecpp_crm_stat is called every CRM time step (dt), not every subcycle time step (dtn).
      ! This is what the original MMF model did (crm_rad%temperature, crm_rad%qv, ...). Do we want to call ecpp_crm_stat
      ! every subcycle time step??? +++mhwang
      call ecpp_crm_stat(ncrms,icrm)
#endif /*ECPP*/

      cwp (:,:,icrm) = 0.
      cwph(:,:,icrm) = 0.
      cwpm(:,:,icrm) = 0.
      cwpl(:,:,icrm) = 0.

      flag_top(:,:,icrm) = .true.

      cltemp(:,:,icrm) = 0.0; cmtemp(:,:,icrm) = 0.0
      chtemp(:,:,icrm) = 0.0; cttemp(:,:,icrm) = 0.0

      do k=1,nzm
        l = plev-k+1
        do j=1,ny
          do i=1,nx
            tmp1 = rho(nz-k,icrm)*adz(nz-k,icrm)*dz(icrm)*(qcl(i,j,nz-k,icrm)+qci(i,j,nz-k,icrm))
            cwp(i,j,icrm) = cwp(i,j,icrm)+tmp1
            cttemp(i,j,icrm) = max(CF3D(i,j,nz-k,icrm), cttemp(i,j,icrm))
            if(cwp(i,j,icrm).gt.cwp_threshold.and.flag_top(i,j,icrm)) then
                crm_output%cldtop(icrm,l) = crm_output%cldtop(icrm,l) + 1
                flag_top(i,j,icrm) = .false.
            endif
            if(pres(nz-k,icrm).ge.700.) then
                cwpl(i,j,icrm) = cwpl(i,j,icrm)+tmp1
                cltemp(i,j,icrm) = max(CF3D(i,j,nz-k,icrm), cltemp(i,j,icrm))
            else if(pres(nz-k,icrm).lt.400.) then
                cwph(i,j,icrm) = cwph(i,j,icrm)+tmp1
                chtemp(i,j,icrm) = max(CF3D(i,j,nz-k,icrm), chtemp(i,j,icrm))
            else
                cwpm(i,j,icrm) = cwpm(i,j,icrm)+tmp1
                cmtemp(i,j,icrm) = max(CF3D(i,j,nz-k,icrm), cmtemp(i,j,icrm))
            endif
            tmp1 = rho(k,icrm)*adz(k,icrm)*dz(icrm)
            if(tmp1*(qcl(i,j,k,icrm)+qci(i,j,k,icrm)).gt.cwp_threshold) then
                 crm_output%cld(icrm,l) = crm_output%cld(icrm,l) + CF3D(i,j,k,icrm)
                 if(w(i,j,k+1,icrm)+w(i,j,k,icrm).gt.2*wmin) then
                   crm_output%mcup (icrm,l) = crm_output%mcup (icrm,l) + rho(k,icrm)*0.5*(w(i,j,k+1,icrm)+w(i,j,k,icrm)) * CF3D(i,j,k,icrm)
                   crm_output%mcuup(icrm,l) = crm_output%mcuup(icrm,l) + rho(k,icrm)*0.5*(w(i,j,k+1,icrm)+w(i,j,k,icrm)) * (1.0 - CF3D(i,j,k,icrm))
                 endif
                 if(w(i,j,k+1,icrm)+w(i,j,k,icrm).lt.-2*wmin) then
                   crm_output%mcdn (icrm,l) = crm_output%mcdn (icrm,l) + rho(k,icrm)*0.5*(w(i,j,k+1,icrm)+w(i,j,k,icrm)) * CF3D(i,j,k,icrm)
                   crm_output%mcudn(icrm,l) = crm_output%mcudn(icrm,l) + rho(k,icrm)*0.5*(w(i,j,k+1,icrm)+w(i,j,k,icrm)) * (1. - CF3D(i,j,k,icrm))
                 endif
            else
                 if(w(i,j,k+1,icrm)+w(i,j,k,icrm).gt.2*wmin) then
                   crm_output%mcuup(icrm,l) = crm_output%mcuup(icrm,l) + rho(k,icrm)*0.5*(w(i,j,k+1,icrm)+w(i,j,k,icrm))
                 endif
                 if(w(i,j,k+1,icrm)+w(i,j,k,icrm).lt.-2*wmin) then
                   crm_output%mcudn(icrm,l) = crm_output%mcudn(icrm,l) + rho(k,icrm)*0.5*(w(i,j,k+1,icrm)+w(i,j,k,icrm))
                 endif
            endif

            !!! only collect radiative inputs during tphysbc() when using SP_CRM_SPLIT
            if ( phys_stage == 1 ) then

              !!! Reduced radiation method allows for fewer radiation calculations
              !!! by collecting statistics and doing radiation over column groups
              i_rad = ceiling( real(i,crm_rknd) * crm_nx_rad_fac )
              j_rad = ceiling( real(j,crm_rknd) * crm_ny_rad_fac )

              crm_rad%temperature  (icrm,i_rad,j_rad,k) = crm_rad%temperature  (icrm,i_rad,j_rad,k) + tabs(i,j,k,icrm)
              crm_rad%qv (icrm,i_rad,j_rad,k) = crm_rad%qv (icrm,i_rad,j_rad,k) + max(real(0.,crm_rknd),qv(i,j,k,icrm))
              crm_rad%qc (icrm,i_rad,j_rad,k) = crm_rad%qc (icrm,i_rad,j_rad,k) + qcl(i,j,k,icrm)
              crm_rad%qi (icrm,i_rad,j_rad,k) = crm_rad%qi (icrm,i_rad,j_rad,k) + qci(i,j,k,icrm)
              crm_rad%cld(icrm,i_rad,j_rad,k) = crm_rad%cld(icrm,i_rad,j_rad,k) + CF3D(i,j,k,icrm)
#ifdef m2005
              crm_rad%nc(icrm,i_rad,j_rad,k) = crm_rad%nc(icrm,i_rad,j_rad,k) + micro_field(i,j,k,incl,icrm)
              crm_rad%ni(icrm,i_rad,j_rad,k) = crm_rad%ni(icrm,i_rad,j_rad,k) + micro_field(i,j,k,inci,icrm)
              crm_rad%qs(icrm,i_rad,j_rad,k) = crm_rad%qs(icrm,i_rad,j_rad,k) + micro_field(i,j,k,iqs,icrm)
              crm_rad%ns(icrm,i_rad,j_rad,k) = crm_rad%ns(icrm,i_rad,j_rad,k) + micro_field(i,j,k,ins,icrm)
#endif
            endif

            crm_output%gliqwp(icrm,l) = crm_output%gliqwp(icrm,l) + qcl(i,j,k,icrm)
            crm_output%gicewp(icrm,l) = crm_output%gicewp(icrm,l) + qci(i,j,k,icrm)
          enddo
        enddo
      enddo

      ! Diagnose mass fluxes to drive CAM's convective transport of tracers.
      ! definition of mass fluxes is taken from Xu et al., 2002, QJRMS.
      do k=1, nzm+1
        l=plev+1-k+1
        do j=1, ny
          do i=1, nx
            if(w(i,j,k,icrm).gt.0.) then
              kx=max(1, k-1)
              qsat = qsatw_crm(tabs(i,j,kx,icrm),pres(kx,icrm))
              if(qcl(i,j,kx,icrm)+qci(i,j,kx,icrm).gt.min(real(1.e-5,crm_rknd),0.01*qsat)) then
                mui_crm(icrm,l) = mui_crm(icrm,l)+rhow(k,icrm)*w(i,j,k,icrm)
              endif
            else if (w(i,j,k,icrm).lt.0.) then
              kx=min(k+1, nzm)
              qsat = qsatw_crm(tabs(i,j,kx,icrm),pres(kx,icrm))
              if(qcl(i,j,kx,icrm)+qci(i,j,kx,icrm).gt.min(real(1.e-5,crm_rknd),0.01*qsat)) then
                mdi_crm(icrm,l) = mdi_crm(icrm,l)+rhow(k,icrm)*w(i,j,k,icrm)
              else if(qpl(i,j,kx,icrm)+qpi(i,j,kx,icrm).gt.1.0e-4) then
                mdi_crm(icrm,l) = mdi_crm(icrm,l)+rhow(k,icrm)*w(i,j,k,icrm)
              endif
            endif
          enddo
        enddo
      enddo

      do j=1,ny
        do i=1,nx
          if(cwp (i,j,icrm).gt.cwp_threshold) crm_output%cltot(icrm) = crm_output%cltot(icrm) + cttemp(i,j,icrm)
          if(cwph(i,j,icrm).gt.cwp_threshold) crm_output%clhgh(icrm) = crm_output%clhgh(icrm) + chtemp(i,j,icrm)
          if(cwpm(i,j,icrm).gt.cwp_threshold) crm_output%clmed(icrm) = crm_output%clmed(icrm) + cmtemp(i,j,icrm)
          if(cwpl(i,j,icrm).gt.cwp_threshold) crm_output%cllow(icrm) = crm_output%cllow(icrm) + cltemp(i,j,icrm)
        enddo
      enddo

    enddo ! icrm
  enddo ! nstep
  !========================================================================================
  !----------------------------------------------------------------------------------------
  ! End main time loop
  !----------------------------------------------------------------------------------------
  !========================================================================================

  do icrm = 1 , ncrms
    !!! only collect radiative inputs during tphysbc()
    if ( phys_stage == 1 ) then

      tmp1 = crm_nx_rad_fac * crm_ny_rad_fac / real(nstop,crm_rknd)

      crm_rad%temperature  (icrm,:,:,:) = crm_rad%temperature  (icrm,:,:,:) * tmp1
      crm_rad%qv (icrm,:,:,:) = crm_rad%qv (icrm,:,:,:) * tmp1
      crm_rad%qc (icrm,:,:,:) = crm_rad%qc (icrm,:,:,:) * tmp1
      crm_rad%qi (icrm,:,:,:) = crm_rad%qi (icrm,:,:,:) * tmp1
      crm_rad%cld(icrm,:,:,:) = crm_rad%cld(icrm,:,:,:) * tmp1
#ifdef m2005
      crm_rad%nc(icrm,:,:,:) = crm_rad%nc(icrm,:,:,:) * tmp1
      crm_rad%ni(icrm,:,:,:) = crm_rad%ni(icrm,:,:,:) * tmp1
      crm_rad%qs(icrm,:,:,:) = crm_rad%qs(icrm,:,:,:) * tmp1
      crm_rad%ns(icrm,:,:,:) = crm_rad%ns(icrm,:,:,:) * tmp1
#endif /* m2005 */

    endif

    ! no CRM tendencies above its top
    tln  (1:ptop-1,icrm) =   crm_input%tl(icrm,1:ptop-1)
    qln  (1:ptop-1,icrm) =   crm_input%ql(icrm,1:ptop-1)
    qccln(1:ptop-1,icrm) = crm_input%qccl(icrm,1:ptop-1)
    qiiln(1:ptop-1,icrm) = crm_input%qiil(icrm,1:ptop-1)
    uln  (1:ptop-1,icrm) =   crm_input%ul(icrm,1:ptop-1)
    vln  (1:ptop-1,icrm) =   crm_input%vl(icrm,1:ptop-1)

    !  Compute tendencies due to CRM:
    tln  (ptop:plev,icrm) = 0.
    qln  (ptop:plev,icrm) = 0.
    qccln(ptop:plev,icrm) = 0.
    qiiln(ptop:plev,icrm) = 0.
    uln  (ptop:plev,icrm) = 0.
    vln  (ptop:plev,icrm) = 0.

#if defined( SP_ESMT )
    uln_esmt(1:ptop-1,icrm)  = crm_input%ul_esmt(icrm,1:ptop-1)
    vln_esmt(1:ptop-1,icrm)  = crm_input%vl_esmt(icrm,1:ptop-1)
    uln_esmt(ptop:plev,icrm) = 0.
    vln_esmt(ptop:plev,icrm) = 0.
#endif /* SP_ESMT */

    colprec=0
    colprecs=0
    do k = 1,nzm
      l = plev-k+1
      do i=1,nx
        do j=1,ny
          colprec = colprec +(qpl(i,j,k,icrm)+qpi(i,j,k,icrm))*crm_input%pdel(icrm,plev-k+1)
          colprecs= colprecs+qpi(i,j,k,icrm)*crm_input%pdel(icrm,plev-k+1)
          tln(l,icrm)  = tln(l,icrm)  +tabs(i,j,k,icrm)
          qln(l,icrm)  = qln(l,icrm)  +qv(i,j,k,icrm)
          qccln(l,icrm)= qccln(l,icrm)+qcl(i,j,k,icrm)
          qiiln(l,icrm)= qiiln(l,icrm)+qci(i,j,k,icrm)
          uln(l,icrm)  = uln(l,icrm)  +u(i,j,k,icrm)
          vln(l,icrm)  = vln(l,icrm)  +v(i,j,k,icrm)

#if defined(SP_ESMT)
          uln_esmt(l,icrm) = uln_esmt(l,icrm)+u_esmt(i,j,k,icrm)
          vln_esmt(l,icrm) = vln_esmt(l,icrm)+v_esmt(i,j,k,icrm)
#endif
        enddo ! j
      enddo ! i
    enddo ! k

    tln  (ptop:plev,icrm) = tln  (ptop:plev,icrm) * factor_xy
    qln  (ptop:plev,icrm) = qln  (ptop:plev,icrm) * factor_xy
    qccln(ptop:plev,icrm) = qccln(ptop:plev,icrm) * factor_xy
    qiiln(ptop:plev,icrm) = qiiln(ptop:plev,icrm) * factor_xy
    uln  (ptop:plev,icrm) = uln  (ptop:plev,icrm) * factor_xy
    vln  (ptop:plev,icrm) = vln  (ptop:plev,icrm) * factor_xy

#if defined(SP_ESMT)
    uln_esmt(ptop:plev,icrm) = uln_esmt(ptop:plev,icrm) * factor_xy
    vln_esmt(ptop:plev,icrm) = vln_esmt(ptop:plev,icrm) * factor_xy

    crm_output%u_tend_esmt(icrm,:) = (uln_esmt(:,icrm) - crm_input%ul_esmt(icrm,:))*icrm_run_time
    crm_output%v_tend_esmt(icrm,:) = (vln_esmt(:,icrm) - crm_input%vl_esmt(icrm,:))*icrm_run_time

    ! don't use tendencies from two top levels,
    crm_output%u_tend_esmt(icrm,ptop:ptop+1) = 0.
    crm_output%v_tend_esmt(icrm,ptop:ptop+1) = 0.
#endif

#if defined(SPMOMTRANS)
    !!! resolved convective momentum transport (CMT) tendencies
    crm_output%ultend(icrm,:) = (uln(:,icrm) - crm_input%ul(icrm,:))*icrm_run_time
    crm_output%vltend(icrm,:) = (vln(:,icrm) - crm_input%vl(icrm,:))*icrm_run_time

    !!! don't use tendencies from two top levels
    crm_output%ultend(icrm,ptop:ptop+1) = 0.
    crm_output%vltend(icrm,ptop:ptop+1) = 0.
#endif /* SPMOMTRANS */

    crm_output%sltend (icrm,:) = cp * (tln  (:,icrm) - crm_input%tl  (icrm,:)) * icrm_run_time
    crm_output%qltend (icrm,:) =      (qln  (:,icrm) - crm_input%ql  (icrm,:)) * icrm_run_time
    crm_output%qcltend(icrm,:) =      (qccln(:,icrm) - crm_input%qccl(icrm,:)) * icrm_run_time
    crm_output%qiltend(icrm,:) =      (qiiln(:,icrm) - crm_input%qiil(icrm,:)) * icrm_run_time
    crm_output%prectend (icrm) = (colprec -crm_output%prectend (icrm))/ggr*factor_xy * icrm_run_time
    crm_output%precstend(icrm) = (colprecs-crm_output%precstend(icrm))/ggr*factor_xy * icrm_run_time

    !!! don't use CRM tendencies from two crm top levels
    !!! radiation tendencies are added back after the CRM call (see crm_physics_tend)
    crm_output%sltend (icrm,ptop:ptop+1) = 0.
    crm_output%qltend (icrm,ptop:ptop+1) = 0.
    crm_output%qcltend(icrm,ptop:ptop+1) = 0.
    crm_output%qiltend(icrm,ptop:ptop+1) = 0.
    !-------------------------------------------------------------
    !
    ! Save the last step to the permanent core:
    crm_state%u_wind  (icrm,1:nx,1:ny,1:nzm) = u   (1:nx,1:ny,1:nzm,icrm)
    crm_state%v_wind  (icrm,1:nx,1:ny,1:nzm) = v   (1:nx,1:ny,1:nzm,icrm)
    crm_state%w_wind  (icrm,1:nx,1:ny,1:nzm) = w   (1:nx,1:ny,1:nzm,icrm)
    crm_state%temperature  (icrm,1:nx,1:ny,1:nzm) = tabs(1:nx,1:ny,1:nzm,icrm)

#ifdef m2005
      crm_state%qt(icrm,1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,1,icrm)
      crm_state%nc(icrm,1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,2,icrm)
      crm_state%qr(icrm,1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,3,icrm)
      crm_state%nr(icrm,1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,4,icrm)
      crm_state%qi(icrm,1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,5,icrm)
      crm_state%ni(icrm,1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,6,icrm)
      crm_state%qs(icrm,1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,7,icrm)
      crm_state%ns(icrm,1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,8,icrm)
      crm_state%qg(icrm,1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,9,icrm)
      crm_state%ng(icrm,1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,10,icrm)
      crm_state%qc(icrm,1:nx,1:ny,1:nzm) = cloudliq(1:nx,1:ny,1:nzm,icrm)
#else
      crm_state%qt(icrm,1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,1,icrm)
      crm_state%qp(icrm,1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,2,icrm)
      crm_state%qn(icrm,1:nx,1:ny,1:nzm) = qn(1:nx,1:ny,1:nzm,icrm)
#endif

    crm_output%tk   (icrm,1:nx,1:ny,1:nzm) = tk  (1:nx, 1:ny, 1:nzm,icrm)
    crm_output%tkh  (icrm,1:nx,1:ny,1:nzm) = tkh (1:nx, 1:ny, 1:nzm,icrm)
#ifdef CLUBB_CRM
    clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  1) = up2       (1:nx, 1:ny, 1:nz )
    clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  2) = vp2       (1:nx, 1:ny, 1:nz )
    clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  3) = wprtp     (1:nx, 1:ny, 1:nz )
    clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  4) = wpthlp    (1:nx, 1:ny, 1:nz )
    clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  5) = wp2       (1:nx, 1:ny, 1:nz )
    clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  6) = wp3       (1:nx, 1:ny, 1:nz )
    clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  7) = rtp2      (1:nx, 1:ny, 1:nz )
    clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  8) = thlp2     (1:nx, 1:ny, 1:nz )
    clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  9) = rtpthlp   (1:nx, 1:ny, 1:nz )
    clubb_buffer(icrm,1:nx, 1:ny, 1:nz , 10) = upwp      (1:nx, 1:ny, 1:nz )
    clubb_buffer(icrm,1:nx, 1:ny, 1:nz , 11) = vpwp      (1:nx, 1:ny, 1:nz )
    clubb_buffer(icrm,1:nx, 1:ny, 1:nz , 12) = cloud_frac(1:nx, 1:ny, 1:nz )
    clubb_buffer(icrm,1:nx, 1:ny, 1:nzm, 13) = t_tndcy   (1:nx, 1:ny, 1:nzm)
    clubb_buffer(icrm,1:nx, 1:ny, 1:nzm, 14) = qc_tndcy  (1:nx, 1:ny, 1:nzm)
    clubb_buffer(icrm,1:nx, 1:ny, 1:nzm, 15) = qv_tndcy  (1:nx, 1:ny, 1:nzm)
    clubb_buffer(icrm,1:nx, 1:ny, 1:nzm, 16) = u_tndcy   (1:nx, 1:ny, 1:nzm)
    clubb_buffer(icrm,1:nx, 1:ny, 1:nzm, 17) = v_tndcy   (1:nx, 1:ny, 1:nzm)

    crm_cld    (icrm,1:nx, 1:ny, 1:nz ) = cloud_frac  (1:nx, 1:ny, 1:nz )
    clubb_tk   (icrm,1:nx, 1:ny, 1:nzm) = tk_clubb    (1:nx, 1:ny, 1:nzm)
    clubb_tkh  (icrm,1:nx, 1:ny, 1:nzm) = tkh_clubb   (1:nx, 1:ny, 1:nzm)
    relvar     (icrm,1:nx, 1:ny, 1:nzm) = relvarg     (1:nx, 1:ny, 1:nzm)
    accre_enhan(icrm,1:nx, 1:ny, 1:nzm) = accre_enhang(1:nx, 1:ny, 1:nzm)
    qclvar     (icrm,1:nx, 1:ny, 1:nzm) = qclvarg     (1:nx, 1:ny, 1:nzm)
#endif /* CLUBB_CRM */

    do k=1,nzm
     do j=1,ny
      do i=1,nx
        crm_output%qcl(icrm,i,j,k) = qcl(i,j,k,icrm)
        crm_output%qci(icrm,i,j,k) = qci(i,j,k,icrm)
        crm_output%qpl(icrm,i,j,k) = qpl(i,j,k,icrm)
        crm_output%qpi(icrm,i,j,k) = qpi(i,j,k,icrm)
#ifdef m2005
        crm_output%wvar(icrm,i,j,k) = wvar (i,j,k,icrm)
        crm_output%aut (icrm,i,j,k) = aut1 (i,j,k,icrm)
        crm_output%acc (icrm,i,j,k) = acc1 (i,j,k,icrm)
        crm_output%evpc(icrm,i,j,k) = evpc1(i,j,k,icrm)
        crm_output%evpr(icrm,i,j,k) = evpr1(i,j,k,icrm)
        crm_output%mlt (icrm,i,j,k) = mlt1 (i,j,k,icrm)
        crm_output%sub (icrm,i,j,k) = sub1 (i,j,k,icrm)
        crm_output%dep (icrm,i,j,k) = dep1 (i,j,k,icrm)
        crm_output%con (icrm,i,j,k) = con1 (i,j,k,icrm)
#endif /* m2005 */
        enddo
      enddo
    enddo
    crm_output%z0m     (icrm) = z0(icrm)
    crm_output%taux(icrm) = taux0(icrm) / dble(nstop)
    crm_output%tauy(icrm) = tauy0(icrm) / dble(nstop)

    !---------------------------------------------------------------
    !  Diagnostics:

    ! hm add 9/7/11, change from GCM-time step avg to end-of-timestep
    do k=1,nzm
      l = plev-k+1
      do j=1,ny
        do i=1,nx
          crm_output%qc_mean(icrm,l) = crm_output%qc_mean(icrm,l) + qcl(i,j,k,icrm)
          crm_output%qi_mean(icrm,l) = crm_output%qi_mean(icrm,l) + qci(i,j,k,icrm)
          crm_output%qr_mean(icrm,l) = crm_output%qr_mean(icrm,l) + qpl(i,j,k,icrm)
#ifdef sam1mom
          omg = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(i,j,k,icrm)-tgrmin)*a_gr))
          crm_output%qg_mean(icrm,l) = crm_output%qg_mean(icrm,l) + qpi(i,j,k,icrm)*omg
          crm_output%qs_mean(icrm,l) = crm_output%qs_mean(icrm,l) + qpi(i,j,k,icrm)*(1.-omg)
#else
          crm_output%qg_mean(icrm,l) = crm_output%qg_mean(icrm,l) + micro_field(i,j,k,iqg,icrm)
          crm_output%qs_mean(icrm,l) = crm_output%qs_mean(icrm,l) + micro_field(i,j,k,iqs,icrm)

          crm_output%nc_mean(icrm,l) = crm_output%nc_mean(icrm,l) + micro_field(i,j,k,incl,icrm)
          crm_output%ni_mean(icrm,l) = crm_output%ni_mean(icrm,l) + micro_field(i,j,k,inci,icrm)
          crm_output%nr_mean(icrm,l) = crm_output%nr_mean(icrm,l) + micro_field(i,j,k,inr,icrm)
          crm_output%ng_mean(icrm,l) = crm_output%ng_mean(icrm,l) + micro_field(i,j,k,ing,icrm)
          crm_output%ns_mean(icrm,l) = crm_output%ns_mean(icrm,l) + micro_field(i,j,k,ins,icrm)
#endif /* sam1mom */
        enddo
      enddo
    enddo

    crm_output%cld   (icrm,:) = min( 1._r8, crm_output%cld   (icrm,:) * factor_xyt )
    crm_output%cldtop(icrm,:) = min( 1._r8, crm_output%cldtop(icrm,:) * factor_xyt )
    crm_output%gicewp(icrm,:) = crm_output%gicewp(icrm,:)*crm_input%pdel(icrm,:)*1000./ggr * factor_xyt
    crm_output%gliqwp(icrm,:) = crm_output%gliqwp(icrm,:)*crm_input%pdel(icrm,:)*1000./ggr * factor_xyt
    crm_output%mcup  (icrm,:) = crm_output%mcup (icrm,:) * factor_xyt
    crm_output%mcdn  (icrm,:) = crm_output%mcdn (icrm,:) * factor_xyt
    crm_output%mcuup (icrm,:) = crm_output%mcuup(icrm,:) * factor_xyt
    crm_output%mcudn (icrm,:) = crm_output%mcudn(icrm,:) * factor_xyt
    crm_output%mctot (icrm,:) = crm_output%mcup(icrm,:) + crm_output%mcdn(icrm,:) + crm_output%mcuup(icrm,:) + crm_output%mcudn(icrm,:)

    crm_output%qc_mean(icrm,:) = crm_output%qc_mean(icrm,:) * factor_xy
    crm_output%qi_mean(icrm,:) = crm_output%qi_mean(icrm,:) * factor_xy
    crm_output%qs_mean(icrm,:) = crm_output%qs_mean(icrm,:) * factor_xy
    crm_output%qg_mean(icrm,:) = crm_output%qg_mean(icrm,:) * factor_xy
    crm_output%qr_mean(icrm,:) = crm_output%qr_mean(icrm,:) * factor_xy
#ifdef m2005
    crm_output%nc_mean(icrm,:) = crm_output%nc_mean(icrm,:) * factor_xy
    crm_output%ni_mean(icrm,:) = crm_output%ni_mean(icrm,:) * factor_xy
    crm_output%ns_mean(icrm,:) = crm_output%ns_mean(icrm,:) * factor_xy
    crm_output%ng_mean(icrm,:) = crm_output%ng_mean(icrm,:) * factor_xy
    crm_output%nr_mean(icrm,:) = crm_output%nr_mean(icrm,:) * factor_xy

    ! hm 8/31/11 new output, gcm-grid- and time-step avg
    ! add loop over i,j do get horizontal avg, and flip vertical array
    do k=1,nzm
      l = plev-k+1
      do j=1,ny
        do i=1,nx
          crm_output%aut_a (icrm,l) = crm_output%aut_a (icrm,l) + aut1a (i,j,k,icrm)
          crm_output%acc_a (icrm,l) = crm_output%acc_a (icrm,l) + acc1a (i,j,k,icrm)
          crm_output%evpc_a(icrm,l) = crm_output%evpc_a(icrm,l) + evpc1a(i,j,k,icrm)
          crm_output%evpr_a(icrm,l) = crm_output%evpr_a(icrm,l) + evpr1a(i,j,k,icrm)
          crm_output%mlt_a (icrm,l) = crm_output%mlt_a (icrm,l) + mlt1a (i,j,k,icrm)
          crm_output%sub_a (icrm,l) = crm_output%sub_a (icrm,l) + sub1a (i,j,k,icrm)
          crm_output%dep_a (icrm,l) = crm_output%dep_a (icrm,l) + dep1a (i,j,k,icrm)
          crm_output%con_a (icrm,l) = crm_output%con_a (icrm,l) + con1a (i,j,k,icrm)
        enddo
      enddo
    enddo

    ! note, rates are divded by dt to get mean rate over step
    crm_output%aut_a (icrm,:) = crm_output%aut_a (icrm,:) * factor_xyt / dt
    crm_output%acc_a (icrm,:) = crm_output%acc_a (icrm,:) * factor_xyt / dt
    crm_output%evpc_a(icrm,:) = crm_output%evpc_a(icrm,:) * factor_xyt / dt
    crm_output%evpr_a(icrm,:) = crm_output%evpr_a(icrm,:) * factor_xyt / dt
    crm_output%mlt_a (icrm,:) = crm_output%mlt_a (icrm,:) * factor_xyt / dt
    crm_output%sub_a (icrm,:) = crm_output%sub_a (icrm,:) * factor_xyt / dt
    crm_output%dep_a (icrm,:) = crm_output%dep_a (icrm,:) * factor_xyt / dt
    crm_output%con_a (icrm,:) = crm_output%con_a (icrm,:) * factor_xyt / dt
#endif /* m2005 */

    crm_output%precc (icrm) = 0.
    crm_output%precl (icrm) = 0.
    crm_output%precsc(icrm) = 0.
    crm_output%precsl(icrm) = 0.
    do j=1,ny
      do i=1,nx
#ifdef sam1mom
        precsfc(i,j,icrm) = precsfc(i,j,icrm)*dz(icrm)/dt/dble(nstop)
        precssfc(i,j,icrm) = precssfc(i,j,icrm)*dz(icrm)/dt/dble(nstop)
#endif /* sam1mom */
#ifdef m2005
        precsfc(i,j,icrm) = precsfc(i,j,icrm)*dz(icrm)/dt/dble(nstop)     !mm/s/dz --> mm/s
        precssfc(i,j,icrm) = precssfc(i,j,icrm)*dz(icrm)/dt/dble(nstop)   !mm/s/dz --> mm/s
#endif /* m2005 */
        if(precsfc(i,j,icrm).gt.10./86400.) then
           crm_output%precc (icrm) = crm_output%precc (icrm) + precsfc(i,j,icrm)
           crm_output%precsc(icrm) = crm_output%precsc(icrm) + precssfc(i,j,icrm)
        else
           crm_output%precl (icrm) = crm_output%precl (icrm) + precsfc(i,j,icrm)
           crm_output%precsl(icrm) = crm_output%precsl(icrm) + precssfc(i,j,icrm)
        endif
      enddo
    enddo
    crm_output%prec_crm(icrm,:,:) = precsfc(:,:,icrm)/1000.           !mm/s --> m/s
    crm_output%precc   (icrm)     = crm_output%precc (icrm)*factor_xy/1000.
    crm_output%precl   (icrm)     = crm_output%precl (icrm)*factor_xy/1000.
    crm_output%precsc  (icrm)     = crm_output%precsc(icrm)*factor_xy/1000.
    crm_output%precsl  (icrm)     = crm_output%precsl(icrm)*factor_xy/1000.

    crm_output%cltot(icrm) = crm_output%cltot(icrm) * factor_xyt
    crm_output%clhgh(icrm) = crm_output%clhgh(icrm) * factor_xyt
    crm_output%clmed(icrm) = crm_output%clmed(icrm) * factor_xyt
    crm_output%cllow(icrm) = crm_output%cllow(icrm) * factor_xyt

    crm_output%jt_crm(icrm) = plev * 1.0
    crm_output%mx_crm(icrm) = 1.0
    do k=1, plev
      crm_output%mu_crm(icrm,k)=0.5*(mui_crm(icrm,k)+mui_crm(icrm,k+1))
      crm_output%md_crm(icrm,k)=0.5*(mdi_crm(icrm,k)+mdi_crm(icrm,k+1))
      crm_output%mu_crm(icrm,k)=crm_output%mu_crm(icrm,k)*ggr/100.          !kg/m2/s --> mb/s
      crm_output%md_crm(icrm,k)=crm_output%md_crm(icrm,k)*ggr/100.          !kg/m2/s --> mb/s
      crm_output%eu_crm(icrm,k) = 0.
      if(mui_crm(icrm,k)-mui_crm(icrm,k+1).gt.0) then
        crm_output%eu_crm(icrm,k)=(mui_crm(icrm,k)-mui_crm(icrm,k+1))*ggr/crm_input%pdel(icrm,k)    !/s
      else
        crm_output%du_crm(icrm,k)=-1.0*(mui_crm(icrm,k)-mui_crm(icrm,k+1))*ggr/crm_input%pdel(icrm,k)   !/s
      endif
      if(mdi_crm(icrm,k+1)-mdi_crm(icrm,k).lt.0) then
        crm_output%ed_crm(icrm,k)=(mdi_crm(icrm,k)-mdi_crm(icrm,k+1))*ggr/crm_input%pdel(icrm,k) ! /s
      else
        dd_crm(icrm,k)=-1.*(mdi_crm(icrm,k)-mdi_crm(icrm,k+1))*ggr/crm_input%pdel(icrm,k)   !/s
      endif
      if(abs(crm_output%mu_crm(icrm,k)).gt.1.0e-15.or.abs(crm_output%md_crm(icrm,k)).gt.1.0e-15) then
        crm_output%jt_crm(icrm) = min(k*1.0_r8, crm_output%jt_crm(icrm))
        crm_output%mx_crm(icrm) = max(k*1.0_r8, crm_output%mx_crm(icrm))
      endif
    enddo

    !-------------------------------------------------------------
    !       Fluxes and other stat:
    !-------------------------------------------------------------
    do k=1,nzm
      u2z = 0.
      v2z = 0.
      w2z = 0.
      do j=1,ny
        do i=1,nx
          u2z = u2z+(u(i,j,k,icrm)-u0(k,icrm))**2
          v2z = v2z+(v(i,j,k,icrm)-v0(k,icrm))**2
          w2z = w2z+0.5*(w(i,j,k+1,icrm)**2+w(i,j,k,icrm)**2)
        enddo
      enddo
      !+++mhwang
      ! mkwsb, mkle, mkadv, mkdiff (also crm_output%flux_u, crm_output%flux_v,icrm) seem not calculted correclty in the spcam3.5 codes.
      ! Only values at the last time step are calculated, but is averaged over the entire GCM
      ! time step.
      !---mhwang

      tmp1 = dz(icrm)/rhow(k,icrm)
      tmp2 = tmp1/dtn                        ! dtn is calculated inside of the icyc loop.
                                             ! It seems wrong to use it here ???? +++mhwang
      mkwsb (k,:,icrm) = mkwsb (k,:,icrm) * tmp1*rhow(k,icrm) * factor_xy/nstop     !kg/m3/s --> kg/m2/s
      mkwle (k,:,icrm) = mkwle (k,:,icrm) * tmp2*rhow(k,icrm) * factor_xy/nstop     !kg/m3   --> kg/m2/s
      mkadv (k,:,icrm) = mkadv (k,:,icrm) * factor_xy*icrm_run_time     ! kg/kg  --> kg/kg/s
      mkdiff(k,:,icrm) = mkdiff(k,:,icrm) * factor_xy*icrm_run_time   ! kg/kg  --> kg/kg/s

      ! qpsrc, qpevp, qpfall in M2005 are calculated in micro_flux.
      qpsrc   (k,icrm) = qpsrc   (k,icrm) * factor_xy*icrm_run_time
      qpevp   (k,icrm) = qpevp   (k,icrm) * factor_xy*icrm_run_time
      qpfall  (k,icrm) = qpfall  (k,icrm) * factor_xy*icrm_run_time   ! kg/kg in M2005 ---> kg/kg/s
      precflux(k,icrm) = precflux(k,icrm) * factor_xy*dz(icrm)/dt/nstop  !kg/m2/dz in M2005 -->kg/m2/s or mm/s (idt_gl=1/dt/nstop)

      l = plev-k+1
      crm_output%flux_u    (icrm,l) = (uwle(k,icrm) + uwsb(k,icrm))*tmp1*factor_xy/nstop
      crm_output%flux_v    (icrm,l) = (vwle(k,icrm) + vwsb(k,icrm))*tmp1*factor_xy/nstop
#ifdef sam1mom
      crm_output%flux_qt   (icrm,l) = mkwle(k,1,icrm) + mkwsb(k,1,icrm)
      crm_output%fluxsgs_qt(icrm,l) = mkwsb(k,1,icrm)
      crm_output%flux_qp   (icrm,l) = mkwle(k,2,icrm) + mkwsb(k,2,icrm)
      crm_output%qt_trans  (icrm,l) = mkadv(k,1,icrm) + mkdiff(k,1,icrm)
      crm_output%qp_trans  (icrm,l) = mkadv(k,2,icrm) + mkdiff(k,2,icrm)
#endif /* sam1mom */
#ifdef m2005
      crm_output%flux_qt   (icrm,l) = mkwle(k,1   ,icrm) + mkwsb(k,1   ,icrm) +  &
                         mkwle(k,iqci,icrm) + mkwsb(k,iqci,icrm)
      crm_output%fluxsgs_qt(icrm,l) = mkwsb(k,1   ,icrm) + mkwsb(k,iqci,icrm)
      crm_output%flux_qp   (icrm,l) = mkwle(k,iqr,icrm) + mkwsb(k,iqr,icrm) +  &
                         mkwle(k,iqs,icrm) + mkwsb(k,iqs,icrm) + mkwle(k,iqg,icrm) + mkwsb(k,iqg,icrm)
      crm_output%qt_trans  (icrm,l) = mkadv (k,1,icrm) + mkadv (k,iqci,icrm) + &
                         mkdiff(k,1,icrm) + mkdiff(k,iqci,icrm)
      crm_output%qp_trans  (icrm,l) = mkadv (k,iqr,icrm) + mkadv (k,iqs,icrm) + mkadv (k,iqg,icrm) + &
                         mkdiff(k,iqr,icrm) + mkdiff(k,iqs,icrm) + mkdiff(k,iqg,icrm)
#endif /* m2005 */
      crm_output%tkesgsz   (icrm,l)= rho(k,icrm)*sum(tke(1:nx,1:ny,k,icrm))*factor_xy
      crm_output%tkez      (icrm,l)= rho(k,icrm)*0.5*(u2z+v2z*YES3D+w2z)*factor_xy + crm_output%tkesgsz(icrm,l)
      crm_output%tkz       (icrm,l) = sum(tk(1:nx, 1:ny, k,icrm)) * factor_xy
      crm_output%precflux      (icrm,l) = precflux(k,icrm)/1000.       !mm/s  -->m/s

      crm_output%qp_fall   (icrm,l) = qpfall(k,icrm)
      crm_output%qp_evp    (icrm,l) = qpevp(k,icrm)
      crm_output%qp_src    (icrm,l) = qpsrc(k,icrm)

      crm_output%qt_ls     (icrm,l) = qtend(k,icrm)
      crm_output%t_ls      (icrm,l) = ttend(k,icrm)
    enddo

#ifdef ECPP
    crm_ecpp_output%abnd         (icrm,:,:,:,:)=0.0
    crm_ecpp_output%abnd_tf      (icrm,:,:,:,:)=0.0
    crm_ecpp_output%massflxbnd   (icrm,:,:,:,:)=0.0
    crm_ecpp_output%acen         (icrm,:,:,:,:)=0.0
    crm_ecpp_output%acen_tf      (icrm,:,:,:,:)=0.0
    crm_ecpp_output%rhcen        (icrm,:,:,:,:)=0.0
    crm_ecpp_output%qcloudcen    (icrm,:,:,:,:)=0.0
    crm_ecpp_output%qicecen      (icrm,:,:,:,:)=0.0
    crm_ecpp_output%qlsinkcen    (icrm,:,:,:,:)=0.0
    crm_ecpp_output%precrcen     (icrm,:,:,:,:)=0.0
    crm_ecpp_output%precsolidcen (icrm,:,:,:,:)=0.0
    crm_ecpp_output%qlsink_bfcen (icrm,:,:,:,:)=0.0
    crm_ecpp_output%qlsink_avgcen(icrm,:,:,:,:)=0.0
    crm_ecpp_output%praincen     (icrm,:,:,:,:)=0.0

    crm_ecpp_output%wupthresh_bnd   (icrm,:)=0.0
    crm_ecpp_output%wdownthresh_bnd (icrm,:)=0.0
    crm_ecpp_output%wwqui_cen       (icrm,:)=0.0
    crm_ecpp_output%wwqui_bnd       (icrm,:)=0.0
    crm_ecpp_output%wwqui_cloudy_cen(icrm,:)=0.0
    crm_ecpp_output%wwqui_cloudy_bnd(icrm,:)=0.0

    ! default is clear, non-precipitating, and quiescent class
    crm_ecpp_output%abnd   (icrm,:,1,1,1)=1.0
    crm_ecpp_output%abnd_tf(icrm,:,1,1,1)=1.0
    crm_ecpp_output%acen   (icrm,:,1,1,1)=1.0
    crm_ecpp_output%acen_tf(icrm,:,1,1,1)=1.0

    do k=1, nzm
      l=plev-k+1
      crm_ecpp_output%acen            (icrm,l,:,:,:) = area_cen_sum        (k,:,1:ncls_ecpp_in,:,icrm)
      crm_ecpp_output%acen_tf         (icrm,l,:,:,:) = area_cen_final      (k,:,1:ncls_ecpp_in,:,icrm)
      crm_ecpp_output%rhcen           (icrm,l,:,:,:) = rh_cen_sum          (k,:,1:ncls_ecpp_in,:,icrm)
      crm_ecpp_output%qcloudcen       (icrm,l,:,:,:) = qcloud_cen_sum      (k,:,1:ncls_ecpp_in,:,icrm)
      crm_ecpp_output%qicecen         (icrm,l,:,:,:) = qice_cen_sum        (k,:,1:ncls_ecpp_in,:,icrm)
      crm_ecpp_output%qlsinkcen       (icrm,l,:,:,:) = qlsink_cen_sum      (k,:,1:ncls_ecpp_in,:,icrm)
      crm_ecpp_output%precrcen        (icrm,l,:,:,:) = precr_cen_sum       (k,:,1:ncls_ecpp_in,:,icrm)
      crm_ecpp_output%precsolidcen    (icrm,l,:,:,:) = precsolid_cen_sum   (k,:,1:ncls_ecpp_in,:,icrm)
      crm_ecpp_output%wwqui_cen       (icrm,l)       = wwqui_cen_sum       (k,icrm)
      crm_ecpp_output%wwqui_cloudy_cen(icrm,l)       = wwqui_cloudy_cen_sum(k,icrm)
      crm_ecpp_output%qlsink_bfcen    (icrm,l,:,:,:) = qlsink_bf_cen_sum   (k,:,1:ncls_ecpp_in,:,icrm)
      crm_ecpp_output%qlsink_avgcen   (icrm,l,:,:,:) = qlsink_avg_cen_sum  (k,:,1:ncls_ecpp_in,:,icrm)
      crm_ecpp_output%praincen        (icrm,l,:,:,:) = prain_cen_sum       (k,:,1:ncls_ecpp_in,:,icrm)
    enddo
    do k=1, nzm+1
      l=plev+1-k+1
      crm_ecpp_output%abnd            (icrm,l,:,:,:) = area_bnd_sum        (k,:,1:ncls_ecpp_in,:,icrm)
      crm_ecpp_output%abnd_tf         (icrm,l,:,:,:) = area_bnd_final      (k,:,1:ncls_ecpp_in,:,icrm)
      crm_ecpp_output%massflxbnd      (icrm,l,:,:,:) = mass_bnd_sum        (k,:,1:ncls_ecpp_in,:,icrm)
      crm_ecpp_output%wupthresh_bnd   (icrm,l)       = wup_thresh          (k,icrm)
      crm_ecpp_output%wdownthresh_bnd (icrm,l)       = wdown_thresh        (k,icrm)
      crm_ecpp_output%wwqui_bnd       (icrm,l)       = wwqui_bnd_sum       (k,icrm)
      crm_ecpp_output%wwqui_cloudy_bnd(icrm,l)       = wwqui_cloudy_bnd_sum(k,icrm)
    enddo
#endif /* ECPP */

    crm_output%timing_factor(icrm) = crm_output%timing_factor(icrm) / nstop

#ifdef CLUBB_CRM
    ! Deallocate CLUBB variables, etc.
    ! -UWM
    if ( doclubb .or. doclubbnoninter ) call clubb_sgs_cleanup( )
#endif
  enddo

#ifdef ECPP
  !!! Deallocate ECPP variables
  call ecpp_crm_cleanup ()
#endif

  deallocate( t00)
  deallocate( fluxbtmp)
  deallocate( fluxttmp)
  deallocate( tln)
  deallocate( qln)
  deallocate( qccln)
  deallocate( qiiln)
  deallocate( uln)
  deallocate( vln)
#if defined(SP_ESMT)
  deallocate( uln_esmt)
  deallocate( vln_esmt)
#endif
  deallocate( cwp)
  deallocate( cwph)
  deallocate( cwpm)
  deallocate( cwpl)
  deallocate( flag_top)
  deallocate( gcolindex)
  deallocate( cltemp)
  deallocate( cmtemp)
  deallocate( chtemp)
  deallocate( cttemp)
#ifdef CLUBB_CRM
  deallocate( rtm_integral_before )
  deallocate( rtm_integral_after )
  deallocate( thlm_integral_before)
  deallocate( thlm_integral_after)
  deallocate( thlm_before)
  deallocate( thlm_after)
  deallocate( rtm_column)
#endif /* CLUBB_CRM */
  deallocate( dd_crm  )
  deallocate( mui_crm )
  deallocate( mdi_crm )

  call deallocate_params()
  call deallocate_grid()
  call deallocate_tracers()
  call deallocate_sgs()
  call deallocate_vars()
  call deallocate_micro()
#ifdef sam1mom
  call deallocate_micro_params()
#endif
#if defined( SP_ESMT )
  call deallocate_scalar_momentum()
#endif

end subroutine crm

end module crm_module
