
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

subroutine crm(lchnk, icol, ncrms, dt_gl, plev, &
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
    use module_ecpp_crm_driver, only: ecpp_crm_stat, ecpp_crm_init, ecpp_crm_cleanup
    use ecppvars              , only: NCLASS_CL, ncls_ecpp_in, NCLASS_PR
#endif /* ECPP */
    use accelerate_crm_mod    , only: use_crm_accel, crm_accel_factor, crm_accel_nstop, accelerate_crm
    use cam_abortutils        , only: endrun
    use time_manager          , only: get_nstep

    implicit none

    !-----------------------------------------------------------------------------------------------
    ! Interface variable declarations
    !-----------------------------------------------------------------------------------------------

    integer , intent(in   ) :: lchnk                            ! chunk identifier (only for lat/lon and random seed)
    integer , intent(in   ) :: ncrms                            ! Number of "vector" GCM columns to push down into CRM for SIMD vectorization / more threading
    integer , intent(in   ) :: plev                             ! number of levels in parent model
    real(r8), intent(in   ) :: dt_gl                            ! global model's time step
    integer , intent(in   ) :: icol                (ncrms)      ! column identifier (only for lat/lon and random seed)
    type(crm_input_type),      intent(in   ) :: crm_input
    type(crm_state_type),      intent(inout) :: crm_state
    type(crm_rad_type), target,intent(inout) :: crm_rad
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
    type(crm_output_type), target,     intent(inout) :: crm_output

    !-----------------------------------------------------------------------------------------------
    ! Local variable declarations
    !-----------------------------------------------------------------------------------------------

    real(r8),       parameter :: umax = 0.5*crm_dx/crm_dt       ! maxumum ampitude of the l.s. wind
    real(r8),       parameter :: wmin = 2.                      ! minimum up/downdraft velocity for stat
    real(crm_rknd), parameter :: cwp_threshold = 0.001          ! threshold for cloud condensate for shaded fraction calculation
    integer,        parameter :: perturb_seed_scale = 1000      ! scaling value for setperturb() seed value (seed = gcol * perturb_seed_scale)
    real(r8)        :: crm_run_time                             ! length of CRM integration
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
    logical :: crm_accel_ceaseflag   ! indicates if accelerate_crm needs to be aborted for remainder of crm call

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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! These pointers are workarounds for OpenACC PGI bugs
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(crm_rknd), pointer :: crm_rad_qrad(:,:,:,:)
    real(crm_rknd), pointer :: crm_output_timing_factor(:)
    real(crm_rknd), pointer :: crm_output_cldtop(:,:)
    real(crm_rknd), pointer :: crm_output_cld(:,:)
    real(crm_rknd), pointer :: crm_output_mcup(:,:)
    real(crm_rknd), pointer :: crm_output_mcuup(:,:)
    real(crm_rknd), pointer :: crm_output_mcdn(:,:)
    real(crm_rknd), pointer :: crm_output_mcudn(:,:)
    real(crm_rknd), pointer :: crm_rad_temperature(:,:,:,:)
    real(crm_rknd), pointer :: crm_rad_qv(:,:,:,:)
    real(crm_rknd), pointer :: crm_rad_qc(:,:,:,:)
    real(crm_rknd), pointer :: crm_rad_qi(:,:,:,:)
    real(crm_rknd), pointer :: crm_rad_cld(:,:,:,:)
    real(crm_rknd), pointer :: crm_output_gliqwp(:,:)
    real(crm_rknd), pointer :: crm_output_gicewp(:,:)
    real(crm_rknd), pointer :: crm_output_cltot(:)
    real(crm_rknd), pointer :: crm_output_clhgh(:)
    real(crm_rknd), pointer :: crm_output_clmed(:)
    real(crm_rknd), pointer :: crm_output_cllow(:)

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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! These pointers are workarounds for OpenACC PGI bugs
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    crm_rad_qrad             => crm_rad%qrad                
    crm_output_timing_factor => crm_output%timing_factor    
    crm_output_cldtop        => crm_output%cldtop           
    crm_output_cld           => crm_output%cld              
    crm_output_mcup          => crm_output%mcup             
    crm_output_mcuup         => crm_output%mcuup            
    crm_output_mcdn          => crm_output%mcdn             
    crm_output_mcudn         => crm_output%mcudn            
    crm_output_mcudn         => crm_output%mcudn            
    crm_rad_temperature      => crm_rad%temperature         
    crm_rad_qv               => crm_rad%qv                  
    crm_rad_qc               => crm_rad%qc                  
    crm_rad_qi               => crm_rad%qi                  
    crm_rad_cld              => crm_rad%cld                 
    crm_output_gliqwp        => crm_output%gliqwp           
    crm_output_gicewp        => crm_output%gicewp           
    crm_output_cltot         => crm_output%cltot            
    crm_output_clhgh         => crm_output%clhgh            
    crm_output_clmed         => crm_output%clmed            
    crm_output_cllow         => crm_output%cllow            

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

  crm_accel_ceaseflag = .false.

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
    fcory(icrm,:) = fcor(icrm)
    fcorzy(icrm,:) = fcorz(icrm)
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
      z(icrm,k) = crm_input%zmid(icrm,plev-k+1) - crm_input%zint(icrm,plev+1)
      zi(k,icrm) = crm_input%zint(icrm,plev-k+2)- crm_input%zint(icrm,plev+1)
      pres(icrm,k) = crm_input%pmid(icrm,plev-k+1)/100.
      presi(icrm,k) = crm_input%pint(icrm,plev-k+2)/100.
      prespot(k,icrm)=(1000./pres(icrm,k))**(rgas/cp)
      bet(icrm,k) = ggr/crm_input%tl(icrm,plev-k+1)
      gamaz(icrm,k)=ggr/cp*z(icrm,k)
    end do ! k
   ! zi(nz,icrm) =  crm_input%zint(plev-nz+2)
    zi(nz,icrm) = crm_input%zint(icrm,plev-nz+2)-crm_input%zint(icrm,plev+1) !+++mhwang, 2012-02-04
    presi(icrm,nz) = crm_input%pint(icrm, plev-nz+2)/100.

    dz(icrm) = 0.5*(z(icrm,1)+z(icrm,2))
    do k=2,nzm
      adzw(icrm,k) = (z(icrm,k)-z(icrm,k-1))/dz(icrm)
    end do
    adzw(icrm,1)  = 1.
    adzw(icrm,nz) = adzw(icrm,nzm)
    !+++mhwang fix the adz bug. (adz needs to be consistent with zi)
    !2012-02-04 Minghuai Wang (minghuai.wang@pnnl.gov)
    do k=1, nzm
      adz(icrm,k)=(zi(k+1,icrm)-zi(k,icrm))/dz(icrm)
    end do

    do k = 1,nzm
      rho(icrm,k) = crm_input%pdel(icrm,plev-k+1)/ggr/(adz(icrm,k)*dz(icrm))
    end do
    do k=2,nzm
    ! rhow(icrm,k) = 0.5*(rho(icrm,k)+rho(icrm,k-1))
    !+++mhwang fix the rhow bug (rhow needes to be consistent with crm_input%pmid)
    !2012-02-04 Minghuai Wang (minghuai.wang@pnnl.gov)
      rhow(icrm,k) = (crm_input%pmid(icrm,plev-k+2)-crm_input%pmid(icrm,plev-k+1))/ggr/(adzw(icrm,k)*dz(icrm))
    end do
    rhow(icrm,1) = 2.*rhow(icrm,2) - rhow(icrm,3)
#ifdef CLUBB_CRM /* Fix extrapolation for 30 point grid */
    if (  2.*rhow(icrm,nzm) - rhow(icrm,nzm-1) > 0. ) then
       rhow(icrm,nz)= 2.*rhow(icrm,nzm) - rhow(icrm,nzm-1)
    else
       rhow(icrm,nz)= sqrt( rhow(icrm,nzm) )
    endif
#else
    rhow(icrm,nz)= 2.*rhow(icrm,nzm) - rhow(icrm,nzm-1)
#endif /*CLUBB_CRM*/

    !  Initialize CRM fields:
    u(icrm,1:nx,1:ny,1:nzm) = crm_state%u_wind(icrm,1:nx,1:ny,1:nzm)
    v(icrm,1:nx,1:ny,1:nzm) = crm_state%v_wind(icrm,1:nx,1:ny,1:nzm)*YES3D
    w(icrm,1:nx,1:ny,1:nzm) = crm_state%w_wind(icrm,1:nx,1:ny,1:nzm)
    tabs(icrm,1:nx,1:ny,1:nzm) = crm_state%temperature(icrm,1:nx,1:ny,1:nzm)

    ! limit the velocity at the very first step:
    if(u(icrm,1,1,1).eq.u(icrm,2,1,1).and.u(icrm,3,1,2).eq.u(icrm,4,1,2)) then
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            u(icrm,i,j,k) = min( umax, max(-umax,u(icrm,i,j,k)) )
            v(icrm,i,j,k) = min( umax, max(-umax,v(icrm,i,j,k)) )*YES3D
          enddo
        enddo
      enddo
    endif


#if defined(SP_ESMT)
    do k=1,nzm
      u_esmt(icrm,:,:,k) = crm_input%ul_esmt(icrm,plev-k+1)
      v_esmt(icrm,:,:,k) = crm_input%vl_esmt(icrm,plev-k+1)
    end do
#endif

      ! Populate microphysics array from crm_state
#ifdef m2005
      micro_field(icrm,1:nx,1:ny,1:nzm,1 )  = crm_state%qt(icrm,1:nx,1:ny,1:nzm)
      micro_field(icrm,1:nx,1:ny,1:nzm,2 )  = crm_state%nc(icrm,1:nx,1:ny,1:nzm)
      micro_field(icrm,1:nx,1:ny,1:nzm,3 )  = crm_state%qr(icrm,1:nx,1:ny,1:nzm)
      micro_field(icrm,1:nx,1:ny,1:nzm,4 )  = crm_state%nr(icrm,1:nx,1:ny,1:nzm)
      micro_field(icrm,1:nx,1:ny,1:nzm,5 )  = crm_state%qi(icrm,1:nx,1:ny,1:nzm)
      micro_field(icrm,1:nx,1:ny,1:nzm,6 )  = crm_state%ni(icrm,1:nx,1:ny,1:nzm)
      micro_field(icrm,1:nx,1:ny,1:nzm,7 )  = crm_state%qs(icrm,1:nx,1:ny,1:nzm)
      micro_field(icrm,1:nx,1:ny,1:nzm,8 )  = crm_state%ns(icrm,1:nx,1:ny,1:nzm)
      micro_field(icrm,1:nx,1:ny,1:nzm,9 )  = crm_state%qg(icrm,1:nx,1:ny,1:nzm)
      micro_field(icrm,1:nx,1:ny,1:nzm,10) = crm_state%ng(icrm,1:nx,1:ny,1:nzm)
      cloudliq(1:nx,1:ny,1:nzm,icrm) = crm_state%qc(icrm,1:nx,1:ny,1:nzm)
#else
      micro_field(icrm,1:nx,1:ny,1:nzm,1) = crm_state%qt(icrm,1:nx,1:ny,1:nzm)
      micro_field(icrm,1:nx,1:ny,1:nzm,2) = crm_state%qp(icrm,1:nx,1:ny,1:nzm)
      qn(icrm,1:nx,1:ny,1:nzm) = crm_state%qn(icrm,1:nx,1:ny,1:nzm)
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
              if( micro_field(icrm,i,j,k,incl).eq.0) micro_field(icrm,i,j,k,incl) = 1.0e6*Nc0/rho(icrm,k)
            endif
          endif
        enddo
      enddo
    enddo
#endif /* m2005 */

    w(icrm,:,:,nz)=0.
    wsub (:,icrm) = 0.      !used in clubb, +++mhwang
    dudt(icrm,1:nx,1:ny,1:nzm,1:3) = 0.
    dvdt(icrm,1:nx,1:ny,1:nzm,1:3) = 0.
    dwdt(icrm,1:nx,1:ny,1:nz,1:3) = 0.
    tke(icrm,1:nx,1:ny,1:nzm) = 0.
    tk(icrm,1:nx,1:ny,1:nzm) = 0.
    tkh(icrm,1:nx,1:ny,1:nzm) = 0.
    p(icrm,1:nx,1:ny,1:nzm) = 0.

    cf3d(icrm,1:nx,1:ny,1:nzm) = 1.
  enddo

  call micro_init(ncrms)

  do icrm = 1 , ncrms
    ! initialize sgs fields
    call sgs_init(ncrms,icrm)

    colprec=0
    colprecs=0
    do k=1,nzm
      u0(icrm,k)=0.
      v0(icrm,k)=0.
      t0(k,icrm)=0.
      t00(k,icrm)=0.
      tabs0(icrm,k)=0.
      q0(k,icrm)=0.
      qv0(icrm,k)=0.
      !+++mhwang these are not initialized ??
      qn0(icrm,k) = 0.0
      qp0(icrm,k) = 0.0
      tke0(k,icrm) = 0.0
      !---mhwang
      do j=1,ny
        do i=1,nx
          t(icrm,i,j,k) = tabs(icrm,i,j,k)+gamaz(icrm,k) &
                    -fac_cond*qcl(icrm,i,j,k)-fac_sub*qci(icrm,i,j,k) &
                    -fac_cond*qpl(icrm,i,j,k)-fac_sub*qpi(icrm,i,j,k)
          colprec=colprec+(qpl(icrm,i,j,k)+qpi(icrm,i,j,k))*crm_input%pdel(icrm,plev-k+1)
          colprecs=colprecs+qpi(icrm,i,j,k)*crm_input%pdel(icrm,plev-k+1)
          u0(icrm,k)=u0(icrm,k)+u(icrm,i,j,k)
          v0(icrm,k)=v0(icrm,k)+v(icrm,i,j,k)
          t0(k,icrm)=t0(k,icrm)+t(icrm,i,j,k)
          t00(k,icrm)=t00(k,icrm)+t(icrm,i,j,k)+fac_cond*qpl(icrm,i,j,k)+fac_sub*qpi(icrm,i,j,k)
          tabs0(icrm,k)=tabs0(icrm,k)+tabs(icrm,i,j,k)
          q0(k,icrm)=q0(k,icrm)+qv(icrm,i,j,k)+qcl(icrm,i,j,k)+qci(icrm,i,j,k)
          qv0(icrm,k) = qv0(icrm,k) + qv(icrm,i,j,k)
          qn0(icrm,k) = qn0(icrm,k) + qcl(icrm,i,j,k) + qci(icrm,i,j,k)
          qp0(icrm,k) = qp0(icrm,k) + qpl(icrm,i,j,k) + qpi(icrm,i,j,k)
          tke0(k,icrm)=tke0(k,icrm)+tke(icrm,i,j,k)
        enddo
      enddo

      u0(icrm,k) = u0(icrm,k) * factor_xy
      v0(icrm,k) = v0(icrm,k) * factor_xy
      t0(k,icrm) = t0(k,icrm) * factor_xy
      t00(k,icrm) = t00(k,icrm) * factor_xy
      tabs0(icrm,k) = tabs0(icrm,k) * factor_xy
      q0(k,icrm) = q0(k,icrm) * factor_xy
      qv0(icrm,k) = qv0(icrm,k) * factor_xy
      qn0(icrm,k) = qn0(icrm,k) * factor_xy
      qp0(icrm,k) = qp0(icrm,k) * factor_xy
      tke0(k,icrm) = tke0(k,icrm) * factor_xy
#ifdef CLUBB_CRM
      ! Update thetav for CLUBB.  This is needed when we have a higher model top
      ! than is in the sounding, because we subsequently use tv0 to initialize
      ! thv_ds_zt/zm, which appear in CLUBB's anelastic buoyancy terms.
      ! -dschanen UWM 11 Feb 2010
      tv0(k,icrm) = tabs0(icrm,k)*prespot(k,icrm)*(1.+epsv*q0(k,icrm))
#endif /* CLUBB_CRM */

      l = plev-k+1
      uln(l,icrm) = min( umax, max(-umax,crm_input%ul(icrm,l)) )
      vln(l,icrm) = min( umax, max(-umax,crm_input%vl(icrm,l)) )*YES3D
      ttend(icrm,k) = (crm_input%tl(icrm,l)+gamaz(icrm,k)- fac_cond*(crm_input%qccl(icrm,l)+crm_input%qiil(icrm,l))-fac_fus*crm_input%qiil(icrm,l)-t00(k,icrm))*idt_gl
      qtend(icrm,k) = (crm_input%ql(icrm,l)+crm_input%qccl(icrm,l)+crm_input%qiil(icrm,l)-q0(k,icrm))*idt_gl
      utend(icrm,k) = (uln(l,icrm)-u0(icrm,k))*idt_gl
      vtend(icrm,k) = (vln(l,icrm)-v0(icrm,k))*idt_gl
      ug0(icrm,k) = uln(l,icrm)
      vg0(icrm,k) = vln(l,icrm)
      tg0(k,icrm) = crm_input%tl(icrm,l)+gamaz(icrm,k)-fac_cond*crm_input%qccl(icrm,l)-fac_sub*crm_input%qiil(icrm,l)
      qg0(k,icrm) = crm_input%ql(icrm,l)+crm_input%qccl(icrm,l)+crm_input%qiil(icrm,l)

    end do ! k

    uhl(icrm) = u0(icrm,1)
    vhl(icrm) = v0(icrm,1)

! estimate roughness length assuming logarithmic profile of velocity near the surface:

    ustar(icrm) = sqrt(crm_input%tau00(icrm)/rho(icrm,1))
    z0(icrm) = z0_est(z(icrm,1),bflx(icrm),wnd(icrm),ustar(icrm))
    z0(icrm) = max(real(0.00001,crm_rknd),min(real(1.,crm_rknd),z0(icrm)))

    crm_output%timing_factor = 0.

    crm_output%prectend(icrm)=colprec
    crm_output%precstend(icrm)=colprecs


#ifdef CLUBB_CRM
    if(doclubb) then
      fluxbu(icrm,:, :) = crm_input%fluxu00(icrm)/rhow(icrm,1)
      fluxbv(icrm,:, :) = crm_input%fluxv00(icrm)/rhow(icrm,1)
      fluxbt(:, :,icrm) = crm_input%fluxt00(icrm)/rhow(icrm,1)
      fluxbq(icrm,:, :) = crm_input%fluxq00(icrm)/rhow(icrm,1)
    else
      fluxbu(icrm,:, :) = 0.
      fluxbv(icrm,:, :) = 0.
      fluxbt(:, :,icrm) = 0.
      fluxbq(icrm,:, :) = 0.
    endif
#else
    fluxbu(icrm,:,:)=0.
    fluxbv(icrm,:,:)=0.
    fluxbt(:,:,icrm)=0.
    fluxbq(icrm,:,:)=0.
#endif /* CLUBB_CRM */
    fluxtu(icrm,:,:)=0.
    fluxtv(icrm,:,:)=0.
    fluxtt  (:,:,icrm)=0.
    fluxtq(icrm,:,:)=0.
    fzero   (:,:,icrm) =0.
    precsfc(icrm,:,:)=0.
    precssfc(icrm,:,:)=0.

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

    uwle(icrm,:)     = 0.
    uwsb(icrm,:)     = 0.
    vwle(icrm,:)     = 0.
    vwsb(icrm,:)     = 0.
    qpsrc(icrm,:)    = 0.
    qpevp(icrm,:)    = 0.
    qpfall(icrm,:)   = 0.
    precflux(icrm,:) = 0.

  enddo
!--------------------------------------------------
#ifdef sam1mom
  if(doprecip) call precip_init(ncrms)
#endif

  do icrm = 1 , ncrms
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
                            latitude(:,:,icrm), longitude(:,:,icrm), z(icrm,:), rho(icrm,:), zi(:,icrm), rhow(icrm,:), tv0(:,icrm), tke(icrm,:,:,:) )
    endif
#endif /* CLUBB_CRM */
  enddo

#ifdef ECPP
  call ecpp_crm_init(ncrms,dt_gl)

  qlsink    = 0.0
  qlsink_bf = 0.0
  prain     = 0.0
  precr     = 0.0
  precsolid = 0.0
#endif /* ECPP */

  nstop = dt_gl/dt
  dt = dt_gl/nstop

  crm_run_time  = dt_gl
  icrm_run_time = 1._r8/crm_run_time

  if (use_crm_accel) then
    call crm_accel_nstop(nstop)  ! reduce nstop by factor of (1 + crm_accel_factor)
  end if


  !$acc enter data copyin(dudt,dvdt,dwdt,misc,adz,bet,tabs0,qv,qv0,qcl,qci,qn0,qpl,qpi,qp0,tabs,t,micro_field,ttend,qtend,utend,vtend,u,u0,v,v0,w,t0,dz,precsfc,precssfc,rho,qifall,tlatqi) async(asyncid)
  !$acc enter data copyin(sstxy,taux0,tauy0,z,z0,fluxbu,fluxbv,bflx,uhl,vhl,adzw,presi,tkelediss,tkesbdiss,tkesbshear,tkesbbuoy,grdf_x,grdf_y,grdf_z,fcory,fcorzy,ug0,vg0,t01,q01,p0,pres,p) async(asyncid)
  !$acc enter data copyin(rhow,uwle,vwle,uwsb,vwsb,w_max,u_max,dt3,cwp,cwph,cwpm,cwpl,flag_top,cltemp,cmtemp,chtemp,cttemp,mkadv,mkwle,sgsadv,sgswle,gamaz,iw_xy,cw_xy,pw_xy,u200_xy,v200_xy) async(asyncid)
  !$acc enter data copyin(usfc_xy,vsfc_xy,w500_xy,swvp_xy,psfc_xy,u850_xy,v850_xy,cloudtopheight,cloudtoptemp,echotopheight,cld_xy,crm_output_timing_factor,crm_rad_qrad,cf3d) async(asyncid)
  !$acc enter data copyin(crm_output_mcudn,crm_output_mcup,crm_output_cld,crm_output_mcdn,crm_output_gliqwp,crm_output_mcuup,crm_rad_qc,crm_rad_cld,crm_rad_qi,crm_rad_temperature) async(asyncid)
  !$acc enter data copyin(crm_rad_qv,crm_output_gicewp,crm_output_cldtop,mdi_crm,mui_crm,crm_output_cltot,crm_output_clhgh,crm_output_clmed,crm_output_cllow,fluxbt,fluxtt,tdiff,twsb,fzero) async(asyncid)
  !$acc enter data copyin(fluxbq,fluxbmk,fluxtq,fluxtmk,sgswsb,mkdiff,mkwsb,qn,qpsrc,qpevp,accrrc,accrsc,accrsi,accrgi,accrgc,coefice,evapg1,evapg2,evapr1,evaps2,evaps1,evapr2) async(asyncid)

  !$acc enter data copyin(sgs_field,sgs_field_diag,tke2,tk2,tk,tke,tkh,twle,tadv,q0,qpfall,tlat,precflux,prec_xy,fluxtu,fluxtv) async(asyncid)

  !========================================================================================
  !----------------------------------------------------------------------------------------
  !   Main time loop
  !----------------------------------------------------------------------------------------
  !========================================================================================
  nstep = 0
  do while (nstep < nstop)
    nstep = nstep + 1

    !$acc parallel loop copy(crm_output_timing_factor) async(asyncid)
    do icrm = 1 , ncrms
      crm_output_timing_factor(icrm) = crm_output_timing_factor(icrm)+1
    enddo

    !------------------------------------------------------------------
    !  Check if the dynamical time step should be decreased
    !  to handle the cases when the flow being locally linearly unstable
    !------------------------------------------------------------------
    call kurant(ncrms)
    !$acc wait(asyncid)

    do icyc=1,ncycle
      icycle = icyc
      dtn = dt/ncycle
      dt3(na) = dtn
      !$acc update device(dt3) async(asyncid)
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
      !$acc parallel loop collapse(4) private(i_rad,j_rad) copy(t) copyin(crm_rad_qrad) async(asyncid)
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              i_rad = ceiling( real(i,crm_rknd) * crm_nx_rad_fac )
              j_rad = ceiling( real(j,crm_rknd) * crm_ny_rad_fac )
              t(icrm,i,j,k) = t(icrm,i,j,k) + crm_rad_qrad(icrm,i_rad,j_rad,k)*dtn
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

      !----------------------------------------------------------
      !	SGS effects on momentum:
      if(dosgs) call sgs_mom(ncrms)

      !-----------------------------------------------------------
      !       Coriolis force:
      if (docoriolis) call coriolis(ncrms)

      !---------------------------------------------------------
      !       compute rhs of the Poisson equation and solve it for pressure.
      call pressure(ncrms)

      !---------------------------------------------------------
      !       find velocity field at n+1/2 timestep needed for advection of scalars:
      !  Note that at the end of the call, the velocities are in nondimensional form.
      call adams(ncrms)

      !----------------------------------------------------------
      !     Update boundaries for all prognostic scalar fields for advection:
      call boundaries(ncrms,2)

      !---------------------------------------------------------
      !      advection of scalars :
      call advect_all_scalars(ncrms)

      !-----------------------------------------------------------
      !    Convert velocity back from nondimensional form:
      call uvw(ncrms)

      !----------------------------------------------------------
      !     Update boundaries for scalars to prepare for SGS effects:
      call boundaries(ncrms,3)

      !---------------------------------------------------------
      !      SGS effects on scalars :
      if (dosgs) call sgs_scalars(ncrms)

      !-----------------------------------------------------------
      !       Calculate PGF for scalar momentum tendency
#if defined( SP_ESMT ) && defined( SP_ESMT_PGF )
      call scalar_momentum_tend(ncrms)
#endif

      !-----------------------------------------------------------
      !       Cloud condensation/evaporation and precipitation processes:
#ifdef CLUBB_CRM
      if(docloud.or.dosmoke.or.doclubb) call micro_proc(ncrms)
#else
      if(docloud.or.dosmoke) call micro_proc(ncrms)
#endif /*CLUBB_CRM*/

      !-----------------------------------------------------------
      !       Apply mean-state acceleration
      if (use_crm_accel .and. .not. crm_accel_ceaseflag) then
        ! Use Jones-Bretherton-Pritchard methodology to accelerate
        ! CRM horizontal mean evolution artificially.
        call accelerate_crm(ncrms, nstep, nstop, crm_accel_ceaseflag)
      endif

      !-----------------------------------------------------------
      !    Compute diagnostics fields:
      call diagnose(ncrms)

      !----------------------------------------------------------
      ! Rotate the dynamic tendency arrays for Adams-bashforth scheme:
      nn=na
      na=nc
      nc=nb
      nb=nn
    enddo ! icycle

#ifdef ECPP
    ! Here ecpp_crm_stat is called every CRM time step (dt), not every subcycle time step (dtn).
    ! This is what the original MMF model did (crm_rad%temperature, crm_rad%qv, ...). Do we want to call ecpp_crm_stat
    ! every subcycle time step??? +++mhwang
    call ecpp_crm_stat(ncrms)
#endif
    !$acc parallel loop collapse(3) copy(cwp,cwph,cwpm,cwpl,flag_top,cltemp,cmtemp,chtemp,cttemp) async(asyncid)
    do icrm = 1 , ncrms
      do j = 1 , ny
        do i = 1 , nx
          cwp (i,j,icrm) = 0.
          cwph(i,j,icrm) = 0.
          cwpm(i,j,icrm) = 0.
          cwpl(i,j,icrm) = 0.

          flag_top(i,j,icrm) = .true.

          cltemp(i,j,icrm) = 0.0; cmtemp(i,j,icrm) = 0.0
          chtemp(i,j,icrm) = 0.0; cttemp(i,j,icrm) = 0.0
        enddo
      enddo
    enddo

    !$acc parallel loop gang vector collapse(3) copyin(cf3d,pres,qci,qv,dz,adz,w,tabs,qcl,rho) &
    !$acc& copy(crm_output_mcudn,crm_output_mcup,cwp,cltemp,cwpl,flag_top,crm_output_cld,cwpm,cttemp,cmtemp,cwph,chtemp,crm_output_mcdn,crm_output_gliqwp,&
    !$acc&      crm_output_mcuup,crm_rad_qc,crm_rad_cld,crm_rad_qi,crm_rad_temperature,crm_rad_qv,crm_output_gicewp,crm_output_cldtop) async(asyncid)
    do icrm = 1 , ncrms
      do j=1,ny
        do i=1,nx
          do k=1,nzm
            l = plev-k+1
            tmp1 = rho(icrm,nz-k)*adz(icrm,nz-k)*dz(icrm)*(qcl(icrm,i,j,nz-k)+qci(icrm,i,j,nz-k))
            cwp(i,j,icrm) = cwp(i,j,icrm)+tmp1
            cttemp(i,j,icrm) = max(cf3d(icrm,i,j,nz-k), cttemp(i,j,icrm))
            if(cwp(i,j,icrm).gt.cwp_threshold.and.flag_top(i,j,icrm)) then
                !$acc atomic update
                crm_output_cldtop(icrm,l) = crm_output_cldtop(icrm,l) + 1
                flag_top(i,j,icrm) = .false.
            endif
            if(pres(icrm,nz-k).ge.700.) then
                cwpl(i,j,icrm) = cwpl(i,j,icrm)+tmp1
                cltemp(i,j,icrm) = max(cf3d(icrm,i,j,nz-k), cltemp(i,j,icrm))
            else if(pres(icrm,nz-k).lt.400.) then
                cwph(i,j,icrm) = cwph(i,j,icrm)+tmp1
                chtemp(i,j,icrm) = max(cf3d(icrm,i,j,nz-k), chtemp(i,j,icrm))
            else
                cwpm(i,j,icrm) = cwpm(i,j,icrm)+tmp1
                cmtemp(i,j,icrm) = max(cf3d(icrm,i,j,nz-k), cmtemp(i,j,icrm))
            endif
            tmp1 = rho(icrm,k)*adz(icrm,k)*dz(icrm)
            if(tmp1*(qcl(icrm,i,j,k)+qci(icrm,i,j,k)).gt.cwp_threshold) then
                 !$acc atomic update
                 crm_output_cld(icrm,l) = crm_output_cld(icrm,l) + cf3d(icrm,i,j,k)
                 if(w(icrm,i,j,k+1)+w(icrm,i,j,k).gt.2*wmin) then
                   tmp = rho(icrm,k)*0.5*(w(icrm,i,j,k+1)+w(icrm,i,j,k)) * cf3d(icrm,i,j,k)
                   !$acc atomic update
                   crm_output_mcup (icrm,l) = crm_output_mcup (icrm,l) + tmp
                   tmp = rho(icrm,k)*0.5*(w(icrm,i,j,k+1)+w(icrm,i,j,k)) * (1.0 - cf3d(icrm,i,j,k))
                   !$acc atomic update
                   crm_output_mcuup(icrm,l) = crm_output_mcuup(icrm,l) + tmp
                 endif
                 if(w(icrm,i,j,k+1)+w(icrm,i,j,k).lt.-2*wmin) then
                   tmp = rho(icrm,k)*0.5*(w(icrm,i,j,k+1)+w(icrm,i,j,k)) * cf3d(icrm,i,j,k)
                   !$acc atomic update
                   crm_output_mcdn (icrm,l) = crm_output_mcdn (icrm,l) + tmp
                   tmp = rho(icrm,k)*0.5*(w(icrm,i,j,k+1)+w(icrm,i,j,k)) * (1. - cf3d(icrm,i,j,k))
                   !$acc atomic update
                   crm_output_mcudn(icrm,l) = crm_output_mcudn(icrm,l) + tmp
                 endif
            else
                 if(w(icrm,i,j,k+1)+w(icrm,i,j,k).gt.2*wmin) then
                   tmp = rho(icrm,k)*0.5*(w(icrm,i,j,k+1)+w(icrm,i,j,k))
                   !$acc atomic update
                   crm_output_mcuup(icrm,l) = crm_output_mcuup(icrm,l) + tmp
                 endif
                 if(w(icrm,i,j,k+1)+w(icrm,i,j,k).lt.-2*wmin) then
                    tmp = rho(icrm,k)*0.5*(w(icrm,i,j,k+1)+w(icrm,i,j,k))
                   !$acc atomic update
                   crm_output_mcudn(icrm,l) = crm_output_mcudn(icrm,l) + tmp
                 endif
            endif


            !!! Reduced radiation method allows for fewer radiation calculations
            !!! by collecting statistics and doing radiation over column groups
            i_rad = ceiling( real(i,crm_rknd) * crm_nx_rad_fac )
            j_rad = ceiling( real(j,crm_rknd) * crm_ny_rad_fac )

            !$acc atomic update
            crm_rad_temperature(icrm,i_rad,j_rad,k) = crm_rad_temperature(icrm,i_rad,j_rad,k) + tabs(icrm,i,j,k)
            tmp = max(real(0.,crm_rknd),qv(icrm,i,j,k))
            !$acc atomic update
            crm_rad_qv         (icrm,i_rad,j_rad,k) = crm_rad_qv         (icrm,i_rad,j_rad,k) + tmp
            !$acc atomic update
            crm_rad_qc         (icrm,i_rad,j_rad,k) = crm_rad_qc         (icrm,i_rad,j_rad,k) + qcl(icrm,i,j,k)
            !$acc atomic update
            crm_rad_qi         (icrm,i_rad,j_rad,k) = crm_rad_qi         (icrm,i_rad,j_rad,k) + qci(icrm,i,j,k)
            !$acc atomic update
            crm_rad_cld        (icrm,i_rad,j_rad,k) = crm_rad_cld        (icrm,i_rad,j_rad,k) + cf3d(icrm,i,j,k)
#ifdef m2005
            !$acc atomic update
            crm_rad%nc         (icrm,i_rad,j_rad,k) = crm_rad%nc         (icrm,i_rad,j_rad,k) + micro_field(icrm,i,j,k,incl)
            !$acc atomic update
            crm_rad%ni         (icrm,i_rad,j_rad,k) = crm_rad%ni         (icrm,i_rad,j_rad,k) + micro_field(icrm,i,j,k,inci)
            !$acc atomic update
            crm_rad%qs         (icrm,i_rad,j_rad,k) = crm_rad%qs         (icrm,i_rad,j_rad,k) + micro_field(icrm,i,j,k,iqs )
            !$acc atomic update
            crm_rad%ns         (icrm,i_rad,j_rad,k) = crm_rad%ns         (icrm,i_rad,j_rad,k) + micro_field(icrm,i,j,k,ins )
#endif
            !$acc atomic update
            crm_output_gliqwp(icrm,l) = crm_output_gliqwp(icrm,l) + qcl(icrm,i,j,k)
            !$acc atomic update
            crm_output_gicewp(icrm,l) = crm_output_gicewp(icrm,l) + qci(icrm,i,j,k)
          enddo
        enddo
      enddo
    enddo

    ! Diagnose mass fluxes to drive CAM's convective transport of tracers.
    ! definition of mass fluxes is taken from Xu et al., 2002, QJRMS.
    !$acc parallel loop collapse(3) copyin(tabs,pres,qcl,qci,qpl,qpi,rhow,w) copy(mdi_crm,mui_crm) async(asyncid)
    do icrm = 1 , ncrms
      do j=1, ny
        do i=1, nx
          do k=1, nzm+1
            l=plev+1-k+1
            if(w(icrm,i,j,k).gt.0.) then
              kx=max(1, k-1)
              qsat = qsatw_crm(tabs(icrm,i,j,kx),pres(icrm,kx))
              if(qcl(icrm,i,j,kx)+qci(icrm,i,j,kx).gt.min(real(1.e-5,crm_rknd),0.01*qsat)) then
                tmp = rhow(icrm,k)*w(icrm,i,j,k)
                !$acc atomic update
                mui_crm(icrm,l) = mui_crm(icrm,l)+tmp
              endif
            else if (w(icrm,i,j,k).lt.0.) then
              kx=min(k+1, nzm)
              qsat = qsatw_crm(tabs(icrm,i,j,kx),pres(icrm,kx))
              if(qcl(icrm,i,j,kx)+qci(icrm,i,j,kx).gt.min(real(1.e-5,crm_rknd),0.01*qsat)) then
                tmp = rhow(icrm,k)*w(icrm,i,j,k)
                !$acc atomic update
                mdi_crm(icrm,l) = mdi_crm(icrm,l)+tmp
              else if(qpl(icrm,i,j,kx)+qpi(icrm,i,j,kx).gt.1.0e-4) then
                tmp = rhow(icrm,k)*w(icrm,i,j,k)
                !$acc atomic update
                mdi_crm(icrm,l) = mdi_crm(icrm,l)+tmp
              endif
            endif
          enddo
        enddo
      enddo
    enddo

    !$acc parallel loop collapse(3) copyin(cwp,cwph,cwpm,cwpl,cttemp,chtemp,cmtemp,cltemp) copy(crm_output_cltot,crm_output_clhgh,crm_output_clmed,crm_output_cllow) async(asyncid)
    do icrm = 1 , ncrms
      do j=1,ny
        do i=1,nx
          if(cwp (i,j,icrm).gt.cwp_threshold) then
            !$acc atomic update
            crm_output_cltot(icrm) = crm_output_cltot(icrm) + cttemp(i,j,icrm)
          endif
          if(cwph(i,j,icrm).gt.cwp_threshold) then
            !$acc atomic update
            crm_output_clhgh(icrm) = crm_output_clhgh(icrm) + chtemp(i,j,icrm)
          endif
          if(cwpm(i,j,icrm).gt.cwp_threshold) then
            !$acc atomic update
            crm_output_clmed(icrm) = crm_output_clmed(icrm) + cmtemp(i,j,icrm)
          endif
          if(cwpl(i,j,icrm).gt.cwp_threshold) then
            !$acc atomic update
            crm_output_cllow(icrm) = crm_output_cllow(icrm) + cltemp(i,j,icrm)
          endif
        enddo
      enddo
    enddo

  enddo ! nstep

  ! for time-averaging crm output statistics
  factor_xyt = factor_xy / real(nstop,crm_rknd) 

  !========================================================================================
  !----------------------------------------------------------------------------------------
  ! End main time loop
  !----------------------------------------------------------------------------------------
  !========================================================================================

  !$acc exit data copyout(sgs_field,sgs_field_diag,tke2,tk2,tk,tke,tkh,twle,tadv,q0,qpfall,tlat,precflux,prec_xy,fluxtu,fluxtv) async(asyncid)

  !$acc exit data copyout(dudt,dvdt,dwdt,misc,adz,bet,tabs0,qv,qv0,qcl,qci,qn0,qpl,qpi,qp0,tabs,t,micro_field,ttend,qtend,utend,vtend,u,u0,v,v0,w,t0,dz,precsfc,precssfc,rho,qifall,tlatqi) async(asyncid)
  !$acc exit data copyout(sstxy,taux0,tauy0,z,z0,fluxbu,fluxbv,bflx,uhl,vhl,adzw,presi,tkelediss,tkesbdiss,tkesbshear,tkesbbuoy,grdf_x,grdf_y,grdf_z,fcory,fcorzy,ug0,vg0,t01,q01,p0,pres,p) async(asyncid)
  !$acc exit data copyout(rhow,uwle,vwle,uwsb,vwsb,w_max,u_max,dt3,cwp,cwph,cwpm,cwpl,flag_top,cltemp,cmtemp,chtemp,cttemp,mkadv,mkwle,sgsadv,sgswle,gamaz,iw_xy,cw_xy,pw_xy,u200_xy,v200_xy) async(asyncid)
  !$acc exit data copyout(usfc_xy,vsfc_xy,w500_xy,swvp_xy,psfc_xy,u850_xy,v850_xy,cloudtopheight,cloudtoptemp,echotopheight,cld_xy,crm_output_timing_factor,crm_rad_qrad,cf3d) async(asyncid)
  !$acc exit data copyout(crm_output_mcudn,crm_output_mcup,crm_output_cld,crm_output_mcdn,crm_output_gliqwp,crm_output_mcuup,crm_rad_qc,crm_rad_cld,crm_rad_qi,crm_rad_temperature) async(asyncid)
  !$acc exit data copyout(crm_rad_qv,crm_output_gicewp,crm_output_cldtop,mdi_crm,mui_crm,crm_output_cltot,crm_output_clhgh,crm_output_clmed,crm_output_cllow,fluxbt,fluxtt,tdiff,twsb,fzero) async(asyncid)
  !$acc exit data copyout(fluxbq,fluxbmk,fluxtq,fluxtmk,sgswsb,mkdiff,mkwsb,qn,qpsrc,qpevp,accrrc,accrsc,accrsi,accrgi,accrgc,coefice,evapg1,evapg2,evapr1,evaps2,evaps1,evapr2) async(asyncid)

  !$acc wait(asyncid)

  do icrm = 1 , ncrms
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
          colprec = colprec +(qpl(icrm,i,j,k)+qpi(icrm,i,j,k))*crm_input%pdel(icrm,plev-k+1)
          colprecs= colprecs+qpi(icrm,i,j,k)*crm_input%pdel(icrm,plev-k+1)
          tln(l,icrm)  = tln(l,icrm)  +tabs(icrm,i,j,k)
          qln(l,icrm)  = qln(l,icrm)  +qv(icrm,i,j,k)
          qccln(l,icrm)= qccln(l,icrm)+qcl(icrm,i,j,k)
          qiiln(l,icrm)= qiiln(l,icrm)+qci(icrm,i,j,k)
          uln(l,icrm)  = uln(l,icrm)  +u(icrm,i,j,k)
          vln(l,icrm)  = vln(l,icrm)  +v(icrm,i,j,k)

#if defined(SP_ESMT)
          uln_esmt(l,icrm) = uln_esmt(l,icrm)+u_esmt(icrm,i,j,k)
          vln_esmt(l,icrm) = vln_esmt(l,icrm)+v_esmt(icrm,i,j,k)
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
    crm_state%u_wind  (icrm,1:nx,1:ny,1:nzm) = u(icrm,1:nx,1:ny,1:nzm)
    crm_state%v_wind  (icrm,1:nx,1:ny,1:nzm) = v(icrm,1:nx,1:ny,1:nzm)
    crm_state%w_wind  (icrm,1:nx,1:ny,1:nzm) = w(icrm,1:nx,1:ny,1:nzm)
    crm_state%temperature  (icrm,1:nx,1:ny,1:nzm) = tabs(icrm,1:nx,1:ny,1:nzm)

#ifdef m2005
      crm_state%qt(icrm,1:nx,1:ny,1:nzm) = micro_field(icrm,1:nx,1:ny,1:nzm,1 )
      crm_state%nc(icrm,1:nx,1:ny,1:nzm) = micro_field(icrm,1:nx,1:ny,1:nzm,2 )
      crm_state%qr(icrm,1:nx,1:ny,1:nzm) = micro_field(icrm,1:nx,1:ny,1:nzm,3 )
      crm_state%nr(icrm,1:nx,1:ny,1:nzm) = micro_field(icrm,1:nx,1:ny,1:nzm,4 )
      crm_state%qi(icrm,1:nx,1:ny,1:nzm) = micro_field(icrm,1:nx,1:ny,1:nzm,5 )
      crm_state%ni(icrm,1:nx,1:ny,1:nzm) = micro_field(icrm,1:nx,1:ny,1:nzm,6 )
      crm_state%qs(icrm,1:nx,1:ny,1:nzm) = micro_field(icrm,1:nx,1:ny,1:nzm,7 )
      crm_state%ns(icrm,1:nx,1:ny,1:nzm) = micro_field(icrm,1:nx,1:ny,1:nzm,8 )
      crm_state%qg(icrm,1:nx,1:ny,1:nzm) = micro_field(icrm,1:nx,1:ny,1:nzm,9 )
      crm_state%ng(icrm,1:nx,1:ny,1:nzm) = micro_field(icrm,1:nx,1:ny,1:nzm,10)
      crm_state%qc(icrm,1:nx,1:ny,1:nzm) = cloudliq(1:nx,1:ny,1:nzm,icrm)
#else
      crm_state%qt(icrm,1:nx,1:ny,1:nzm) = micro_field(icrm,1:nx,1:ny,1:nzm,1)
      crm_state%qp(icrm,1:nx,1:ny,1:nzm) = micro_field(icrm,1:nx,1:ny,1:nzm,2)
      crm_state%qn(icrm,1:nx,1:ny,1:nzm) = qn(icrm,1:nx,1:ny,1:nzm)
#endif

    crm_output%tk   (icrm,1:nx,1:ny,1:nzm) = tk(icrm,1:nx, 1:ny, 1:nzm)
    crm_output%tkh  (icrm,1:nx,1:ny,1:nzm) = tkh(icrm,1:nx, 1:ny, 1:nzm)
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
        crm_output%qcl(icrm,i,j,k) = qcl(icrm,i,j,k)
        crm_output%qci(icrm,i,j,k) = qci(icrm,i,j,k)
        crm_output%qpl(icrm,i,j,k) = qpl(icrm,i,j,k)
        crm_output%qpi(icrm,i,j,k) = qpi(icrm,i,j,k)
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
          crm_output%qc_mean(icrm,l) = crm_output%qc_mean(icrm,l) + qcl(icrm,i,j,k)
          crm_output%qi_mean(icrm,l) = crm_output%qi_mean(icrm,l) + qci(icrm,i,j,k)
          crm_output%qr_mean(icrm,l) = crm_output%qr_mean(icrm,l) + qpl(icrm,i,j,k)
#ifdef sam1mom
          omg = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(icrm,i,j,k)-tgrmin)*a_gr))
          crm_output%qg_mean(icrm,l) = crm_output%qg_mean(icrm,l) + qpi(icrm,i,j,k)*omg
          crm_output%qs_mean(icrm,l) = crm_output%qs_mean(icrm,l) + qpi(icrm,i,j,k)*(1.-omg)
#else
          crm_output%qg_mean(icrm,l) = crm_output%qg_mean(icrm,l) + micro_field(icrm,i,j,k,iqg)
          crm_output%qs_mean(icrm,l) = crm_output%qs_mean(icrm,l) + micro_field(icrm,i,j,k,iqs)

          crm_output%nc_mean(icrm,l) = crm_output%nc_mean(icrm,l) + micro_field(icrm,i,j,k,incl)
          crm_output%ni_mean(icrm,l) = crm_output%ni_mean(icrm,l) + micro_field(icrm,i,j,k,inci)
          crm_output%nr_mean(icrm,l) = crm_output%nr_mean(icrm,l) + micro_field(icrm,i,j,k,inr )
          crm_output%ng_mean(icrm,l) = crm_output%ng_mean(icrm,l) + micro_field(icrm,i,j,k,ing )
          crm_output%ns_mean(icrm,l) = crm_output%ns_mean(icrm,l) + micro_field(icrm,i,j,k,ins )
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
        precsfc(icrm,i,j) = precsfc(icrm,i,j)*dz(icrm)/dt/dble(nstop)
        precssfc(icrm,i,j) = precssfc(icrm,i,j)*dz(icrm)/dt/dble(nstop)
#endif /* sam1mom */
#ifdef m2005
        precsfc(icrm,i,j) = precsfc(icrm,i,j)*dz(icrm)/dt/dble(nstop)     !mm/s/dz --> mm/s
        precssfc(icrm,i,j) = precssfc(icrm,i,j)*dz(icrm)/dt/dble(nstop)   !mm/s/dz --> mm/s
#endif /* m2005 */
        if(precsfc(icrm,i,j).gt.10./86400.) then
           crm_output%precc (icrm) = crm_output%precc (icrm) + precsfc(icrm,i,j)
           crm_output%precsc(icrm) = crm_output%precsc(icrm) + precssfc(icrm,i,j)
        else
           crm_output%precl (icrm) = crm_output%precl (icrm) + precsfc(icrm,i,j)
           crm_output%precsl(icrm) = crm_output%precsl(icrm) + precssfc(icrm,i,j)
        endif
      enddo
    enddo
    crm_output%prec_crm(icrm,:,:) = precsfc(icrm,:,:)/1000.           !mm/s --> m/s
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
          u2z = u2z+(u(icrm,i,j,k)-u0(icrm,k))**2
          v2z = v2z+(v(icrm,i,j,k)-v0(icrm,k))**2
          w2z = w2z+0.5*(w(icrm,i,j,k+1)**2+w(icrm,i,j,k)**2)
        enddo
      enddo
      !+++mhwang
      ! mkwsb, mkle, mkadv, mkdiff (also crm_output%flux_u, crm_output%flux_v,icrm) seem not calculted correclty in the spcam3.5 codes.
      ! Only values at the last time step are calculated, but is averaged over the entire GCM
      ! time step.
      !---mhwang

      tmp1 = dz(icrm)/rhow(icrm,k)
      tmp2 = tmp1/dtn                        ! dtn is calculated inside of the icyc loop.
                                             ! It seems wrong to use it here ???? +++mhwang
      mkwsb(icrm,k,:) = mkwsb(icrm,k,:) * tmp1*rhow(icrm,k) * factor_xy/nstop     !kg/m3/s --> kg/m2/s
      mkwle(icrm,k,:) = mkwle(icrm,k,:) * tmp2*rhow(icrm,k) * factor_xy/nstop     !kg/m3   --> kg/m2/s
      mkadv(icrm,k,:) = mkadv(icrm,k,:) * factor_xy*icrm_run_time     ! kg/kg  --> kg/kg/s
      mkdiff(icrm,k,:) = mkdiff(icrm,k,:) * factor_xy*icrm_run_time   ! kg/kg  --> kg/kg/s

      ! qpsrc, qpevp, qpfall in M2005 are calculated in micro_flux.
      qpsrc(icrm,k) = qpsrc(icrm,k) * factor_xy*icrm_run_time
      qpevp(icrm,k) = qpevp(icrm,k) * factor_xy*icrm_run_time
      qpfall(icrm,k) = qpfall(icrm,k) * factor_xy*icrm_run_time   ! kg/kg in M2005 ---> kg/kg/s
      precflux(icrm,k) = precflux(icrm,k) * factor_xy*dz(icrm)/dt/nstop  !kg/m2/dz in M2005 -->kg/m2/s or mm/s (idt_gl=1/dt/nstop)

      l = plev-k+1
      crm_output%flux_u    (icrm,l) = (uwle(icrm,k) + uwsb(icrm,k))*tmp1*factor_xy/nstop
      crm_output%flux_v    (icrm,l) = (vwle(icrm,k) + vwsb(icrm,k))*tmp1*factor_xy/nstop
#ifdef sam1mom
      crm_output%flux_qt   (icrm,l) = mkwle(icrm,k,1) + mkwsb(icrm,k,1)
      crm_output%fluxsgs_qt(icrm,l) = mkwsb(icrm,k,1)
      crm_output%flux_qp   (icrm,l) = mkwle(icrm,k,2) + mkwsb(icrm,k,2)
      crm_output%qt_trans  (icrm,l) = mkadv(icrm,k,1) + mkdiff(icrm,k,1)
      crm_output%qp_trans  (icrm,l) = mkadv(icrm,k,2) + mkdiff(icrm,k,2)
#endif /* sam1mom */
#ifdef m2005
      crm_output%flux_qt   (icrm,l) = mkwle(icrm,k,1   ) + mkwsb(icrm,k,1   ) +  &
                         mkwle(icrm,k,iqci) + mkwsb(icrm,k,iqci)
      crm_output%fluxsgs_qt(icrm,l) = mkwsb(icrm,k,1   ) + mkwsb(icrm,k,iqci)
      crm_output%flux_qp   (icrm,l) = mkwle(icrm,k,iqr) + mkwsb(icrm,k,iqr) +  &
                         mkwle(icrm,k,iqs) + mkwsb(icrm,k,iqs) + mkwle(icrm,k,iqg) + mkwsb(icrm,k,iqg)
      crm_output%qt_trans  (icrm,l) = mkadv(icrm,k,1) + mkadv(icrm,k,iqci) + &
                         mkdiff(icrm,k,1) + mkdiff(icrm,k,iqci)
      crm_output%qp_trans  (icrm,l) = mkadv(icrm,k,iqr) + mkadv(icrm,k,iqs) + mkadv(icrm,k,iqg) + &
                         mkdiff(icrm,k,iqr) + mkdiff(icrm,k,iqs) + mkdiff(icrm,k,iqg)
#endif /* m2005 */
      crm_output%tkesgsz   (icrm,l)= rho(icrm,k)*sum(tke(icrm,1:nx,1:ny,k))*factor_xy
      crm_output%tkez      (icrm,l)= rho(icrm,k)*0.5*(u2z+v2z*YES3D+w2z)*factor_xy + crm_output%tkesgsz(icrm,l)
      crm_output%tkz       (icrm,l) = sum(tk(icrm,1:nx, 1:ny, k)) * factor_xy
      crm_output%precflux      (icrm,l) = precflux(icrm,k)/1000.       !mm/s  -->m/s

      crm_output%qp_fall   (icrm,l) = qpfall(icrm,k)
      crm_output%qp_evp    (icrm,l) = qpevp(icrm,k)
      crm_output%qp_src    (icrm,l) = qpsrc(icrm,k)

      crm_output%qt_ls     (icrm,l) = qtend(icrm,k)
      crm_output%t_ls      (icrm,l) = ttend(icrm,k)
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
