module vars
  use grid
  use params, only: crm_rknd
#ifdef CRM
#ifdef MODAL_AERO
  use modal_aero_data,   only: ntot_amode
#endif
#endif

  use openacc_pool
  implicit none
  !--------------------------------------------------------------------
  ! prognostic variables:

  real(crm_rknd), pointer :: u  (:,:,:,:) !REDIM ! x-wind
  real(crm_rknd), pointer :: v  (:,:,:,:) !REDIM ! y-wind
  real(crm_rknd), pointer :: w  (:,:,:,:) !REDIM ! z-wind
  real(crm_rknd), pointer :: t  (:,:,:,:) !REDIM ! liquid/ice water static energy

  !--------------------------------------------------------------------
  ! diagnostic variables:

  real(crm_rknd), pointer :: p   (:,:,:,:) !REDIM       ! perturbation pressure (from Poison eq)
  real(crm_rknd), pointer :: tabs(:,:,:,:) !REDIM        ! temperature
  real(crm_rknd), pointer :: qv  (:,:,:,:) !REDIM       ! water vapor
  real(crm_rknd), pointer :: qcl (:,:,:,:) !REDIM       ! liquid water  (condensate)
  real(crm_rknd), pointer :: qpl (:,:,:,:) !REDIM       ! liquid water  (precipitation)
  real(crm_rknd), pointer :: qci (:,:,:,:) !REDIM       ! ice water  (condensate)
  real(crm_rknd), pointer :: qpi (:,:,:,:) !REDIM       ! ice water  (precipitation)

  real(crm_rknd), pointer :: tke2(:,:,:,:) !REDIM      ! SGS TKE
  real(crm_rknd), pointer :: tk2 (:,:,:,:) !REDIM    ! SGS eddyviscosity

  !--------------------------------------------------------------------
  ! time-tendencies for prognostic variables

  real(crm_rknd), pointer :: dudt(:,:,:,:,:) !REDIM
  real(crm_rknd), pointer :: dvdt(:,:,:,:,:) !REDIM
  real(crm_rknd), pointer :: dwdt(:,:,:,:,:) !REDIM

  !----------------------------------------------------------------
  ! Temporary storage array:

  real(crm_rknd), pointer :: misc(:,:,:,:) !REDIM
  !------------------------------------------------------------------
  ! fluxes at the top and bottom of the domain:

  real(crm_rknd), pointer :: fluxbu  (:,:,:) !REDIM
  real(crm_rknd), pointer :: fluxbv  (:,:,:) !REDIM
  real(crm_rknd), pointer :: fluxbt  (:,:,:) !REDIM
  real(crm_rknd), pointer :: fluxbq  (:,:,:) !REDIM
  real(crm_rknd), pointer :: fluxtu  (:,:,:) !REDIM
  real(crm_rknd), pointer :: fluxtv  (:,:,:) !REDIM
  real(crm_rknd), pointer :: fluxtt  (:,:,:) !REDIM
  real(crm_rknd), pointer :: fluxtq  (:,:,:) !REDIM
  real(crm_rknd), pointer :: fzero   (:,:,:) !REDIM
  real(crm_rknd), pointer :: precsfc (:,:,:) !REDIM ! surface precip. rate
  real(crm_rknd), pointer :: precssfc(:,:,:) !REDIM ! surface ice precip. rate

  !-----------------------------------------------------------------
  ! profiles

  real(crm_rknd), pointer :: t0   (:,:) !REDIM
  real(crm_rknd), pointer :: q0   (:,:) !REDIM
  real(crm_rknd), pointer :: qv0  (:,:) !REDIM
  real(crm_rknd), pointer :: tabs0(:,:) !REDIM
  real(crm_rknd), pointer :: tl0  (:,:) !REDIM
  real(crm_rknd), pointer :: tv0  (:,:) !REDIM
  real(crm_rknd), pointer :: u0   (:,:) !REDIM
  real(crm_rknd), pointer :: v0   (:,:) !REDIM
  real(crm_rknd), pointer :: tg0  (:,:) !REDIM
  real(crm_rknd), pointer :: qg0  (:,:) !REDIM
  real(crm_rknd), pointer :: ug0  (:,:) !REDIM
  real(crm_rknd), pointer :: vg0  (:,:) !REDIM
  real(crm_rknd), pointer :: p0   (:,:) !REDIM
  real(crm_rknd), pointer :: tke0 (:,:) !REDIM
  real(crm_rknd), pointer :: t01  (:,:) !REDIM
  real(crm_rknd), pointer :: q01  (:,:) !REDIM
  real(crm_rknd), pointer :: qp0  (:,:) !REDIM
  real(crm_rknd), pointer :: qn0  (:,:) !REDIM

  !-----------------------------------------------------------------
  ! reference vertical profiles:

  real(crm_rknd), pointer :: prespot(:,:) !REDIM ! (1000./pres)**R/cp
  real(crm_rknd), pointer :: rho    (:,:) !REDIM ! air density at pressure levels,kg/m3
  real(crm_rknd), pointer :: rhow   (:,:) !REDIM ! air density at vertical velocity levels,kg/m3
  real(crm_rknd), pointer :: bet    (:,:) !REDIM ! = ggr/tv0
  real(crm_rknd), pointer :: gamaz  (:,:) !REDIM ! ggr/cp*z
  real(crm_rknd), pointer :: wsub   (:,:) !REDIM ! Large-scale subsidence velocity,m/s
  real(crm_rknd), pointer :: qtend  (:,:) !REDIM ! Large-scale tendency for total water
  real(crm_rknd), pointer :: ttend  (:,:) !REDIM ! Large-scale tendency for temp.
  real(crm_rknd), pointer :: utend  (:,:) !REDIM ! Large-scale tendency for u
  real(crm_rknd), pointer :: vtend  (:,:) !REDIM ! Large-scale tendency for v

  !---------------------------------------------------------------------
  !  Horizontally varying stuff (as a function of xy)
  !
  real(crm_rknd), pointer :: sstxy    (:,:,:) !REDIM  !  surface temperature xy-distribution
  real(crm_rknd), pointer :: fcory    (:,:)   !REDIM  !  Coriolis parameter xy-distribution
  real(crm_rknd), pointer :: fcorzy   (:,:)   !REDIM  !  z-Coriolis parameter xy-distribution
  real(crm_rknd), pointer :: latitude (:,:,:) !REDIM  ! latitude (icrm,degrees)
  real(crm_rknd), pointer :: longitude(:,:,:) !REDIM  ! longitude(icrm,degrees)
  real(crm_rknd), pointer :: prec_xy  (:,:,:) !REDIM  ! mean precip. rate for outout
  real(crm_rknd), pointer :: pw_xy    (:,:,:) !REDIM  ! precipitable water
  real(crm_rknd), pointer :: cw_xy    (:,:,:) !REDIM  ! cloud water path
  real(crm_rknd), pointer :: iw_xy    (:,:,:) !REDIM  ! ice water path
  real(crm_rknd), pointer :: cld_xy   (:,:,:) !REDIM  ! cloud frequency
  real(crm_rknd), pointer :: u200_xy  (:,:,:) !REDIM  ! u-wind at 200 mb
  real(crm_rknd), pointer :: usfc_xy  (:,:,:) !REDIM  ! u-wind at at the surface
  real(crm_rknd), pointer :: v200_xy  (:,:,:) !REDIM  ! v-wind at 200 mb
  real(crm_rknd), pointer :: vsfc_xy  (:,:,:) !REDIM  ! v-wind at the surface
  real(crm_rknd), pointer :: w500_xy  (:,:,:) !REDIM  ! w at 500 mb

  !----------------------------------------------------------------------
  !  Vertical profiles of quantities sampled for statitistics purposes:

  real(crm_rknd), pointer :: twle     (:,:) !REDIM
  real(crm_rknd), pointer :: twsb     (:,:) !REDIM
  real(crm_rknd), pointer :: precflux (:,:) !REDIM
  real(crm_rknd), pointer :: uwle     (:,:) !REDIM
  real(crm_rknd), pointer :: uwsb     (:,:) !REDIM
  real(crm_rknd), pointer :: vwle     (:,:) !REDIM
  real(crm_rknd), pointer :: vwsb     (:,:) !REDIM
  real(crm_rknd), pointer :: tkelediss(:,:) !REDIM
  real(crm_rknd), pointer :: t2leadv  (:,:) !REDIM
  real(crm_rknd), pointer :: t2legrad (:,:) !REDIM
  real(crm_rknd), pointer :: t2lediff (:,:) !REDIM
  real(crm_rknd), pointer :: t2lediss (:,:) !REDIM
  real(crm_rknd), pointer :: twleadv  (:,:) !REDIM
  real(crm_rknd), pointer :: twlediff (:,:) !REDIM
  real(crm_rknd), pointer :: tadv     (:,:) !REDIM
  real(crm_rknd), pointer :: tdiff    (:,:) !REDIM
  real(crm_rknd), pointer :: tlat     (:,:) !REDIM
  real(crm_rknd), pointer :: tlatqi   (:,:) !REDIM
  real(crm_rknd), pointer :: qifall   (:,:) !REDIM
  real(crm_rknd), pointer :: qpfall   (:,:) !REDIM
  real(crm_rknd), pointer :: w_max    (:)   !REDIM
  real(crm_rknd), pointer :: u_max    (:)   !REDIM


  ! register functions:


  !real(crm_rknd), external :: esatw_crm,esati_crm,dtesatw_crm,dtesati_crm
  !real(crm_rknd), external :: qsatw_crm,qsati_crm,dtqsatw_crm,dtqsati_crm
  !integer, external :: lenstr, bytes_in_rec

  ! energy conservation diagnostics:

  real(8), pointer :: total_water_before (:) !REDIM
  real(8), pointer :: total_water_after  (:) !REDIM
  real(8), pointer :: total_water_evap   (:) !REDIM
  real(8), pointer :: total_water_prec   (:) !REDIM
  real(8), pointer :: total_water_ls     (:) !REDIM
  real(8), pointer :: total_water_clubb  (:) !REDIM
  real(8), pointer :: total_energy_before(:) !REDIM
  real(8), pointer :: total_energy_after (:) !REDIM
  real(8), pointer :: total_energy_evap  (:) !REDIM
  real(8), pointer :: total_energy_prec  (:) !REDIM
  real(8), pointer :: total_energy_ls    (:) !REDIM
  real(8), pointer :: total_energy_clubb (:) !REDIM
  real(8), pointer :: total_energy_rad   (:) !REDIM
  real(8), pointer :: qtotmicro          (:,:) !REDIM  ! total water for water conservation test in microphysics +++mhwang


  !===========================================================================
  ! UW ADDITIONS

  ! conditional average statistics, subsumes cloud_factor, core_factor, coredn_factor
  real(crm_rknd), pointer :: CF3D(:,:,:,:) !REDIM  ! Cloud fraction
  ! =1.0 when there is no fractional cloudiness scheme
  ! = cloud fraction produced by fractioal cloudiness scheme when avaiable

  ! 850 mbar horizontal winds
  real(crm_rknd), pointer :: u850_xy(:,:,:) !REDIM ! zonal velocity at 850 mb
  real(crm_rknd), pointer :: v850_xy(:,:,:) !REDIM ! meridional velocity at 850 mb

  ! Surface pressure
  real(crm_rknd), pointer :: psfc_xy(:,:,:) !REDIM ! pressure (in millibar) at lowest grid point

  ! Saturated water vapor path, useful for computing column relative humidity
  real(crm_rknd), pointer :: swvp_xy(:,:,:) !REDIM  ! saturated water vapor path (wrt water)

  ! Cloud and echo top heights, and cloud top temperature (instantaneous)
  real(crm_rknd), pointer :: cloudtopheight(:,:,:) !REDIM
  real(crm_rknd), pointer :: echotopheight (:,:,:) !REDIM
  real(crm_rknd), pointer :: cloudtoptemp  (:,:,:) !REDIM

  ! END UW ADDITIONS
#if (defined CRM && defined MODAL_AERO)
  real(crm_rknd), pointer :: naer (:,:,:) !REDIM    ! Aerosol number concentration [/m3]
  real(crm_rknd), pointer :: vaer (:,:,:) !REDIM    ! aerosol volume concentration [m3/m3]
  real(crm_rknd), pointer :: hgaer(:,:,:) !REDIM    ! hygroscopicity of aerosol mode
#endif

integer :: ncondavg
integer :: icondavg_cld
integer :: icondavg_cor
integer :: icondavg_cordn
integer :: icondavg_satdn
integer :: icondavg_satup
integer :: icondavg_env
real(crm_rknd), pointer :: condavg_factor(:,:) ! replaces cloud_factor, core_factor
real(crm_rknd), pointer :: condavg_mask(:,:,:,:) ! indicator array for various conditional averages
character(LEN=8), pointer :: condavgname(:) ! array of short names
character(LEN=25), pointer :: condavglongname(:) ! array of long names

contains

  subroutine allocate_vars(ncrms)
    use openacc_utils
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) :: zero
    ! allocate( u (ncrms,dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm) ) ! x-wind
    ! allocate( v (ncrms,dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm) ) ! y-wind
    ! allocate( w (ncrms,dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz ) ) ! z-wind
    ! allocate( t (ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ) ! liquid/ice water static energy
    ! allocate( p   (ncrms,0:nx, (1-YES3D):ny, nzm) )     ! perturbation pressure (from Poison eq)
    ! allocate( tke2(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) )   ! SGS TKE
    ! allocate( tk2 (ncrms,0:nxp1, (1-YES3D):nyp1, nzm) ) ! SGS eddyviscosity
    ! allocate( sstxy    (ncrms,0:nx,(1-YES3D):ny) ) !  surface temperature xy-distribution
    ! allocate( fcory    (ncrms,0:ny) )             !  Coriolis parameter xy-distribution
    ! allocate( tabs               (ncrms,nx, ny, nzm) )                 ! temperature
    ! allocate( qv                 (ncrms,nx, ny, nzm) )                ! water vapor
    ! allocate( qcl                (ncrms,nx, ny, nzm) )                ! liquid water  (condensate)
    ! allocate( qpl                (ncrms,nx, ny, nzm) )                ! liquid water  (precipitation)
    ! allocate( qci                (ncrms,nx, ny, nzm) )                ! ice water  (condensate)
    ! allocate( qpi                (ncrms,nx, ny, nzm) )                ! ice water  (precipitation)
    ! allocate( dudt               (ncrms,nxp1, ny, nzm, 3) )
    ! allocate( dvdt               (ncrms,nx, nyp1, nzm, 3) )
    ! allocate( dwdt               (ncrms,nx, ny, nz,  3) )
    ! allocate( misc               (ncrms,nx, ny, nz) )
    ! allocate( fluxbu             (ncrms,nx, ny) )
    ! allocate( fluxbv             (ncrms,nx, ny) )
    ! allocate( fluxbt             (ncrms,nx, ny) )
    ! allocate( fluxbq             (ncrms,nx, ny) )
    ! allocate( fluxtu             (ncrms,nx, ny) )
    ! allocate( fluxtv             (ncrms,nx, ny) )
    ! allocate( fluxtt             (ncrms,nx, ny) )
    ! allocate( fluxtq             (ncrms,nx, ny) )
    ! allocate( fzero              (ncrms,nx, ny) )
    ! allocate( precsfc            (ncrms,nx,ny) ) ! surface precip. rate
    ! allocate( precssfc           (ncrms,nx,ny) ) ! surface ice precip. rate
    ! allocate( t0                 (ncrms,nzm) )
    ! allocate( q0                 (ncrms,nzm) )
    ! allocate( qv0                (ncrms,nzm) )
    ! allocate( tabs0              (ncrms,nzm) )
    ! allocate( tl0                (ncrms,nzm) )
    ! allocate( tv0                (ncrms,nzm) )
    ! allocate( u0                 (ncrms,nzm) )
    ! allocate( v0                 (ncrms,nzm) )
    ! allocate( tg0                (ncrms,nzm) )
    ! allocate( qg0                (ncrms,nzm) )
    ! allocate( ug0                (ncrms,nzm) )
    ! allocate( vg0                (ncrms,nzm) )
    ! allocate( p0                 (ncrms,nzm) )
    ! allocate( tke0               (ncrms,nzm) )
    ! allocate( t01                (ncrms,nzm) )
    ! allocate( q01                (ncrms,nzm) )
    ! allocate( qp0                (ncrms,nzm) )
    ! allocate( qn0                (ncrms,nzm) )
    ! allocate( prespot            (ncrms,nzm) ) ! (1000./pres)**R/cp
    ! allocate( rho                (ncrms,nzm) ) ! air density at pressure levels,kg/m3
    ! allocate( rhow               (ncrms,nz ) ) ! air density at vertical velocity levels,kg/m3
    ! allocate( bet                (ncrms,nzm) ) ! = ggr/tv0
    ! allocate( gamaz              (ncrms,nzm) ) ! ggr/cp*z
    ! allocate( wsub               (ncrms,nz ) ) ! Large-scale subsidence velocity,m/s
    ! allocate( qtend              (ncrms,nzm) ) ! Large-scale tendency for total water
    ! allocate( ttend              (ncrms,nzm) ) ! Large-scale tendency for temp.
    ! allocate( utend              (ncrms,nzm) ) ! Large-scale tendency for u
    ! allocate( vtend              (ncrms,nzm) ) ! Large-scale tendency for v
    ! allocate( fcorzy             (ncrms,ny) )               !  z-Coriolis parameter xy-distribution
    ! allocate( latitude           (ncrms,nx,ny) )            ! latitude (icrm,degrees)
    ! allocate( longitude          (ncrms,nx,ny) )            ! longitude(icrm,degrees)
    ! allocate( prec_xy            (ncrms,nx,ny) )            ! mean precip. rate for outout
    ! allocate( pw_xy              (ncrms,nx,ny) )            ! precipitable water
    ! allocate( cw_xy              (ncrms,nx,ny) )            ! cloud water path
    ! allocate( iw_xy              (ncrms,nx,ny) )            ! ice water path
    ! allocate( cld_xy             (ncrms,nx,ny) )            ! cloud frequency
    ! allocate( u200_xy            (ncrms,nx,ny) )            ! u-wind at 200 mb
    ! allocate( usfc_xy            (ncrms,nx,ny) )            ! u-wind at at the surface
    ! allocate( v200_xy            (ncrms,nx,ny) )            ! v-wind at 200 mb
    ! allocate( vsfc_xy            (ncrms,nx,ny) )            ! v-wind at the surface
    ! allocate( w500_xy            (ncrms,nx,ny) )            ! w at 500 mb
    ! allocate( twle               (ncrms,nz) )
    ! allocate( twsb               (ncrms,nz) )
    ! allocate( precflux           (ncrms,nz) )
    ! allocate( uwle               (ncrms,nz) )
    ! allocate( uwsb               (ncrms,nz) )
    ! allocate( vwle               (ncrms,nz) )
    ! allocate( vwsb               (ncrms,nz) )
    ! allocate( tkelediss          (ncrms,nz) )
    ! allocate( t2leadv            (ncrms,nz) )
    ! allocate( t2legrad           (ncrms,nz) )
    ! allocate( t2lediff           (ncrms,nz) )
    ! allocate( t2lediss           (ncrms,nz) )
    ! allocate( twleadv            (ncrms,nz) )
    ! allocate( twlediff           (ncrms,nz) )
    ! allocate( tadv               (ncrms,nz) )
    ! allocate( tdiff              (ncrms,nz) )
    ! allocate( tlat               (ncrms,nz) )
    ! allocate( tlatqi             (ncrms,nz) )
    ! allocate( qifall             (ncrms,nz) )
    ! allocate( qpfall             (ncrms,nz) )
    ! allocate( w_max              (ncrms) )
    ! allocate( u_max              (ncrms) )
    ! allocate( total_water_before (ncrms) )
    ! allocate( total_water_after  (ncrms) )
    ! allocate( total_water_evap   (ncrms) )
    ! allocate( total_water_prec   (ncrms) )
    ! allocate( total_water_ls     (ncrms) )
    ! allocate( total_water_clubb  (ncrms) )
    ! allocate( total_energy_before(ncrms) )
    ! allocate( total_energy_after (ncrms) )
    ! allocate( total_energy_evap  (ncrms) )
    ! allocate( total_energy_prec  (ncrms) )
    ! allocate( total_energy_ls    (ncrms) )
    ! allocate( total_energy_clubb (ncrms) )
    ! allocate( total_energy_rad   (ncrms) )
    ! allocate( qtotmicro          (ncrms,5) )
    ! allocate( CF3D               (ncrms,1:nx, 1:ny, 1:nzm) )
    ! allocate( u850_xy            (ncrms,nx,ny) )
    ! allocate( v850_xy            (ncrms,nx,ny) )
    ! allocate( psfc_xy            (ncrms,nx,ny) )
    ! allocate( swvp_xy            (ncrms,nx,ny) )
    ! allocate( cloudtopheight     (ncrms,nx,ny) )
    ! allocate( echotopheight      (ncrms,nx,ny) )
    ! allocate( cloudtoptemp       (ncrms,nx,ny) )
#if (defined CRM && defined MODAL_AERO)
    ! allocate( naer          (ncrms,nzm, ntot_amode) )
    ! allocate( vaer          (ncrms,nzm, ntot_amode) )
    ! allocate( hgaer         (ncrms,nzm, ntot_amode) )
#endif

    call pool_push( u     , (/1,dimx1_u,dimy1_u,1/),(/ncrms,dimx2_u,dimy2_u,nzm/))
    call pool_push( v     , (/1,dimx1_v,dimy1_v,1/),(/ncrms,dimx2_v,dimy2_v,nzm/))
    call pool_push( w     , (/1,dimx1_w,dimy1_w,1/),(/ncrms,dimx2_w,dimy2_w,nz /))
    call pool_push( t     , (/1,dimx1_s,dimy1_s,1/),(/ncrms,dimx2_s,dimy2_s,nzm/))
    call pool_push( tke2  , (/1,dimx1_s,dimy1_s,1/),(/ncrms,dimx2_s,dimy2_s,nzm/))
    call pool_push( p     , (/1,0,(1-YES3D),1/),(/ncrms,nx  ,ny  ,nzm/))
    call pool_push( tk2   , (/1,0,(1-YES3D),1/),(/ncrms,nxp1,nyp1,nzm/))
    call pool_push( sstxy , (/1,0,(1-YES3D)/),(/ncrms,nx,ny/))
    call pool_push( fcory , (/1,0/),(/ncrms,ny/))
    call pool_push( tabs                , (/ncrms,nx, ny, nzm/) )
    call pool_push( qv                  , (/ncrms,nx, ny, nzm/) )
    call pool_push( qcl                 , (/ncrms,nx, ny, nzm/) )
    call pool_push( qpl                 , (/ncrms,nx, ny, nzm/) )
    call pool_push( qci                 , (/ncrms,nx, ny, nzm/) )
    call pool_push( qpi                 , (/ncrms,nx, ny, nzm/) )
    call pool_push( dudt                , (/ncrms,nxp1, ny, nzm, 3/) )
    call pool_push( dvdt                , (/ncrms,nx, nyp1, nzm, 3/) )
    call pool_push( dwdt                , (/ncrms,nx, ny, nz,  3/) )
    call pool_push( misc                , (/ncrms,nx, ny, nz/) )
    call pool_push( fluxbu              , (/ncrms,nx, ny/) )
    call pool_push( fluxbv              , (/ncrms,nx, ny/) )
    call pool_push( fluxbt              , (/ncrms,nx, ny/) )
    call pool_push( fluxbq              , (/ncrms,nx, ny/) )
    call pool_push( fluxtu              , (/ncrms,nx, ny/) )
    call pool_push( fluxtv              , (/ncrms,nx, ny/) )
    call pool_push( fluxtt              , (/ncrms,nx, ny/) )
    call pool_push( fluxtq              , (/ncrms,nx, ny/) )
    call pool_push( fzero               , (/ncrms,nx, ny/) )
    call pool_push( precsfc             , (/ncrms,nx,ny/) )
    call pool_push( precssfc            , (/ncrms,nx,ny/) )
    call pool_push( t0                  , (/ncrms,nzm/) )
    call pool_push( q0                  , (/ncrms,nzm/) )
    call pool_push( qv0                 , (/ncrms,nzm/) )
    call pool_push( tabs0               , (/ncrms,nzm/) )
    call pool_push( tl0                 , (/ncrms,nzm/) )
    call pool_push( tv0                 , (/ncrms,nzm/) )
    call pool_push( u0                  , (/ncrms,nzm/) )
    call pool_push( v0                  , (/ncrms,nzm/) )
    call pool_push( tg0                 , (/ncrms,nzm/) )
    call pool_push( qg0                 , (/ncrms,nzm/) )
    call pool_push( ug0                 , (/ncrms,nzm/) )
    call pool_push( vg0                 , (/ncrms,nzm/) )
    call pool_push( p0                  , (/ncrms,nzm/) )
    call pool_push( tke0                , (/ncrms,nzm/) )
    call pool_push( t01                 , (/ncrms,nzm/) )
    call pool_push( q01                 , (/ncrms,nzm/) )
    call pool_push( qp0                 , (/ncrms,nzm/) )
    call pool_push( qn0                 , (/ncrms,nzm/) )
    call pool_push( prespot             , (/ncrms,nzm/) )
    call pool_push( rho                 , (/ncrms,nzm/) )
    call pool_push( rhow                , (/ncrms,nz /) )
    call pool_push( bet                 , (/ncrms,nzm/) )
    call pool_push( gamaz               , (/ncrms,nzm/) )
    call pool_push( wsub                , (/ncrms,nz /) )
    call pool_push( qtend               , (/ncrms,nzm/) )
    call pool_push( ttend               , (/ncrms,nzm/) )
    call pool_push( utend               , (/ncrms,nzm/) )
    call pool_push( vtend               , (/ncrms,nzm/) )
    call pool_push( fcorzy              , (/ncrms,ny/) )
    call pool_push( latitude            , (/ncrms,nx,ny/) )
    call pool_push( longitude           , (/ncrms,nx,ny/) )
    call pool_push( prec_xy             , (/ncrms,nx,ny/) )
    call pool_push( pw_xy               , (/ncrms,nx,ny/) )
    call pool_push( cw_xy               , (/ncrms,nx,ny/) )
    call pool_push( iw_xy               , (/ncrms,nx,ny/) )
    call pool_push( cld_xy              , (/ncrms,nx,ny/) )
    call pool_push( u200_xy             , (/ncrms,nx,ny/) )
    call pool_push( usfc_xy             , (/ncrms,nx,ny/) )
    call pool_push( v200_xy             , (/ncrms,nx,ny/) )
    call pool_push( vsfc_xy             , (/ncrms,nx,ny/) )
    call pool_push( w500_xy             , (/ncrms,nx,ny/) )
    call pool_push( twle                , (/ncrms,nz/) )
    call pool_push( twsb                , (/ncrms,nz/) )
    call pool_push( precflux            , (/ncrms,nz/) )
    call pool_push( uwle                , (/ncrms,nz/) )
    call pool_push( uwsb                , (/ncrms,nz/) )
    call pool_push( vwle                , (/ncrms,nz/) )
    call pool_push( vwsb                , (/ncrms,nz/) )
    call pool_push( tkelediss           , (/ncrms,nz/) )
    call pool_push( t2leadv             , (/ncrms,nz/) )
    call pool_push( t2legrad            , (/ncrms,nz/) )
    call pool_push( t2lediff            , (/ncrms,nz/) )
    call pool_push( t2lediss            , (/ncrms,nz/) )
    call pool_push( twleadv             , (/ncrms,nz/) )
    call pool_push( twlediff            , (/ncrms,nz/) )
    call pool_push( tadv                , (/ncrms,nz/) )
    call pool_push( tdiff               , (/ncrms,nz/) )
    call pool_push( tlat                , (/ncrms,nz/) )
    call pool_push( tlatqi              , (/ncrms,nz/) )
    call pool_push( qifall              , (/ncrms,nz/) )
    call pool_push( qpfall              , (/ncrms,nz/) )
    call pool_push( w_max               , (/ncrms/) )
    call pool_push( u_max               , (/ncrms/) )
    call pool_push( total_water_before  , (/ncrms/) )
    call pool_push( total_water_after   , (/ncrms/) )
    call pool_push( total_water_evap    , (/ncrms/) )
    call pool_push( total_water_prec    , (/ncrms/) )
    call pool_push( total_water_ls      , (/ncrms/) )
    call pool_push( total_water_clubb   , (/ncrms/) )
    call pool_push( total_energy_before , (/ncrms/) )
    call pool_push( total_energy_after  , (/ncrms/) )
    call pool_push( total_energy_evap   , (/ncrms/) )
    call pool_push( total_energy_prec   , (/ncrms/) )
    call pool_push( total_energy_ls     , (/ncrms/) )
    call pool_push( total_energy_clubb  , (/ncrms/) )
    call pool_push( total_energy_rad    , (/ncrms/) )
    call pool_push( qtotmicro           , (/ncrms,5/) )
    call pool_push( CF3D                , (/ncrms,nx, ny, nzm/) )
    call pool_push( u850_xy             , (/ncrms,nx,ny/) )
    call pool_push( v850_xy             , (/ncrms,nx,ny/) )
    call pool_push( psfc_xy             , (/ncrms,nx,ny/) )
    call pool_push( swvp_xy             , (/ncrms,nx,ny/) )
    call pool_push( cloudtopheight      , (/ncrms,nx,ny/) )
    call pool_push( echotopheight       , (/ncrms,nx,ny/) )
    call pool_push( cloudtoptemp        , (/ncrms,nx,ny/) )
#if (defined CRM && defined MODAL_AERO)
    call pool_push( naer  , (/ncrms,nzm, ntot_amode/) )
    call pool_push( vaer  , (/ncrms,nzm, ntot_amode/) )
    call pool_push( hgaer , (/ncrms,nzm, ntot_amode/) )
#endif

    zero = 0.

    call memzero_crm_rknd( u                    , product(shape(u                   )) )
    call memzero_crm_rknd( v                    , product(shape(v                   )) )
    call memzero_crm_rknd( w                    , product(shape(w                   )) )
    call memzero_crm_rknd( t                    , product(shape(t                   )) )
    call memzero_crm_rknd( p                    , product(shape(p                   )) )
    call memzero_crm_rknd( tabs                 , product(shape(tabs                )) )
    call memzero_crm_rknd( qv                   , product(shape(qv                  )) )
    call memzero_crm_rknd( qcl                  , product(shape(qcl                 )) )
    call memzero_crm_rknd( qpl                  , product(shape(qpl                 )) )
    call memzero_crm_rknd( qci                  , product(shape(qci                 )) )
    call memzero_crm_rknd( qpi                  , product(shape(qpi                 )) )
    call memzero_crm_rknd( tke2                 , product(shape(tke2                )) )
    call memzero_crm_rknd( tk2                  , product(shape(tk2                 )) )
    call memzero_crm_rknd( dudt                 , product(shape(dudt                )) )
    call memzero_crm_rknd( dvdt                 , product(shape(dvdt                )) )
    call memzero_crm_rknd( dwdt                 , product(shape(dwdt                )) )
    call memzero_crm_rknd( misc                 , product(shape(misc                )) )
    call memzero_crm_rknd( fluxbu               , product(shape(fluxbu              )) )
    call memzero_crm_rknd( fluxbv               , product(shape(fluxbv              )) )
    call memzero_crm_rknd( fluxbt               , product(shape(fluxbt              )) )
    call memzero_crm_rknd( fluxbq               , product(shape(fluxbq              )) )
    call memzero_crm_rknd( fluxtu               , product(shape(fluxtu              )) )
    call memzero_crm_rknd( fluxtv               , product(shape(fluxtv              )) )
    call memzero_crm_rknd( fluxtt               , product(shape(fluxtt              )) )
    call memzero_crm_rknd( fluxtq               , product(shape(fluxtq              )) )
    call memzero_crm_rknd( fzero                , product(shape(fzero               )) )
    call memzero_crm_rknd( precsfc              , product(shape(precsfc             )) )
    call memzero_crm_rknd( precssfc             , product(shape(precssfc            )) )
    call memzero_crm_rknd( t0                   , product(shape(t0                  )) )
    call memzero_crm_rknd( q0                   , product(shape(q0                  )) )
    call memzero_crm_rknd( qv0                  , product(shape(qv0                 )) )
    call memzero_crm_rknd( tabs0                , product(shape(tabs0               )) )
    call memzero_crm_rknd( tl0                  , product(shape(tl0                 )) )
    call memzero_crm_rknd( tv0                  , product(shape(tv0                 )) )
    call memzero_crm_rknd( u0                   , product(shape(u0                  )) )
    call memzero_crm_rknd( v0                   , product(shape(v0                  )) )
    call memzero_crm_rknd( tg0                  , product(shape(tg0                 )) )
    call memzero_crm_rknd( qg0                  , product(shape(qg0                 )) )
    call memzero_crm_rknd( ug0                  , product(shape(ug0                 )) )
    call memzero_crm_rknd( vg0                  , product(shape(vg0                 )) )
    call memzero_crm_rknd( p0                   , product(shape(p0                  )) )
    call memzero_crm_rknd( tke0                 , product(shape(tke0                )) )
    call memzero_crm_rknd( t01                  , product(shape(t01                 )) )
    call memzero_crm_rknd( q01                  , product(shape(q01                 )) )
    call memzero_crm_rknd( qp0                  , product(shape(qp0                 )) )
    call memzero_crm_rknd( qn0                  , product(shape(qn0                 )) )
    call memzero_crm_rknd( prespot              , product(shape(prespot             )) )
    call memzero_crm_rknd( rho                  , product(shape(rho                 )) )
    call memzero_crm_rknd( rhow                 , product(shape(rhow                )) )
    call memzero_crm_rknd( bet                  , product(shape(bet                 )) )
    call memzero_crm_rknd( gamaz                , product(shape(gamaz               )) )
    call memzero_crm_rknd( wsub                 , product(shape(wsub                )) )
    call memzero_crm_rknd( qtend                , product(shape(qtend               )) )
    call memzero_crm_rknd( ttend                , product(shape(ttend               )) )
    call memzero_crm_rknd( utend                , product(shape(utend               )) )
    call memzero_crm_rknd( vtend                , product(shape(vtend               )) )
    call memzero_crm_rknd( sstxy                , product(shape(sstxy               )) )
    call memzero_crm_rknd( fcory                , product(shape(fcory               )) )
    call memzero_crm_rknd( fcorzy               , product(shape(fcorzy              )) )
    call memzero_crm_rknd( latitude             , product(shape(latitude            )) )
    call memzero_crm_rknd( longitude            , product(shape(longitude           )) )
    call memzero_crm_rknd( prec_xy              , product(shape(prec_xy             )) )
    call memzero_crm_rknd( pw_xy                , product(shape(pw_xy               )) )
    call memzero_crm_rknd( cw_xy                , product(shape(cw_xy               )) )
    call memzero_crm_rknd( iw_xy                , product(shape(iw_xy               )) )
    call memzero_crm_rknd( cld_xy               , product(shape(cld_xy              )) )
    call memzero_crm_rknd( u200_xy              , product(shape(u200_xy             )) )
    call memzero_crm_rknd( usfc_xy              , product(shape(usfc_xy             )) )
    call memzero_crm_rknd( v200_xy              , product(shape(v200_xy             )) )
    call memzero_crm_rknd( vsfc_xy              , product(shape(vsfc_xy             )) )
    call memzero_crm_rknd( w500_xy              , product(shape(w500_xy             )) )
    call memzero_crm_rknd( twle                 , product(shape(twle                )) )
    call memzero_crm_rknd( twsb                 , product(shape(twsb                )) )
    call memzero_crm_rknd( precflux             , product(shape(precflux            )) )
    call memzero_crm_rknd( uwle                 , product(shape(uwle                )) )
    call memzero_crm_rknd( uwsb                 , product(shape(uwsb                )) )
    call memzero_crm_rknd( vwle                 , product(shape(vwle                )) )
    call memzero_crm_rknd( vwsb                 , product(shape(vwsb                )) )
    call memzero_crm_rknd( tkelediss            , product(shape(tkelediss           )) )
    call memzero_crm_rknd( t2leadv              , product(shape(t2leadv             )) )
    call memzero_crm_rknd( t2legrad             , product(shape(t2legrad            )) )
    call memzero_crm_rknd( t2lediff             , product(shape(t2lediff            )) )
    call memzero_crm_rknd( t2lediss             , product(shape(t2lediss            )) )
    call memzero_crm_rknd( twleadv              , product(shape(twleadv             )) )
    call memzero_crm_rknd( twlediff             , product(shape(twlediff            )) )
    call memzero_crm_rknd( tadv                 , product(shape(tadv                )) )
    call memzero_crm_rknd( tdiff                , product(shape(tdiff               )) )
    call memzero_crm_rknd( tlat                 , product(shape(tlat                )) )
    call memzero_crm_rknd( tlatqi               , product(shape(tlatqi              )) )
    call memzero_crm_rknd( qifall               , product(shape(qifall              )) )
    call memzero_crm_rknd( qpfall               , product(shape(qpfall              )) )
    call memzero_crm_rknd( w_max                , product(shape(w_max               )) )
    call memzero_crm_rknd( u_max                , product(shape(u_max               )) )
    call memzero_real8   ( total_water_before   , product(shape(total_water_before  )) )
    call memzero_real8   ( total_water_after    , product(shape(total_water_after   )) )
    call memzero_real8   ( total_water_evap     , product(shape(total_water_evap    )) )
    call memzero_real8   ( total_water_prec     , product(shape(total_water_prec    )) )
    call memzero_real8   ( total_water_ls       , product(shape(total_water_ls      )) )
    call memzero_real8   ( total_water_clubb    , product(shape(total_water_clubb   )) )
    call memzero_real8   ( total_energy_before  , product(shape(total_energy_before )) )
    call memzero_real8   ( total_energy_after   , product(shape(total_energy_after  )) )
    call memzero_real8   ( total_energy_evap    , product(shape(total_energy_evap   )) )
    call memzero_real8   ( total_energy_prec    , product(shape(total_energy_prec   )) )
    call memzero_real8   ( total_energy_ls      , product(shape(total_energy_ls     )) )
    call memzero_real8   ( total_energy_clubb   , product(shape(total_energy_clubb  )) )
    call memzero_real8   ( total_energy_rad     , product(shape(total_energy_rad    )) )
    call memzero_crm_rknd( qtotmicro            , product(shape(qtotmicro           )) )
    call memzero_crm_rknd( CF3D                 , product(shape(CF3D                )) )
    call memzero_crm_rknd( u850_xy              , product(shape(u850_xy             )) )
    call memzero_crm_rknd( v850_xy              , product(shape(v850_xy             )) )
    call memzero_crm_rknd( psfc_xy              , product(shape(psfc_xy             )) )
    call memzero_crm_rknd( swvp_xy              , product(shape(swvp_xy             )) )
    call memzero_crm_rknd( cloudtopheight       , product(shape(cloudtopheight      )) )
    call memzero_crm_rknd( echotopheight        , product(shape(echotopheight       )) )
    call memzero_crm_rknd( cloudtoptemp         , product(shape(cloudtoptemp        )) )
#if (defined CRM && defined MODAL_AERO)
    call memzero_crm_rknd( naer  , product(shape(naer )) )
    call memzero_crm_rknd( vaer  , product(shape(vaer )) )
    call memzero_crm_rknd( hgaer , product(shape(hgaer)) )
#endif
  end subroutine allocate_vars


  subroutine deallocate_vars
    implicit none
    ! deallocate( u ) ! x-wind
    ! deallocate( v ) ! y-wind
    ! deallocate( w ) ! z-wind
    ! deallocate( t ) ! liquid/ice water static energy
    ! deallocate( p    )   ! perturbation pressure (from Poison eq)
    ! deallocate( tabs )   ! temperature
    ! deallocate( qv   )   ! water vapor
    ! deallocate( qcl  )   ! liquid water  (condensate)
    ! deallocate( qpl  )   ! liquid water  (precipitation)
    ! deallocate( qci  )   ! ice water  (condensate)
    ! deallocate( qpi  )   ! ice water  (precipitation)
    ! deallocate( tke2 )   ! SGS TKE
    ! deallocate( tk2  )   ! SGS eddyviscosity
    ! deallocate( dudt    )
    ! deallocate( dvdt    )
    ! deallocate( dwdt    )
    ! deallocate( misc    )
    ! deallocate( fluxbu  )
    ! deallocate( fluxbv  )
    ! deallocate( fluxbt  )
    ! deallocate( fluxbq  )
    ! deallocate( fluxtu  )
    ! deallocate( fluxtv  )
    ! deallocate( fluxtt  )
    ! deallocate( fluxtq  )
    ! deallocate( fzero   )
    ! deallocate( precsfc ) ! surface precip. rate
    ! deallocate( precssfc ) ! surface ice precip. rate
    ! deallocate( t0    )
    ! deallocate( q0    )
    ! deallocate( qv0   )
    ! deallocate( tabs0 )
    ! deallocate( tl0   )
    ! deallocate( tv0   )
    ! deallocate( u0    )
    ! deallocate( v0    )
    ! deallocate( tg0   )
    ! deallocate( qg0   )
    ! deallocate( ug0   )
    ! deallocate( vg0   )
    ! deallocate( p0    )
    ! deallocate( tke0  )
    ! deallocate( t01   )
    ! deallocate( q01   )
    ! deallocate( qp0   )
    ! deallocate( qn0   )
    ! deallocate( prespot   ) ! (1000./pres)**R/cp
    ! deallocate( rho       ) ! air density at pressure levels,kg/m3
    ! deallocate( rhow      ) ! air density at vertical velocity levels,kg/m3
    ! deallocate( bet       ) ! = ggr/tv0
    ! deallocate( gamaz     ) ! ggr/cp*z
    ! deallocate( wsub      ) ! Large-scale subsidence velocity,m/s
    ! deallocate( qtend     ) ! Large-scale tendency for total water
    ! deallocate( ttend     ) ! Large-scale tendency for temp.
    ! deallocate( utend     ) ! Large-scale tendency for u
    ! deallocate( vtend     ) ! Large-scale tendency for v
    ! deallocate( sstxy     ) !  surface temperature xy-distribution
    ! deallocate( fcory     )            !  Coriolis parameter xy-distribution
    ! deallocate( fcorzy    )            !  z-Coriolis parameter xy-distribution
    ! deallocate( latitude  )            ! latitude (icrm,degrees)
    ! deallocate( longitude )            ! longitude(icrm,degrees)
    ! deallocate( prec_xy   )            ! mean precip. rate for outout
    ! deallocate( pw_xy     )            ! precipitable water
    ! deallocate( cw_xy     )            ! cloud water path
    ! deallocate( iw_xy     )            ! ice water path
    ! deallocate( cld_xy    )            ! cloud frequency
    ! deallocate( u200_xy   )            ! u-wind at 200 mb
    ! deallocate( usfc_xy   )            ! u-wind at at the surface
    ! deallocate( v200_xy   )            ! v-wind at 200 mb
    ! deallocate( vsfc_xy   )            ! v-wind at the surface
    ! deallocate( w500_xy   )            ! w at 500 mb
    ! deallocate( twle      )
    ! deallocate( twsb      )
    ! deallocate( precflux  )
    ! deallocate( uwle      )
    ! deallocate( uwsb      )
    ! deallocate( vwle      )
    ! deallocate( vwsb      )
    ! deallocate( tkelediss )
    ! deallocate( t2leadv   )
    ! deallocate( t2legrad  )
    ! deallocate( t2lediff  )
    ! deallocate( t2lediss  )
    ! deallocate( twleadv   )
    ! deallocate( twlediff  )
    ! deallocate( tadv      )
    ! deallocate( tdiff     )
    ! deallocate( tlat      )
    ! deallocate( tlatqi    )
    ! deallocate( qifall    )
    ! deallocate( qpfall    )
    ! deallocate( w_max     )
    ! deallocate( u_max     )
    ! deallocate( total_water_before  )
    ! deallocate( total_water_after   )
    ! deallocate( total_water_evap    )
    ! deallocate( total_water_prec    )
    ! deallocate( total_water_ls      )
    ! deallocate( total_water_clubb   )
    ! deallocate( total_energy_before )
    ! deallocate( total_energy_after  )
    ! deallocate( total_energy_evap   )
    ! deallocate( total_energy_prec   )
    ! deallocate( total_energy_ls     )
    ! deallocate( total_energy_clubb  )
    ! deallocate( total_energy_rad    )
    ! deallocate( qtotmicro )
    ! deallocate( CF3D           )
    ! deallocate( u850_xy        )
    ! deallocate( v850_xy        )
    ! deallocate( psfc_xy        )
    ! deallocate( swvp_xy        )
    ! deallocate( cloudtopheight )
    ! deallocate( echotopheight  )
    ! deallocate( cloudtoptemp   )
    call pool_pop_multiple(115)
#if (defined CRM && defined MODAL_AERO)
    ! deallocate( naer           )
    ! deallocate( vaer           )
    ! deallocate( hgaer          )
    call pool_pop_multiple(3)
#endif
  end subroutine deallocate_vars

end module vars
