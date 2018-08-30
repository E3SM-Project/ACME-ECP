module vars
  use grid
  use params, only: crm_rknd
#ifdef CRM
#ifdef MODAL_AERO
  use modal_aero_data,   only: ntot_amode
#endif
#endif

  implicit none
  !--------------------------------------------------------------------
  ! prognostic variables:

  real(crm_rknd), allocatable :: u   (:,:,:) ! x-wind
  real(crm_rknd), allocatable :: v   (:,:,:) ! y-wind
  real(crm_rknd), allocatable :: w   (:,:,:) ! z-wind
  real(crm_rknd), allocatable :: t   (:,:,:) ! liquid/ice water static energy

  !--------------------------------------------------------------------
  ! diagnostic variables:

  real(crm_rknd), allocatable :: p       (:,:,:,:)     ! perturbation pressure (from Poison eq)
  real(crm_rknd), allocatable :: tabs    (:,:,:,:)                 ! temperature
  real(crm_rknd), allocatable :: qv      (:,:,:,:)                ! water vapor
  real(crm_rknd), allocatable :: qcl     (:,:,:,:)                ! liquid water  (condensate)
  real(crm_rknd), allocatable :: qpl     (:,:,:,:)                ! liquid water  (precipitation)
  real(crm_rknd), allocatable :: qci     (:,:,:,:)                ! ice water  (condensate)
  real(crm_rknd), allocatable :: qpi     (:,:,:,:)                ! ice water  (precipitation)
  real(crm_rknd), allocatable :: tke2    (:,:,:,:)   ! SGS TKE
  real(crm_rknd), allocatable :: tk2     (:,:,:,:) ! SGS eddyviscosity

  !--------------------------------------------------------------------
  ! time-tendencies for prognostic variables

  real(crm_rknd), allocatable :: dudt   (:,:,:,:,:)
  real(crm_rknd), allocatable :: dvdt   (:,:,:,:,:)
  real(crm_rknd), allocatable :: dwdt   (:,:,:,:,:)

  !----------------------------------------------------------------
  ! Temporary storage array:

  real(crm_rknd), allocatable :: misc(:,:,:)
  !------------------------------------------------------------------
  ! fluxes at the top and bottom of the domain:

  real(crm_rknd), allocatable :: fluxbu  (:,:)
  real(crm_rknd), allocatable :: fluxbv  (:,:)
  real(crm_rknd), allocatable :: fluxbt  (:,:)
  real(crm_rknd), allocatable :: fluxbq  (:,:)
  real(crm_rknd), allocatable :: fluxtu  (:,:)
  real(crm_rknd), allocatable :: fluxtv  (:,:)
  real(crm_rknd), allocatable :: fluxtt  (:,:)
  real(crm_rknd), allocatable :: fluxtq  (:,:)
  real(crm_rknd), allocatable :: fzero   (:,:)
  real(crm_rknd), allocatable :: precsfc (:,:) ! surface precip. rate
  real(crm_rknd), allocatable :: precssfc(:,:) ! surface ice precip. rate

  !-----------------------------------------------------------------
  ! profiles

  real(crm_rknd), allocatable :: t0   (:)
  real(crm_rknd), allocatable :: q0   (:)
  real(crm_rknd), allocatable :: qv0  (:)
  real(crm_rknd), allocatable :: tabs0(:)
  real(crm_rknd), allocatable :: tl0  (:)
  real(crm_rknd), allocatable :: tv0  (:)
  real(crm_rknd), allocatable :: u0   (:)
  real(crm_rknd), allocatable :: v0   (:)
  real(crm_rknd), allocatable :: tg0  (:)
  real(crm_rknd), allocatable :: qg0  (:)
  real(crm_rknd), allocatable :: ug0  (:)
  real(crm_rknd), allocatable :: vg0  (:)
  real(crm_rknd), allocatable :: p0   (:)
  real(crm_rknd), allocatable :: tke0 (:)
  real(crm_rknd), allocatable :: t01  (:)
  real(crm_rknd), allocatable :: q01  (:)
  real(crm_rknd), allocatable :: qp0  (:)
  real(crm_rknd), allocatable :: qn0  (:)
  !----------------------------------------------------------------
  ! "observed" (read from snd file) surface characteristics

  real(crm_rknd)  sstobs, lhobs, shobs
  !----------------------------------------------------------------
  !  Domain top stuff:

  real(crm_rknd)   gamt0    ! gradient of t() at the top,K/m
  real(crm_rknd)   gamq0    ! gradient of q() at the top,g/g/m

  !-----------------------------------------------------------------
  ! reference vertical profiles:

  real(crm_rknd), allocatable :: prespot(:)  ! (1000./pres)**R/cp
  real(crm_rknd), allocatable :: rho    (:)   ! air density at pressure levels,kg/m3
  real(crm_rknd), allocatable :: rhow   (:)   ! air density at vertical velocity levels,kg/m3
  real(crm_rknd), allocatable :: bet    (:)   ! = ggr/tv0
  real(crm_rknd), allocatable :: gamaz  (:) ! ggr/cp*z
  real(crm_rknd), allocatable :: wsub   (:)   ! Large-scale subsidence velocity,m/s
  real(crm_rknd), allocatable :: qtend  (:) ! Large-scale tendency for total water
  real(crm_rknd), allocatable :: ttend  (:) ! Large-scale tendency for temp.
  real(crm_rknd), allocatable :: utend  (:) ! Large-scale tendency for u
  real(crm_rknd), allocatable :: vtend  (:) ! Large-scale tendency for v

  !---------------------------------------------------------------------
  ! Large-scale and surface forcing:

  integer nlsf  ! number of large-scale forcing profiles
  integer nrfc  ! number of radiative forcing profiles
  integer nsfc  ! number of surface forcing profiles
  integer nsnd  ! number of observed soundings
  integer nzlsf ! number of large-scale forcing profiles
  integer nzrfc ! number of radiative forcing profiles
  integer nzsnd ! number of observed soundings

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! MRN: Already previously allocated. I'm leaving these alone
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(crm_rknd), allocatable :: dqls(:,:) ! Large-scale tendency for total water
  real(crm_rknd), allocatable :: dtls(:,:) ! Large-scale tendency for temp.
  real(crm_rknd), allocatable :: ugls(:,:) ! Large-scale wind in X-direction
  real(crm_rknd), allocatable :: vgls(:,:) ! Large-scale wind in Y-direction
  real(crm_rknd), allocatable :: wgls(:,:) ! Large-scale subsidence velocity,m/s
  real(crm_rknd), allocatable :: pres0ls(:)! Surface pressure, mb
  real(crm_rknd), allocatable :: zls(:,:)  ! Height
  real(crm_rknd), allocatable :: pls(:,:)  ! Pressure
  real(crm_rknd), allocatable :: dayls(:)  ! Large-scale forcing arrays time (days)
  real(crm_rknd), allocatable :: dtrfc(:,:)! Radiative tendency for pot. temp.
  real(crm_rknd), allocatable :: dayrfc(:) ! Radiative forcing arrays time (days)
  real(crm_rknd), allocatable :: prfc(:,:) ! Pressure/Height
  real(crm_rknd), allocatable :: sstsfc(:) ! SSTs
  real(crm_rknd), allocatable :: shsfc(:)   ! Sensible heat flux,W/m2
  real(crm_rknd), allocatable :: lhsfc(:)  ! Latent heat flux,W/m2
  real(crm_rknd), allocatable :: tausfc(:) ! Surface drag,m2/s2
  real(crm_rknd), allocatable :: daysfc(:) ! Surface forcing arrays time (days)
  real(crm_rknd), allocatable :: usnd(:,:) ! Observed zonal wind
  real(crm_rknd), allocatable :: vsnd(:,:) ! Observed meriod wind
  real(crm_rknd), allocatable :: tsnd(:,:) ! Observed Abs. temperature
  real(crm_rknd), allocatable :: qsnd(:,:) ! Observed Moisture
  real(crm_rknd), allocatable :: zsnd(:,:) ! Height
  real(crm_rknd), allocatable :: psnd(:,:) ! Pressure
  real(crm_rknd), allocatable :: daysnd(:) ! number of sounding samples
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !---------------------------------------------------------------------
  !  Horizontally varying stuff (as a function of xy)
  !
  real(crm_rknd), allocatable :: sstxy    (:,:) !  surface temperature xy-distribution
  real(crm_rknd), allocatable :: fcory    (:)      !  Coriolis parameter xy-distribution
  real(crm_rknd), allocatable :: fcorzy   (:)      !  z-Coriolis parameter xy-distribution
  real(crm_rknd), allocatable :: latitude (:,:)      ! latitude (degrees)
  real(crm_rknd), allocatable :: longitude(:,:)      ! longitude(degrees)
  real(crm_rknd), allocatable :: prec_xy  (:,:) ! mean precip. rate for outout
  real(crm_rknd), allocatable :: shf_xy   (:,:) ! mean precip. rate for outout
  real(crm_rknd), allocatable :: lhf_xy   (:,:) ! mean precip. rate for outout
  real(crm_rknd), allocatable :: lwns_xy  (:,:) ! mean net lw at SFC
  real(crm_rknd), allocatable :: swns_xy  (:,:) ! mean net sw at SFC
  real(crm_rknd), allocatable :: lwnsc_xy (:,:) ! clear-sky mean net lw at SFC
  real(crm_rknd), allocatable :: swnsc_xy (:,:) ! clear-sky mean net sw at SFC
  real(crm_rknd), allocatable :: lwnt_xy  (:,:) ! mean net lw at TOA
  real(crm_rknd), allocatable :: swnt_xy  (:,:) ! mean net sw at TOA
  real(crm_rknd), allocatable :: lwntc_xy (:,:) ! clear-sky mean net lw at TOA
  real(crm_rknd), allocatable :: swntc_xy (:,:) ! clear-sky mean net sw at TOA
  real(crm_rknd), allocatable :: solin_xy (:,:) ! solar TOA insolation
  real(crm_rknd), allocatable :: pw_xy    (:,:)   ! precipitable water
  real(crm_rknd), allocatable :: cw_xy    (:,:)   ! cloud water path
  real(crm_rknd), allocatable :: iw_xy    (:,:)   ! ice water path
  real(crm_rknd), allocatable :: cld_xy   (:,:)   ! cloud frequency
  real(crm_rknd), allocatable :: u200_xy  (:,:) ! u-wind at 200 mb
  real(crm_rknd), allocatable :: usfc_xy  (:,:) ! u-wind at at the surface
  real(crm_rknd), allocatable :: v200_xy  (:,:) ! v-wind at 200 mb
  real(crm_rknd), allocatable :: vsfc_xy  (:,:) ! v-wind at the surface
  real(crm_rknd), allocatable :: w500_xy  (:,:) ! w at 500 mb
  real(crm_rknd), allocatable :: qocean_xy(:,:) ! ocean cooling in W/m2

  !----------------------------------------------------------------------
  ! Vertical profiles of quantities sampled for statitistics purposes:

  real(crm_rknd) :: w_max
  real(crm_rknd) :: u_max
  real(crm_rknd) :: s_acld
  real(crm_rknd) :: s_acldcold
  real(crm_rknd) :: s_ar
  real(crm_rknd) :: s_arthr
  real(crm_rknd) :: s_sst
  real(crm_rknd) :: s_acldl
  real(crm_rknd) :: s_acldm
  real(crm_rknd) :: s_acldh
  real(crm_rknd) :: ncmn
  real(crm_rknd) :: nrmn
  real(crm_rknd) :: z_inv
  real(crm_rknd) :: z_cb
  real(crm_rknd) :: z_ct
  real(crm_rknd) :: z_cbmn
  real(crm_rknd) :: z_ctmn
  real(crm_rknd) :: z2_inv
  real(crm_rknd) :: z2_cb
  real(crm_rknd) :: z2_ct
  real(crm_rknd) :: cwpmean
  real(crm_rknd) :: cwp2
  real(crm_rknd) :: precmean
  real(crm_rknd) :: prec2
  real(crm_rknd) :: precmax
  real(crm_rknd) :: nrainy
  real(crm_rknd) :: ncloudy
  real(crm_rknd) :: s_acldisccp
  real(crm_rknd) :: s_acldlisccp
  real(crm_rknd) :: s_acldmisccp
  real(crm_rknd) :: s_acldhisccp
  real(crm_rknd) :: s_ptopisccp
  real(crm_rknd) :: s_acldmodis
  real(crm_rknd) :: s_acldlmodis
  real(crm_rknd) :: s_acldmmodis
  real(crm_rknd) :: s_acldhmodis
  real(crm_rknd) :: s_ptopmodis
  real(crm_rknd) :: s_acldmisr
  real(crm_rknd) :: s_ztopmisr
  real(crm_rknd) :: s_relmodis
  real(crm_rknd) :: s_reimodis
  real(crm_rknd) :: s_lwpmodis
  real(crm_rknd) :: s_iwpmodis
  real(crm_rknd) :: s_tbisccp
  real(crm_rknd) :: s_tbclrisccp
  real(crm_rknd) :: s_acldliqmodis
  real(crm_rknd) :: s_acldicemodis
  real(crm_rknd) :: s_cldtauisccp
  real(crm_rknd) :: s_cldtaumodis
  real(crm_rknd) :: s_cldtaulmodis
  real(crm_rknd) :: s_cldtauimodis
  real(crm_rknd) :: s_cldalbisccp
  real(crm_rknd) :: s_flns
  real(crm_rknd) :: s_flnt
  real(crm_rknd) :: s_flntoa
  real(crm_rknd) :: s_flnsc
  real(crm_rknd) :: s_flntoac
  real(crm_rknd) :: s_flds
  real(crm_rknd) :: s_fsns
  real(crm_rknd) :: s_fsnt
  real(crm_rknd) :: s_fsntoa
  real(crm_rknd) :: s_fsnsc
  real(crm_rknd) :: s_fsntoac
  real(crm_rknd) :: s_fsds
  real(crm_rknd) :: s_solin
  real(crm_rknd), allocatable :: twle(:)
  real(crm_rknd), allocatable :: twsb(:)
  real(crm_rknd), allocatable :: precflux(:)
  real(crm_rknd), allocatable :: uwle(:)
  real(crm_rknd), allocatable :: uwsb(:)
  real(crm_rknd), allocatable :: vwle(:)
  real(crm_rknd), allocatable :: vwsb(:)
  real(crm_rknd), allocatable :: radlwup(:)
  real(crm_rknd), allocatable :: radlwdn(:)
  real(crm_rknd), allocatable :: radswup(:)
  real(crm_rknd), allocatable :: radswdn(:)
  real(crm_rknd), allocatable :: radqrlw(:)
  real(crm_rknd), allocatable :: radqrsw(:)
  real(crm_rknd), allocatable :: tkeleadv(:)
  real(crm_rknd), allocatable :: tkelepress(:)
  real(crm_rknd), allocatable :: tkelediss(:)
  real(crm_rknd), allocatable :: tkelediff(:)
  real(crm_rknd), allocatable :: tkelebuoy(:)
  real(crm_rknd), allocatable :: t2leadv(:)
  real(crm_rknd), allocatable :: t2legrad(:)
  real(crm_rknd), allocatable :: t2lediff(:)
  real(crm_rknd), allocatable :: t2leprec(:)
  real(crm_rknd), allocatable :: t2lediss(:)
  real(crm_rknd), allocatable :: q2leadv(:)
  real(crm_rknd), allocatable :: q2legrad(:)
  real(crm_rknd), allocatable :: q2lediff(:)
  real(crm_rknd), allocatable :: q2leprec(:)
  real(crm_rknd), allocatable :: q2lediss(:)
  real(crm_rknd), allocatable :: twleadv(:)
  real(crm_rknd), allocatable :: twlediff(:)
  real(crm_rknd), allocatable :: twlepres(:)
  real(crm_rknd), allocatable :: twlebuoy(:)
  real(crm_rknd), allocatable :: twleprec(:)
  real(crm_rknd), allocatable :: qwleadv(:)
  real(crm_rknd), allocatable :: qwlediff(:)
  real(crm_rknd), allocatable :: qwlepres(:)
  real(crm_rknd), allocatable :: qwlebuoy(:)
  real(crm_rknd), allocatable :: qwleprec(:)
  real(crm_rknd), allocatable :: momleadv(:,:)
  real(crm_rknd), allocatable :: momlepress(:,:)
  real(crm_rknd), allocatable :: momlebuoy(:,:)
  real(crm_rknd), allocatable :: momlediff(:,:)
  real(crm_rknd), allocatable :: tadv(:)
  real(crm_rknd), allocatable :: tdiff(:)
  real(crm_rknd), allocatable :: tlat(:)
  real(crm_rknd), allocatable :: tlatqi(:)
  real(crm_rknd), allocatable :: qifall(:)
  real(crm_rknd), allocatable :: qpfall(:)
  real(crm_rknd), allocatable :: tdiff_xy(:)
  real(crm_rknd), allocatable :: tdiff_z(:)
  real(crm_rknd), allocatable :: ttest0(:)
  real(crm_rknd), allocatable :: ttest1(:)
  real(crm_rknd), allocatable :: ttest2(:,:)  !+++mhwang test

  ! energy conservation diagnostics:

  real(8) total_water_before, total_water_after
  real(8) total_water_evap, total_water_prec, total_water_ls
  !#ifdef CLUBB_CRM
  real(8) total_water_clubb
  real(8) total_energy_before, total_energy_after
  real(8) total_energy_evap, total_energy_prec, total_energy_ls
  real(8) total_energy_clubb, total_energy_rad
  !#endif
  real(8), allocatable :: qtotmicro(:)  ! total water for water conservation test in microphysics +++mhwang

  !===========================================================================
  ! UW ADDITIONS

  ! conditional average statistics, subsumes cloud_factor, core_factor, coredn_factor
  integer :: ncondavg, icondavg_cld, icondavg_cor, icondavg_cordn, &
  icondavg_satdn, icondavg_satup, icondavg_env

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! MRN: Already allocated, I'm leaving them alone
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(crm_rknd), allocatable :: condavg_factor(:,:) ! replaces cloud_factor, core_factor
  real(crm_rknd), allocatable :: condavg_mask(:,:,:,:) ! indicator array for various conditional averages
  character(LEN=8), allocatable :: condavgname(:) ! array of short names
  character(LEN=25), allocatable :: condavglongname(:) ! array of long names
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(crm_rknd), allocatable :: qlsvadv(:) ! Large-scale vertical advection tendency for total water
  real(crm_rknd), allocatable :: tlsvadv(:) ! Large-scale vertical advection tendency for temperature
  real(crm_rknd), allocatable :: ulsvadv(:) ! Large-scale vertical advection tendency for zonal velocity
  real(crm_rknd), allocatable :: vlsvadv(:) ! Large-scale vertical advection tendency for meridional velocity
  real(crm_rknd), allocatable :: qnudge(:) ! Nudging of horiz.-averaged total water profile
  real(crm_rknd), allocatable :: tnudge(:) ! Nudging of horiz.-averaged temperature profile
  real(crm_rknd), allocatable :: unudge(:) ! Nudging of horiz.-averaged zonal velocity
  real(crm_rknd), allocatable :: vnudge(:) ! Nudging of horiz.-averaged meridional velocity
  real(crm_rknd), allocatable :: qstor(:) ! Storage of horiz.-averaged total water profile
  real(crm_rknd), allocatable :: tstor(:) ! Storage of horiz.-averaged temperature profile
  real(crm_rknd), allocatable :: ustor(:) ! Storage of horiz.-averaged zonal velocity
  real(crm_rknd), allocatable :: vstor(:) ! Storage of horiz.-averaged meridional velocity
  real(crm_rknd), allocatable :: qtostor(:) ! Storage of horiz.-averaged total water profile (vapor + liquid)
  real(crm_rknd), allocatable :: utendcor(:) ! coriolis acceleration of zonal velocity
  real(crm_rknd), allocatable :: vtendcor(:) ! coriolis acceleration of meridional velocity
  real(crm_rknd), allocatable :: CF3D(:,:,:)  ! Cloud fraction
  ! =1.0 when there is no fractional cloudiness scheme
  ! = cloud fraction produced by fractioal cloudiness scheme when avaiable

  ! 850 mbar horizontal winds
  real(crm_rknd), allocatable :: u850_xy(:,:) ! zonal velocity at 850 mb
  real(crm_rknd), allocatable :: v850_xy(:,:) ! meridional velocity at 850 mb

  ! Surface pressure
  real(crm_rknd), allocatable :: psfc_xy(:,:) ! pressure (in millibar) at lowest grid point

  ! Saturated water vapor path, useful for computing column relative humidity
  real(crm_rknd), allocatable :: swvp_xy(:,:)  ! saturated water vapor path (wrt water)

  ! Cloud and echo top heights, and cloud top temperature (instantaneous)
  real(crm_rknd), allocatable :: cloudtopheight(:,:)
  real(crm_rknd), allocatable :: echotopheight (:,:)
  real(crm_rknd), allocatable :: cloudtoptemp  (:,:)

  ! END UW ADDITIONS
  !===========================================================================
#if (defined CRM && defined MODAL_AERO)
  real(crm_rknd), allocatable :: naer (:,:)     ! Aerosol number concentration [/m3]
  real(crm_rknd), allocatable :: vaer (:,:)     ! aerosol volume concentration [m3/m3]
  real(crm_rknd), allocatable :: hgaer(:,:)    ! hygroscopicity of aerosol mode
#endif


contains


  subroutine allocate_vars(ncrms)
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) :: zero
    allocate( u   (dimx1_u:dimx2_u,dimy1_u:dimy2_u,nzm)  )
    allocate( v   (dimx1_v:dimx2_v,dimy1_v:dimy2_v,nzm)  )
    allocate( w   (dimx1_w:dimx2_w,dimy1_w:dimy2_w,nz )  )
    allocate( t   (dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm)  )
    allocate( p       (0:nx, (1-YES3D):ny, nzm, ncrms)      )
    allocate( tabs    (nx, ny, nzm, ncrms)                  )
    allocate( qv      (nx, ny, nzm, ncrms)                 )
    allocate( qcl     (nx, ny, nzm, ncrms)                 )
    allocate( qpl     (nx, ny, nzm, ncrms)                 )
    allocate( qci     (nx, ny, nzm, ncrms)                 )
    allocate( qpi     (nx, ny, nzm, ncrms)                 )
    allocate( tke2(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, ncrms)    )
    allocate( tk2  (0:nxp1, (1-YES3D):nyp1, nzm, ncrms)  )
    allocate( dudt   (nxp1, ny, nzm, 3, ncrms) )
    allocate( dvdt   (nx, nyp1, nzm, 3, ncrms) )
    allocate( dwdt   (nx, ny  , nz,  3, ncrms) )
    allocate( misc(nx, ny, nz) )
    allocate( fluxbu  (nx,ny) )
    allocate( fluxbv  (nx,ny) )
    allocate( fluxbt  (nx,ny) )
    allocate( fluxbq  (nx,ny) )
    allocate( fluxtu  (nx,ny) )
    allocate( fluxtv  (nx,ny) )
    allocate( fluxtt  (nx,ny) )
    allocate( fluxtq  (nx,ny) )
    allocate( fzero   (nx,ny) )
    allocate( precsfc (nx,ny)  )
    allocate( precssfc(nx,ny)  )
    allocate( t0   (nzm) )
    allocate( q0   (nzm) )
    allocate( qv0  (nzm) )
    allocate( tabs0(nzm) )
    allocate( tl0  (nzm) )
    allocate( tv0  (nzm) )
    allocate( u0   (nzm) )
    allocate( v0   (nzm) )
    allocate( tg0  (nzm) )
    allocate( qg0  (nzm) )
    allocate( ug0  (nzm) )
    allocate( vg0  (nzm) )
    allocate( p0   (nzm) )
    allocate( tke0 (nzm) )
    allocate( t01  (nzm) )
    allocate( q01  (nzm) )
    allocate( qp0  (nzm) )
    allocate( qn0  (nzm) )
    allocate( prespot(nzm)   )
    allocate( rho    (nzm)     )
    allocate( rhow   (nz )    )
    allocate( bet    (nzm)     )
    allocate( gamaz  (nzm)  )
    allocate( wsub   (nz )    )
    allocate( qtend  (nzm)  )
    allocate( ttend  (nzm)  )
    allocate( utend  (nzm)  )
    allocate( vtend  (nzm)  )
    allocate( sstxy    (0:nx,(1-YES3D):ny)   )
    allocate( fcory    (0:ny)       )
    allocate( fcorzy   (ny)       )
    allocate( latitude (nx,ny)        )
    allocate( longitude(nx,ny)        )
    allocate( prec_xy  (nx,ny)  )
    allocate( shf_xy   (nx,ny)  )
    allocate( lhf_xy   (nx,ny)  )
    allocate( lwns_xy  (nx,ny)  )
    allocate( swns_xy  (nx,ny)  )
    allocate( lwnsc_xy (nx,ny)  )
    allocate( swnsc_xy (nx,ny)  )
    allocate( lwnt_xy  (nx,ny)  )
    allocate( swnt_xy  (nx,ny)  )
    allocate( lwntc_xy (nx,ny)  )
    allocate( swntc_xy (nx,ny)  )
    allocate( solin_xy (nx,ny)  )
    allocate( pw_xy    (nx,ny)    )
    allocate( cw_xy    (nx,ny)    )
    allocate( iw_xy    (nx,ny)    )
    allocate( cld_xy   (nx,ny)    )
    allocate( u200_xy  (nx,ny)  )
    allocate( usfc_xy  (nx,ny)  )
    allocate( v200_xy  (nx,ny)  )
    allocate( vsfc_xy  (nx,ny)  )
    allocate( w500_xy  (nx,ny)  )
    allocate( qocean_xy(nx,ny)  )
    allocate( twle(nz) )
    allocate( twsb(nz) )
    allocate( precflux(nz) )
    allocate( uwle(nz) )
    allocate( uwsb(nz) )
    allocate( vwle(nz) )
    allocate( vwsb(nz) )
    allocate( radlwup(nz) )
    allocate( radlwdn(nz) )
    allocate( radswup(nz) )
    allocate( radswdn(nz) )
    allocate( radqrlw(nz) )
    allocate( radqrsw(nz) )
    allocate( tkeleadv(nz) )
    allocate( tkelepress(nz) )
    allocate( tkelediss(nz) )
    allocate( tkelediff(nz) )
    allocate( tkelebuoy(nz) )
    allocate( t2leadv(nz) )
    allocate( t2legrad(nz) )
    allocate( t2lediff(nz) )
    allocate( t2leprec(nz) )
    allocate( t2lediss(nz) )
    allocate( q2leadv(nz) )
    allocate( q2legrad(nz) )
    allocate( q2lediff(nz) )
    allocate( q2leprec(nz) )
    allocate( q2lediss(nz) )
    allocate( twleadv(nz) )
    allocate( twlediff(nz) )
    allocate( twlepres(nz) )
    allocate( twlebuoy(nz) )
    allocate( twleprec(nz) )
    allocate( qwleadv(nz) )
    allocate( qwlediff(nz) )
    allocate( qwlepres(nz) )
    allocate( qwlebuoy(nz) )
    allocate( qwleprec(nz) )
    allocate( momleadv(nz,3) )
    allocate( momlepress(nz,3) )
    allocate( momlebuoy(nz,3) )
    allocate( momlediff(nz,3) )
    allocate( tadv(nz) )
    allocate( tdiff(nz) )
    allocate( tlat(nz) )
    allocate( tlatqi(nz) )
    allocate( qifall(nz) )
    allocate( qpfall(nz) )
    allocate( tdiff_xy(nz) )
    allocate( tdiff_z(nz) )
    allocate( ttest0(nzm) )
    allocate( ttest1(nz) )
    allocate( ttest2(nz,10)   )
    allocate( qtotmicro(5)   )
    allocate( qlsvadv(nzm)  )
    allocate( tlsvadv(nzm)  )
    allocate( ulsvadv(nzm)  )
    allocate( vlsvadv(nzm)  )
    allocate( qnudge(nzm)  )
    allocate( tnudge(nzm)  )
    allocate( unudge(nzm)  )
    allocate( vnudge(nzm)  )
    allocate( qstor(nzm)  )
    allocate( tstor(nzm)  )
    allocate( ustor(nzm)  )
    allocate( vstor(nzm)  )
    allocate( qtostor(nzm)  )
    allocate( utendcor(nzm)  )
    allocate( vtendcor(nzm)  )
    allocate( CF3D(1:nx, 1:ny, 1:nzm)   )
    allocate( u850_xy(nx,ny)  )
    allocate( v850_xy(nx,ny)  )
    allocate( psfc_xy(nx,ny)  )
    allocate( swvp_xy(nx,ny)   )
    allocate( cloudtopheight(nx,ny) )
    allocate( echotopheight (nx,ny) )
    allocate( cloudtoptemp  (nx,ny) )
#if (defined CRM && defined MODAL_AERO)
    allocate( naer (nzm, ntot_amode) )
    allocate( vaer (nzm, ntot_amode) )
    allocate( hgaer(nzm, ntot_amode) )
#endif

    zero = 0

    u = zero
    v = zero
    w = zero
    t = zero
    p = zero
    tabs = zero
    qv = zero
    qcl = zero
    qpl = zero
    qci = zero
    qpi = zero
    tke2 = zero
    tk2 = zero
    dudt = zero
    dvdt = zero
    dwdt = zero
    misc = zero
    fluxbu = zero
    fluxbv = zero
    fluxbt = zero
    fluxbq = zero
    fluxtu = zero
    fluxtv = zero
    fluxtt = zero
    fluxtq = zero
    fzero = zero
    precsfc = zero
    precssfc = zero
    t0 = zero
    q0 = zero
    qv0 = zero
    tabs0 = zero
    tl0 = zero
    tv0 = zero
    u0 = zero
    v0 = zero
    tg0 = zero
    qg0 = zero
    ug0 = zero
    vg0 = zero
    p0 = zero
    tke0 = zero
    t01 = zero
    q01 = zero
    qp0 = zero
    qn0 = zero
    prespot = zero
    rho = zero
    rhow = zero
    bet = zero
    gamaz = zero
    wsub = zero
    qtend = zero
    ttend = zero
    utend = zero
    vtend = zero
    sstxy = zero
    fcory = zero
    fcorzy = zero
    latitude = zero
    longitude = zero
    prec_xy = zero
    shf_xy = zero
    lhf_xy = zero
    lwns_xy = zero
    swns_xy = zero
    lwnsc_xy = zero
    swnsc_xy = zero
    lwnt_xy = zero
    swnt_xy = zero
    lwntc_xy = zero
    swntc_xy = zero
    solin_xy = zero
    pw_xy = zero
    cw_xy = zero
    iw_xy = zero
    cld_xy = zero
    u200_xy = zero
    usfc_xy = zero
    v200_xy = zero
    vsfc_xy = zero
    w500_xy = zero
    qocean_xy = zero
    twle = zero
    twsb = zero
    precflux = zero
    uwle = zero
    uwsb = zero
    vwle = zero
    vwsb = zero
    radlwup = zero
    radlwdn = zero
    radswup = zero
    radswdn = zero
    radqrlw = zero
    radqrsw = zero
    tkeleadv = zero
    tkelepress = zero
    tkelediss = zero
    tkelediff = zero
    tkelebuoy = zero
    t2leadv = zero
    t2legrad = zero
    t2lediff = zero
    t2leprec = zero
    t2lediss = zero
    q2leadv = zero
    q2legrad = zero
    q2lediff = zero
    q2leprec = zero
    q2lediss = zero
    twleadv = zero
    twlediff = zero
    twlepres = zero
    twlebuoy = zero
    twleprec = zero
    qwleadv = zero
    qwlediff = zero
    qwlepres = zero
    qwlebuoy = zero
    qwleprec = zero
    momleadv = zero
    momlepress = zero
    momlebuoy = zero
    momlediff = zero
    tadv = zero
    tdiff = zero
    tlat = zero
    tlatqi = zero
    qifall = zero
    qpfall = zero
    tdiff_xy = zero
    tdiff_z = zero
    ttest0 = zero
    ttest1 = zero
    ttest2 = zero
    qtotmicro = zero
    qlsvadv = zero
    tlsvadv = zero
    ulsvadv = zero
    vlsvadv = zero
    qnudge = zero
    tnudge = zero
    unudge = zero
    vnudge = zero
    qstor = zero
    tstor = zero
    ustor = zero
    vstor = zero
    qtostor = zero
    utendcor = zero
    vtendcor = zero
    CF3D = zero
    u850_xy = zero
    v850_xy = zero
    psfc_xy = zero
    swvp_xy = zero
    cloudtopheight = zero
    echotopheight = zero
    cloudtoptemp = zero
#if (defined CRM && defined MODAL_AERO)
    naer = zero
    vaer = zero
    hgaer = zero
#endif
  end subroutine allocate_vars


  subroutine deallocate_vars()
    implicit none
    deallocate( u )
    deallocate( v )
    deallocate( w )
    deallocate( t )
    deallocate( p )
    deallocate( tabs )
    deallocate( qv )
    deallocate( qcl )
    deallocate( qpl )
    deallocate( qci )
    deallocate( qpi )
    deallocate( tke2 )
    deallocate( tk2 )
    deallocate( dudt )
    deallocate( dvdt )
    deallocate( dwdt )
    deallocate( misc )
    deallocate( fluxbu )
    deallocate( fluxbv )
    deallocate( fluxbt )
    deallocate( fluxbq )
    deallocate( fluxtu )
    deallocate( fluxtv )
    deallocate( fluxtt )
    deallocate( fluxtq )
    deallocate( fzero )
    deallocate( precsfc )
    deallocate( precssfc )
    deallocate( t0 )
    deallocate( q0 )
    deallocate( qv0 )
    deallocate( tabs0 )
    deallocate( tl0 )
    deallocate( tv0 )
    deallocate( u0 )
    deallocate( v0 )
    deallocate( tg0 )
    deallocate( qg0 )
    deallocate( ug0 )
    deallocate( vg0 )
    deallocate( p0 )
    deallocate( tke0 )
    deallocate( t01 )
    deallocate( q01 )
    deallocate( qp0 )
    deallocate( qn0 )
    deallocate( prespot )
    deallocate( rho )
    deallocate( rhow )
    deallocate( bet )
    deallocate( gamaz )
    deallocate( wsub )
    deallocate( qtend )
    deallocate( ttend )
    deallocate( utend )
    deallocate( vtend )
    deallocate( sstxy )
    deallocate( fcory )
    deallocate( fcorzy )
    deallocate( latitude )
    deallocate( longitude )
    deallocate( prec_xy )
    deallocate( shf_xy )
    deallocate( lhf_xy )
    deallocate( lwns_xy )
    deallocate( swns_xy )
    deallocate( lwnsc_xy )
    deallocate( swnsc_xy )
    deallocate( lwnt_xy )
    deallocate( swnt_xy )
    deallocate( lwntc_xy )
    deallocate( swntc_xy )
    deallocate( solin_xy )
    deallocate( pw_xy )
    deallocate( cw_xy )
    deallocate( iw_xy )
    deallocate( cld_xy )
    deallocate( u200_xy )
    deallocate( usfc_xy )
    deallocate( v200_xy )
    deallocate( vsfc_xy )
    deallocate( w500_xy )
    deallocate( qocean_xy )
    deallocate( twle )
    deallocate( twsb )
    deallocate( precflux )
    deallocate( uwle )
    deallocate( uwsb )
    deallocate( vwle )
    deallocate( vwsb )
    deallocate( radlwup )
    deallocate( radlwdn )
    deallocate( radswup )
    deallocate( radswdn )
    deallocate( radqrlw )
    deallocate( radqrsw )
    deallocate( tkeleadv )
    deallocate( tkelepress )
    deallocate( tkelediss )
    deallocate( tkelediff )
    deallocate( tkelebuoy )
    deallocate( t2leadv )
    deallocate( t2legrad )
    deallocate( t2lediff )
    deallocate( t2leprec )
    deallocate( t2lediss )
    deallocate( q2leadv )
    deallocate( q2legrad )
    deallocate( q2lediff )
    deallocate( q2leprec )
    deallocate( q2lediss )
    deallocate( twleadv )
    deallocate( twlediff )
    deallocate( twlepres )
    deallocate( twlebuoy )
    deallocate( twleprec )
    deallocate( qwleadv )
    deallocate( qwlediff )
    deallocate( qwlepres )
    deallocate( qwlebuoy )
    deallocate( qwleprec )
    deallocate( momleadv )
    deallocate( momlepress )
    deallocate( momlebuoy )
    deallocate( momlediff )
    deallocate( tadv )
    deallocate( tdiff )
    deallocate( tlat )
    deallocate( tlatqi )
    deallocate( qifall )
    deallocate( qpfall )
    deallocate( tdiff_xy )
    deallocate( tdiff_z )
    deallocate( ttest0 )
    deallocate( ttest1 )
    deallocate( ttest2 )
    deallocate( qtotmicro )
    deallocate( qlsvadv )
    deallocate( tlsvadv )
    deallocate( ulsvadv )
    deallocate( vlsvadv )
    deallocate( qnudge )
    deallocate( tnudge )
    deallocate( unudge )
    deallocate( vnudge )
    deallocate( qstor )
    deallocate( tstor )
    deallocate( ustor )
    deallocate( vstor )
    deallocate( qtostor )
    deallocate( utendcor )
    deallocate( vtendcor )
    deallocate( CF3D )
    deallocate( u850_xy )
    deallocate( v850_xy )
    deallocate( psfc_xy )
    deallocate( swvp_xy )
    deallocate( cloudtopheight )
    deallocate( echotopheight )
    deallocate( cloudtoptemp )
#if (defined CRM && defined MODAL_AERO)
    deallocate( naer )
    deallocate( vaer )
    deallocate( hgaer  )
#endif
end subroutine deallocate_vars


end module vars
