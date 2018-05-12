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

  real(crm_rknd) :: u   (dimx1_u:dimx2_u,dimy1_u:dimy2_u,nzm) ! x-wind
  real(crm_rknd) :: v   (dimx1_v:dimx2_v,dimy1_v:dimy2_v,nzm) ! y-wind
  real(crm_rknd) :: w   (dimx1_w:dimx2_w,dimy1_w:dimy2_w,nz ) ! z-wind
  real(crm_rknd) :: t   (dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm) ! liquid/ice water static energy

  !--------------------------------------------------------------------
  ! diagnostic variables:

  real(crm_rknd) :: p       (0:nx, (1-YES3D):ny, nzm)     ! perturbation pressure (from Poison eq)
  real(crm_rknd) :: tabs    (nx, ny, nzm)                 ! temperature
  real(crm_rknd) :: qv      (nx, ny, nzm)                ! water vapor
  real(crm_rknd) :: qcl     (nx, ny, nzm)                ! liquid water  (condensate)
  real(crm_rknd) :: qpl     (nx, ny, nzm)                ! liquid water  (precipitation)
  real(crm_rknd) :: qci     (nx, ny, nzm)                ! ice water  (condensate)
  real(crm_rknd) :: qpi     (nx, ny, nzm)                ! ice water  (precipitation)
  real(crm_rknd) :: tke2(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)   ! SGS TKE
  real(crm_rknd) :: tk2  (0:nxp1, (1-YES3D):nyp1, nzm) ! SGS eddyviscosity

  !--------------------------------------------------------------------
  ! time-tendencies for prognostic variables

  real(crm_rknd) :: dudt   (nxp1, ny, nzm, 3)
  real(crm_rknd) :: dvdt   (nx, nyp1, nzm, 3)
  real(crm_rknd) :: dwdt   (nx, ny  , nz,  3)

  !----------------------------------------------------------------
  ! Temporary storage array:

  real(crm_rknd) :: misc(nx, ny, nz)
  !------------------------------------------------------------------
  ! fluxes at the top and bottom of the domain:

  real(crm_rknd) :: fluxbu  (nx,ny)
  real(crm_rknd) :: fluxbv  (nx,ny)
  real(crm_rknd) :: fluxbt  (nx,ny)
  real(crm_rknd) :: fluxbq  (nx,ny)
  real(crm_rknd) :: fluxtu  (nx,ny)
  real(crm_rknd) :: fluxtv  (nx,ny)
  real(crm_rknd) :: fluxtt  (nx,ny)
  real(crm_rknd) :: fluxtq  (nx,ny)
  real(crm_rknd) :: fzero   (nx,ny)
  real(crm_rknd) :: precsfc (nx,ny) ! surface precip. rate
  real(crm_rknd) :: precssfc(nx,ny) ! surface ice precip. rate

  !-----------------------------------------------------------------
  ! profiles

  real(crm_rknd) :: t0   (nzm)
  real(crm_rknd) :: q0   (nzm)
  real(crm_rknd) :: qv0  (nzm)
  real(crm_rknd) :: tabs0(nzm)
  real(crm_rknd) :: tl0  (nzm)
  real(crm_rknd) :: tv0  (nzm)
  real(crm_rknd) :: u0   (nzm)
  real(crm_rknd) :: v0   (nzm)
  real(crm_rknd) :: tg0  (nzm)
  real(crm_rknd) :: qg0  (nzm)
  real(crm_rknd) :: ug0  (nzm)
  real(crm_rknd) :: vg0  (nzm)
  real(crm_rknd) :: p0   (nzm)
  real(crm_rknd) :: tke0 (nzm)
  real(crm_rknd) :: t01  (nzm)
  real(crm_rknd) :: q01  (nzm)
  real(crm_rknd) :: qp0  (nzm)
  real(crm_rknd) :: qn0  (nzm)
  !----------------------------------------------------------------
  ! "observed" (read from snd file) surface characteristics

  real(crm_rknd)  sstobs, lhobs, shobs
  !----------------------------------------------------------------
  !  Domain top stuff:

  real(crm_rknd)   gamt0    ! gradient of t() at the top,K/m
  real(crm_rknd)   gamq0    ! gradient of q() at the top,g/g/m

  !-----------------------------------------------------------------
  ! reference vertical profiles:

  real(crm_rknd) :: prespot(nzm)  ! (1000./pres)**R/cp
  real(crm_rknd) :: rho    (nzm)	  ! air density at pressure levels,kg/m3
  real(crm_rknd) :: rhow   (nz )   ! air density at vertical velocity levels,kg/m3
  real(crm_rknd) :: bet    (nzm)	  ! = ggr/tv0
  real(crm_rknd) :: gamaz  (nzm) ! ggr/cp*z
  real(crm_rknd) :: wsub   (nz )   ! Large-scale subsidence velocity,m/s
  real(crm_rknd) :: qtend  (nzm) ! Large-scale tendency for total water
  real(crm_rknd) :: ttend  (nzm) ! Large-scale tendency for temp.
  real(crm_rknd) :: utend  (nzm) ! Large-scale tendency for u
  real(crm_rknd) :: vtend  (nzm) ! Large-scale tendency for v

  !---------------------------------------------------------------------
  ! Large-scale and surface forcing:

  integer nlsf	! number of large-scale forcing profiles
  integer nrfc	! number of radiative forcing profiles
  integer nsfc	! number of surface forcing profiles
  integer nsnd	! number of observed soundings
  integer nzlsf	! number of large-scale forcing profiles
  integer nzrfc	! number of radiative forcing profiles
  integer nzsnd	! number of observed soundings

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
  real(crm_rknd) sstxy    (0:nx,(1-YES3D):ny)	!  surface temperature xy-distribution
  real(crm_rknd) fcory    (0:ny)      !  Coriolis parameter xy-distribution
  real(crm_rknd) fcorzy   (ny)      !  z-Coriolis parameter xy-distribution
  real(crm_rknd) latitude (nx,ny)	     ! latitude (degrees)
  real(crm_rknd) longitude(nx,ny)	     ! longitude(degrees)
  real(crm_rknd) prec_xy  (nx,ny) ! mean precip. rate for outout
  real(crm_rknd) shf_xy   (nx,ny) ! mean precip. rate for outout
  real(crm_rknd) lhf_xy   (nx,ny) ! mean precip. rate for outout
  real(crm_rknd) lwns_xy  (nx,ny) ! mean net lw at SFC
  real(crm_rknd) swns_xy  (nx,ny) ! mean net sw at SFC
  real(crm_rknd) lwnsc_xy (nx,ny) ! clear-sky mean net lw at SFC
  real(crm_rknd) swnsc_xy (nx,ny) ! clear-sky mean net sw at SFC
  real(crm_rknd) lwnt_xy  (nx,ny) ! mean net lw at TOA
  real(crm_rknd) swnt_xy  (nx,ny) ! mean net sw at TOA
  real(crm_rknd) lwntc_xy (nx,ny) ! clear-sky mean net lw at TOA
  real(crm_rknd) swntc_xy (nx,ny) ! clear-sky mean net sw at TOA
  real(crm_rknd) solin_xy (nx,ny) ! solar TOA insolation
  real(crm_rknd) pw_xy    (nx,ny)   ! precipitable water
  real(crm_rknd) cw_xy    (nx,ny)   ! cloud water path
  real(crm_rknd) iw_xy    (nx,ny)   ! ice water path
  real(crm_rknd) cld_xy   (nx,ny)   ! cloud frequency
  real(crm_rknd) u200_xy  (nx,ny) ! u-wind at 200 mb
  real(crm_rknd) usfc_xy  (nx,ny) ! u-wind at at the surface
  real(crm_rknd) v200_xy  (nx,ny) ! v-wind at 200 mb
  real(crm_rknd) vsfc_xy  (nx,ny) ! v-wind at the surface
  real(crm_rknd) w500_xy  (nx,ny) ! w at 500 mb
  real(crm_rknd) qocean_xy(nx,ny) ! ocean cooling in W/m2

  !----------------------------------------------------------------------
  !	Vertical profiles of quantities sampled for statitistics purposes:

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
  real(crm_rknd) :: twle(nz)
  real(crm_rknd) :: twsb(nz)
  real(crm_rknd) :: precflux(nz)
  real(crm_rknd) :: uwle(nz)
  real(crm_rknd) :: uwsb(nz)
  real(crm_rknd) :: vwle(nz)
  real(crm_rknd) :: vwsb(nz)
  real(crm_rknd) :: radlwup(nz)
  real(crm_rknd) :: radlwdn(nz)
  real(crm_rknd) :: radswup(nz)
  real(crm_rknd) :: radswdn(nz)
  real(crm_rknd) :: radqrlw(nz)
  real(crm_rknd) :: radqrsw(nz)
  real(crm_rknd) :: tkeleadv(nz)
  real(crm_rknd) :: tkelepress(nz)
  real(crm_rknd) :: tkelediss(nz)
  real(crm_rknd) :: tkelediff(nz)
  real(crm_rknd) :: tkelebuoy(nz)
  real(crm_rknd) :: t2leadv(nz)
  real(crm_rknd) :: t2legrad(nz)
  real(crm_rknd) :: t2lediff(nz)
  real(crm_rknd) :: t2leprec(nz)
  real(crm_rknd) :: t2lediss(nz)
  real(crm_rknd) :: q2leadv(nz)
  real(crm_rknd) :: q2legrad(nz)
  real(crm_rknd) :: q2lediff(nz)
  real(crm_rknd) :: q2leprec(nz)
  real(crm_rknd) :: q2lediss(nz)
  real(crm_rknd) :: twleadv(nz)
  real(crm_rknd) :: twlediff(nz)
  real(crm_rknd) :: twlepres(nz)
  real(crm_rknd) :: twlebuoy(nz)
  real(crm_rknd) :: twleprec(nz)
  real(crm_rknd) :: qwleadv(nz)
  real(crm_rknd) :: qwlediff(nz)
  real(crm_rknd) :: qwlepres(nz)
  real(crm_rknd) :: qwlebuoy(nz)
  real(crm_rknd) :: qwleprec(nz)
  real(crm_rknd) :: momleadv(nz,3)
  real(crm_rknd) :: momlepress(nz,3)
  real(crm_rknd) :: momlebuoy(nz,3)
  real(crm_rknd) :: momlediff(nz,3)
  real(crm_rknd) :: tadv(nz)
  real(crm_rknd) :: tdiff(nz)
  real(crm_rknd) :: tlat(nz)
  real(crm_rknd) :: tlatqi(nz)
  real(crm_rknd) :: qifall(nz)
  real(crm_rknd) :: qpfall(nz)
  real(crm_rknd) :: tdiff_xy(nz)
  real(crm_rknd) :: tdiff_z(nz)
  real(crm_rknd) :: ttest0(nzm)
  real(crm_rknd) :: ttest1(nz)
  real(crm_rknd) :: ttest2(nz,10)  !+++mhwang test

  ! energy conservation diagnostics:

  real(8) total_water_before, total_water_after
  real(8) total_water_evap, total_water_prec, total_water_ls
  !#ifdef CLUBB_CRM
  real(8) total_water_clubb
  real(8) total_energy_before, total_energy_after
  real(8) total_energy_evap, total_energy_prec, total_energy_ls
  real(8) total_energy_clubb, total_energy_rad
  !#endif
  real(8) :: qtotmicro(5)  ! total water for water conservation test in microphysics +++mhwang

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

  real(crm_rknd) :: qlsvadv(nzm) ! Large-scale vertical advection tendency for total water
  real(crm_rknd) :: tlsvadv(nzm) ! Large-scale vertical advection tendency for temperature
  real(crm_rknd) :: ulsvadv(nzm) ! Large-scale vertical advection tendency for zonal velocity
  real(crm_rknd) :: vlsvadv(nzm) ! Large-scale vertical advection tendency for meridional velocity
  real(crm_rknd) :: qnudge(nzm) ! Nudging of horiz.-averaged total water profile
  real(crm_rknd) :: tnudge(nzm) ! Nudging of horiz.-averaged temperature profile
  real(crm_rknd) :: unudge(nzm) ! Nudging of horiz.-averaged zonal velocity
  real(crm_rknd) :: vnudge(nzm) ! Nudging of horiz.-averaged meridional velocity
  real(crm_rknd) :: qstor(nzm) ! Storage of horiz.-averaged total water profile
  real(crm_rknd) :: tstor(nzm) ! Storage of horiz.-averaged temperature profile
  real(crm_rknd) :: ustor(nzm) ! Storage of horiz.-averaged zonal velocity
  real(crm_rknd) :: vstor(nzm) ! Storage of horiz.-averaged meridional velocity
  real(crm_rknd) :: qtostor(nzm) ! Storage of horiz.-averaged total water profile (vapor + liquid)
  real(crm_rknd) :: utendcor(nzm) ! coriolis acceleration of zonal velocity
  real(crm_rknd) :: vtendcor(nzm) ! coriolis acceleration of meridional velocity
  real(crm_rknd) :: CF3D(1:nx, 1:ny, 1:nzm)  ! Cloud fraction
  ! =1.0 when there is no fractional cloudiness scheme
  ! = cloud fraction produced by fractioal cloudiness scheme when avaiable

  ! 850 mbar horizontal winds
  real(crm_rknd) :: u850_xy(nx,ny) ! zonal velocity at 850 mb
  real(crm_rknd) :: v850_xy(nx,ny) ! meridional velocity at 850 mb

  ! Surface pressure
  real(crm_rknd) :: psfc_xy(nx,ny) ! pressure (in millibar) at lowest grid point

  ! Saturated water vapor path, useful for computing column relative humidity
  real(crm_rknd) :: swvp_xy(nx,ny)  ! saturated water vapor path (wrt water)

  ! Cloud and echo top heights, and cloud top temperature (instantaneous)
  real(crm_rknd) :: cloudtopheight(nx,ny)
  real(crm_rknd) :: echotopheight (nx,ny)
  real(crm_rknd) :: cloudtoptemp  (nx,ny)

  ! END UW ADDITIONS
  !===========================================================================
  ! Initial bubble parameters. Activated when perturb_type = 2
  real(crm_rknd) bubble_x0
  real(crm_rknd) bubble_y0
  real(crm_rknd) bubble_z0
  real(crm_rknd) bubble_radius_hor
  real(crm_rknd) bubble_radius_ver
  real(crm_rknd) bubble_dtemp
  real(crm_rknd) bubble_dq
#if (defined CRM && defined MODAL_AERO)
  real(crm_rknd) :: naer (nzm, ntot_amode)     ! Aerosol number concentration [/m3]
  real(crm_rknd) :: vaer (nzm, ntot_amode)     ! aerosol volume concentration [m3/m3]
  real(crm_rknd) :: hgaer(nzm, ntot_amode)    ! hygroscopicity of aerosol mode
#endif

end module vars
