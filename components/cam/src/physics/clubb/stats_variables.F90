!-------------------------------------------------------------------------------
! $Id$
!-------------------------------------------------------------------------------

! Description:
!   Holds pointers and other variables for statistics to be written to 
!   GrADS files and netCDF files.
!-------------------------------------------------------------------------------
module stats_variables


  use stats_type, only:  &
      stats ! Type

  use clubb_precision, only:  & 
      core_rknd ! Variable(s)

  implicit none

  private ! Set Default Scope

  ! Sampling and output frequencies
  real( kind = core_rknd ), public :: &
    stats_tsamp = 0._core_rknd, & ! Sampling interval   [s]
    stats_tout  = 0._core_rknd ! Output interval     [s]


  logical, public ::  &
    l_stats            = .false., & ! Main flag to turn statistics on/off
    l_output_rad_files = .false., & ! Flag to turn off radiation statistics output
    l_netcdf           = .false., & ! Output to NetCDF format
    l_grads            = .false., & ! Output to GrADS format
    l_silhs_out        = .false., & ! Output SILHS files (stats_lh_zt and stats_lh_sfc)
    l_allow_small_stats_tout = .false. ! Do not stop if output timestep is too low for
                      ! requested format, e.g. l_grads = .true. and
                      ! stats_tout < 60.0


  logical, public :: & 
    l_stats_samp = .false., & ! Sample flag for current time step
    l_stats_last = .false.    ! Last time step of output period


  character(len=200), public ::  & 
    fname_zt     = '', & ! Name of the stats file for thermodynamic grid fields
    fname_lh_zt  = '', & ! Name of the stats file for LH variables on the stats_zt grid
    fname_lh_sfc = '', & ! Name of the stats file for LH variables on the stats_zt grid
    fname_zm     = '', & ! Name of the stats file for momentum grid fields
    fname_rad_zt = '', & ! Name of the stats file for the stats_zt radiation grid fields
    fname_rad_zm = '', & ! Name of the stats file for the stats_zm radiation grid fields
    fname_sfc    = ''    ! Name of the stats file for surface only fields


!       Indices for statistics in stats_zt file

  integer, public :: & 
     ithlm = 0, & 
     ithvm = 0, & 
     irtm = 0, & 
     ircm = 0, &
     irvm = 0, & 
     ium = 0, & 
     ivm = 0, & 
     iwm_zt = 0, &
     iwm_zm = 0, &
     ium_ref = 0,&
     ivm_ref = 0, & 
     iug = 0, & 
     ivg = 0, & 
     icloud_frac = 0, &
     iice_supersat_frac = 0, &
     ircm_in_layer = 0, &
     ircm_in_cloud = 0, &
     icloud_cover = 0, &
     ip_in_Pa = 0, & 
     iexner = 0, & 
     irho_ds_zt = 0, &
     ithv_ds_zt = 0, &
     iLscale = 0, & 
     iwp3 = 0, & 
     ithlp3 = 0, &
     irtp3 = 0, &
     iwpthlp2 = 0, & 
     iwp2thlp = 0, & 
     iwprtp2 = 0, & 
     iwp2rtp = 0, &
     iSkw_zt = 0, &
     iSkthl_zt = 0, &
     iSkrt_zt = 0, &
     ircm_supersat_adj = 0

  integer, public :: & 
     iLscale_up = 0, & 
     iLscale_down = 0, & 
     iLscale_pert_1 = 0, & 
     iLscale_pert_2 = 0, & 
     itau_zt = 0, & 
     iKh_zt = 0, & 
     iwp2thvp = 0, & 
     iwp2rcp = 0, & 
     iwprtpthlp = 0, &
     irc_coef = 0, &
     isigma_sqd_w_zt = 0, & 
     irho = 0

  integer, public :: &
     itau_zm_simp = 0


  integer, dimension(:), allocatable, public :: & 
     ihm_1, &
     ihm_2

  integer, public :: & 
     iprecip_frac = 0, &
     iprecip_frac_1 = 0, &
     iprecip_frac_2 = 0, &
     iNcnm = 0 

  integer, dimension(:), allocatable, public :: &
     imu_hm_1,         &
     imu_hm_2,         &
     imu_hm_1_n,       &
     imu_hm_2_n,       &
     isigma_hm_1,      &
     isigma_hm_2,      &
     isigma_hm_1_n,    &
     isigma_hm_2_n,    &
     icorr_w_hm_1,     &
     icorr_w_hm_2,     &
     icorr_chi_hm_1,   &
     icorr_chi_hm_2,   &
     icorr_eta_hm_1,   &
     icorr_eta_hm_2,   &
     icorr_Ncn_hm_1,   &
     icorr_Ncn_hm_2,   &
     icorr_w_hm_1_n,   &
     icorr_w_hm_2_n,   &
     icorr_chi_hm_1_n, &
     icorr_chi_hm_2_n, &
     icorr_eta_hm_1_n, &
     icorr_eta_hm_2_n, &
     icorr_Ncn_hm_1_n, &
     icorr_Ncn_hm_2_n

  integer, dimension(:,:), allocatable, public :: &
     icorr_hmx_hmy_1,   &
     icorr_hmx_hmy_2,   &
     icorr_hmx_hmy_1_n, &
     icorr_hmx_hmy_2_n

  integer, public :: &
     imu_Ncn_1 = 0,      &
     imu_Ncn_2 = 0,      &
     imu_Ncn_1_n = 0,    &
     imu_Ncn_2_n = 0,    &
     isigma_Ncn_1 = 0,   &
     isigma_Ncn_2 = 0,   &
     isigma_Ncn_1_n = 0, &
     isigma_Ncn_2_n = 0

  integer, public :: &
     icorr_w_chi_1_ca = 0, &
     icorr_w_chi_2_ca = 0, &
     icorr_w_eta_1_ca = 0, &
     icorr_w_eta_2_ca = 0, &
     icorr_w_Ncn_1 = 0,  &
     icorr_w_Ncn_2 = 0,  &
     icorr_chi_eta_1_ca = 0, &
     icorr_chi_eta_2_ca = 0, &
     icorr_chi_Ncn_1 = 0,  &
     icorr_chi_Ncn_2 = 0,  &
     icorr_eta_Ncn_1 = 0,  &
     icorr_eta_Ncn_2 = 0

  integer, public :: &
     icorr_w_Ncn_1_n = 0, &
     icorr_w_Ncn_2_n = 0, &
     icorr_chi_Ncn_1_n = 0, &
     icorr_chi_Ncn_2_n = 0, &
     icorr_eta_Ncn_1_n = 0, &
     icorr_eta_Ncn_2_n = 0

  integer, dimension(:), allocatable, public :: &
    isilhs_variance_category, &
    ilh_samp_frac_category


  integer, public :: & 
     iNcm = 0,             & ! Brian
     iNccnm = 0,           & 
     iNc_in_cloud = 0,     &
     iNc_activated = 0,    &
     isnowslope = 0,       & ! Adam Smith, 22 April 2008
     ised_rcm = 0,         & ! Brian
     irsat = 0,            & ! Brian
     irsati = 0,           & 
     irrm = 0,          & ! Brian
     im_vol_rad_rain = 0,  & ! Brian
     im_vol_rad_cloud = 0, & ! COAMPS only. dschanen 6 Dec 2006
     iprecip_rate_zt = 0,    & ! Brian
     iAKm = 0,             & ! analytic Kessler.  Vince Larson 22 May 2005 
     ilh_AKm = 0,          & ! LH Kessler.  Vince Larson  22 May 2005
     iradht = 0,           & ! Radiative heating.
     iradht_LW = 0,        & !   "           "   Long-wave component
     iradht_SW = 0,        & !   "           "   Short-wave component
     irel_humidity = 0

  integer, public :: & 
     iAKstd = 0,     &
     iAKstd_cld = 0, & 
     iAKm_rcm = 0, & 
     iAKm_rcc = 0


  integer, public :: & 
   irfrzm = 0

  ! Skewness functions on stats_zt grid
  integer, public :: &
    iC11_Skw_fnc = 0


  integer, public :: &
    icloud_frac_zm = 0, &
    iice_supersat_frac_zm = 0, &
    ircm_zm = 0, &
    irtm_zm = 0, &
    ithlm_zm = 0


  integer, public :: &
    iw_1_zm = 0, &
    iw_2_zm = 0, &
    ivarnce_w_1_zm = 0, &
    ivarnce_w_2_zm = 0, &
    imixt_frac_zm = 0


  integer, public :: &
    ilh_rcm_avg = 0, &
    ik_lh_start = 0


  integer, public :: & 
     iNrm = 0,       & ! Rain droplet number concentration
     iNim = 0,       & ! Ice number concentration
     iNsm = 0,    & ! Snow number concentration
     iNgm = 0    ! Graupel number concentration

  integer, public :: & 
     iT_in_K      ! Absolute temperature

  integer, public :: &
    ieff_rad_cloud = 0, &
    ieff_rad_ice = 0, &
    ieff_rad_snow = 0, &
    ieff_rad_rain = 0, &
    ieff_rad_graupel = 0


  integer, public :: & 
    irsm = 0, &
    irgm = 0, & 
    irim = 0, & 
    idiam = 0,           & ! Diameter of ice crystal           [m]
    imass_ice_cryst = 0, & ! Mass of a single ice crystal      [kg]
    ircm_icedfs = 0,     & ! Change in liquid water due to ice [kg/kg/s]
    iu_T_cm = 0         ! Fallspeed of ice crystal in cm/s  [cm s^{-1}]



  ! thlm/rtm budget terms
  integer, public :: & 
    irtm_bt = 0,         & ! rtm total time tendency
    irtm_ma = 0,         & ! rtm mean advect. term
    irtm_ta = 0,         & ! rtm turb. advect. term
    irtm_forcing = 0,    & ! rtm large scale forcing term
    irtm_mc = 0,         & ! rtm change from microphysics
    irtm_sdmp = 0,       & ! rtm change from sponge damping
    irvm_mc = 0,         & ! rvm change from microphysics
    ircm_mc = 0,         & ! rcm change from microphysics
    ircm_sd_mg_morr = 0, & ! rcm sedimentation tendency
    irtm_mfl = 0,        & ! rtm change due to monotonic flux limiter
    irtm_tacl = 0,       & ! rtm correction from turbulent advection (wprtp) clipping
    irtm_cl = 0,         & ! rtm clipping term
    irtm_pd = 0,         & ! thlm postive definite adj term
    ithlm_bt = 0,        & ! thlm total time tendency
    ithlm_ma = 0,        & ! thlm mean advect. term
    ithlm_ta = 0,        & ! thlm turb. advect. term
    ithlm_forcing = 0,   & ! thlm large scale forcing term
    ithlm_sdmp = 0,      & ! thlm change from sponge damping
    ithlm_mc = 0,        & ! thlm change from microphysics
    ithlm_mfl = 0,       & ! thlm change due to monotonic flux limiter
    ithlm_tacl = 0,      & ! thlm correction from turbulent advection (wpthlp) clipping
    ithlm_cl = 0           ! thlm clipping term


  !monatonic flux limiter diagnostic terms
  integer, public :: &
    ithlm_mfl_min = 0, &
    ithlm_mfl_max = 0, &
    iwpthlp_entermfl = 0, &
    iwpthlp_exit_mfl = 0, &
    iwpthlp_mfl_min = 0, &
    iwpthlp_mfl_max = 0, &
    irtm_mfl_min = 0, &
    irtm_mfl_max = 0, &
    iwprtp_enter_mfl = 0, &
    iwprtp_exit_mfl = 0, &
    iwprtp_mfl_min = 0, &
    iwprtp_mfl_max = 0, &
    ithlm_enter_mfl = 0, &
    ithlm_exit_mfl = 0, &
    ithlm_old = 0, &
    ithlm_without_ta = 0, &
    irtm_enter_mfl = 0, &
    irtm_exit_mfl = 0, &
    irtm_old = 0, &
    irtm_without_ta = 0


  integer, public :: & 
     iwp3_bt  = 0, & 
     iwp3_ma  = 0, & 
     iwp3_ta  = 0, & 
     iwp3_tp  = 0, & 
     iwp3_ac  = 0, & 
     iwp3_bp1 = 0, & 
     iwp3_bp2 = 0, & 
     iwp3_pr1 = 0, & 
     iwp3_pr2 = 0, & 
     iwp3_pr3 = 0, &
     iwp3_dp1 = 0, &
     iwp3_sdmp = 0, &
     iwp3_cl  = 0, &
     iwp3_splat = 0


  integer, public :: &
    irtp3_bt  = 0, &
    irtp3_tp  = 0, &
    irtp3_ac  = 0, &
    irtp3_dp  = 0, &
    ithlp3_bt = 0, &
    ithlp3_tp = 0, &
    ithlp3_ac = 0, &
    ithlp3_dp = 0


  ! Rain mixing ratio budgets
  integer, public :: & 
     irrm_bt = 0, &
     irrm_ma = 0, &
     irrm_ta = 0, &
     irrm_sd = 0, &
     irrm_ts = 0, &
     irrm_sd_morr = 0, &
     irrm_cond = 0, &
     irrm_auto = 0, &
     irrm_accr = 0, &
     irrm_cond_adj = 0, &
     irrm_src_adj = 0, &
     irrm_mc_nonadj = 0, &
     irrm_mc = 0, &
     irrm_hf = 0, &
     irrm_wvhf = 0, &
     irrm_cl = 0


  integer, public :: &
     iNrm_bt = 0, &
     iNrm_ma = 0, &
     iNrm_ta = 0, &
     iNrm_sd = 0, &
     iNrm_ts = 0, &
     iNrm_cond = 0, &
     iNrm_auto = 0, &
     iNrm_cond_adj = 0, &
     iNrm_src_adj = 0, &
     iNrm_mc = 0, &
     iNrm_cl = 0



  ! Snow/Ice/Graupel mixing ratio budgets
  integer, public :: & 
     irsm_bt = 0, & 
     irsm_ma = 0, & 
     irsm_sd = 0, & 
     irsm_sd_morr = 0, &
     irsm_ta = 0, & 
     irsm_mc = 0, & 
     irsm_hf = 0, &
     irsm_wvhf = 0, &
     irsm_cl = 0, &
     irsm_sd_morr_int = 0


  integer, public :: & 
     irgm_bt = 0, & 
     irgm_ma = 0, & 
     irgm_sd = 0, & 
     irgm_sd_morr = 0, &
     irgm_ta = 0, & 
     irgm_mc = 0, & 
     irgm_hf = 0, &
     irgm_wvhf = 0, &
     irgm_cl = 0


  integer, public :: & 
     irim_bt = 0, & 
     irim_ma = 0, & 
     irim_sd = 0, & 
     irim_sd_mg_morr = 0, &
     irim_ta = 0, & 
     irim_mc = 0, & 
     irim_hf = 0, &
     irim_wvhf = 0, &
     irim_cl = 0


  integer, public :: &
    iNsm_bt = 0,  &
    iNsm_ma = 0,  &
    iNsm_sd = 0,  &
    iNsm_ta = 0,  &
    iNsm_mc = 0,  &
    iNsm_cl = 0


  integer, public :: &
    iNgm_bt = 0, &
    iNgm_ma = 0, &
    iNgm_sd = 0, &
    iNgm_ta = 0, &
    iNgm_mc = 0, &
    iNgm_cl = 0


  integer, public :: &
    iNim_bt = 0, &
    iNim_ma = 0, &
    iNim_sd = 0, &
    iNim_ta = 0, &
    iNim_mc = 0, &
    iNim_cl = 0


  integer, public :: &
    iNcm_bt = 0, &
    iNcm_ma = 0, &
    iNcm_ta = 0, &
    iNcm_mc = 0, &
    iNcm_cl = 0, &
    iNcm_act = 0


  ! Covariances between w, r_t, theta_l and KK microphysics tendencies.
  ! Additionally, covariances between r_r and N_r and KK rain drop mean
  ! volume radius.  These are all calculated on thermodynamic grid levels.
  integer, public :: &
    iw_KK_evap_covar_zt = 0,   & ! Covariance of w and KK evaporation tendency.
    irt_KK_evap_covar_zt = 0,  & ! Covariance of r_t and KK evaporation tendency.
    ithl_KK_evap_covar_zt = 0, & ! Covariance of theta_l and KK evap. tendency.
    iw_KK_auto_covar_zt = 0,   & ! Covariance of w and KK autoconversion tendency.
    irt_KK_auto_covar_zt = 0,  & ! Covariance of r_t and KK autoconversion tendency.
    ithl_KK_auto_covar_zt = 0, & ! Covariance of theta_l and KK autoconv. tendency.
    iw_KK_accr_covar_zt = 0,   & ! Covariance of w and KK accretion tendency.
    irt_KK_accr_covar_zt = 0,  & ! Covariance of r_t and KK accretion tendency.
    ithl_KK_accr_covar_zt = 0, & ! Covariance of theta_l and KK accretion tendency.
    irr_KK_mvr_covar_zt = 0,   & ! Covariance of r_r and KK mean volume radius.
    iNr_KK_mvr_covar_zt = 0,   & ! Covariance of N_r and KK mean volume radius.
    iKK_mvr_variance_zt = 0      ! Variance of KK rain drop mean volume radius.


  ! Wind budgets
  integer, public :: & 
     ivm_bt = 0, & 
     ivm_ma = 0, & 
     ivm_ta = 0, & 
     ivm_gf = 0, & 
     ivm_cf = 0, &
     ivm_f = 0, &
     ivm_sdmp = 0, &
     ivm_ndg = 0, &
     ivm_mfl = 0


  integer, public :: & 
     ium_bt = 0, & 
     ium_ma = 0, & 
     ium_ta = 0, & 
     ium_gf = 0, & 
     ium_cf = 0, & 
     ium_f = 0, &
     ium_sdmp = 0, &
     ium_ndg = 0, &
     ium_mfl = 0



  ! PDF parameters
  integer, public :: & 
     imixt_frac = 0, & 
     iw_1 = 0, & 
     iw_2 = 0, & 
     ivarnce_w_1 = 0, & 
     ivarnce_w_2 = 0, & 
     ithl_1 = 0, & 
     ithl_2 = 0, & 
     ivarnce_thl_1 = 0, & 
     ivarnce_thl_2 = 0, & 
     irt_1 = 0, & 
     irt_2 = 0, & 
     ivarnce_rt_1 = 0, & 
     ivarnce_rt_2 = 0, & 
     irc_1 = 0, & 
     irc_2 = 0, & 
     irsatl_1 = 0, & 
     irsatl_2 = 0, & 
     icloud_frac_1 = 0, & 
     icloud_frac_2 = 0

  integer, public :: & 
     ichi_1 = 0, &
     ichi_2 = 0, &
     istdev_chi_1 = 0, & 
     istdev_chi_2 = 0, &
     ichip2 = 0, &
     istdev_eta_1 = 0, &
     istdev_eta_2 = 0, &
     icovar_chi_eta_1 = 0, &
     icovar_chi_eta_2 = 0, &
     icorr_w_chi_1 = 0, &
     icorr_w_chi_2 = 0, &
     icorr_w_eta_1 = 0, &
     icorr_w_eta_2 = 0, &
     icorr_chi_eta_1 = 0, &
     icorr_chi_eta_2 = 0, &
     icorr_w_rt_1 = 0, &
     icorr_w_rt_2 = 0, &
     icorr_w_thl_1 = 0, &
     icorr_w_thl_2 = 0, &
     icorr_rt_thl_1 = 0, &
     icorr_rt_thl_2 = 0, &
     icrt_1 = 0, &
     icrt_2 = 0, &
     icthl_1 = 0, &
     icthl_2 = 0

  integer, public :: &
    iF_w = 0, &
    iF_rt = 0, &
    iF_thl = 0, &
    imin_F_w = 0, &
    imax_F_w = 0, &
    imin_F_rt = 0, &
    imax_F_rt = 0, &
    imin_F_thl = 0, &
    imax_F_thl = 0


  integer, public :: &
    icoef_wprtp2_implicit = 0, &
    iterm_wprtp2_explicit = 0, &
    icoef_wpthlp2_implicit = 0, &
    iterm_wpthlp2_explicit = 0, &
    icoef_wprtpthlp_implicit = 0, &
    iterm_wprtpthlp_explicit = 0, &
    icoef_wp2rtp_implicit = 0, &
    iterm_wp2rtp_explicit = 0, &
    icoef_wp2thlp_implicit = 0, &
    iterm_wp2thlp_explicit = 0


  integer, public :: & 
    iwp2_zt = 0, & 
    ithlp2_zt = 0, & 
    iwpthlp_zt = 0, & 
    iwprtp_zt = 0, & 
    irtp2_zt = 0, & 
    irtpthlp_zt = 0, &
    iup2_zt = 0, &
    ivp2_zt = 0, &
    iupwp_zt = 0, &
    ivpwp_zt = 0


  integer, dimension(:), allocatable, public :: &
    iwp2hmp


  integer, dimension(:), allocatable, public :: &
    ihydrometp2,  &
    iwphydrometp, &
    irtphmp,      &
    ithlphmp


  integer, dimension(:,:), allocatable, public :: &
    ihmxphmyp


  integer, dimension(:), allocatable, public :: &
    ihmp2_zt


  integer, public :: &
    ichi = 0

  integer, target, allocatable, dimension(:), public :: & 
    isclrm,   & ! Passive scalar mean (1)
    isclrm_f    ! Passive scalar forcing (1)

! Used to calculate clear-sky radiative fluxes.
  integer, public :: &
    ifulwcl = 0, ifdlwcl = 0, ifdswcl = 0, ifuswcl = 0


  integer, target, allocatable, dimension(:), public :: & 
    iedsclrm,   & ! Eddy-diff. scalar term (1)
    iedsclrm_f    ! Eddy-diffusivity scalar forcing (1)


  integer, public :: &
    ilh_thlm_mc = 0,      & ! Latin hypercube estimate of thlm_mc
    ilh_rvm_mc = 0,       & ! Latin hypercube estimate of rvm_mc
    ilh_rcm_mc = 0,       & ! Latin hypercube estimate of rcm_mc
    ilh_Ncm_mc = 0,       & ! Latin hypercube estimate of Ncm_mc
    ilh_rrm_mc = 0,    & ! Latin hypercube estimate of rrm_mc
    ilh_Nrm_mc = 0,       & ! Latin hypercube estimate of Nrm_mc
    ilh_rsm_mc = 0,    & ! Latin hypercube estimate of rsm_mc
    ilh_Nsm_mc = 0,    & ! Latin hypercube estimate of Nsm_mc
    ilh_rgm_mc = 0, & ! Latin hypercube estimate of rgm_mc
    ilh_Ngm_mc = 0, & ! Latin hypercube estimate of Ngm_mc
    ilh_rim_mc = 0,     & ! Latin hypercube estimate of rim_mc
    ilh_Nim_mc = 0          ! Latin hypercube estimate of Nim_mc

  integer, public :: &
    ilh_rrm_auto = 0, & ! Latin hypercube estimate of autoconversion
    ilh_rrm_accr = 0, & ! Latin hypercube estimate of accretion
    ilh_rrm_evap = 0, & ! Latin hypercube estimate of evaporation
    ilh_Nrm_auto    = 0, & ! Latin hypercube estimate of Nrm autoconversion
    ilh_Nrm_cond    = 0, & ! Latin hypercube estimate of Nrm evaporation
    ilh_m_vol_rad_rain = 0, &
    ilh_rrm_mc_nonadj = 0


  integer, public :: &
    ilh_rrm_src_adj  = 0, & ! Latin hypercube estimate of source adjustment (KK only!)
    ilh_rrm_cond_adj = 0, & ! Latin hypercube estimate of evap adjustment (KK only!)
    ilh_Nrm_src_adj     = 0, & ! Latin hypercube estimate of Nrm source adjustmet (KK only!)
    ilh_Nrm_cond_adj    = 0    ! Latin hypercube estimate of Nrm evap adjustment (KK only!)

  integer, public :: &
    ilh_Vrr = 0, & ! Latin hypercube estimate of rrm sedimentation velocity
    ilh_VNr = 0    ! Latin hypercube estimate of Nrm sedimentation velocity

  integer, public :: &
    ilh_rrm = 0, &
    ilh_Nrm = 0, &
    ilh_rim = 0, &
    ilh_Nim = 0, &
    ilh_rsm = 0, &
    ilh_Nsm = 0, &
    ilh_rgm = 0, &
    ilh_Ngm = 0, &
    ilh_thlm = 0, &
    ilh_rcm = 0, &
    ilh_Ncm = 0, &
    ilh_Ncnm = 0, &
    ilh_rvm = 0, &
    ilh_wm = 0, &
    ilh_cloud_frac = 0, &
    ilh_chi = 0, &
    ilh_eta = 0, &
    ilh_precip_frac = 0, &
    ilh_mixt_frac = 0


  integer, public :: &
    ilh_cloud_frac_unweighted  = 0,  &
    ilh_precip_frac_unweighted = 0,  &
    ilh_mixt_frac_unweighted   = 0


  integer, public :: &
    ilh_wp2_zt = 0, &
    ilh_Nrp2_zt = 0, &
    ilh_Ncnp2_zt = 0, &
    ilh_Ncp2_zt = 0, &
    ilh_rcp2_zt = 0, &
    ilh_rtp2_zt = 0, &
    ilh_thlp2_zt = 0, &
    ilh_rrp2_zt = 0, &
    ilh_chip2 = 0 ! Eric Raut

  ! SILHS covariance estimate indicies
  integer, public :: &
    ilh_rtp2_mc    = 0, &
    ilh_thlp2_mc   = 0, &
    ilh_wprtp_mc   = 0, &
    ilh_wpthlp_mc  = 0, &
    ilh_rtpthlp_mc = 0


  ! Indices for Morrison budgets
  integer, public :: &
    iPSMLT = 0, &
    iEVPMS = 0, &
    iPRACS = 0, &
    iEVPMG = 0, &
    iPRACG = 0, &
    iPGMLT = 0, &
    iMNUCCC = 0, &
    iPSACWS = 0, &
    iPSACWI = 0, &
    iQMULTS = 0, &
    iQMULTG = 0, &
    iPSACWG = 0, &
    iPGSACW = 0, &
    iPRD = 0, &
    iPRCI = 0, &
    iPRAI = 0, &
    iQMULTR = 0, &
    iQMULTRG = 0, &
    iMNUCCD = 0, &
    iPRACI = 0, &
    iPRACIS = 0, &
    iEPRD = 0, &
    iMNUCCR = 0, &
    iPIACR = 0, &
    iPIACRS = 0, &
    iPGRACS = 0, &
    iPRDS = 0, &
    iEPRDS = 0, &
    iPSACR = 0, &
    iPRDG = 0, &
    iEPRDG = 0


  ! More indices for Morrison budgets!!
  integer, public :: &
    iNGSTEN = 0, &
    iNRSTEN = 0, &
    iNISTEN = 0, &
    iNSSTEN = 0, &
    iNCSTEN = 0, &
    iNPRC1 = 0,  &
    iNRAGG = 0,  &
    iNPRACG = 0, &
    iNSUBR = 0,  &
    iNSMLTR = 0, &
    iNGMLTR = 0, &
    iNPRACS = 0, &
    iNNUCCR = 0, &
    iNIACR = 0,  &
    iNIACRS = 0, &
    iNGRACS = 0, &
    iNSMLTS = 0, &
    iNSAGG = 0,  &
    iNPRCI = 0, &
    iNSCNG = 0, &
    iNSUBS = 0, &
    iPRC = 0, &
    iPRA = 0, &
    iPRE = 0


  ! More indices for Morrison budgets!!
  integer, public :: &
    iPCC = 0, & 
    iNNUCCC = 0, & 
    iNPSACWS = 0, &
    iNPRA = 0, & 
    iNPRC = 0, & 
    iNPSACWI = 0, &
    iNPSACWG = 0, &
    iNPRAI = 0, &
    iNMULTS = 0, & 
    iNMULTG = 0, & 
    iNMULTR = 0, & 
    iNMULTRG = 0, & 
    iNNUCCD = 0, & 
    iNSUBI = 0, & 
    iNGMLTG = 0, &
    iNSUBG = 0, &
    iNACT = 0

  integer, public :: &
    iSIZEFIX_NR = 0, &
    iSIZEFIX_NC = 0, &
    iSIZEFIX_NI = 0, &
    iSIZEFIX_NS = 0, &
    iSIZEFIX_NG = 0, &
    iNEGFIX_NR = 0, &
    iNEGFIX_NC = 0, &
    iNEGFIX_NI = 0, &
    iNEGFIX_NS = 0, &
    iNEGFIX_NG = 0, &
    iNIM_MORR_CL = 0, &
    iQC_INST = 0, &
    iQR_INST = 0, &
    iQI_INST = 0, &
    iQS_INST = 0, &
    iQG_INST = 0, &
    iNC_INST = 0, &
    iNR_INST = 0, &
    iNI_INST = 0, &
    iNS_INST = 0, &
    iNG_INST = 0, &
    iT_in_K_mc = 0, &
    ihl_on_Cp_residual = 0, &
    iqto_residual = 0


  ! Indices for statistics in stats_zm file
  integer, public :: & 
     iwp2 = 0, & 
     irtp2 = 0, & 
     ithlp2 = 0, & 
     irtpthlp = 0, & 
     iwprtp = 0, & 
     iwpthlp = 0, & 
     iwp4 = 0, & 
     iwpthvp = 0, & 
     irtpthvp = 0, & 
     ithlpthvp = 0, & 
     itau_zm = 0, & 
     iKh_zm = 0, & 
     iwprcp = 0, & 
     irc_coef_zm = 0, &
     ithlprcp = 0, & 
     irtprcp = 0, & 
     ircp2 = 0, & 
     iupwp = 0, & 
     ivpwp = 0, &
     iupthlp = 0, &
     iuprtp = 0, &
     ivpthlp = 0, &
     ivprtp = 0, &
     iupthvp = 0, &
     iuprcp = 0, &
     ivpthvp = 0, &
     ivprcp = 0, &
     iSkw_zm = 0, &
     iSkthl_zm = 0, &
     iSkrt_zm = 0

  integer, public :: &
     irho_zm = 0, & 
     isigma_sqd_w = 0, &
     irho_ds_zm = 0, &
     ithv_ds_zm = 0, &
     iem = 0, & 
     ishear = 0, & ! Brian
     imean_w_up = 0, &
     imean_w_down = 0, &
     iFrad = 0, & 
     iFrad_LW = 0,   & ! Brian
     iFrad_SW = 0,   & ! Brian
     iFrad_LW_up = 0,   & 
     iFrad_SW_up = 0,   & 
     iFrad_LW_down = 0,   & 
     iFrad_SW_down = 0,   & 
     iFprec = 0,          & ! Brian
     iFcsed = 0             ! Brian

   ! Stability correction applied to Kh_N2_zm (diffusion on rtm and thlm)
   integer, public :: &
     istability_correction = 0 ! schemena


  integer, dimension(:), allocatable, public :: &
    iK_hm

  ! Skewness Functions on stats_zm grid
  integer, public :: &
    igamma_Skw_fnc = 0,  &
    iC6rt_Skw_fnc = 0,   &
    iC6thl_Skw_fnc = 0,  &
    iC7_Skw_fnc = 0,     &
    iC1_Skw_fnc = 0,     &
    ibrunt_vaisala_freq_sqd = 0, &
    iRichardson_num = 0, &
    ishear_sqd = 0


  integer, public :: &
    icoef_wp4_implicit = 0


  ! Covariance of w and cloud droplet concentration, < w'N_c' >
  integer, public :: &
    iwpNcp = 0


  ! Sedimentation velocities
  integer, public :: & 
    iVNr = 0,    &
    iVrr = 0,    &
    iVNc = 0,    &
    iVrc = 0,    &
    iVNs = 0, &
    iVrs = 0, &
    iVNi = 0,  &
    iVri = 0,  &
    iVrg = 0


  ! Covariance of sedimentation velocity and hydrometeor, <V_xx'x_x'>.
  integer, public :: &
    iVrrprrp = 0,         &
    iVNrpNrp = 0,         &
    iVrrprrp_expcalc = 0, &
    iVNrpNrp_expcalc = 0


  integer, public :: & 
     iwp2_bt = 0, & 
     iwp2_ma = 0, & 
     iwp2_ta = 0, & 
     iwp2_ac = 0, & 
     iwp2_bp = 0, & 
     iwp2_pr1 = 0, & 
     iwp2_pr2 = 0, & 
     iwp2_pr3 = 0, & 
     iwp2_dp1 = 0, & 
     iwp2_dp2 = 0, &
     iwp2_sdmp = 0, &
     iwp2_pd = 0, & 
     iwp2_cl = 0, &
     iwp2_sf = 0, &
     iwp2_splat = 0


  integer, public :: & 
     iwprtp_bt = 0,      & 
     iwprtp_ma = 0,      & 
     iwprtp_ta = 0,      & 
     iwprtp_tp = 0,      & 
     iwprtp_ac = 0,      & 
     iwprtp_bp = 0,      & 
     iwprtp_pr1 = 0,     & 
     iwprtp_pr2 = 0,     & 
     iwprtp_pr3 = 0,     & 
     iwprtp_dp1 = 0,     &
     iwprtp_mfl = 0,     &
     iwprtp_cl = 0,      & 
     iwprtp_sicl = 0,    & 
     iwprtp_pd = 0,      &
     iwprtp_forcing = 0, &
     iwprtp_mc = 0


  integer, public :: & 
     iwpthlp_bt = 0,      & 
     iwpthlp_ma = 0,      & 
     iwpthlp_ta = 0,      & 
     iwpthlp_tp = 0,      & 
     iwpthlp_ac = 0,      & 
     iwpthlp_bp = 0,      & 
     iwpthlp_pr1 = 0,     & 
     iwpthlp_pr2 = 0,     & 
     iwpthlp_pr3 = 0,     & 
     iwpthlp_dp1 = 0,     &
     iwpthlp_mfl = 0,     &
     iwpthlp_cl = 0,      & 
     iwpthlp_sicl = 0,    &
     iwpthlp_forcing = 0, &
     iwpthlp_mc = 0


  integer, public :: & 
     iupwp_bt = 0,  &
     iupwp_ma = 0,  &
     iupwp_ta = 0,  &
     iupwp_tp = 0,  &
     iupwp_ac = 0,  &
     iupwp_bp = 0,  &
     iupwp_pr1 = 0, &
     iupwp_pr2 = 0, &
     iupwp_pr3 = 0, &
     iupwp_pr4 = 0, &
     iupwp_dp1 = 0, &
     iupwp_mfl = 0, &
     iupwp_cl = 0


  integer, public :: & 
     ivpwp_bt = 0,  &
     ivpwp_ma = 0,  &
     ivpwp_ta = 0,  &
     ivpwp_tp = 0,  &
     ivpwp_ac = 0,  &
     ivpwp_bp = 0,  &
     ivpwp_pr1 = 0, &
     ivpwp_pr2 = 0, &
     ivpwp_pr3 = 0, &
     ivpwp_pr4 = 0, &
     ivpwp_dp1 = 0, &
     ivpwp_mfl = 0, &
     ivpwp_cl = 0


!    Dr. Golaz's new variance budget terms
!    qt was changed to rt to avoid confusion

  integer, public :: & 
     irtp2_bt = 0,      & 
     irtp2_ma = 0,      & 
     irtp2_ta = 0,      & 
     irtp2_tp = 0,      & 
     irtp2_dp1 = 0,     & 
     irtp2_dp2 = 0,     & 
     irtp2_pd = 0,      & 
     irtp2_cl = 0,      &
     irtp2_sf = 0,      &
     irtp2_forcing = 0, &
     irtp2_mc = 0
     

  integer, public :: & 
     ithlp2_bt = 0,      & 
     ithlp2_ma = 0,      & 
     ithlp2_ta = 0,      & 
     ithlp2_tp = 0,      & 
     ithlp2_dp1 = 0,     & 
     ithlp2_dp2 = 0,     & 
     ithlp2_pd = 0,      & 
     ithlp2_cl = 0,      &
     ithlp2_sf = 0,      &
     ithlp2_forcing = 0, &
     ithlp2_mc = 0


  integer, public :: & 
    irtpthlp_bt = 0,      & 
    irtpthlp_ma = 0,      & 
    irtpthlp_ta = 0,      & 
    irtpthlp_tp1 = 0,     & 
    irtpthlp_tp2 = 0,     & 
    irtpthlp_dp1 = 0,     & 
    irtpthlp_dp2 = 0,     & 
    irtpthlp_cl = 0,      &
    irtpthlp_sf = 0,      &
    irtpthlp_forcing = 0, &
    irtpthlp_mc = 0


  integer, public :: & 
    iup2 = 0, & 
    ivp2 = 0


  integer, public :: & 
    iup2_bt = 0, & 
    iup2_ta = 0, & 
    iup2_tp = 0, & 
    iup2_ma = 0, & 
    iup2_dp1 = 0, & 
    iup2_dp2 = 0, & 
    iup2_pr1 = 0, & 
    iup2_pr2 = 0, & 
    iup2_sdmp = 0, & 
    iup2_pd = 0, & 
    iup2_cl = 0, &
    iup2_sf = 0, &
    iup2_splat = 0, &
    ivp2_bt = 0, & 
    ivp2_ta = 0, & 
    ivp2_tp = 0, & 
    ivp2_ma = 0, & 
    ivp2_dp1 = 0, & 
    ivp2_dp2 = 0, & 
    ivp2_pr1 = 0, & 
    ivp2_pr2 = 0, & 
    ivp2_sdmp = 0, & 
    ivp2_pd = 0, & 
    ivp2_cl = 0, &
    ivp2_sf = 0, &
    ivp2_splat = 0


!       Passive scalars.  Note that floating point roundoff may make
!       mathematically equivalent variables different values.
  integer,target, allocatable, dimension(:), public :: & 
    isclrprtp,           & ! sclr'(1)rt'     / rt'^2
    isclrp2,             & ! sclr'(1)^2      / rt'^2
    isclrpthvp,          & ! sclr'(1)th_v'   / rt'th_v' 
    isclrpthlp,          & ! sclr'(1)th_l'   / rt'th_l' 
    isclrprcp,           & ! sclr'(1)rc'     / rt'rc'
    iwpsclrp,            & ! w'slcr'(1)      / w'rt'
    iwp2sclrp,           & ! w'^2 sclr'(1)   / w'^2 rt'
    iwpsclrp2,           & ! w'sclr'(1)^2    / w'rt'^2
    iwpsclrprtp,         & ! w'sclr'(1)rt'   / w'rt'^2
    iwpsclrpthlp           ! w'sclr'(1)th_l' / w'rt'th_l'


  integer, target, allocatable, dimension(:), public :: & 
     iwpedsclrp ! eddy sclr'(1)w'


  ! Indices for statistics in stats_rad_zt file
  integer, public :: &
    iT_in_K_rad = 0, &
    ircil_rad = 0, &
    io3l_rad = 0, &
    irsm_rad = 0, &
    ircm_in_cloud_rad = 0, &
    icloud_frac_rad = 0, & 
    iice_supersat_frac_rad = 0, &
    iradht_rad = 0, &
    iradht_LW_rad = 0, &
    iradht_SW_rad = 0, &
    ip_in_mb_rad = 0, &
    isp_humidity_rad = 0


  ! Indices for statistics in stats_rad_zm file
  integer, public :: &
    iFrad_LW_rad = 0, &
    iFrad_SW_rad = 0, &
    iFrad_SW_up_rad = 0, &
    iFrad_LW_up_rad = 0, &
    iFrad_SW_down_rad = 0, &
    iFrad_LW_down_rad = 0


  ! Indices for statistics in stats_sfc file

  integer, public :: & 
    iustar = 0, &
    isoil_heat_flux = 0,&
    iveg_T_in_K = 0,&
    isfc_soil_T_in_K = 0, &
    ideep_soil_T_in_K = 0,& 
    ilh = 0, & 
    ish = 0, & 
    icc = 0, & 
    ilwp = 0, &
    ivwp = 0, &        ! nielsenb
    iiwp = 0, &        ! nielsenb
    iswp = 0, &        ! nielsenb
    irwp = 0, &
    iz_cloud_base = 0, & 
    iz_inversion = 0, & 
    iprecip_rate_sfc = 0,    &    ! Brian
    irain_flux_sfc = 0,   &    ! Brian
    irrm_sfc = 0, & ! Brian
    iwpthlp_sfc = 0, &
    iprecip_frac_tol = 0

  integer, public :: &
    iwprtp_sfc = 0, &
    iupwp_sfc = 0, &
    ivpwp_sfc = 0, &
    ithlm_vert_avg = 0, &
    irtm_vert_avg = 0, &
    ium_vert_avg = 0, &
    ivm_vert_avg = 0, &
    iwp2_vert_avg = 0, & ! nielsenb
    iup2_vert_avg = 0, &
    ivp2_vert_avg = 0, &
    irtp2_vert_avg = 0, &
    ithlp2_vert_avg = 0, &
    iT_sfc         ! kcwhite

  integer, public :: & 
    iwp23_matrix_condt_num = 0, & 
    irtm_matrix_condt_num = 0, & 
    ithlm_matrix_condt_num = 0, & 
    irtp2_matrix_condt_num = 0, & 
    ithlp2_matrix_condt_num = 0, & 
    irtpthlp_matrix_condt_num = 0, & 
    iup2_vp2_matrix_condt_num = 0, & 
    iwindm_matrix_condt_num = 0

  integer, public :: & 
    imorr_snow_rate = 0


  integer, public :: &
    irtm_spur_src = 0,    &
    ithlm_spur_src = 0


  integer, public :: &
    iSkw_velocity = 0, & ! Skewness velocity
    iwp3_zm = 0, &
    ithlp3_zm = 0, &
    irtp3_zm = 0, &
    ia3_coef = 0, &
    ia3_coef_zt = 0

  integer, public :: &
    iwp3_on_wp2 = 0, &  ! w'^3 / w'^2 [m/s]
    iwp3_on_wp2_zt = 0  ! w'^3 / w'^2 [m/s]

  integer, public :: & 
    ilh_morr_snow_rate = 0

  integer, public :: & 
    ilh_vwp = 0, &
    ilh_lwp = 0, &
    ilh_sample_weights_sum = 0, &
    ilh_sample_weights_avg = 0


  integer, public :: &
    icloud_frac_refined = 0, &
    ircm_refined = 0

  integer, public :: &
    irtp2_from_chi = 0


  ! Variables that contains all the statistics

  type (stats), target, public :: stats_zt,   &    ! stats_zt grid
                                  stats_zm,   &    ! stats_zm grid
                                  stats_lh_zt,  &  ! stats_lh_zt grid
                                  stats_lh_sfc,  & ! stats_lh_sfc grid
                                  stats_rad_zt,  & ! stats_rad_zt grid
                                  stats_rad_zm,  & ! stats_rad_zm grid
                                  stats_sfc        ! stats_sfc


  ! Scratch space

  real( kind = core_rknd ), dimension(:), allocatable, public :: &
    ztscr01, ztscr02, ztscr03, & 
    ztscr04, ztscr05, ztscr06, & 
    ztscr07, ztscr08, ztscr09, & 
    ztscr10, ztscr11, ztscr12, & 
    ztscr13, ztscr14, ztscr15, & 
    ztscr16, ztscr17, ztscr18, &
    ztscr19, ztscr20, ztscr21


  real( kind = core_rknd ), dimension(:), allocatable, public :: &
    zmscr01, zmscr02, zmscr03, &
    zmscr04, zmscr05, zmscr06, & 
    zmscr07, zmscr08, zmscr09, & 
    zmscr10, zmscr11, zmscr12, & 
    zmscr13, zmscr14, zmscr15, &
    zmscr16, zmscr17


end module stats_variables
