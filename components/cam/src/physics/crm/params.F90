module params

  ! use grid, only: nzm
#ifdef CLUBB_CRM
  ! Use the CLUBB values for these constants for consistency
  use constants_clubb, only: Cp_clubb => Cp, grav_clubb => grav, Lv_clubb => Lv, Lf_clubb => Lf, &
  Ls_clubb => Ls, Rv_clubb => Rv, Rd_clubb => Rd, pi_clubb => pi
#else

#ifdef CRM
  use shr_const_mod, only: shr_const_rdair, shr_const_cpdair, shr_const_latvap, &
  shr_const_latice, shr_const_latsub, shr_const_rgas, &
  shr_const_mwwv, shr_const_stebol, shr_const_tkfrz, &
  shr_const_mwdair, shr_const_g, shr_const_karman, &
  shr_const_rhofw
#endif /*CRM*/

#endif

  implicit none

#ifdef CRM_SINGLE_PRECISION
  integer, parameter :: crm_rknd = selected_real_kind( 6) ! 4 byte real
#else
  ! whannah - default precision of real - kind(1.d0)
  integer, parameter :: crm_rknd = selected_real_kind(12) ! 8 byte real
#endif

  !   Constants:

#ifdef CLUBB_CRM
  ! Define Cp, ggr, etc. in module constants_clubb
  real(crm_rknd), parameter :: cp    = real( Cp_clubb   ,crm_rknd)
  real(crm_rknd), parameter :: ggr   = real( grav_clubb ,crm_rknd)
  real(crm_rknd), parameter :: lcond = real( Lv_clubb   ,crm_rknd)
  real(crm_rknd), parameter :: lfus  = real( Lf_clubb   ,crm_rknd)
  real(crm_rknd), parameter :: lsub  = real( Ls_clubb   ,crm_rknd)
  real(crm_rknd), parameter :: rv    = real( Rv_clubb   ,crm_rknd)
  real(crm_rknd), parameter :: rgas  = real( Rd_clubb   ,crm_rknd)
#else
#ifndef CRM
  real(crm_rknd), parameter :: cp    = 1004.          ! Specific heat of air, J/kg/K
  real(crm_rknd), parameter :: ggr   = 9.81           ! Gravity acceleration, m/s2
  real(crm_rknd), parameter :: lcond = 2.5104e+06     ! Latent heat of condensation, J/kg
  real(crm_rknd), parameter :: lfus  = 0.3336e+06   ! Latent heat of fusion, J/kg
  real(crm_rknd), parameter :: lsub  = 2.8440e+06     ! Latent heat of sublimation, J/kg
  real(crm_rknd), parameter :: rv    = 461.           ! Gas constant for water vapor, J/kg/K
  real(crm_rknd), parameter :: rgas  = 287.           ! Gas constant for dry air, J/kg/K
#else
  real(crm_rknd), parameter :: cp    = real( shr_const_cpdair ,crm_rknd)
  real(crm_rknd), parameter :: ggr   = real( shr_const_g      ,crm_rknd)
  real(crm_rknd), parameter :: lcond = real( shr_const_latvap ,crm_rknd)
  real(crm_rknd), parameter :: lfus  = real( shr_const_latice ,crm_rknd)
  real(crm_rknd), parameter :: lsub  = real( lcond + lfus     ,crm_rknd)
  real(crm_rknd), parameter :: rgas  = real( shr_const_rdair  ,crm_rknd)
  real(crm_rknd), parameter :: rv    = real( shr_const_rgas/shr_const_mwwv ,crm_rknd)
#endif
#endif
  real(crm_rknd), parameter :: diffelq = 2.21e-05     ! Diffusivity of water vapor, m2/s
  real(crm_rknd), parameter :: therco = 2.40e-02      ! Thermal conductivity of air, J/m/s/K
  real(crm_rknd), parameter :: muelq = 1.717e-05      ! Dynamic viscosity of air

  real(crm_rknd), parameter :: fac_cond = lcond/cp
  real(crm_rknd), parameter :: fac_fus  = lfus/cp
  real(crm_rknd), parameter :: fac_sub  = lsub/cp

#ifdef CLUBB_CRM
  real(crm_rknd), parameter ::  pi =  real( pi_clubb ,crm_rknd)
#else
  real(crm_rknd), parameter ::  pi = 3.141592653589793
#endif



  !
  ! internally set parameters:

  real(crm_rknd)   epsv     ! = (1-eps)/eps, where eps= Rv/Ra, or =0. if dosmoke=.true.
  logical:: dosubsidence = .false.
  real(crm_rknd), allocatable :: fcorz(:)      ! Vertical Coriolis parameter

  !----------------------------------------------
  ! Parameters set by PARAMETERS namelist:
  ! Initialized to default values.
  !----------------------------------------------

  real(crm_rknd):: ug = 0.        ! Velocity of the Domain's drift in x direction
  real(crm_rknd):: vg = 0.        ! Velocity of the Domain's drift in y direction
  real(crm_rknd), allocatable :: fcor(:)  ! Coriolis parameter
  real(crm_rknd), allocatable :: longitude0(:)    ! latitude of the domain's center
  real(crm_rknd), allocatable :: latitude0 (:)    ! longitude of the domain's center

  real(crm_rknd)::   z0     =0.035  ! roughness length
  real(crm_rknd)::   soil_wetness =1.! wetness coeff for soil (from 0 to 1.)
  integer:: ocean_type =0 ! type of SST forcing
  logical:: cem =.false.    ! flag for Cloud Ensemble Model
  logical:: les =.false.    ! flag for Large-Eddy Simulation
  logical:: ocean =.false.  ! flag indicating that surface is water
  logical:: land =.false.   ! flag indicating that surface is land
  logical:: sfc_flx_fxd =.false. ! surface sensible flux is fixed
  logical:: sfc_tau_fxd =.false.! surface drag is fixed

  real(crm_rknd):: timelargescale =0. ! time to start large-scale forcing

  ! nudging boundaries (between z1 and z2, where z2 > z1):
  real(crm_rknd):: nudging_uv_z1 =-1., nudging_uv_z2 = 1000000.
  real(crm_rknd):: nudging_t_z1 =-1., nudging_t_z2 = 1000000.
  real(crm_rknd):: nudging_q_z1 =-1., nudging_q_z2 = 1000000.
  real(crm_rknd):: tauls = 99999999.    ! nudging-to-large-scaler-profile time-scale
  real(crm_rknd):: tautqls = 99999999.! nudging-to-large-scaler-profile time-scale for scalars

  logical:: dodamping = .false.
  logical:: doupperbound = .false.
  logical:: docloud = .false.
  logical:: doclubb = .false. ! Enabled the CLUBB parameterization (interactively)
  logical:: doclubb_sfc_fluxes = .false. ! Apply the surface fluxes within the CLUBB code rather than SAM
  logical:: doclubbnoninter = .false. ! Enable the CLUBB parameterization (non-interactively)
  logical:: docam_sfc_fluxes = .false.   ! Apply the surface fluxes within CAM
  logical:: doprecip = .false.
  logical:: dolongwave = .false.
  logical:: doshortwave = .false.
  logical:: dosgs = .false.
  logical:: docoriolis = .false.
  logical:: docoriolisz = .false.
  logical:: dofplane = .true.
  logical:: dosurface = .false.
  logical:: dolargescale = .false.
  logical:: doradforcing = .false.
  logical:: dosfcforcing = .false.
  logical:: doradsimple = .false.
  logical:: donudging_uv = .false.
  logical:: donudging_tq = .false.
  logical:: donudging_t = .false.
  logical:: donudging_q = .false.
  logical:: doensemble = .false.
  logical:: dowallx = .false.
  logical:: dowally = .false.
  logical:: docolumn = .false.
  logical:: docup = .false.
  logical:: doperpetual = .false.
  logical:: doseasons = .false.
  logical:: doradhomo = .false.
  logical:: dosfchomo = .false.
  logical:: dossthomo = .false.
  logical:: dodynamicocean = .false.
  logical:: dosolarconstant = .false.
  logical:: dotracers = .false.
  logical:: dosmoke = .false.
  logical:: notracegases = .false.

  ! Specify solar constant and zenith angle for perpetual insolation.
  ! Based onn Tompkins and Graig (1998)
  ! Note that if doperpetual=.true. and dosolarconstant=.false.
  ! the insolation will be set to the daily-averaged value on day0.
  real(crm_rknd):: solar_constant = 685. ! solar constant (in W/m2)
  real(crm_rknd):: zenith_angle = 51.7   ! zenith angle (in degrees)

  integer:: nensemble =0   ! the number of subensemble set of perturbations
  integer:: perturb_type  = 0 ! type of initial noise in setperturb()
  integer:: nclubb = 1 ! SAM timesteps per CLUBB timestep
  ! Initial bubble parameters. Activated when perturb_type = 2
  real(crm_rknd):: bubble_x0 = 0.
  real(crm_rknd):: bubble_y0 = 0.
  real(crm_rknd):: bubble_z0 = 0.
  real(crm_rknd):: bubble_radius_hor = 0.
  real(crm_rknd):: bubble_radius_ver = 0.
  real(crm_rknd):: bubble_dtemp = 0.
  real(crm_rknd):: bubble_dq = 0.

  real(crm_rknd) uhl      ! current large-scale velocity in x near sfc
  real(crm_rknd) vhl      ! current large-scale velocity in y near sfc
  real(crm_rknd) taux0    ! surface stress in x, m2/s2
  real(crm_rknd) tauy0    ! surface stress in y, m2/s2


contains

  
  subroutine allocate_params(ncrms)
    implicit none
    integer, intent(in) :: ncrms
    allocate(fcor (ncrms))
    allocate(fcorz(ncrms))
    allocate(longitude0(ncrms))
    allocate(latitude0 (ncrms))

    fcor  = 0
    fcorz = 0
    longitude0 = 0
    latitude0  = 0
  end subroutine allocate_params

  
  subroutine deallocate_params()
    implicit none
    deallocate(fcor )
    deallocate(fcorz)
    deallocate(longitude0)
    deallocate(latitude0 )
  end subroutine deallocate_params


end module params
