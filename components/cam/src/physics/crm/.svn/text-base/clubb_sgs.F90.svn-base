!-------------------------------------------------------------------------------
! $Id: clubb_sgs.F90 831 2012-02-14 19:34:37Z dschanen@uwm.edu $
module clubb_sgs
#ifdef CLUBB_CRM
! Description:
! Contains function and subroutines for interfacing with the UW Milwaukee
! Single-Column Model.
!-------------------------------------------------------------------------------

  use clubb_core, only: &
    setup_clubb_core, advance_clubb_core, &
    cleanup_clubb_core

  use clubb_precision, only: &
    time_precision, & ! Constant(s)
    core_rknd

  use domain, only: &
    nsubdomains_x, &
    nsubdomains_y

  implicit none

  private

  public :: clubb_sgs_setup, advance_clubb_sgs, clubb_sgs_cleanup, &
    apply_clubb_sgs_tndcy, apply_clubb_sgs_tndcy_scalar, apply_clubb_sgs_tndcy_mom, t2thetal

  logical, private :: lstats_clubb

  integer, dimension(nsubdomains_x*nsubdomains_y), private :: &
    sample_nodes, x_samp, y_samp

  integer, private :: x_samp_node, y_samp_node

  contains
!-------------------------------------------------------------------------------
  subroutine clubb_sgs_setup( dt_clubb, latitude, longitude, z, rho, zi, rhow, tv0, tke )

! Description:
!   Initialize UWM CLUBB.

! References:
!   None
!-------------------------------------------------------------------------------

    ! From the CLUBB directory
    use error_code, only: &
      clubb_no_error, set_clubb_debug_level ! Subroutines

    use parameter_indices, only: &
      nparams ! Constant

    use constants_clubb, only: &
      em_min, w_tol_sqd, rt_tol, thl_tol, zero_threshold, & ! Constants
      fstderr, fstdout


    use grid_class, only: &
      zm2zt, zt2zm, & ! Functions
      gr ! Derived type

    ! These are only needed if we're using a passive scalar
    use array_index, only: &
      iisclr_rt, iisclr_thl, iisclr_CO2, &    ! [kg/kg]/[K]/[1e6 mol/mol]
      iiedsclr_rt, iiedsclr_thl, iiedsclr_CO2 ! "    "

    use parameters_tunable, only: &
      read_parameters ! Subroutine

    use stats_subs, only: &
      stats_init ! Subroutine

    use model_flags, only: &
      l_use_boussinesq, & ! Variables
      l_tke_aniso

    ! From the SAM directory
    use grid, only: rank, nx, ny, nz, nzm, dx, dy, time, case, caseid, &
      nrestart, dimx1_s, dimx2_s, dimy1_s, dimy2_s, ntracers, doclubb_sfc_fluxes ! Variable(s)

    use params, only: lcond, cp ! Constants

    use clubbvars, only: &
      upwp,   &! u'w'.                 [m^2/s^2]
      vpwp,   &! u'w'.                 [m^2/s^2]
      up2,    &! u'^2                  [m^2/s^2]
      vp2,    &! v'^2                  [m^2/s^2]
      wprtp,  &! w' r_t'.              [(m kg)/(s kg)]
      wpthlp, &! w' th_l'.             [(m K)/s]
      wprcp,  &! w' r_c'               [(kg/kg) m/s]
      wp2,    &! w'^2.                 [m^2/s^2]
      rtp2,   &! r_t'^2.               [(kg/kg)^2]
      thlp2,  &! th_l'^2.              [K^2]
      rtpthlp,&! r_t' th_l'.           [(kg K)/kg]
      wp3      ! w'^3.                 [m^3/s^3]

    use clubbvars, only: &
      tracer_tndcy, & ! Time tendency of the SAM set of tracers
      t_tndcy,  & ! CLUBB contribution to moist static energy  [K/s]
      qc_tndcy, & ! CLUBB contribution to liquid water         [kg/kg/s]
      qv_tndcy, & ! CLUBB contribution to vapor water          [kg/kg/s]
      u_tndcy,  & ! CLUBB contribution to x-wind               [m/s^2]
      v_tndcy     ! CLUBB contribution to y-wind               [m/s^2]

    use clubbvars, only: &
      sclrp2,      & ! Passive scalar variance.       [{units vary}^2]
      sclrpthlp,   & ! Passive scalar covariance.     [{units vary}^2]
      sclrprtp,    & ! Passive scalar covariance.     [{units vary}^2]
      wpsclrp        ! w'sclr'                        [units vary m/s]

    use clubbvars, only: &
      rho_ds_zm,       & ! Dry, static density on momentum levels      [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermodynamic levels [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density on momentum levels [m^3/kg]
      invrs_rho_ds_zt, & ! Inv. dry, static density on thermo. levels  [m^3/kg]
      thv_ds_zm,       & ! Dry, base-state theta_v on momentum levels  [K]
      thv_ds_zt          ! Dry, base-state theta_v on thermo. levels   [K]

    use clubbvars, only: &
      sclr_tol,   & ! Tolerance on high-order scalars
      edsclr_dim, & ! Number of eddy-diffusivity scalars
      sclr_dim      ! Numer of high-order scalars

    use clubbvars, only: &
      tndcy_precision ! Precision of CLUBB's contribution to the tendencies of mean variables

#ifdef CRM
    use clubbvars, only: lrestart_clubb
#endif
              

    implicit none

    ! Constant parameters
    logical, parameter :: &
      l_uv_nudge       = .false.,  & ! Use u/v nudging (not used)
      l_implemented    = .true.      ! Implemented in a host model (always true)

    integer, parameter :: &
      grid_type    = 2, &  ! The 2 option specifies stretched thermodynamic levels
      hydromet_dim = 0, &  ! The hydromet array in SAM-CLUBB is currently 0 elements
      iunit = 50           ! Fortran I/O unit

    character(len=6), parameter :: &
      saturation_equation = "flatau" ! Flatau polynomial approximation for SVP

    real(kind=core_rknd), parameter :: &
      theta0   = 300._core_rknd, &! Reference temperature                     [K]
      ts_nudge = 86400._time_precision ! Time scale for u/v nudging (not used)     [s]

    ! Input Variables
    real(kind=time_precision), intent(in) :: &
      dt_clubb ! SAM-CLUBB subcycled model timestep   [s]

    real, dimension(nx, ny), intent(in) :: &
      latitude,  & ! Latitudes for SAM's dynamical core  [degrees_N]
      longitude    ! Longitudes for SAM's dynamical core [degrees_E]

    real, dimension(nzm), intent(in) :: &
      z,   & ! Thermodynamic/Scalar grid in SAM         [m]
      rho    ! Thermodynamic/Scalar density in SAM      [kg/m^3]

    real, dimension(nz), intent(in) :: &
      zi, & ! Momentum/Vertical Velocity grid in SAM    [m]
      rhow  ! Momentum/Vertical Velocity density in SAM [kg/m^3]

    real, dimension(nzm), intent(in) :: &
      tv0 ! Virtual potential temperature from SAM      [K]

    real, dimension(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm), intent(in) :: &
      tke   ! SGS TKE                   [m^2/s]

    ! Local Variables
    real(kind=core_rknd), dimension(nparams) :: &
      clubb_params ! These adjustable CLUBB parameters (C1, C2 ...)

    ! 1D variables with ghost points at the lowest level
    real(kind=core_rknd), dimension(nz) :: &
      zt, & ! Thermodynamic grid        [m]
      zm, & ! Momentum grid             [m]
      em    ! Turbulent kinetic energy  [-]

    logical :: l_stats  ! Stats enabled (T/F)

    real(kind=time_precision) :: &
      stats_tsamp, & ! Sampling interval for a single column of CLUBB data  [s]
      stats_tout     ! Output interval for a single column of CLUBB data    [s]

    character(len=10)  :: stats_fmt     ! Format of stats output (netCDF/GrADS)
    character(len=250) :: fname_prefix  ! Prefix for stats filename

    ! Horizontal grid spacings (i.e., dx and dy), used for computing Lscale_max
    real(kind=core_rknd) :: host_dx, host_dy ! [m]

    real(kind=core_rknd), dimension(1) :: &
      rlat, rlon ! Latitude and Longitude for stats [degrees]

    integer :: err_code ! Code for when CLUBB fails

    integer :: i, j, ig, jg ! Indices

    integer :: ilen  ! Length of a string

    logical :: l_host_applies_sfc_fluxes ! Whether the host model applies the surface fluxes

    namelist /stats_setting/ l_stats, stats_fmt, stats_tsamp, stats_tout, &
      sample_nodes, x_samp, y_samp

!-------------------------------------------------------------------------------
!  SAM uses an Arakawa C type grid for the 3D quantities.  The UWM SCM has an
!  additional `ghost' point on the lowest pressure/thermodynamic level.
!  i.e.
!
!  SAM vert. vel. grid          UWM SCM moment. grid
!
!      Dimension    Elevation         Dimension       Elevation 
!  . . . (nz  ) . . zi(nz  )    . . . (gr%nzmax  ) . . gr%zm(gr%nzmax  ) . . .
!  . . . (nz-1) . . zi(nz-1)    . . . (gr%nzmax-1) . . gr%zm(gr%nzmax-1) . . .
!           |          |                   |                    |
!  . . . (1   ) . . zi(1   )    . . . (1        ) . . gr%zm(1      )   . . .
!
!  In SAM the lowest grid point on the vertical velocity levels (or `interface' 
!  levels) is always 0 meters.  The UWM SCM supports an arbitrary starting
!  point for the momentum grid, but this code assumes 0 meters.
!
!  SAM pressure grid            UWM SCM thermo. grid
!
!      Dimension    Elevation         Dimension       Elevation 
!  . . . (nz-1) . . z(nz-1)     . . . (gr%nzmax  ) . . gr%zt(gr%nzmax  ) . . .
!  . . . (nz-2) . . z(nz-2)     . . . (gr%nzmax-1) . . gr%zt(gr%nzmax-1) . . .
!           |          |                   |                    |
!  . . . (1   ) . . z(1   )     . . . (2        ) . . gr%zt(2        ) . . .
!  / / /  N/A   / /   N/A       / / / (1        ) / / gr%zt(1        ) / / /
!
!  Note that the lowest SCM point is below ground.
!-------------------------------------------------------------------------------

    !----- Begin Code -----

    ! Set the ghost point to be the distance between the first interface level,
    ! which is always zero, and the first pressure level.
    zt(1)    = real( -z(1), kind=core_rknd ) ! [m]
    ! All other pressure levels are determined by the host model
    zt(2:nz) = real( z(1:nzm), kind=core_rknd ) ! [m]

    zm = real( zi, kind=core_rknd )

    ! Set the SCM parameters (C2, etc. ) based on default values
    !call read_parameters( -99, "", clubb_params )

    ! Set the SCM parameters (C2, etc. ) based on a namelist
#ifdef CRM
    ! Set the SCM parameters (C2, etc. ) based on default values 
    call read_parameters( -99, "", clubb_params )
#else
    ! Set the SCM parameters (C2, etc. ) based on a namelist
    call read_parameters( iunit, "CLUBB_PARAMETERS/tunable_parameters.in", clubb_params )
#endif

    ! Set the debug level.  Level 2 has additional computational expense since
    ! it checks the array variables in CLUBB for invalid values.
    call set_clubb_debug_level( 0 )

    host_dx = real( dx, kind=core_rknd )
    host_dy = real( dy, kind=core_rknd )

    ! These are for emulating total water or thetal for testing purposes
    iisclr_rt  = -1
    iisclr_thl = -1
    iisclr_CO2 = -1

    iiedsclr_rt  = -1
    iiedsclr_thl = -1
    iiedsclr_CO2 = -1

    ! Sanity check
    if ( sclr_dim > 0 .and. edsclr_dim > 0 ) then
      write(fstderr,*) "Only one scalar scheme can be enabled at one time"
      call task_abort()
    end if

    ! This is the tolerance on total water in the CLUBB SCM
    ! Other tracers will need this value set according to their order of 
    ! magnitude and the units they are in. Keep in mind that the variable
    ! sclrp2 will be clipped to a minimum value of sclr_tol^2
    sclr_tol(1:sclr_dim) = 1.e-8_core_rknd ! total water is in kg/kg

    ! Determine whether clubb is applying the surface flux or the host model 
    ! from the namelist variable doclubb_sfc_fluxes
    l_host_applies_sfc_fluxes = .not. doclubb_sfc_fluxes

    call setup_clubb_core     &
         ( nz, theta0, ts_nudge, & ! In
           hydromet_dim,  sclr_dim, &  !  In
           sclr_tol, edsclr_dim, clubb_params, & ! In
           l_host_applies_sfc_fluxes, & ! In
           l_uv_nudge, saturation_equation,  & ! In
           l_implemented, grid_type, zm(2), zm(1), zm(nz), & ! In
           zm(1:nz), zt(1:nz), & ! In
           host_dx, host_dy, zm(1), & ! In
           err_code )

    if ( err_code /= CLUBB_no_error ) then
      write(fstderr,*) "Initialization of CLUBB failed"
      call task_abort()
    end if

    ! Initialize stats_setting
    l_stats = .false.
    stats_fmt = "grads"
    stats_tsamp = 60._time_precision
    stats_tout = 60._time_precision
    sample_nodes(:) = -1 ! Which nodes are outputting a column
    x_samp(:) = -1       ! Which x point for the nth node
    y_samp(:) = -1       ! Which y point for the nth node

    ! Figure out which node and points we're sampling
#ifdef CRM
      l_stats      = .false.
      stats_fmt   = 'grads'
      stats_tsamp = 60.0
      stats_tout  = 60.0
      sample_nodes(:) = -1
      x_samp(:)     = 1
      y_samp(:)     = 1
      lstats_clubb = .false.
#else
    open(unit=10, file="clubb_stats")
    read(unit=10, nml=stats_setting)
    close(unit=10)

    if ( is_a_sample_node( rank ) .and. l_stats ) then

      ! Determine and save the local x and y to write to be written to disk
      call get_sample_points( rank, x_samp_node, y_samp_node )

      ! Figure out the position on the global grid
      call task_rank_to_index( rank, ig, jg )

      ! The filename follows the following format:
      ! case_caseid_x<point on global grid>_y<point on global grid>_<grid>
      ! e.g. (variables in single quotes)
      ! 'BOMEX'_'64x64x75_scm_LES'_x000'1'_y00'10'_'zt'
      fname_prefix = trim( case )//"_"//trim( caseid )
      ilen = len( trim( fname_prefix ) )
      fname_prefix = trim( fname_prefix )//"_x0000_y0000"
      write(unit=fname_prefix(ilen+3:ilen+6),fmt='(i4.4)') ig+x_samp_node
      write(unit=fname_prefix(ilen+9:ilen+12),fmt='(i4.4)') jg+y_samp_node
      rlat = real( latitude(x_samp_node,y_samp_node), kind=core_rknd )
      rlon = real( longitude(x_samp_node,y_samp_node), kind=core_rknd )

      ! Use a bogus date, since SAM does not track the year, and it would require
      ! some work to covert the `day' variable to MMDD format
      call stats_init( iunit, fname_prefix, "./OUT_STAT/", l_stats, &
                       stats_fmt, stats_tsamp, stats_tout, "clubb_stats", &
                       nz, zt, zm, nz, zt, nz, zm, 1, 4, 1900, &
                       rlat, rlon, &
                       time, dt_clubb )

      ! If CLUBB stats are on for this node, toggle a flag in this module
      write(fstdout,*) "CLUBB stats enabled"
      lstats_clubb = .true.
    else
      lstats_clubb = .false.
      x_samp_node = -1
      y_samp_node = -1
    end if
#endif /*CRM*/

#ifdef CRM
!+++mhwang, 2012-02-06 (Minghuai.Wang@pnnl.gov)
!    rho_ds_zm, rho_ds_zt, thv_ds_zt, thv_ds_zm, invrs_rho_ds_zm, invrs_rho_ds_zt are needed
!    to be copied from those from the GCM at the beginning of each GCM time step. 
     if (lrestart_clubb) then
     ! Set variables for the use of the anelastic equation set in CLUBB.
     ! Set the value of dry, static, base-state density.
      rho_ds_zm(:) = rhow(:)
      rho_ds_zt(2:nz) = rho(1:nzm)
      rho_ds_zt(1) = LIN_EXT( rho_ds_zt(3), rho_ds_zt(2), gr%zt(3), gr%zt(2), gr%zt(1) )
     ! Set the value of dry, base-state theta_v.
      thv_ds_zt(2:nz) = tv0(1:nzm)
      thv_ds_zt(1) = tv0(1)
      thv_ds_zm(:) = zt2zm( thv_ds_zt )

      ! Set the value of inverse dry, static, base-state density based on the
      ! value of dry, static, base-state density.
      invrs_rho_ds_zm(:) = 1.0 / rho_ds_zm(:)
      invrs_rho_ds_zt(:) = 1.0 / rho_ds_zt(:)
     end if
#endif /*CRM*/

    ! If this is restart run, just return at this point and do not re-initialize
    !  any variables as we would a run starting from the beginning.

#ifndef CRM
    if ( nrestart /= 0 ) return
#else
   if  (lrestart_clubb ) return
#endif

    if ( sclr_dim > 0 ) then
      sclrp2    = 0._core_rknd
      sclrprtp  = 0._core_rknd
      sclrpthlp = 0._core_rknd
      wpsclrp   = 0._core_rknd
    end if

    ! Initialize CLUBB's tendencies to 0
    t_tndcy  = 0._tndcy_precision
    qc_tndcy = 0._tndcy_precision
    qv_tndcy = 0._tndcy_precision
    u_tndcy  = 0._tndcy_precision
    v_tndcy  = 0._tndcy_precision

    if ( ntracers > 0 ) then
      tracer_tndcy = 0._tndcy_precision
    end if

    ! SAM's dynamical core is anelastic, so l_use_boussineq should probably be
    ! set to false generally, as it is by default in the CLUBB SCM.
    if ( l_use_boussinesq ) then
      rho_ds_zm(:) = 1._core_rknd
      rho_ds_zt(:) = 1._core_rknd
      ! Set the value of dry, base-state theta_v.
      thv_ds_zm(:) = theta0
      thv_ds_zt(:) = theta0
    else 
      ! Set variables for the use of the anelastic equation set in CLUBB.
      ! Set the value of dry, static, base-state density.
      rho_ds_zm(:) = real( rhow(:), kind=core_rknd )
      rho_ds_zt(2:nz) = real( rho(1:nzm), kind=core_rknd )
      rho_ds_zt(1) = LIN_EXT( rho_ds_zt(3), rho_ds_zt(2), gr%zt(3), gr%zt(2), gr%zt(1) )
      ! Set the value of dry, base-state theta_v.
      thv_ds_zt(2:nz) = real( tv0(1:nzm), kind=core_rknd )
      thv_ds_zt(1) = real( tv0(1), kind=core_rknd )
      thv_ds_zm(:) = zt2zm( thv_ds_zt )
    end if
    ! Set the value of inverse dry, static, base-state density based on the
    ! value of dry, static, base-state density.
    invrs_rho_ds_zm(:) = 1.0_core_rknd / rho_ds_zm(:)
    invrs_rho_ds_zt(:) = 1.0_core_rknd / rho_ds_zt(:)

    ! Determine the initial value of some variables as in WRF-HOC

    wprtp(:,:,:)       = 0._core_rknd ! w'rt'
    wpthlp(:,:,:)      = 0._core_rknd ! w'thl'
    wprcp(:,:,:)       = 0._core_rknd ! w'rc'
    wp3(:,:,:)         = 0._core_rknd ! w'^3
    wp2(:,:,:)         = w_tol_sqd    ! w'^2
    up2(:,:,:)         = w_tol_sqd    ! u'^2
    vp2(:,:,:)         = w_tol_sqd    ! v'^2
    rtp2(:,:,:)        = rt_tol**2    ! rt'^2
    thlp2(:,:,:)       = thl_tol**2   ! thl'^2
    rtpthlp(:,:,:)     = 0._core_rknd ! rt'thl'
    upwp(:,:,:)        = 0._core_rknd ! u'w'
    vpwp(:,:,:)        = 0._core_rknd ! v'w'

    do i=1, nx, 1
      do j=1, ny, 1

        ! Extrapolate intial SGS TKE and use it to compute wp2
        ! This value is going to depend on initial noise and whether 
        ! Smagorinksy diffusion is enabled
        em(2:nz) = real( tke(i,j,1:nzm), kind=core_rknd )
        em(1)    = LIN_EXT( em(3), em(2), gr%zt(3), gr%zt(2), gr%zt(1) )
        em(1:nz) = max( zt2zm( em(1:nz) ), em_min )

        !!!! Initialize w'^2 based on initial SGS TKE !!!!

        if ( l_tke_aniso ) then

          ! SGS TKE:  em = (1/2) * ( w'^2 + u'^2 + v'^2 )
          ! Evenly divide SGS TKE into its component
          ! contributions (w'^2, u'^2, and v'^2).

          wp2(i,j,1:nz) = (2._core_rknd/3._core_rknd) * em(1:nz)
          up2(i,j,1:nz) = (2._core_rknd/3._core_rknd) * em(1:nz)
          vp2(i,j,1:nz) = (2._core_rknd/3._core_rknd) * em(1:nz)

        else

          ! Assume isotropy for initialization of wp2
          ! SGS TKE:  em = (3/2) * w'^2

          wp2(i,j,1:nz) = (2._core_rknd/3._core_rknd) * em(1:nz)

        end if

      end do ! j=1..ny
    end do ! i=1..nx

    return
  end subroutine clubb_sgs_setup

!-------------------------------------------------------------------------------
  subroutine advance_clubb_sgs( dt_clubb, time_initial, time_current,         &
                                rho, rhow, wsub, u, v, w, qpl, qci, qpi, &
                                t, qv, qcl )

! Description:
!   Advance Cloud Layer Unified By Binormals one timestep.

! References:
!   ``A PDF-Based Model for Boundary Layer Clouds. Part I:
!     Method and Model Description'' Golaz, et al. (2002)
!   JAS, Vol. 59, pp. 3540--3551.
!-------------------------------------------------------------------------------

  ! From SAM
    use grid, only: &
      nx, ny, nxp1, nyp1, nz, nzm,&! Local grid dimensions
      nx_gl, ny_gl,   &! Global grid dimensions
      dimx1_s, dimx2_s, dimy1_s, dimy2_s,& ! Scalars dimensions
      dimx1_u, dimx2_u, dimy1_u, dimy2_u,& ! U wind dimensions
      dimx1_v, dimx2_v, dimy1_v, dimy2_v,& ! V wind dimensions
      dimx1_w, dimx2_w, dimy1_w, dimy2_w,& ! W wind dimensions
      YES3D, rank, pres, doclubb, dompi, & 
      ntracers 

    use params, only: cp, lfus, lsub, &
      ug, vg ! ug and vg are scalars, not arrays

    use vars, only: fcory, fluxbt, fluxbq, fluxbu, fluxbv, gamaz, prespot

    use microphysics, only: nmicro_fields

    use clubbvars, only: &
      upwp,        &! u'w'.                 [m^2/s^2]
      vpwp,        &! u'w'.                 [m^2/s^2]
      up2,         &! u'^2                  [m^2/s^2]
      vp2,         &! v'^2                  [m^2/s^2]
      wprtp,       &! w' r_t'.              [(m kg)/(s kg)]
      wpthlp,      &! w' th_l'.             [(m K)/s]
      wprcp,       &! w' r_c'.              [(kg/kg) m/s]
      wp2,         &! w'^2.                 [m^2/s^2]
      rtp2,        &! r_t'^2.               [(kg/kg)^2]
      thlp2,       &! th_l'^2.              [K^2]
      rtpthlp,     &! r_t' th_l'.           [(kg K)/kg]
      rcm,         &! Cloud water           [kg/kg]
      cloud_frac,  &! Cloud Fraction.       [-]
      rcm_in_layer,&! rcm in cloud layer    [kg/kg]
      cloud_cover, &! Cloud Cover           [-]
      wp3,         &! w'^3.                 [m^3/s^3]
      um,          &! x-wind                [m/s]
      vm            ! y-wind                [m/s]

    use clubbvars, only: &
      sclrp2,      & ! Passive scalar variance.       [{units vary}^2]
      sclrpthlp,   & ! Passive scalar covariance.     [{units vary}^2]
      sclrprtp,    & ! Passive scalar covariance.     [{units vary}^2]
      wpsclrp        ! w'sclr'                        [units vary m/s]

    use clubbvars, only: &
      u_tndcy,&  ! CLUBB contribution to the x wind
      v_tndcy,&  ! CLUBB contribution to the y wind
      qv_tndcy,& ! CLUBB contribution to vapor water mixing ratio
      qc_tndcy,& ! CLUBB contribution to liquid water mixing ratio
      t_tndcy    ! CLUBB contribution to moist static energy

    use clubbvars, only: &
      tracer_tndcy ! CLUBB contribution to a set of tracers

    use clubbvars, only: &
      rho_ds_zm,       & ! Dry, static density on momentum levels      [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermodynamic levels [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density on momentum levels [m^3/kg]
      invrs_rho_ds_zt, & ! Inv. dry, static density on thermo. levels  [m^3/kg]
      thv_ds_zm,       & ! Dry, base-state theta_v on momentum levels  [K]
      thv_ds_zt          ! Dry, base-state theta_v on thermo. levels   [K]

    use clubbvars, only: &
      sclr_dim, & ! Constant(s)
      edsclr_dim

    use clubbvars, only: &
      tndcy_precision ! Constant(s)

#ifndef CRM
    use tracers, only: &
#else
    use crmtracers, only: &
#endif
      fluxbtr, &
      tracer

    ! From CLUBB
    use error_code, only: &
      clubb_no_error, &    ! Constant
      clubb_at_least_debug_level ! Function

    use grid_class, only: &
      zm2zt, zt2zm, & ! Functions
      gr ! Derived type

    use stats_variables, only: &
      l_stats, l_stats_samp ! Logicals

    use stats_subs, only: &
      stats_begin_timestep, stats_end_timestep ! Subroutines

    use pdf_parameter_module, only: &
      pdf_parameter

    use constants_clubb, only: &
      fstderr ! Constant

    implicit none

    ! Parameters
    logical, parameter :: & 
      l_implemented = .true., & ! CLUBB is implemented in a host model, so this is true
      l_advect      = .false. ! Whether to advect around the high-order moments

    real(kind=core_rknd), parameter, dimension(nz) :: &
      zero = 0.0_core_rknd ! Field of zeros

    ! Input
    real(kind=time_precision), intent(in) :: &
      dt_clubb ! Timestep size for CLUBB [s]

    real(kind=time_precision), intent(in) :: &
      time_initial, time_current ! Initial and current time [s]

    real, intent(in), dimension(nzm) :: &
      rho  ! Air density        [kg/m^3]

    real, intent(in), dimension(nz) :: &
      wsub,&! Imposed vertical velocity       [m/s]
      rhow  ! Density on vert velocity grid   [kg/m^3]

    real, intent(in), dimension(dimx1_u:dimx2_u,dimy1_u:dimy2_u,nzm) :: &
      u  ! u wind    [m/s]

    real, intent(in), dimension(dimx1_v:dimx2_v,dimy1_v:dimy2_v,nzm) :: &
      v  ! v wind    [m/s]

    real, intent(in), dimension(dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz ) :: &
      w ! Vertical wind  [m/s]

    real, intent(in), dimension(nx,ny,nzm) :: &
      qpl,& ! Liquid water mixing ratio (precipitation) [kg/kg]
      qci,& ! Cloud ice water mixing ratio              [kg/kg]
      qpi   ! Snow + graupel mixing ratio (precip)      [kg/kg]

    real, intent(in), dimension(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm) :: &
      t     ! Moist static energy           [K]

    real, intent(in), dimension(nx,ny,nzm) :: &
      qv, & ! Water vapor mixing ratio                  [kg/kg]
      qcl   ! Liquid water mixing ratio (condensate)    [kg/kg]

    ! Local Variables
    real(kind=core_rknd) :: &
      wpthlp_sfc, &! w' theta_l' at surface    [(m K)/s]
      wprtp_sfc,  &! w' r_t' at surface        [(kg m)/( kg s)]
      upwp_sfc,   &! u'w' at surface           [m^2/s^2]
      vpwp_sfc     ! v'w' at surface           [m^2/s^2]

    real(kind=core_rknd), dimension(nz) :: &
      thlm,    &! Liquid water potential temperature (theta_l)  [K]
      rtm,     &! Total water mixing ratio                      [kg/kg] 
      p_in_Pa, &! Pressure                                      [Pa] 
      rho_zt,  &! Density on pressure levels                    [kg/m^3]
      rho_zm,  &! Density on momentum levels                    [kg/m^3]
      exner,   &! Exner function                                [-]
      wm_zm,   &! Imposed subs. + perturbation w on vertical vel. levels    [m/s]
      wm_zt     ! Imposed subs. + perturbation w on pressure levels         [m/s]

    real, dimension(nz) :: &
      dum     ! Dummy array for advection

    real(kind=core_rknd), allocatable, dimension(:,:) :: &
      sclrm,          & ! Array for high order passive scalars
      sclrm_forcing,  & ! Large-scale forcing array for passive scalars
      edsclrm,        & ! Array for eddy passive scalars
      edsclrm_forcing   ! Large-scale forcing array for eddy passive scalars

    real(kind=core_rknd), allocatable, dimension(:) :: &
      wpedsclrp_sfc, & ! Array for passive scalar surface flux
      wpsclrp_sfc      ! Array for high order scalar surface flux

    ! Thermo grid versions of variables on the momentum grid
    real, dimension(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nz) :: &
      wp2_zt, rtp2_zt, thlp2_zt, rtpthlp_zt, &
      wprtp_zt, wpthlp_zt, up2_zt, vp2_zt, &
      um_r4, vm_r4, um_old, vm_old ! wind arrays

    real(kind=tndcy_precision), dimension(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nz) :: &
      um_change, vm_change ! Change in u/v      [m/s^2]

    type(pdf_parameter), allocatable, dimension(:) :: &
      pdf_params ! PDF parameters       [units vary]

    integer :: err_code

    ! Array indices
    integer :: i, j, k, ig, jg, ip1, jp1, jm1, indx

!-------------------------------------------------------------------------------

    !----- Begin Code -----

#ifndef CRM
    call t_startf('advance_clubb') ! For timing
#endif

    ! Initialize err_code to CLUBB_no_error.  In the event of the singular 
    ! matrix, etc. the variable will be set to the appropriate error code 
    ! within advance_clubb_core
    err_code = CLUBB_no_error

    ! Density is in correct units
    rho_zt(2:nz) = real( rho(1:nzm), kind=core_rknd )
    rho_zt(1) = LIN_EXT( rho_zt(3), rho_zt(2), gr%zt(3), gr%zt(2), gr%zt(1) )

    rho_zm(1:nz) = real( rhow(1:nz), kind=core_rknd )

    ! Compute and extrapolate Exner function
    exner(2:nz) = 1.0_core_rknd / real( prespot(1:nzm), kind=core_rknd )
    exner(1)    = 1.0_core_rknd / LIN_EXT( exner(3), exner(2), gr%zt(3), gr%zt(2), gr%zt(1) )

    ! Allocate passive scalar arrays
    allocate( wpsclrp_sfc(sclr_dim), sclrm(nz,sclr_dim), &
              sclrm_forcing(nz,sclr_dim) )
    allocate( wpedsclrp_sfc(edsclr_dim), edsclrm(nz,edsclr_dim), &
              edsclrm_forcing(nz,edsclr_dim) )

    ! Allocate variables for the PDF closure scheme
    allocate( pdf_params(1:nz) )
 
    um_r4 = 0.0
    vm_r4 = 0.0
    do i = 1, nx, 1
      do j = 1, ny, 1

        ip1 = min( nxp1, i+1 )  ! This is redundant, but we include it for safety
        jp1 = min( nyp1, j+1 )  ! This prevents an array out of bounds error
                                !   for dvdt in a 2D simulation

        ! Average u-wind (east-west wind) to scalar points.
        um_r4(i,j,2:nz) = 0.5 * ( u(i,j,1:nzm) + u(ip1,j,1:nzm) ) + ug
        um_r4(i,j,1) = um_r4(i,j,2)

        ! Average v-wind (north-south wind) to scalar points.
        vm_r4(i,j,2:nz) = 0.5 * ( v(i,j,1:nzm) + v(i,jp1,1:nzm) ) + vg
        vm_r4(i,j,1) = vm_r4(i,j,2)
      end do
    end do

    ! Adjust the ghost points to allow for interpolation back on to 
    !   the u & v grid points
#ifndef CRM
    if ( dompi ) then
      call task_exchange( um_r4(:,:,2:nz), dimx1_s, dimx2_s, dimy1_s, dimy2_s, &
                          nzm, 3,3,3,3, ntracers+nmicro_fields+19)
      call task_exchange( vm_r4(:,:,2:nz), dimx1_s, dimx2_s, dimy1_s, dimy2_s, &
                          nzm, 3,3,3,3, ntracers+nmicro_fields+20)
    else
#endif /*CRM*/
      call bound_exchange( um_r4(:,:,2:nz), dimx1_s, dimx2_s, dimy1_s, dimy2_s, &
                           nzm, 3,3,3,3, ntracers+nmicro_fields+19)
      call bound_exchange( vm_r4(:,:,2:nz), dimx1_s, dimx2_s, dimy1_s, dimy2_s, &
                           nzm, 3,3,3,3, ntracers+nmicro_fields+20)
#ifndef CRM
    end if
#endif /*CRM*/

    ! Preserve value of u and v to calculate total change from CLUBB
    um_old = um_r4
    vm_old = vm_r4

    ! Copy the SAM precision values into CLUBB precision arrays
    um = real( um_r4, kind=core_rknd )
    vm = real( vm_r4, kind=core_rknd )

    do i=1, nx, 1

      do j=1, ny, 1
        
        ! Sample from a single column
        if ( is_a_sample_node( rank ) .and. i == x_samp_node .and. j == y_samp_node & 
            .and. lstats_clubb ) then
          l_stats = .true.
          call stats_begin_timestep( time_current-time_initial+dt_clubb )
        else
          l_stats = .false.
          l_stats_samp = .false.
        end if

        ! The 2-D flux arrays are already in the correct units
        wprtp_sfc  = real( fluxbq(i,j), kind=core_rknd )  ! [m kg/kg s] 
        wpthlp_sfc = real( fluxbt(i,j), kind=core_rknd )  ! [m K/s]
! Vince Larson set sfc momentum flux constant, as a temporary band-aid.  
! 25 Feb 2008.
        ! These are set for the purposes of computing sfc_var, but this value is
        ! not applied to the value of u and v in SAM.
        upwp_sfc = real( fluxbu(i,j), kind=core_rknd )
        vpwp_sfc = real( fluxbv(i,j), kind=core_rknd )
! End of Vince Larson's change

        ! Set the surface flux of the two scalar types to the tracer flux at the
        ! bottom of the domain, and set edsclrm to the tracer
        do indx = 1, edsclr_dim, 1
          wpedsclrp_sfc(indx) = real( fluxbtr(i,j,indx), kind=core_rknd )
          edsclrm(2:nz,indx)  = real( tracer(i,j,1:nzm,indx), kind=core_rknd )
          edsclrm(1,indx) = real( LIN_EXT( edsclrm(3,indx), edsclrm(2,indx), &
                                     gr%zt(3), gr%zt(2), gr%zt(1) ), kind=core_rknd )

          edsclrm_forcing(1:nz,indx) = 0.0_core_rknd
        end do

        do indx = 1, sclr_dim, 1
          wpsclrp_sfc(indx) = real( fluxbtr(i,j,indx), kind=core_rknd )
          sclrm(2:nz,indx)  = real( tracer(i,j,1:nzm,indx), kind=core_rknd )
          sclrm(1,indx) = LIN_EXT( sclrm(3,indx), sclrm(2,indx), &
                                 gr%zt(3), gr%zt(2), gr%zt(1) )
          sclrm_forcing(1:nz,indx) = 0.0_core_rknd
        end do

  
        ! Check for negative values of water vapor being fed from SAM into CLUBB
        if ( clubb_at_least_debug_level( 2 ) ) then
          do k=1,nzm
            if ( qv(i,j,k) < 0. ) then
              write(fstderr,*) 'SAM has fed into CLUBB negative rv at grid point i,j,k =', &
                i, j, k
            end if
          end do         

          ! Check for negative values of cloud water being fed from SAM into CLUBB
          do k=1,nzm
            if ( qcl(i,j,k)  < 0. ) then
              write(fstderr,*) 'SAM has fed into CLUBB negative qcl at grid point i,j.k =', &
                i, j, k
            end if
          end do
        end if ! clubb_at_least_debug_level( 2 )

        ! Total water. Since the SCM does not account for ice, we sum only the
        ! non-precipitating liquid and vapor

        ! Total water is the sum of non-precipitating liquid + vapor
        rtm(2:nz) = real( qv(i,j,1:nzm) + qcl(i,j,1:nzm), kind=core_rknd )
        rtm(1)    = rtm(2)

        ! Cloud water is total non-precipitating liquid
        rcm(i,j,2:nz) = real( qcl(i,j,1:nzm), kind=core_rknd )
        rcm(i,j,1)    = 0.0_core_rknd ! No below ground cloud water

        ! Note: t is moist static energy, which is not quite the same as liquid
        ! potential temperature.
        thlm(2:nz) = t2thetal( t(i,j,1:nzm), gamaz(1:nzm), &
                               qcl(i,j,1:nzm), qpl(i,j,1:nzm), &
                               qci(i,j,1:nzm), qpi(i,j,1:nzm), &
                               prespot(1:nzm) )
        thlm(1)    = thlm(2)

        ! The w variable requires no extrapolation

    ! Vince Larson added option for l_advect = .true. .  13 Mar 2008.
    !   SAM's subroutine 'subsidence' imposes wsub on t, q, u, and v.
    !   SAM advects all means using u, v, w.
    !   When implemented in a host model, CLUBB imposes wm_zm/wm_zt on higher-order
    ! moments but not means. 
    !   (l_advect=.true.) advects all higher-order moments using u, v, w.
         if ( l_advect )  then
            wm_zt(1) = 0._core_rknd
            wm_zt(2:nz) = real( wsub(1:nzm), kind=core_rknd ) ! Use this if l_advect = .true.
            wm_zm = zt2zm( wm_zt )
         else ! l_advect = .false.
            ! Higher-order moments are advected vertically but not horizontally.
            ! In principle, this could lead to undesirable accumulation.
            wm_zt(1) = 0._core_rknd ! Set ghost point to 0.
            wm_zt(2:nz) = real( wsub(1:nzm), kind=core_rknd ) ! wsub is on the t-levels
            wm_zm(1:nz) = zt2zm( wm_zt ) ! Interpolate imposed subsidence to m-levels

            ! Resolved vertical velocity is on the momentum levels
            wm_zm(1:nz) = wm_zm(1:nz) + real( w(i,j,1:nz), kind=core_rknd ) 
            ! Interpolate resolved w to t-levels
            wm_zt(1:nz) = wm_zt + zm2zt( real( w(i,j,1:nz), kind=core_rknd ) ) 
         end if
    ! End Vince Larson's commenting

        ! Add in pressure perturbation, extrapolate, & convert from mb to Pa.
        ! Vince Larson of UWM removed perturbation pressure to avoid 
        !      negative pressure at domain top in ARM9707.  22 Dec 2007.
    !   pr(2:nz) = 100. *  ( pres(1:nzm) + p(i,j,1:nzm) )
    !   pr(1)    = 100. * LIN_EXT( pres(2)+p(i,j,2), pres(1)+p(i,j,1), &
    !                              gr%zt(3), gr%zt(2), gr%zt(1) )
        P_in_Pa(2:nz) = 100._core_rknd *  real( pres(1:nzm), kind=core_rknd )
        P_in_Pa(1)    = 100._core_rknd * LIN_EXT( P_in_Pa(3), P_in_Pa(2), &
                                             gr%zt(3), gr%zt(2), gr%zt(1) )

        !  End Vince Larson's change.  

        ! Call the single column model, CLUBB
        call advance_clubb_core &
             ( l_implemented, dt_clubb, real( fcory(j), kind=core_rknd ), gr%zm(1), &
               zero(:), zero(:), zero(:), zero(:), &
               sclrm_forcing, edsclrm_forcing, wm_zm(:), wm_zt(:), &
               wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, &
               wpsclrp_sfc, wpedsclrp_sfc, & 
               P_in_Pa(:), rho_zm(:), rho_zt(:), exner(:), &
               rho_ds_zm(:), rho_ds_zt(:), invrs_rho_ds_zm(:), &
               invrs_rho_ds_zt(:), thv_ds_zm(:), thv_ds_zt(:), &
               um(i,j,:), vm(i,j,:), upwp(i,j,:), vpwp(i,j,:), up2(i,j,:), vp2(i,j,:), &
               thlm(:), rtm(:), wprtp(i,j,:), wpthlp(i,j,:), &
               wp2(i,j,:), wp3(i,j,:), rtp2(i,j,:), thlp2(i,j,:), rtpthlp(i,j,:), &
               sclrm, sclrp2(i,j,:,:), sclrprtp(i,j,:,:), sclrpthlp(i,j,:,:), &
               wpsclrp(i,j,:,:), edsclrm, err_code, &
               rcm(i,j,:), wprcp(i,j,:), cloud_frac(i,j,:), &
               rcm_in_layer(i,j,:), cloud_cover(i,j,:), pdf_params )

        ! Sample stats from a single column
        if ( is_a_sample_node( rank ) .and. i == x_samp_node .and. j == y_samp_node ) then
          call stats_end_timestep( )
        end if

        ! Check if a critical error has occured within the CLUBB model
        if ( err_code /= clubb_no_error ) then
          call task_rank_to_index( rank, ig, jg )
          write(fstderr,*) "Task #:", rank
          write(fstderr,*) "Single-column model failed at: ", "nx=", i, ";", "ny=", j, ";"
          write(fstderr,*) "x global=", i+ig, ";", "y global=", j+jg, ";"
          call task_abort( )
        end if

        ! If we're not doing a doclubbnoninter run, then we feed the results back
        !   into the 3D SAM model arrays
        if ( doclubb ) then

          ! Check for negative values of water vapor
          if ( clubb_at_least_debug_level( 2 ) ) then
            do k=1,nz
              if ( ( rtm(k) - rcm(i,j,k) ) < 0._core_rknd ) then
                write(fstderr,*) 'CLUBB has produced negative rvm at grid level k=', k
              end if
            end do
          end if ! clubb_at_least_debug_level( 2 )

          ! Re-compute vapor for total water and liquid from CLUBB
          !qv(i,j,1:nzm) = rtm(2:nz) - rcm(i,j,2:nz)
          qv_tndcy(i,j,1:nzm) = &
            ( rtm(2:nz) - rcm(i,j,2:nz) - real( qv(i,j,1:nzm), kind=core_rknd ) ) / dt_clubb

          if ( clubb_at_least_debug_level( 2 ) ) then
            ! Check for negative values of cloud water
            do k=1,nz
              if ( rcm(i,j,k)  < 0._core_rknd ) then
                write(fstderr,*) 'CLUBB has produced negative rcm at grid level k=', k
              end if
            end do
          end if ! clubb_at_least_debug_level( 2 )

          ! Re-compute qcl based on new rcm 
          !qcl(i,j,1:nzm) = rcm(i,j,2:nz)
          ! Compute tendency of total water due to CLUBB
          qc_tndcy(i,j,1:nzm) = ( rcm(i,j,2:nz) - real( qcl(i,j,1:nzm), kind=core_rknd ) ) &
                              / dt_clubb

          ! Compute moist static energy based on new thetal
!         t(i,j,1:nzm) = thetal2t( thlm(2:nz), gamaz(1:nzm), &
!                                  qcl(i,j,1:nzm), qpl(i,j,1:nzm), &
!                                  qci(i,j,1:nzm), qpi(i,j,1:nzm), &
!                                  prespot(1:nzm) )

          ! Compute tendency of moist static energy due to CLUBB
          ! Note that this formula assumes qci/qpl/qpi won't change rapidly in
          ! the time between successive clubb calls in order to avoid calling 
          ! thetal2t on at every SAM timestep -dschanen 27 Oct 08
          t_tndcy(i,j,1:nzm) = &
            ( thetal2t( thlm(2:nz), gamaz(1:nzm), rcm(i,j,2:nz), &
                      qpl(i,j,1:nzm), qci(i,j,1:nzm), qpi(i,j,1:nzm), prespot(1:nzm) ) &
              - real( t(i,j,1:nzm), kind=core_rknd ) )  / dt_clubb

          do indx = 1, edsclr_dim 
            tracer_tndcy(i,j,1:nzm,indx) = &
              ( edsclrm(2:nz,indx) - real( tracer(i,j,1:nzm,indx), kind=core_rknd ) ) &
               / dt_clubb
          end do

          do indx = 1, sclr_dim
            tracer_tndcy(i,j,1:nzm,indx) = &
              ( sclrm(2:nz,indx) - real( tracer(i,j,1:nzm,indx), kind=core_rknd ) ) / dt_clubb
          end do

        end if ! doclubb

      end do ! j

    end do   ! i

    ! De-allocate temporary arrays.  This is just in case the compiler isn't
    ! 100% Fortran 95 compliant and doesn't de-allocate this memory when it
    ! leaves the scope of advance_clubb_sgs
    deallocate( wpsclrp_sfc, sclrm )
    deallocate( wpedsclrp_sfc, edsclrm )
    deallocate( pdf_params )

    ! Copy back the value from the CLUBB precision um and vm
    um_r4 = real( um )
    vm_r4 = real( vm )

    if ( doclubb ) then

      ! Adjust the ghost points to allow for interpolation back onto the u & v grid
#ifndef CRM
      if ( dompi ) then
        call task_exchange( um_r4(:,:,2:nz), dimx1_s, dimx2_s, dimy1_s, dimy2_s, &
                            nzm, 3,3,3,3, ntracers+nmicro_fields+19)
        call task_exchange( vm_r4(:,:,2:nz), dimx1_s, dimx2_s, dimy1_s, dimy2_s, &
                            nzm, 3,3,3,3, ntracers+nmicro_fields+20)
      else
#endif /*CRM*/
        call bound_exchange( um_r4(:,:,2:nz), dimx1_s, dimx2_s, dimy1_s, dimy2_s, &
                             nzm, 3,3,3,3, ntracers+nmicro_fields+19)
        call bound_exchange( vm_r4(:,:,2:nz), dimx1_s, dimx2_s, dimy1_s, dimy2_s, &
                             nzm, 3,3,3,3, ntracers+nmicro_fields+20)
#ifndef CRM
      end if
#endif

      ! Compute the total change in u due to the CLUBB part of the code
      um_change = real( um_r4 - um_old, kind=tndcy_precision ) / dt_clubb
      vm_change = real( vm_r4 - vm_old, kind=tndcy_precision ) / dt_clubb

      ! Average the contributions of CLUBB to the wind back on to the u and v grid
      ! This has shown to make the model unstable at fine horizontal resolution.
      ! To interpolate across subdomain boundaries requires that we 
      !  transfer information using MPI (via task_exchange).
      do i=1, nx, 1
        do j=1, ny, 1
          jm1 = max( dimy1_s, j-1 ) ! For the 2D case vm wind

          ! The horiztontal grid in SAM is always evenly spaced, so we just use
          ! 0.5 *( x(n-1)+x(n) ) to interpolate back to the u,v point on the Arakawa C grid
          u_tndcy(i,j,1:nzm) = &
            0.4_tndcy_precision * & ! This is a made up coefficient to reduce numerical instability
            0.5_tndcy_precision * &
            real( um_change(i,j,2:nz) + um_change(i-1,j,2:nz), kind=tndcy_precision )
          v_tndcy(i,j,1:nzm) = &
            0.4_tndcy_precision * & ! This is a made up coefficient to reduce numerical instability
            0.5_tndcy_precision * &
            real( vm_change(i,j,2:nz) + vm_change(i,jm1,2:nz), kind=tndcy_precision )

        end do ! j
      
      end do ! i

    end if ! doclubb
  

! Vince Larson attempted to advect higher-order moments horizontally.  
!     26 Feb 2008.

! Horizontal advection of higher-order moments.

! The following method has the drawback of requiring two interpolations, 
!    which unnecesarily smooths the fields in the vertical.
! In preparation for advection, interpolate to thermodynamic (scalar) vertical gridpoints.
!   (wp3 is already on the thermodynamic gridpoints.)


!print*, 'Before advection, wp2(nx,ny,:) =', wp2(nx,ny,:)
! For now we default to not doing this, because the interpolation seems to cause
! and artificial rise in fields such as moisture at a coarse model resolution.
! -dschanen 29 Apr 2008
  if ( l_advect ) then

    do i=1, nx, 1
      do j=1, ny, 1

        wp2_zt(i,j,:)     = real( zm2zt( wp2(i,j,:) ) )
        up2_zt(i,j,:)     = real( zm2zt( up2(i,j,:) ) )
        vp2_zt(i,j,:)     = real( zm2zt( vp2(i,j,:) ) )
        rtp2_zt(i,j,:)    = real( zm2zt( rtp2(i,j,:) ) )
        thlp2_zt(i,j,:)   = real( zm2zt( thlp2(i,j,:) ) )
        rtpthlp_zt(i,j,:) = real( zm2zt( rtpthlp(i,j,:) ) )
        wprtp_zt(i,j,:)   = real( zm2zt( wprtp(i,j,:) ) )
        wpthlp_zt(i,j,:)  = real( zm2zt( wpthlp(i,j,:) ) )

      end do ! j
    end do   ! i

#ifndef CRM
    if ( dompi ) then
  
      call task_exchange( wp2_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+10 )
      call task_exchange( rtp2_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+11 )
      call task_exchange( thlp2_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+12 )
      call task_exchange( rtpthlp_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+13 )
      call task_exchange( wprtp_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+14 )
      call task_exchange( wpthlp_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+15 )
      call task_exchange( wp3(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+16 )
      call task_exchange( up2_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+17 )
      call task_exchange( vp2_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+18 )
    else
#endif /*CRM*/

      call bound_exchange( wp2_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+10 )
      call bound_exchange( rtp2_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+11 )
      call bound_exchange( thlp2_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+12 )
      call bound_exchange( rtpthlp_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+13 )
      call bound_exchange( wprtp_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+14 )
      call bound_exchange( wpthlp_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+15 )
      call bound_exchange( wp3(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+16 )
      call bound_exchange( up2_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+17 )
      call bound_exchange( vp2_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),   &
        dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm, 3,3,3,3, &
        ntracers+nmicro_fields+18 )

#ifndef CRM
    end if
#endif

    ! Now call the standard SAM advection subroutine for scalars
    call advect_scalar( wp2_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz), &
                        dum(1:nz), dum(1:nz), dum(1:nzm), &
                        dum(1:nzm), dum(1:nzm), .false. )

    call advect_scalar( wp3(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz),  &
                        dum(1:nz), dum(1:nz), dum(1:nzm), &
                        dum(1:nzm), dum(1:nzm), .false. )

    call advect_scalar( rtp2_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz), &
                        dum(1:nz), dum(1:nz), dum(1:nzm), &
                        dum(1:nzm), dum(1:nzm), .false. )

    call advect_scalar( thlp2_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz), &
                        dum(1:nz), dum(1:nz), dum(1:nzm), &
                        dum(1:nzm), dum(1:nzm), .false. )

    call advect_scalar( rtpthlp_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz), &
                        dum(1:nz), dum(1:nz), dum(1:nzm), &
                        dum(1:nzm), dum(1:nzm), .false. )

    call advect_scalar( wprtp_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz), &
                        dum(1:nz), dum(1:nz), dum(1:nzm), &
                        dum(1:nzm), dum(1:nzm), .false. )

    call advect_scalar( wpthlp_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz), &
                        dum(1:nz), dum(1:nz), dum(1:nzm), &
                        dum(1:nzm), dum(1:nzm), .false. )

    call advect_scalar( up2_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz), &
                        dum(1:nz), dum(1:nz), dum(1:nzm), &
                        dum(1:nzm), dum(1:nzm), .false. )

    call advect_scalar( vp2_zt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,2:nz), &
                        dum(1:nz), dum(1:nz), dum(1:nzm), &
                        dum(1:nzm), dum(1:nzm), .false. )

!print*, 'After advect, wp2_zt(dimx2_s,dimy2_s,:) =', wp2_zt(dimx2_s,dimy2_s,:)
!
!do i=dimx1_s, dimx2_s, 1
!  do j=dimy1_s, dimy2_s, 1
!    if ( any ( rtp2_zt(i,j,:) < 0.0 ) ) then
!    print*, 'After advect, rtp2_zt at ', i, j, " = ",  rtp2_zt(i,j,:)
!    end if
!  end do ! i
!end do ! j
! Now interpolate back to momentum gridpoints.
!   (wp3 is already on the thermodynamic gridpoints.)
! do i=dimx1_s, dimx2_s, 1
!   do j=dimy1_s, dimy2_s, 1
    do i=1, nx, 1
      do j=1, ny, 1
     
        wp2(i,j,:)     = zt2zm( real( wp2_zt(i,j,:), kind=core_rknd ) )
        up2(i,j,:)     = zt2zm( real( up2_zt(i,j,:), kind=core_rknd ) )
        vp2(i,j,:)     = zt2zm( real( vp2_zt(i,j,:), kind=core_rknd ) )
        rtp2(i,j,:)    = zt2zm( real( rtp2_zt(i,j,:), kind=core_rknd ) )
        thlp2(i,j,:)   = zt2zm( real( thlp2_zt(i,j,:), kind=core_rknd ) )
        rtpthlp(i,j,:) = zt2zm( real( rtpthlp_zt(i,j,:), kind=core_rknd ) )
        wprtp(i,j,:)   = zt2zm( real( wprtp_zt(i,j,:), kind=core_rknd ) )
        wpthlp(i,j,:)  = zt2zm( real( wpthlp_zt(i,j,:), kind=core_rknd ) )

      end do ! j
    end do   ! i

    ! Clip variances where the top point is negative
    where ( wp2(:,:,nz) < 0._core_rknd ) wp2(:,:,nz) = 0._core_rknd
    where ( up2(:,:,nz) < 0._core_rknd ) up2(:,:,nz) = 0._core_rknd
    where ( vp2(:,:,nz) < 0._core_rknd ) vp2(:,:,nz) = 0._core_rknd
    where ( rtp2(:,:,nz) < 0._core_rknd ) rtp2(:,:,nz) = 0._core_rknd
    where ( thlp2(:,:,nz) < 0._core_rknd ) thlp2(:,:,nz) = 0._core_rknd

    ! Clip variances where the bottom point is negative
    where ( wp2(:,:,1) < 0._core_rknd ) wp2(:,:,1) = 0._core_rknd
    where ( up2(:,:,1) < 0._core_rknd ) up2(:,:,1) = 0._core_rknd
    where ( vp2(:,:,1) < 0._core_rknd ) vp2(:,:,1) = 0._core_rknd
    where ( rtp2(:,:,1) < 0._core_rknd ) rtp2(:,:,1) = 0._core_rknd
    where ( thlp2(:,:,1) < 0._core_rknd ) thlp2(:,:,1) = 0._core_rknd


!do i=1, nx, 1
!  do j=1, ny, 1
!    if ( any ( rtp2(i,j,:) < 0.0 ) ) then
!    print*, 'After interp, rtp2 at ', i, j, " = ",  rtp2(i,j,:)
!    end if
!  end do ! i
!end do ! j
!
!print*, 'After interp back, wp2(nx,ny,:) =', wp2(nx,ny,:)
!! End of Vince Larson's changes.
    end if ! ladvect

#ifndef CRM
    call t_stopf('advance_clubb') ! For timing
#endif

    return
  end subroutine advance_clubb_sgs

!-------------------------------------------------------------------------------
  subroutine apply_clubb_sgs_tndcy( dt, t, qv, qcl, dudt, dvdt )

    use grid, only: &
      nx, nxp1, ny, nyp1, dimx1_s, dimx2_s, dimy1_s, dimy2_s, nz, nzm, na, &
      rank

    use domain, only: &
      ntracers

#ifndef CRM
    use tracers, only: &
#else
    use crmtracers, only: &
#endif 
      tracer

    use clubbvars, only: &
      u_tndcy, & ! CLUBB contribution to the x wind
      v_tndcy, & ! CLUBB contribution to the y wind
      t_tndcy, & ! CLUBB contribution to moist static energy
      qc_tndcy,& ! CLUBB contribution to liquid water mixing ratio
      qv_tndcy   ! CLUBB contribution to vapor water mixing ratio

    use clubbvars, only: &
      tracer_tndcy

    use clubbvars, only: &
      sclr_dim, & ! Constant(s)
      edsclr_dim

    use clubbvars, only: &
     rho_ds_zt, & ! Variable(s)
     rho_ds_zm

    use error_code, only: clubb_at_least_debug_level

    use fill_holes, only: fill_holes_driver

    implicit none

    intrinsic :: any

    ! In variables
    real(kind=time_precision), intent(in) :: &
      dt ! Timestep [s]

    ! In/Out variables
    real, intent(inout), dimension(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm) :: &
      t     ! Moist static energy           [K]

    real, intent(inout), dimension(nx,ny,nzm) :: &
      qv, & ! Water vapor mixing ratio                  [kg/kg]
      qcl   ! Liquid water mixing ratio (condensate)    [kg/kg]

    real, intent(inout), dimension(nxp1,ny,nzm,3) :: &
      dudt ! u wind tendency [m/s^2]

    real, intent(inout), dimension(nx,nyp1,nzm,3) :: &
      dvdt ! v wind tendency [m/s^2]

    ! Local Variables
    real(kind=core_rknd), dimension(nz) :: tmpqv, tmpqcl

    real(kind=core_rknd) :: threshold ! Threshold on clipping [units vary]

    integer :: i, j, ig, jg

    ! --- Begin Code --- 

#ifndef CRM
    call t_startf('apply_clubb_sgs_tndcy') ! For timing
#endif

    ! Since dudt/dvdt are already time tendencies, we just add the contribution 
    ! to the existing SAM contribution
    dudt(1:nx,1:ny,1:nzm,na) = dudt(1:nx,1:ny,1:nzm,na) + real( u_tndcy(1:nx,1:ny,1:nzm) )
    dvdt(1:nx,1:ny,1:nzm,na) = dvdt(1:nx,1:ny,1:nzm,na) + real( v_tndcy(1:nx,1:ny,1:nzm) )

    tmpqv  = 0.0_core_rknd
    tmpqcl = 0.0_core_rknd

    ! Add clubb tendency to qv, qc, t, and tracers
    do i = 1, nx, 1
      do j = 1, ny, 1

        t(i,j,1:nzm) = t(i,j,1:nzm) + real( dt*t_tndcy(i,j,1:nzm) )

        tmpqv(2:nz)  = real( qv(i,j,1:nzm), kind=core_rknd )  + dt*qv_tndcy(i,j,1:nzm)
        tmpqcl(2:nz) = real( qcl(i,j,1:nzm), kind=core_rknd ) + dt*qc_tndcy(i,j,1:nzm)

        if ( edsclr_dim > 0 .or. sclr_dim > 0 ) then
          tracer(i,j,1:nzm,1:ntracers) = tracer(i,j,1:nzm,1:ntracers) &
            + real( dt*tracer_tndcy(i,j,1:nzm,1:ntracers) )
        end if

        ! Apply hole-filling scheme to qv as needed
        threshold = 0._core_rknd
        if ( any( tmpqv(2:nz) < threshold ) ) then

          ! CLUBB's tendency in this column will produce a negative vapor water,
          ! so we apply hole-filling
          if ( clubb_at_least_debug_level( 1 ) ) then
            call task_rank_to_index( rank, ig, jg )
            write(0,*) "Task #:", rank
            write(0,*) "Applying hole-filling scheme to vapor water mixing ratio at:", &
              "nx=", i, ";", "ny=", j, ";"
            write(0,*) "x global=", i+ig, ";", "y global=", j+jg, ";"
          end if

          call fill_holes_driver( 2, threshold, "zt", rho_ds_zt, rho_ds_zm, tmpqv )

        end if

        ! Update qv
        qv(i,j,1:nzm) = real( tmpqv(2:nz) )

        threshold = 0._core_rknd
        ! Apply hole-filling scheme to qcl as needed
        if ( any( tmpqcl(2:nz) < threshold ) ) then

          ! CLUBB's tendency in this column will produce a negative cloud water,
          ! so we apply hole-filling
          if ( clubb_at_least_debug_level( 1 ) ) then
            call task_rank_to_index( rank, ig, jg )
            write(0,*) "Task #:", rank
            write(0,*) "Applying hole-filling scheme to cloud water mixing ratio at:", &
              "nx=", i, ";", "ny=", j, ";"
            write(0,*) "x global=", i+ig, ";", "y global=", j+jg, ";"
          end if

          call fill_holes_driver( 2, threshold, "zt", rho_ds_zt, rho_ds_zm, tmpqcl )

        end if

        ! Update qcl
        qcl(i,j,1:nzm) = real( tmpqcl(2:nz) )

      end do ! j = 1, ny
    end do ! i = 1, nx

#ifndef CRM
    call t_stopf('apply_clubb_sgs_tndcy') ! For timing
#endif

    return
  end subroutine apply_clubb_sgs_tndcy

!-------------------------------------------------------------------------------
  subroutine apply_clubb_sgs_tndcy_mom( dt, dudt, dvdt )

    use grid, only: &
      nx, nxp1, ny, nyp1, dimx1_s, dimx2_s, dimy1_s, dimy2_s, nz, nzm, na, &
      rank

    use clubbvars, only: &
      u_tndcy, & ! CLUBB contribution to the x wind
      v_tndcy  ! CLUBB contribution to the y wind

    implicit none

    intrinsic :: any

    ! In variables
    real(kind=time_precision), intent(in) :: &
      dt ! Timestep [s]

    real, intent(inout), dimension(nxp1,ny,nzm,3) :: &
      dudt ! u wind tendency [m/s^2]

    real, intent(inout), dimension(nx,nyp1,nzm,3) :: &
      dvdt ! v wind tendency [m/s^2]

!   local variables
    integer :: i, j, ig, jg

    ! --- Begin Code --- 

#ifndef CRM
    call t_startf('apply_clubb_sgs_tndcy_mom') ! For timing
#endif

    ! Since dudt/dvdt are already time tendencies, we just add the contribution 
    ! to the existing SAM contribution
    dudt(1:nx,1:ny,1:nzm,na) = dudt(1:nx,1:ny,1:nzm,na) + real( u_tndcy(1:nx,1:ny,1:nzm) )
    dvdt(1:nx,1:ny,1:nzm,na) = dvdt(1:nx,1:ny,1:nzm,na) + real( v_tndcy(1:nx,1:ny,1:nzm) )

#ifndef CRM
    call t_stopf('apply_clubb_sgs_tndcy_mom') ! For timing
#endif

    return
  end subroutine apply_clubb_sgs_tndcy_mom

!-------------------------------------------------------------------------------
  subroutine apply_clubb_sgs_tndcy_scalar( dt, t, qv, qcl)

    use grid, only: &
      nx, nxp1, ny, nyp1, dimx1_s, dimx2_s, dimy1_s, dimy2_s, nz, nzm, na, &
      rank, adz, dz, doclubb_sfc_fluxes

    use vars, only: rho

    use domain, only: &
      ntracers

#ifndef CRM
    use tracers, only: &
#else
    use crmtracers, only: &
#endif 
      tracer

    use clubbvars, only: &
      t_tndcy, & ! CLUBB contribution to moist static energy
      qc_tndcy,& ! CLUBB contribution to liquid water mixing ratio
      qv_tndcy   ! CLUBB contribution to vapor water mixing ratio

    use clubbvars, only: &
      tracer_tndcy

    use clubbvars, only: &
      sclr_dim, & ! Constant(s)
      edsclr_dim

    use clubbvars, only: &
     rho_ds_zt, & ! Variable(s)
     rho_ds_zm

    use error_code, only: clubb_at_least_debug_level

    use fill_holes, only: fill_holes_driver

    implicit none

    intrinsic :: any

    ! In variables
    real(kind=time_precision), intent(in) :: &
      dt ! Timestep [s]

    ! In/Out variables
    real, intent(inout), dimension(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm) :: &
      t     ! Moist static energy           [K]

    real, intent(inout), dimension(nx,ny,nzm) :: &
      qv, & ! Water vapor mixing ratio                  [kg/kg]
      qcl   ! Liquid water mixing ratio (condensate)    [kg/kg]

    ! Local Variables
    real(kind=core_rknd), dimension(nz) :: tmpqv, tmpqcl

    real(kind=core_rknd) :: threshold ! Threshold on clipping [units vary]

    real(kind=core_rknd), dimension(2) :: t_total 

    integer :: i, j, ig, jg, k

    ! --- Begin Code --- 

#ifndef CRM
    call t_startf('apply_clubb_sgs_tndcy_scalar') ! For timing
#endif

    tmpqv  = 0.0_core_rknd
    tmpqcl = 0.0_core_rknd

    ! Add clubb tendency to qv, qc, t, and tracers
    do i = 1, nx, 1
      do j = 1, ny, 1

! add energy conservation check and fix for CLUBB
! Minghuai Wang, 2012-06
        t_total = 0.0_core_rknd
        do k=1, nzm
          t_total(1) = t_total(1) +real(t(i,j,k)*rho(k)*adz(k)*dz, kind=core_rknd)
          t(i,j,k) = t(i,j,k) + real( dt*t_tndcy(i,j,k) )
          t_total(2) = t_total(2) +real(t(i,j,k)*rho(k)*adz(k)*dz, kind=core_rknd)
        end do
        if(abs(t_total(2)-t_total(1))/t_total(1).gt.5.0e-7) then
          write(0, *) 'energy conervation issue in clubb', i,j,   & 
             abs(t_total(2)-t_total(1))/t_total(1),  t_total(1)
        end if
        if(.not.doclubb_sfc_fluxes) then
           t(i,j,1:nzm) = t(i,j,1:nzm) * real(t_total(1)/t_total(2))
        else 
           write(0, *) 'need add surface fluxes in energy conservation fix'
           stop
        end if

!        t(i,j,1:nzm) = t(i,j,1:nzm) + real( dt*t_tndcy(i,j,1:nzm) )

        tmpqv(2:nz)  = real( qv(i,j,1:nzm), kind=core_rknd )  + dt*qv_tndcy(i,j,1:nzm)
        tmpqcl(2:nz) = real( qcl(i,j,1:nzm), kind=core_rknd ) + dt*qc_tndcy(i,j,1:nzm)

        if ( edsclr_dim > 0 .or. sclr_dim > 0 ) then
          tracer(i,j,1:nzm,1:ntracers) = tracer(i,j,1:nzm,1:ntracers) &
            + real( dt*tracer_tndcy(i,j,1:nzm,1:ntracers) )
        end if

        ! Apply hole-filling scheme to qv as needed
        threshold = 0._core_rknd
        if ( any( tmpqv(2:nz) < threshold ) ) then

          ! CLUBB's tendency in this column will produce a negative vapor water,
          ! so we apply hole-filling
          if ( clubb_at_least_debug_level( 1 ) ) then
            call task_rank_to_index( rank, ig, jg )
            write(0,*) "Task #:", rank
            write(0,*) "Applying hole-filling scheme to vapor water mixing ratio at:", &
              "nx=", i, ";", "ny=", j, ";"
            write(0,*) "x global=", i+ig, ";", "y global=", j+jg, ";"
          end if

          call fill_holes_driver( 2, threshold, "zt", rho_ds_zt, rho_ds_zm, tmpqv )

        end if

        ! Update qv
        qv(i,j,1:nzm) = real( tmpqv(2:nz) )

        threshold = 0._core_rknd
        ! Apply hole-filling scheme to qcl as needed
        if ( any( tmpqcl(2:nz) < threshold ) ) then

          ! CLUBB's tendency in this column will produce a negative cloud water,
          ! so we apply hole-filling
          if ( clubb_at_least_debug_level( 1 ) ) then
            call task_rank_to_index( rank, ig, jg )
            write(0,*) "Task #:", rank
            write(0,*) "Applying hole-filling scheme to cloud water mixing ratio at:", &
              "nx=", i, ";", "ny=", j, ";"
            write(0,*) "x global=", i+ig, ";", "y global=", j+jg, ";"
          end if

          call fill_holes_driver( 2, threshold, "zt", rho_ds_zt, rho_ds_zm, tmpqcl )

        end if

        ! Update qcl
        qcl(i,j,1:nzm) = real( tmpqcl(2:nz) )

      end do ! j = 1, ny
    end do ! i = 1, nx

#ifndef CRM
    call t_stopf('apply_clubb_sgs_tndcy_scalar') ! For timing
#endif

    return
  end subroutine apply_clubb_sgs_tndcy_scalar

!-------------------------------------------------------------------------------
  subroutine clubb_sgs_cleanup( )
!   Description:
!     De-allocate memory and exit.
!-------------------------------------------------------------------------------
    use grid, only: rank

    use stats_subs, only: stats_finalize

    implicit none

    !----- Begin Code -----

    call cleanup_clubb_core( .true. )

    if ( is_a_sample_node( rank ) ) then
      call stats_finalize( )
    end if

    return
  end subroutine clubb_sgs_cleanup

!-------------------------------------------------------------------------------
  elemental function t2thetal( t, gamaz, qcl, qpl, qci, qpi, prespot ) &
    result( thl )
! Description:
!   Convert moist static energy into the liquid potential temperature 
!   used in CLUBB.
!-------------------------------------------------------------------------------
    use params, only: &
      fac_cond, & ! Variables
      fac_sub

    implicit none

    ! Input variables
    real, intent(in) :: &
      t,            & ! Moist static energy                 [K]
      gamaz,        & ! grav/Cp*z                           [m]
      qcl,          & ! Cloud water mixing ration           [kg/kg]
      qpl,          & ! Rain water mixing ratio (liquid)    [kg/kg]
      qci,          & ! Cloud water mixing ratio (ice)      [kg/kg]
      qpi,          & ! Snow+Graupel mixing ratio           [kg/kg]
      prespot         ! Exner^-1                            [-]

    ! Result
    real(kind=core_rknd) :: thl     ! Liquid pot. temperature       [K]

    real :: tabs    ! Absolute temp.      [K]

    !----- Begin Code -----

    ! Compute absolute temperature from t 
    ! Formula comes from module diagnose.
    tabs =  t - gamaz + fac_cond * ( qcl + qpl ) + fac_sub * ( qci + qpi )

    ! Compute thetal (don't include ice because CLUBB doesn't) 
    thl = real( prespot * ( tabs - fac_cond * qcl ), kind=core_rknd )

    return
  end function t2thetal

!-------------------------------------------------------------------------------
  elemental function thetal2t( thl, gamaz, qcl, qpl, qci, qpi, prespot ) &
    result( t )

! Description:
!   Convert liquid potential temperature into moist static energy.
! References:
!   None
!-------------------------------------------------------------------------------
  use params, only: &
    fac_cond, & ! Variables
    fac_sub

  implicit none

  ! Input Variables
  real(kind=core_rknd), intent(in) :: &
    thl,          & ! Liquid potential temperature        [K]
    qcl             ! Cloud water mixing ration           [kg/kg]

  real, intent(in) :: &
    gamaz,        & ! grav/Cp*z                           [m]
    qpl,          & ! Rain water mixing ratio (liquid)    [kg/kg]
    qci,          & ! Cloud water mixing ratio (ice)      [kg/kg]
    qpi,          & ! Snow+Graupel mixing ratio           [kg/kg]
    prespot         ! Exner^-1                            [-]

    ! Result
    real(kind=core_rknd) :: t ! Moist static energy         [K]

    real(kind=core_rknd) :: &
      tabs,   & ! Absolute temp.      [K]
      theta     ! Pot. temp.          [K]

    !----- Begin Code -----

    ! Compute absolute temperature from thl
    ! Use fac_cond since CLUBB's thl does not account for ice
    theta = thl + real( prespot * fac_cond, kind=core_rknd ) * qcl
    tabs  = theta / real( prespot, kind=core_rknd )
    ! Compute moist static energy 
    ! Formula comes from module diagnose
    t = tabs + real( gamaz, kind=core_rknd ) &
             - real( fac_cond, kind=core_rknd ) * ( qcl + real( qpl, kind=core_rknd ) ) &
             - real( fac_sub * ( qci + qpi ), kind=core_rknd )

    return
  end function thetal2t

!-------------------------------------------------------------------------------
  FUNCTION LIN_EXT( var_high, var_low, height_high, height_low, height_ext )

! Author: Brian M. Griffin,  UW Milwaukee

! References: None

! Description:
! This function computes a linear extension of the value of variable.
! Given two known values of a variable at two height values, the value
! of that variable at a height outside of those two height levels 
! (rather than a height between those two height levels) is computed.
!
! Here is a diagram:
!
!  -------------------------------- Height to be extended to; linear extension
!
!  ################################ Height high, know variable value
!
!
!
!  ################################ Height low, know variable value
!
!
!
!  -------------------------------- Height to be extended to; linear extension
!
!
! FORMULA:
!
! variable(@ Height extension) =
!
! [ (variable(@ Height high) - variable(@ Height low)) / (Height high - Height low) ]
! * (Height extension - Height high)  +  variable(@ Height high)
!-------------------------------------------------------------------------------

    IMPLICIT NONE

    ! Input Variables
    REAL(kind=core_rknd), INTENT(IN):: var_high
    REAL(kind=core_rknd), INTENT(IN):: var_low
    REAL(kind=core_rknd), INTENT(IN):: height_high
    REAL(kind=core_rknd), INTENT(IN):: height_low
    REAL(kind=core_rknd), INTENT(IN):: height_ext

    ! Output Variable
    REAL(kind=core_rknd):: lin_ext

    !----- Begin Code -----

    lin_ext = ( var_high - var_low ) / ( height_high - height_low ) &
         * ( height_ext - height_high ) + var_high

    RETURN
  END FUNCTION LIN_EXT

  !-----------------------------------------------------------------------------
  logical function is_a_sample_node( rank )

  ! Description:
  !   Determine if we're output single-columns stats from this node.
  ! References:
  !   None
  !-----------------------------------------------------------------------------

    implicit none

    ! External
    intrinsic :: any, spread, size

    ! Input Variable
    integer, intent(in) :: rank

    integer :: iter

    ! ---- Begin Code ----

    ! Initialize
    is_a_sample_node = .false.

    ! Determine if we're sampling a column of stats from this node
    do iter = 1, size( sample_nodes )
      if ( sample_nodes(iter) == rank ) then
        is_a_sample_node = .true.
        exit
      end if
    end do

    return
  end function is_a_sample_node
 !-----------------------------------------------------------------------------
  subroutine get_sample_points( rank, i, j )

  ! Description:
  !   Output the local x and y location to be output for this particular node.
  ! 
  ! References:
  !   None
  !-----------------------------------------------------------------------------

    implicit none

    ! Input Variable
    integer, intent(in) :: rank

    ! Output Variables
    integer, intent(out) :: i, j

    integer :: iter

    ! ---- Begin Code ----

    i = -1
    j = -1
    do iter = 1, size( sample_nodes )
      if ( sample_nodes(iter) == rank ) then
        i = x_samp(iter); j = y_samp(iter)
        exit
      end if
    end do

    return
  end subroutine get_sample_points

#endif /*CLUBB_CRM*/
end module clubb_sgs
