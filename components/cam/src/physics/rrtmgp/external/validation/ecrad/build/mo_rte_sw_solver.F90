module mo_sw_solver
  !
  ! RTE modules
  !
  use mo_rte_kind,           only: wp, wl
  use mo_optical_props,      only: ty_optical_props_arry, ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
  !
  ! ECRAD modules
  !
  use parkind1,                 only : jprb ! Working precision
  use radiation_config, only         : config_type
  use radiation_single_level, only   : single_level_type
  use radiation_thermodynamics, only : thermodynamics_type
  use radiation_cloud, only          : cloud_type
  use radiation_flux, only           : flux_type
  use radiation_homogeneous_sw, only : solver_homogeneous_sw


  implicit none
  private

  public :: sw_solver

contains
  ! ---------------------------------------------------------------
  !
  ! Solution to the radiative transfer equation assuming internal emission
  !
  function sw_solver(ncol, nlay, ngpt, top_is_1,           &
                     atmos, mu0, sfc_alb_dir, sfc_alb_dif, &
                     inc_flux, flux_up, flux_dn, flux_dir, inc_flux_dif)
    integer,                         intent( in) :: ncol, nlay, ngpt !< Number of columns, layers, g-points
    logical,                         intent( in) :: top_is_1         !< True if arrays are indexed top to bottom.
    class(ty_optical_props_arry),    intent( in) :: atmos            ! Optical properties of the atmosphere
    real(wp), dimension(ncol),       intent( in) :: mu0              !< cosine of solar zenith angle
    real(wp), dimension(ncol,ngpt),  intent( in) :: sfc_alb_dir, sfc_alb_dif
                                                                     !< surface albedo for direct and diffuse radiation
    real(wp), dimension(ncol,ngpt),  intent( in) :: inc_flux          !< direct beam incident flux at top-of-atmosphere [W/m2]
    real(wp), dimension(ncol,nlay+1,ngpt), &
                                     intent(inout) :: flux_up, flux_dn, &  ! Fluxes [W/m2]
                                                      flux_dir             ! Downward direct
                                                                         ! Top level (= merge(1, nlay+1, top_is_1)
                                                                         ! must contain incident flux boundary condition
    real(wp), dimension(ncol,ngpt), optional, &
                                    intent( in) :: inc_flux_dif     !< diffuse incident flux at top-of-atmosphere [W/m2]
    character(len=128)                           :: sw_solver

    ! --------------------------------------------------
    ! Derived types for the inputs to ECRAD
    type(config_type)         :: config
    type(single_level_type)   :: single_level   ! need lw_emissivity, sw_albedo, sw_albedo_direct, cos_sza
    type(thermodynamics_type) :: thermodynamics ! not used
    type(cloud_type)          :: cloud          ! need only cloud%fraction
    ! Derived type containing outputs from the radiation scheme
    type(flux_type)           :: flux           !

    real(jprb), dimension(:,:,:), allocatable :: od, ssa, g, od_cloud, ssa_cloud, g_cloud, planck_hl
    real(jprb), dimension(:,:  ), allocatable :: incoming_sw

    integer :: icol, ilay, igpt, istartcol, iendcol
    ! --------------------------------------------------
    ! Calculations shared by all solvers ...
    !  ... are there any? We should leave anything that depends on spec_cfg elsewhere so this is pure RT.
    sw_solver = ""
    if(.not. top_is_1) then
      sw_solver = "sw_solver: atmosphere has to be ordered top to bottom (top_is_1 = .false.) for ECRAD"
      return
    end if

    !
    ! Incident flux
    !
    flux_dir(:,1,:) = inc_flux(:,:) * spread(mu0, 2, ngpt)
    allocate(incoming_sw(ngpt, ncol))
    incoming_sw = transpose(inc_flux)

    istartcol = 1; iendcol = ncol
    !
    ! ECRAD types
    ! thermo type is unused

    ! cloud
    allocate(cloud%fraction(ncol,nlay))
    cloud%fraction(1:ncol,1:nlay) = 0._wp

    ! single-level
    allocate(single_level%cos_sza(ncol), &
             single_level%sw_albedo(ncol,ngpt), &
             single_level%sw_albedo_direct(ncol,ngpt))
    single_level%cos_sza(1:ncol)  = mu0
    single_level%sw_albedo        = sfc_alb_dif
    single_level%sw_albedo_direct = sfc_alb_dir

    ! config type
    config%do_clear = .true.
    config%do_sw    = .true.
    config%do_lw    = .false.
    config%do_save_gpoint_flux = .true.
    config%do_save_spectral_flux = .true.
    config%n_spec_sw             = ngpt
    config%n_bands_sw            = ngpt ! Would normally be the number of bands but we don't have access
    config%n_g_sw                = ngpt
    config%i_band_from_reordered_g_sw = [(igpt, igpt = 1, ngpt)]
    config%i_albedo_from_band_sw      = [(igpt, igpt = 1, ngpt)]

    call flux%allocate(config, 1, ncol, nlay)

    !
    !  Optical property arrays are dimensioned ngpt, nlev, ncol
    !     although emissivity and albedo are ncol, nbnd
    !
    ! Tranpose arrays from ncol, nlay, ngpt to ngpt, nlay, ncol
    allocate(od (ngpt, nlay, ncol),  od_cloud(ngpt, nlay, ncol), &
             ssa(ngpt, nlay, ncol), ssa_cloud(ngpt, nlay, ncol), &
             g  (ngpt, nlay, ncol),   g_cloud(ngpt, nlay, ncol))
    select type(atmos)
      type is (ty_optical_props_1scl)
        config%do_lw_aerosol_scattering = .false.
        config%do_lw_cloud_scattering   = .false.
        do igpt = 1, ngpt
          do ilay = 1, nlay
            do icol = 1, ncol
              od(igpt, ilay, icol) = atmos%tau(icol, ilay, igpt)
              ssa(igpt, ilay, icol) = 0._wp
              g  (igpt, ilay, icol) = 0._wp
              od_cloud (igpt, ilay, icol) = 0._wp
              ssa_cloud(igpt, ilay, icol) = 0._wp
              g_cloud  (igpt, ilay, icol) = 0._wp
            end do
          end do
        end do
      type is (ty_optical_props_2str)
        config%do_lw_aerosol_scattering = .true.
        config%do_lw_cloud_scattering   = .true.
        do igpt = 1, ngpt
          do ilay = 1, nlay
            do icol = 1, ncol
              od (igpt, ilay, icol) = atmos%tau(icol, ilay, igpt)
              ssa(igpt, ilay, icol) = atmos%ssa(icol, ilay, igpt)
              g  (igpt, ilay, icol) = atmos%g  (icol, ilay, igpt)
              od_cloud (igpt, ilay, icol) = 0._wp
              ssa_cloud(igpt, ilay, icol) = 0._wp
              g_cloud  (igpt, ilay, icol) = 0._wp
            end do
          end do
        end do
      class default
        call stop_on_err("Don't understand this type of atmosphere")
    end select

    call solver_homogeneous_sw(nlay,istartcol,iendcol, &
         &  config, single_level, thermodynamics, cloud, &
         &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, incoming_sw, &
         &  flux)
    !
    ! Homogeneous solver reports only broadband fluxes summed across g-points
    !   Broadband sums will be correct below; band fluxes will be all messed up
    !
    flux_up(1:ncol,1:nlay+1,1) = flux%sw_up; flux_up(1:ncol,1:nlay+1,2:ngpt) = 0._wp
    flux_dn(1:ncol,1:nlay+1,1) = flux%sw_dn; flux_dn(1:ncol,1:nlay+1,2:ngpt) = 0._wp
    flux_dir(1:ncol,1:nlay+1,1) = flux%sw_dn_direct; flux_dir(1:ncol,1:nlay+1,2:ngpt) = 0._wp

  end function sw_solver
  ! ---------------------------------------------------------------
end module mo_sw_solver
