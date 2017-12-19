! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2017,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!

!
! This module provides an interface to RRTMGP for a common use case --
!   users want to start from gas concentrations, pressures, and temperatures,
!   and compute clear-sky (aerosol plus gases) and all-sky fluxes.
! The routines here have the same names as those in mo_rrtmgp_[ls]w; normally users
!   will use either this module or the underling modules, but not both
!
module mo_rrtmgp_clr_all_sky
  use mo_rte_kind,   only: wp
  use mo_spectral_disc, only: ty_spectral_disc
  use mo_gas_optics, &
                        only: ty_gas_optics_specification
  use mo_gas_concentrations, &
                        only: ty_gas_concs
  use mo_optical_props, only: ty_optical_props, &
                              ty_optical_props_arry, ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
  use mo_fluxes,        only: ty_fluxes
  use mo_rte_lw,     only: base_rte_lw_init => rte_lw_init, &
                              base_rte_lw      => rte_lw
  use mo_rte_sw,     only: base_rte_sw_init => rte_sw_init, &
                              base_rte_sw      => rte_sw
  implicit none
  private

  public :: rte_lw_init, rte_lw, &
            rte_sw_init, rte_sw

  !
  ! Configuration information
  !
  integer :: n_lw_streams = 0 ! Number of streams: 0 means no scattering
                            ! Require even, positive value when setting

  integer :: n_sw_streams = 2 ! Number of streams: 0 means no scattering
                            ! Require even, positive value when setting

contains
  ! --------------------------------------------------
  ! Initialization functions
  ! --------------------------------------------------
  function rte_lw_init(nlwstreams, nangles) result(error_msg)
    integer,           optional, intent( in) :: nlwstreams ! Scattering/no scattering
    integer,           optional, intent( in) :: nangles    ! number of quadrature angles for
                                                           ! no-scattering calculation
    character(len=128)                       :: error_msg

    error_msg = base_rte_lw_init(nlwstreams, nangles)
    if(len_trim(error_msg) == 0 .and. present(nlwstreams)) n_lw_streams = nlwstreams

  end function rte_lw_init
  ! --------------------------------------------------
  function rte_sw_init(nswstreams) result(error_msg)
    integer,           optional, intent( in) :: nswstreams
    character(len=128)                       :: error_msg

    error_msg = base_rte_sw_init(nswstreams)
    if(len_trim(error_msg) == 0 .and. present(nswstreams)) n_sw_streams = nswstreams
  end function rte_sw_init
  ! --------------------------------------------------
  !
  ! Interfaces using clear (gas + aerosol) and all-sky categories, starting from
  !   pressures, temperatures, and gas amounts for the gas contribution
  !
  ! --------------------------------------------------
  function rte_lw(k_dist, gas_concs, p_lay, t_lay, p_lev, &
                     t_sfc, sfc_emis, cloud_props,           &
                     allsky_fluxes, clrsky_fluxes,           &
                     aer_props, col_dry, t_lev, inc_flux) result(error_msg)
    type(ty_gas_optics_specification), intent(in   ) :: k_dist       !< derived type with spectral information
    type(ty_gas_concs),                intent(in   ) :: gas_concs    !< derived type encapsulating gas concentrations
    real(wp), dimension(:,:),          intent(in   ) :: p_lay, t_lay !< pressure [Pa], temperature [K] at layer centers (ncol,nlay)
    real(wp), dimension(:,:),          intent(in   ) :: p_lev        !< pressure at levels/interfaces [Pa] (ncol,nlay+1)
    real(wp), dimension(:),            intent(in   ) :: t_sfc     !< surface temperature           [K]  (ncol)
    real(wp), dimension(:,:),          intent(in   ) :: sfc_emis  !< emissivity at surface         []   (nband, ncol)
    class(ty_optical_props),           intent(in   ) :: cloud_props !< cloud optical properties (ncol,nlay,ngpt)
    class(ty_fluxes),                  intent(inout) :: allsky_fluxes, clrsky_fluxes

    ! Optional inputs
    class(ty_optical_props),  &
              optional,       intent(in ) :: aer_props   !< aerosol optical properties
    real(wp), dimension(:,:), &
              optional,       intent(in ) :: col_dry !< Molecular number density (ncol, nlay)
    real(wp), dimension(:,:), target, &
              optional,       intent(in ) :: t_lev     !< temperature at levels [K] (ncol, nlay+1)
    real(wp), dimension(:,:), target, &
              optional,       intent(in ) :: inc_flux   !< incident flux at domain top [W/m2] (ncol, ngpts)

    character(len=128)                    :: error_msg
    ! --------------------------------
    ! Local variables
    !
    class(ty_optical_props_arry),   allocatable :: optical_props

    real(wp), dimension(:,:,:),     allocatable :: lay_src, lev_src_inc, lev_src_dec
    real(wp), dimension(:,:),       allocatable :: sfc_src
    integer :: ncol, nlay, ngpt, nband
    logical :: top_at_1
    ! --------------------------------
    ! Problem sizes
    !
    error_msg = ""

    ncol  = size(p_lay, 1)
    nlay  = size(p_lay, 2)
    ngpt  = k_dist%get_ngpt()
    nband = k_dist%get_nband()

    top_at_1 = p_lay(1, 1) < p_lay(1, nlay)

    ! ------------------------------------------------------------------------------------
    !  Error checking
    !
    if(present(aer_props)) then
      if(any([aer_props%get_ncol(), &
              aer_props%get_nlay()] /= [ncol, nlay])) &
        error_msg = "rrtmpg_lw: aerosol properties inconsistently sized"
      if(.not. any(aer_props%get_ngpt() /= [ngpt, nband])) &
        error_msg = "rrtmpg_lw: aerosol properties inconsistently sized"
    end if

    if(present(t_lev)) then
      if(any([size(t_lev, 1), &
              size(t_lev, 2)] /= [ncol, nlay+1])) &
        error_msg = "rrtmpg_lw: t_lev inconsistently sized"
    end if

    if(present(inc_flux)) then
      if(any([size(inc_flux, 1), &
              size(inc_flux, 2)] /= [ncol, ngpt])) &
        error_msg = "rrtmpg_lw: incident flux inconsistently sized"
    end if
    if(len_trim(error_msg) > 0) return

    ! ------------------------------------------------------------------------------------
    ! Optical properties arrays
    !
    select case(n_lw_streams)
      case(0)
        allocate(ty_optical_props_1scl::optical_props)
      case(2)
        allocate(ty_optical_props_2str::optical_props)
      case default
        allocate(ty_optical_props_nstr::optical_props)
    end select

    select type (optical_props)
      class is (ty_optical_props_1scl) ! No scattering
        error_msg = optical_props%init_1scl(ncol, nlay, ngpt)
      class is (ty_optical_props_2str)
        error_msg = optical_props%init_2str(ncol, nlay, ngpt)
      class is (ty_optical_props_nstr)
        error_msg = optical_props%init_nstr(n_lw_streams/2, ncol, nlay, ngpt)
    end select
    if (error_msg /= '') return

    allocate(lay_src    (ncol, nlay, ngpt), &
             lev_src_inc(ncol, nlay, ngpt), &
             lev_src_dec(ncol, nlay, ngpt), &
             sfc_src    (ncol,       ngpt))
    ! ------------------------------------------------------------------------------------
    ! Clear skies
    !
    ! Gas optical depth -- pressure need to be expressed as Pa
    !
    error_msg = k_dist%gas_optics(p_lay, p_lev, t_lay, t_sfc, gas_concs,  &
                                  optical_props,                                          &
                                  lay_src, lev_src_inc, lev_src_dec, sfc_src,             &
                                  col_dry, t_lev)
    if (error_msg /= '') return
    ! ----------------------------------------------------
    ! Clear sky is gases + aerosols (if they're supplied)
    !
    if(present(aer_props)) then
      if(aer_props%get_ngpt() == ngpt) then
        error_msg = aer_props%increment(optical_props)
      else
        error_msg = aer_props%increment(optical_props, k_dist%get_band_lims_gpoint())
      end if
    end if
    if(error_msg /= '') return

    error_msg = base_rte_lw(optical_props, top_at_1, k_dist, &
                               lay_src, lev_src_inc, lev_src_dec, sfc_emis, sfc_src, &
                               clrsky_fluxes,                   &
                               inc_flux)
    if(error_msg /= '') return
    ! ------------------------------------------------------------------------------------
    ! All-sky fluxes = clear skies + clouds
    !
    if(cloud_props%get_ngpt() == ngpt) then
      error_msg = cloud_props%increment(optical_props)
    else
      error_msg = cloud_props%increment(optical_props, k_dist%get_band_lims_gpoint())
    end if
    if(error_msg /= '') return

    error_msg = base_rte_lw(optical_props, top_at_1, k_dist, &
                               lay_src, lev_src_inc, lev_src_dec, sfc_emis, sfc_src, &
                               allsky_fluxes,                   &
                               inc_flux)

  end function rte_lw
  ! --------------------------------------------------
  ! --------------------------------------------------
  function rte_sw(k_dist, gas_concs, p_lay, t_lay, p_lev, &
                                 mu0, sfc_alb_dir, sfc_alb_dif, cloud_props, &
                                 allsky_fluxes, clrsky_fluxes,           &
                                 aer_props, col_dry, inc_flux, tsi_scaling) result(error_msg)
    type(ty_gas_optics_specification), intent(in   ) :: k_dist       !< derived type with spectral information
    type(ty_gas_concs),                intent(in   ) :: gas_concs    !< derived type encapsulating gas concentrations
    real(wp), dimension(:,:),          intent(in   ) :: p_lay, t_lay !< pressure [Pa], temperature [K] at layer centers (ncol,nlay)
    real(wp), dimension(:,:),          intent(in   ) :: p_lev        !< pressure at levels/interfaces [Pa] (ncol,nlay+1)
    real(wp), dimension(:  ),          intent(in   ) :: mu0          !< cosine of solar zenith angle
    real(wp), dimension(:,:),          intent(in   ) :: sfc_alb_dir, sfc_alb_dif
                                                        !  surface albedo for direct and diffuse radiation (band, col)
    class(ty_optical_props),           intent(in   ) :: cloud_props !< cloud optical properties (ncol,nlay,ngpt)
    class(ty_fluxes),                  intent(inout) :: allsky_fluxes, clrsky_fluxes

    ! Optional inputs
    class(ty_optical_props),   target, &
              optional,       intent(in ) :: aer_props   !< aerosol optical properties
    real(wp), dimension(:,:), &
              optional,       intent(in ) :: col_dry, &  !< Molecular number density (ncol, nlay)
                                             inc_flux    !< incident flux at domain top [W/m2] (ncol, ngpts)
    real(wp), optional,       intent(in ) :: tsi_scaling !< Optional scaling for total solar irradiance

    character(len=128)                    :: error_msg
    ! --------------------------------
    ! Local variables
    !
    class(ty_optical_props_arry), allocatable :: optical_props
    real(wp), dimension(:,:),     allocatable :: toa_flux
    integer :: ncol, nlay, ngpt, nband
    logical :: top_at_1
    ! --------------------------------
    ! Problem sizes
    !
    error_msg = ""

    ncol  = size(p_lay, 1)
    nlay  = size(p_lay, 2)
    ngpt  = k_dist%get_ngpt()
    nband = k_dist%get_nband()

    top_at_1 = p_lay(1, 1) < p_lay(1, nlay)

    ! ------------------------------------------------------------------------------------
    !  Error checking
    !
    if(present(inc_flux) .and. present(tsi_scaling)) &
      error_msg = "rrtmpg_sw: only one of inc_flux, tsi_scaling may be supplied."

    if(present(aer_props)) then
      if(any([aer_props%get_ncol(), &
              aer_props%get_nlay()] /= [ncol, nlay])) &
        error_msg = "rrtmpg_sw: aerosol properties inconsistently sized"
      if(.not. any(aer_props%get_ngpt() /= [ngpt, nband])) &
        error_msg = "rrtmpg_sw: aerosol properties inconsistently sized"
    end if

    if(present(tsi_scaling)) then
      if(tsi_scaling <= 0._wp) &
        error_msg = "rrtmpg_sw: tsi_scaling is < 0"
    end if

    if(present(inc_flux)) then
      if(any([size(inc_flux, 1), &
              size(inc_flux, 2)] /= [ncol, ngpt])) &
        error_msg = "rrtmpg_sw: incident flux inconsistently sized"
    end if
    if(len_trim(error_msg) > 0) return

    ! ------------------------------------------------------------------------------------
    ! Optical properties arrays
    !
    select case(n_sw_streams)
      case(0)
        allocate(ty_optical_props_1scl::optical_props)
      case(2)
        allocate(ty_optical_props_2str::optical_props)
      case default
        allocate(ty_optical_props_nstr::optical_props)
    end select

    select type (optical_props)
      class is (ty_optical_props_1scl) ! No scattering
        error_msg = optical_props%init_1scl(ncol, nlay, ngpt)
      class is (ty_optical_props_2str)
        error_msg = optical_props%init_2str(ncol, nlay, ngpt)
      class is (ty_optical_props_nstr)
        error_msg = optical_props%init_nstr(n_sw_streams/2, ncol, nlay, ngpt)
    end select
    if (error_msg /= '') return

    allocate(toa_flux(ncol, ngpt))
    ! ------------------------------------------------------------------------------------
    ! Clear skies
    !
    ! Gas optical depth -- pressure need to be expressed as Pa
    !
    error_msg = k_dist%gas_optics(p_lay, p_lev, t_lay, gas_concs,  &
                                  optical_props, toa_flux,                          &
                                  col_dry)
    if (error_msg /= '') return
    !
    ! If users have supplied an incident flux, use that
    !
    if(present(inc_flux))    toa_flux(:,:) = inc_flux(:,:)
    if(present(tsi_scaling)) toa_flux(:,:) = toa_flux(:,:) * tsi_scaling
    ! ----------------------------------------------------
    ! Clear sky is gases + aerosols (if they're supplied)
    !
    if(present(aer_props)) then
      if(aer_props%get_ngpt() == ngpt) then
        error_msg = aer_props%increment(optical_props)
      else
        error_msg = aer_props%increment(optical_props, k_dist%get_band_lims_gpoint())
      end if
    end if
    if(error_msg /= '') return

    error_msg = base_rte_sw(optical_props, top_at_1, k_dist, &
                               mu0, toa_flux,                   &
                               sfc_alb_dir, sfc_alb_dif,        &
                               clrsky_fluxes)

    if(error_msg /= '') return
    ! ------------------------------------------------------------------------------------
    ! All-sky fluxes = clear skies + clouds
    !
    if(cloud_props%get_ngpt() == ngpt) then
      error_msg = cloud_props%increment(optical_props)
    else
      error_msg = cloud_props%increment(optical_props, k_dist%get_band_lims_gpoint())
    end if
    if(error_msg /= '') return

    error_msg = base_rte_sw(optical_props, top_at_1, k_dist, &
                               mu0, toa_flux,                   &
                               sfc_alb_dir, sfc_alb_dif,        &
                               allsky_fluxes)

  end function rte_sw

end module mo_rrtmgp_clr_all_sky
