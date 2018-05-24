! Module: mo_lw_solver

! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2015,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!
! Description:  RT solver for rte_lw.  Performs 2-stream no-scattering, (future) 2-stream scattering, and
! (future) n-stream scattering.

module mo_lw_solver
  !
  ! ECRAD modules
  !
  use parkind1,                 only : jprb ! Working precision
  use radiation_config, only         : config_type
  use radiation_single_level, only   : single_level_type
  use radiation_thermodynamics, only : thermodynamics_type
  use radiation_cloud, only          : cloud_type
  use radiation_flux, only           : flux_type
  use radiation_homogeneous_lw, only : solver_homogeneous_lw
  !
  ! RTE modules.
  !
  use mo_rte_kind,              only: wp, wl
  use mo_optical_props,         only: ty_optical_props_arry, ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
  implicit none
  private

  !
  ! RTE
  !
  real(wp), parameter :: pi = acos(-1._wp)
  public :: lw_solver_init, lw_solver
contains
  ! ---------------------------------------------------------------
  function lw_solver_init(n_angles) result(error_msg)
    integer, intent(in) :: n_angles
    character(len=128)  :: error_msg

    error_msg = ""

  end function lw_solver_init
  ! ---------------------------------------------------------------
  !
  ! Solution to the radiative transfer equation assuming internal emission
  !
  function lw_solver(ncol, nlay, ngpt, top_is_1, &
                     atmos, lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, &
                     flux_up, flux_dn,           &
                     inc_flux, Ds) result(error_msg)
    integer,                               intent( in) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
    logical,                               intent( in) :: top_is_1       ! True if arrays are indexed top to bottom.
    class(ty_optical_props_arry),          intent( in) :: atmos          ! Optical properties of the atmosphere
    real(wp), dimension(ncol,nlay,  ngpt), intent( in) :: lay_source     ! Planck source at layer average temperature
                                                                         ! [W/m2] (ncol, nlay, ngpt)
    real(wp), dimension(ncol,nlay,ngpt), intent( in) :: lev_source_inc, lev_source_dec
                                      ! Planck source at layer edge for radiation in increasing/decreasing ilay direction
                                      ! Includes spectral weighting that accounts for state-dependent frequency to g-space mapping
                                      ! lev_source_dec applies the mapping in layer i to the Planck function at layer i
                                      ! lev_source_inc applies the mapping in layer i to the Planck function at layer i+1
                                      ! [W/m2] (ncol, nlay, ngpt)
    real(wp), dimension(ncol,       ngpt), intent( in) :: sfc_emis         ! Surface emissivity      []
    real(wp), dimension(ncol,       ngpt), intent( in) :: sfc_src          ! Surface source function [W/m2]
    real(wp), dimension(ncol,nlay+1,ngpt), intent(out) :: flux_up, flux_dn ! Fluxes [W/m2]
                                                                           ! Top level (= merge(1, nlay+1, top_is_1)
                                                                           ! must contain incident flux boundary condition
    real(wp), dimension(ncol,       ngpt), optional, &                     ! incident diffuse flux [W/m2]
                                           intent( in) :: inc_flux
    real(wp), dimension(ncol,       ngpt), optional, &                     ! "User"-supplied integration secants
                                           intent( in) :: Ds               ! One per column / g-point
    character(len=128)                                 :: error_msg

    ! ----------------------------------------------------------------------------
    ! Derived types for the inputs to ECRAD
    type(config_type)         :: config
    type(single_level_type)   :: single_level   ! need lw_emissivity, sw_albedo, sw_albedo_direct, cos_sza
    type(thermodynamics_type) :: thermodynamics ! not used
    type(cloud_type)          :: cloud          ! need only cloud%fraction
    ! Derived type containing outputs from the radiation scheme
    type(flux_type)           :: flux           !

    real(jprb), dimension(:,:,:), allocatable :: od, ssa, g, od_cloud, ssa_cloud, g_cloud, planck_hl
    real(jprb), dimension(:,:  ), allocatable :: planck_surf
    ! ----------------------------------------------------------------------------
    integer               :: istartcol, iendcol ! range of columns to process
    integer               :: icol, ilay, igpt
    real(wp), parameter :: pi = acos(-1._wp)
    real(wp), dimension(:,:), allocatable :: lev_source
    ! ----------------------------------------------------------------------------
    error_msg = ""
    if(.not. top_is_1) then
      error_msg = "lw_solver: atmosphere has to be ordered top to bottom (top_is_1 = .false.) for ECRAD"
      return
    end if
    !
    ! Source functions
    !
    allocate(planck_surf(ngpt,ncol), planck_hl(ngpt,nlay+1,ncol), lev_source(ncol, nlay+1))
    do igpt = 1, ngpt
      call lw_combine_sources(ncol, nlay, top_is_1, &
                              lev_source_inc(:,:,igpt), lev_source_dec(:,:,igpt), &
                              lev_source)
      do icol = 1, ncol
        do ilay = 1, nlay + 1
          planck_hl(igpt, ilay, icol) = lev_source(icol, ilay) * pi
        end do
      end do
    end do
    planck_surf = transpose(sfc_src) * pi

    istartcol = 1; iendcol = ncol
    !
    ! ECRAD types
    ! thermo type is unused

    ! cloud
    allocate(cloud%fraction(ncol,nlay))
    cloud%fraction(1:ncol,1:nlay) = 0._wp

    ! single-level
    allocate(single_level%lw_emissivity(ncol,ngpt))
    single_level%lw_emissivity = sfc_emis(1:ncol,1:ngpt)

    ! config type
    config%do_clear = .true.
    config%do_sw    = .false.
    config%do_lw    = .true.
    config%do_lw_derivatives = .false.
    config%n_bands_lw           = ngpt ! Would normally be the number of bands but we don't have access
    config%n_g_lw               = ngpt
    config%n_g_lw_if_scattering = ngpt
    config%i_band_from_reordered_g_lw = [(igpt, igpt = 1, ngpt)]
    config%i_emiss_from_band_lw       = [(igpt, igpt = 1, ngpt)]

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
              od (igpt, ilay, icol) = atmos%tau(icol, ilay, igpt)
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

    call solver_homogeneous_lw(nlay,istartcol,iendcol, &
        &  config, single_level, thermodynamics, cloud, &
        &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, planck_hl, planck_surf, &
        &  flux)
    !
    ! Homogeneous solver reports only broadband fluxes summed across g-points
    !   Broadband sums will be correct below; band fluxes will be all messed up
    !
    flux_up(1:ncol,1:nlay+1,1) = flux%lw_up; flux_up(1:ncol,1:nlay+1,2:ngpt) = 0._wp
    flux_dn(1:ncol,1:nlay+1,1) = flux%lw_dn; flux_dn(1:ncol,1:nlay+1,2:ngpt) = 0._wp

   end function lw_solver

   subroutine lw_combine_sources(ncol, nlay, top_is_1, &
                                 lev_src_inc, lev_src_dec, lev_source) bind(C)
     integer,                           intent(in ) :: ncol, nlay
     logical,                           intent(in ) :: top_is_1
     real(wp), dimension(ncol, nlay  ), intent(in ) :: lev_src_inc, lev_src_dec
     real(wp), dimension(ncol, nlay+1), intent(out) :: lev_source

     integer :: icol, ilay
     ! ---------------------------------------------------------------
     ilay = 1
     do icol = 1,ncol
       lev_source(icol, ilay) =        lev_src_dec(icol, ilay)
     end do
     do ilay = 2, nlay
       do icol = 1,ncol
         lev_source(icol, ilay) = sqrt(lev_src_dec(icol, ilay) * &
                                       lev_src_inc(icol, ilay-1))
       end do
     end do
     ilay = nlay+1
     do icol = 1,ncol
       lev_source(icol, ilay) =        lev_src_inc(icol, ilay-1)
     end do

   end subroutine lw_combine_sources

end module mo_lw_solver
