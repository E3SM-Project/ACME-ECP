! This code is part of Radiative Transfer for Energetics (RTE)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2015-2016,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!
! Description:  RT solver for rte_sw. Provides a no-scattering (direct-beam only)
!   and two-stream/adding calculations.

module mo_sw_solver
  use mo_rte_kind,           only: wp, wl
  use mo_optical_props,      only: ty_optical_props_arry, &
                                   ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
  use mo_rte_solver_kernels, only: apply_BC, sw_solver_noscat, sw_solver_2stream
  implicit none
  private

  public :: sw_solver

contains
  ! ---------------------------------------------------------------
  !
  ! Solution to the radiative transfer equation assuming external direction and possibly diffuse sources
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
    real(wp), dimension(ncol,ngpt), intent( in) :: inc_flux          !< direct beam incident flux at top-of-atmosphere [W/m2]
    real(wp), dimension(ncol,nlay+1,ngpt), &
                                     intent(out) :: flux_up, flux_dn, &  ! Fluxes [W/m2]
                                                    flux_dir             ! Downward direct
                                                                         ! Top level (= merge(1, nlay+1, top_is_1)
                                                                         ! must contain incident flux boundary condition
    real(wp), dimension(ncol,ngpt), optional, &
                                     intent( in) :: inc_flux_dif     !< diffuse incident flux at top-of-atmosphere [W/m2]
    character(len=128)                           :: sw_solver

    ! --------------------------------------------------
    sw_solver = ""
    !
    ! Apply boundary conditions
    !   On input flux_dn is the diffuse component; the last action in each solver is to add
    !   direct and diffuse to represent the total, consistent with the LW
    !
    call apply_BC(ncol, nlay, ngpt, logical(top_is_1, wl),   inc_flux, mu0, flux_dir)
    if(present(inc_flux_dif)) then
      call apply_BC(ncol, nlay, ngpt, logical(top_is_1, wl), inc_flux_dif,  flux_dn )
    else
      call apply_BC(ncol, nlay, ngpt, logical(top_is_1, wl),                flux_dn )
    end if

    select type (atmos)
      class is (ty_optical_props_1scl)
        !
        ! Direct beam only
        !
        call sw_solver_noscat(ncol, nlay, ngpt, logical(top_is_1, wl), &
                              atmos%tau, mu0,                          &
                              flux_dir)
        ! No diffuse flux
        flux_up = 0._wp
        flux_dn = 0._wp
        sw_solver = ""
      class is (ty_optical_props_2str)
        !
        ! two-stream calculation with scattering
        !
        call sw_solver_2stream(ncol, nlay, ngpt, logical(top_is_1, wl), &
                               atmos%tau, atmos%ssa, atmos%g, mu0, &
                               sfc_alb_dir, sfc_alb_dif,           &
                               flux_up, flux_dn, flux_dir)
        sw_solver = ""
      class is (ty_optical_props_nstr)
        !
        ! n-stream calculation
        !
        ! not yet implemented so fail
        !
        sw_solver = 'sw_solver(...ty_optical_props_nstr...) not yet implemented'
    end select
  end function sw_solver
  ! ---------------------------------------------------------------
end module mo_sw_solver
