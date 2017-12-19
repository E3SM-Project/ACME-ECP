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
  use mo_rte_kind,           only: wp, wl
  use mo_optical_props,         only: ty_optical_props_arry, ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
  use mo_rte_solver_kernels, only: lw_solver_noscat, lw_solver_noscat_quad
  implicit none
  private

  public :: lw_solver_init, lw_solver
  ! -------------------------------------------------------------------------------------
  integer,    private :: n_quad_angs = 1
  integer,  parameter :: max_gauss_pts = 4
  real(wp), parameter,                         &  
    dimension(max_gauss_pts, max_gauss_pts) :: & 
      !
      ! Weights and angle secants for first order (k=1) Gaussian quadrature. 
      !   Values from Table 2, Clough et al, 1992, doi:10.1029/92JD01419
      !   after Abramowitz & Stegun 1972, page 921 
      !
      gauss_Ds  = RESHAPE([1.66_wp,               0._wp,         0._wp,         0._wp, &  ! Diffusivity angle, not Gaussian angle  
                           1.18350343_wp, 2.81649655_wp,         0._wp,         0._wp, & 
                           1.09719858_wp, 1.69338507_wp, 4.70941630_wp,         0._wp, & 
                           1.06056257_wp, 1.38282560_wp, 2.40148179_wp, 7.15513024_wp], & 
                          [max_gauss_pts, max_gauss_pts]),              &
      gauss_wts = RESHAPE([0.5_wp,          0._wp,           0._wp,           0._wp, & 
                           0.3180413817_wp, 0.1819586183_wp, 0._wp,           0._wp, & 
                           0.2009319137_wp, 0.2292411064_wp, 0.0698269799_wp, 0._wp, & 
                           0.1355069134_wp, 0.2034645680_wp, 0.1298475476_wp, 0.0311809710_wp], & 
                           [max_gauss_pts, max_gauss_pts])
contains
  ! ---------------------------------------------------------------
  function lw_solver_init(n_angles) result(error_msg) 
    integer, intent(in) :: n_angles
    character(len=128)  :: error_msg
    
    error_msg = "" 
    if(n_angles <= 0)            error_msg = "lw_solver_init: n_angles provided is less than 0"
    if(n_angles > max_gauss_pts) error_msg = "lw_solver_init: n_angles provided greater than max allowed"
    if(error_msg == "") n_quad_angs = n_angles
   
  end function lw_solver_init
  ! ---------------------------------------------------------------
  !
  ! Solution to the radiative transfer equation assuming internal emission
  !
  function lw_solver(ncol, nlay, ngpt, top_is_1, &
                     atmos, lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, & 
                     flux_up, flux_dn,           & 
                     Ds) result(error_msg) 
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
    real(wp), dimension(ncol,       ngpt), optional, &                     ! "User"-supplied integration secants    
                                           intent( in) :: Ds               ! One per column / g-point 
                                                                               
    character(len=128)                                 :: error_msg
    ! --------------------------------------------------
    error_msg = ""
    select type (atmos)
      class is (ty_optical_props_1scl)
        !
        ! No scattering two-stream calculation
        !
        if(present(Ds)) then 
          !
          ! "Users" have supplied quadrature angle cosines per column, g-point 
          !  
          if(any(Ds < 1._wp)) then 
            error_msg = "Supplied quadrature secants out of bounds" 
            return 
          end if 
          call lw_solver_noscat(ncol, nlay, ngpt, logical(top_is_1, wl), & 
                                Ds, gauss_wts(1,1),                      &
                                atmos%tau, lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, &
                                flux_up, flux_dn)
        else 
          !
          ! First-order Gaussian quadrature  
          !  
          call lw_solver_noscat_quad(ncol, nlay, ngpt, logical(top_is_1, wl), & 
                                n_quad_angs, gauss_Ds(1:n_quad_angs,n_quad_angs), gauss_wts(1:n_quad_angs,n_quad_angs), & 
                                atmos%tau, lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src,  & 
                                flux_up, flux_dn)

        end if
      class is (ty_optical_props_2str)
        !
        ! two-stream calculation with scattering
        !
        error_msg = 'lw_solver(...ty_optical_props_2str...) not yet implemented'
      class is (ty_optical_props_nstr)
        !
        ! n-stream calculation
        !
        error_msg = 'lw_solver(...ty_optical_props_nstr...) not yet implemented'
    end select

  end function lw_solver
  ! ---------------------------------------------------------------
end module mo_lw_solver
