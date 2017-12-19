! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2016-2017,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!
! Description:  Numeric calculations for radiative transfer solvers.

module mo_rte_solver_kernels
  use,  intrinsic :: iso_c_binding
  use mo_rte_kind, only: wp, wl
  implicit none
  private
  public :: lw_solver_noscat,  lw_solver_noscat_quad, &
            sw_solver_noscat, sw_solver_2stream

  ! These routines don't really need to be visible but making them so is useful for testing.
  public :: two_stream, &
            adding_sw

  real(wp), parameter :: pi = acos(-1._wp)
contains
  ! ---------------------------------------------------------------
  !
  ! Description: LW transport, no scattering, multi-angle quadrature
  !   Users provide a set of weights and quadrature angles
  !
  ! ---------------------------------------------------------------
  subroutine lw_solver_noscat_quad(ncol, nlay, ngpt, top_is_1, nmus, Ds, weights, &
                                   tau, lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, flux_up, flux_dn) &
                                   bind (C)
    integer,                    intent( in) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
    logical(wl),                intent( in) :: top_is_1
    integer,                    intent( in) :: nmus          ! number of quadrature angles
    real(wp), dimension(nmus),  intent( in) :: Ds, weights  ! quadrature secants, weights
    real(wp), dimension(ncol,nlay,  ngpt), intent( in) :: tau          ! Absorption optical thickness []
    real(wp), dimension(ncol,nlay,  ngpt), intent( in) :: lay_source   ! Planck source at layer average temperature [W/m2]
    real(wp), dimension(ncol,nlay,  ngpt), intent( in) :: lev_source_inc
                                        ! Planck source at layer edge for radiation in increasing ilay direction [W/m2]
                                        ! Includes spectral weighting that accounts for state-dependent frequency to g-space mapping
    real(wp), dimension(ncol,nlay,  ngpt), intent( in) :: lev_source_dec
                                               ! Planck source at layer edge for radiation in decreasing ilay direction [W/m2]
    real(wp), dimension(ncol,       ngpt), intent( in) :: sfc_emis         ! Surface emissivity      []
    real(wp), dimension(ncol,       ngpt), intent( in) :: sfc_src          ! Surface source function [W/m2]
    real(wp), dimension(ncol,nlay+1,ngpt), intent(out) :: flux_up, flux_dn ! Radiances [W/m2-str]
                                                                           ! Top level must contain incident flux boundary condition
    ! Local variables
    real(wp), dimension(ncol,nlay+1,ngpt) :: radn_dn, radn_up ! Fluxes per quad angle
    real(wp), dimension(ncol,       ngpt) :: Ds_ncol

    integer :: imu, top_level
    ! ------------------------------------
    Ds_ncol(:,:) = Ds(1)
    call lw_solver_noscat(ncol, nlay, ngpt, &
                          top_is_1, Ds_ncol, weights(1), tau, &
                          lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, &
                          flux_up, flux_dn)

    top_level = MERGE(1, nlay+1, top_is_1)
    radn_dn(:,top_level,:) = flux_dn(:, top_level, :) ! Flux boundary condition

    do imu = 2, nmus
      Ds_ncol(:,:) = Ds(imu)
      call lw_solver_noscat(ncol, nlay, ngpt, &
                            top_is_1, Ds_ncol, weights(imu), tau, &
                            lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, &
                            radn_up, radn_dn)
      flux_up(:,:,:) = flux_up(:,:,:) + radn_up(:,:,:)
      flux_dn(:,:,:) = flux_dn(:,:,:) + radn_dn(:,:,:)
    end do
  end subroutine lw_solver_noscat_quad
  ! ---------------------------------------------------------------
  !
  ! Description: LW transport, no scattering, mu (cosine of integration angle) specified by column
  !   Does radiation calculation at user-supplied angles; converts radiances to flux
  !   using user-supplied weights
  !
  ! ---------------------------------------------------------------
  subroutine lw_solver_noscat(ncol, nlay, ngpt, top_is_1, D, weight,                             &
                              tau, lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, &
                              radn_up, radn_dn) bind (C)
    integer,                    intent( in) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
    logical(wl),                intent( in) :: top_is_1
    real(wp), dimension(ncol,       ngpt), intent( in) :: D            ! secant of propagation angle  []
    real(wp),                              intent( in) :: weight       ! quadrature weight
    real(wp), dimension(ncol,nlay,  ngpt), intent( in) :: tau          ! Absorption optical thickness []
    real(wp), dimension(ncol,nlay,  ngpt), intent( in) :: lay_source   ! Planck source at layer average temperature [W/m2]
    real(wp), dimension(ncol,nlay,  ngpt), target, &
                                           intent( in) :: lev_source_inc, lev_source_dec
                                                                       ! Planck source at layer edge for radiation in increasing/decreasing ilay direction
                                                                       ! lev_source_dec applies the mapping in layer i to the Planck function at layer i
                                                                       ! lev_source_inc applies the mapping in layer i to the Planck function at layer i+1
    real(wp), dimension(ncol,       ngpt), intent( in) :: sfc_emis         ! Surface emissivity      []
    real(wp), dimension(ncol,       ngpt), intent( in) :: sfc_src          ! Surface source function [W/m2]
    real(wp), dimension(ncol,nlay+1,ngpt), intent(out) :: radn_up, radn_dn ! Radiances [W/m2-str]
                                                                           ! Top level must contain incident flux boundary condition

    ! Local variables
    integer :: sfc_lev
    real(wp), dimension(ncol,nlay  ) :: tau_loc, &  ! path length (tau/mu)
                                        trans       ! transmissivity  = exp(-tau)
    real(wp), dimension(ncol,nlay  ) :: source_dn, source_up

    real(wp), dimension(:,:,:), pointer :: lev_source_up, lev_source_dn ! Mapping increasing/decreasing indicies to up/down

    integer :: ilev, igpt, top_level
    ! ------------------------------------
    ! Which way is up?
    ! Level Planck sources for upward and downward radiation
    ! When top_is_1, lev_source_up => lev_source_dec
    !                lev_source_dn => lev_source_inc, and vice-versa

    ! Transport is for intensity
    !   convert flux at top of domain to intensity assuming azimuthal isotropy
    !
    top_level = MERGE(1, nlay+1, top_is_1)
    if(top_is_1) then
      lev_source_up => lev_source_dec
      lev_source_dn => lev_source_inc
    else
      lev_source_up => lev_source_inc
      lev_source_dn => lev_source_dec
    end if

    radn_dn(:,top_level,:) = radn_dn(:,top_level,:)/(2._wp * pi * weight)
    do igpt = 1, ngpt
      do ilev = 1, nlay
        tau_loc(:,ilev) = tau(:,ilev,igpt)*D(:,igpt)
        trans  (:,ilev) = exp(-tau_loc(:,ilev))
      end do
      !
      ! Compute radiation emitted by each layer in each direction
      !
      call lw_source_noscat(ncol, nlay, &
                            lay_source(:,:,igpt), lev_source_up(:,:,igpt), lev_source_dn(:,:,igpt), &
                            tau_loc, trans, source_dn, source_up)

      ! Indexing into arrays for upward and downward propagation depends on the vertical
      !   orientation of the arrays (whether the domain top is at the first or last index)
      ! There are explicit routines for bottom-up and top-up transport so compilers will have no trouble optimizing them.
      if(top_is_1) then
        call lw_transport_noscat_top1(ncol, nlay,  &
                                 tau_loc, trans, source_dn, source_up, &
                                 sfc_emis(:,igpt), sfc_src(:,igpt), radn_up(:,:,igpt), radn_dn(:,:,igpt))
      else
        call lw_transport_noscat_topN(ncol, nlay,  &
                                 tau_loc, trans, source_dn, source_up, &
                                 sfc_emis(:,igpt), sfc_src(:,igpt), radn_up(:,:,igpt), radn_dn(:,:,igpt))
      end if
      !
      ! Convert intensity to flux assuming azimuthal isotropy and quadrature weight
      !
      radn_dn(:,:,igpt) = 2._wp * pi * weight * radn_dn(:,:,igpt)
      radn_up(:,:,igpt) = 2._wp * pi * weight * radn_up(:,:,igpt)
    end do  ! g point loop

  end subroutine lw_solver_noscat
! ---------------------------------------------------------------
  !
  ! Compute LW source function for upward and downward emission at levels using linear-in-tau assumption
  ! See Clough et al., 1992, doi: 10.1029/92JD01419, Eq 13
  !
  ! ---------------------------------------------------------------
  subroutine lw_source_noscat(ncol, nlay, lay_source, lev_source_up, lev_source_dn, tau, trans, source_dn, source_up) &
                            bind (C)
    integer,                         intent(in) :: ncol, nlay
    real(wp), dimension(ncol, nlay), intent(in) :: lay_source, & ! Planck source at layer center
                                                   lev_source_up, & ! Planck source at levels (layer edges),
                                                   lev_source_dn, & !   increasing/decreasing layer index
                                                   tau,        & ! Optical path (tau/mu)
                                                   trans         ! Transmissivity (exp(-tau))

    real(wp), dimension(ncol, nlay), intent(out):: source_dn, source_up
                                                                   ! Source function at layer edges
                                                                   ! Down at the bottom of the layer, up at the top

    integer             :: icol, ilay
    real(wp)            :: fact
    real(wp), parameter :: tau_thresh = sqrt(epsilon(tau)), one_third = 1._wp/3._wp
    ! ---------------------------------------------------------------
    do ilay = 1, nlay
      do icol = 1, ncol
      !
      ! Weighting factor. Use 2nd order series expansion when rounding error (~tau^2)
      !   is of order epsilon (smallest difference from 1. in working precision)
      !   Thanks to Peter Blossey
      !
      fact = merge((1._wp - trans(icol,ilay))/tau(icol,ilay) - trans(icol,ilay), &
                   tau(icol, ilay) * ( 0.5_wp -  tau(icol, ilay) * one_third), &
                   tau(icol, ilay) > tau_thresh)
      !
      ! Equation below is developed in Clough et al., 1992, doi:10.1029/92JD01419, Eq 13
      !
      source_dn(icol,ilay) = (1._wp - trans(icol,ilay)) * lev_source_dn(icol,ilay) + &
                              2._wp * fact * (lay_source(icol,ilay) - lev_source_dn(icol,ilay))
      source_up(icol,ilay) = (1._wp - trans(icol,ilay)) * lev_source_up(icol,ilay  ) + &
                              2._wp * fact * (lay_source(icol,ilay) - lev_source_up(icol,ilay))
      end do
    end do
  end subroutine lw_source_noscat
  ! ---------------------------------------------------------------
  !
  ! Longwave no-scattering transport when vertical index 1 is the top of the domain
  !
  ! ---------------------------------------------------------------
  subroutine lw_transport_noscat_top1(ncol, nlay, &
                                      tau, trans, source_dn, source_up, sfc_emis, sfc_src, &
                                      radn_up, radn_dn) bind (C)
    integer,                          intent( in) :: ncol, nlay ! Number of columns, layers, g-points
    real(wp), dimension(ncol,nlay  ), intent( in) :: tau, &     ! Absorption optical thickness, pre-divided by mu []
                                                     trans      ! transmissivity = exp(-tau)
    real(wp), dimension(ncol,nlay  ), intent( in) :: source_dn, &
                                                     source_up  ! Diffuse radiation emitted by the layer
    real(wp), dimension(ncol       ), intent( in) :: sfc_emis         ! Surface emissivity      []
    real(wp), dimension(ncol       ), intent( in) :: sfc_src          ! Surface source function [W/m2]
    real(wp), dimension(ncol,nlay+1), intent(out) :: radn_up, radn_dn ! Radiances [W/m2-str]
                                                                      ! Top level must contain incident flux boundary condition
    ! Local variables
    integer :: ilev
    ! ---------------------------------------------------
    ! Downward propagation
    ! For the flux at this level, what was the previous level, and which layer has the
    !   radiation just passed through?
    ! layer index = level index - 1
    ! previous level is up (-1)
    do ilev = 2, nlay+1
      radn_dn(:,ilev) = trans(:,ilev-1)*radn_dn(:,ilev-1) + source_dn(:,ilev-1)
    end do

    ! Surface reflection and emission
    radn_up(:,nlay+1) = radn_dn(:,nlay+1)*(1._wp-sfc_emis(:)) + sfc_src(:)*sfc_emis(:)

    ! Upward propagation
    ! layer index = level index
    ! previous level is down (+1)
    do ilev = nlay, 1, -1
      radn_up(:,ilev) = trans(:,ilev  )*radn_up(:,ilev+1) + source_up(:,ilev)
    end do
  end subroutine lw_transport_noscat_top1
! ---------------------------------------------------------------
  !
  ! Longwave no-scattering transport when vertical index N is the top of the domain
  !
  subroutine lw_transport_noscat_topN(ncol, nlay,  &
                                      tau, trans, source_dn, source_up, sfc_emis, sfc_src, &
                                      radn_up, radn_dn) bind (C)
    integer,                          intent( in) :: ncol, nlay ! Number of columns, layers, g-points
    real(wp), dimension(ncol,nlay  ), intent( in) :: tau, &     ! Absorption optical thickness, pre-divided by mu []
                                                     trans      ! transmissivity = exp(-tau)
    real(wp), dimension(ncol,nlay  ), intent( in) :: source_dn, &
                                                     source_up  ! Diffuse radiation emitted by the layer
    real(wp), dimension(ncol       ), intent( in) :: sfc_emis         ! Surface emissivity      []
    real(wp), dimension(ncol       ), intent( in) :: sfc_src          ! Surface source function [W/m2]
    real(wp), dimension(ncol,nlay+1), intent(out) :: radn_up, radn_dn ! Radiances [W/m2-str]
                                                                      ! Top level must contain incident flux boundary condition
    ! Local variables
    integer :: ilev
    ! ---------------------------------------------------
    ! layer index = level index
    ! previous level is up (+1)
    do ilev = nlay, 1, -1
      radn_dn(:,ilev) = trans(:,ilev  ) * radn_dn(:,ilev+1) + source_dn(:,ilev)
    end do

    ! Surface reflection and emission
    radn_up(:,1) = radn_dn(:,1)*(1._wp - sfc_emis(:)) + sfc_src(:)*sfc_emis(:)

    ! Upward propagation
    ! layer index = level index - 1
    ! previous level is down (-1)
    do ilev = 2, nlay+1
      radn_up(:,ilev) = trans(:,ilev-1) * radn_up(:,ilev-1) +  source_up(:,ilev-1)
    end do
  end subroutine lw_transport_noscat_topN
! ---------------------------------------------------------------
!   Shortwave kernels
! ---------------------------------------------------------------
  pure subroutine sw_solver_noscat(ncol, nlay, ngpt, &
                              top_is_1, tau, mu0, flux_dir) bind (C)
    integer,                    intent( in) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
    logical(wl),                intent( in) :: top_is_1
    real(wp), dimension(ncol,nlay,  ngpt), intent( in) :: tau          ! Absorption optical thickness []
    real(wp), dimension(ncol            ), intent( in) :: mu0          ! cosine of solar zenith angle
    real(wp), dimension(ncol,nlay+1,ngpt), intent(out) :: flux_dir     ! Direct-beam flux, spectral [W/m2]
                                                                       ! Top level must contain incident flux boundary condition
    integer :: icol, ilev, igpt
    real(wp) :: mu0_inv(ncol)

    ! ------------------------------------
    mu0_inv(1:ncol) = 1._wp/mu0(1:ncol)
    ! Indexing into arrays for upward and downward propagation depends on the vertical
    !   orientation of the arrays (whether the domain top is at the first or last index)
    ! We write the loops out explicitly so compilers will have no trouble optimizing them.

    ! Downward propagation
    if(top_is_1) then
      ! For the flux at this level, what was the previous level, and which layer has the
      !   radiation just passed through?
      ! layer index = level index - 1
      ! previous level is up (-1)
      do igpt = 1, ngpt
        do ilev = 2, nlay+1
          flux_dir(:,ilev,igpt) = flux_dir(:,ilev-1,igpt) * exp(-tau(:,ilev,igpt)*mu0_inv(:))
        end do
      end do
    else
      ! layer index = level index
      ! previous level is up (+1)
      do igpt = 1, ngpt
        do ilev = nlay, 1, -1
          flux_dir(:,ilev,igpt) = flux_dir(:,ilev+1,igpt) * exp(-tau(:,ilev,igpt)*mu0_inv(:))
        end do
      end do
    end if
  end subroutine sw_solver_noscat
  ! ---------------------------------------------------------------
   subroutine sw_solver_2stream (ncol, nlay, ngpt, top_is_1, &
                                 tau, ssa, g, mu0,           &
                                 sfc_alb_dir, sfc_alb_dif,   &
                                 flux_up, flux_dn, flux_dir) bind (C)
    integer,                    intent( in) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
    logical(wl),                intent( in) :: top_is_1
    real(wp), dimension(ncol,nlay,  ngpt), intent( in) :: tau, &  ! Optical thickness,
                                                          ssa, &  ! single-scattering albedo,
                                                          g       ! asymmetry parameter []
    real(wp), dimension(ncol            ), intent( in) :: mu0     ! cosine of solar zenith angle
    real(wp), dimension(ncol,       ngpt), intent( in) :: sfc_alb_dir, sfc_alb_dif
                                                                  ! Spectral albedo of surface to direct and diffuse radiation
    real(wp), dimension(ncol,nlay+1,ngpt), &
                                           intent(out) :: flux_up, flux_dn, &  ! Fluxes [W/m2]
                                                          flux_dir             ! Downward direct
                                                                               ! Top level (= merge(1, nlay+1, top_is_1)
                                                                               ! must contain incident flux boundary condition
    integer :: igpt
    real(wp), dimension(ncol,nlay) :: Rdif, Tdif, Rdir, Tdir, Tnoscat
    ! ------------------------------------
    !
    ! It's not clear how to be most efficient here. The two-stream calculation is atomic
    !   but the transport (adding) is not. Combining the loops reduces the memory footprint
    !   for the transmission and reflection arrays by a factor of ngpt, normally a lot,
    !   but then the size for the two-stream problem is greatly reduced.
    !
    do igpt = 1, ngpt
      !
      ! Compute cell properties
      !
      call two_stream(ncol, nlay, mu0,                                &
                      tau (:,:,igpt), ssa (:,:,igpt), g   (:,:,igpt), &
                      Rdif, Tdif, Rdir, Tdir, Tnoscat)
      call adding_sw(ncol, nlay, top_is_1,                     &
                     Rdif, Tdif, Rdir, Tdir, Tnoscat,          &
                     sfc_alb_dif(:,igpt), sfc_alb_dir(:,igpt), &
                     flux_up(:,:,igpt), flux_dn(:,:,igpt), flux_dir(:,:,igpt))
      !
      ! adding_sw breaks out direct and diffuse; flux_dn is total to be consistent with LW.
      !
      flux_dn(:,:,igpt) = flux_dn(:,:,igpt) + flux_dir(:,:,igpt)
    end do

  end subroutine sw_solver_2stream

! ---------------------------------------------------------------
! Two-stream solutions to direct and diffuse reflectance and transmittance for a layer
!    with optical depth tau, single scattering albedo w0, and asymmetery parameter g.
!
! Equations are developed in Meador and Weaver, 1980,
!    doi:10.1175/1520-0469(1980)037<0630:TSATRT>2.0.CO;2
!
  pure subroutine two_stream(ncol, nlay, mu0, tau, w0, g, &
                             Rdif, Tdif, Rdir, Tdir, Tnoscat) bind (C)
    integer,                        intent(in)  :: ncol, nlay
    real(wp), dimension(ncol),      intent(in)  :: mu0
    real(wp), dimension(ncol,nlay), intent(in)  :: tau, w0, g
    real(wp), dimension(ncol,nlay), intent(out) :: Rdif, Tdif, Rdir, Tdir, Tnoscat

    ! -----------------------
    integer  :: i, j

    ! Variables used in Meador and Weaver
    real(wp) :: gamma1(ncol), gamma2(ncol), gamma3(ncol), gamma4(ncol)
    real(wp) :: alpha1(ncol), alpha2(ncol), k(ncol)

    ! Ancillary variables
    real(wp) :: RT_term(ncol)
    real(wp) :: exp_minusktau(ncol), exp_minus2ktau(ncol)
    real(WP) :: k_mu, k_gamma3, k_gamma4
    real(wp) :: mu0_inv(ncol)
    ! ---------------------------------
    mu0_inv(1:ncol) = 1._wp/mu0(1:ncol)
    do j = 1, nlay
      do i = 1, ncol
        ! Zdunkowski Practical Improved Flux Method "PIFM"
        !  (Zdunkowski et al., 1980;  Contributions to Atmospheric Physics 53, 147-66)
        !
        gamma1(i)= (8._wp - w0(i,j) * (5._wp + 3._wp * g(i,j))) * .25_wp
        gamma2(i)=  3._wp *(w0(i,j) * (1._wp -         g(i,j))) * .25_wp
        gamma3(i)= (2._wp - 3._wp * mu0(i) *           g(i,j) ) * .25_wp
        gamma4(i)=  1._wp - gamma3(i)

        alpha1(i) = gamma1(i) * gamma4(i) + gamma2(i) * gamma3(i)           ! Eq. 16
        alpha2(i) = gamma1(i) * gamma3(i) + gamma2(i) * gamma4(i)           ! Eq. 17
      end do

      ! Written to encourage vectorization of exponential, square root
      ! Eq 18;  k = SQRT(gamma1**2 - gamma2**2), limited below to avoid div by 0.
      !   k = 0 for isotropic, conservative scattering; this lower limit on k
      !   gives relative error with respect to conservative solution
      !   of < 0.1% in Rdif down to tau = 10^-9
      k(1:ncol) = sqrt(max((gamma1(1:ncol) - gamma2(1:ncol)) * &
                           (gamma1(1:ncol) + gamma2(1:ncol)),  &
                           1.e-12_wp))
      exp_minusktau(1:ncol) = exp(-tau(1:ncol,j)*k(1:ncol))

      !
      ! Diffuse reflection and transmission
      !
      do i = 1, ncol
        exp_minus2ktau(i) = exp_minusktau(i) * exp_minusktau(i)

        ! Refactored to avoid rounding errors when k, gamma1 are of very different magnitudes
        RT_term(i) = 1._wp / (k     (i) * (1._wp + exp_minus2ktau(i))  + &
                              gamma1(i) * (1._wp - exp_minus2ktau(i)) )

        ! Equation 25
        Rdif(i,j) = RT_term(i) * gamma2(i) * (1._wp - exp_minus2ktau(i))

        ! Equation 26
        Tdif(i,j) = RT_term(i) * 2._wp * k(i) * exp_minusktau(i)
      end do

      !
      ! Transmittance of direct, unscattered beam. Also used below
      !
      Tnoscat(1:ncol,j) = exp(-tau(1:ncol,j)*mu0_inv(1:ncol))

      !
      ! Direct reflect and transmission
      !
      do i = 1, ncol
        k_mu     = k(i) * mu0(i)
        k_gamma3 = k(i) * gamma3(i)
        k_gamma4 = k(i) * gamma4(i)

        !
        ! Equation 14, multiplying top and bottom by exp(-k*tau)
        !   and rearranging to avoid div by 0.
        !
        RT_term(i) =  w0(i,j) * RT_term(i)/merge(1._wp - k_mu*k_mu, &
                                                 epsilon(1._wp),    &
                                                 abs(1._wp - k_mu*k_mu) >= epsilon(1._wp))

        Rdir(i,j) = RT_term(i)  *                                        &
            ((1._wp - k_mu) * (alpha2(i) + k_gamma3)                     - &
             (1._wp + k_mu) * (alpha2(i) - k_gamma3) * exp_minus2ktau(i) - &
             2.0_wp * (k_gamma3 - alpha2(i) * k_mu)  * exp_minusktau (i) * Tnoscat(i,j))

        !
        ! Equation 15, multiplying top and bottom by exp(-k*tau),
        !   multiplying through by exp(-tau/mu0) to
        !   prefer underflow to overflow
        !
        Tdir(i,j) = Tnoscat(i,j) - &
                  RT_term(i) * ((1._wp + k_mu) * (alpha1(i) + k_gamma4)                     * Tnoscat(i,j) - &
                                (1._wp - k_mu) * (alpha1(i) - k_gamma4) * exp_minus2ktau(i) * Tnoscat(i,j) - &
                                2.0_wp * (k_gamma4 + alpha1(i) * k_mu)  * exp_minusktau (i))

        Tdir(i,j) = Tdir(i,j) - Tnoscat(i,j)
      end do
    end do
  end subroutine two_stream
! ---------------------------------------------------------------
! Transport of solar radiation through a vertically layered atmosphere.
!   Equations are after Shonk and Hogan 2008, doi:10.1175/2007JCLI1940.1 (SH08)
!
 ! Indexing into arrays for upward and downward propagation depends on the vertical
 !   orientation of the arrays (whether the domain top is at the first or last index)
 ! We write the loops out explicitly so compilers will have no trouble optimizing them.

  pure subroutine adding_sw(ncol, nlay, top_is_1,            &
                            rdif, tdif, Rdir, Tdir, Tnoscat, &
                            sfc_alb_dif, sfc_alb_dir,        &
                            flux_up, flux_dn_dif, flux_dn_dir) bind (C)
    integer,                          intent(in   ) :: ncol, nlay
    logical(wl),                      intent(in   ) :: top_is_1
    real(wp), dimension(ncol,nlay  ), intent(in   ) :: rdif, tdif, Rdir, Tdir, Tnoscat
    real(wp), dimension(ncol       ), intent(in   ) :: sfc_alb_dif, sfc_alb_dir
    real(wp), dimension(ncol,nlay+1), intent(  out) :: flux_up
    ! intent(inout) because top layer includes incident flux
    real(wp), dimension(ncol,nlay+1), intent(inout) :: flux_dn_dif, flux_dn_dir
    ! ------------------
    integer :: ilev
    real(wp), dimension(ncol,nlay+1) :: albedo, &  ! reflectivity to diffuse radiation below this level
                                                   ! alpha in SH08
                                        dir_src    ! source of diffuse radiation from direct radiation
                                                   ! G in SH08
    real(wp), dimension(ncol,nlay  ) :: denom      ! beta in SH08
    ! ------------------
    !
    ! Direct beam flux
    !
    if(top_is_1) then
      ! For the flux at this level, what was the previous level, and which layer has the
      !   radiation just passed through?
      ! layer index = level index - 1
      ! previous level is up (-1)
      do ilev = 2, nlay+1
        flux_dn_dir(:,ilev) = flux_dn_dir(:,ilev-1) * Tnoscat(:,ilev-1)
      end do
    else
      ! layer index = level index
      ! previous level is up (+1)
      do ilev = nlay, 1, -1
        flux_dn_dir(:,ilev) = flux_dn_dir(:,ilev+1) * Tnoscat(:,ilev)
      end do
    end if

    if(top_is_1) then
      ! layer index = level index - 1
      ! previous level is up (-1)
      ilev = nlay + 1
      ! Albedo of lowest level is the surface albedo...
      albedo(:,ilev)  = sfc_alb_dif(:)
      ! ... and source of diffuse radiation is direct beam
      dir_src(:,ilev) = sfc_alb_dir(:) * flux_dn_dir(:,ilev) ! Should this be weighted by cos(mu)?

      !
      ! From bottom to top of atmosphere --
      !   compute albedo and contribution of direct beam to diffuse radiation
      !
      do ilev = nlay, 1, -1
        denom(:, ilev) = 1._wp/(1._wp - rdif(:,ilev)*albedo(:,ilev+1))                 ! Eq 10
        albedo(:,ilev) = rdif(:,ilev) + &
                         tdif(:,ilev)*tdif(:,ilev) * albedo(:,ilev+1) * denom(:,ilev) ! Equation 9
        !
        ! Equation 11 -- source is "scattering of the direct solar beam into the diffuse components"
        !   i.e. reflected flux from this layer (Rdir*flux_dn_dir) and
        !   transmitted through the layer (Tdir*flux_dn_dir) and reflected from layers below (albedo)
        !
        dir_src(:,ilev) = Rdir(:,ilev) * flux_dn_dir(:,ilev) + &
                          tdif(:,ilev) * denom(:,ilev) *       &
                            (dir_src(:,ilev+1) + albedo(:,ilev+1)*Tdir(:,ilev)*flux_dn_dir(:,ilev))
      end do

      ! Eq 12, at the top of the domain upwelling diffuse is due to ...
      ilev = 1
      flux_up(:,ilev) = flux_dn_dif(:,ilev) * albedo(:,ilev) + & ! ... reflection of incident diffuse and
                        dir_src(:,ilev)                          ! scattering by the direct beam below

      !
      ! From the top of the atmosphere downward -- compute fluxes
      !
      do ilev = 2, nlay+1
        ! ... The equation doesn't have the denominator in it but it seems like it should
        flux_dn_dif(:,ilev) = (tdif(:,ilev-1)*flux_dn_dif(:,ilev-1) + &  ! Equation 13
                               rdif(:,ilev-1)*dir_src(:,ilev) +       &
                               Tdir(:,ilev-1)*flux_dn_dir(:,ilev-1)) * denom(:,ilev-1)
        flux_up(:,ilev) = flux_dn_dif(:,ilev) * albedo(:,ilev) + & ! Equation 12
                          dir_src(:,ilev)
      end do
    else
      ! layer index = level index
      ! previous level is up (+1)
      ilev = 1
      ! Albedo of lowest level is the surface albedo...
      albedo(:,ilev)  = sfc_alb_dif(:)
      ! ... and source of diffuse radiation is direct beam
      dir_src(:,ilev) = sfc_alb_dir(:) * flux_dn_dir(:,ilev) ! Should this be weighted by cos(mu)?

      !
      ! From bottom to top of atmosphere --
      !   compute albedo and contribution of direct beam to diffuse radiation
      !
      do ilev = 1, nlay
        denom(:, ilev  ) = 1._wp/(1._wp - rdif(:,ilev)*albedo(:,ilev))                ! Eq 10
        albedo(:,ilev+1) = rdif(:,ilev) + &
                           tdif(:,ilev)*tdif(:,ilev) * albedo(:,ilev) * denom(:,ilev) ! Equation 9
        !
        ! Equation 11 -- source is "scattering of the direct solar beam into the diffuse components"
        !   i.e. reflected flux from this layer (Rdir*flux_dn_dir) and
        !   transmitted through the layer (Tdir*flux_dn_dir) and reflected from layers below (albedo)
        !
        dir_src(:,ilev+1) = Rdir(:,ilev) * flux_dn_dir(:,ilev+1) + &
                            tdif(:,ilev) * denom(:,ilev) *       &
                              (dir_src(:,ilev) + albedo(:,ilev)*Tdir(:,ilev)*flux_dn_dir(:,ilev+1))
      end do

      ! Eq 12, at the top of the domain upwelling diffuse is due to ...
      ilev = nlay+1
      flux_up(:,ilev) = flux_dn_dif(:,ilev) * albedo(:,ilev) + & ! ... reflection of incident diffuse and
                        dir_src(:,ilev)                          ! scattering by the direct beam below

      !
      ! From the top of the atmosphere downward -- compute fluxes
      !
      do ilev = nlay, 1, -1
        ! ... The equation doesn't have the denominator in it but it seems like it should
        flux_dn_dif(:,ilev) = (tdif(:,ilev)*flux_dn_dif(:,ilev+1) + &  ! Equation 13
                               rdif(:,ilev)*dir_src(:,ilev) +       &
                               Tdir(:,ilev)*flux_dn_dir(:,ilev+1)) * denom(:,ilev)
        flux_up(:,ilev) = flux_dn_dif(:,ilev) * albedo(:,ilev) + & ! Equation 12
                          dir_src(:,ilev)
      end do
    end if

  end subroutine adding_sw
end module mo_rte_solver_kernels
