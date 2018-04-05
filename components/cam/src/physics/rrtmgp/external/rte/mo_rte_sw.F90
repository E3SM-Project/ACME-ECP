! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2015-2017,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!

!
! Output is via type ty_fluxes defined in this module
!   (this means the argument list to rte_sw() need not change).
!   Fluxes by g-point are available within the code. Most applications will expect some
!   reduced version of this detail e.g. broadband flux profiles.
! This implementation provides broadband fluxes but users are free to add others in
!   rte_sw_opt()
!
!
module mo_rte_sw
  use mo_rte_kind,   only: wp
  use mo_spectral_disc, only: ty_spectral_disc
  use mo_optical_props, only: ty_optical_props, &
                              ty_optical_props_arry, ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
  use mo_fluxes,        only: ty_fluxes
  use mo_sw_solver,     only: sw_solver
  implicit none
  private

  public :: rte_sw_init, rte_sw

  !
  ! Configuration information
  !
  integer :: nstreams = 2 ! Number of streams: 0 means no scattering
                          ! Require even, positive value when setting

contains
  ! --------------------------------------------------
  ! Initialization function
  ! --------------------------------------------------
  function rte_sw_init(nswstreams) result(error_msg)
    integer,           optional, intent( in) :: nswstreams
    character(len=128)                       :: error_msg

    error_msg = ""
    if(present(nswstreams)) then
      if(nswstreams >= 0) then  ! Check for an even number?
        nstreams = nswstreams
      else
        error_msg = "rte_sw_init: nswstreams provided is less than 0"
      end if
    end if
  end function rte_sw_init

  ! --------------------------------------------------
  ! Public interfaces to rte_sw()
  ! --------------------------------------------------
  !
  ! Interface using source functions and vectors of optical properties as inputs; vectors of fluxes as outputs.
  !
  ! --------------------------------------------------
  function rte_sw(optical_props, top_at_1, k_dist, &
                     mu0, inc_flux,                   &
                     sfc_alb_dir, sfc_alb_dif,        &
                     fluxes) result(error_msg)
    class(ty_optical_props_arry), intent(in   ) :: optical_props   ! Array of ty_optical_props. This type is abstract
                                                                   ! and needs to be made concrete, either as an array
                                                                   ! (class ty_optical_props_arry) or in some user-defined way
    logical,                      intent(in   ) :: top_at_1        ! Is the top of the domain at index 1?
                                                                 ! (if not, ordering is bottom-to-top)
    real(wp), dimension(:),       intent(in   ) :: mu0             ! cosine of solar zenith angle (ncol)
    real(wp), dimension(:,:),     intent(in   ) :: inc_flux,    &  ! incident flux at top of domain [W/m2] (ncol, ngpt)
                                                   sfc_alb_dir, &  ! surface albedo for direct and
                                                   sfc_alb_dif     ! diffuse radiation (nband, ncol)
    class(ty_spectral_disc),      intent(in   ) :: k_dist          ! derived type with spectral information
    class(ty_fluxes),             intent(inout) :: fluxes          ! Array of ty_fluxes. Default computes broadband fluxes at all levels
                                                                   !   if output arrays are defined. Can be extended per user desires.
    character(len=128)                          :: error_msg       ! If empty, calculation was successful
    ! --------------------------------
    !
    ! Local variables
    !
    integer :: ncol, nlay, ngpt, nband
    integer :: top_lev, icol
    integer :: alloc_stat, k, i

    real(wp), dimension(optical_props%get_ncol(),   &
                        optical_props%get_nlay()+1, &
                        optical_props%get_ngpt()) :: gpt_flux_up, gpt_flux_dn, gpt_flux_dir
    real(wp), dimension(optical_props%get_ncol(), &
                        optical_props%get_ngpt()) :: sfc_alb_dir_gpt, sfc_alb_dif_gpt
    ! ------------------------------------------------------------------------------------
    ncol  = optical_props%get_ncol()
    nlay  = optical_props%get_nlay()
    ngpt  = optical_props%get_ngpt()
    nband = k_dist%get_nband()
    error_msg = ""

    ! ------------------------------------------------------------------------------------
    !
    ! Error checking -- consistency of sizes and validity of values
    !
    ! --------------------------------
    if(.not. fluxes%are_desired()) then
      error_msg = "rte_sw: no space allocated for fluxes"
      return
    end if

    !
    ! k-distribution and optical properties are consistent
    !
    if(k_dist%get_ngpt() /= ngpt) then
      error_msg = "rte_sw: number of gpoints inconsistent between optical_props, k_dist"
      return
    end if
    !
    ! Sizes and values of input arrays
    !
    if(     size(mu0, 1)                                /=  ncol        ) &
      error_msg = "rte_sw: mu0 inconsistently sized"
    if(any(mu0 <= 0._wp)) &
      error_msg = "rte_sw: one or more mu0 <= 0"
    if(any([size(inc_flux, 1),    size(inc_flux, 2)]    /= [ncol, ngpt])) &
      error_msg = "rte_sw: inc_flux inconsistently sized"
    if(any(inc_flux <  0._wp)) &
      error_msg = "rte_sw: one or more inc_flux < 0"

    if(any([size(sfc_alb_dir, 1), size(sfc_alb_dir, 2)] /= [nband, ncol])) &
      error_msg = "rte_sw: sfc_alb_dir inconsistently sized"
    if(any(sfc_alb_dir < 0._wp .or. sfc_alb_dir > 1._wp)) &
      error_msg = "rte_sw: sfc_alb_dir out of bounds [0,1]"
    if(any([size(sfc_alb_dif, 1), size(sfc_alb_dif, 2)] /= [nband, ncol])) &
      error_msg = "rte_sw: sfc_alb_dif inconsistently sized"
    if(any(sfc_alb_dif < 0._wp .or. sfc_alb_dif > 1._wp)) &
      error_msg = "rte_sw: sfc_alb_dif out of bounds [0,1]"

    if(len_trim(error_msg) > 0) return

    !
    ! Ensure values of tau, ssa, and g are reasonable
    !
    error_msg =  optical_props%validate()
    if(len_trim(error_msg) > 0) then
      if(len_trim(optical_props%get_name()) > 0) &
        error_msg = 'rte_sw: ' //  trim(error_msg)
      return
    end if

    ! ------------------------------------------------------------------------------------
    ! Lower boundary condition -- expand surface albedos by band to gpoints
    !   and switch dimension ordering
    do icol = 1, ncol
      sfc_alb_dir_gpt(icol, 1:ngpt) = k_dist%expand(sfc_alb_dir(:,icol))
    end do
    do icol = 1, ncol
      sfc_alb_dif_gpt(icol, 1:ngpt) = k_dist%expand(sfc_alb_dif(:,icol))
    end do

    ! ------------------------------------------------------------------------------------
    !
    ! Apply boundary conditions
    top_lev = MERGE(1, nlay+1, top_at_1)
    ! Should allow for diffuse upper boundary condition for generality
    gpt_flux_dn(:,top_lev,:) = 0._wp
    gpt_flux_dir(:,top_lev,:) = inc_flux(:,:) * spread(mu0(:), dim = 2, ncopies = ngpt)

    !
    ! Compute the radiative transfer...
    !
    error_msg = sw_solver(ncol, nlay, ngpt, top_at_1,                     &
                          optical_props, mu0, sfc_alb_dir_gpt, sfc_alb_dif_gpt, &
                          gpt_flux_up, gpt_flux_dn, gpt_flux_dir)
    if (error_msg /= '') return
    !
    ! ...and reduce spectral fluxes to desired output quantities
    !
    error_msg = fluxes%reduce(gpt_flux_up, gpt_flux_dn, k_dist, top_at_1, gpt_flux_dir)
  end function rte_sw
end module mo_rte_sw
