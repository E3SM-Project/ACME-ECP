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

! Output is via type ty_fluxes defined in this module
!   (this means the argument list to rte_lw() need not change).
!   Fluxes by g-point are available within the code. Most applications will expect some
!   reduced version of this detail e.g. broadband flux profiles.
! This implementation provides broadband fluxes but users are free to add others in
!   rte_lw_opt()
!
!
module mo_rte_lw
  use mo_rte_kind,   only: wp
  use mo_spectral_disc, only: ty_spectral_disc
  use mo_optical_props, only: ty_optical_props, &
                              ty_optical_props_arry, ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
  use mo_fluxes,        only: ty_fluxes
  use mo_lw_solver,     only: lw_solver_init, lw_solver
  implicit none
  private

  public :: rte_lw_init, rte_lw

  !
  ! Configuration information
  !
  integer :: nstreams = 0 ! Number of streams: 0 means no scattering
                          ! Require even, positive value when setting

contains
  ! --------------------------------------------------
  ! Initialization function
  ! --------------------------------------------------
  function rte_lw_init(nlwstreams, nangles) result(error_msg)
    integer,           optional, intent( in) :: nlwstreams ! Scattering/no scattering
    integer,           optional, intent( in) :: nangles    ! number of quadrature angles for
                                                           ! no-scattering calculation
    character(len=128)                       :: error_msg

    error_msg = ""
    if(present(nlwstreams)) then
      if(nlwstreams >= 0) then  ! Check for an even number?
        nstreams = nlwstreams
      else
        error_msg = "rte_lw_init: nlwstreams provided is less than 0"
        return
      end if
    end if
    if(present(nangles)) then
      error_msg = lw_solver_init(n_angles = nangles)
      if(error_msg /= "") return
    end if
  end function rte_lw_init

  ! --------------------------------------------------
  !
  ! Interface using only optical properties and source functions as inputs; fluxes as outputs.
  !
  ! --------------------------------------------------
  function rte_lw(optical_props, top_at_1, k_dist, &
                     lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_source, &
                     fluxes, &
                     inc_flux) result(error_msg)
    class(ty_optical_props_arry), intent(in   ) :: optical_props     ! Array of ty_optical_props. This type is abstract
                                                                     ! and needs to be made concrete, either as an array
                                                                     ! (class ty_optical_props_arry) or in some user-defined way
    logical,                      intent(in   ) :: top_at_1          ! Is the top of the domain at index 1?
                                                                     ! (if not, ordering is bottom-to-top)
    real(wp), dimension(:,:,:),   intent(in   ) :: lay_source, &     ! Planck source at layer average temperature
                                                                     ! [W/m2] (ncol, nlay, ngpt)
                                                   lev_source_inc, & ! Planck source at layer edge for radiation,
                                                   lev_source_dec    ! [W/m2] (ncol, nlay+1, ngpt)
                                                                     ! in increasing/decreasing ilay direction
                                                                     ! Includes spectral weighting that accounts for state-dependent
                                                                     ! frequency to g-space mapping
    real(wp), dimension(:,:),     intent(in   ) :: sfc_emis    ! emissivity at surface [] (nband, ncol)
    real(wp), dimension(:,:),     intent(in   ) :: sfc_source     ! source function at surface [W/m2] (ncol, ngpt)
    class(ty_spectral_disc),      intent(in   ) :: k_dist      ! derived type with spectral information
    class(ty_fluxes),             intent(inout) :: fluxes      ! Array of ty_fluxes. Default computes broadband fluxes at all levels
                                                               !   if output arrays are defined. Can be extended per user desires.
    real(wp), dimension(:,:),   &
              target, optional, intent(in   ) :: inc_flux    ! incident flux at domain top [W/m2] (ncol, ngpts)
    character(len=128)                        :: error_msg   ! If empty, calculation was successful
    ! --------------------------------
    !
    ! Local variables
    !
    integer :: ncol, nlay, ngpt, nband, ncomp
    integer :: top_lev
    integer :: comp, icol
    real(wp), dimension(optical_props%get_ncol(),   &
                        optical_props%get_nlay()+1, &
                        optical_props%get_ngpt()) :: gpt_flux_up, gpt_flux_dn
    real(wp), dimension(optical_props%get_ncol(), &
                        optical_props%get_ngpt()) :: sfc_emis_gpt
    ! ------------------------------------------------------------------------------------
    !
    ! Error checking
    !   if inc_flux is present it has the right dimensions, is positive definite
    !
    ! --------------------------------
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
      error_msg = "rte_lw: no space allocated for fluxes"
      return
    end if

    !
    ! k-distribution and optical properties are consistent
    !
    if(k_dist%get_ngpt() /= ngpt) then
      error_msg = "rte_sw_alt: number of gpoints inconsistent between optical_props, k_dist"
      return
    end if

    !
    ! Source functions
    !
    if(any([size(lay_source,1), size(lay_source,2), size(lay_source,3)]            /= [ncol, nlay, ngpt])) &
      error_msg = "rte_lw: lay_source inconsistently sized"
    if(any(lay_source < 0._wp)) &
      error_msg = "rte_lw: lay_source has values < 0"
    if(len_trim(error_msg) > 0) return
    if(any([size(lev_source_inc,1), size(lev_source_inc,2), size(lev_source_inc,3)] /= [ncol, nlay, ngpt])) &
      error_msg = "rte_lw: lev_source_inc inconsistently sized"
    if(any(lev_source_inc < 0._wp)) &
      error_msg = "rte_lw: lev_source_inc has values < 0"
    if(len_trim(error_msg) > 0) return
    if(any([size(lev_source_dec,1), size(lev_source_dec,2), size(lev_source_dec,3)] /= [ncol, nlay, ngpt])) &
      error_msg = "rte_lw: lev_source_dec inconsistently sized"
    if(any(lev_source_dec < 0._wp)) &
      error_msg = "rte_lw: lev_source_dec has values < 0"
    if(len_trim(error_msg) > 0) return

    !
    ! Surface source, emissivity
    !
    if(any([size(sfc_source, 1), size(sfc_source, 2)] /= [ncol, ngpt])) &
      error_msg = "rte_lw: sfc_source inconsistently sized"
    if(any(sfc_source < 0._wp )) &
      error_msg = "rte_lw: sfc_source has values < 0"

    if(any([size(sfc_emis,1), size(sfc_emis,2)] /= [nband, ncol])) &
      error_msg = "rte_lw: sfc_emis inconsistently sized"
    if(any(sfc_emis < 0._wp .or. sfc_emis > 1._wp)) &
      error_msg = "rte_lw: sfc_emis has values < 0 or > 1"
    if(len_trim(error_msg) > 0) return

    !
    ! Incident flux, if present
    !
    if(present(inc_flux)) then
      if(any([size(inc_flux,1), size(inc_flux,2)] /= [ncol, ngpt])) &
        error_msg = "rte_lw: inc_flux inconsistently sized"
      if(any(inc_flux < 0._wp)) &
        error_msg = "rte_lw: inc_flux has values < 0"
    end if

    !
    ! Ensure values of tau, ssa, and g are reasonable
    !
    error_msg =  optical_props%validate()
    if(len_trim(error_msg) > 0) then
      error_msg = "rte_lw: " // trim(error_msg)
      return
    end if

    ! ------------------------------------------------------------------------------------
    ! Apply boundary conditions
    !   Upper boundary condition
    top_lev = MERGE(1, nlay+1, top_at_1)
    if(present(inc_flux)) then
      gpt_flux_dn(:,top_lev,:) = inc_flux(:,:)
    else
      gpt_flux_dn(:,top_lev,:) = 0._wp
    end if

    !    Lower boundary condition -- expand surface emissivity by band to gpoints
    do icol = 1, ncol
      sfc_emis_gpt(icol, 1:ngpt) = k_dist%expand(sfc_emis(:,icol))
    end do

    !
    ! Compute the radiative transfer...
    !
    error_msg = lw_solver(ncol, nlay, ngpt, top_at_1,                                &
                          optical_props, lay_source, lev_source_inc, lev_source_dec, &
                          sfc_emis_gpt, sfc_source, gpt_flux_up, gpt_flux_dn)
    if (error_msg /= '') return
    !
    ! ...and reduce spectral fluxes to desired output quantities
    !
    error_msg = fluxes%reduce(gpt_flux_up, gpt_flux_dn, k_dist, top_at_1)
  end function rte_lw
end module mo_rte_lw
