! Module: mo_gas_optics

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
! Description:  Specifies all properties related to the k-distributions.  This includes each band's spectral properties and
! absorbing gases.

module mo_gas_optics
  use mo_rte_kind,        only: wp, wl
  use mo_rrtmgp_constants,   only: avogad, m_dry, m_h2o, grav
  use mo_spectral_disc,      only: ty_spectral_disc
  use mo_gas_optics_kernels, only: interpolation, gas_optical_depths_major, &
                                   gas_optical_depths_minor, &
                                   gas_optical_depths_rayleigh, source
  use mo_util_string,        only : lower_case, string_in_array, string_loc_in_array
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_optical_props,      only: ty_optical_props_arry, ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
  use mo_util_reorder
  implicit none
  private

  real(wp), parameter :: pi = acos(-1._wp)

  ! -----------------------------------------------------------------------------------
  type, extends(ty_spectral_disc), public :: ty_gas_optics_specification
    private
    character(32), &
              dimension(:),   allocatable :: gas_names  ! gas names

    integer,  dimension(:,:), allocatable :: flavor        ! major species pair; (2,nflav)
    integer,  dimension(:,:), allocatable :: gpoint_flavor ! flavor = gpoint_flavor(lower or upper atmosphere, g-point)

    ! -----------------------------------------------------------------------------------
    ! Temperature and pressure interpolation grids
    real(wp), dimension(:),  allocatable :: press_ref,  press_ref_log, temp_ref
    ! Volume mixing ratios for reference atmosphere; vmr_ref(lower or upper atmosphere, gas, temp)
    real(wp), dimension(:,:,:), allocatable :: vmr_ref

    real(wp) :: press_ref_min, press_ref_max, &  ! min and max pressure of interpolation grids
                temp_ref_min,  temp_ref_max      ! min and max temperature
    real(wp) :: press_ref_log_delta, & ! difference in ln pressure between consecutive reference levels
                temp_ref_delta,      & ! temperature difference between consecutive reference levels
                press_ref_trop_log     ! log of reference pressure separating the lower and upper atmosphere
    ! -----------------------------------------------------------------------------------
    ! Absorption coefficients
      ! ----- major gas (also referred to as "key species")absorption coefficients ; kmajor(g-point,eta,pressure,temperature)
    real(wp), dimension(:,:,:,:), allocatable :: kmajor
    ! ----- minor species
      ! stored absorption coefficients due to minor absorbing gases in lower/upper part of atmosphere;
      ! kminor_lower(contributor,eta,temperature)
    real(wp), dimension(:,:,:), allocatable :: kminor_lower, kminor_upper
    ! Name of the absorbing gas & unique identifying name (n_minor)
    character(len=256), dimension(:), allocatable :: gas_minor, &
      identifier_minor
    ! Description of each minor gas contribution, separately for upper and lower atmospheres
  	!
  	! Arrays are dimensioned with n_minor_lower, n_minor_upper
  	!
  	! Name of the absorbing gas listed once for each band in which it is active (n_minor)
    character(len=256), dimension(:), allocatable :: minor_gases_lower, &
      minor_gases_upper
  	! Starting and ending g-points for each minor gas in each band (2, n_minor)
    integer, dimension(:,:), allocatable :: minor_limits_gpt_lower, &
      minor_limits_gpt_upper
    ! Does the minor gas absorption coefficient scale with density? (n_minor)
    logical, dimension(:), allocatable :: minor_scales_with_density_lower
    logical, dimension(:), allocatable :: minor_scales_with_density_upper
  	! Does the minor gas absorption coefficient scale with the amount of another gas? (n_minor)
    character(len=256), dimension(:), allocatable :: scaling_gas_lower, &
      scaling_gas_upper
  	! If the minor gas absorption coefficient depends on another gas, does it depend on the concentration itself or
  	!   the concentration of all gases besides the scaling gas? (n_minor)
    logical, dimension(:), allocatable :: scale_by_complement_lower, &
      scale_by_complement_upper
    integer, dimension(:), allocatable :: kminor_start_lower, kminor_start_upper
    ! -----------------------------------------------------------------------------------
    ! ----- Rayleigh scattering
      ! stored scattering coefficients due to molecules in atmosphere;
      ! krayl(g-point,eta,temperature,upper/lower atmosphere)
    real(wp), dimension(:,:,:,:), allocatable :: krayl
    ! -----------------------------------------------------------------------------------
    ! Planck function spectral mapping
    !   Allocated only when gas optics object is internal-source
    !
    real(wp), dimension(:,:,:,:), allocatable :: planck_frac   ! stored fraction of Planck irradiance in band for given g-point
                                                               ! planck_frac(eta,temperature,pressure,g-point)
    real(wp), dimension(:,:),     allocatable :: totplnk       ! integrated Planck irradiance by band; (reference temperatures,band)
    real(wp)                                  :: totplnk_delta ! temperature steps in totplnk
    ! -----------------------------------------------------------------------------------
    ! Solar source function spectral mapping
    !   Allocated only when gas optics object is external-source
    !
    real(wp), dimension(:), allocatable :: solar_src ! incoming solar irradiance(g-point)
    ! -----------------------------------------------------------------------------------
    ! Ancillary
    ! -----------------------------------------------------------------------------------
    ! Index into %gas_names -- is this a key species in any band?
    logical, dimension(:), allocatable :: is_key
    ! -----------------------------------------------------------------------------------

  contains
    ! Type-bound procedures
    ! Public procedures
    ! public interface
    generic,   public :: init       => init_int,       init_ext
    generic,   public :: gas_optics => gas_optics_int, gas_optics_ext
    procedure, public :: is_internal_source_present
    procedure, public :: is_external_source_present
    procedure, public :: get_ngas
    procedure, public :: get_gases
    procedure, public :: get_press_ref_min
    procedure, public :: get_press_ref_max
    procedure, public :: get_temp_ref_min
    procedure, public :: get_temp_ref_max
    ! Internal procedures
    procedure, private :: init_int
    procedure, private :: init_ext
    procedure, private :: gas_optics_int
    procedure, private :: gas_optics_ext
    procedure, private :: check_key_species_present
    procedure, private :: get_minor_list
    procedure, private :: get_nflav
    procedure, private :: get_nlay_ref
    procedure, private :: get_neta
    procedure, private :: compute_gas_tau_core
  end type
  ! -----------------------------------------------------------------------------------
  public :: get_col_dry ! Utility function, not type-bound

  interface check_range
    module procedure check_range_1D, check_range_2D, check_range_3D
  end interface check_range

  interface check_extent
    module procedure check_extent_1D, check_extent_2D, check_extent_3D
  end interface check_extent
contains
  ! --------------------------------------------------------------------------------------
  !
  ! Public procedures
  !
  ! --------------------------------------------------------------------------------------
  !
  ! Two functions to define array sizes needed by gas_optics()
  !
  pure function get_ngas(this)
    ! return the number of gases registered in the spectral configuration
    class(ty_gas_optics_specification), intent(in) :: this
    integer                                        :: get_ngas

    get_ngas = size(this%gas_names)
  end function get_ngas
  !--------------------------------------------------------------------------------------------------------------------
  pure function get_nflav(this)
    ! return the number of distinct major gas pairs in the spectral bands (referred to as
    ! "flavors" - all bands have a flavor even if there is one or no major gas)
    class(ty_gas_optics_specification), intent(in) :: this
    integer                                        :: get_nflav

    get_nflav = size(this%flavor,dim=2)
  end function get_nflav
  !--------------------------------------------------------------------------------------------------------------------

  ! Compute gas optical depth and, optionally, Planck source functions,
  !  given temperature, pressure, and composition
  function gas_optics_int(this,                                   &
                      play, plev, tlay, tsfc, gas_desc,           & ! mandatory inputs
                      optical_props,                              & ! mandatory outputs
                      lay_src, lev_src_inc, lev_src_dec, sfc_src, & ! internal-source specific outputs
                      col_dry, tlev)                              & ! optional inputs
                      result(error_msg)
    ! inputs
    class(ty_gas_optics_specification), intent(in) :: this
    real(wp), dimension(:,:), intent(in   ) :: play, &   ! layer pressures [Pa, mb]; (ncol,nlay)
                                               plev, &   ! level pressures [Pa, mb]; (ncol,nlay+1)
                                               tlay      ! layer temperatures [K]; (ncol,nlay)
    real(wp), dimension(:),   intent(in   ) :: tsfc      ! surface skin temperatures [K]; (ncol)
    type(ty_gas_concs),       intent(in   ) :: gas_desc  ! Gas volume mixing ratios
    ! output
    class(ty_optical_props_arry),  &
                              intent(inout) :: optical_props
    character(len=128)                      :: error_msg
    ! internal source functions (LW only)
    ! These include spectral weighting that accounts for state-dependent frequency to k-distribution mapping
    ! [W/m2]
    real(wp), dimension(:,:,:), intent(  out) :: lay_src, &  ! source for average layer temperature; (ncol,nlay,ngpt)
                                                 lev_src_inc, lev_src_dec
                                                             ! level source radiances in increasing/decreasing
                                                             ! ilay direction (ncol,nlay+1,ngpt)
    real(wp), dimension(:,:),   intent(  out) :: sfc_src     ! surface Planck source; (ncol,ngpts)
    ! Optional inputs
    real(wp), dimension(:,:),   intent(in   ), &
                           optional, target :: col_dry, &  ! Column dry amount; dim(ncol,nlay)
                                               tlev        ! level temperatures [K]l (ncol,nlay+1)
    ! ----------------------------------------------------------
    ! Local variables
    ! Interpolation coefficients to save for use in source function
    integer,  dimension(size(play,dim=1), size(play,dim=2)) :: jtemp, jpress
    logical,  dimension(size(play,dim=1), size(play,dim=2)) :: tropo
    real(wp), dimension(2,2,2,this%get_nflav(),size(play,dim=1), size(play,dim=2)) :: fmajor
    integer,  dimension(2,    this%get_nflav(),size(play,dim=1), size(play,dim=2)) :: jeta
    ! ----------------------------------------------------------
    ! dimensions - determined from problem size
    integer :: ncol, nlay ! number of columns, layers
    ! dimensions - provided by k-distribution
    integer :: ngpt, nband, ngas, nflav ! Number of g-points, bands, gas, gas "flavors" (major species combinations)

    ! index
    integer :: icol, ilay, igpt, igas

    ! Variables for temperature at layer edges
    ! [K] (ncol, nlay+1)
    real(wp), dimension(size(play,dim=1),size(play,dim=2)+1), target  :: tlev_arr
    real(wp), dimension(:,:),                                 pointer :: tlev_wk => NULL()

    ! ----------------------------------------------------------
    ! Code starts
    !
    error_msg = ""
    error_msg = compute_gas_taus(this,                       &
                                 play, plev, tlay, gas_desc, &
                                 optical_props,              &
                                 jtemp, jpress, jeta, tropo, fmajor, &
                                 col_dry)
    if(error_msg  /= '') return

    ! init from array dimensions
    ncol = size(play,dim=1)
    nlay = size(play,dim=2)
    ngpt = this%get_ngpt()
    nband = this%get_nband()
    ngas = this%get_ngas()
    nflav = this%get_nflav()

    !
    ! Planck source function
    !   Check input data sizes and values
    !
    error_msg = check_extent(tsfc, ncol, 'tsfc')
    if(error_msg  /= '') return
    error_msg = check_range(tsfc, this%temp_ref_min,  this%temp_ref_max,  'tsfc')
    if(error_msg  /= '') return
    if(present(tlev)) then
      error_msg = check_extent(tlev, ncol, nlay+1, 'tlev')
      if(error_msg  /= '') return
      error_msg = check_range(tlev, this%temp_ref_min, this%temp_ref_max, 'tlev')
      if(error_msg  /= '') return
    end if

    !
    !   output extents
    !
    error_msg = check_extent(sfc_src,     ncol,       ngpt, 'sfc_src')
    if(error_msg  /= '') return
    error_msg = check_extent(lay_src,     ncol, nlay, ngpt, 'lay_src')
    if(error_msg  /= '') return
    error_msg = check_extent(lev_src_inc, ncol, nlay, ngpt, 'lev_src_inc')
    if(error_msg  /= '') return
    error_msg = check_extent(lev_src_dec, ncol, nlay, ngpt, 'lev_src_dec')
    if(error_msg  /= '') return

    !
    ! Source function needs temperature at interfaces/levels and at layer centers
    !
    if (present(tlev)) then
      !   Users might have provided these
      tlev_wk => tlev
    else
       tlev_wk => tlev_arr
       !
       ! Interpolate temperature to levels if not provided
       !   Interpolation and extrapolation at boundaries is weighted by pressure
       !
       do icol = 1, ncol
         tlev_arr(icol,1) = tlay(icol,1) &
                           + (plev(icol,1)-play(icol,1))*(tlay(icol,2)-tlay(icol,1))  &
              &                                           / (play(icol,2)-play(icol,1))
       end do
       do ilay = 2, nlay
         do icol = 1, ncol
           tlev_arr(icol,ilay) = (play(icol,ilay-1)*tlay(icol,ilay-1)*(plev(icol,ilay  )-play(icol,ilay)) &
                                +  play(icol,ilay  )*tlay(icol,ilay  )*(play(icol,ilay-1)-plev(icol,ilay))) /  &
                                  (plev(icol,ilay)*(play(icol,ilay-1) - play(icol,ilay)))
         end do
       end do
       do icol = 1, ncol
         tlev_arr(icol,nlay+1) = tlay(icol,nlay)                                                             &
                                + (plev(icol,nlay+1)-play(icol,nlay))*(tlay(icol,nlay)-tlay(icol,nlay-1))  &
                                                                      / (play(icol,nlay)-play(icol,nlay-1))
       end do
     end if

!   Get internal source functions at layers and levels, which depend on mapping from spectral space that creates k-distribution.
    call source(ncol, nlay, ngpt, nband, ngas, nflav, &
                tlay, tlev_wk, tsfc, merge(1,nlay,play(1,1) > play(1,nlay)), &
                fmajor, jeta, tropo, jtemp, jpress,                    &
                this%get_gpoint_bands(), this%planck_frac, this%temp_ref_min,    &
                this%totplnk_delta, this%totplnk, this%gpoint_flavor,  &
                sfc_src, lay_src, lev_src_inc, lev_src_dec)


  end function gas_optics_int
  !------------------------------------------------------------------------------------------

  ! Compute gas optical depth
  !  given temperature, pressure, and composition
  function gas_optics_ext(this,        &
    play, plev, tlay, gas_desc,        & ! mandatory inputs
    optical_props, toa_src,            & ! mandatory outputs
    col_dry) result(error_msg)           ! optional input

    class(ty_gas_optics_specification), intent(in) :: this
    real(wp), dimension(:,:), intent(in   ) :: play, &   ! layer pressures [Pa, mb]; (ncol,nlay)
                                               plev, &   ! level pressures [Pa, mb]; (ncol,nlay+1)
                                               tlay      ! layer temperatures [K]; (ncol,nlay)
    type(ty_gas_concs),       intent(in   ) :: gas_desc  ! Gas volume mixing ratios
    ! output
    class(ty_optical_props_arry),  &
                              intent(inout) :: optical_props
    real(wp), dimension(:,:), intent(  out) :: toa_src     ! Incoming solar irradiance(ncol,ngpt)
    character(len=128)                      :: error_msg

    ! Optional inputs
    real(wp), dimension(:,:), intent(in   ), &
                           optional, target :: col_dry ! Column dry amount; dim(ncol,nlay)
    ! ----------------------------------------------------------
    ! Local variables
    ! Interpolation coefficients
    integer,  dimension(size(play,dim=1), size(play,dim=2)) :: jtemp, jpress
    logical,  dimension(size(play,dim=1), size(play,dim=2)) :: tropo
    real(wp), dimension(2,2,2,this%get_nflav(),size(play,dim=1), size(play,dim=2)) :: fmajor
    integer,  dimension(2,    this%get_nflav(),size(play,dim=1), size(play,dim=2)) :: jeta
    ! ----------------------------------------------------------
    integer :: ncol, ngpt

    ! ----------------------------------------------------------
    ! Code starts
    !
    error_msg = ""
    error_msg = compute_gas_taus(this,                       &
                                 play, plev, tlay, gas_desc, &
                                 optical_props,              &
                                 jtemp, jpress, jeta, tropo, fmajor, &
                                 col_dry)
    if(error_msg  /= '') return

    ncol = size(play,dim=1)
    ngpt = this%get_ngpt()
    error_msg = check_extent(toa_src,     ncol,         ngpt, 'toa_src')
    if(error_msg  /= '') return
    toa_src(:,:) = spread(this%solar_src(:), dim=1, ncopies=ncol)

  end function gas_optics_ext
  !------------------------------------------------------------------------------------------
  !
  ! Returns optical properties and interpolation coefficients
  !
  function compute_gas_taus(this,                       &
                            play, plev, tlay, gas_desc, &
                            optical_props,              &
                            jtemp, jpress, jeta, tropo, fmajor, &
                            col_dry) result(error_msg)

    class(ty_gas_optics_specification), &
                                      intent(in   ) :: this
    real(wp), dimension(:,:),         intent(in   ) :: play, &   ! layer pressures [Pa, mb]; (ncol,nlay)
                                                       plev, &   ! level pressures [Pa, mb]; (ncol,nlay+1)
                                                       tlay      ! layer temperatures [K]; (ncol,nlay)
    type(ty_gas_concs),               intent(in   ) :: gas_desc  ! Gas volume mixing ratios

    class(ty_optical_props_arry),     intent(inout) :: optical_props
    ! Interpolation coefficients for use in internal source function
    integer,  dimension(:,:),         intent(  out) :: jtemp, jpress
    integer,  dimension(:,:,:,:),     intent(  out) :: jeta
    logical,  dimension(:,:),         intent(  out) :: tropo
    real(wp), dimension(:,:,:,:,:,:), intent(  out) :: fmajor
    character(len=128)                            :: error_msg

    ! Optional inputs
    real(wp), dimension(:,:), intent(in   ), &
                           optional, target :: col_dry ! Column dry amount; dim(ncol,nlay)
    ! ----------------------------------------------------------
    ! Local variables
    ! gas amounts
    real(wp), dimension(size(play,dim=1), size(play,dim=2))                  :: one_vmr ! a single volume mixing ratio, (ncol, nlay)

    real(wp), dimension(size(optical_props%tau,dim=3), &
                        size(optical_props%tau,dim=2), &
                        size(optical_props%tau,dim=1)) :: tau  ! optical depth; (ngpt, nlay, ncol)
    real(wp), dimension(size(optical_props%tau,dim=3), &
                        size(optical_props%tau,dim=2), &
                        size(optical_props%tau,dim=1)) :: tau_rayleigh ! optical depth; (ngpt, nlay, ncol)

    logical, dimension(this%get_ngas())              :: gas_is_present  ! Is the concentration of each gas known to the
                                                                        !   k-distribution available in the set of concentrations?
    integer, dimension(2, this%get_nflav())          :: flavor          ! indices of two major absorbing gases per band, referring to
                                                                        ! gas-which-are-present

    real(wp),          dimension(:,:,:), allocatable :: vmr             ! volume mixing ratios; (nlay,ncol,ngas)
    character(len=32), dimension(:),     allocatable :: terse_gas_names ! The gases known to the k-distribution with
                                                                        ! concentrations present
    character(len=32), dimension(:), allocatable     :: minor_gas_list
    integer                                          :: imnr
    ! ----------------------------------------------------------
    ! dimensions - determined from problem size
    integer :: ncol, nlay ! number of columns, layers
    ! dimensions - provided by k-distribution
    integer :: ngpt, nband, ngas, nflav ! Number of g-points, bands, gas, gas "flavors" (major species combinations)
    ! index
    integer :: igas, iband, idx_h2o

    ! Number of molecules per cm^2
    real(wp), dimension(size(play,dim=1), size(play,dim=2)), target  :: col_dry_arr
    real(wp), dimension(:,:),                                pointer :: col_dry_wk => NULL()
    ! ----------------------------------------------------------
    ! Code starts
    !
    error_msg = ''
    ! Check for initialization
    if (.not. this%is_initialized()) then
      error_msg = 'ERROR: spectral configuration not loaded'
      return
    end if
    !
    ! Check for presence of key species in ty_gas_concs; return error if any key species are not present
    !
    error_msg = this%check_key_species_present(gas_desc)
    if (error_msg /= '') return
    ! init from array dimensions
    ncol = size(play,dim=1)
    nlay = size(play,dim=2)
    ngpt = this%get_ngpt()
    nband = this%get_nband()
    ngas = this%get_ngas()
    nflav = this%get_nflav()

    !
    ! Check input data sizes and values
    !
    error_msg = check_extent(play, ncol, nlay,   'play')
    if(error_msg  /= '') return
    error_msg = check_extent(plev, ncol, nlay+1, 'plev')
    if(error_msg  /= '') return
    error_msg = check_extent(tlay, ncol, nlay,   'tlay')
    if(error_msg  /= '') return
    error_msg = check_range(play, this%press_ref_min,this%press_ref_max, 'play')
    if(error_msg  /= '') return
    error_msg = check_range(plev, this%press_ref_min, this%press_ref_max, 'plev')
    if(error_msg  /= '') return
    error_msg = check_range(tlay, this%temp_ref_min,  this%temp_ref_max,  'tlay')
    if(error_msg  /= '') return
    if(present(col_dry)) then
      error_msg = check_extent(col_dry, ncol, nlay, 'col_dry')
      if(error_msg  /= '') return
      error_msg = check_range(col_dry, 0._wp, huge(col_dry), 'col_dry')
      if(error_msg  /= '') return
    end if

    !
    ! Code to be replaced when gas optics calculations are more thoroughly kernel-ized
    !
    allocate(vmr(ncol, nlay, ngas))
    do igas = 1, ngas
      ! Get vmr only for gases provided in ty_gas_concs
      if (any (lower_case(this%gas_names(igas)) == gas_desc%gas_name(:))) then
         error_msg = gas_desc%get_vmr(this%gas_names(igas),one_vmr)
         if (error_msg /= '') return
         vmr(:,:,igas) = one_vmr
      else
      ! Temporarily set missing gas amounts to zero; may not be needed when missing
      ! gases skipped in calculation
         vmr(:,:,igas) = 0._wp
      end if
    end do

    !
    ! Construct arrays of mixing ratios using only those gases that are both known to the
    !   k-distribution and have concentrations available
    !   Revise flavors and minor gases activities
    !
    if(.false.) then
      !
      ! Which gases from the k-distribution are present in the set of gas_concentrations?
      !
      do igas = 1, this%get_ngas()
        gas_is_present(igas) = string_in_array(this%gas_names(igas), gas_desc%gas_name)
      end do
      ngas = count(gas_is_present)
      allocate(terse_gas_names(ngas), vmr(ncol,nlay,ngas))
      terse_gas_names(:) = pack(this%gas_names, mask=gas_is_present)
      !
      ! Expand volume mixing ratio of available gases to 3D fields
      !
      do igas = 1, ngas
        error_msg = gas_desc%get_vmr(terse_gas_names(igas), one_vmr)
        if (error_msg /= '') return
        vmr(:,:,igas) = one_vmr
      end do
      !
      ! Revise mappings into concentration arrays
      !
      do iband = 1, this%get_nflav()
        flavor(1, iband) = string_loc_in_array(this%gas_names(this%flavor(1, iband)), terse_gas_names)
        flavor(2, iband) = string_loc_in_array(this%gas_names(this%flavor(2, iband)), terse_gas_names)
      end do
    end if
    !
    !
    !

    ! Compute dry air column amounts (number of molecule per cm^2) if user hasn't provided them
    if (present(col_dry)) then
      col_dry_wk => col_dry
    else
      ! terse_gas_names never gets allocated or assigned because the above block
      ! of code is wrapped in an if (.false.), so we get gas names from
      ! this%gas_names, which seems to be what corresponds with the construction
      ! of the vmr array
      idx_h2o = string_loc_in_array('h2o', this%gas_names)
      if (idx_h2o > 0) then
         col_dry_arr = get_col_dry(vmr(:,:,idx_h2o), plev, tlay) ! dry air column amounts computation
         col_dry_wk => col_dry_arr
      else
         error_msg = 'compute_gas_taus: h2o not found in gas list.'
         return
      end if
    end if

    ! Make list of minor gases that are defined in specification and have available concentrations
    ! Includes key species that are also considered minor at some g-points
    minor_gas_list = this%get_minor_list(gas_desc, ngas, this%gas_names)

    ! Compute gas optical depths.
    error_msg = this%compute_gas_tau_core(play, tlay, vmr, col_dry_wk, minor_gas_list, &
                                     ncol, nlay, ngpt, nband, ngas, nflav, &
                                     tau, tau_rayleigh, &
                                     fmajor, jeta, tropo, jtemp, jpress)
    if (error_msg /= '') return

    ! Combine optical depths and reorder for radiative transfer solver.
    call combine_and_reorder(tau, tau_rayleigh, allocated(this%krayl), optical_props)

  end function compute_gas_taus
  !------------------------------------------------------------------------------------------
  !
  ! This function should be made into a kernel.
  !
  function compute_gas_tau_core(this,            &
    play, tlay, vmr, col_dry, minor_gas_list, & !  inputs
    ncol, nlay, ngpt, nband, ngas, nflav,   &
    tau, tau_rayleigh,                      & ! mandatory outputs
    fmajor_out, jeta_out, tropo_out, jtemp_out, jpress_out) result(error_msg)

    class(ty_gas_optics_specification), intent(in) :: this

    ! dimensions
    integer, intent(in) :: ncol  ! Number of columns
    integer, intent(in) :: nlay  ! Number of layers
    integer, intent(in) :: ngpt  ! Number of gpts
    integer, intent(in) :: nband ! Number of bands
    integer, intent(in) :: ngas  ! Number of gases
    integer, intent(in) :: nflav ! Number of gas flavors

    real(wp), dimension(ncol,nlay  ), intent(in) :: play   ! Layer pressures [Pa, mb]
    real(wp), dimension(ncol,nlay  ), intent(in) :: tlay   ! Layer temperatures [K]
    real(wp), dimension(ncol,nlay  ,ngas), &
                                      intent(in) :: vmr ! volume mixing ratios
    real(wp), dimension(ncol,nlay  ), intent(in) :: col_dry ! Column amount of dry air
    ! List of minor gases to be processed
    character(len=32), dimension(:), intent(in)  :: minor_gas_list

    ! output
    real(wp), dimension(ngpt,nlay,ncol), intent(out) :: tau          ! gas absorption optical depth
    real(wp), dimension(ngpt,nlay,ncol), intent(out) :: tau_rayleigh ! Rayleigh scattering optical depth

    real(wp), dimension(2,2,2,nflav,ncol,nlay), optional, intent(out) :: fmajor_out
    integer,  dimension(2,    nflav,ncol,nlay), optional, intent(out) :: jeta_out
    logical,  dimension(            ncol,nlay), optional, intent(out) :: tropo_out
    integer,  dimension(            ncol,nlay), optional, intent(out) :: jtemp_out
    integer,  dimension(            ncol,nlay), optional, intent(out) :: jpress_out

    ! result
    character(len=128) :: error_msg
    ! ----------------------------------------------------------
    ! Local variables
    ! index
    integer :: igas
    integer,  dimension(ngpt) :: gpt_flv_lower, gpt_flv_upper
    integer :: idx_h2o ! index of some gases
    ! Planck fractions
    ! gas amounts
    real(wp), dimension(ncol,nlay,0:ngas) :: col_gas ! column amounts for each gas
    integer, dimension(ngas)  :: idx_gas_list      ! Index of minor gases to be processed

    ! temperature variables
    integer,  dimension(ncol,nlay) :: jtemp ! interpolation index for temperature
    ! pressure variables
    integer,  dimension(ncol,nlay) :: jpress ! interpolation index for pressure
    logical,  dimension(ncol,nlay) :: tropo ! true for lower atmosphere; false for upper atmosphere
    integer, dimension(ncol,2) :: itropo_lower ! layer boundaries of lower atmosphere
    integer, dimension(ncol,2) :: itropo_upper ! layer boundaries of upper atmosphere
    integer, dimension(2,     nflav,ncol,nlay) :: jeta ! interpolation index for binary species parameter (eta)
                                                     ! index(1) : reference temperature level
                                                     ! index(2) : flavor
                                                     ! index(3) : layer

    real(wp), dimension(2,    nflav,ncol,nlay) :: col_mix ! combination of major species's column amounts
                                                         ! index(1) : reference temperature level
                                                         ! index(2) : flavor
                                                         ! index(3) : layer

    real(wp), dimension(2,2,2,nflav,ncol,nlay) :: fmajor ! interpolation fractions for major species
                                                            ! index(1) : reference eta level (temperature dependent)
                                                            ! index(2) : reference pressure level
                                                            ! index(3) : reference temperature level
                                                            ! index(4) : flavor
                                                            ! index(5) : layer

    real(wp), dimension(2,2,  nflav,ncol,nlay) :: fminor ! interpolation fractions for minor species
                                                          ! index(1) : reference eta level (temperature dependent)
                                                          ! index(2) : reference temperature level
                                                          ! index(3) : flavor
                                                          ! index(4) : layer


    ! ----------------------------------------------------------
    ! Code starts
    !

    error_msg = ''
    ! set up minor gases
    idx_h2o = string_loc_in_array('h2o', this%gas_names)
    do igas = 1, size(minor_gas_list)
      idx_gas_list(igas) = string_loc_in_array(minor_gas_list(igas), this%gas_names)
    end do

    ! compute column gas amounts
    col_gas(:,:,0) = col_dry(:,:)
    do igas = 1, ngas
      col_gas(:,:,igas) = vmr(:,:,igas) * col_dry(:,:)
    end do

    gpt_flv_lower = this%gpoint_flavor(1,:)
    gpt_flv_upper = this%gpoint_flavor(2,:)
    tau(:,:,:) = 0._wp
    ! ---- calculate gas optical depths ----
    call interpolation( &
      ncol,nlay,nflav,this%get_neta(), & ! dimensions
      this%flavor,this%press_ref_log,this%temp_ref,this%press_ref_log_delta,this%temp_ref_min, & ! inputs from object
      this%temp_ref_delta, this%press_ref_trop_log,this%vmr_ref,this%get_nlay_ref(), &
      play,tlay,col_gas, & ! local input
      jtemp,fmajor,fminor,col_mix,tropo,itropo_lower,itropo_upper,jeta,jpress) ! output
    call gas_optical_depths_major( & ! optical depths from major abosrbing gases
      ncol,nlay,ngpt,nflav, & ! dimensions
      this%gpoint_flavor,this%kmajor, & ! inputs from object
      col_mix,fmajor,&
      jeta,tropo,jtemp,jpress, & ! local input
      tau)
    call gas_optical_depths_minor( & !optical depths from minor gases in lower atmosphere, includes h2o continuum
      ncol,nlay,ngpt,ngas,nflav, & ! dimensions
      idx_h2o,&
      gpt_flv_lower, & ! inputs from object
      this%gas_names, &
      this%gas_minor, this%identifier_minor, &
      this%kminor_lower, &
      this%minor_gases_lower, &
      this%minor_limits_gpt_lower, &
      this%minor_scales_with_density_lower, &
      this%scaling_gas_lower, &
      this%scale_by_complement_lower, &
      this%kminor_start_lower, &
      play,tlay, &
      col_gas,idx_gas_list, &
      fminor,jeta,itropo_lower,jtemp, & ! local input
      tau)
    call gas_optical_depths_minor( & !optical depths from minor gases in upper atmosphere, includes h2o continuum
      ncol,nlay,ngpt,ngas,nflav, & ! dimensions
      idx_h2o,&
      gpt_flv_upper, & ! inputs from object
      this%gas_names, &
      this%gas_minor, this%identifier_minor, &
      this%kminor_upper, &
      this%minor_gases_upper, &
      this%minor_limits_gpt_upper, &
      this%minor_scales_with_density_upper, &
      this%scaling_gas_upper, &
      this%scale_by_complement_upper, &
      this%kminor_start_upper, &
      play,tlay, &
      col_gas,idx_gas_list, &
      fminor,jeta,itropo_upper,jtemp, & ! local input
      tau)
    if (allocated(this%krayl)) then
      call gas_optical_depths_rayleigh( & !Rayleigh scattering optical depths
        ncol,nlay,ngpt,ngas,nflav, & ! dimensions
        this%gpoint_flavor,this%krayl, & ! inputs from object
        idx_h2o,play,tlay,col_dry,col_gas,&
        fminor,jeta,tropo,jtemp, & ! local input
        tau_rayleigh)
    end if

    ! This is an internal function -- we can assume that all or none of these are present
    if(present(fmajor_out)) then
      fmajor_out = fmajor
      jeta_out   = jeta
      tropo_out  = tropo
      jtemp_out  = jtemp
      jpress_out = jpress
    end if

  end function compute_gas_tau_core
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Initialization
  !
  !--------------------------------------------------------------------------------------------------------------------
  ! Initialize object based on data read from netCDF file however the user desires.
  !  Rayleigh scattering tables may or may not be present; this is indicated with allocation status
  ! This interface is for the internal-sources object -- includes Plank functions and fractions
  !
  function init_int(this, available_gases, gas_names, key_species, &
                    band2gpt, band_lims_wavenum,            &
                    press_ref, press_ref_trop, temp_ref, &
                    temp_ref_p, temp_ref_t, vmr_ref,     &
                    kmajor, kminor_lower, kminor_upper, &
                    gas_minor,identifier_minor,&
                    minor_gases_lower, minor_gases_upper, &
                    minor_limits_gpt_lower, minor_limits_gpt_upper, &
                    minor_scales_with_density_lower, &
                    minor_scales_with_density_upper, &
                    scaling_gas_lower, scaling_gas_upper, &
                    scale_by_complement_lower, &
                    scale_by_complement_upper, &
                    kminor_start_lower, &
                    kminor_start_upper, &
                    totplnk, planck_frac, rayl_lower, rayl_upper) result(err_message)
    class(ty_gas_optics_specification), intent(inout) :: this
    class(ty_gas_concs),                intent(in   ) :: available_gases ! Which gases does the host model have available?
    character(len=*), dimension(:), intent(in) :: gas_names
    integer,  dimension(:,:,:),   intent(in) :: key_species
    integer,  dimension(:,:),     intent(in) :: band2gpt
    real(wp), dimension(:,:),     intent(in) :: band_lims_wavenum
    real(wp), dimension(:),       intent(in) :: press_ref, temp_ref
    real(wp),                     intent(in) :: press_ref_trop, temp_ref_p, temp_ref_t
    real(wp), dimension(:,:,:),   intent(in) :: vmr_ref
    real(wp), dimension(:,:,:,:), intent(in) :: kmajor
    real(wp), dimension(:,:,:),   intent(in) :: kminor_lower, kminor_upper
    real(wp), dimension(:,:),     intent(in) :: totplnk
    real(wp), dimension(:,:,:,:), intent(in) :: planck_frac
    real(wp), dimension(:,:,:),   intent(in), &
                                 allocatable :: rayl_lower, rayl_upper
    character(len=256), dimension(:), &
                                  intent(in) :: gas_minor,identifier_minor
    character(len=256), dimension(:), &
                                  intent(in) :: minor_gases_lower, &
                                                minor_gases_upper
    integer,  dimension(:,:),     intent(in) :: &
                                                minor_limits_gpt_lower, &
                                                minor_limits_gpt_upper
    logical,  dimension(:),       intent(in) :: &
                                                minor_scales_with_density_lower, &
                                                minor_scales_with_density_upper
    character(len=256), dimension(:),intent(in) :: &
                                                scaling_gas_lower, &
                                                scaling_gas_upper

    logical, dimension(:), intent(in) :: &
                                                scale_by_complement_lower,&
                                                scale_by_complement_upper
    integer, dimension(:), intent(in) :: &
                                                kminor_start_lower,&
                                                kminor_start_upper
    character(len = 128) err_message
    ! ----
    err_message = init_abs_coeffs(this, &
                                  gas_names, key_species,    &
                                  band2gpt, band_lims_wavenum, &
                                  press_ref, temp_ref,       &
                                  press_ref_trop, temp_ref_p, temp_ref_t, &
                                  vmr_ref,                   &
                                  kmajor, kminor_lower, kminor_upper, &
                                  gas_minor,identifier_minor,&
                                  minor_gases_lower, minor_gases_upper, &
                                  minor_limits_gpt_lower, &
                                  minor_limits_gpt_upper, &
                                  minor_scales_with_density_lower, &
                                  minor_scales_with_density_upper, &
                                  scaling_gas_lower, scaling_gas_upper, &
                                  scale_by_complement_lower, &
                                  scale_by_complement_upper, &
                                  kminor_start_lower, &
                                  kminor_start_upper, &
                                  rayl_lower, rayl_upper)
    ! Planck function tables
    !
    this%totplnk = totplnk
    this%planck_frac = planck_frac
    ! Temperature steps for Planck function interpolation
    !   Assumes that temperature minimum and max are the same for the absorption coefficient grid and the
    !   Planck grid and the Planck grid is equally spaced
    this%totplnk_delta =  (this%temp_ref_max-this%temp_ref_min) / (size(this%totplnk,dim=1)-1)
  end function init_int

  !--------------------------------------------------------------------------------------------------------------------
  ! Initialize object based on data read from netCDF file however the user desires.
  !  Rayleigh scattering tables may or may not be present; this is indicated with allocation status
  ! This interface is for the external-sources object -- includes TOA source function table
  !
  function init_ext(this, available_gases, gas_names, key_species,        &
                    band2gpt, band_lims_wavenum,           &
                    press_ref, press_ref_trop, temp_ref, &
                    temp_ref_p, temp_ref_t, vmr_ref,     &
                    kmajor, kminor_lower, kminor_upper, &
                    gas_minor,identifier_minor, &
                    minor_gases_lower, minor_gases_upper, &
                    minor_limits_gpt_lower, minor_limits_gpt_upper, &
                    minor_scales_with_density_lower, &
                    minor_scales_with_density_upper, &
                    scaling_gas_lower, scaling_gas_upper, &
                    scale_by_complement_lower, &
                    scale_by_complement_upper, &
                    kminor_start_lower, &
                    kminor_start_upper, &
                    solar_src, rayl_lower, rayl_upper)  result(err_message)
    class(ty_gas_optics_specification), intent(inout) :: this
    class(ty_gas_concs),                intent(in   ) :: available_gases ! Which gases does the host model have available?
    character(len=*), &
              dimension(:),       intent(in) :: gas_names
    integer,  dimension(:,:,:),   intent(in) :: key_species
    integer,  dimension(:,:),     intent(in) :: band2gpt
    real(wp), dimension(:,:),     intent(in) :: band_lims_wavenum
    real(wp), dimension(:),       intent(in) :: press_ref, temp_ref
    real(wp),                     intent(in) :: press_ref_trop, temp_ref_p, temp_ref_t
    real(wp), dimension(:,:,:),   intent(in) :: vmr_ref
    real(wp), dimension(:,:,:,:), intent(in) :: kmajor
    real(wp), dimension(:,:,:),   intent(in) :: kminor_lower, kminor_upper
    character(len=256), dimension(:), &
                                  intent(in) :: gas_minor, &
                                                identifier_minor
    character(len=256), dimension(:), &
                                  intent(in) :: minor_gases_lower, &
                                                minor_gases_upper
    integer,  dimension(:,:),     intent(in) :: &
                                                minor_limits_gpt_lower, &
                                                minor_limits_gpt_upper
    logical,  dimension(:),       intent(in) :: &
                                                minor_scales_with_density_lower, &
                                                minor_scales_with_density_upper
    character(len=256), dimension(:),intent(in) :: &
                                                scaling_gas_lower, &
                                                scaling_gas_upper
    logical,  dimension(:),       intent(in) :: &
                                                scale_by_complement_lower, &
                                                scale_by_complement_upper
    integer,  dimension(:),       intent(in) :: &
                                                kminor_start_lower, &
                                                kminor_start_upper
    real(wp), dimension(:),       intent(in), allocatable :: solar_src
                                                            ! allocatable status to change when solar source is present in file
    real(wp), dimension(:,:,:), intent(in), allocatable :: rayl_lower, rayl_upper
    character(len = 128) err_message
    ! ----
    err_message = init_abs_coeffs(this, &
                                  gas_names, key_species,    &
                                  band2gpt, band_lims_wavenum, &
                                  press_ref, temp_ref,       &
                                  press_ref_trop, temp_ref_p, temp_ref_t, &
                                  vmr_ref,                   &
                                  kmajor, kminor_lower, kminor_upper, &
                                  gas_minor,identifier_minor, &
                                  minor_gases_lower, minor_gases_upper, &
                                  minor_limits_gpt_lower, &
                                  minor_limits_gpt_upper, &
                                  minor_scales_with_density_lower, &
                                  minor_scales_with_density_upper, &
                                  scaling_gas_lower, scaling_gas_upper, &
                                  scale_by_complement_lower, &
                                  scale_by_complement_upper, &
                                  kminor_start_lower, &
                                  kminor_start_upper, &
                                  rayl_lower, rayl_upper)
    !
    ! Solar source table init
    !
    this%solar_src = solar_src

  end function init_ext
  !--------------------------------------------------------------------------------------------------------------------
  ! Initialize absorption coefficient arrays,
  !   including Rayleigh scattering tables if provided (allocated)
  !
  function init_abs_coeffs(this, &
                           gas_names, key_species,    &
                           band2gpt, band_lims_wavenum, &
                           press_ref, temp_ref,       &
                           press_ref_trop, temp_ref_p, temp_ref_t, &
                           vmr_ref,                   &
                           kmajor, kminor_lower, kminor_upper, &
                           gas_minor,identifier_minor,&
                           minor_gases_lower, minor_gases_upper, &
                           minor_limits_gpt_lower, &
                           minor_limits_gpt_upper, &
                           minor_scales_with_density_lower, &
                           minor_scales_with_density_upper, &
                           scaling_gas_lower, scaling_gas_upper, &
                           scale_by_complement_lower, &
                           scale_by_complement_upper, &
                           kminor_start_lower, &
                           kminor_start_upper, &
                           rayl_lower, rayl_upper) result(err_message)
    class(ty_gas_optics_specification), intent(inout) :: this
    character(len=*), &
              dimension(:),       intent(in) :: gas_names
    integer,  dimension(:,:,:),   intent(in) :: key_species
    integer,  dimension(:,:),     intent(in) :: band2gpt
    real(wp), dimension(:,:),     intent(in) :: band_lims_wavenum
    real(wp), dimension(:),       intent(in) :: press_ref, temp_ref
    real(wp),                     intent(in) :: press_ref_trop, temp_ref_p, temp_ref_t
    real(wp), dimension(:,:,:),   intent(in) :: vmr_ref
    real(wp), dimension(:,:,:,:), intent(in) :: kmajor
    real(wp), dimension(:,:,:),   intent(in) :: kminor_lower, kminor_upper
    character(len=256), dimension(:), &
                                  intent(in) :: gas_minor, &
                                                identifier_minor
    character(len=256), dimension(:), &
                                  intent(in) :: minor_gases_lower, &
                                                minor_gases_upper
    integer,  dimension(:,:),     intent(in) :: &
                                                minor_limits_gpt_lower, &
                                                minor_limits_gpt_upper
    logical,  dimension(:),       intent(in) :: &
                                                minor_scales_with_density_lower, &
                                                minor_scales_with_density_upper
    character(len=256), dimension(:),intent(in) :: &
                                                scaling_gas_lower, &
                                                scaling_gas_upper
    logical,  dimension(:), intent(in) :: &
                                                scale_by_complement_lower, &
                                                scale_by_complement_upper
    integer,  dimension(:), intent(in) :: &
                                                kminor_start_lower, &
                                                kminor_start_upper
    real(wp), dimension(:,:,:),   intent(in), &
                                 allocatable :: rayl_lower, rayl_upper
    real(wp), dimension(:,:,:), allocatable  :: vmr_ref_tmp
    character(len=128)                       :: err_message
    ! --------------------------------------
    err_message = this%ty_spectral_disc%init(band2gpt, band_lims_wavenum)
    if(len_trim(err_message) /= 0) return
    ! Assignment
    !   includes allocation
    if (allocated(vmr_ref_tmp)) deallocate(vmr_ref_tmp)
    allocate(vmr_ref_tmp(size(vmr_ref,dim=1),0:size(vmr_ref,dim=2), &
      size(vmr_ref,dim=3)))
    vmr_ref_tmp = vmr_ref

    this%gas_names = gas_names
    this%press_ref = press_ref
    this%temp_ref  = temp_ref
    this%vmr_ref   = vmr_ref_tmp
    this%kmajor       = kmajor
    this%kminor_lower = kminor_lower
    this%kminor_upper = kminor_upper
    this%gas_minor = gas_minor
    this%identifier_minor = identifier_minor
    this%minor_gases_lower = minor_gases_lower
    this%minor_gases_upper = minor_gases_upper
    this%minor_limits_gpt_lower = minor_limits_gpt_lower
    this%minor_limits_gpt_upper = minor_limits_gpt_upper
    this%minor_scales_with_density_lower = minor_scales_with_density_lower
    this%minor_scales_with_density_upper = minor_scales_with_density_upper
    this%scaling_gas_lower = scaling_gas_lower
    this%scaling_gas_upper = scaling_gas_upper
    this%scale_by_complement_lower = scale_by_complement_lower
    this%scale_by_complement_upper = scale_by_complement_upper
    this%kminor_start_lower = kminor_start_lower
    this%kminor_start_upper = kminor_start_upper
    if(allocated(rayl_lower) .neqv. allocated(rayl_upper)) then
      err_message = "rayl_lower and rayl_upper must have the same allocation status"
      return
    end if
    if (allocated(rayl_lower)) then
      allocate(this%krayl(size(rayl_lower,dim=1),size(rayl_lower,dim=2),size(rayl_lower,dim=3),2))
      this%krayl(:,:,:,1) = rayl_lower
      this%krayl(:,:,:,2) = rayl_upper
    end if

    ! ---- post processing ----
    ! Incoming coefficients file has units of Pa
    this%press_ref(:) = this%press_ref(:)

    ! creates log reference pressure
    allocate(this%press_ref_log(size(this%press_ref)))
    this%press_ref_log(:) = log(this%press_ref(:))

    ! log scale of reference pressure
    this%press_ref_trop_log = log(press_ref_trop)

    ! create flavor list
    call create_flavor(key_species, this%flavor)
    ! create gpoint_flavor list
    call create_gpoint_flavor(key_species, this%get_gpoint_bands(), this%flavor, this%gpoint_flavor)

    ! minimum, maximum reference temperature, pressure -- assumes low-to-high ordering
    !   for T, high-to-low ordering for p
    this%temp_ref_min  = this%temp_ref (1)
    this%temp_ref_max  = this%temp_ref (size(this%temp_ref))
    this%press_ref_min = this%press_ref(size(this%press_ref))
    this%press_ref_max = this%press_ref(1)

    ! creates press_ref_log, temp_ref_delta
    this%press_ref_log_delta = (log(this%press_ref_min)-log(this%press_ref_max))/(size(this%press_ref)-1)
    this%temp_ref_delta      = (this%temp_ref_max-this%temp_ref_min)/(size(this%temp_ref)-1)

    !
    ! Which species are key in one or more bands?
    !   this%flavor is an index into this%gas_names
    !
    if (allocated(this%is_key)) deallocate(this%is_key) ! Shouldn't ever happen...
    allocate(this%is_key(this%get_ngas()))
    this%is_key(:) = .False.
    this%is_key(pack(this%flavor(:,:), this%flavor(:,:) /= 0)) = .true.

  end function init_abs_coeffs

  !------------------------------------------------------------------------------------------
  !
  ! Ensure that every key gas required by the k-distribution is
  !    present in the gas concentration object
  !
  function check_key_species_present(this, gas_desc) result(error_msg)
    class(ty_gas_optics_specification), intent(in) :: this
    class(ty_gas_concs),                intent(in) :: gas_desc
    character(len=128)                             :: error_msg

    ! Local variables
    character(len=32), dimension(count(this%is_key(:)  )) :: key_gas_names
    integer                                               :: igas
    ! --------------------------------------
    error_msg = ""
    key_gas_names = pack(this%gas_names, mask=this%is_key)
    do igas = 1, size(key_gas_names)
      if(.not. string_in_array(key_gas_names(igas), gas_desc%gas_name)) &
        error_msg = ' ' // trim(lower_case(key_gas_names(igas))) // trim(error_msg)
    end do
    if(len_trim(error_msg) > 0) error_msg = "gas_optics: required gases" // trim(error_msg) // " are not provided"

  end function check_key_species_present

  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Function to define names of key and minor gases to be used by gas_optics().
  ! The final list gases includes those that are defined in gas_optics_specification
  ! and are provided in ty_gas_concs.
  !
  function get_minor_list(this, gas_desc, ngas, names_spec)
    class(ty_gas_optics_specification), intent(in)       :: this
    class(ty_gas_concs), intent(in)                      :: gas_desc
    integer, intent(in)                                  :: ngas
    character(32), dimension(ngas), intent(in)           :: names_spec

    ! List of minor gases to be used in gas_optics()
    character(len=32), dimension(:), allocatable         :: get_minor_list
    ! Logical flag for minor species in specification (T = minor; F = not minor)
    logical, dimension(size(names_spec))                 :: gas_is_present
    integer                                              :: igas, icnt

    if (allocated(get_minor_list)) deallocate(get_minor_list)
    do igas = 1, this%get_ngas()
      gas_is_present(igas) = string_in_array(names_spec(igas), gas_desc%gas_name)
    end do
    icnt = count(gas_is_present)
    allocate(get_minor_list(icnt))
    get_minor_list(:) = pack(this%gas_names, mask=gas_is_present)


  end function get_minor_list

  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Inquiry functions
  !
  !--------------------------------------------------------------------------------------------------------------------

  ! return true if initialized for internal sources, false otherwise
  pure function is_internal_source_present(this)
    class(ty_gas_optics_specification), intent(in) :: this
    logical                                        :: is_internal_source_present
    is_internal_source_present = allocated(this%totplnk).and.allocated(this%planck_frac)
  end function is_internal_source_present
  !--------------------------------------------------------------------------------------------------------------------

  ! return true if initialized for external sources, false otherwise
  pure function is_external_source_present(this)
    class(ty_gas_optics_specification), intent(in) :: this
    logical                                        :: is_external_source_present
    is_external_source_present = allocated(this%solar_src)
  end function is_external_source_present

  !--------------------------------------------------------------------------------------------------------------------
  ! return the gas names
  pure function get_gases(this)
    class(ty_gas_optics_specification), intent(in) :: this
    character(32), dimension(this%get_ngas())     :: get_gases

    get_gases = this%gas_names
  end function get_gases
  !--------------------------------------------------------------------------------------------------------------------
  ! return the minimum pressure on the interpolation grids
  pure function get_press_ref_min(this)
    class(ty_gas_optics_specification), intent(in) :: this
    real(wp)                                       :: get_press_ref_min

    get_press_ref_min = this%press_ref_min
  end function get_press_ref_min

  !--------------------------------------------------------------------------------------------------------------------
  ! return the maximum pressure on the interpolation grids
  pure function get_press_ref_max(this)
    class(ty_gas_optics_specification), intent(in) :: this
    real(wp)                                       :: get_press_ref_max

    get_press_ref_max = this%press_ref_max
  end function get_press_ref_max

  !--------------------------------------------------------------------------------------------------------------------
  ! return the minimum temparature on the interpolation grids
  pure function get_temp_ref_min(this)
    class(ty_gas_optics_specification), intent(in) :: this
    real(wp)                                       :: get_temp_ref_min

    get_temp_ref_min = this%temp_ref_min
  end function get_temp_ref_min

  !--------------------------------------------------------------------------------------------------------------------
  ! return the maximum temparature on the interpolation grids
  pure function get_temp_ref_max(this)
    class(ty_gas_optics_specification), intent(in) :: this
    real(wp)                                       :: get_temp_ref_max

    get_temp_ref_max = this%temp_ref_max
  end function get_temp_ref_max


  !--------------------------------------------------------------------------------------------------------------------
  ! --- gas optical depth calculations
  !--------------------------------------------------------------------------------------------------------------------
  ! Utility function, provided for user convenience
  ! computes column amounts of dry air using hydrostatic equation
  function get_col_dry(vmr_h2o, plev, tlay, latitude) result(col_dry)
    ! input
    real(wp), dimension(:,:), intent(in) :: vmr_h2o  ! volume mixing ratio of all gases excluding water; (ncol,nlay)
    real(wp), dimension(:,:), intent(in) :: plev     ! Layer boundary pressures [Pa, mb] (ncol,nlay+1)
    real(wp), dimension(:,:), intent(in) :: tlay     ! Layer temperatures [K] (ncol,nlay)
    real(wp), dimension(:),   optional, &
                              intent(in) :: latitude ! Latitude [degrees] (ncol)
    ! output
    real(wp), dimension(size(tlay,dim=1),size(tlay,dim=2)) :: col_dry ! Column dry amount (ncol,nlay)
    ! ------------------------------------------------
    ! first and second term of Helmert formula
    real(wp), parameter :: helmert1 = 9.80665_wp
    real(wp), parameter :: helmert2 = 0.02586_wp
    ! local variables
    real(wp), dimension(size(tlay,dim=1)                 ) :: g0 ! (ncol)
    real(wp), dimension(size(tlay,dim=1),size(tlay,dim=2)) :: delta_plev ! (ncol,nlay)
    real(wp), dimension(size(tlay,dim=1),size(tlay,dim=2)) :: m_air ! average mass of air; (ncol,nlay)
    integer :: nlev, nlay
    ! ------------------------------------------------
    nlay = size(tlay, dim=2)
    nlev = size(plev, dim=2)

    if(present(latitude)) then
      g0(:) = helmert1 - helmert2 * cos(2.0_wp * pi * latitude(:) / 180.0_wp) ! acceleration due to gravity [m/s^2]
    else
      g0(:) = grav
    end if
    delta_plev(:,:) = abs(plev(:,1:nlev-1) - plev(:,2:nlev))

    ! Get average mass of air
    m_air(:,:) = (m_dry+m_h2o*vmr_h2o(:,:))/(1.+vmr_h2o(:,:))

    ! Hydrostatic equation
    col_dry(:,:) = 10._wp*delta_plev(:,:)*avogad/(1000._wp*m_air(:,:)*100._wp*spread(g0(:),dim=2,ncopies=nlay))
    col_dry(:,:) = col_dry(:,:)/(1._wp+vmr_h2o(:,:))
  end function get_col_dry
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Internal procedures
  !
  !--------------------------------------------------------------------------------------------------------------------
  pure function rewrite_key_species_pair(key_species_pair)
    ! (0,0) becomes (2,2) -- because absorption coefficients for these g-points will be 0.
    integer, dimension(2) :: rewrite_key_species_pair
    integer, dimension(2), intent(in) :: key_species_pair
    rewrite_key_species_pair = key_species_pair
    if (all(key_species_pair(:).eq.(/0,0/))) then
      rewrite_key_species_pair(:) = (/2,2/)
    end if
  end function

  ! ---------------------------------------------------------------------------------------
  ! true is key_species_pair exists in key_species_list
  pure function key_species_pair_exists(key_species_list, key_species_pair)
    logical :: key_species_pair_exists
    integer, dimension(:,:), intent(in) :: key_species_list
    integer, dimension(2), intent(in) :: key_species_pair
    integer :: i
    do i=1,size(key_species_list,dim=2)
      if (all(key_species_list(:,i).eq.key_species_pair(:))) then
        key_species_pair_exists = .true.
        return
      end if
    end do
    key_species_pair_exists = .false.
  end function key_species_pair_exists

  ! ---------------------------------------------------------------------------------------
  ! create flavor list --
  !   an unordered array of extent (2,:) containing all possible pairs of key species
  !   used in either upper or lower atmos
  !
  subroutine create_flavor(key_species, flavor)
    integer, dimension(:,:,:), intent(in) :: key_species
    integer, dimension(:,:), allocatable, intent(out) :: flavor
    integer, dimension(2,size(key_species,3)*2) :: key_species_list

    integer :: ibnd, iatm, i, iflavor
    ! prepare list of key_species
    i = 1
    do ibnd=1,size(key_species,3)
      do iatm=1,size(key_species,1)
        key_species_list(:,i) = key_species(:,iatm,ibnd)
        i = i + 1
      end do
    end do
    ! rewrite single key_species pairs
    do i=1,size(key_species_list,2)
        key_species_list(:,i) = rewrite_key_species_pair(key_species_list(:,i))
    end do
    ! count unique key species pairs
    iflavor = 0
    do i=1,size(key_species_list,2)
      if (.not.key_species_pair_exists(key_species_list(:,1:i-1),key_species_list(:,i))) then
        iflavor = iflavor + 1
      end if
    end do
    ! fill flavors
    allocate(flavor(2,iflavor))
    iflavor = 0
    do i=1,size(key_species_list,2)
      if (.not.key_species_pair_exists(key_species_list(:,1:i-1),key_species_list(:,i))) then
        iflavor = iflavor + 1
        flavor(:,iflavor) = key_species_list(:,i)
      end if
    end do
  end subroutine create_flavor
! ---------------------------------------------------------------------------------------

  ! returns flavor index; -1 if not found
  pure function key_species_pair2flavor(flavor, key_species_pair)
    integer :: key_species_pair2flavor
    integer, dimension(:,:), intent(in) :: flavor
    integer, dimension(2), intent(in) :: key_species_pair
    integer :: iflav
    do iflav=1,size(flavor,2)
      if (all(key_species_pair(:).eq.flavor(:,iflav))) then
        key_species_pair2flavor = iflav
        return
      end if
    end do
    key_species_pair2flavor = -1
  end function key_species_pair2flavor

  ! ---------------------------------------------------------------------------------------
  ! create gpoint_flavor list
  !   a map pointing from each g-point to the corresponding entry in the "flavor list"
  !
  subroutine create_gpoint_flavor(key_species, gpt2band, flavor, gpoint_flavor)
    integer, dimension(:,:,:), intent(in) :: key_species
    integer, dimension(:), intent(in) :: gpt2band
    integer, dimension(:,:), intent(in) :: flavor
    integer, dimension(:,:), intent(out), allocatable :: gpoint_flavor
    integer :: ngpt, igpt, iatm
    ngpt = size(gpt2band)
    allocate(gpoint_flavor(2,ngpt))
    do igpt=1,ngpt
      do iatm=1,2
        gpoint_flavor(iatm,igpt) = key_species_pair2flavor( &
          flavor, &
          rewrite_key_species_pair(key_species(:,iatm,gpt2band(igpt))) &
        )
      end do
    end do
  end subroutine create_gpoint_flavor

 !--------------------------------------------------------------------------------------------------------------------
 !
 ! Utility function to combine optical depths from gas absorption and Rayleigh scattering
 !   (and reorder them for convenience, while we're at it)
 !
 subroutine combine_and_reorder(tau, tau_rayleigh, has_rayleigh, optical_props)
    real(wp), dimension(:,:,:),   intent(in) :: tau
    real(wp), dimension(:,:,:),   intent(in) :: tau_rayleigh
    logical,                      intent(in) :: has_rayleigh
    class(ty_optical_props_arry), intent(inout) :: optical_props

    integer :: icol, ilay, igpt, ncol, nlay, ngpt

    ncol = size(tau, 3)
    nlay = size(tau, 2)
    ngpt = size(tau, 1)

    if (.not. has_rayleigh) then
      ! index reorder (ngpt, nlay, ncol) -> (ncol,nlay,gpt)
      optical_props%tau = reorder123x321(tau)
    else
      ! combine optical depth and rayleigh scattering
      select type(optical_props)
        type is (ty_optical_props_1scl)
          ! User is asking for absorption optical depth
          optical_props%tau = reorder123x321(tau)
        type is (ty_optical_props_2str)
          do icol = 1, ncol
            do ilay = 1, nlay
              do igpt = 1, ngpt
                optical_props%tau(icol,ilay,igpt) = tau(igpt,ilay,icol) + tau_rayleigh(igpt,ilay,icol)
                optical_props%ssa(icol,ilay,igpt) = tau_rayleigh(igpt,ilay,icol) / optical_props%tau(icol,ilay,igpt)
              end do
            end do
          end do
          optical_props%g = 0._wp
        type is (ty_optical_props_nstr) ! We ought to be able to combine this with above
          do icol = 1, ncol
            do ilay = 1, nlay
              do igpt = 1, ngpt
                optical_props%tau(icol,ilay,igpt) = tau(igpt,ilay,icol) + tau_rayleigh(igpt,ilay,icol)
                optical_props%ssa(icol,ilay,igpt) = tau_rayleigh(igpt,ilay,icol) / optical_props%tau(icol,ilay,igpt)
              end do
            end do
          end do
          optical_props%p = 0._wp
          optical_props%p(2,:,:,:) = 0.1_wp
      end select

    end if
  end subroutine combine_and_reorder

  !--------------------------------------------------------------------------------------------------------------------
  ! return the number of reference pressure layers
  pure function get_nlay_ref(this)
    class(ty_gas_optics_specification), intent(in) :: this
    integer                                        :: get_nlay_ref

    get_nlay_ref = size(this%kmajor,dim=3)
  end function get_nlay_ref

  !--------------------------------------------------------------------------------------------------------------------
  ! return eta dimension
  pure function get_neta(this)
    class(ty_gas_optics_specification), intent(in) :: this
    integer                                        :: get_neta

    get_neta = size(this%kmajor,dim=2)
  end function

  !--------------------------------------------------------------------------------------------------------------------
  ! Generic procedures for checking sizes, limits
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Extents
  !
  function check_extent_1d(array, nx, label)
    real(wp), dimension(:),     intent(in) :: array
    integer,                    intent(in) :: nx
    character(len=*),           intent(in) :: label
    character(len=128)                     :: check_extent_1d

    check_extent_1d = ""
    if(size(array,1) /= nx) &
      check_extent_1d = trim(label) // ' has incorrect size.'
  end function check_extent_1d
  ! --------------------------------------------------------------------------------------
  function check_extent_2d(array, nx, ny, label)
    real(wp), dimension(:,:),   intent(in) :: array
    integer,                    intent(in) :: nx, ny
    character(len=*),           intent(in) :: label
    character(len=128)                     :: check_extent_2d

    check_extent_2d = ""
    if(size(array,1) /= nx .or. size(array,2) /= ny) &
      check_extent_2d = trim(label) // ' has incorrect size.'
  end function check_extent_2d
  ! --------------------------------------------------------------------------------------
  function check_extent_3d(array, nx, ny, nz, label)
    real(wp), dimension(:,:,:), intent(in) :: array
    integer,                    intent(in) :: nx, ny, nz
    character(len=*),           intent(in) :: label
    character(len=128)                     :: check_extent_3d

    check_extent_3d = ""
    if(size(array,1) /= nx .or. size(array,2) /= ny .or. size(array,3) /= nz) &
      check_extent_3d = trim(label) // ' has incorrect size.'
  end function check_extent_3d
  ! --------------------------------------------------------------------------------------
  !
  ! Values
  !
  ! --------------------------------------------------------------------------------------
  function check_range_1D(val, minV, maxV, label)
    real(wp), dimension(:),     intent(in) :: val
    real(wp),                   intent(in) :: minV, maxV
    character(len=*),           intent(in) :: label
    character(len=128)                     :: check_range_1D

    check_range_1D = ""
    if(any(val < minV) .or. any(val > maxV)) &
      check_range_1D = trim(label) // ' values out of range.'
  end function check_range_1D
  ! --------------------------------------------------------------------------------------
  function check_range_2D(val, minV, maxV, label)
    real(wp), dimension(:,:),   intent(in) :: val
    real(wp),                   intent(in) :: minV, maxV
    character(len=*),           intent(in) :: label
    character(len=128)                     :: check_range_2D

    check_range_2D = ""
    if(any(val < minV) .or. any(val > maxV)) &
      check_range_2D = trim(label) // ' values out of range.'
  end function check_range_2D
  ! --------------------------------------------------------------------------------------
  function check_range_3D(val, minV, maxV, label)
    real(wp), dimension(:,:,:), intent(in) :: val
    real(wp),                   intent(in) :: minV, maxV
    character(len=*),           intent(in) :: label
    character(len=128)                     :: check_range_3D

    check_range_3D = ""
    if(any(val < minV) .or. any(val > maxV)) &
      check_range_3D = trim(label) // ' values out of range.'
  end function check_range_3D
  !------------------------------------------------------------------------------------------


end module mo_gas_optics
