! This code is part of Radiative Transfer for Energetics (RTE)
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
! Description:  Describes spectral discretization for a k-distribution, including the wavenumber limits
!   of each band (spectral region) and the mapping between g-points and bands

module mo_spectral_disc
  use mo_rte_kind,        only: wp
  implicit none
  private

  ! -----------------------------------------------------------------------------------
  type, public :: ty_spectral_disc
    private
    integer,  dimension(:,:), allocatable :: band2gpt       ! (begin g-point, end g-point) = band2gpt(2,band)
    integer,  dimension(:),   allocatable :: gpt2band       ! band = gpt2band(g-point)
    real(wp), dimension(:,:), allocatable :: band_lims_wvn  ! (upper and lower wavenumber by band) = band_lims_wvn(2,band)
  contains
    generic,   public :: init => init_disc
    procedure, public :: init_disc
    generic,   public :: is_initialized => is_initialized_desc
    procedure, public :: is_initialized_desc
    procedure, public :: get_nband
    procedure, public :: get_ngpt
    procedure, public :: get_gpoint_bands
    procedure, public :: convert_band2gpt
    procedure, public :: convert_gpt2band
    procedure, public :: get_band_lims_gpoint
    procedure, public :: get_band_lims_wavenumber
    procedure, public :: get_band_lims_wavelength
    procedure, public :: expand
  end type ty_spectral_disc

contains
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Initialization
  !
  function init_disc(this, band_lims_gpt, band_lims_wvn) result(err_message)
    class(ty_spectral_disc),  intent(inout) :: this
    integer,  dimension(:,:), intent(in   ) :: band_lims_gpt
    real(wp), dimension(:,:), intent(in   ) :: band_lims_wvn
    character(len = 128)                    :: err_message

    integer :: iband
    ! -------------------------
    !
    ! Error checking -- are the arrays the size we expect, contain positive values?
    !
    err_message = ""
    if(size(band_lims_gpt,1) /= 2 .or. size(band_lims_wvn,1) /= 2) &
      err_message = "ty_spectral_disc%init(): band_lims_gpt or band_lims_wvn 1st dim wrong size"
    if(size(band_lims_gpt,2) /=        size(band_lims_wvn,2)     ) &
      err_message = "ty_spectral_disc%init(): band_lims_gpt or band_lims_wvn sized inconsistently"
    if(any(band_lims_gpt < 1 .or.      band_lims_wvn < 1) ) &
      err_message = "ty_spectral_disc%init(): band_lims_gpt or band_lims_wvn have values < 1"
    if(len_trim(err_message) > 0) return

    ! Assignment, includes allocation
    this%band2gpt      = band_lims_gpt
    this%band_lims_wvn = band_lims_wvn

    !
    ! Make a map between g-points and bands
    !   Efficient only when g-point indexes start at 1 and are contiguous.
    !
    allocate(this%gpt2band(maxval(band_lims_gpt)))
    do iband=1,size(band_lims_gpt,dim=2)
      this%gpt2band(band_lims_gpt(1,iband):band_lims_gpt(2,iband)) = iband
    end do
  end function init_disc
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! return true if initialized, false otherwise
  !
  pure function is_initialized_desc(this)
    class(ty_spectral_disc), intent(in) :: this
    logical                             :: is_initialized_desc

    is_initialized_desc = allocated(this%band2gpt)
  end function is_initialized_desc

  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Number of bands
  !
  pure function get_nband(this)
    class(ty_spectral_disc), intent(in) :: this
    integer                                        :: get_nband

    get_nband = size(this%band2gpt,dim=2)
  end function get_nband

  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Number of g-points
  !
  pure function get_ngpt(this)
    class(ty_spectral_disc), intent(in) :: this
    integer                                        :: get_ngpt

    get_ngpt = maxval(this%band2gpt)
  end function get_ngpt

  !--------------------------------------------------------------------------------------------------------------------
  !
  ! The first and last g-point of all bands at once
  ! dimension (2, nbands)
  !
  pure function get_band_lims_gpoint(this)
    class(ty_spectral_disc), intent(in) :: this
    integer, dimension(size(this%band2gpt,dim=1), size(this%band2gpt,dim=2)) &
                                                   :: get_band_lims_gpoint

    get_band_lims_gpoint = this%band2gpt
  end function get_band_lims_gpoint

  !--------------------------------------------------------------------------------------------------------------------
  !
  ! First and last g-point of a specific band
  !
  pure function convert_band2gpt(this, band)
    class(ty_spectral_disc), intent(in) :: this
    integer,                            intent(in) :: band
    integer, dimension(2)                          :: convert_band2gpt

    convert_band2gpt(:) = this%band2gpt(:,band)
  end function convert_band2gpt

  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Lower and upper wavenumber of all bands
  ! (upper and lower wavenumber by band) = band_lims_wvn(2,band)
  !
  pure function get_band_lims_wavenumber(this)
    class(ty_spectral_disc), intent(in) :: this
    real(wp), dimension(size(this%band_lims_wvn,1), size(this%band_lims_wvn,2)) &
                                                   :: get_band_lims_wavenumber

    get_band_lims_wavenumber(:,:) = this%band_lims_wvn(:,:)
  end function get_band_lims_wavenumber

  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Lower and upper wavelength of all bands
  !
  pure function get_band_lims_wavelength(this)
    class(ty_spectral_disc), intent(in) :: this
    real(wp), dimension(size(this%band_lims_wvn,1), size(this%band_lims_wvn,2)) &
                                                   :: get_band_lims_wavelength

    get_band_lims_wavelength(:,:) = 1._wp/this%band_lims_wvn(:,:)
  end function get_band_lims_wavelength

  !--------------------------------------------------------------------------------------------------------------------
  ! Bands for all the g-points at once
  ! dimension (ngpt)
  !
  pure function get_gpoint_bands(this)
    class(ty_spectral_disc), intent(in) :: this
    integer, dimension(size(this%gpt2band,dim=1)) &
                                                   :: get_gpoint_bands

    get_gpoint_bands(:) = this%gpt2band(:)
  end function get_gpoint_bands

  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Band associated with a specific g-point
  !
  pure function convert_gpt2band(this, gpt)
    class(ty_spectral_disc), intent(in) :: this
    integer,                            intent(in) :: gpt
    integer                                        :: convert_gpt2band

    convert_gpt2band = this%gpt2band(gpt)
  end function convert_gpt2band

  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Expand an array of dimension arr_in(nband) to dimension arr_out(ngpt)
  !
  pure function expand(this, arr_in) result(arr_out)
    class(ty_spectral_disc), intent(in) :: this
    real(wp), dimension(:), intent(in) :: arr_in ! (nband)
    real(wp), dimension(size(this%gpt2band)) :: arr_out

    integer :: iband

    do iband=1,this%get_nband()
      arr_out(this%band2gpt(1,iband):this%band2gpt(2,iband)) = arr_in(iband)
    end do
  end function expand
  !--------------------------------------------------------------------------------------------------------------------

end module mo_spectral_disc
