! This code is part of Radiative Transfer for Energetics (RTE)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2015-2018,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!
! Description:  Sets up arrays needed for optical properties.

module mo_optical_props
  use mo_rte_kind,              only: wp
  use mo_optical_props_kernels, only: &
        increment_1scalar_by_1scalar, increment_1scalar_by_2stream, increment_1scalar_by_nstream, &
        increment_2stream_by_1scalar, increment_2stream_by_2stream, increment_2stream_by_nstream, &
        increment_nstream_by_1scalar, increment_nstream_by_2stream, increment_nstream_by_nstream, &
        inc_1scalar_by_1scalar_bybnd, inc_1scalar_by_2stream_bybnd, inc_1scalar_by_nstream_bybnd, &
        inc_2stream_by_1scalar_bybnd, inc_2stream_by_2stream_bybnd, inc_2stream_by_nstream_bybnd, &
        inc_nstream_by_1scalar_bybnd, inc_nstream_by_2stream_bybnd, inc_nstream_by_nstream_bybnd, &
        delta_scale_2str_kernel, &
        any_vals_less_than, any_vals_outside
  implicit none
  integer, parameter :: name_len = 32
  !----------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------
  !
  ! Base class for optical properties
  !   Describes the spectral discretization including the wavenumber limits
  !   of each band (spectral region) and the mapping between g-points and bands
  !
  type, public :: ty_optical_props
    integer,  dimension(:,:), allocatable :: band2gpt       ! (begin g-point, end g-point) = band2gpt(2,band)
    integer,  dimension(:),   allocatable :: gpt2band       ! band = gpt2band(g-point)
    real(wp), dimension(:,:), allocatable :: band_lims_wvn  ! (upper and lower wavenumber by band) = band_lims_wvn(2,band)
    character(len=name_len)               :: name = ""
  contains
    generic,   public  :: init => init_base, init_base_from_copy
    procedure, private :: init_base
    procedure, private :: init_base_from_copy
    procedure, public  :: is_initialized => is_initialized_base
    procedure, private :: is_initialized_base
    procedure, public  :: finalize => finalize_base
    procedure, private :: finalize_base
    procedure, public  :: get_nband
    procedure, public  :: get_ngpt
    procedure, public  :: get_gpoint_bands
    procedure, public  :: convert_band2gpt
    procedure, public  :: convert_gpt2band
    procedure, public  :: get_band_lims_gpoint
    procedure, public  :: get_band_lims_wavenumber
    procedure, public  :: get_band_lims_wavelength
    procedure, public  :: expand
    procedure, public  :: set_name
    procedure, public  :: get_name
  end type
  !----------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------
  !
  ! Optical properties as arrays, normally dimensioned ncol, nlay, ngpt/nbnd
  !   The abstract base class for arrays defines what procedures will be available
  !   The optical depth field is also part of the abstract base class, since
  !    any representation of values as arrays needs an optical depth field
  !
  type, extends(ty_optical_props), abstract, public :: ty_optical_props_arry
    real(wp), dimension(:,:,:), allocatable :: tau ! optical depth (ncol, nlay, ngpt)
  contains
    procedure, public  :: get_ncol
    procedure, public  :: get_nlay
    !
    ! Routines to increment another set of values
    !
    procedure, private :: increment_bygpoint
    procedure, private :: increment_byband
    generic,   public  :: increment => increment_bygpoint, increment_byband

    !
    ! Deferred procedures -- each must be implemented in each child class with
    !   arguments following the abstract interface (defined below)
    !
    procedure(validate_abstract),     deferred, public  :: validate
    procedure(delta_scale_abstract),  deferred, public  :: delta_scale
    procedure(subset_range_abstract), deferred, public  :: get_subset
  end type
  !
  ! Interfaces for the methods to be implemented
  !
  abstract interface
    !
    ! Validation function looks only at internal data
    !
    function validate_abstract(this) result(err_message)
      import ty_optical_props_arry
      class(ty_optical_props_arry),  intent(in) :: this
      character(len=128)  :: err_message
    end function validate_abstract

    !
    ! Delta-scaling
    !
    function delta_scale_abstract(this, for) result(err_message)
      import ty_optical_props_arry
      import wp
      class(ty_optical_props_arry),  intent(inout) :: this
      real(wp), dimension(:,:,:), optional, &
                                     intent(in   ) :: for
      ! Forward scattering fraction; g**2 if not provided
      character(len=128)  :: err_message
    end function delta_scale_abstract

    !
    ! Subsetting -- currently there are only routines with start col and count
    !
    function subset_range_abstract(full, start, n, subset) result(err_message)
      import ty_optical_props_arry
      class(ty_optical_props_arry), intent(inout) :: full
      integer,                      intent(in   ) :: start, n
      class(ty_optical_props_arry), intent(inout) :: subset
      character(128)                              :: err_message
    end function subset_range_abstract
  end interface
  !----------------------------------------------------------------------------------------

  !   ty_optical_props_arry  includes only (extinction) optical depth
  !   Class two-stream adds arrays for single scattering albedo ssa and
  !     asymmetry parameter needed in two-stream methods
  !   Class n-stream adds arrays for single scattering albedo ssa and
  !     phase function moments (index 1 = g) for use with discrete ordinate methods
  !
  type, extends(ty_optical_props_arry) :: ty_optical_props_1scl
  contains
    procedure, public  :: validate => validate_1scalar
    procedure, public  :: get_subset => subset_1scl_range
    procedure, public  :: delta_scale => delta_scale_1scl

    procedure, private :: alloc_only_1scl
    procedure, private :: init_and_alloc_1scl
    generic,   public  :: alloc_1scl => alloc_only_1scl, init_and_alloc_1scl
  end type

  ! --- 2 stream ------------------------------------------------------------------------
  type, extends(ty_optical_props_arry) :: ty_optical_props_2str
    real(wp), dimension(:,:,:), allocatable :: ssa ! single-scattering albedo (ncol, nlay, ngpt)
    real(wp), dimension(:,:,:), allocatable :: g   ! asymmetry parameter (ncol, nlay, ngpt)
  contains
    procedure, public  :: validate => validate_2stream
    procedure, public  :: get_subset => subset_2str_range
    procedure, public  :: delta_scale => delta_scale_2str

    procedure, private :: alloc_only_2str
    procedure, private :: init_and_alloc_2str
    generic,   public  :: alloc_2str => alloc_only_2str, init_and_alloc_2str
  end type

  ! --- n stream ------------------------------------------------------------------------
  type, extends(ty_optical_props_arry) :: ty_optical_props_nstr
    real(wp), dimension(:,:,:),   allocatable :: ssa ! single-scattering albedo (ncol, nlay, ngpt)
    real(wp), dimension(:,:,:,:), allocatable :: p   ! phase-function moments (nmom, ncol, nlay, ngpt)
  contains
    procedure, public :: validate => validate_nstream
    procedure, public :: get_subset => subset_nstr_range
    procedure, public :: delta_scale => delta_scale_nstr
    procedure, public :: get_nmom

    procedure, private :: alloc_only_nstr
    procedure, private :: init_and_alloc_nstr
    generic,   public  :: alloc_nstr => alloc_only_nstr, init_and_alloc_nstr
  end type

contains
  ! ------------------------------------------------------------------------------------------
  !
  !  Routines for the base class: initialization, validity checking, finalization
  !
  ! ------------------------------------------------------------------------------------------
  !
  ! Base class: Initialization
  !
  function init_base(this, band_lims_gpt, band_lims_wvn, name) result(err_message)
    class(ty_optical_props),    intent(inout) :: this
    integer,  dimension(:,:),   intent(in   ) :: band_lims_gpt
    real(wp), dimension(:,:),   intent(in   ) :: band_lims_wvn
    character(len=*), optional, intent(in   ) :: name
    character(len = 128)                      :: err_message

    integer :: iband
    ! -------------------------
    !
    ! Error checking -- are the arrays the size we expect, contain positive values?
    !
    err_message = ""
    if(size(band_lims_gpt,1) /= 2 .or. size(band_lims_wvn,1) /= 2) &
      err_message = "optical_props%init(): band_lims_gpt or band_lims_wvn 1st dim wrong size"
    if(size(band_lims_gpt,2) /=        size(band_lims_wvn,2)     ) &
      err_message = "optical_props%init(): band_lims_gpt or band_lims_wvn sized inconsistently"
    if(any(band_lims_gpt < 1 .or.      band_lims_wvn < 1) ) &
      err_message = "optical_props%init(): band_lims_gpt or band_lims_wvn have values < 1"
    if(len_trim(err_message) > 0) return

    ! Assignment
    allocate(this%band2gpt     (2,size(band_lims_gpt,2)), &
             this%band_lims_wvn(2,size(band_lims_gpt,2)))
    this%band2gpt      = band_lims_gpt
    this%band_lims_wvn = band_lims_wvn
    if(present(name)) this%name = trim(name)

    !
    ! Make a map between g-points and bands
    !   Efficient only when g-point indexes start at 1 and are contiguous.
    !
    allocate(this%gpt2band(maxval(band_lims_gpt)))
    do iband=1,size(band_lims_gpt,dim=2)
      this%gpt2band(band_lims_gpt(1,iband):band_lims_gpt(2,iband)) = iband
    end do
  end function init_base
  !--------------------------------------------------------------------------------------------------------------------
  function init_base_from_copy(this, spectral_desc) result(err_message)
    class(ty_optical_props),    intent(inout) :: this
    class(ty_optical_props),    intent(in   ) :: spectral_desc
    character(len = 128)                      :: err_message

    if(.not. spectral_desc%is_initialized()) then
      err_message = "optical_props%init(): can't initialize based on un-initialized input"
      return
    else
      err_message = this%init(spectral_desc%get_band_lims_gpoint(), &
                              spectral_desc%get_band_lims_wavenumber())
      call this%set_name(spectral_desc%get_name())
    end if
  end function init_base_from_copy
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Base class: return true if initialized, false otherwise
  !
  pure function is_initialized_base(this)
    class(ty_optical_props), intent(in) :: this
    logical                             :: is_initialized_base

    is_initialized_base = allocated(this%band2gpt)
  end function is_initialized_base
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Base class: finalize (deallocate memory)
  !
  subroutine finalize_base(this)
    class(ty_optical_props),    intent(inout) :: this

    if(allocated(this%band2gpt)) deallocate(this%band2gpt)
    if(allocated(this%gpt2band)) deallocate(this%gpt2band)
    if(allocated(this%band_lims_wvn)) &
                                 deallocate(this%band_lims_wvn)
    this%name = ""
  end subroutine finalize_base
  ! ------------------------------------------------------------------------------------------
  !
  !  Routines for array classes: initialization, allocation, and finalization
  !
  ! ------------------------------------------------------------------------------------------
  ! Straight allocation routines
  !
  ! --- 1 scalar ------------------------------------------------------------------------
  function alloc_only_1scl(this, ncol, nlay) result(err_message)
    class(ty_optical_props_1scl) :: this
    integer,          intent(in) :: ncol, nlay
    character(len=128)           :: err_message

    err_message = ""
    if(.not. this%is_initialized()) then
      err_message = "optical_props%alloc: spectral discretization hasn't been provided"
      return
    end if
    if(any([ncol, nlay, this%get_ngpt()] <= 0)) then
      err_message = "optical_props%alloc: must provide positive extents for ncol, nlay, ngpt"
    else
      if(allocated(this%tau)) deallocate(this%tau)
      allocate(this%tau(ncol,nlay,this%get_ngpt()))
    end if
  end function alloc_only_1scl

  ! --- 2 stream ------------------------------------------------------------------------
  function alloc_only_2str(this, ncol, nlay) result(err_message)
    class(ty_optical_props_2str)    :: this
    integer,             intent(in) :: ncol, nlay
    character(len=128)              :: err_message

    err_message = ""
    if(any([ncol, nlay, this%get_ngpt()] <= 0)) then
      err_message = "optical_props%alloc: must provide positive extents for ncol, nlay, ngpt"
    else
      if(allocated(this%tau)) deallocate(this%tau)
      allocate(this%tau(ncol,nlay,this%get_ngpt()))
    end if
    if(allocated(this%ssa)) deallocate(this%ssa)
    if(allocated(this%g  )) deallocate(this%g  )
    if(err_message == "") allocate(this%ssa(ncol,nlay,this%get_ngpt()), this%g(ncol,nlay,this%get_ngpt()))
  end function alloc_only_2str

  ! --- n stream ------------------------------------------------------------------------
  function alloc_only_nstr(this, nmom, ncol, nlay) result(err_message)
    class(ty_optical_props_nstr)    :: this
    integer,             intent(in) :: nmom ! number of moments
    integer,             intent(in) :: ncol, nlay
    character(len=128)              :: err_message

    err_message = ""
    if(any([ncol, nlay, this%get_ngpt()] <= 0)) then
      err_message = "optical_props%alloc: must provide positive extents for ncol, nlay, ngpt"
    else
      if(allocated(this%tau)) deallocate(this%tau)
      allocate(this%tau(ncol,nlay,this%get_ngpt()))
    end if
    if(allocated(this%ssa)) deallocate(this%ssa)
    if(allocated(this%p  )) deallocate(this%p  )
    if(err_message == "") allocate(this%ssa(ncol,nlay,this%get_ngpt()), this%p(nmom,ncol,nlay,this%get_ngpt()))
  end function alloc_only_nstr
  ! ------------------------------------------------------------------------------------------
  ! Combined allocation/initialization routines
  !
  ! ---------------------------------------------------------------------------
  function init_and_alloc_1scl(this, spectral_desc, ncol, nlay) result(err_message)
    class(ty_optical_props_1scl)             :: this
    class(ty_optical_props     ), intent(in) :: spectral_desc
    integer,                      intent(in) :: ncol, nlay
    character(len=128)                       :: err_message

    err_message = ""
    if(this%ty_optical_props%is_initialized()) call this%ty_optical_props%finalize()
    err_message = this%ty_optical_props%init(spectral_desc%get_band_lims_gpoint(), &
                                             spectral_desc%get_band_lims_wavenumber())
    if(err_message /= "") return
    err_message = this%alloc_1scl(ncol, nlay)
  end function init_and_alloc_1scl
  ! ---------------------------------------------------------------------------
  function init_and_alloc_2str(this, spectral_desc, ncol, nlay) result(err_message)
    class(ty_optical_props_2str)             :: this
    class(ty_optical_props     ), intent(in) :: spectral_desc
    integer,                      intent(in) :: ncol, nlay
    character(len=128)                       :: err_message

    err_message = ""
    if(this%ty_optical_props%is_initialized()) call this%ty_optical_props%finalize()
    err_message = this%ty_optical_props%init(spectral_desc%get_band_lims_gpoint(), &
                                             spectral_desc%get_band_lims_wavenumber())
    if(err_message /= "") return
    err_message = this%alloc_2str(ncol, nlay)
  end function init_and_alloc_2str
  ! ---------------------------------------------------------------------------
  function init_and_alloc_nstr(this, spectral_desc, nmom, ncol, nlay) result(err_message)
    class(ty_optical_props_nstr)             :: this
    class(ty_optical_props     ), intent(in) :: spectral_desc
    integer,                      intent(in) :: nmom, ncol, nlay
    character(len=128)                       :: err_message

    err_message = ""
    if(this%ty_optical_props%is_initialized()) call this%ty_optical_props%finalize()
    err_message = this%ty_optical_props%init(spectral_desc%get_band_lims_gpoint(), &
                                             spectral_desc%get_band_lims_wavenumber())
    if(err_message /= "") return
    err_message = this%alloc_nstr(nmom, ncol, nlay)
  end function init_and_alloc_nstr
  ! ------------------------------------------------------------------------------------------
  !
  !  Routines for array classes: delta-scaling, validation (ensuring all values can be used )
  !
  ! ------------------------------------------------------------------------------------------
  ! --- delta scaling
  ! ------------------------------------------------------------------------------------------
  function delta_scale_1scl(this, for) result(err_message)
    class(ty_optical_props_1scl), intent(inout) :: this
    real(wp), dimension(:,:,:), optional, &
                                  intent(in   ) :: for
    character(128)                              :: err_message
    !
    ! Nothing to do
    !
    err_message = ""
  end function delta_scale_1scl
  ! ------------------------------------------------------------------------------------------
  function delta_scale_2str(this, for) result(err_message)
    class(ty_optical_props_2str), intent(inout) :: this
    real(wp), dimension(:,:,:), optional, &
                                  intent(in   ) :: for
    ! Forward scattering fraction; g**2 if not provided
    character(128)                              :: err_message

    integer :: ncol, nlay, ngpt
    integer :: icol, ilay, igpt
    ! --------------------------------
    ncol = this%get_ncol()
    nlay = this%get_nlay()
    ngpt = this%get_ngpt()
    err_message = ""

    if(present(for)) then
      if(any([size(for, 1), size(for, 2), size(for, 3)] /= [ncol, nlay, ngpt])) then
        err_message = "delta_scale: dimension of 'for' don't match optical properties arrays"
        return
      end if
      if(any(for < 0._wp .or. for > 1._wp)) then
        err_message = "delta_scale: values of 'for' out of bounds [0,1]"
        return
      end if
      call delta_scale_2str_kernel(ncol, nlay, ngpt, this%tau, this%ssa, this%g, for)
    else
      call delta_scale_2str_kernel(ncol, nlay, ngpt, this%tau, this%ssa, this%g)
    end if

  end function delta_scale_2str
  ! ------------------------------------------------------------------------------------------
  function delta_scale_nstr(this, for) result(err_message)
    class(ty_optical_props_nstr), intent(inout) :: this
    real(wp), dimension(:,:,:), optional, &
                                 intent(in   ) :: for
    character(128)                             :: err_message

    err_message = 'delta_scale_nstr: Not yet implemented'
  end function delta_scale_nstr
  ! ------------------------------------------------------------------------------------------
  ! --- Validation
  ! ------------------------------------------------------------------------------------------
  function validate_1scalar(this) result(err_message)
    class(ty_optical_props_1scl), intent(in) :: this
    character(len=128)                       :: err_message

    err_message = ''
    if(.not. allocated(this%tau)) then
      err_message = "validate: tau not allocated/initialized"
      return
    end if
    if(any_vals_less_than(size(this%tau, 1), size(this%tau, 2), size(this%tau, 3), this%tau, 0._wp)) &
      err_message = "validate: tau values out of range"
    if(len_trim(err_message) > 0 .and. len_trim(this%get_name()) > 0) &
      err_message = trim(this%get_name()) // ': ' // trim(err_message)

  end function validate_1scalar

  ! ------------------------------------------------------------------------------------------
  function validate_2stream(this) result(err_message)
    class(ty_optical_props_2str), intent(in) :: this
    character(len=128)                       :: err_message

    integer :: varSizes(3)

    err_message = ''
    !
    ! Array allocation status, sizing
    !
    if(.not. all([allocated(this%tau), allocated(this%ssa), allocated(this%g)])) then
      err_message = "validate: arrays not allocated/initialized"
      return
    end if
    varSizes =   [size(this%tau, 1), size(this%tau, 2), size(this%tau, 3)]
    if(.not. all([size(this%ssa, 1), size(this%ssa, 2), size(this%ssa, 3)] == varSizes) .or. &
       .not. all([size(this%g,   1), size(this%g,   2), size(this%g,   3)] == varSizes))     &
    err_message = "validate: arrays not sized consistently"
    !
    ! Valid values
    !
    if(any_vals_less_than(varSizes(1), varSizes(2), varSizes(3), this%tau,  0._wp)) &
      err_message = "validate: tau values out of range"
    if(any_vals_outside  (varSizes(1), varSizes(2), varSizes(3), this%ssa,  0._wp, 1._wp)) &
      err_message = "validate: ssa values out of range"
    if(any_vals_outside  (varSizes(1), varSizes(2), varSizes(3), this%g  , -1._wp, 1._wp)) &
      err_message = "validate: g values out of range"

    if(len_trim(err_message) > 0 .and. len_trim(this%get_name()) > 0) &
      err_message = trim(this%get_name()) // ': ' // trim(err_message)

  end function validate_2stream

  ! ------------------------------------------------------------------------------------------
  function validate_nstream(this) result(err_message)
    class(ty_optical_props_nstr), intent(in) :: this
    character(len=128)                       :: err_message

    integer :: varSizes(3)

    err_message = ''
    !
    ! Array allocation status, sizing
    !
    if(.not. all([allocated(this%tau), allocated(this%ssa), allocated(this%p)])) then
      err_message = "validate: arrays not allocated/initialized"
      return
    end if
    varSizes =   [size(this%tau, 1), size(this%tau, 2), size(this%tau, 3)]
    if(.not. all([size(this%ssa, 1), size(this%ssa, 2), size(this%ssa, 3)] == varSizes) .or. &
       .not. all([size(this%p,   2), size(this%p,   3), size(this%p,   4)] == varSizes))     &
    err_message = "validate: arrays not sized consistently"
    !
    ! Valid values
    !
    if(any_vals_less_than(varSizes(1), varSizes(2), varSizes(3), this%tau,  0._wp)) &
      err_message = "validate: tau values out of range"
    if(any_vals_outside  (varSizes(1), varSizes(2), varSizes(3), this%ssa,  0._wp, 1._wp)) &
      err_message = "validate: ssa values out of range"
    if(any_vals_outside  (varSizes(1), varSizes(2), varSizes(3), this%p(2,:,:,:),  &
                                                                           -1._wp, 1._wp)) &
      err_message = "validate: p(2,:,:,:)  = g values out of range"

    if(len_trim(err_message) > 0 .and. len_trim(this%get_name()) > 0) &
        err_message = trim(this%get_name()) // ': ' // trim(err_message)
  end function validate_nstream

  ! ------------------------------------------------------------------------------------------
  !
  !  Routines for array classes: subsetting of optical properties arrays along x (col) direction
  !
  ! Allocate class, then arrays; copy. Could probably be more efficient if
  !   classes used pointers internally.
  !
  ! This set takes start position and number as scalars
  ! ------------------------------------------------------------------------------------------

  function subset_1scl_range(full, start, n, subset) result(err_message)
    class(ty_optical_props_1scl), intent(inout) :: full
    integer,                      intent(in   ) :: start, n
    class(ty_optical_props_arry), intent(inout) :: subset
    character(128)                              :: err_message

    integer :: nlay, ngpt, nmom

    err_message = ""
    if(.not. full%is_initialized()) then
      err_message = "optical_props%subset: Asking for a subset of uninitialized data"
      return
    end if
    nlay = full%get_nlay()
    ngpt = full%get_ngpt()
    if(start < 1 .or. start + n-1 > full%get_ncol()) &
       err_message = "optical_props%subset: Asking for columns outside range"
    if(err_message /= "") return

    if(subset%is_initialized()) call subset%finalize()
    err_message = subset%init(full%band2gpt, full%band_lims_wvn, full%name)
    ! Seems like the deallocation statements should be needed under Fortran 2003
    !   but Intel compiler doesn't run without them
    if(allocated(subset%tau)) deallocate(subset%tau)
    select type (subset)
      class is (ty_optical_props_1scl)
        err_message = subset%alloc_1scl(n, nlay)
        if(err_message /= "") return
      class is (ty_optical_props_2str)
        if(allocated(subset%ssa)) deallocate(subset%ssa)
        if(allocated(subset%g  )) deallocate(subset%g  )
        err_message = subset%alloc_2str(n, nlay)
        if(err_message /= "") return
        subset%ssa(1:n,:,:) = 0._wp
        subset%g  (1:n,:,:) = 0._wp
      class is (ty_optical_props_nstr)
        if(allocated(subset%ssa)) deallocate(subset%ssa)
        if(allocated(subset%p  )) then
          nmom = subset%get_nmom()
          deallocate(subset%p  )
        else
          nmom = 1
        end if
        err_message = subset%alloc_nstr(nmom, n, nlay)
        if(err_message /= "") return
        subset%ssa(1:n,:,:) = 0._wp
        subset%p(:,1:n,:,:) = 0._wp
    end select
    subset%tau(1:n,:,:) = full%tau(start:start+n-1,:,:)

  end function subset_1scl_range
  ! ------------------------------------------------------------------------------------------
  function subset_2str_range(full, start, n, subset) result(err_message)
    class(ty_optical_props_2str), intent(inout) :: full
    integer,                      intent(in   ) :: start, n
    class(ty_optical_props_arry), intent(inout) :: subset
    character(128)                              :: err_message

    integer :: nlay, ngpt, nmom

    err_message = ""
    if(.not. full%is_initialized()) then
      err_message = "optical_props%subset: Asking for a subset of uninitialized data"
      return
    end if
    nlay = full%get_nlay()
    ngpt = full%get_ngpt()
    if(start < 1 .or. start + n-1 > full%get_ncol()) &
       err_message = "optical_props%subset: Asking for columns outside range"
    if(err_message /= "") return

    if(subset%is_initialized()) call subset%finalize()
    err_message = subset%init(full%band2gpt, full%band_lims_wvn, full%name)
    select type (subset)
      class is (ty_optical_props_1scl)
        err_message = subset%alloc_1scl(n, nlay)
        if(err_message /= "") return
        subset%tau(1:n,:,:) = full%tau(start:start+n-1,:,:) * &
                     (1._wp - full%ssa(start:start+n-1,:,:))
      class is (ty_optical_props_2str)
        if(allocated(subset%ssa)) deallocate(subset%ssa)
        if(allocated(subset%g  )) deallocate(subset%g  )
        err_message = subset%alloc_2str(n, nlay)
        if(err_message /= "") return
        subset%tau(1:n,:,:) = full%tau(start:start+n-1,:,:)
        subset%ssa(1:n,:,:) = full%ssa(start:start+n-1,:,:)
        subset%g  (1:n,:,:) = full%g  (start:start+n-1,:,:)
      class is (ty_optical_props_nstr)
        if(allocated(subset%ssa)) deallocate(subset%ssa)
        if(allocated(subset%p  )) then
          nmom = subset%get_nmom()
          deallocate(subset%p  )
        else
          nmom = 1
        end if
        err_message = subset%alloc_nstr(nmom, n, nlay)
        if(err_message /= "") return
        subset%tau(1:n,:,:) = full%tau(start:start+n-1,:,:)
        subset%ssa(1:n,:,:) = full%ssa(start:start+n-1,:,:)
        subset%p(1,1:n,:,:) = full%g  (start:start+n-1,:,:)
        subset%p(2:,:, :,:) = 0._wp
    end select

  end function subset_2str_range
  ! ------------------------------------------------------------------------------------------
  function subset_nstr_range(full, start, n, subset) result(err_message)
    class(ty_optical_props_nstr), intent(inout) :: full
    integer,                      intent(in   ) :: start, n
    class(ty_optical_props_arry), intent(inout) :: subset
    character(128)                              :: err_message

    integer :: nlay, ngpt, nmom

    err_message = ""
    if(.not. full%is_initialized()) then
      err_message = "optical_props%subset: Asking for a subset of uninitialized data"
      return
    end if
    nlay = full%get_nlay()
    ngpt = full%get_ngpt()
    if(start < 1 .or. start + n-1 > full%get_ncol()) &
       err_message = "optical_props%subset: Asking for columns outside range"
    if(err_message /= "") return

    if(subset%is_initialized()) call subset%finalize()
    err_message = subset%init(full%band2gpt, full%band_lims_wvn, full%name)
    if(allocated(subset%tau)) deallocate(subset%tau)
    select type (subset)
      class is (ty_optical_props_1scl)
        err_message = subset%alloc_1scl(n, nlay)
        if(err_message /= "") return
        subset%tau(1:n,:,:) = full%tau(start:start+n-1,:,:) * &
                     (1._wp - full%ssa(start:start+n-1,:,:))
      class is (ty_optical_props_2str)
        if(allocated(subset%ssa)) deallocate(subset%ssa)
        if(allocated(subset%g  )) deallocate(subset%g  )
        err_message = subset%alloc_2str(n, nlay)
        if(err_message /= "") return
        subset%tau(1:n,:,:) = full%tau(start:start+n-1,:,:)
        subset%ssa(1:n,:,:) = full%ssa(start:start+n-1,:,:)
        subset%g  (1:n,:,:) = full%p(1,start:start+n-1,:,:)
      class is (ty_optical_props_nstr)
        if(allocated(subset%ssa)) deallocate(subset%ssa)
        if(allocated(subset%p  )) deallocate(subset%p  )
        err_message = subset%alloc_nstr(nmom, n, nlay)
        if(err_message /= "") return
        subset%ssa(1:n,:,:) = full%ssa(start:start+n-1,:,:)
        subset%p(:,1:n,:,:) = full%p(:,start:start+n-1,:,:)
    end select

  end function subset_nstr_range

  ! ------------------------------------------------------------------------------------------
  !
  !  Routines for array classes: incrementing
  !
  !   The increment function, called as op1%increment(op2) or op1%increment(op2, band_limits),
  !   adds the values of op1 to op2, changing op2 and leaving op1 untouched. This lets us make
  !   functions that implement the same behavior.
  !
  !  It might make sense to break this big function up in type-bound procedures for each
  !    variant of ty_optical_props_arry
  ! -----------------------------------------------------------------------------------------
  function increment_bygpoint(op1, op2)result(err_message)
    class(ty_optical_props_arry), intent(in   ) :: op1
    class(ty_optical_props_arry), intent(inout) :: op2
    character(128)                              :: err_message
    ! -----
    integer :: ncol, nlay, ngpt, nmom1
    ! -----
    err_message = ""
    ncol = op1%get_ncol()
    nlay = op1%get_nlay()
    ngpt = op1%get_ngpt()

    !
    ! Rudimentary error checking -- users are responsible for ensuring consistency of
    !   array sizes within an object of ty_optical props
    !
    if(any([op2%get_ncol(), op2%get_nlay(), op2%get_ngpt()] /= [ncol, nlay, ngpt])) then
      err_message = "ty_optical_props%increment: optical properties objects are inconsistently sized"
      return
    end if

    select type (op2)
      class is (ty_optical_props_arry)
        select type (op1)
         class is (ty_optical_props_arry)
           call increment_1scalar_by_1scalar(ncol, nlay, ngpt, &
                                             op2%tau,          &
                                             op1%tau)
         class is (ty_optical_props_2str)
           call increment_1scalar_by_2stream(ncol, nlay, ngpt, &
                                             op2%tau,          &
                                             op1%tau, op1%ssa)

         class is (ty_optical_props_nstr)
           call increment_1scalar_by_nstream(ncol, nlay, ngpt, &
                                             op2%tau,          &
                                             op1%tau, op1%ssa)
        end select

    class is (ty_optical_props_2str)
      select type (op1)
        class is (ty_optical_props_arry)
          call increment_2stream_by_1scalar(ncol, nlay, ngpt,   &
                                            op2%tau, op2%ssa,&
                                            op1%tau)
        class is (ty_optical_props_2str)
          call increment_2stream_by_2stream(ncol, nlay, ngpt,        &
                                            op2%tau, op2%ssa, op2%g, &
                                            op1%tau, op1%ssa, op1%g)
        class is (ty_optical_props_nstr)
          call increment_2stream_by_nstream(ncol, nlay, ngpt, op1%get_nmom(), &
                                            op2%tau, op2%ssa, op2%g, &
                                            op1%tau, op1%ssa, op1%p)
      end select

    class is (ty_optical_props_nstr)
      nmom1 = op2%get_nmom()
      select type (op1)
        class is (ty_optical_props_arry)
          call increment_nstream_by_1scalar(ncol, nlay, ngpt, &
                                            op2%tau, op2%ssa, &
                                            op1%tau)
        class is (ty_optical_props_2str)
          call increment_nstream_by_2stream(ncol, nlay, ngpt, nmom1, &
                                            op2%tau, op2%ssa, op2%p, &
                                            op1%tau, op1%ssa, op1%g)
        class is (ty_optical_props_nstr)
          call increment_nstream_by_nstream(ncol, nlay, ngpt, nmom1, op1%get_nmom(), &
                                            op2%tau, op2%ssa, op2%p, &
                                            op1%tau, op1%ssa, op1%p)
      end select
    end select
  end function increment_bygpoint
  ! -----------------------------------------------------------------------------------------
  function increment_byband(op1, op2, gpt_lims) result(err_message)
    class(ty_optical_props_arry), intent(in   ) :: op1
    class(ty_optical_props_arry), intent(inout) :: op2
    integer, dimension(:,:),      intent(in   ) :: gpt_lims  ! (begin g-point, end g-point) = gpt_lims(2,band)
    character(128)                              :: err_message

    ! -----
    integer :: ncol, nlay, ngpt, nbnd, nmom1
    ! -----
    ncol = op2%get_ncol()
    nlay = op2%get_nlay()
    ngpt = op2%get_ngpt()
    nbnd = size(gpt_lims, 2)

    err_message = ""
    !
    ! Rudimentary error checking -- users are responsible for ensuring consistency of
    !   array sizes within an object of ty_optical props
    !
    if(any([op1%get_ncol(), op1%get_nlay()] /= [ncol, nlay])) &
      err_message = "ty_optical_props%increment_by: optical properties objects are inconsistently sized"
    if(op1%get_ngpt() /= nbnd)                           &
      err_message = "ty_optical_props%increment_by: " // &
                    "number of bands not consistent between optical properties objects, g-point limits"
    if(minval(gpt_lims) < 1 .or. maxval(gpt_lims) > ngpt) &
      err_message = "ty_optical_props%increment_by: " // &
                    "band limits not consistent with number of gpoints"
    if(err_message /= "") return

    select type (op2)
      class is (ty_optical_props_arry)
        select type (op1)
          class is (ty_optical_props_arry)
            call inc_1scalar_by_1scalar_bybnd(ncol, nlay, ngpt, &
                                              op2%tau,          &
                                              op1%tau,          &
                                              nbnd, gpt_lims)
          class is (ty_optical_props_2str)
            call inc_1scalar_by_2stream_bybnd(ncol, nlay, ngpt, &
                                              op2%tau,          &
                                              op1%tau, op1%ssa, &
                                              nbnd, gpt_lims)
          class is (ty_optical_props_nstr)
            call inc_1scalar_by_nstream_bybnd(ncol, nlay, ngpt, &
                                              op2%tau,          &
                                              op1%tau, op1%ssa, &
                                              nbnd, gpt_lims)
        end select

      class is (ty_optical_props_2str)
        select type (op1)
          class is (ty_optical_props_arry)
            call inc_2stream_by_1scalar_bybnd(ncol, nlay, ngpt, &
                                              op2%tau, op2%ssa, &
                                              op1%tau,          &
                                              nbnd, gpt_lims)
          class is (ty_optical_props_2str)
            call inc_2stream_by_2stream_bybnd(ncol, nlay, ngpt,        &
                                              op2%tau, op2%ssa, op2%g, &
                                              op1%tau, op1%ssa, op1%g, &
                                              nbnd, gpt_lims)
          class is (ty_optical_props_nstr)
            call inc_2stream_by_nstream_bybnd(ncol, nlay, ngpt, op1%get_nmom(), &
                                              op2%tau, op2%ssa, op2%g, &
                                              op1%tau, op1%ssa, op1%p, &
                                              nbnd, gpt_lims)
        end select

      class is (ty_optical_props_nstr)
        nmom1 = op2%get_nmom()
        select type (op1)
          class is (ty_optical_props_arry)
            call inc_nstream_by_1scalar_bybnd(ncol, nlay, ngpt, &
                                              op2%tau, op2%ssa, &
                                              op1%tau,          &
                                              nbnd, gpt_lims)
          class is (ty_optical_props_2str)
            call inc_nstream_by_2stream_bybnd(ncol, nlay, ngpt, nmom1, &
                                              op2%tau, op2%ssa, op2%p, &
                                              op1%tau, op1%ssa, op1%g, &
                                              nbnd, gpt_lims)
          class is (ty_optical_props_nstr)
            call inc_nstream_by_nstream_bybnd(ncol, nlay, ngpt, nmom1, op1%get_nmom(), &
                                              op2%tau, op2%ssa, op2%p, &
                                              op1%tau, op1%ssa, op1%p, &
                                              nbnd, gpt_lims)
        end select
    end select
  end function increment_byband
  ! ------------------------------------------------------------------------------------------
  !
  !  Routines for array classes: problem sizes
  !
  ! ------------------------------------------------------------------------------------------
  pure function get_arry_extent(this, dim)
    class(ty_optical_props_arry), intent(in   ) :: this
    integer,                      intent(in   ) :: dim
    integer                                     :: get_arry_extent

    if(allocated(this%tau)) then
      get_arry_extent = size(this%tau, dim)
    else
      get_arry_extent = 0
    end if
  end function get_arry_extent
  ! ------------------------------------------------------------------------------------------
  pure function get_ncol(this)
    class(ty_optical_props_arry), intent(in   ) :: this
    integer                                     :: get_ncol

    get_ncol = get_arry_extent(this, 1)
  end function get_ncol
  ! ------------------------------------------------------------------------------------------
  pure function get_nlay(this)
    class(ty_optical_props_arry), intent(in   ) :: this
    integer                                     :: get_nlay

    get_nlay = get_arry_extent(this, 2)
  end function get_nlay
  ! ------------------------------------------------------------------------------------------
  pure function get_nmom(this)
    class(ty_optical_props_nstr), intent(in   ) :: this
    integer                                     :: get_nmom

    if(allocated(this%p)) then
      get_nmom = size(this%p, 1)
    else
      get_nmom = 0
    end if
  end function get_nmom
  ! ------------------------------------------------------------------------------------------
  !
  !  Routines for base class: spectral discretization
  !
  ! ------------------------------------------------------------------------------------------
  !
  ! Number of bands
  !
  pure function get_nband(this)
    class(ty_optical_props), intent(in) :: this
    integer                             :: get_nband

    if(this%is_initialized()) then
      get_nband = size(this%band2gpt,dim=2)
    else
      get_nband = 0
    end if
  end function get_nband

  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Number of g-points
  !
  pure function get_ngpt(this)
    class(ty_optical_props), intent(in) :: this
    integer                             :: get_ngpt

    if(this%is_initialized()) then
      get_ngpt = maxval(this%band2gpt)
    else
      get_ngpt = 0
    end if
  end function get_ngpt

  !--------------------------------------------------------------------------------------------------------------------
  !
  ! The first and last g-point of all bands at once
  ! dimension (2, nbands)
  !
  pure function get_band_lims_gpoint(this)
    class(ty_optical_props), intent(in) :: this
    integer, dimension(size(this%band2gpt,dim=1), size(this%band2gpt,dim=2)) &
                                        :: get_band_lims_gpoint

    get_band_lims_gpoint = this%band2gpt
  end function get_band_lims_gpoint

  !--------------------------------------------------------------------------------------------------------------------
  !
  ! First and last g-point of a specific band
  !
  pure function convert_band2gpt(this, band)
    class(ty_optical_props), intent(in) :: this
    integer,                 intent(in) :: band
    integer, dimension(2)               :: convert_band2gpt

    if(this%is_initialized()) then
      convert_band2gpt(:) = this%band2gpt(:,band)
    else
      convert_band2gpt(:) = 0
    end if
  end function convert_band2gpt

  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Lower and upper wavenumber of all bands
  ! (upper and lower wavenumber by band) = band_lims_wvn(2,band)
  !
  pure function get_band_lims_wavenumber(this)
    class(ty_optical_props), intent(in) :: this
    real(wp), dimension(size(this%band_lims_wvn,1), size(this%band_lims_wvn,2)) &
                                        :: get_band_lims_wavenumber

    if(this%is_initialized()) then
      get_band_lims_wavenumber(:,:) = this%band_lims_wvn(:,:)
    else
      get_band_lims_wavenumber(:,:) = 0._wp
    end if
  end function get_band_lims_wavenumber

  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Lower and upper wavelength of all bands
  !
  pure function get_band_lims_wavelength(this)
    class(ty_optical_props), intent(in) :: this
    real(wp), dimension(size(this%band_lims_wvn,1), size(this%band_lims_wvn,2)) &
                                        :: get_band_lims_wavelength

    if(this%is_initialized()) then
      get_band_lims_wavelength(:,:) = 1._wp/this%band_lims_wvn(:,:)
    else
      get_band_lims_wavelength(:,:) = 0._wp
    end if
  end function get_band_lims_wavelength

  !--------------------------------------------------------------------------------------------------------------------
  ! Bands for all the g-points at once
  ! dimension (ngpt)
  !
  pure function get_gpoint_bands(this)
    class(ty_optical_props), intent(in) :: this
    integer, dimension(size(this%gpt2band,dim=1)) &
                                        :: get_gpoint_bands

    if(this%is_initialized()) then
      get_gpoint_bands(:) = this%gpt2band(:)
    else
      get_gpoint_bands(:) = 0
    end if
  end function get_gpoint_bands

  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Band associated with a specific g-point
  !
  pure function convert_gpt2band(this, gpt)
    class(ty_optical_props), intent(in) :: this
    integer,                            intent(in) :: gpt
    integer                             :: convert_gpt2band

    if(this%is_initialized()) then
      convert_gpt2band = this%gpt2band(gpt)
    else
      convert_gpt2band = 0
    end if
  end function convert_gpt2band

  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Expand an array of dimension arr_in(nband) to dimension arr_out(ngpt)
  !
  pure function expand(this, arr_in) result(arr_out)
    class(ty_optical_props), intent(in) :: this
    real(wp), dimension(:),  intent(in) :: arr_in ! (nband)
    real(wp), dimension(size(this%gpt2band)) :: arr_out

    integer :: iband

    do iband=1,this%get_nband()
      arr_out(this%band2gpt(1,iband):this%band2gpt(2,iband)) = arr_in(iband)
    end do
  end function expand
  !--------------------------------------------------------------------------------------------------------------------
  ! --- Setting/getting the name
  ! ------------------------------------------------------------------------------------------
  subroutine set_name(this, name)
    class(ty_optical_props),  intent(inout) :: this
    character(len=*),         intent(in   ) :: name

    this%name = trim(name)
  end subroutine set_name
  ! --------------------------------------------------------
  function get_name(this)
    class(ty_optical_props),  intent(in   ) :: this
    character(len=name_len)                 :: get_name

      get_name = trim(this%name)
  end function get_name
  ! ------------------------------------------------------------------------------------------

end module mo_optical_props
