! Module: mo_optical_props

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
! Description:  Sets up arrays needed for optical properties.

module mo_optical_props
  use mo_rte_kind,           only: wp
  use mo_optical_props_kernels, only: &
        increment_1scalar_by_1scalar, increment_1scalar_by_2stream, increment_1scalar_by_nstream, &
        increment_2stream_by_1scalar, increment_2stream_by_2stream, increment_2stream_by_nstream, &
        increment_nstream_by_1scalar, increment_nstream_by_2stream, increment_nstream_by_nstream, & 
        inc_1scalar_by_1scalar_bybnd, inc_1scalar_by_2stream_bybnd, inc_1scalar_by_nstream_bybnd, &
        inc_2stream_by_1scalar_bybnd, inc_2stream_by_2stream_bybnd, inc_2stream_by_nstream_bybnd, &
        inc_nstream_by_1scalar_bybnd, inc_nstream_by_2stream_bybnd, inc_nstream_by_nstream_bybnd
  implicit none
  integer, parameter :: name_len = 32
  !----------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------
  !
  ! Abstract base class for optical properties 
  !   Variables of this type can't be declared directly. 
  !   Instead, the type needs to be extended and the procedures implemented.   
  !
  type, abstract, public :: ty_optical_props
    !
    ! Everyone needs a name
    !
    character(len=name_len)                 :: name = "" 
  contains
    !
    ! Procedures to set and retrieve the name 
    !
    procedure, public :: set_name
    procedure, public :: get_name
    
    !
    ! Deferred procedures. Concrete types need to implement procedures that follow the interfaces below 
    !
    ! Functions to report domain size 
    !   ngpt is the number of spectral points; for some instances it might be bands 
    ! 
    procedure(get_scalar_int_abstract), deferred, & 
               public :: get_ncol
    procedure(get_scalar_int_abstract), deferred, & 
               public :: get_nlay
    procedure(get_scalar_int_abstract), deferred, & 
               public :: get_ngpt 
    !
    ! Functions to add optical properties to existing arrays. 
    !   Concrete types need to implement incrementing both one-to-one (by gpoint) and 
    !   many-to-one (by band); these could do nothing but return error strings 
    !   
    procedure(increment_bygpt_abstract ), deferred, & 
               private :: increment_bygpt 
    procedure(increment_byband_abstract), deferred, & 
               private :: increment_byband  
    !
    !  Incrementing functions have two distinct interfaces. 
    !   This refers to them both with a single name (overloading) 
    ! 
    generic,   public  :: increment => increment_bygpt, increment_byband
  end type
  abstract interface
    !
    ! Abstract interface to get single scalar values -- used for ncol, nlay, ngpt 
    !
    function get_scalar_int_abstract(this) 
      import ty_optical_props
      class(ty_optical_props),  intent(in   ) :: this
      integer                                 :: get_scalar_int_abstract
    end function get_scalar_int_abstract     
  end interface 
  !----------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------
  !
  ! Optical properties as arrays, normally dimensioned ncol, nlay, ngpt/nbnd 
  !
  !   Abstract class ty_optical_props_arry implements procedures from the base class and 
  !   defines new procedures that concrete classes will implement. 
  ! 
  !   Any representation of values as arrays needs an optical depth field so that's defined 
  !   in the abstract base class too. 
  !
  type, extends(ty_optical_props), abstract, public :: ty_optical_props_arry
    real(wp), dimension(:,:,:), allocatable :: tau ! optical depth (ncol, nlay, ngpt)
  contains
    !
    ! Implementation of deferred procedures from abstract base class. 
    !
    procedure, private :: get_ncol_arry 
    procedure, private :: get_nlay_arry
    procedure, private :: get_ngpt_arry
    procedure, public  :: get_ncol => get_ncol_arry
    procedure, public  :: get_nlay => get_nlay_arry
    procedure, public  :: get_ngpt => get_ngpt_arry
    
    procedure, public  :: increment_bygpt  => increment_arry
    procedure, public  :: increment_byband => increment_arry_byband

    !   Should array classes implement init ? 
    
    !
    ! Deferred procedures -- each must be implemented in each child class with 
    !   arguments following the abstract interface (defined below)
    !
    procedure(validate_abstract),     deferred, public  :: validate
    procedure(delta_scale_abstract),  deferred, public  :: delta_scale
    procedure(subset_range_abstract), deferred, public  :: get_subset
    !
    ! increment_by routines; deprecated 
    !   op1%increment_by(op2) changes the values of op1 
    ! 
    procedure, private :: increment_gpt_by
    procedure, private :: increment_band_by
    generic,   public  :: increment_by => increment_gpt_by, increment_band_by
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
      real(wp), dimension(:,:,:), target, optional, &
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

    !
    ! Incrementing -- any object extending ty_optical_props need to be able to increment arrays. 
    !
    function increment_bygpt_abstract(op1, op2) result(err_message)
      import ty_optical_props, ty_optical_props_arry
      class(ty_optical_props),       intent(in   ) :: op1
      class(ty_optical_props_arry),  intent(inout) :: op2
      character(len=128)                           :: err_message
    end function increment_bygpt_abstract
    ! ------------------
    function increment_byband_abstract(op1, op2, gpt_lims) result(err_message)
      import ty_optical_props, ty_optical_props_arry
      class(ty_optical_props),       intent(in   ) :: op1
      class(ty_optical_props_arry),  intent(inout) :: op2
      integer, dimension(:,:),       intent(in   ) :: gpt_lims  ! (begin g-point, end g-point) = gpt_lims(2,band)
      character(len=128)                           :: err_message
    end function increment_byband_abstract

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

    procedure, public  :: init_1scl
  end type

  ! --- 2 stream ------------------------------------------------------------------------
  type, extends(ty_optical_props_arry) :: ty_optical_props_2str
    real(wp), dimension(:,:,:), allocatable :: ssa ! single-scattering albedo (ncol, nlay, ngpt)
    real(wp), dimension(:,:,:), allocatable :: g   ! asymmetry parameter (ncol, nlay, ngpt)
  contains
    procedure, public  :: validate => validate_2stream
    procedure, public  :: get_subset => subset_2str_range
    procedure, public  :: delta_scale => delta_scale_2str

    procedure, public  :: init_2str
  end type

  ! --- n stream ------------------------------------------------------------------------
  type, extends(ty_optical_props_arry) :: ty_optical_props_nstr
    real(wp), dimension(:,:,:), allocatable :: ssa ! single-scattering albedo (ncol, nlay, ngpt)
    real(wp), dimension(:,:,:,:), allocatable :: p ! phase-function moments (nmom, ncol, nlay, ngpt)
  contains
    procedure, public  :: validate => validate_nstream
    procedure, public  :: get_subset => subset_nstr_range
    procedure, public  :: delta_scale => delta_scale_nstr

    procedure, public  :: init_nstr
  end type

contains
  ! ------------------------------------------------------------------------------------------
  ! --- Initialization routines 
  ! ------------------------------------------------------------------------------------------
  ! --- 1 scalar ------------------------------------------------------------------------
  function init_1scl(this, ncol, nlay, ngpt, name) result(err_message)
    class(ty_optical_props_1scl)    :: this
    integer,             intent(in) :: ncol, nlay, ngpt
    character(len=*), optional, intent(in) :: name
    character(len=128)              :: err_message

    err_message = "" 
    if(any([ncol, nlay, ngpt] <= 0)) then 
      err_message = "optical_props%init: must provide positive extents for ncol, nlay, ngpt"
    else
      if(allocated(this%tau)) deallocate(this%tau) 
      allocate(this%tau(ncol,nlay,ngpt))
    end if 
    if(present(name)) this%name = name
  end function init_1scl

  ! --- 2 stream ------------------------------------------------------------------------
  function init_2str(this, ncol, nlay, ngpt, name) result(err_message)
    class(ty_optical_props_2str)    :: this
    integer,             intent(in) :: ncol, nlay, ngpt
    character(len=*), optional, intent(in) :: name
    character(len=128)              :: err_message

    err_message = "" 
    if(any([ncol, nlay, ngpt] <= 0)) then 
      err_message = "optical_props%init: must provide positive extents for ncol, nlay, ngpt"
    else
      if(allocated(this%tau)) deallocate(this%tau) 
      allocate(this%tau(ncol,nlay,ngpt))
    end if 
    if(present(name)) this%name = name
    if(allocated(this%ssa)) deallocate(this%ssa) 
    if(allocated(this%g  )) deallocate(this%g  ) 
    if(err_message == "") allocate(this%ssa(ncol,nlay,ngpt), this%g(ncol,nlay,ngpt))
  end function init_2str

  ! --- n stream ------------------------------------------------------------------------
  function init_nstr(this, nmom, ncol, nlay, ngpt, name) result(err_message)
    class(ty_optical_props_nstr)    :: this
    integer, intent(in)             :: nmom ! number of moments
    integer,             intent(in) :: ncol, nlay, ngpt
    character(len=*), optional, intent(in) :: name
    character(len=128)              :: err_message

    err_message = "" 
    if(any([ncol, nlay, ngpt] <= 0)) then 
      err_message = "optical_props%init: must provide positive extents for ncol, nlay, ngpt"
    else
      if(allocated(this%tau)) deallocate(this%tau) 
      allocate(this%tau(ncol,nlay,ngpt))
    end if 
    if(present(name)) this%name = name
    if(allocated(this%ssa)) deallocate(this%ssa) 
    if(allocated(this%p  )) deallocate(this%p  ) 
    if(err_message == "") allocate(this%ssa(ncol,nlay,ngpt), this%p(nmom,ncol,nlay,ngpt))
  end function init_nstr
  ! ------------------------------------------------------------------------------------------
  ! --- delta scaling ------------------------------------------------------------------------
  ! ------------------------------------------------------------------------------------------
  function delta_scale_1scl(this, for) result(err_message)
    class(ty_optical_props_1scl), intent(inout) :: this
    real(wp), dimension(:,:,:), target, optional, &
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
    real(wp), dimension(:,:,:), target, optional, &
                                  intent(in   ) :: for
    ! Forward scattering fraction; g**2 if not provided
    character(128)                              :: err_message

    integer :: ncol, nlay, ngpt
    integer :: icol, ilay, igpt
    real(wp), dimension(size(this%tau,dim=1), size(this%tau,dim=2), size(this%tau,dim=3)), target :: f_loc
    real(wp), dimension(:,:,:), pointer :: f
    real(wp) :: wf ! Temporary -- should it be vector?
    ! --------------------------------
    ncol = size(this%tau,dim=1)
    nlay = size(this%tau,dim=2)
    ngpt = size(this%tau,dim=3)
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
      f => for
    else
      f_loc(:,:,:) = this%g(:,:,:) * this%g(:,:,:)
      f => f_loc
    end if

    do igpt = 1, ngpt
      do ilay = 1, nlay
        do icol = 1, ncol
          wf = this%ssa(icol,ilay,igpt) * f(icol,ilay,igpt)
          this%tau(icol,ilay,igpt) = (1._wp - wf) * this%tau(icol,ilay,igpt)
          this%ssa(icol,ilay,igpt) = (this%ssa(icol,ilay,igpt) - wf) / (1.0_wp - wf)
          this%g  (icol,ilay,igpt) = (this%g  (icol,ilay,igpt) - f(icol,ilay,igpt)) / &
                                       (1._wp - f(icol,ilay,igpt))
        end do
      end do
    end do

  end function delta_scale_2str
  ! ------------------------------------------------------------------------------------------

  function delta_scale_nstr(this, for) result(err_message)
    class(ty_optical_props_nstr), intent(inout) :: this
    real(wp), dimension(:,:,:), target, optional, &
                                 intent(in   ) :: for
    ! Forward scattering fraction; g**2 if not provided
    character(128)                             :: err_message

    err_message = 'delta_scale_nstr: Not yet implemented'
  end function delta_scale_nstr
  ! ------------------------------------------------------------------------------------------
  ! --- Validate values ----------------------------------------------------------------------
  ! ------------------------------------------------------------------------------------------
  function validate_1scalar(this) result(err_message)
    class(ty_optical_props_1scl), intent(in) :: this
    character(len=128)                       :: err_message
    
    err_message = ''
    if(any(this%tau <  0._wp)) & 
      err_message = "validate: tau values out of range" 
    if(len_trim(err_message) > 0 .and. len_trim(this%name) > 0) &
      err_message = trim(this%name) // ': ' // trim(err_message) 
    
  end function validate_1scalar
  
  ! ------------------------------------------------------------------------------------------
  function validate_2stream(this) result(err_message)
    class(ty_optical_props_2str), intent(in) :: this
    character(len=128)                       :: err_message
    
    err_message = ''
    if(any(this%tau <  0._wp)) & 
      err_message = "validate: tau values out of range" 
    if(any(this%ssa <  0._wp) .or. any(this%ssa > 1._wp)) & 
      err_message = "validate: ssa values out of range" 
    if(any(this%g   < -1._wp) .or. any(this%g   > 1._wp)) & 
      err_message = "validate: g values out of range" 
    if(len_trim(err_message) > 0 .and. len_trim(this%name) > 0) &
      err_message = trim(this%name) // ': ' // trim(err_message) 
    
  end function validate_2stream

  ! ------------------------------------------------------------------------------------------
  function validate_nstream(this) result(err_message)
    class(ty_optical_props_nstr), intent(in) :: this
    character(len=128)                       :: err_message
    
    err_message = ''
    if(any(this%tau <  0._wp)) & 
      err_message = "validate: tau values out of range" 
    if(any(this%ssa <  0._wp) .or. any(this%ssa > 1._wp)) & 
      err_message = "validate: ssa values out of range" 
    if(any(this%p(2,:,:,:) < -1._wp) .or. & 
       any(this%p(2,:,:,:) > 1._wp))      & 
      err_message = "validate: p(2,:,:,:)  = g values out of range" 
    if(len_trim(err_message) > 0 .and. len_trim(this%name) > 0) &
      err_message = trim(this%name) // ': ' // trim(err_message) 
  end function validate_nstream

  ! ------------------------------------------------------------------------------------------
  ! --- Return subsets of optical properties arrays along x (col) direction
  ! ------------------------------------------------------------------------------------------
  
  !
  ! Allocate class, then arrays; copy. Could probably be more efficient if 
  !   classes used pointers internally. 
  !
  
  ! ------------------------------------------------------------------------------------------
  ! This set takes start position and number as scalars
  ! ------------------------------------------------------------------------------------------
  
  function subset_1scl_range(full, start, n, subset) result(err_message)
    class(ty_optical_props_1scl), intent(inout) :: full
    integer,                      intent(in   ) :: start, n
    class(ty_optical_props_arry), intent(inout) :: subset 
    character(128)                              :: err_message

    integer :: nlay, ngpt, nmom

    err_message = "" 
    nlay = size(full%tau,dim=2)
    ngpt = size(full%tau,dim=3)
    if(start < 1 .or. start + n-1 > size(full%tau, 1)) & 
       err_message = "optical_props%subset: Asking for columns outside range" 
    if(err_message /= "") return 
    
    ! Seems like the deallocation statements should be needed under Fortran 2003    
    !   but Intel compiler doesn't run without them
    call subset%set_name(full%get_name())
    if(allocated(subset%tau)) deallocate(subset%tau)
    select type (subset)
      class is (ty_optical_props_1scl)  
        err_message = subset%init_1scl(n, nlay, ngpt)
        if(err_message /= "") return 
      class is (ty_optical_props_2str)     
        if(allocated(subset%ssa)) deallocate(subset%ssa)
        if(allocated(subset%g  )) deallocate(subset%g  )
        err_message = subset%init_2str(n, nlay, ngpt)
        if(err_message /= "") return 
        subset%ssa(1:n,:,:) = 0._wp
        subset%g  (1:n,:,:) = 0._wp
      class is (ty_optical_props_nstr) 
        if(allocated(subset%ssa)) deallocate(subset%ssa)
        if(allocated(subset%p  )) then
          nmom = size(subset%p,1)
          deallocate(subset%p  )
        else
          nmom = 1
        end if
        err_message = subset%init_nstr(nmom, n, nlay, ngpt)
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
    nlay = size(full%tau,dim=2)
    ngpt = size(full%tau,dim=3)
    if(start < 1 .or. start + n-1 > size(full%tau, 1)) & 
       err_message = "optical_props%subset: Asking for columns outside range" 
    if(err_message /= "") return 
    
    call subset%set_name(full%get_name())
    if(allocated(subset%tau)) deallocate(subset%tau)
    select type (subset)
      class is (ty_optical_props_1scl)  
        err_message = subset%init_1scl(n, nlay, ngpt)
        if(err_message /= "") return 
        subset%tau(1:n,:,:) = full%tau(start:start+n-1,:,:) * & 
                     (1._wp - full%ssa(start:start+n-1,:,:))
      class is (ty_optical_props_2str)     
        if(allocated(subset%ssa)) deallocate(subset%ssa)
        if(allocated(subset%g  )) deallocate(subset%g  )
        err_message = subset%init_2str(n, nlay, ngpt)
        if(err_message /= "") return 
        subset%tau(1:n,:,:) = full%tau(start:start+n-1,:,:) 
        subset%ssa(1:n,:,:) = full%ssa(start:start+n-1,:,:)
        subset%g  (1:n,:,:) = full%g  (start:start+n-1,:,:)
      class is (ty_optical_props_nstr) 
        if(allocated(subset%ssa)) deallocate(subset%ssa)
        if(allocated(subset%p  )) then
          nmom = size(subset%p,1)
          deallocate(subset%p  )
        else
          nmom = 1
        end if
        err_message = subset%init_nstr(nmom, n, nlay, ngpt)
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
    nlay = size(full%tau,dim=2)
    ngpt = size(full%tau,dim=3)
    nmom = size(full%p  ,dim=1) 
    if(start < 1 .or. start + n-1 > size(full%tau, 1)) & 
       err_message = "optical_props%subset: Asking for columns outside range" 
    if(err_message /= "") return 
    
    call subset%set_name(full%get_name())
    if(allocated(subset%tau)) deallocate(subset%tau)
    select type (subset)
      class is (ty_optical_props_1scl)  
        err_message = subset%init_1scl(n, nlay, ngpt)
        if(err_message /= "") return 
        subset%tau(1:n,:,:) = full%tau(start:start+n-1,:,:) * & 
                     (1._wp - full%ssa(start:start+n-1,:,:))
      class is (ty_optical_props_2str)     
        if(allocated(subset%ssa)) deallocate(subset%ssa)
        if(allocated(subset%g  )) deallocate(subset%g  )
        err_message = subset%init_2str(n, nlay, ngpt)
        if(err_message /= "") return 
        subset%tau(1:n,:,:) = full%tau(start:start+n-1,:,:) 
        subset%ssa(1:n,:,:) = full%ssa(start:start+n-1,:,:)
        subset%g  (1:n,:,:) = full%p(1,start:start+n-1,:,:)
      class is (ty_optical_props_nstr) 
        if(allocated(subset%ssa)) deallocate(subset%ssa)
        if(allocated(subset%p  )) deallocate(subset%p  )
        err_message = subset%init_nstr(nmom, n, nlay, ngpt)
        if(err_message /= "") return 
        subset%ssa(1:n,:,:) = full%ssa(start:start+n-1,:,:)
        subset%p(:,1:n,:,:) = full%p(:,start:start+n-1,:,:)
    end select 

  end function subset_nstr_range
  
  ! -----------------------------------------------------------------------------------------
  ! --- increment: add optical properties of op1 to op2 
  !   The increment function, called as op1%increment(op2) or op1%increment(op2, band_limits), 
  !   adds the values of op1 to op2, changing op2 and leaving op1 untouched. This lets us make 
  !   functions that implement the same behavior. 
  !  
  !  It might make sense to break this big function up in type-bound procedures for each 
  !    variant of ty_optical_props_arry 
  ! -----------------------------------------------------------------------------------------
  function increment_arry(op1, op2)result(err_message)
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
           ncol = size(op2%tau,1) 
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
          call increment_2stream_by_nstream(ncol, nlay, ngpt, size(op1%p, 1), &
                                            op2%tau, op2%ssa, op2%g, &
                                            op1%tau, op1%ssa, op1%p)
      end select
      
    class is (ty_optical_props_nstr)
      nmom1 = size(op2%p, 1) 
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
          call increment_nstream_by_nstream(ncol, nlay, ngpt, nmom1, size(op1%p, 1), &
                                            op2%tau, op2%ssa, op2%p, &
                                            op1%tau, op1%ssa, op1%p)
      end select
    end select
  end function increment_arry
  ! -----------------------------------------------------------------------------------------
  function increment_arry_byband(op1, op2, gpt_lims) result(err_message)
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
    if(any([op1%get_ncol(), op1%get_nlay()] /= [ncol, nlay])) then
      err_message = "ty_optical_props%increment_by: " // &
                    "optical properties objects are inconsistently sized"
      return
    end if
    if(op1%get_ngpt() /= nbnd) then
      err_message = "ty_optical_props%increment_by: " // & 
                    "number of bands not consistent between " // &
                    "optical properties objects, g-point limits" 
      return
    end if
    if(minval(gpt_lims) < 1 .or. maxval(gpt_lims) > ngpt) then 
      err_message = "ty_optical_props%increment_by: " // & 
                    "band limits not consistent with number of gpoints" 
      return
    end if
    if(err_message /= "") return 

    ! Validate
    err_message = op1%validate()
    if (err_message /= "") then
       err_message = err_message // " (before incrementing)"
       return
    end if
    err_message = op2%validate()
    if (err_message /= "") then
       err_message = err_message // " (before incrementing)"
       return
    end if
        
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
            call inc_2stream_by_nstream_bybnd(ncol, nlay, ngpt, size(op1%p, 1), &
                                              op2%tau, op2%ssa, op2%g, &
                                              op1%tau, op1%ssa, op1%p, & 
                                              nbnd, gpt_lims)
        end select
      
      class is (ty_optical_props_nstr)
        nmom1 = size(op2%p, 1) 
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
            call inc_nstream_by_nstream_bybnd(ncol, nlay, ngpt, nmom1, size(op1%p, 1), &
                                              op2%tau, op2%ssa, op2%p, &
                                              op1%tau, op1%ssa, op1%p, & 
                                              nbnd, gpt_lims)
        end select
    end select

    ! Validate
    err_message = op1%validate()
    if (err_message /= "") then
       err_message = err_message // " (after incrementing)"
       return
    end if
    err_message = op2%validate()
    if (err_message /= "") then
       err_message = err_message // " (after incrementing)"
       return
    end if
 
  end function increment_arry_byband
  ! ------------------------------------------------------------------------------------------
  ! Incrementing -- adding optical properties together 
  !   Originally we used increment_by, so that op1%increment_by(op1) changed the values of op1 
  !   and left the values of op2 untouched. These are the increment_by routines. increment_by_band 
  !   also expands values defined on bands to values defined on g-points (really, the extent of 
  !   the third array dimension). 
  ! ------------------------------------------------------------------------------------------
  function increment_gpt_by(op1, op2)result(err_message)
    class(ty_optical_props_arry), intent(inout) :: op1
    class(ty_optical_props     ), intent(in   ) :: op2
    character(128)                         :: err_message
    ! ----- 
    integer :: ncol, nlay, ngpt, nmom1
    ! ----- 
    err_message = op2%increment(op1) 
    
  end function increment_gpt_by
  ! -----------------------------------------------------------------------------------------
  function increment_band_by(op1, op2, gpt_lims) result(err_message)
    class(ty_optical_props_arry), intent(inout) :: op1
    class(ty_optical_props     ), intent(in   ) :: op2
    integer, dimension(:,:),      intent(in   ) :: gpt_lims  ! (begin g-point, end g-point) = gpt_lims(2,band)
    character(128)                              :: err_message

    ! ----- 
    integer :: ncol, nlay, ngpt, nbnd, nmom1
    ! -----
    err_message = op2%increment(op1, gpt_lims) 
  end function increment_band_by
  ! ------------------------------------------------------------------------------------------
  ! Functions for problem sizes  
  !   Dimensions is (ncol, nlay, ngpt) although ngpt might also refer to the number of bands
  !   For arrays these are the array extents. 
  ! ------------------------------------------------------------------------------------------
  function get_arry_extent(this, dim) 
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
  function get_ncol_arry(this) 
    class(ty_optical_props_arry), intent(in   ) :: this
    integer                                     :: get_ncol_arry

    get_ncol_arry = get_arry_extent(this, 1) 
  end function get_ncol_arry
  ! ------------------------------------------------------------------------------------------
  function get_nlay_arry(this) 
    class(ty_optical_props_arry), intent(in   ) :: this
    integer                                     :: get_nlay_arry
    
    get_nlay_arry = get_arry_extent(this, 2) 
  end function get_nlay_arry
  ! ------------------------------------------------------------------------------------------
  function get_ngpt_arry(this) 
    class(ty_optical_props_arry), intent(in   ) :: this
    integer                                     :: get_ngpt_arry
    
    get_ngpt_arry = get_arry_extent(this, 3) 
  end function get_ngpt_arry
  ! ------------------------------------------------------------------------------------------
  ! ------------------------------------------------------------------------------------------
  !
  ! Functions for the base class 
  !
  ! ------------------------------------------------------------------------------------------
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
    character(len=name_len)                      :: get_name
      
      get_name = trim(this%name)
  end function get_name
  ! ------------------------------------------------------------------------------------------
end module mo_optical_props
