module radiation
!---------------------------------------------------------------------------------
! Purpose:
!
! CAM interface to RRTMGP
!
! Revision history:
! May 2004, D. B. Coleman, Initial version of interface module.
! Jul 2004, B. Eaton,      Use interfaces from new shortwave, longwave, and ozone modules.
! Feb 2005, B. Eaton,      Add namelist variables and control of when calcs are done.
! May 2008, Mike Iacono    Initial version for RRTMG
! Jun 2009, Minghuai Wang, The MMF treatment is added to the subroutine of 
!                          radiation_tend. These modifications are based on the 
!                          spcam3.5, which was developled by Marat Khairoutdinov. 
!                          The spcam3.5 only have one radiation package (camrt). 
!                          See comments in radiation_tend for the details. 
! Jul 2009, Minghuai Wang, For the Morrison's two momenent microphysics in SAM, 
!                          droplet and ice crystal effective radius used in the 
!                          radiation code are calculated at each CRM column by 
!                          calling m2005_effradius
! Oct 2009, Minghuai Wang, CRM-scale aerosol water is used to calculate aerosol 
!                          optical depth
! Nov 2010, J. Kay,        Add COSP simulator calls
! Nov 2017, B. Hillman,    Modify extensively for use with RRTMGP
!---------------------------------------------------------------------------------

! E3SM-specific modules that are used throughout this module. An effort was made
! to keep imports as local as possible, so we only load a few of these at the
! module (rather than the subroutine) level.
use shr_kind_mod,     only: r8=>shr_kind_r8, cl=>shr_kind_cl
use ppgrid,           only: pcols, pver, pverp, begchunk, endchunk
use cam_abortutils,   only: endrun
use scamMod,          only: scm_crm_mode, single_column, swrad_off
use rad_constituents, only: N_DIAG

! RRTMGP gas optics object to store coefficient information. This is imported
! here so that we can make the k_dist objects module data and only load them
! once.
use mo_gas_optics, only: ty_gas_optics_specification

! PIO libraries needed for reading netcdf files
use cam_pio_utils, only: cam_pio_openfile
use pio, only: file_desc_t, var_desc_t,               &
               PIO_NOERR, PIO_INTERNAL_ERROR,         &
               pio_seterrorhandling, PIO_BCAST_ERROR, &
               pio_inq_dimlen, pio_inq_dimid,         &
               pio_inq_varid, pio_def_var,            &
               pio_put_var, pio_get_var,              &
               PIO_NOWRITE, pio_closefile
use netcdf

implicit none
private
save

! Public routines provided by this module.
! TODO: radiation_defaultopts, radiation_setops, and radiation_printops exist
! only because the radiation namelist has traditionally been read at the driver
! level by runtime_opts. I am not sure this is the best solution, and in fact
! CESM seems to have gone away from this practice altogether, opting instead for
! individual modules to be responsible for reading their own namelists. This
! might be a better practice going forward, and would simplify the logic here
! somewhat. To do this, we would need to use the radiation_readnl routine
! (copied below), and add a line in runtime_opts that calls this routine, and
! then remove the code in runtime_opts that reads the radiation namelist
! variables from the CAM namelist. The calls to radiation_defaultopts,
! radiation_setops, and radiation_printops would also then need to be removed
! from runtime_ops.
public :: &
   radiation_register,    &! registers radiation physics buffer fields
   radiation_defaultopts, &! set default values of namelist variables in runtime_opts
   radiation_setopts,     &! set namelist values from runtime_opts
   radiation_printopts,   &! print namelist values to log
   radiation_nextsw_cday, &! calendar day of next radiation calculation
   radiation_do,          &! query which radiation calcs are done this timestep
   radiation_init,        &! calls radini
   radiation_readnl,      &! read radiation namelist
   radiation_tend          ! moved from radctl.F90

! Counter variables for use with the CFMIP Observation Simulator Package (COSP).
! TODO: This seems like somewhat of an awkward way of determining when to run
! COSP, and it might make more sense to implement a more elegant routine to call
! that determines whether or not to run COSP at a given timestep (similar to the
! radiation_do routine in this module).
integer, public, allocatable :: cosp_cnt(:)       ! counter for cosp
integer, public              :: cosp_cnt_init = 0 ! initial value for cosp counter

! Private module data
integer :: qrs_idx      = 0 
integer :: qrl_idx      = 0 
integer :: su_idx       = 0 
integer :: sd_idx       = 0 
integer :: lu_idx       = 0 
integer :: ld_idx       = 0 
integer :: cldfsnow_idx = 0 
integer :: cld_idx      = 0 
integer :: concld_idx   = 0
integer :: rel_idx      = 0
integer :: rel_fn_idx   = 0 
integer :: rei_idx      = 0
integer :: dei_idx      = 0

! Declare namelist variables as module data. This also sets default values for
! namelist variables.
integer :: iradsw = -1  ! freq. of shortwave radiation calc in time steps (positive)
                        ! or hours (negative).
integer :: iradlw = -1  ! frequency of longwave rad. calc. in time steps (positive)
                        ! or hours (negative).
integer :: irad_always = 0  ! Specifies length of time in timesteps (positive)
                            ! or hours (negative) SW/LW radiation will be
                            ! run continuously from the start of an
                            ! initial or restart run

! The spectralflux flag determines whether or not spectral (per band) fluxes are
! calculated. If true, upward and downward fluxes are calculated per band,
! otherwise just broadband "shortwave" and "longwave" fluxes are calculated.
! TODO: while it seems that setting spectralflux = .true. add spectral fluxes to
! the physics buffer, I do not see where these fields are added to the history
! buffer. It might be that the output of these fields are just handled by the
! physics buffer output routines (they are declared with "global" scope, so
! written on restart files at least), but it would be good to include them on
! the history buffer too so that we can annotate them with a description of what
! they are, rather than expecting the user to know what they are from the pbuf
! fields.
logical :: spectralflux  = .false.  ! calculate fluxes (up and down) per band.

! Flag to indicate whether or not to use the radiation timestep for solar zenith
! angle calculations. If true, use the radiation timestep for all solar zenith 
! angle (cosz) calculations.
! TODO: How does this differ if value is .false.?
logical :: use_rad_dt_cosz  = .false. 

! Model data that is not controlled by namelist fields specifically follows
! below.

character(len=4) :: diag(0:N_DIAG) =(/'    ','_d1 ','_d2 ','_d3 ','_d4 ','_d5 ', &
                                      '_d6 ','_d7 ','_d8 ','_d9 ','_d10'/)

! Output diagnostic brightness temperatures at the top of the
! atmosphere for 7 TOVS/HIRS channels (2,4,6,8,10,11,12) and 4 TOVS/MSU 
! channels (1,2,3,4).
! TODO: where are these options set?
logical :: dohirs = .false. 
integer :: ihirsfq = 1      ! frequency (timesteps) of brightness temperature calcs

! time step to use for the shr_orb_cosz calculation, if use_rad_dt_cosz set to true
! TODO: where is this set, and what is shr_orb_cosz? Alternative solar zenith
! angle calculation? What is the other behavior?
real(r8) :: dt_avg = 0.0_r8  

! k-distribution coefficients. These will be populated by reading from the
! RRTMGP coefficients files, specified by coefficients_file_sw and
! coefficients_file_lw in the radiation namelist. They exist as module data
! because we only want to load those files once.
type(ty_gas_optics_specification) :: k_dist_sw, k_dist_lw

! k-distribution coefficients files to read from. These are set via namelist
! variables.
character(len=cl) :: coefficients_file_sw, coefficients_file_lw

! Number of shortwave and longwave bands in use by the RRTMGP radiation code.
! This information will be stored in the k_dist_sw and k_dist_lw objects and may
! be retrieved using the k_dist_sw%get_nbands() and k_dist_lw%get_nbands()
! methods, but I think we need to save these as private module data so that we
! can automatically allocate arrays later in subroutine headers, i.e.:
!
!     real(r8) :: cld_tau(pcols,pver,nswbands)
!
! and so forth. Previously some of this existed in radconstants.F90, but I do
! not think we need to use that.
integer :: nswbands, nlwbands

! These are used in the set_gases routines to convert between mmr and vmr
real(r8), parameter :: amdw = 1.607793_r8    ! Molecular weight of dry air / water vapor
real(r8), parameter :: amdc = 0.658114_r8    ! Molecular weight of dry air / carbon dioxide
real(r8), parameter :: amdo = 0.603428_r8    ! Molecular weight of dry air / ozone
real(r8), parameter :: amdm = 1.805423_r8    ! Molecular weight of dry air / methane
real(r8), parameter :: amdn = 0.658090_r8    ! Molecular weight of dry air / nitrous oxide
real(r8), parameter :: amdo2 = 0.905140_r8   ! Molecular weight of dry air / oxygen
real(r8), parameter :: amdc1 = 0.210852_r8   ! Molecular weight of dry air / CFC11
real(r8), parameter :: amdc2 = 0.239546_r8   ! Molecular weight of dry air / CFC12

! Indices for copying data between cam and rrtmgp arrays
! Assume the rrtmgp vertical index goes bottom to top of atm
integer :: ktopcamm ! cam index of top layer
integer :: ktopradm ! rrtmgp index of layer corresponding to ktopcamm
integer :: ktopcami ! cam index of top interface
integer :: ktopradi ! rrtmgp index of interface corresponding to ktopcami

type, private :: gas_optics_coefficients
   character(32), dimension(:), allocatable :: gas_names
   character(256), dimension(:), allocatable :: gas_minor
   character(256), dimension(:), allocatable :: identifier_minor
   character(256), dimension(:), allocatable :: minor_gases_lower
   character(256), dimension(:), allocatable :: minor_gases_upper
   character(256), dimension(:), allocatable :: scaling_gas_lower
   character(256), dimension(:), allocatable :: scaling_gas_upper
   integer,  dimension(:,:,:), allocatable :: key_species
   integer,  dimension(:,:), allocatable :: minor_limits_gpt_lower
   integer,  dimension(:,:), allocatable :: minor_limits_gpt_upper
   logical,  dimension(:), allocatable :: minor_scales_with_density_lower
   logical,  dimension(:), allocatable :: minor_scales_with_density_upper
   logical,  dimension(:), allocatable :: scale_by_complement_lower
   logical,  dimension(:), allocatable :: scale_by_complement_upper
   integer,  dimension(:), allocatable :: kminor_start_lower
   integer,  dimension(:), allocatable :: kminor_start_upper
   integer,  dimension(:,:), allocatable :: bnd_limits_gpt
   real(r8), dimension(:,:), allocatable :: bnd_limits_wavenumber
   real(r8), dimension(:), allocatable :: press_ref, temp_ref
   real(r8) :: press_ref_trop
   real(r8) :: absorption_coefficient_ref_T
   real(r8) :: absorption_coefficient_ref_P
   real(r8), dimension(:,:,:), allocatable :: vmr_ref
   real(r8), dimension(:,:,:,:), allocatable :: kmajor
   real(r8), dimension(:,:,:), allocatable :: kminor_lower, kminor_upper


   real(r8), dimension(:,:), allocatable :: totplnk
   real(r8), dimension(:,:,:,:), allocatable :: planck_fraction
   real(r8), dimension(:), allocatable :: solar_src
   real(r8), dimension(:,:,:), allocatable :: rayl_lower, rayl_upper

contains
   procedure, public :: load_from_file=>coefficients_load_from_file
   procedure, public :: free_allocated_memory=>coefficients_free_allocated_memory

end type gas_optics_coefficients

type(gas_optics_coefficients) :: coefficients_sw, coefficients_lw

! give this module a name
character(len=*), parameter :: module_name = 'radiation'

interface read_field
   module procedure read_real_scalar, &
                    read_real_1d, read_real_2d, read_real_3d, read_real_4d, &
                    read_integer_1d, read_integer_2d, read_integer_3d, read_integer_4d, &
                    read_logical_1d, read_character_1d
end interface read_field

!===============================================================================

contains

function coefficients_free_allocated_memory(this) result(error_message)
   class(gas_optics_coefficients), intent(inout) :: this
   character(len=128) :: error_message

   ! Naively deallocate; TODO: better to check if allocated like below
   deallocate(this%gas_names, this%key_species, this%bnd_limits_gpt, &
              this%bnd_limits_wavenumber, this%press_ref, &
              this%temp_ref, &
              this%vmr_ref, this%kmajor, &
              this%kminor_lower, this%kminor_upper)

   ! These may or may not be allocated
   if (allocated(this%totplnk)) deallocate(this%totplnk)
   if (allocated(this%planck_fraction)) deallocate(this%planck_fraction)
   if (allocated(this%solar_src)) deallocate(this%solar_src)
   if (allocated(this%rayl_lower)) deallocate(this%rayl_lower)
   if (allocated(this%rayl_upper)) deallocate(this%rayl_upper)

   ! Return empty string for no errors
   error_message = ""

end function coefficients_free_allocated_memory

!===============================================================================

subroutine radiation_readnl(nlfile, dtime_in)
!-------------------------------------------------------------------------------
! Purpose: Read radiation_nl namelist group.
!-------------------------------------------------------------------------------

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
#ifdef SPMD
   !use mpishorthand,    only: mpicom, mpilog, mpiint, mpichar
#endif
   use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_integer, mpi_logical, &
                              mpi_character, masterproc
   use time_manager, only: get_step_size
   use cam_logfile, only: iulog

   ! File containing namelist input
   character(len=*), intent(in) :: nlfile
   integer, intent(in), optional :: dtime_in

   ! Local variables
   integer :: unitn, ierr
   integer :: dtime  ! timestep size
   character(len=*), parameter :: subroutine_name = 'radiation_readnl'

   character(len=cl) :: rrtmgp_coefficients_file_lw, rrtmgp_coefficients_file_sw
   integer :: rrtmgp_iradsw, rrtmgp_iradlw, rrtmgp_irad_always
   logical :: rrtmgp_use_rad_dt_cosz, rrtmgp_spectralflux

   ! Variables defined in namelist
   ! TODO: why are these prefaced with rrtmgp instead of radiation?
   namelist /radiation_nl/ rrtmgp_coefficients_file_lw, rrtmgp_coefficients_file_sw, &
                        rrtmgp_iradsw, rrtmgp_iradlw, rrtmgp_irad_always, &
                        rrtmgp_use_rad_dt_cosz, rrtmgp_spectralflux

   ! Read the namelist, only if called from master process
   ! TODO: better documentation and cleaner logic here?
   if (masterproc) then
      unitn = getunit()
      open(unitn, file=trim(nlfile), status='old')
      call find_group_name(unitn, 'radiation_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, radiation_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subroutine_name // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

   ! Broadcast namelist variables
   call mpibcast(rrtmgp_coefficients_file_lw, cl, mpi_character, mstrid, mpicom, ierr)
   call mpibcast(rrtmgp_coefficients_file_sw, cl, mpi_character, mstrid, mpicom, ierr)
   call mpibcast(rrtmgp_iradsw, 1, mpi_integer, mstrid, mpicom, ierr)
   call mpibcast(rrtmgp_iradlw, 1, mpi_integer, mstrid, mpicom, ierr)
   call mpibcast(rrtmgp_irad_always, 1, mpi_integer, mstrid, mpicom, ierr)
   call mpibcast(rrtmgp_use_rad_dt_cosz, 1, mpi_logical, mstrid, mpicom, ierr)
   call mpibcast(rrtmgp_spectralflux, 1, mpi_logical, mstrid, mpicom, ierr)

   ! Set module data
   coefficients_file_lw = rrtmgp_coefficients_file_lw
   coefficients_file_sw = rrtmgp_coefficients_file_sw
   iradsw          = rrtmgp_iradsw
   iradlw          = rrtmgp_iradlw
   irad_always     = rrtmgp_irad_always
   use_rad_dt_cosz = rrtmgp_use_rad_dt_cosz
   spectralflux    = rrtmgp_spectralflux

   ! Convert iradsw, iradlw and irad_always from hours to timesteps if necessary
   if (present(dtime_in)) then
      dtime = dtime_in
   else
      dtime  = get_step_size()
   end if
   if (iradsw      < 0) iradsw      = nint((-iradsw     *3600._r8)/dtime)
   if (iradlw      < 0) iradlw      = nint((-iradlw     *3600._r8)/dtime)
   if (irad_always < 0) irad_always = nint((-irad_always*3600._r8)/dtime)

   !-----------------------------------------------------------------------
   ! Print runtime options to log.
   !-----------------------------------------------------------------------

   if (masterproc) then
      write(iulog,*) 'RRTMGP radiation scheme parameters:'
      write(iulog,10) trim(coefficients_file_lw), trim(coefficients_file_sw), iradsw, iradlw, &
                      irad_always, use_rad_dt_cosz, spectralflux
   end if

10 format('  LW coefficents file: ',                                a/, &
          '  SW coefficents file: ',                                a/, &
          '  Frequency (timesteps) of Shortwave Radiation calc:  ',i5/, &
          '  Frequency (timesteps) of Longwave Radiation calc:   ',i5/, &
          '  SW/LW calc done every timestep for first N steps. N=',i5/, &
          '  Use average zenith angle:                           ',l5/, &
          '  Output spectrally resolved fluxes:                  ',l5/)

end subroutine radiation_readnl

subroutine radiation_register()

   !----------------------------------------------------------------------------
   ! Register radiation fields in the physics buffer
   !----------------------------------------------------------------------------

   use physics_buffer, only: pbuf_add_field, dtype_r8
   !use radconstants,   only: nswbands, nlwbands

   integer :: idx  ! dummy index for adding fields to physics buffer

   ! Heating rate profiles; QRS is the shortwave radiative heating rate, and QRL
   ! is the longwave radiative heating rate
   ! TODO: Do QRS and QRL need to be set to "global" scope on the physics
   ! buffer? Doing so forces them to be written to restart files, but is this
   ! needed? It does not look like their values get read anywhere in this
   ! module, just overwritten, so they may not need to be written to restarts.
   call pbuf_add_field('QRS', 'global', dtype_r8, (/pcols,pver/), qrs_idx)
   call pbuf_add_field('QRL', 'global', dtype_r8, (/pcols,pver/), qrl_idx)

   ! Chemistry interface needs shortwave down at surface
   ! TODO: why add the others? Are these needed?
   call pbuf_add_field('FSDS', 'global', dtype_r8, (/pcols/), idx)
   call pbuf_add_field('FSNS', 'global', dtype_r8, (/pcols/), idx)
   call pbuf_add_field('FSNT', 'global', dtype_r8, (/pcols/), idx)
   call pbuf_add_field('FLNS', 'global', dtype_r8, (/pcols/), idx)
   call pbuf_add_field('FLNT', 'global', dtype_r8, (/pcols/), idx)

   ! If the namelist has been configured for preserving the spectral fluxes, then create
   ! physics buffer variables to store the results. These are fluxes per
   ! spectral band, as follows:
   !     SU: shortwave up
   !     SD: shortwave down
   !     LU: longwave up
   !     LD: longwave down
   ! TODO: Do these need to be "global"? Why not "physpkg"?
   if (spectralflux) then
      call pbuf_add_field('SU', 'global', dtype_r8, (/pcols,pverp,nswbands/), su_idx)
      call pbuf_add_field('SD', 'global', dtype_r8, (/pcols,pverp,nswbands/), sd_idx)
      call pbuf_add_field('LU', 'global', dtype_r8, (/pcols,pverp,nlwbands/), lu_idx)
      call pbuf_add_field('LD', 'global', dtype_r8, (/pcols,pverp,nlwbands/), ld_idx)
   end if

end subroutine radiation_register

!===============================================================================

subroutine radiation_defaultopts(iradsw_out, iradlw_out, iradae_out, &
                                 irad_always_out, spectralflux_out, &
                                 use_rad_dt_cosz_out, &
                                 coefficients_file_sw_out, &
                                 coefficients_file_lw_out)

   !---------------------------------------------------------------------------- 
   ! Purpose: Return default runtime options
   ! TODO: remove this routine
   !----------------------------------------------------------------------------

   integer, intent(out), optional :: iradsw_out
   integer, intent(out), optional :: iradlw_out
   integer, intent(out), optional :: iradae_out
   integer, intent(out), optional :: irad_always_out
   logical, intent(out), optional :: spectralflux_out
   logical, intent(out), optional :: use_rad_dt_cosz_out
   character(len=*), intent(out), optional :: coefficients_file_sw_out
   character(len=*), intent(out), optional :: coefficients_file_lw_out

   if (present(iradsw_out)) iradsw_out = iradsw
   if (present(iradlw_out)) iradlw_out = iradlw
   if (present(iradae_out)) iradae_out = -999
   if (present(irad_always_out)) irad_always_out = irad_always
   if (present(spectralflux_out)) spectralflux_out = spectralflux
   if (present(use_rad_dt_cosz_out)) use_rad_dt_cosz_out = use_rad_dt_cosz
   if (present(coefficients_file_sw_out)) coefficients_file_sw_out = coefficients_file_sw
   if (present(coefficients_file_lw_out)) coefficients_file_lw_out = coefficients_file_lw

end subroutine radiation_defaultopts

!===============================================================================

subroutine radiation_setopts(dtime, nhtfrq, iradsw_in, iradlw_in, iradae_in, &
                             irad_always_in, spectralflux_in, &
                             use_rad_dt_cosz_in, &
                             coefficients_file_sw_in, coefficients_file_lw_in)
                          
   !----------------------------------------------------------------------------
   ! Purpose: Set runtime options
   !
   ! History: Copied from RRTMG implementation
   !
   ! *** NOTE *** This routine needs information about dtime (init by dycore) 
   !              and nhtfrq (init by history) to do its checking.  Being called
   !              from runtime_opts provides these values possibly before they
   !              have been set in the modules responsible for them.
   !
   ! TODO: remove this routine; not needed if module reads its own namelist.
   !----------------------------------------------------------------------------

   ! Arguments
   integer, intent(in)           :: dtime           ! timestep size (s)
   integer, intent(in)           :: nhtfrq          ! output frequency of primary history file
   integer, intent(in), optional :: iradsw_in
   integer, intent(in), optional :: iradlw_in
   integer, intent(in), optional :: iradae_in
   integer, intent(in), optional :: irad_always_in
   logical, intent(in), optional :: spectralflux_in
   logical, intent(in), optional :: use_rad_dt_cosz_in
   character(len=*), intent(in), optional :: coefficients_file_sw_in
   character(len=*), intent(in), optional :: coefficients_file_lw_in
   
   ! Local namespace
   integer :: iradae   ! not used by RRTMGP

   ! Check for optional arguments and set corresponding module data if present
   if (present(iradsw_in)) iradsw = iradsw_in
   if (present(iradlw_in)) iradlw = iradlw_in
   if (present(iradae_in)) iradae = iradae_in
   if (present(irad_always_in)) irad_always = irad_always_in
   if (present(spectralflux_in)) spectralflux = spectralflux_in
   if (present(use_rad_dt_cosz_in)) use_rad_dt_cosz = use_rad_dt_cosz_in
   if (present(coefficients_file_sw_in)) coefficients_file_sw = coefficients_file_sw_in
   if (present(coefficients_file_lw_in)) coefficients_file_lw = coefficients_file_lw_in

   ! Convert iradsw, iradlw and irad_always from hours to timesteps if necessary
   if (iradsw      < 0) iradsw      = nint((-iradsw     *3600._r8)/dtime)
   if (iradlw      < 0) iradlw      = nint((-iradlw     *3600._r8)/dtime)
   if (irad_always < 0) irad_always = nint((-irad_always*3600._r8)/dtime)

   ! Has user specified iradae?
   if (iradae /= -999) then
      call endrun('radiation_setopts: iradae not used by RRTMGP.')
   end if

end subroutine radiation_setopts

!===============================================================================

subroutine radiation_printopts

   !----------------------------------------------------------------------- 
   ! Purpose: Print runtime options to log.
   !
   ! History: Copied from RRTMG implementation.
   !
   ! TODO: remove this routine; not needed if radiation reads its own namelist.
   !-----------------------------------------------------------------------

   use cam_logfile, only: iulog

   if(irad_always /= 0) write(iulog,10) irad_always
   write(iulog,20) iradsw,iradlw
10 format(' Execute SW/LW radiation continuously for the first ',i5, ' timestep(s) of this run')
20 format(' Frequency of Shortwave Radiation calc. (IRADSW)    ',i5/, &
          ' Frequency of Longwave Radiation calc. (IRADLW)     ',i5)

end subroutine radiation_printopts

!================================================================================================

function radiation_do(op, timestep)

   !---------------------------------------------------------------------------- 
   ! Purpose: Returns true if the specified operation is done this timestep.
   !
   ! History: Copied from RRTMG implementation.
   !----------------------------------------------------------------------------

   use time_manager, only: get_nstep
   use phys_control, only: phys_getopts

   ! Arguments
   character(len=*), intent(in)  :: op             ! name of operation
   integer, intent(in), optional :: timestep       ! model timestep
   logical                       :: radiation_do   ! return value

   ! Local variables
   integer :: nstep             ! current timestep number
   logical :: use_SPCAM

   if (present(timestep)) then
      nstep = timestep
   else
      nstep = get_nstep()
   end if

   ! When the physics package is superparameterization, we do radiation at every
   ! timestep (why?). Otherwise we observe the namelist options that set the
   ! radiation timestep.
   !
   ! TODO: Do we need to require radiation to be done at every timestep for
   ! superparameterized runs? Can we experiment with NOT doing radiation every
   ! timestep? Move this to build/configure step to enable this.
   call phys_getopts(use_SPCAM_out=use_SPCAM)
   if (use_SPCAM) then
      radiation_do = .true. 
   else
      select case (op)
         case ('sw') ! do a shortwave heating calc this timestep?
            radiation_do = nstep == 0 .or. iradsw == 1                     &
                          .or. (mod(nstep-1,iradsw) == 0 .and. nstep /= 1) &
                          .or. nstep <= irad_always
         case ('lw') ! do a longwave heating calc this timestep?
            radiation_do = nstep == 0 .or. iradlw == 1                     &
                          .or. (mod(nstep-1,iradlw) == 0 .and. nstep /= 1) &
                          .or. nstep <= irad_always
         case ('aeres') ! write absorptivity/emissivity to restart file this timestep?
            ! for RRTMG there is no abs/ems restart file
            radiation_do = .false.
         case default
            call endrun('radiation_do: unknown operation:'//op)
      end select
   end if
end function radiation_do

!================================================================================================

real(r8) function radiation_nextsw_cday()
     
   !----------------------------------------------------------------------- 
   ! Purpose: Returns calendar day of next sw radiation calculation
   !
   ! History: Copied from RRTMG implementation.
   !-----------------------------------------------------------------------

   use time_manager, only: get_curr_calday, get_nstep, get_step_size

   ! Local variables
   integer :: nstep      ! timestep counter
   logical :: dosw       ! true => do shosrtwave calc   
   integer :: offset     ! offset for calendar day calculation
   integer :: dTime      ! integer timestep size
   real(r8):: calday     ! calendar day of 
   !-----------------------------------------------------------------------

   radiation_nextsw_cday = -1._r8
   dosw   = .false.
   nstep  = get_nstep()
   dtime  = get_step_size()
   offset = 0
   do while (.not. dosw)
      nstep = nstep + 1
      offset = offset + dtime
      if (radiation_do('sw', nstep)) then
         radiation_nextsw_cday = get_curr_calday(offset=offset) 
         dosw = .true.
      end if
   end do
   if(radiation_nextsw_cday == -1._r8) then
      call endrun('error in radiation_nextsw_cday')
   end if
        
end function radiation_nextsw_cday

!================================================================================================

subroutine radiation_init()
!-------------------------------------------------------------------------------
! Purpose: Initialize the radiation parameterization and add fields to the 
! history buffer
! 
! History: Copied from RRTMG implemenation.
!-------------------------------------------------------------------------------
   use physics_buffer,     only: pbuf_get_index
   use cam_history,        only: addfld, horiz_only, add_default
   use constituents,       only: cnst_get_ind
   use physconst,          only: gravit, cpair, epsilo, stebol, &
                                 pstd, mwdry, mwco2, mwo3
   use phys_control,       only: phys_getopts
   use rad_constituents,   only: N_DIAG, rad_cnst_get_call_list, rad_cnst_get_info
   use cospsimulator_intr, only: docosp, cospsimulator_intr_init
   use hirsbt,             only: hirsbt_init
   use hirsbtpar,          only: hirsname, msuname
   use modal_aer_opt,      only: modal_aer_opt_init
   use time_manager,       only: get_nstep, get_step_size, is_first_restart_step
   use radiation_data,     only: init_rad_data

   ! RRTMGP modules
   use mo_rrtmgp_clr_all_sky, only: &
         rrtmgp_sw_init=>rte_sw_init, &
         rrtmgp_lw_init=>rte_lw_init
   use mo_load_coefficients, only: rrtmgp_load_coefficients=>load_and_init
   use mo_gas_concentrations, only: ty_gas_concs

   integer :: icall, nmodes
   logical :: active_calls(0:N_DIAG)
   integer :: nstep                       ! current timestep number
   logical :: history_amwg                ! output the variables used by the AMWG diag package
   logical :: history_vdiag               ! output the variables used by the AMWG variability diag package
   logical :: history_budget              ! output tendencies and state variables for CAM4
                                          ! temperature, water vapor, cloud ice and cloud
                                          ! liquid budgets.
   integer :: history_budget_histfile_num ! output history file number for budget fields
   integer :: err
   integer :: dtime  ! time step

   logical :: use_SPCAM  ! SPCAM flag

   character(len=128) :: error_message

   ! Dummy ty_gas_concs object for call to set_gases; does not seem to be used
   ! in that method at all.
   type(ty_gas_concs) :: available_gases

   !-----------------------------------------------------------------------

   ! Initialize output fields for offline driver.
   ! TODO: do we need to keep this functionality? Where is the offline driver?
   ! Do we need to write a new offline driver for RRTMGP?
   call init_rad_data() 

   ! Read gas optics coefficients from file; NOTE: available_gases is a dummy
   ! argument, and is not actually used by the initialization routine. If it was
   ! we would have some more trouble here, because we would need to use the
   ! rad_constituents interface to initialize the gas concentrations object
   ! properly via the set_gas_concentrations routine added to this module.
   call rrtmgp_load_coefficients(k_dist_sw, coefficients_file_sw, available_gases)
   call rrtmgp_load_coefficients(k_dist_lw, coefficients_file_lw, available_gases)

   ! TODO: check that things got loaded properly?


   ! Get number of bands used in shortwave and longwave and set module data
   ! appropriately so that these sizes can be used to allocate array sizes.
   nswbands = k_dist_sw%get_nband()
   nlwbands = k_dist_lw%get_nband()

   ! Initialize the shortwave and longwave drivers
   call handle_rrtmgp_error(rrtmgp_sw_init())
   call handle_rrtmgp_error(rrtmgp_lw_init())

   ! Set the radiation timestep for cosz calculations if requested using the 
   ! adjusted iradsw value from radiation
   if (use_rad_dt_cosz)  then
      dtime  = get_step_size()
      dt_avg = iradsw*dtime
   end if

   call phys_getopts(history_amwg_out = history_amwg,    &
                     history_vdiag_out = history_vdiag,   &
                     history_budget_out = history_budget,  &
                     history_budget_histfile_num_out = history_budget_histfile_num)
   
   ! Determine whether modal aerosols are affecting the climate, and if so
   ! then initialize the modal aerosol optics module
   call rad_cnst_get_info(0, nmodes=nmodes)
   if (nmodes > 0) call modal_aer_opt_init()

   call hirsbt_init()

   ! "irad_always" is number of time steps to execute radiation continuously from start of
   ! initial OR restart run
   nstep = get_nstep()
   if (irad_always > 0) then
      irad_always = irad_always + nstep
   end if

   ! Initialize the satellite simulator package (COSP). Should this be moved to
   ! a higher level? This should probably not be specific to a given radiation
   ! package.
   if (docosp) then
      ! Initialization for the simulator package.
      call cospsimulator_intr_init
   
      ! Allocate and initialize the counter that is used to determine when to
      ! call the simulator package.
      allocate(cosp_cnt(begchunk:endchunk))
      if (is_first_restart_step()) then
         cosp_cnt(begchunk:endchunk) = cosp_cnt_init
      else
         cosp_cnt(begchunk:endchunk) = 0     
      end if
   end if

   ! Shortwave radiation
   call addfld('TOT_CLD_VISTAU', (/ 'lev' /), 'A',   '1', 'Total gbx cloud extinction visible sw optical depth', &
                                                      sampling_seq='rad_lwsw', flag_xyfill=.true.)
   call addfld('TOT_ICLD_VISTAU', (/ 'lev' /), 'A',  '1', 'Total in-cloud extinction visible sw optical depth', &
                                                      sampling_seq='rad_lwsw', flag_xyfill=.true.)
   call addfld('LIQ_ICLD_VISTAU', (/ 'lev' /), 'A',  '1', 'Liquid in-cloud extinction visible sw optical depth', &
                                                      sampling_seq='rad_lwsw', flag_xyfill=.true.)
   call addfld('ICE_ICLD_VISTAU', (/ 'lev' /), 'A',  '1', 'Ice in-cloud extinction visible sw optical depth', &
                                                      sampling_seq='rad_lwsw', flag_xyfill=.true.)

   call add_default('TOT_CLD_VISTAU',  1, ' ')
   call add_default('TOT_ICLD_VISTAU', 1, ' ')

   ! get list of active radiation calls
   call rad_cnst_get_call_list(active_calls)

   ! TODO: what is active_calls and diag?
   do icall = 0,N_DIAG

      if (active_calls(icall)) then
         call addfld('SOLIN'//diag(icall), horiz_only, 'A',   'W/m2', &
                     'Solar insolation', &
                     sampling_seq='rad_lwsw')
         call addfld('SOLL'//diag(icall),  horiz_only, 'A',   'W/m2', &
                     'Solar downward near infrared direct  to surface', &
                     sampling_seq='rad_lwsw')
         call addfld('SOLS'//diag(icall),  horiz_only, 'A',   'W/m2', &
                     'Solar downward visible direct  to surface', &
                     sampling_seq='rad_lwsw')
         call addfld('SOLLD'//diag(icall), horiz_only, 'A',   'W/m2', &
                     'Solar downward near infrared diffuse to surface', &
                     sampling_seq='rad_lwsw')
         call addfld('SOLSD'//diag(icall), horiz_only, 'A',   'W/m2', &
                     'Solar downward visible diffuse to surface', &
                     sampling_seq='rad_lwsw')
         call addfld('QRS'//diag(icall),   (/ 'lev' /), 'A',  'K/s', &
                     'Solar heating rate', &
                     sampling_seq='rad_lwsw')
         call addfld('QRSC'//diag(icall),   (/ 'lev' /), 'A', 'K/s', &
                     'Clearsky solar heating rate', &
                     sampling_seq='rad_lwsw')
         call addfld('FSNS'//diag(icall),   horiz_only,  'A', 'W/m2', &
                     'Net solar flux at surface', &
                     sampling_seq='rad_lwsw')
         call addfld('FSNT'//diag(icall),   horiz_only,  'A', 'W/m2', &
                     'Net solar flux at top of model', &
                     sampling_seq='rad_lwsw')
         call addfld('FSNTOA'//diag(icall), horiz_only,  'A', 'W/m2', &
                     'Net solar flux at top of atmosphere', &
                     sampling_seq='rad_lwsw')
         call addfld('FSUTOA'//diag(icall), horiz_only,  'A',  'W/m2', &
                     'Upwelling solar flux at top of atmosphere', &
                     sampling_seq='rad_lwsw')
         call addfld('FSNTOAC'//diag(icall), horiz_only, 'A',  'W/m2', &
                     'Clearsky net solar flux at top of atmosphere', &
                     sampling_seq='rad_lwsw')
         call addfld('FSUTOAC'//diag(icall), horiz_only, 'A',  'W/m2', &
                     'Clearsky upwelling solar flux at top of atmosphere', &
                     sampling_seq='rad_lwsw')
         call addfld('FSN200'//diag(icall), horiz_only,  'A',  'W/m2', &
                     'Net shortwave flux at 200 mb', &
                     sampling_seq='rad_lwsw')
         call addfld('FSN200C'//diag(icall), horiz_only, 'A',  'W/m2', &
                     'Clearsky net shortwave flux at 200 mb', &
                     sampling_seq='rad_lwsw')
         call addfld('FSNTC'//diag(icall), horiz_only,   'A',  'W/m2', &
                     'Clearsky net solar flux at top of model', &
                     sampling_seq='rad_lwsw')
         call addfld('FSNSC'//diag(icall), horiz_only,   'A',  'W/m2', &
                     'Clearsky net solar flux at surface', &
                     sampling_seq='rad_lwsw')
         call addfld('FSDSC'//diag(icall), horiz_only,   'A',  'W/m2', &
                     'Clearsky downwelling solar flux at surface', &
                     sampling_seq='rad_lwsw')
         call addfld('FSDS'//diag(icall), horiz_only,    'A',  'W/m2', &
                     'Downwelling solar flux at surface', &
                     sampling_seq='rad_lwsw')
         call addfld('FUS'//diag(icall),  (/ 'ilev' /),  'I',  'W/m2', &
                     'Shortwave upward flux')
         call addfld('FDS'//diag(icall),  (/ 'ilev' /),  'I',  'W/m2', &
                     'Shortwave downward flux')
         call addfld('FUSC'//diag(icall),  (/ 'ilev' /), 'I',  'W/m2', &
                     'Shortwave clear-sky upward flux')
         call addfld('FDSC'//diag(icall),  (/ 'ilev' /), 'I',  'W/m2', &
                     'Shortwave clear-sky downward flux')
         call addfld('FSNIRTOA'//diag(icall), horiz_only, 'A', 'W/m2', &
                     'Net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere', &
                     sampling_seq='rad_lwsw')
         call addfld('FSNRTOAC'//diag(icall), horiz_only, 'A', 'W/m2', &
                     'Clearsky net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere', &
                     sampling_seq='rad_lwsw')
         call addfld('FSNRTOAS'//diag(icall), horiz_only, 'A', 'W/m2', &
                     'Net near-infrared flux (>= 0.7 microns) at top of atmosphere', &
                     sampling_seq='rad_lwsw')
         call addfld('SWCF'//diag(icall),     horiz_only, 'A', 'W/m2', &
                     'Shortwave cloud forcing', &
                     sampling_seq='rad_lwsw')

         if (history_amwg) then
            call add_default('SOLIN'//diag(icall),   1, ' ')
            call add_default('QRS'//diag(icall),     1, ' ')
            call add_default('FSNS'//diag(icall),    1, ' ')
            call add_default('FSNT'//diag(icall),    1, ' ')
            call add_default('FSNTOA'//diag(icall),  1, ' ')
            call add_default('FSUTOA'//diag(icall),  1, ' ')
            call add_default('FSNTOAC'//diag(icall), 1, ' ')
            call add_default('FSUTOAC'//diag(icall), 1, ' ')
            call add_default('FSNTC'//diag(icall),   1, ' ')
            call add_default('FSNSC'//diag(icall),   1, ' ')
            call add_default('FSDSC'//diag(icall),   1, ' ')
            call add_default('FSDS'//diag(icall),    1, ' ')
            call add_default('SWCF'//diag(icall),    1, ' ')
         end if  ! history_amwg
      end if  ! active_calls(icall)
   end do  ! icall = 0,N_DIAG

   if (single_column .and. scm_crm_mode) then
      call add_default ('FUS     ', 1, ' ')
      call add_default ('FUSC    ', 1, ' ')
      call add_default ('FDS     ', 1, ' ')
      call add_default ('FDSC    ', 1, ' ')
   end if

   ! Longwave radiation
   do icall = 0,N_DIAG
      if (active_calls(icall)) then
         call addfld('QRL'//diag(icall),     (/'lev'/),  'A', 'K/s',  &
                     'Longwave heating rate',                         &
                     sampling_seq='rad_lwsw')
         call addfld('QRLC'//diag(icall),    (/'lev'/),  'A', 'K/s',  &
                     'Clearsky longwave heating rate',                &
                     sampling_seq='rad_lwsw')
         call addfld('FLDS'//diag(icall),    horiz_only, 'A', 'W/m2', &
                     'Downwelling longwave flux at surface',          &
                     sampling_seq='rad_lwsw')
         call addfld('FLDSC'//diag(icall),   horiz_only, 'A', 'W/m2', &
                     'Clearsky Downwelling longwave flux at surface', &
                     sampling_seq='rad_lwsw')
         call addfld('FLNS'//diag(icall),    horiz_only, 'A', 'W/m2', &
                     'Net longwave flux at surface',                  &
                     sampling_seq='rad_lwsw')
         call addfld('FLNT'//diag(icall),    horiz_only, 'A', 'W/m2', &
                     'Net longwave flux at top of model',             &
                     sampling_seq='rad_lwsw')
         call addfld('FLUT'//diag(icall),    horiz_only, 'A', 'W/m2', &
                     'Upwelling longwave flux at top of model',       &
                     sampling_seq='rad_lwsw')
         call addfld('FLUTC'//diag(icall),   horiz_only, 'A', 'W/m2', &
                     'Clearsky upwelling longwave flux at top of model', &
                     sampling_seq='rad_lwsw')
         call addfld('FLNTC'//diag(icall),   horiz_only, 'A', 'W/m2', &
                     'Clearsky net longwave flux at top of model',    &
                     sampling_seq='rad_lwsw')
         call addfld('LWCF'//diag(icall),    horiz_only, 'A', 'W/m2', &
                     'Longwave cloud forcing',                        &
                     sampling_seq='rad_lwsw')
         call addfld('FLN200'//diag(icall),  horiz_only, 'A', 'W/m2', &
                     'Net longwave flux at 200 mb',                   &
                     sampling_seq='rad_lwsw')
         call addfld('FLN200C'//diag(icall), horiz_only, 'A', 'W/m2', &
                     'Clearsky net longwave flux at 200 mb',          &
                     sampling_seq='rad_lwsw')
         call addfld('FLNSC'//diag(icall),   horiz_only, 'A', 'W/m2', &
                     'Clearsky net longwave flux at surface',         &
                     sampling_seq='rad_lwsw')
         call addfld('FUL'//diag(icall),     (/'ilev'/), 'I', 'W/m2', &
                     'Longwave upward flux')
         call addfld('FDL'//diag(icall),     (/'ilev'/), 'I', 'W/m2', &
                     'Longwave downward flux')
         call addfld('FULC'//diag(icall),    (/'ilev'/), 'I', 'W/m2', &
                     'Longwave clear-sky upward flux')
         call addfld('FDLC'//diag(icall),    (/'ilev'/), 'I', 'W/m2', &
                     'Longwave clear-sky downward flux')

         if (history_amwg) then
            call add_default('QRL'//diag(icall),   1, ' ')
            call add_default('FLNS'//diag(icall),  1, ' ')
            call add_default('FLDS'//diag(icall),  1, ' ')
            call add_default('FLNT'//diag(icall),  1, ' ')
            call add_default('FLUT'//diag(icall),  1, ' ')
            call add_default('FLUTC'//diag(icall), 1, ' ')
            call add_default('FLNTC'//diag(icall), 1, ' ')
            call add_default('FLNSC'//diag(icall), 1, ' ')
            call add_default('LWCF'//diag(icall),  1, ' ')
         endif  ! history_amwg
      end if
   end do

   call addfld('EMIS', (/ 'lev' /), 'A', '1', 'Cloud longwave emissivity')

   ! Add default fields for single column mode
   if (single_column.and.scm_crm_mode) then
      call add_default ('FUL     ', 1, ' ')
      call add_default ('FULC    ', 1, ' ')
      call add_default ('FDL     ', 1, ' ')
      call add_default ('FDLC    ', 1, ' ')
   endif

   ! HIRS/MSU diagnostic brightness temperatures
   if (dohirs) then
      call addfld (hirsname(1),horiz_only,'A','K','HIRS CH2 infra-red brightness temperature')
      call addfld (hirsname(2),horiz_only,'A','K','HIRS CH4 infra-red brightness temperature')
      call addfld (hirsname(3),horiz_only,'A','K','HIRS CH6 infra-red brightness temperature')
      call addfld (hirsname(4),horiz_only,'A','K','HIRS CH8 infra-red brightness temperature')
      call addfld (hirsname(5),horiz_only,'A','K','HIRS CH10 infra-red brightness temperature')
      call addfld (hirsname(6),horiz_only,'A','K','HIRS CH11 infra-red brightness temperature')
      call addfld (hirsname(7),horiz_only,'A','K','HIRS CH12 infra-red brightness temperature')
      call addfld (msuname(1),horiz_only,'A','K','MSU CH1 microwave brightness temperature')
      call addfld (msuname(2),horiz_only,'A','K','MSU CH2 microwave brightness temperature')
      call addfld (msuname(3),horiz_only,'A','K','MSU CH3 microwave brightness temperature')
      call addfld (msuname(4),horiz_only,'A','K','MSU CH4 microwave brightness temperature')

      ! Add HIRS/MSU fields to default history files
      call add_default (hirsname(1), 1, ' ')
      call add_default (hirsname(2), 1, ' ')
      call add_default (hirsname(3), 1, ' ')
      call add_default (hirsname(4), 1, ' ')
      call add_default (hirsname(5), 1, ' ')
      call add_default (hirsname(6), 1, ' ')
      call add_default (hirsname(7), 1, ' ')
      call add_default (msuname(1), 1, ' ')
      call add_default (msuname(2), 1, ' ')
      call add_default (msuname(3), 1, ' ')
      call add_default (msuname(4), 1, ' ')
   end if

   ! Heating rate needed for d(theta)/dt computation
   call addfld ('HR',(/ 'lev' /), 'A','K/s','Heating rate needed for d(theta)/dt computation')

   if (history_budget .and. history_budget_histfile_num > 1) then
      call add_default ('QRL     ', history_budget_histfile_num, ' ')
      call add_default ('QRS     ', history_budget_histfile_num, ' ')
   end if

   if (history_vdiag) then
      call add_default('FLUT', 2, ' ')
      call add_default('FLUT', 3, ' ')
   end if

   if (history_budget .and. history_budget_histfile_num > 1) then
      call add_default ('QRL     ', history_budget_histfile_num, ' ')
      call add_default ('QRS     ', history_budget_histfile_num, ' ')
   end if

   if (history_vdiag) then
      call add_default('FLUT', 2, ' ')
      call add_default('FLUT', 3, ' ')
   end if 

   cldfsnow_idx = pbuf_get_index('CLDFSNOW',errcode=err)
   cld_idx      = pbuf_get_index('CLD')
   concld_idx   = pbuf_get_index('CONCLD')
   rel_idx      = pbuf_get_index('REL')
   rei_idx      = pbuf_get_index('REI')
   dei_idx      = pbuf_get_index('DEI')

   if (cldfsnow_idx > 0) then
      call addfld ('CLDFSNOW',(/ 'lev' /),'I','1','CLDFSNOW',flag_xyfill=.true.)
      call addfld('SNOW_ICLD_VISTAU', (/ 'lev' /), 'A', '1', &
                  'Snow in-cloud extinction visible sw optical depth', &
                  sampling_seq='rad_lwsw', flag_xyfill=.true.)
   endif

end subroutine radiation_init

!===============================================================================
  
subroutine radiation_tend(state, ptend, pbuf, cam_out, cam_in, net_flx)
!-------------------------------------------------------------------------------
! 
! Purpose: 
! Driver for radiation computation.
! 
! Method: 
! Radiation uses cgs units, so conversions must be done from
! model fields to radiation fields.
!
! Revision history:
! May 2004    D.B. Coleman     Merge of code from radctl.F90 and parts of tphysbc.F90.
! 2004-08-09  B. Eaton         Add pointer variables for constituents.
! 2004-08-24  B. Eaton         Access O3 and GHG constituents from chem_get_cnst.
! 2004-08-30  B. Eaton         Replace chem_get_cnst by rad_constituent_get.
! 2007-11-05  M. Iacono        Install rrtmg_lw and sw as radiation model.
! 2007-12-27  M. Iacono        Modify to use CAM cloud optical properties with rrtmg.
! 2009-06     Minghuai Wang,   add treatments for MMF CAM
!              These modifications are based on the spcam3.5, which was developled
!              by Marat Khairoutdinov (mkhairoutdin@ms.cc.sunysb.edu). The spcam3.5 
!              only have one radiation package (camrt). 
!              Short wave and long wave radiation codes are called for every column
!              of the CRM domain. CRM fields are named as *_crm, and domain-averaged fields are
!              named as *_m. The domain-averaged fields and CRM fields are outputed at the end
!              of the CRM domain loop (last_column=.true.).
!              Several variables in state are updated with those from CRM output
!              (liquid water, qc_rad; ice water, qi_rad; water vapor, qv_rad; 
!              and temperature, t_rad). Several variables in pbuf are also updated 
!              with those in CRM domain (cld, cicewp, cliqwp, csnowp, cldfsnow).  
!              At the end of the radiation calculation, state and those in pbuf are
!              restored to the old values. 
!              Finally, a new cloud simulator are called, which takes account of cloud fileds 
!              from the CRM domain. 
! 2009-07-13, Minghuai Wang: MMF CAM
!             When Morrison's two momenent microphysics is used in SAM, droplet and ice crystal effective radius
!             used in this radiation code are calcualted at each CRM column by calling m2005_effradius
! 2009-10-21, Minghuai Wang: MMF CAM
!             CRM-scale aerosol water is used to calculate aerosol optical depth
! 2017-11     Ben Hillman: Adapted for use with RRTMGP
!----------------------------------------------------------------------------------------------

   ! Performance module needed for timing functions
   use perf_mod,        only: t_startf, t_stopf

   ! CAM derived types; needed for surface exchange fields, physics state, and
   ! tendency fields
   use camsrfexch,      only: cam_out_t, cam_in_t
   use physics_types,   only: physics_state, physics_ptend

   ! Utilities for interacting with the physics buffer
   use physics_buffer,  only: physics_buffer_desc, pbuf_get_field, &
                              pbuf_get_index

   ! For calculating radiative heating tendencies?
   use radheat,         only: radheat_tend

   use physconst,       only: cpair, cappa

   use rad_constituents, only: N_DIAG, rad_cnst_get_call_list

   ! These are needed for calculating solar zenith angle
   use phys_grid,       only: get_rlat_all_p, get_rlon_all_p
   use time_manager,    only: get_curr_calday
   use orbit,           only: zenith

   ! RRTMGP radiation drivers
   use mo_rrtmgp_clr_all_sky, only: rrtmgp_sw=>rte_sw, rrtmgp_lw=>rte_lw
   use mo_gas_concentrations, only: ty_gas_concs
   use mo_optical_props, only: ty_optical_props, ty_optical_props_1scl
   use mo_fluxes_byband, only:  ty_fluxes_byband

   !use cloud_rad_props,  only: get_ice_optics_sw, get_liquid_optics_sw, &
   !                            liquid_cloud_get_rad_props_lw, &
   !                            ice_cloud_get_rad_props_lw, cloud_rad_props_get_lw, &
   !                            snow_cloud_get_rad_props_lw, get_snow_optics_sw
   !use slingo,           only: slingo_liq_get_rad_props_lw, slingo_liq_optics_sw
   !use ebert_curry,      only: ec_ice_optics_sw, ec_ice_get_rad_props_lw
   !use rad_solar_var,    only: get_variability

   ! Arguments

   ! Physics state variables
   type(physics_state), intent(in), target :: state

   ! Heating tendencies calculated in this subroutine
   type(physics_ptend), intent(out) :: ptend

   ! Fields from other parameterizations that persist across timesteps
   type(physics_buffer_desc), pointer :: pbuf(:)

   ! Surface fluxes?
   type(cam_out_t), intent(inout) :: cam_out
   type(cam_in_t), intent(in) :: cam_in

   ! Net flux calculated in this routine; used to check energy conservation in
   ! the physics package driver?
   real(r8), intent(inout) :: net_flx(pcols)

   ! ---------------------------------------------------------------------------
   ! Local variables
   logical :: dosw, dolw
   integer :: nstep  ! current timestep number

   real(r8) :: heating_rate(pcols,pver)  ! Net heating rate

!   ! Combined cloud radiative parameters (liquid, ice, and possibly snow). Note
!   ! that these are "in cloud" not "in cell"
!   real(r8) :: c_cld_tau    (nswbands,pcols,pver) ! cloud extinction optical depth
!   real(r8) :: c_cld_tau_w  (nswbands,pcols,pver) ! cloud single scattering albedo * tau
!   real(r8) :: c_cld_tau_w_g(nswbands,pcols,pver) ! cloud assymetry parameter * w * tau
!   real(r8) :: c_cld_tau_w_f(nswbands,pcols,pver) ! cloud forward scattered fraction * w * tau
!   real(r8) :: c_cld_lw_abs (nlwbands,pcols,pver) ! cloud absorption optics depth (LW)
!
!   ! Cloud radiative parameters (liquid and ice cloud). Note that these are are
!   ! "in cloud" not "in cell"
!   real(r8) :: cld_tau    (nswbands,pcols,pver) ! cloud extinction optical depth
!   real(r8) :: cld_tau_w  (nswbands,pcols,pver) ! cloud single scattering albedo * tau
!   real(r8) :: cld_tau_w_g(nswbands,pcols,pver) ! cloud assymetry parameter * w * tau
!   real(r8) :: cld_tau_w_f(nswbands,pcols,pver) ! cloud forward scattered fraction * w * tau
!   real(r8) :: cld_lw_abs (nlwbands,pcols,pver) ! cloud absorption optics depth (LW)
!
!   ! Ice cloud radiative parameters. Note that these are "in cloud" not "in cell".
!   real(r8) :: ice_tau    (nswbands,pcols,pver) ! ice extinction optical depth
!   real(r8) :: ice_tau_w  (nswbands,pcols,pver) ! ice single scattering albedo * tau
!   real(r8) :: ice_tau_w_g(nswbands,pcols,pver) ! ice assymetry parameter * tau * w
!   real(r8) :: ice_tau_w_f(nswbands,pcols,pver) ! ice forward scattered fraction * tau * w
!   real(r8) :: ice_lw_abs (nlwbands,pcols,pver) ! ice absorption optics depth (LW)
!
!   ! "Snow" cloud radiative parameters. Note that these are "in cloud" not "in cell".
!   real(r8) :: snow_tau    (nswbands,pcols,pver) ! snow extinction optical depth
!   real(r8) :: snow_tau_w  (nswbands,pcols,pver) ! snow single scattering albedo * tau
!   real(r8) :: snow_tau_w_g(nswbands,pcols,pver) ! snow assymetry parameter * tau * w
!   real(r8) :: snow_tau_w_f(nswbands,pcols,pver) ! snow forward scattered fraction * tau * w
!   real(r8) :: snow_lw_abs (nlwbands,pcols,pver) ! snow absorption optics depth (LW)
!
!   ! Liquid cloud radiative parameters. Note that these are "in cloud" not "in cell".
!   real(r8) :: liq_tau    (nswbands,pcols,pver) ! liquid extinction optical depth
!   real(r8) :: liq_tau_w  (nswbands,pcols,pver) ! liquid single scattering albedo * tau
!   real(r8) :: liq_tau_w_g(nswbands,pcols,pver) ! liquid assymetry parameter * tau * w
!   real(r8) :: liq_tau_w_f(nswbands,pcols,pver) ! liquid forward scattered fraction * tau * w
!   real(r8) :: liq_lw_abs (nlwbands,pcols,pver) ! liquid absorption optics depth (LW)
!
!   ! Aerosol radiative properties; why the heck are these ordered differently
!   ! from the cloud radiative properties? Is this a carry-over from old RRTMG
!   ! that we no longer need to perpetuate?
!   real(r8) :: aer_tau    (pcols,0:pver,nswbands) ! aerosol extinction optical depth
!   real(r8) :: aer_tau_w  (pcols,0:pver,nswbands) ! aerosol single scattering albedo * tau
!   real(r8) :: aer_tau_w_g(pcols,0:pver,nswbands) ! aerosol assymetry parameter * w * tau
!   real(r8) :: aer_tau_w_f(pcols,0:pver,nswbands) ! aerosol forward scattered fraction * w * tau
!   real(r8) :: aer_lw_abs (pcols,pver,nlwbands)   ! aerosol absorption optics depth (LW)

   ! Pointers to heating rates on physics buffer
   real(r8), pointer :: qrs(:,:) => null()  ! shortwave radiative heating rate 
   real(r8), pointer :: qrl(:,:) => null()  ! longwave  radiative heating rate 

   ! Variables used to calculate solar zenith angle for the current timestep.
   real(r8) :: calday         ! current calendar day
   real(r8) :: clat(pcols)    ! current latitudes(radians)
   real(r8) :: clon(pcols)    ! current longitudes(radians)
   real(r8) :: coszrs(pcols)  ! Cosine solar zenith angle

   ! Flag to carry (QRS,QRL)*dp across time steps. 
   ! TODO: what does this mean?
   logical :: conserve_energy = .true.

   ! Loop variables; i is used to loop over columns, and k is used to loop over
   ! levels
   integer :: i, k
   integer :: icall

   integer :: ncol, nlay

   ! Indices to highest and lowest CAM levels with valid pressure values to be
   ! used in the radiative transfer parameterization
   integer :: ktop, kbot

   ! Gathered indicies of day and night columns 
   ! chunk_column_index = day_indices(daylight_column_index)
   integer :: nday     ! Number of daylight columns
   integer, allocatable :: day_indices(:)   ! Indicies of daylight coumns

   ! Set name for this subroutine for log files
   character(*), parameter :: subroutine_name = 'radiation_tend'

   ! For loops over diagnostic calls (TODO: what does this mean?)
   logical :: active_calls(0:N_DIAG)

   real(r8), allocatable :: surface_emissivity(:,:)

   ! State fields that are passed into RRTMGP. Some of these may need to
   ! modified from what exist in the physics_state object, i.e. to clip
   ! temperatures to make sure they are within the valid range.
   real(r8), allocatable :: t_rad(:,:)  ! Temperature
   real(r8), allocatable :: pmid_rad(:,:)  ! Pressure at level midpoints
   real(r8), allocatable :: pint_rad(:,:)  ! Pressure at level interfaces
   
   ! RRTMGP types
   type(ty_gas_concs) :: gas_concentrations_sw, gas_concentrations_lw
   type(ty_optical_props_1scl) :: aerosol_optics_lw
   type(ty_optical_props_1scl) :: cloud_optics_lw
   type(ty_fluxes_byband) :: flux_lw_allsky, flux_lw_clearsky

   !----------------------------------------------------------------------
   
   ! Number of physics columns in this "chunk"; used in multiple places
   ! throughout this subroutine, so set once for convenience
   ncol = state%ncol

   ! Figure out which levels to use; only use those levels between the min and
   ! max of press_ref read from the coefficients file. Below code would do this
   ! on a column by column basis, but there are smarter ways of doing this, and
   ! I do not think we want to look over columns here, but rather push the loop
   ! over columns down further into the kernels. I think the implementation that
   ! Brian Eaton at NCAR used was to look at the reference pressure pref. We
   ! could also look at the min and max over all columns in this chunk. That
   ! might be the most transparent thing to do.

   ! Find highest level where all pressures are greater than minimum
   do k = 1,pverp
      if (all(state%pint(:ncol,k) > k_dist_lw%get_press_ref_min())) then
         ktop = k
         exit
      end if
   end do

   ! Find lowest level where all pressures are less than maximum. Starting from
   ! the bottom-most model level (pverp), take the first level where all
   ! pressures are above maximum.
   do k = pverp,1,-1
      if (all(state%pint(:ncol,k) < k_dist_lw%get_press_ref_max())) then
         kbot = k
         exit
      end if
   end do
         
   ! Size of grid for radiation; note kbot > ktop because CAM grid goes from TOA
   ! to surface. RRTMGP grid will go from surface to TOA internally, but this is
   ! handled automatically so we do not worry about this here.
   nlay = kbot - ktop

   ! Sanity check on nlay
   if (nlay <= 0) then
      call endrun(module_name // ': no valid pressure levels?')
   end if
   
   ! Set pointers to heating rates stored on physics buffer. These will be
   ! modified in this routine.
   call pbuf_get_field(pbuf, pbuf_get_index('QRS'), qrs)
   call pbuf_get_field(pbuf, pbuf_get_index('QRL'), qrl)
  
   ! Check if we are supposed to do shortwave and longwave stuff at this
   ! timestep, and if we are then we begin setting optical properties and then
   ! doing the radiative transfer separately for shortwave and longwave.
   dosw = radiation_do('sw')
   dolw = radiation_do('lw')
   if (dosw .or. dolw) then

      ! Loop over diagnostic calls (TODO: more documentation on what this
      ! means)
      ! get list of active radiation calls
      call rad_cnst_get_call_list(active_calls)

      ! Do shortwave stuff...
      if (dosw) then

         ! Get cosine solar zenith angle for current time step. If the swrad_off flag
         ! is set, meaning we should not do SW radiation, then we just set coszrs to
         ! zero everywhere.
         if (swrad_off) then
            coszrs(:) = 0._r8
         else
            calday = get_curr_calday()
            call get_rlat_all_p(state%lchnk, ncol, clat)
            call get_rlon_all_p(state%lchnk, ncol, clon)
            call zenith(calday, clat, clon, coszrs, ncol, dt_avg)
         end if

         ! Gather night/day column indices for subsetting SW inputs; we only want to
         ! do the shortwave radiative transfer during the daytime to save
         ! computational cost (and because RRTMGP will fail for cosine solar zenith
         ! angles less than or equal to zero)
         allocate(day_indices(ncol))
         nday = 0
         do i = 1,ncol
            if (coszrs(i) > 0.0_r8) then
               nday = nday + 1
               day_indices(nday) = i
            else
            end if
         end do

         ! State variables only as large as nday
         allocate(t_rad(nday,play), pmid_rad(nday,play), pint_rad(nday,play+1))
         do i = 1,nday
            t_rad(i,:play) = state%t(day_indices(i),ktop:kbot)
            pmid_rad(i,:play) = state%pmid(day_indices(i),ktop:kbot)
            pint_rad(i,:play+1) = state%pmid(day_indices(i),ktop:kbot+1)
         end do

         ! Do shortwave cloud optics calculations
         call t_startf('shortwave cloud optics')
         !call set_cloud_optics_sw(state, pbuf, cloud_sw)
         call t_stopf('shortwave cloud optics')

         ! Get shortwave aerosol optics
         call t_startf('shortwave aerosol optics')
         !call set_aerosol_optics_sw(state, pbuf, aerosol_sw)
         call t_stopf('shortwave aerosol optics')

         ! Subset optical properties to get only daytime columns
         !call subset_daytime_optics()

         ! The climate (icall==0) calculation must occur last, so we loop
         ! backwards.
         do icall = N_DIAG,0,-1
            if (active_calls(icall)) then

               ! Set gas concentrations (I believe the active gases may change
               ! for different values of icall, which is why we do this within
               ! the loop).
               call set_gas_concentrations(icall, state, pbuf, pver, &
                                           gas_concentrations_sw, &
                                           n_day_indices=nday, &
                                           day_indices=day_indices)

               ! Do shortwave radiative transfer calculations
               call t_startf('shortwave radiation calculations')
               !errmsg = rrtmgp_sw(k_dist_sw, gas_concentrations_sw, &
               !                   pmid_day, t_day, pint_day, &
               !                   coszrs_day, alb_dir, alb_dif, cloud_sw, &
               !                   allsky_fluxes_sw, clrsky_fluxes_sw, &
               !                   aer_props=aerosol_sw)
               call t_stopf('shortwave radiation calculations')

               ! Map RRTMGP outputs to CAM outputs
               ! TODO
            end if
         end do

         ! Clean up
         deallocate(t_rad, pmid_rad, pint_rad, day_indices)

      end if  ! dosw

      ! Do longwave stuff...
      if (dolw) then

         ! initialize aerosol optics object
         ! NOTE: we can initialize the aerosol optical properties by ngpt or by
         ! nlwbands. We initialize by nlwbands here, because that is what the
         ! internal CAM aer_rad_props routines expect. There is logic
         ! encapsulated in RRTMGP to figure out how the optical properties are
         ! split up.
         call handle_rrtmgp_error(aerosol_optics_lw%init_1scl( &
            ncol, nlay, nlwbands, &
            name='longwave aerosol optics' &
         ))

         ! initialize cloud optics object
         call handle_rrtmgp_error(cloud_optics_lw%init_1scl( &
            ncol, nlay, k_dist_lw%get_ngpt(), &
            name='longwave cloud optics' &
         ))

         ! Manually set cloud optical depth to zero since we do not have a
         ! function to set appropriately yet. TODO: move this to the
         ! "set_cloud_optics_lw" routine.
         cloud_optics_lw%tau(:,:,:) = 0.0
         call handle_rrtmgp_error(cloud_optics_lw%validate())

         ! Set surface emissivity to 1 here
         ! TODO: set this more intelligently?
         allocate(surface_emissivity(nlwbands,ncol))
         surface_emissivity(:,:) = 1.0_r8

         ! Allocate longwave outputs; why is this not part of the
         ! ty_fluxes_byband object?
         ! NOTE: fluxes defined at interfaces, so initialize to have vertical
         ! dimension pver+1
         call initialize_rrtmgp_fluxes( &
            ncol, nlay+1, nlwbands, &
            flux_lw_allsky &
         )
         call initialize_rrtmgp_fluxes( &
            ncol, nlay+1, nlwbands, &
            flux_lw_clearsky &
         )
            
         ! Do longwave cloud optics calculations
         !call t_startf('longwave cloud optics')
         !call set_cloud_optics_lw(state, pbuf, cloud_lw)
         !call t_stopf('longwave cloud optics')

         ! Allocate arrays to copy physics state fields to (need to do this
         ! because we may need to modify them, and we do not want to modify
         ! the physics state from within parameterizations)
         allocate(t_rad(ncol,nlay), pmid_rad(ncol,nlay), pint_rad(ncol,nlay))

         ! Loop over diagnostic calls (what does this mean?)
         do icall = N_DIAG,0,-1
            if (active_calls(icall)) then

               ! Set gas concentrations (I believe the active gases may change
               ! for different values of icall, which is why we do this within
               ! the loop).
               call t_startf('longwave gas concentrations')
               call set_gas_concentrations(icall, state, pbuf, pver, &
                                           gas_concentrations_lw)
               call t_stopf('longwave gas concentrations')

               ! Get longwave aerosol optics
               call t_startf('longwave aerosol optics')
               call set_aerosol_optics_lw(icall, state, pbuf, aerosol_optics_lw)
               call t_stopf('longwave aerosol optics')

               ! Populate state variables
               t_rad(:ncol,:nlay) = state%t(:ncol,ktop:kbot)
               pmid_rad(:ncol,:nlay) = state%pmid(:ncol,ktop:kbot)
               pint_rad(:ncol,:nlay+1) = state%pint(:ncol,ktop:kbot+1)

               ! Clip temperature values in case they are out of range
               where (t_rad <  k_dist_lw%get_temp_ref_min()) 
                  t_rad = k_dist_lw%get_temp_ref_min()
               endwhere
               where (t_rad > k_dist_lw%get_temp_ref_max())
                  t_rad = k_dist_lw%get_temp_ref_max()
               endwhere

               ! Do longwave radiative transfer calculations
               call t_startf('longwave radiation calculations')
               call handle_rrtmgp_error(rrtmgp_lw( &
                  k_dist_lw, gas_concentrations_lw, &
                  pmid_rad(:ncol,:nlay), t_rad(:ncol,:nlay), &
                  pint_rad(:ncol,:nlay+1), cam_in%ts(:ncol), &
                  surface_emissivity(:nlwbands,:ncol), &
                  cloud_optics_lw, &
                  flux_lw_allsky, flux_lw_clearsky &
                  !aer_props=aerosol_optics_lw &
               ))
               call t_stopf('longwave radiation calculations')

               deallocate(t_rad)

               ! Calculate longwave heating rate
               call calculate_heating_rate(ktop, kbot, state, flux_lw_allsky, qrl)

               ! Write outputs to history tapes
               call radiation_output_lw( &
                  icall, ktop, kbot, state, flux_lw_allsky, flux_lw_clearsky, qrl &
               )

            end if  ! active calls
         end do  ! loop over diagnostic calls
               
      end if  ! dolw

   else  ! dosw .or. dolw
      ! convert radiative heating rates from Q*dp to Q for energy conservation
      if (conserve_energy) then
         do k = 1,pver
            do i = 1,ncol
               qrs(i,k) = qrs(i,k)/state%pdel(i,k)
               qrl(i,k) = qrl(i,k)/state%pdel(i,k)
            end do
         end do
      end if  
   end if  ! dosw .or. dolw

   ! Compute net radiative heating tendency
   call t_startf('radheat_tend')
   !call radheat_tend(state, pbuf,  ptend, qrl, qrs, fsns, &
   !                  fsnt, flns, flnt, cam_in%asdir, net_flx)
   call t_stopf('radheat_tend')

   ! Compute net heating rate for dtheta/dt
   ! TODO: how is this different than above?
   call t_startf('heating_rate')
   do k=1,pver
      do i = 1,ncol
         heating_rate(i,k) = (qrs(i,k) + qrl(i,k))/cpair * (1.e5_r8/state%pmid(i,k))**cappa
      end do
   end do
   call t_stopf('heating_rate')

   ! convert radiative heating rates to Q*dp for energy conservation
   ! TODO: this gets converted back here? what is going on here?
   if (conserve_energy) then
      do k = 1,pver
         do i = 1,ncol
            qrs(i,k) = qrs(i,k)*state%pdel(i,k)
            qrl(i,k) = qrl(i,k)*state%pdel(i,k)
         end do
      end do
   end if

   ! copy net sw flux to cam_out structure
   !cam_out%netsw(:ncol) = fsns(:ncol)

end subroutine radiation_tend


subroutine calculate_heating_rate(ktop, kbot, state, fluxes, heating_rate)

   use mo_fluxes_byband, only: ty_fluxes_byband
   use physics_types, only: physics_state
   use physconst, only: gravit

   integer, intent(in) :: ktop, kbot
   type(physics_state), intent(in) :: state
   type(ty_fluxes_byband), intent(in) :: fluxes
   real(r8), intent(inout) :: heating_rate(:,:)

   integer :: i, k_cam, k_rad

   ! Loop over levels and calculate heating rates; note that the fluxes *should*
   ! be defined at interfaces, so the loop over pver and grabbing the current
   ! and next value of k should be safe.
   k_rad = 1
   do k_cam = ktop,kbot
      do i = 1,state%ncol
         heating_rate(i,k_cam) = ( &
            fluxes%flux_net(i,k_rad+1) - fluxes%flux_net(i,k_rad) &
         ) * gravit / state%pdel(i,k_cam)
      end do
      k_rad = k_rad + 1
   end do

end subroutine calculate_heating_rate
      

! Send longwave fluxes and heating rates to history buffer
subroutine radiation_output_lw(icall, ktop, kbot, state, flux_allsky, flux_clearsky, qrl)
   use physconst, only: cpair
   use physics_types, only: physics_state
   use mo_fluxes_byband, only: ty_fluxes_byband
   use cam_history, only: outfld
   
   integer, intent(in) :: icall
   integer, intent(in) :: ktop, kbot  ! indices to limits of CAM grid
   type(physics_state), intent(in) :: state
   type(ty_fluxes_byband), intent(in) :: flux_allsky
   type(ty_fluxes_byband), intent(in) :: flux_clearsky
   real(r8), intent(in) :: qrl(:,:)

   real(r8) :: flut(pcols), flnt(pcols), flns(pcols), flds(pcols)
   real(r8) :: flntc(pcols), flnsc(pcols)
   integer :: ncol

   ncol = state%ncol
   
   flnt(:ncol) = flux_allsky%flux_up(:ncol,ktop) - flux_allsky%flux_dn(:ncol,ktop)
   flns(:ncol) = flux_allsky%flux_up(:ncol,kbot+1) - flux_allsky%flux_dn(:ncol,kbot+1)
   flut(:ncol) = flux_allsky%flux_up(:ncol,ktop)
   flds(:ncol) = flux_allsky%flux_dn(:ncol,kbot+1)
   
   ! Heating rates
   call outfld('QRL'//diag(icall), qrl(:ncol,:pver)/cpair, ncol, state%lchnk)

   ! All-sky fluxes
   call outfld('FLNT'//diag(icall), flnt, pcols, state%lchnk)
   call outfld('FLNS'//diag(icall), flns, pcols, state%lchnk)
   call outfld('FLUT'//diag(icall), flut, pcols, state%lchnk)
   call outfld('FLDS'//diag(icall), flds, pcols, state%lchnk)

   ! Clear-sky fluxes
   call outfld('FLNSC'//diag(icall), flnsc, pcols, state%lchnk)

end subroutine radiation_output_lw


! For some reason the RRTMGP flux objects do not include initialization
! routines, but rather expect the user to associate the individual fluxes (which
! are pointers) with appropriate targets. Rather, this routine treats those
! pointers as allocatable members and allocates space for them. TODO: is this
! appropriate use?
subroutine initialize_rrtmgp_fluxes(ncol, nlevels, nbands, fluxes)

   use mo_fluxes_byband, only: ty_fluxes_byband
   integer, intent(in) :: ncol, nlevels, nbands
   type(ty_fluxes_byband), intent(inout) :: fluxes

   ! Allocate flux arrays; fluxes defined at interfaces, so dimension nlevels+1
   allocate(fluxes%flux_up(ncol, nlevels))
   allocate(fluxes%flux_dn(ncol, nlevels))
   allocate(fluxes%flux_net(ncol, nlevels))
   allocate(fluxes%bnd_flux_up(ncol, nlevels, nbands))
   allocate(fluxes%bnd_flux_dn(ncol, nlevels, nbands))
   allocate(fluxes%bnd_flux_net(ncol, nlevels, nbands))

   ! Direct fluxes (no scattering?)
   !allocate(fluxes%flux_dn_dir(ncol, nlevels))
   !allocate(fluxes%bnd_flux_dn_dir(ncol, nlevels, nbands))

end subroutine initialize_rrtmgp_fluxes


! Populate a RRTMGP optical properties object with CAM-computed aerosol
! absorption optical depth.
subroutine set_aerosol_optics_lw(icall, state, pbuf, optical_properties)

   use physics_types, only: physics_state
   use physics_buffer, only: physics_buffer_desc
   use mo_optical_props, only: ty_optical_props_1scl
   use aer_rad_props, only: aer_rad_props_lw

   integer, intent(in) :: icall
   type(physics_state), intent(in) :: state
   type(physics_buffer_desc), pointer :: pbuf(:)
   type(ty_optical_props_1scl), intent(inout) :: optical_properties
   real(r8) :: absorption_tau(pcols,pver,nlwbands)
   character(len=*), parameter :: subroutine_name = 'set_aerosol_optics_lw'

   ! Make sure aerosol optical properties have been initialized
   if (.not. allocated(optical_properties%tau)) then
      call endrun(subroutine_name // ': optical_properties not initialized.')
      !call handle_rrtmgp_error(optical_properties%init_1scl( &
      !   state%ncol, pver+1, nlwbands, 'aerosol'  &
      !))
   end if 

   ! Get aerosol absorption optical depth from CAM routine
   absorption_tau(:,:,:) = 0.0
   call aer_rad_props_lw(icall, state, pbuf, absorption_tau)

   ! Make sure all optical depths are within range
   if (any(absorption_tau < 0)) then
      call endrun(subroutine_name // ': aerosol absorption_tau has negative values.')
   end if

   ! Populate the RRTMGP optical properties object with CAM absorption optical
   ! depth
   optical_properties%tau(:,:,:) = 0.0
   optical_properties%tau(:state%ncol,:pver,:) = absorption_tau(:state%ncol,:pver,:)

   ! Validate
   if (any(optical_properties%tau < 0)) then
      call endrun(subroutine_name // ': aerosol optical_properties%tau has negative values.')
   end if
end subroutine


subroutine set_cloud_optics_lw(state, pbuf, optical_properties)

   use physics_types, only: physics_state
   use physics_buffer, only: physics_buffer_desc
   use mo_optical_props, only: ty_optical_props_1scl

   type(physics_state), intent(in) :: state
   type(physics_buffer_desc), pointer :: pbuf(:)
   type(ty_optical_props_1scl), intent(inout) :: optical_properties
   real(r8) :: absorption_tau(pcols,pver,nlwbands)

   ! Make sure optical properties have been initialized
   if (.not. allocated(optical_properties%tau)) then
      call handle_rrtmgp_error(optical_properties%init_1scl( &
         state%ncol, pver, nlwbands, 'cloud'  &
      ))
   end if 

   ! Get optical properties using given cloud optics scheme
   absorption_tau(:,:,:) = 0.0

   ! Do subcolumn sampling for McICA

   ! Set optical_properties
   optical_properties%tau(:,:,:) = 0.0

end subroutine set_cloud_optics_lw


subroutine set_gas_concentrations(icall, state, pbuf, nlayers, &
                                  gas_concentrations, &
                                  n_day_indices, day_indices)

   use physics_types, only: physics_state
   use physics_buffer, only: physics_buffer_desc
   use rad_constituents, only: rad_cnst_get_gas
   use mo_gas_concentrations, only: ty_gas_concs

   integer, intent(in) :: icall
   type(physics_state), intent(in) :: state
   type(physics_buffer_desc), pointer :: pbuf(:)
   integer, intent(in) :: nlayers
   type(ty_gas_concs), intent(inout) :: gas_concentrations
   integer, intent(in), optional :: n_day_indices
   integer, intent(in), optional :: day_indices(:)

   ! Local variables
   real(r8), allocatable :: vol_mix_ratio(:,:)
   real(r8), pointer :: mass_mix_ratio(:,:)
   character(len=3), dimension(8) :: gas_species = (/ &
      'H2O', 'CO2', 'O3 ', 'N2O', &
      'CO ', 'CH4', 'O2 ', 'N2 ' &
   /)
   real(r8), dimension(8) :: mol_weight_gas = (/ &
      18.02, 44.01, 48.00, 44.013, &
      28.01, 16.04, 15.999, 28.0134 &
   /)  ! g/mol
   real(r8), parameter :: mol_weight_air = 28.97  ! g/mol
                                    
   ! Defaults for gases that are not available (TODO: is this still accurate?)
   real(r8), parameter :: co_vol_mix_ratio = 1.0e-7_r8
   real(r8), parameter :: n2_vol_mix_ratio = 0.7906_r8

   ! Loop index for loop over gas species
   integer :: igas, iday
   

   ! Allocate array to hold volume mixing ratios translated from CAM mass mixing
   ! ratios
   if (present(n_day_indices)) then
      allocate(vol_mix_ratio(n_day_indices,nlayers))
   else
      allocate(vol_mix_ratio(state%ncol,nlayers))
   end if

   ! For each gas species needed for RRTMGP, read the mass mixing ratio from the
   ! CAM rad_constituents interface, convert to volume mixing ratios, and
   ! subset for daytime-only indices if needed.
   do igas = 1,size(gas_species)

      select case(trim(gas_species(igas)))

         case('CO')

            ! CO not available, use default
            vol_mix_ratio(:,:) = co_vol_mix_ratio

         case('N2')

            ! N2 not available, use default
            vol_mix_ratio(:,:) = n2_vol_mix_ratio

         case('H2O')

            ! Water vapor is represented as specific humidity in CAM, so we
            ! need to handle water a little differently
            call rad_cnst_get_gas(icall, trim(gas_species(igas)), state, pbuf, &
                                  mass_mix_ratio)
            
            ! Convert to volume mixing ratio by multiplying by the ratio of
            ! molecular weight of dry air to molecular weight of gas. Note that
            ! first specific humidity (held in the mass_mix_ratio array read
            ! from rad_constituents) is converted to an actual mass mixing
            ! ratio.
            if (present(n_day_indices) .and. present(day_indices)) then
               do iday = 1,n_day_indices
                  vol_mix_ratio(iday,:) = mass_mix_ratio(day_indices(iday),:) / ( &
                     1._r8 - mass_mix_ratio(day_indices(iday),:) &
                  ) * mol_weight_air / mol_weight_gas(igas)
               end do
            else
               vol_mix_ratio(:,:) = mass_mix_ratio(:,:) / ( &
                  1._r8 - mass_mix_ratio(:,:) &
               ) * mol_weight_air / mol_weight_gas(igas)
            end if


         case DEFAULT

            ! Get mass mixing ratio from the rad_constituents interface
            call rad_cnst_get_gas(icall, trim(gas_species(igas)), state, pbuf, &
                                  mass_mix_ratio)

            ! Convert to volume mixing ratio by multiplying by the ratio of
            ! molecular weight of dry air to molecular weight of gas
            if (present(n_day_indices) .and. present(day_indices)) then
               do iday = 1,n_day_indices
                  vol_mix_ratio(iday,:) = mass_mix_ratio(day_indices(iday),:) &
                     * mol_weight_air / mol_weight_gas(igas)
               end do         
            else
               vol_mix_ratio(:,:) = mass_mix_ratio(:,:) &
                  * mol_weight_air / mol_weight_gas(igas)
            end if

      end select

      ! Make sure we do not have any negative volume mixing ratios
      if (any(vol_mix_ratio < 0)) then
         call endrun('set_gas_concentrations: negative gas volume mixing ratios')
      end if

      ! Populate the RRTMGP gas concentration object with values for this gas
      call handle_rrtmgp_error( &
         gas_concentrations%set_vmr(trim(gas_species(igas)),vol_mix_ratio) & 
      )
      
   end do

   ! Free memory allocated to the volume mixing ratio array
   deallocate(vol_mix_ratio)
   
end subroutine set_gas_concentrations


subroutine handle_rrtmgp_error(error_message, stop_on_error)
   character(len=*), intent(in) :: error_message
   logical, intent(in), optional :: stop_on_error
   logical :: stop_on_error_local = .true.

   ! Allow passing of an optional flag to not stop the run if an error is
   ! encountered. This allows this subroutine to be used when inquiring if a
   ! variable exists without failing.
   if (present(stop_on_error)) then
      stop_on_error_local = stop_on_error
   else
      stop_on_error_local = .true.
   end if

   ! If we encounter an error, fail if we require success. Otherwise do
   ! nothing and return silently.
   if (len(trim(error_message)) > 0) then
      if (stop_on_error_local) then
         call endrun(module_name // ' RRTMGP failure: ' // error_message)
      end if
   end if
end subroutine handle_rrtmgp_error


subroutine handle_nf90_error(error_status, stop_on_error, message)
   integer, intent(in) :: error_status
   logical, intent(in), optional :: stop_on_error
   character(len=*), intent(in), optional :: message
   logical :: stop_on_error_local = .true.

   ! Allow passing of an optional flag to not stop the run if an error is
   ! encountered. This allows this subroutine to be used when inquiring if a
   ! variable exists without failing.
   if (present(stop_on_error)) then
      stop_on_error_local = stop_on_error
   else
      stop_on_error_local = .true.
   end if

   ! If we encounter an error, fail if we require success. Otherwise do
   ! nothing and return silently.
   if (error_status /= NF90_NOERR) then
      if (stop_on_error_local) then
         if (present(message)) then
            call endrun(module_name // ' netCDF failure: ' // message)
         else
            call endrun(module_name // ' netCDF failure.')
         end if
      else
         if (present(message)) then
            print *, module_name // ' netCDF failure: ' // message
         else
            print *, module_name // ' netCDF failure.'
         end if
      end if
   end if
end subroutine handle_nf90_error


subroutine set_cloud_optics_sw(state, pbuf, cloud_sw)
   
   use physics_types, only: physics_state
   use physics_buffer, only: physics_buffer_desc
   use mo_optical_props, only: ty_optical_props_2str

   type(physics_state), intent(in) :: state
   type(physics_buffer_desc), pointer :: pbuf
   type(ty_optical_props_2str), intent(out) :: cloud_sw

   !if (use_SPCAM) then
   !   ! get cloud optics from resolved fields
   !   call get_spcam_optics(state, pbuf, cloud_sw)
   !else
   !   !call get_optics_sw(state, pbuf, cld_tau, cld_tau_w, &
   !                      cld_tau_w_g, cld_tau_w_f)

   !   ! Do MCICA sampling
   !   call do_mcica_cloud_sampling(cld_tau, cld_tau_w, cld_tau_w_g, &
   !                                cld_tau_w_f, cloud_sw)
   !end if

end subroutine set_cloud_optics_sw

subroutine output_radiation_diagnostics()
   real(r8) solin    (pcols)       ! Solar incident flux
   real(r8) fsntoa   (pcols)       ! Net solar flux at TOA
   real(r8) fsutoa   (pcols)       ! Upwelling solar flux at TOA
   real(r8) fsntoac  (pcols)       ! Clear sky net solar flux at TOA
   real(r8) fsnirt   (pcols)       ! Near-IR flux absorbed at toa
   real(r8) fsnrtc   (pcols)       ! Clear sky near-IR flux absorbed at toa
   real(r8) fsnirtsq (pcols)       ! Near-IR flux absorbed at toa >= 0.7 microns
   real(r8) fsntc    (pcols)       ! Clear sky total column abs solar flux
   real(r8) fsnsc    (pcols)       ! Clear sky surface abs solar flux
   real(r8) fsdsc    (pcols)       ! Clear sky surface downwelling solar flux
   real(r8) flut     (pcols)       ! Upward flux at top of model
   real(r8) lwcf     (pcols)       ! longwave cloud forcing
   real(r8) swcf     (pcols)       ! shortwave cloud forcing
   real(r8) flutc    (pcols)       ! Upward Clear Sky flux at top of model
   real(r8) flntc    (pcols)       ! Clear sky lw flux at model top
   real(r8) flnsc    (pcols)       ! Clear sky lw flux at srf (up-down)
   real(r8) fldsc    (pcols)       ! Clear sky lw flux at srf (down)
   real(r8) fln200   (pcols)       ! net longwave flux interpolated to 200 mb
   real(r8) fln200c  (pcols)       ! net clearsky longwave flux interpolated to 200 mb
   real(r8) fns      (pcols,pverp) ! net shortwave flux
   real(r8) fcns     (pcols,pverp) ! net clear-sky shortwave flux
   real(r8) fsn200   (pcols)       ! fns interpolated to 200 mb
   real(r8) fsn200c  (pcols)       ! fcns interpolated to 200 mb
   real(r8) fnl      (pcols,pverp) ! net longwave flux
   real(r8) fcnl     (pcols,pverp) ! net clear-sky longwave flux
   real(r8) qtot
   real(r8) factor_xy
   real(r8) trad       (pcols,pver)
   real(r8) qvrad      (pcols,pver)
   real(r8) fice       (pcols,pver)
   real(r8) aod400     (pcols)   ! AOD at 400nm at CRM grids
   real(r8) aod700     (pcols)   ! AOD at 700nm at CRM grids

   integer :: nct_tot_icld_vistau(pcols,pver) ! the number of CRM columns that has in-cloud visible sw optical depth 
   integer :: nct_liq_icld_vistau(pcols,pver) ! the number of CRM column that has liq in-cloud visible sw optical depth 
   integer :: nct_ice_icld_vistau(pcols,pver) ! the number of CRM column that has ice in-cloud visible sw optical depth 
   integer :: nct_snow_icld_vistau(pcols,pver) ! the number of CRM column that has snow in-cloud visible sw optical depth 
 
end subroutine output_radiation_diagnostics

!===============================================================================
subroutine rrtmgp_set_gases(icall, state, pbuf, nlay, gas_concs)

   ! The gases in the LW coefficients file are:
   ! H2O, CO2, O3, N2O, CO, CH4, O2, N2

   ! The memory management for the gas_concs object is internal.  The arrays passed to it
   ! are copied to the internally allocated memory.  Each call to the set_vmr method checks
   ! whether the gas already has memory allocated, and if it does that memory is deallocated
   ! and new memory is allocated.

   use physics_types, only: physics_state
   use physics_buffer, only: physics_buffer_desc
   use rad_constituents, only: rad_cnst_get_gas
   use mo_gas_concentrations, only: ty_gas_concs

   ! arguments
   integer,                     intent(in)    :: icall      ! index of climate/diagnostic radiation call
   type(physics_state), target, intent(in)    :: state
   type(physics_buffer_desc),   pointer       :: pbuf(:)
   integer,                     intent(in)    :: nlay
   type(ty_gas_concs),          intent(inout) :: gas_concs

   ! local variables
   integer :: ncol

   real(r8), pointer     :: gas_mmr(:,:)
   real(r8), allocatable :: gas_vmr(:,:)

   character(len=128)          :: errmsg
   character(len=*), parameter :: sub = 'rrtmgp_set_gases'
   !--------------------------------------------------------------------------------

   ncol = state%ncol

   ! allocate array to pass to RRTMGP gas description object
   allocate(gas_vmr(ncol,nlay))

   ! Access gas mmr from CAM data structures.  Subset and convert mmr -> vmr.
   ! If an extra layer is used copy the mixing ratios from CAM's top model layer.

   ! H20
   call rad_cnst_get_gas(icall, 'H2O', state, pbuf, gas_mmr)
   ! water vapor represented as specific humidity in CAM
   gas_vmr(:,:ktopradm) = (gas_mmr(:ncol,pver:ktopcamm:-1) / &
                       (1._r8 - gas_mmr(:ncol,pver:ktopcamm:-1))) * amdw
   if (nlay == pverp) gas_vmr(:,nlay) = gas_vmr(:,ktopradm)

   errmsg = gas_concs%set_vmr('H2O', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting H2O: '//trim(errmsg))

   ! CO2
   call rad_cnst_get_gas(icall, 'CO2', state, pbuf, gas_mmr)
   gas_vmr(:,:ktopradm) = gas_mmr(:ncol,pver:ktopcamm:-1)*amdc
   if (nlay == pverp) gas_vmr(:,nlay) = gas_vmr(:,ktopradm) 

   errmsg = gas_concs%set_vmr('CO2', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting CO2: '//trim(errmsg))

   ! O3
   call rad_cnst_get_gas(icall, 'O3',  state, pbuf, gas_mmr)
   gas_vmr(:,:ktopradm) = gas_mmr(:ncol,pver:ktopcamm:-1)*amdo
   if (nlay == pverp) gas_vmr(:,nlay) = gas_vmr(:,ktopradm) 

   errmsg = gas_concs%set_vmr('O3', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting O3: '//trim(errmsg))

   ! N2O
   call rad_cnst_get_gas(icall, 'N2O', state, pbuf, gas_mmr)
   gas_vmr(:,:ktopradm) = gas_mmr(:ncol,pver:ktopcamm:-1)*amdn
   if (nlay == pverp) gas_vmr(:,nlay) = gas_vmr(:,ktopradm) 

   errmsg = gas_concs%set_vmr('N2O', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting N2O: '//trim(errmsg))

   ! CO not available
   gas_vmr  = 1.e-7_r8

   errmsg = gas_concs%set_vmr('CO', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting CO: '//trim(errmsg))

   ! CH4
   call rad_cnst_get_gas(icall, 'CH4', state, pbuf, gas_mmr)
   gas_vmr(:,:ktopradm) = gas_mmr(:ncol,pver:ktopcamm:-1)*amdm
   if (nlay == pverp) gas_vmr(:,nlay) = gas_vmr(:,ktopradm) 

   errmsg = gas_concs%set_vmr('CH4', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting CH4: '//trim(errmsg))

   ! O2
   call rad_cnst_get_gas(icall, 'O2',  state, pbuf, gas_mmr)
   gas_vmr(:,:ktopradm) = gas_mmr(:ncol,pver:ktopcamm:-1)*amdo2
   if (nlay == pverp) gas_vmr(:,nlay) = gas_vmr(:,ktopradm) 

   errmsg = gas_concs%set_vmr('O2', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting O2: '//trim(errmsg))

   ! N2 not available
   gas_vmr  = 0.7906_r8

   errmsg = gas_concs%set_vmr('N2', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting N2'//trim(errmsg))

   ! CFCs not used?
   !call rad_cnst_get_gas(icall, 'CFC11', state, pbuf, gas_mmr)
   !call rad_cnst_get_gas(icall, 'CFC12', state, pbuf, gas_mmr)

   deallocate(gas_vmr)

end subroutine rrtmgp_set_gases


subroutine rrtmgp_set_gases_lw(icall, state, pbuf, nlay, gas_concs)

   ! The gases in the LW coefficients file are:
   ! H2O, CO2, O3, N2O, CO, CH4, O2, N2

   ! The memory management for the gas_concs object is internal.  The arrays passed to it
   ! are copied to the internally allocated memory.  Each call to the set_vmr method checks
   ! whether the gas already has memory allocated, and if it does that memory is deallocated
   ! and new memory is allocated.

   use physics_types, only: physics_state
   use physics_buffer, only: physics_buffer_desc
   use rad_constituents, only: rad_cnst_get_gas
   use mo_gas_concentrations, only: ty_gas_concs

   ! arguments
   integer,                     intent(in)    :: icall      ! index of climate/diagnostic radiation call
   type(physics_state), target, intent(in)    :: state
   type(physics_buffer_desc),   pointer       :: pbuf(:)
   integer,                     intent(in)    :: nlay
   type(ty_gas_concs),          intent(inout) :: gas_concs

   ! local variables
   integer :: ncol

   real(r8), pointer     :: gas_mmr(:,:)
   real(r8), allocatable :: gas_vmr(:,:)

   character(len=128)          :: errmsg
   character(len=*), parameter :: sub = 'rrtmgp_set_gases_lw'
   !--------------------------------------------------------------------------------

   ncol = state%ncol

   ! allocate array to pass to RRTMGP gas description object
   allocate(gas_vmr(ncol,nlay))

   ! Access gas mmr from CAM data structures.  Subset and convert mmr -> vmr.
   ! If an extra layer is used copy the mixing ratios from CAM's top model layer.

   ! H20
   call rad_cnst_get_gas(icall, 'H2O', state, pbuf, gas_mmr)
   ! water vapor represented as specific humidity in CAM
   gas_vmr(:,:ktopradm) = (gas_mmr(:ncol,pver:ktopcamm:-1) / &
                       (1._r8 - gas_mmr(:ncol,pver:ktopcamm:-1))) * amdw
   if (nlay == pverp) gas_vmr(:,nlay) = gas_vmr(:,ktopradm)

   errmsg = gas_concs%set_vmr('H2O', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting H2O: '//trim(errmsg))

   ! CO2
   call rad_cnst_get_gas(icall, 'CO2', state, pbuf, gas_mmr)
   gas_vmr(:,:ktopradm) = gas_mmr(:ncol,pver:ktopcamm:-1)*amdc
   if (nlay == pverp) gas_vmr(:,nlay) = gas_vmr(:,ktopradm) 

   errmsg = gas_concs%set_vmr('CO2', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting CO2: '//trim(errmsg))

   ! O3
   call rad_cnst_get_gas(icall, 'O3',  state, pbuf, gas_mmr)
   gas_vmr(:,:ktopradm) = gas_mmr(:ncol,pver:ktopcamm:-1)*amdo
   if (nlay == pverp) gas_vmr(:,nlay) = gas_vmr(:,ktopradm) 

   errmsg = gas_concs%set_vmr('O3', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting O3: '//trim(errmsg))

   ! N2O
   call rad_cnst_get_gas(icall, 'N2O', state, pbuf, gas_mmr)
   gas_vmr(:,:ktopradm) = gas_mmr(:ncol,pver:ktopcamm:-1)*amdn
   if (nlay == pverp) gas_vmr(:,nlay) = gas_vmr(:,ktopradm) 

   errmsg = gas_concs%set_vmr('N2O', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting N2O: '//trim(errmsg))

   ! CO not available
   gas_vmr  = 1.e-7_r8

   errmsg = gas_concs%set_vmr('CO', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting CO: '//trim(errmsg))

   ! CH4
   call rad_cnst_get_gas(icall, 'CH4', state, pbuf, gas_mmr)
   gas_vmr(:,:ktopradm) = gas_mmr(:ncol,pver:ktopcamm:-1)*amdm
   if (nlay == pverp) gas_vmr(:,nlay) = gas_vmr(:,ktopradm) 

   errmsg = gas_concs%set_vmr('CH4', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting CH4: '//trim(errmsg))

   ! O2
   call rad_cnst_get_gas(icall, 'O2',  state, pbuf, gas_mmr)
   gas_vmr(:,:ktopradm) = gas_mmr(:ncol,pver:ktopcamm:-1)*amdo2
   if (nlay == pverp) gas_vmr(:,nlay) = gas_vmr(:,ktopradm) 

   errmsg = gas_concs%set_vmr('O2', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting O2: '//trim(errmsg))

   ! N2 not available
   gas_vmr  = 0.7906_r8

   errmsg = gas_concs%set_vmr('N2', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting N2'//trim(errmsg))

   ! CFCs not used?
   !call rad_cnst_get_gas(icall, 'CFC11', state, pbuf, gas_mmr)
   !call rad_cnst_get_gas(icall, 'CFC12', state, pbuf, gas_mmr)

   deallocate(gas_vmr)

end subroutine rrtmgp_set_gases_lw

!==================================================================================================

subroutine rrtmgp_set_gases_sw( &
   icall, state, pbuf, nlay, nday, &
   idxday, gas_concs)

   ! The gases in the LW coefficients file are:
   ! H2O, CO2, O3, N2O, CO, CH4, O2, N2

   use physics_types, only: physics_state
   use physics_buffer, only: physics_buffer_desc
   use rad_constituents, only: rad_cnst_get_gas
   use mo_gas_concentrations, only: ty_gas_concs

   ! arguments
   integer,                     intent(in)    :: icall      ! index of climate/diagnostic radiation call
   type(physics_state), target, intent(in)    :: state
   type(physics_buffer_desc),   pointer       :: pbuf(:)
   integer,                     intent(in)    :: nlay
   integer,                     intent(in)    :: nday
   integer,                     intent(in)    :: idxday(:)
   type(ty_gas_concs),          intent(inout) :: gas_concs

   ! local variables
   integer :: ncol, i

   real(r8), pointer     :: gas_mmr(:,:)
   real(r8), allocatable :: gas_vmr(:,:)

   character(len=128)          :: errmsg
   character(len=*), parameter :: sub = 'rrtmgp_set_gases'
   !--------------------------------------------------------------------------------

   ncol = state%ncol

   ! allocate array to pass to RRTMGP gas description object
   allocate(gas_vmr(nday,nlay))

   ! Access gas mmr from CAM data structures.  Subset and convert mmr -> vmr.
   ! If an extra layer is used copy the mixing ratios from CAM's top model layer.

   ! H20
   call rad_cnst_get_gas(icall, 'H2O', state, pbuf, gas_mmr)
   do i = 1, nday
      ! water vapor represented as specific humidity in CAM
      gas_vmr(i,:ktopradm) = (gas_mmr(idxday(i),pver:ktopcamm:-1) / &
                             (1._r8 - gas_mmr(idxday(i),pver:ktopcamm:-1))) * amdw
   end do
   if (nlay == pverp) gas_vmr(:,nlay) = gas_vmr(:,ktopradm) 

   errmsg = gas_concs%set_vmr('H2O', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting H2O: '//trim(errmsg))

   ! CO2
   call rad_cnst_get_gas(icall, 'CO2', state, pbuf, gas_mmr)
   do i = 1, nday
      gas_vmr(i,:ktopradm) = gas_mmr(idxday(i),pver:ktopcamm:-1)*amdc
   end do
   if (nlay == pverp) gas_vmr(:,nlay) = gas_vmr(:,ktopradm) 

   errmsg = gas_concs%set_vmr('CO2', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting CO2: '//trim(errmsg))

   ! O3
   call rad_cnst_get_gas(icall, 'O3',  state, pbuf, gas_mmr)
   do i = 1, nday
      gas_vmr(i,:ktopradm) = gas_mmr(idxday(i),pver:ktopcamm:-1)*amdo

   end do
   if (nlay == pverp) gas_vmr(:,nlay) = gas_vmr(:,ktopradm) 

   errmsg = gas_concs%set_vmr('O3', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting O3: '//trim(errmsg))

   ! N2O
   call rad_cnst_get_gas(icall, 'N2O', state, pbuf, gas_mmr)
   do i = 1, nday
      gas_vmr(i,:ktopradm) = gas_mmr(idxday(i),pver:ktopcamm:-1)*amdn
   end do
   if (nlay == pverp) gas_vmr(:,nlay) = gas_vmr(:,ktopradm) 

   errmsg = gas_concs%set_vmr('N2O', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting N2O: '//trim(errmsg))

   ! CO not available
   gas_vmr = 1.e-7_r8

   errmsg = gas_concs%set_vmr('CO', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting CO: '//trim(errmsg))

   ! CH4
   call rad_cnst_get_gas(icall, 'CH4', state, pbuf, gas_mmr)
   do i = 1, nday
      gas_vmr(i,:ktopradm) = gas_mmr(idxday(i),pver:ktopcamm:-1)*amdm
   end do
   if (nlay == pverp) gas_vmr(:,nlay) = gas_vmr(:,ktopradm) 

   errmsg = gas_concs%set_vmr('CH4', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting CH4: '//trim(errmsg))

   ! O2
   call rad_cnst_get_gas(icall, 'O2',  state, pbuf, gas_mmr)
   do i = 1, nday
      gas_vmr(i,:ktopradm) = gas_mmr(idxday(i),pver:ktopcamm:-1)*amdo2
   end do
   if (nlay == pverp) gas_vmr(:,nlay) = gas_vmr(:,ktopradm) 

   errmsg = gas_concs%set_vmr('O2', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting O2: '//trim(errmsg))


   ! N2 not available
   gas_vmr = 0.7906_r8

   errmsg = gas_concs%set_vmr('N2', gas_vmr)
   if (len_trim(errmsg) > 0) call endrun(sub//': error setting N2'//trim(errmsg))

   ! CFCs not used?
   !call rad_cnst_get_gas(icall, 'CFC11', state, pbuf, gas_mmr)
   !call rad_cnst_get_gas(icall, 'CFC12', state, pbuf, gas_mmr)

   deallocate(gas_vmr)

end subroutine rrtmgp_set_gases_sw


function get_dimension_length_pio(file_handle, dimension_name, stop_on_error) result(length)

   type(file_desc_t), intent(in) :: file_handle
   character(len=*), intent(in) :: dimension_name
   logical, intent(in), optional :: stop_on_error
   integer :: length

   integer :: ierr, dimension_id
   character(len=*), parameter :: sub = 'get_dimension_length'
   character(len=128) :: error_message
   logical :: stop_on_error_local = .true.

   if (present(stop_on_error)) then
      stop_on_error_local = stop_on_error
   else
      stop_on_error_local = .true.
   end if

   ! Try to read dimension size
   ierr = pio_inq_dimid(file_handle, dimension_name, dimension_id)

   ! If unsuccessful...
   if (ierr /= PIO_NOERR) then
      error_message = sub//': '//dimension_name//' not found'
      if (stop_on_error_local) call endrun(error_message)
   else
      ierr = pio_inq_dimlen(file_handle, dimension_id, length)
   end if
end function get_dimension_length_pio


subroutine coefficients_load_from_file(this, filename)

   ! Read data from coefficients file.

   use cam_pio_utils, only: cam_pio_openfile
   use pio, only: file_desc_t, var_desc_t,               &
                  PIO_NOERR, PIO_INTERNAL_ERROR,         &
                  pio_seterrorhandling, PIO_BCAST_ERROR, &
                  pio_inq_dimlen, pio_inq_dimid,         &
                  pio_inq_varid, pio_def_var,            &
                  pio_put_var, pio_get_var,              &
                  PIO_NOWRITE, pio_closefile
   use ioFileMod, only: getfil

   ! arguments
   class(gas_optics_coefficients) :: this
   character(len=*), intent(in) :: filename

   ! local variables
   integer :: file_handle    ! file handle
   character(len=256) :: locfn ! path to actual file used

   ! File dimensions
   integer ::            &
      absorber,          &
      atmos_layer,       &
      bnd,               &
      pressure,          &
      temperature,       &
      key_absorber,    &
      pressure_interp,   &
      mixing_fraction,   &
      gpt,               &
      contributors_lower, contributors_upper, &
      minor_absorber_intervals_lower, &
      minor_absorber_intervals_upper, &
      temperature_Planck, &
      pair

   
   integer :: i, j, k
   integer :: did, vid
   integer :: ierr, varid

   character(len=128) :: error_msg

   character(len=*), parameter :: sub = 'coefs_init'
   !----------------------------------------------------------------------------

   ! Open file
   print *, 'Attempt to load file ', filename
   call getfil(filename, locfn, 0)
   !call cam_pio_openfile(file_handle, locfn, PIO_NOWRITE)
   !call pio_seterrorhandling(file_handle, PIO_BCAST_ERROR)
   call handle_nf90_error(nf90_open(locfn, NF90_NOWRITE, file_handle))

   ! Get dimensions
   absorber = get_dimension_length(file_handle, 'absorber')
   atmos_layer = get_dimension_length(file_handle, 'atmos_layer')
   bnd = get_dimension_length(file_handle, 'bnd')
   pressure = get_dimension_length(file_handle, 'pressure')
   temperature = get_dimension_length(file_handle, 'temperature')
   key_absorber = get_dimension_length(file_handle, 'key_absorber')
   pressure_interp = get_dimension_length(file_handle, 'pressure_interp')
   mixing_fraction = get_dimension_length(file_handle, 'mixing_fraction')
   gpt = get_dimension_length(file_handle, 'gpt')
   pair = get_dimension_length(file_handle, 'pair')
   minor_absorber_intervals_lower = get_dimension_length(file_handle, 'minor_absorber_intervals_lower')
   minor_absorber_intervals_upper = get_dimension_length(file_handle, 'minor_absorber_intervals_upper')
   contributors_lower = get_dimension_length(file_handle, 'contributors_lower')
   contributors_upper = get_dimension_length(file_handle, 'contributors_upper')
   temperature_Planck = get_dimension_length(file_handle, 'temperature_Planck', stop_on_error=.false.)

   ! Get variables

   ! names of absorbing gases
   allocate(this%gas_names(absorber))
   call read_field(file_handle, 'gas_names', this%gas_names)

   ! key species pair for each band
   allocate(this%key_species(2,atmos_layer,bnd))
   call read_field(file_handle, 'key_species', this%key_species)

   ! beginning and ending gpoint for each band
   allocate(this%bnd_limits_gpt(2,bnd))
   call read_field(file_handle, 'bnd_limits_gpt', this%bnd_limits_gpt)

   ! beginning and ending wavenumber for each band
   allocate(this%bnd_limits_wavenumber(2,bnd))
   call read_field(file_handle, 'bnd_limits_wavenumber', this%bnd_limits_wavenumber)

   ! pressures [hPa] for reference atmosphere; press_ref(# reference layers)
   allocate(this%press_ref(pressure))
   call read_field(file_handle, 'press_ref', this%press_ref)

   ! reference pressure separating the lower and upper atmosphere
   call read_field(file_handle, 'press_ref_trop', this%press_ref_trop)

   ! temperatures [K] for reference atmosphere; temp_ref(# reference layers)
   allocate(this%temp_ref(temperature))
   call read_field(file_handle, 'temp_ref', this%temp_ref)

   ! standard spectroscopic reference temperature [K]
   call read_field(file_handle, 'absorption_coefficient_ref_T', this%absorption_coefficient_ref_T) 

   ! standard spectroscopic reference pressure [hPa]
   call read_field(file_handle, 'absorption_coefficient_ref_P', this%absorption_coefficient_ref_P)

   ! volume mixing ratios for reference atmosphere
   allocate(this%vmr_ref(atmos_layer,key_absorber,temperature))
   call read_field(file_handle, 'vmr_ref', this%vmr_ref)

   ! absorption coefficients due to major absorbing gases
   allocate(this%kmajor(gpt,mixing_fraction,pressure_interp,temperature))
   call read_field(file_handle, 'kmajor', this%kmajor)

   ! new stuff...
   allocate(this%minor_limits_gpt_lower(pair,minor_absorber_intervals_lower))
   call read_field( &
      file_handle, 'minor_limits_gpt_lower', &
      this%minor_limits_gpt_lower &
   )

   allocate(this%minor_limits_gpt_upper(pair,minor_absorber_intervals_upper))
   call read_field( &
      file_handle, 'minor_limits_gpt_upper', &
      this%minor_limits_gpt_upper &
   )

   ! Read logicals...need to use netCDF for this
   allocate(this%minor_scales_with_density_lower(minor_absorber_intervals_lower))
   call read_field( &
      file_handle, 'minor_scales_with_density_lower', &
      this%minor_scales_with_density_lower &
   )

   allocate(this%minor_scales_with_density_upper(minor_absorber_intervals_upper))
   call read_field( &
      file_handle, 'minor_scales_with_density_upper', &
      this%minor_scales_with_density_upper &
   )

   allocate(this%scale_by_complement_lower(minor_absorber_intervals_lower))
   call read_field( &
      file_handle, 'scale_by_complement_lower', &
      this%scale_by_complement_lower &
   )

   allocate(this%scale_by_complement_upper(minor_absorber_intervals_upper))
   call read_field( &
      file_handle, 'scale_by_complement_upper', &
      this%scale_by_complement_upper &
   )

   ! absorption coefficients due to minor absorbing gases in lower part of atmosphere
   allocate(this%kminor_lower(contributors_lower,mixing_fraction,temperature))
   call read_field( &
      file_handle, 'kminor_lower', &
      this%kminor_lower &
   )

   ! absorption coefficients due to minor absorbing gases in upper part of atmosphere
   allocate(this%kminor_upper(contributors_upper,mixing_fraction,temperature))
   call read_field( &
      file_handle, 'kminor_upper', &
      this%kminor_upper &
   )

   ! new stuff...
   allocate(this%kminor_start_lower(minor_absorber_intervals_lower))
   call read_field( &
      file_handle, 'kminor_start_lower', &
      this%kminor_start_lower &
   )

   allocate(this%kminor_start_upper(minor_absorber_intervals_upper))
   call read_field( &
      file_handle, 'kminor_start_upper', &
      this%kminor_start_upper &
   )

   allocate(this%scaling_gas_lower(minor_absorber_intervals_lower))
   call read_field( &
      file_handle, 'scaling_gas_lower', &
      this%scaling_gas_lower &
   )

   allocate(this%scaling_gas_upper(minor_absorber_intervals_upper))
   call read_field( &
      file_handle, 'scaling_gas_upper', &
      this%scaling_gas_upper &
   )

   ! integrated Planck function by band
   if (var_exists(file_handle, 'totplnk')) then
      allocate(this%totplnk(temperature_Planck,bnd))
      call read_field( &
         file_handle, 'totplnk', &
         this%totplnk &
      )
   end if

   ! Planck fractions
   if (var_exists(file_handle, 'planck_fraction')) then
      allocate(this%planck_fraction(gpt,mixing_fraction,pressure_interp,temperature))
      call read_field( &
         file_handle, 'planck_fraction', &
         this%planck_fraction &
      )
   end if

   ! solar_src
   if (var_exists(file_handle, 'solar_src')) then
      allocate(this%solar_src(gpt))
      call read_field( &
         file_handle, 'solar_src', &
         this%solar_src &
      )
   end if

   ! rayleigh scattering contribution in lower part of atmosphere
   if (var_exists(file_handle, 'rayl_lower')) then
      allocate(this%rayl_lower(gpt,mixing_fraction,temperature))
      call read_field( &
         file_handle, 'rayl_lower', &
         this%rayl_lower &
      )
   end if

   ! rayleigh scattering contribution in upper part of atmosphere
   if (var_exists(file_handle, 'rayl_upper')) then
      allocate(this%rayl_upper(gpt,mixing_fraction,temperature))
      call read_field( &
         file_handle, 'rayl_upper', &
         this%rayl_upper &
      )
   end if

   ! Close file
   call handle_nf90_error(nf90_close(file_handle))

contains

      subroutine handle_pio_error(error_status, stop_on_error, message)
      integer, intent(in) :: error_status
      logical, intent(in), optional :: stop_on_error
      character(len=*), intent(in), optional :: message
      logical :: stop_on_error_local = .true.

      ! Allow passing of an optional flag to not stop the run if an error is
      ! encountered. This allows this subroutine to be used when inquiring if a
      ! variable exists without failing.
      if (present(stop_on_error)) then
         stop_on_error_local = stop_on_error
      else
         stop_on_error_local = .true.
      end if

      ! If we encounter an error, fail if we require success. Otherwise do
      ! nothing and return silently.
      if (error_status /= PIO_NOERR) then
         if (stop_on_error_local) then
            if (present(message)) then
               call endrun(module_name // ' PIO failure: ' // message)
            else
               call endrun(module_name // ' PIO failure.')
            end if
         end if
      end if
   end subroutine handle_pio_error

end subroutine coefficients_load_from_file

logical function var_exists(file_handle, var_name)

   integer, intent(in) :: file_handle
   character(len=*), intent(in) :: var_name
   integer :: varid, error_status

   ! try to find variable in file
   error_status = nf90_inq_varid(file_handle, var_name, varid)

   ! If variable ID query was successful, error_status will be PIO_NOERR
   ! and we set the existence flags to .true., otherwise we set the
   ! existence flags to .false.
   if (error_status == NF90_NOERR) then
      var_exists = .true.
   else
      var_exists = .false.
   end if

end function var_exists


function get_dimension_length(file_handle, dimension_name, stop_on_error) result(length)

   integer, intent(in) :: file_handle
   character(len=*), intent(in) :: dimension_name
   logical, intent(in), optional :: stop_on_error
   integer :: length

   integer :: ierr, dimension_id
   character(len=*), parameter :: sub = 'get_dimension_length'
   character(len=128) :: error_message
   logical :: stop_on_error_local = .true.

   if (present(stop_on_error)) then
      stop_on_error_local = stop_on_error
   else
      stop_on_error_local = .true.
   end if

   ! Try to read dimension size
   ierr = nf90_inq_dimid(file_handle, dimension_name, dimension_id)

   ! If unsuccessful...
   if (ierr == NF90_NOERR) then
      ierr = nf90_inquire_dimension(file_handle, dimension_id, len=length)
   else
      error_message = sub//': '//dimension_name//' not found'
      if (stop_on_error_local) call endrun(error_message)
   end if
end function get_dimension_length


subroutine read_integer_1d(file_handle, var_name, var_data)

   integer, intent(in) :: file_handle
   character(len=*), intent(in) :: var_name
   integer, intent(out) :: var_data(:)
   integer :: varid

   ! Get variable id
   call handle_nf90_error(nf90_inq_varid(file_handle, var_name, varid))

   ! Read variable data
   call handle_nf90_error(nf90_get_var(file_handle, varid, var_data))

end subroutine read_integer_1d

subroutine read_integer_2d(file_handle, var_name, var_data)

   integer, intent(in) :: file_handle
   character(len=*), intent(in) :: var_name
   integer, intent(out) :: var_data(:,:)
   integer :: varid

   ! Get variable id
   call handle_nf90_error(nf90_inq_varid(file_handle, var_name, varid))

   ! Read variable data
   call handle_nf90_error(nf90_get_var(file_handle, varid, var_data))

end subroutine read_integer_2d


subroutine read_integer_3d(file_handle, var_name, var_data)

   integer, intent(in) :: file_handle
   character(len=*), intent(in) :: var_name
   integer, intent(out) :: var_data(:,:,:)
   integer :: varid

   ! Get variable id
   call handle_nf90_error(nf90_inq_varid(file_handle, var_name, varid))

   ! Read variable data
   call handle_nf90_error(nf90_get_var(file_handle, varid, var_data))

end subroutine read_integer_3d


subroutine read_integer_4d(file_handle, var_name, var_data)

   integer, intent(in) :: file_handle
   character(len=*), intent(in) :: var_name
   integer, intent(out) :: var_data(:,:,:,:)
   integer :: varid

   ! Get variable id
   call handle_nf90_error(nf90_inq_varid(file_handle, var_name, varid))

   ! Read variable data
   call handle_nf90_error(nf90_get_var(file_handle, varid, var_data))

end subroutine read_integer_4d


subroutine read_real_1d(file_handle, var_name, var_data)

   integer, intent(in) :: file_handle
   character(len=*), intent(in) :: var_name
   real(r8), intent(out) :: var_data(:)
   integer :: varid

   ! Get variable id
   call handle_nf90_error(nf90_inq_varid(file_handle, var_name, varid))

   ! Read variable data
   call handle_nf90_error(nf90_get_var(file_handle, varid, var_data))

end subroutine read_real_1d

subroutine read_real_2d(file_handle, var_name, var_data)

   integer, intent(in) :: file_handle
   character(len=*), intent(in) :: var_name
   real(r8), intent(out) :: var_data(:,:)
   integer :: varid

   ! Get variable id
   call handle_nf90_error(nf90_inq_varid(file_handle, var_name, varid))

   ! Read variable data
   call handle_nf90_error(nf90_get_var(file_handle, varid, var_data))

end subroutine read_real_2d


subroutine read_real_3d(file_handle, var_name, var_data)

   integer, intent(in) :: file_handle
   character(len=*), intent(in) :: var_name
   real(r8), intent(out) :: var_data(:,:,:)
   integer :: varid

   ! Get variable id
   call handle_nf90_error(nf90_inq_varid(file_handle, var_name, varid))

   ! Read variable data
   call handle_nf90_error(nf90_get_var(file_handle, varid, var_data))

end subroutine read_real_3d


subroutine read_real_4d(file_handle, var_name, var_data)

   integer, intent(in) :: file_handle
   character(len=*), intent(in) :: var_name
   real(r8), intent(out) :: var_data(:,:,:,:)
   integer :: varid

   ! Get variable id
   call handle_nf90_error(nf90_inq_varid(file_handle, var_name, varid))

   ! Read variable data
   call handle_nf90_error(nf90_get_var(file_handle, varid, var_data))

end subroutine read_real_4d


subroutine read_logical_1d(file_handle, var_name, var_data)

   integer, intent(in) :: file_handle
   character(len=*), intent(in) :: var_name
   logical, intent(out) :: var_data(:)
   integer, allocatable :: var_data_tmp(:)
   integer :: varid

   ! NOTE:
   ! To read logicals, we need to read them as integers first, and then populate
   ! the logical array properly, assuming 0 for false and 1 for true.

   ! Get variable id
   call handle_nf90_error(nf90_inq_varid(file_handle, var_name, varid))

   ! Allocate array to hold temporary integer data
   allocate(var_data_tmp(size(var_data,1)))

   ! Read variable data to temporary integer array
   call handle_nf90_error(nf90_get_var(file_handle, varid, var_data_tmp))

   ! Populate output logical array based on values of integer temporary array
   where (var_data_tmp == 1)
      var_data = .true.
   elsewhere
      var_data = .false.
   endwhere

   ! Deallocate temporary array
   deallocate(var_data_tmp)

end subroutine read_logical_1d


subroutine read_character_1d(file_handle, var_name, var_data)

   integer, intent(in) :: file_handle
   character(len=*), intent(in) :: var_name
   character(len=*), intent(out) :: var_data(:)
   integer :: varid

   ! Get variable id
   call handle_nf90_error(nf90_inq_varid(file_handle, var_name, varid))

   ! Read variable data
   call handle_nf90_error(nf90_get_var(file_handle, varid, var_data))

end subroutine read_character_1d


subroutine read_real_scalar(file_handle, var_name, var_data)

   integer, intent(in) :: file_handle
   character(len=*), intent(in) :: var_name
   real(r8) :: var_data
   integer :: varid

   ! Get variable id
   call handle_nf90_error(nf90_inq_varid(file_handle, var_name, varid))

   ! Read variable data
   call handle_nf90_error(nf90_get_var(file_handle, varid, var_data))

end subroutine read_real_scalar


end module radiation
