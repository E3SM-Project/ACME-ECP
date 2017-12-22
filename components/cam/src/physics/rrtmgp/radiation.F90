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
   character(32), dimension(:),  allocatable :: gas_names
   integer,  dimension(:,:,:),   allocatable :: key_species
   integer,  dimension(:,:),     allocatable :: bnd_limits_gpt
   real(r8), dimension(:,:),     allocatable :: bnd_limits_wavenumber
   real(r8), dimension(:),       allocatable :: press_ref, temp_ref
   real(r8)                                  :: press_ref_trop
   real(r8)                                  :: absorption_coefficient_ref_T
   real(r8)                                  :: absorption_coefficient_ref_P
   real(r8), dimension(:,:,:),   allocatable :: vmr_ref
   real(r8), dimension(:,:,:,:), allocatable :: kmajor
   real(r8), dimension(:,:,:),   allocatable :: wv_self, wv_for
   real(r8), dimension(:,:,:,:), allocatable :: kminor_lower, kminor_upper
   real(r8), dimension(:,:),     allocatable :: totplnk
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

!===============================================================================

contains

function coefficients_free_allocated_memory(this) result(error_message)
   class(gas_optics_coefficients), intent(inout) :: this
   character(len=128) :: error_message

   ! Naively deallocate; TODO: better to check if allocated like below
   deallocate(this%gas_names, this%key_species, this%bnd_limits_gpt, &
              this%bnd_limits_wavenumber, this%press_ref, &
              this%temp_ref, &
              this%vmr_ref, this%kmajor, this%wv_self, this%wv_for, &
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

subroutine radiation_readnl(nlfile)
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

   ! Local variables
   integer :: unitn, ierr
   integer :: dtime  ! timestep size
   character(len=*), parameter :: subroutine_name = 'radiation_readnl'

   character(len=cl) :: rrtmgp_coefficients_file_lw, rrtmgp_coefficients_file_sw
   integer :: rrtmgp_iradsw, rrtmgp_iradlw, rrtmgp_irad_always
   logical :: rrtmgp_use_rad_dt_cosz, rrtmgp_spectralflux

   ! Variables defined in namelist
   ! TODO: why are these prefaced with rrtmgp instead of radiation?
   namelist /rrtmgp_nl/ rrtmgp_coefficients_file_lw, rrtmgp_coefficients_file_sw, &
                        rrtmgp_iradsw, rrtmgp_iradlw, rrtmgp_irad_always, &
                        rrtmgp_use_rad_dt_cosz, rrtmgp_spectralflux

   ! Read the namelist, only if called from master process
   ! TODO: better documentation and cleaner logic here?
   if (masterproc) then
      unitn = getunit()
      open(unitn, file=trim(nlfile), status='old')
      call find_group_name(unitn, 'rrtmgp_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, rrtmgp_nl, iostat=ierr)
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
   dtime  = get_step_size()
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
                                 use_rad_dt_cosz_out)

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

   if (present(iradsw_out)) iradsw_out = iradsw
   if (present(iradlw_out)) iradlw_out = iradlw
   if (present(iradae_out)) iradae_out = -999
   if (present(irad_always_out)) irad_always_out = irad_always
   if (present(spectralflux_out)) spectralflux_out = spectralflux
   if (present(use_rad_dt_cosz_out)) use_rad_dt_cosz_out = use_rad_dt_cosz

end subroutine radiation_defaultopts

!===============================================================================

subroutine radiation_setopts(dtime, nhtfrq, iradsw_in, iradlw_in, iradae_in, &
                             irad_always_in, spectralflux_in, &
                             use_rad_dt_cosz_in)
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
   
   ! Local namespace
   integer :: iradae   ! not used by RRTMGP

   ! Check for optional arguments and set corresponding module data if present
   if (present(iradsw_in)) iradsw = iradsw_in
   if (present(iradlw_in)) iradlw = iradlw_in
   if (present(iradae_in)) iradae = iradae_in
   if (present(irad_always_in)) irad_always = irad_always_in
   if (present(spectralflux_in)) spectralflux = spectralflux_in
   if (present(use_rad_dt_cosz_in)) use_rad_dt_cosz = use_rad_dt_cosz_in

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

   type(ty_gas_concs) :: gas_concentrations_sw, gas_concentrations_lw

   !-----------------------------------------------------------------------

   ! Initialize output fields for offline driver.
   ! TODO: do we need to keep this functionality? Where is the offline driver?
   ! Do we need to write a new offline driver for RRTMGP?
   call init_rad_data() 

   ! Initialize gas concentrations. This is needed here because the
   ! initialization routine for the kdist objects need information about what
   ! gases are available. Generally, the _set_gases_ routines here grab gas
   ! concentrations for arbitrary values of icall (the climate/diagnostic calls;
   ! i.e., various loops over the radiation for different diagnostics?), but we
   ! just set this equal to 0 here to get the gases that affect the climate.
   ! TODO: it might be better to instead implement a function that just reads
   ! the coefficients files, and stores the fields in a derived type, and then
   ! to use these routines below just in the radiation_tend subroutine for each
   ! loop over the climate/diagnostic icall just before setting up the rrtmgp
   ! objects for the radiative transfer calls. 
   !icall = 0
   !call rrtmgp_set_gases(icall, state, pbuf, pver, gas_concentrations_sw)
   !call rrtmgp_set_gases(icall, state, pbuf, pver, gas_concentrations_lw)

   ! Read gas optics coefficients from file, but do *not* initialize k_dist
   ! objects.
   call coefficients_lw%load_from_file(coefficients_file_lw)
   call coefficients_sw%load_from_file(coefficients_file_sw)

   ! Get number of bands used in shortwave and longwave and set module data
   ! appropriately so that these sizes can be used to allocate array sizes.
   !nswbands = k_dist_sw%get_nband()
   !nlwbands = k_dist_lw%get_nband()

   ! Initialize the shortwave and longwave drivers
   error_message = rrtmgp_sw_init()
   error_message = rrtmgp_lw_init()

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

   ! These are needed for calculating solar zenith angle
   use phys_grid,       only: get_rlat_all_p, get_rlon_all_p
   use time_manager,    only: get_curr_calday
   use orbit,           only: zenith

   ! RRTMGP radiation drivers
   use mo_rrtmgp_clr_all_sky, only: rrtmgp_sw=>rte_sw, rrtmgp_lw=>rte_lw

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
   integer i, k

   ! Gathered indicies of day and night columns 
   ! chunk_column_index = daytime_indices(daylight_column_index)
   integer :: n_daytime_columns     ! Number of daylight columns
   integer :: n_nighttime_columns   ! Number of night columns
   integer, dimension(pcols) :: daytime_indices   ! Indicies of daylight coumns
   integer, dimension(pcols) :: nighttime_indices ! Indicies of night coumns

   ! Set name for this subroutine for log files
   character(*), parameter :: subroutine_name = 'radiation_tend'

   !----------------------------------------------------------------------
   
   ! Set pointers to heating rates stored on physics buffer
   call pbuf_get_field(pbuf, pbuf_get_index('QRS'), qrs)
   call pbuf_get_field(pbuf, pbuf_get_index('QRL'), qrl)
  
   ! Cosine solar zenith angle for current time step
   if (swrad_off) then
      coszrs(:) = 0._r8
   else
      calday = get_curr_calday()
      call get_rlat_all_p(state%lchnk, state%ncol, clat)
      call get_rlon_all_p(state%lchnk, state%ncol, clon)
      call zenith(calday, clat, clon, coszrs, state%ncol, dt_avg)
   end if

   ! Gather night/day column indices for subsetting SW inputs; we only want to
   ! do the shortwave radiative transfer during the daytime to save
   ! computational cost (and because RRTMGP will fail for cosine solar zenith
   ! angles less than or equal to zero)
   n_daytime_columns = 0
   n_nighttime_columns = 0
   do i = 1, state%ncol
      if (coszrs(i) > 0.0_r8) then
         n_daytime_columns = n_daytime_columns + 1
         daytime_indices(n_daytime_columns) = i
      else
         n_nighttime_columns = n_nighttime_columns + 1
         nighttime_indices(n_nighttime_columns) = i
      end if
   end do

   ! Check if we are supposed to do shortwave and longwave stuff at this
   ! timestep, and if we are then we begin setting optical properties and then
   ! doing the radiative transfer separately for shortwave and longwave.
   dosw = radiation_do('sw')
   dolw = radiation_do('lw')
   if (dosw .or. dolw) then

      ! initialize RRTMGP state
      ! TODO

      ! Read k-distribution coefficients from file and initialize k-distribution
      ! objects.
      !call rrtmgp_load_coefficients(k_dist_sw, coefficients_file_sw, &
      !                              gas_concentrations_sw)
      !call rrtmgp_load_coefficients(k_dist_lw, coefficients_file_lw, &
      !                              gas_concentrations_lw)

      ! Do shortwave stuff...
      if (dosw) then
         ! Do shortwave cloud optics calculations
         call t_startf('shortwave cloud optics')
         !call set_cloud_optics_sw(state, pbuf, cloud_sw)
         call t_stopf('shortwave cloud optics')

         ! Get shortwave gas optics
         call t_startf('shortwave gas concentrations')
         !call set_gas_optics_sw(state, pbuf, gas_concentrations_sw)
         call t_stopf('shortwave gas concentrations')

         ! Get shortwave aerosol optics
         call t_startf('shortwave aerosol optics')
         !call set_aerosol_optics_sw(state, pbuf, aerosol_sw)
         call t_stopf('shortwave aerosol optics')

         ! Subset optical properties to get only daytime columns
         !call subset_daytime_optics()

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
         
      end if  ! dosw

      ! Do longwave stuff...
      if (dolw) then

         ! Do longwave cloud optics calculations
         call t_startf('longwave cloud optics')
         !call set_cloud_optics_lw(state, pbuf, cloud_lw)
         call t_stopf('longwave cloud optics')

         ! Get longwave gas optics
         call t_startf('longwave gas concentrations')
         !call set_gas_optics_lw(state, pbuf, gas_concentrations_lw)
         call t_stopf('longwave gas concentrations')

         ! Get longwave aerosol optics
         call t_startf('longwave aerosol optics')
         !call set_aerosol_optics_lw(state, pbuf, aerosol_lw)
         call t_stopf('longwave aerosol optics')

         ! Do longwave radiative transfer calculations
         call t_startf('longwave radiation calculations')
         !errmsg = rrtmgp_lw(kdist_lw, gas_concentrations_lw, &
         !                   pmid_rad, t_day, pint_rad, &
         !                   t_sfc, emis_sfc, cloud_lw, flw, flwc, &
         !                   aer_props=aerosol_lw)
         call t_stopf('longwave radiation calculations')

         ! Map RRTMGP outputs to CAM outputs
         
      end if  ! dolw

   else  ! dosw .or. dolw
      ! convert radiative heating rates from Q*dp to Q for energy conservation
      if (conserve_energy) then
         do k = 1,pver
            do i = 1,state%ncol
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
      do i = 1,state%ncol
         heating_rate(i,k) = (qrs(i,k) + qrl(i,k))/cpair * (1.e5_r8/state%pmid(i,k))**cappa
      end do
   end do
   call t_stopf('heating_rate')

   ! convert radiative heating rates to Q*dp for energy conservation
   ! TODO: this gets converted back here? what is going on here?
   if (conserve_energy) then
      do k = 1,pver
         do i = 1,state%ncol
            qrs(i,k) = qrs(i,k)*state%pdel(i,k)
            qrl(i,k) = qrl(i,k)*state%pdel(i,k)
         end do
      end do
   end if

   ! copy net sw flux to cam_out structure
   !cam_out%netsw(:state%ncol) = fsns(:state%ncol)

end subroutine radiation_tend


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


#if 0
type gas_coefficients_type

end type
subroutine read_coefficients(coefficients_file, coefficients)
   ! My own subroutine to read coefficients *without* initializing a
   ! ty_gas_optics_specification object, because such needs information about
   ! the gases that may differ with icall.
   character(len=*), intent(in) :: coefficients_file
   type(gas_coefficients_type), intent(in) :: coefficients


end subroutine
#endif


#if 0
subroutine initialize_gas_optics(coefficients_file, gas_optics)
   ! Purpose: Read data from coefficients file.

   character(len=*), intent(in) :: coefficients_file
   type(ty_gas_optics_specification), intent(out) :: gas_optics

   if (allocated(totplnk) .and. allocated(planck_frac)) then
      error_msg = kdist%init( &
         gas_names, key_species, &
         band2gpt, band_lims_wavenum, &
         press_ref, press_ref_trop, temp_ref, &
         temp_ref_p, temp_ref_t, vmr_ref, &
         kmajor, selfrefin, forrefin, kminor_lower, kminor_upper, &
         totplnk, planck_frac, rayl_lower, rayl_upper)
   else if (allocated(solar_src)) then
      error_msg = kdist%init( &
         gas_names, key_species, &
         band2gpt, band_lims_wavenum, &
         press_ref, press_ref_trop, temp_ref, &
         temp_ref_p, temp_ref_t, vmr_ref, &
         kmajor, selfrefin, forrefin, kminor_lower, kminor_upper, &
         solar_src, rayl_lower, rayl_upper)
   else
      call endrun('must supply either totplnk and planck_frac, or solar_src')
   end if

end subroutine initialize_gas_optics
#endif


#if 0
subroutine coefs_init(coefs_file, kdist)
   use cam_pio_utils, only: cam_pio_openfile
   use pio, only: file_desc_t, var_desc_t,               &
                  PIO_NOERR, PIO_INTERNAL_ERROR,         &
                  pio_seterrorhandling, PIO_BCAST_ERROR, &
                  pio_inq_dimlen, pio_inq_dimid,         &
                  pio_inq_varid, pio_def_var,            &
                  pio_put_var, pio_get_var,              &
                  PIO_NOWRITE, pio_closefile
   use ioFileMod, only: getfil


   ! Read data from coefficients file.  Initialize the kdist object.

   ! arguments
   character(len=*),                  intent(in)  :: coefs_file
   type(ty_gas_optics_specification), intent(out) :: kdist

   ! local variables
   type(file_desc_t)  :: file_handle    ! pio file handle
   character(len=256) :: locfn ! path to actual file used

   ! File dimensions
   integer ::            &
      absorber,          &
      atmos_layer,       &
      bnd,               &
      pressure,          &
      temperature,       &
      major_absorber,    &
      pressure_interp,   &
      mixing_fraction,   &
      gpt,               &
      temperature_Planck
   
   integer :: i, j, k
   integer :: did, vid
   integer :: ierr

   character(32), dimension(:),  allocatable :: gas_names
   integer,  dimension(:,:,:),   allocatable :: key_species
   integer,  dimension(:,:),     allocatable :: band2gpt
   real(r8), dimension(:,:),     allocatable :: band_lim_wavenum
   real(r8), dimension(:),       allocatable :: press_ref, temp_ref
   real(r8)                                  :: press_ref_trop, temp_ref_t, temp_ref_p
   real(r8), dimension(:,:,:),   allocatable :: vmr_ref
   real(r8), dimension(:,:,:,:), allocatable :: kmajor
   real(r8), dimension(:,:,:),   allocatable :: selfrefin, forrefin
   real(r8), dimension(:,:,:,:), allocatable :: kminor_lower, kminor_upper
   real(r8), dimension(:,:),     allocatable :: totplnk
   real(r8), dimension(:,:,:,:), allocatable :: planck_frac
   real(r8), dimension(:),       allocatable :: solar_src
   real(r8), dimension(:,:,:),   allocatable :: rayl_lower, rayl_upper

   character(len=128) :: error_msg

   character(len=*), parameter :: sub = 'coefs_init'
   !----------------------------------------------------------------------------

   ! Open file
   call getfil(coefs_file, locfn, 0)
   call cam_pio_openfile(fh, locfn, PIO_NOWRITE)

   call pio_seterrorhandling(fh, PIO_BCAST_ERROR)

   ! Get dimensions and check for consistency with parameter values

   ierr = pio_inq_dimid(fh, 'absorber', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': absorber not found')
   ierr = pio_inq_dimlen(fh, did, absorber)

   ierr = pio_inq_dimid(fh, 'atmos_layer', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': atmos_layer not found')
   ierr = pio_inq_dimlen(fh, did, atmos_layer)

   ierr = pio_inq_dimid(fh, 'bnd', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': bnd not found')
   ierr = pio_inq_dimlen(fh, did, bnd)

   ierr = pio_inq_dimid(fh, 'pressure', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': pressure not found')
   ierr = pio_inq_dimlen(fh, did, pressure)

   ierr = pio_inq_dimid(fh, 'temperature', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': temperature not found')
   ierr = pio_inq_dimlen(fh, did, temperature)

   ierr = pio_inq_dimid(fh, 'major_absorber', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': major_absorber not found')
   ierr = pio_inq_dimlen(fh, did, major_absorber)

   ierr = pio_inq_dimid(fh, 'pressure_interp', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': pressure_interp not found')
   ierr = pio_inq_dimlen(fh, did, pressure_interp)

   ierr = pio_inq_dimid(fh, 'mixing_fraction', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': mixing_fraction not found')
   ierr = pio_inq_dimlen(fh, did, mixing_fraction)

   ierr = pio_inq_dimid(fh, 'gpt', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': gpt not found')
   ierr = pio_inq_dimlen(fh, did, gpt)

   temperature_Planck = 0
   ierr = pio_inq_dimid(fh, 'temperature_Planck', did)
   if (ierr == PIO_NOERR) then
      ierr = pio_inq_dimlen(fh, did, temperature_Planck)
   end if
   
   ! Get variables

   ! names of absorbing gases
   allocate(gas_names(absorber))
   ierr = pio_inq_varid(fh, 'gas_names', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': gas_names not found')
   ierr = pio_get_var(fh, vid, gas_names)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading gas_names')

   ! key species pair for each band
   allocate(key_species(2,atmos_layer,bnd))
   ierr = pio_inq_varid(fh, 'key_species', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': key_species not found')
   ierr = pio_get_var(fh, vid, key_species)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading key_species')

   ! beginning and ending gpoint for each band
   allocate(band2gpt(2,bnd))
   ierr = pio_inq_varid(fh, 'bnd_limits_gpt', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': bnd_limits_gpt not found')
   ierr = pio_get_var(fh, vid, band2gpt)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading bnd_limits_gpt')

   ! beginning and ending wavenumber for each band
   allocate(band_lims_wavenum(2,bnd))
   ierr = pio_inq_varid(fh, 'bnd_limits_wavenumber', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': bnd_limits_wavenumber not found')
   ierr = pio_get_var(fh, vid, band_lims_wavenum)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading bnd_limits_wavenumber')

   ! pressures [hPa] for reference atmosphere; press_ref(# reference layers)
   allocate(press_ref(pressure))
   ierr = pio_inq_varid(fh, 'press_ref', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': press_ref not found')
   ierr = pio_get_var(fh, vid, press_ref)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading press_ref')

   ! reference pressure separating the lower and upper atmosphere
   ierr = pio_inq_varid(fh, 'press_ref_trop', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': press_ref_trop not found')
   ierr = pio_get_var(fh, vid, press_ref_trop)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading press_ref_trop')

   ! temperatures [K] for reference atmosphere; temp_ref(# reference layers)
   allocate(temp_ref(temperature))
   ierr = pio_inq_varid(fh, 'temp_ref', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': temp_ref not found')
   ierr = pio_get_var(fh, vid, temp_ref)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading temp_ref')

   ! standard spectroscopic reference temperature [K]
   ierr = pio_inq_varid(fh, 'absorption_coefficient_ref_T', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': absorption_coefficient_ref_T not found')
   ierr = pio_get_var(fh, vid, temp_ref_t)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading absorption_coefficient_ref_T')

   ! standard spectroscopic reference pressure [hPa]
   ierr = pio_inq_varid(fh, 'absorption_coefficient_ref_P', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': absorption_coefficient_ref_P not found')
   ierr = pio_get_var(fh, vid, temp_ref_p)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading absorption_coefficient_ref_P')

   ! volume mixing ratios for reference atmosphere
   allocate(vmr_ref(atmos_layer,major_absorber,temperature))
   ierr = pio_inq_varid(fh, 'vmr_ref', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': vmr_ref not found')
   ierr = pio_get_var(fh, vid, vmr_ref)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading vmr_ref')

   ! absorption coefficients due to major absorbing gases
   allocate(kmajor(gpt,mixing_fraction,pressure_interp,temperature))
   ierr = pio_inq_varid(fh, 'kmajor', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': kmajor not found')
   ierr = pio_get_var(fh, vid, kmajor)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading kmajor')

   ! absorption coefficients due to water vapor self continuum
   allocate(selfrefin(gpt,mixing_fraction,temperature))
   ierr = pio_inq_varid(fh, 'wv_self', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': wv_self not found')
   ierr = pio_get_var(fh, vid, selfrefin)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading wv_self')

   ! absorption coefficients due to water vapor foreign continuum
   allocate(forrefin(gpt,mixing_fraction,temperature))
   ierr = pio_inq_varid(fh, 'wv_for', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': wv_for not found')
   ierr = pio_get_var(fh, vid, forrefin)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading wv_for')

   ! absorption coefficients due to minor absorbing gases in lower part of atmosphere
   allocate(kminor_lower(absorber,gpt,mixing_fraction,temperature))
   ierr = pio_inq_varid(fh, 'kminor_lower', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': kminor_lower not found')
   ierr = pio_get_var(fh, vid, kminor_lower)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading kminor_lower')

   ! absorption coefficients due to minor absorbing gases in upper part of atmosphere
   allocate(kminor_upper(absorber,gpt,mixing_fraction,temperature))
   ierr = pio_inq_varid(fh, 'kminor_upper', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': kminor_upper not found')
   ierr = pio_get_var(fh, vid, kminor_upper)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading kminor_upper')

   ! integrated Planck function by band
   ierr = pio_inq_varid(fh, 'totplnk', vid)
   if (ierr == PIO_NOERR) then
      allocate(totplnk(temperature_Planck,bnd))
      ierr = pio_get_var(fh, vid, totplnk)
      if (ierr /= PIO_NOERR) call endrun(sub//': error reading totplnk')
   end if

   ! Planck fractions
   ierr = pio_inq_varid(fh, 'plank_fraction', vid)
   if (ierr == PIO_NOERR) then
      allocate(planck_frac(gpt,mixing_fraction,pressure_interp,temperature))
      ierr = pio_get_var(fh, vid, planck_frac)
      if (ierr /= PIO_NOERR) call endrun(sub//': error reading plank_fraction')
   end if

   ! solar_src
   ierr = pio_inq_varid(fh, 'solar_source', vid)
   if (ierr == PIO_NOERR) then
      allocate(solar_src(gpt))
      ierr = pio_get_var(fh, vid, solar_src)
      if (ierr /= PIO_NOERR) call endrun(sub//': error reading solar_source')
   end if

   ! rayleigh scattering contribution in lower part of atmosphere
   ierr = pio_inq_varid(fh, 'rayl_lower', vid)
   if (ierr == PIO_NOERR) then
      allocate(rayl_lower(gpt,mixing_fraction,temperature))
      ierr = pio_get_var(fh, vid, rayl_lower)
      if (ierr /= PIO_NOERR) call endrun(sub//': error reading rayl_lower')
   end if

   ! rayleigh scattering contribution in upper part of atmosphere
   ierr = pio_inq_varid(fh, 'rayl_upper', vid)
   if (ierr == PIO_NOERR) then
      allocate(rayl_upper(gpt,mixing_fraction,temperature))
      ierr = pio_get_var(fh, vid, rayl_upper)
      if (ierr /= PIO_NOERR) call endrun(sub//': error reading rayl_upper')
   end if

   ! Close file
   call pio_closefile(fh)

   if (allocated(totplnk) .and. allocated(planck_frac)) then
      error_msg = kdist%init( &
         gas_names, key_species, &
         band2gpt, band_lims_wavenum, &
         press_ref, press_ref_trop, temp_ref, &
         temp_ref_p, temp_ref_t, vmr_ref, &
         kmajor, selfrefin, forrefin, kminor_lower, kminor_upper, &
         totplnk, planck_frac, rayl_lower, rayl_upper)
   else if (allocated(solar_src)) then
      error_msg = kdist%init( &
         gas_names, key_species, &
         band2gpt, band_lims_wavenum, &
         press_ref, press_ref_trop, temp_ref, &
         temp_ref_p, temp_ref_t, vmr_ref, &
         kmajor, selfrefin, forrefin, kminor_lower, kminor_upper, &
         solar_src, rayl_lower, rayl_upper)
   else
      error_msg = 'must supply either totplnk and planck_frac, or solar_src'
   end if
   if (error_msg /= ' ') call endrun(sub//': ERROR: '//trim(error_msg))

   deallocate( &
      gas_names, key_species,       &
      band2gpt, band_lims_wavenum,    &
      press_ref, temp_ref, vmr_ref, &
      kmajor, selfrefin, forrefin,  &
      kminor_lower, kminor_upper)
   if (allocated(totplnk))     deallocate(totplnk)
   if (allocated(planck_frac)) deallocate(planck_frac)
   if (allocated(solar_src))   deallocate(solar_src)
   if (allocated(rayl_lower))  deallocate(rayl_lower)
   if (allocated(rayl_upper))  deallocate(rayl_upper)


end subroutine coefs_init
#endif

function get_dimension_length(file_handle, dimension_name, stop_on_error) result(length)

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
end function get_dimension_length


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
   type(file_desc_t)  :: file_handle    ! pio file handle
   character(len=256) :: locfn ! path to actual file used

   ! File dimensions
   integer ::            &
      absorber,          &
      atmos_layer,       &
      bnd,               &
      pressure,          &
      temperature,       &
      major_absorber,    &
      pressure_interp,   &
      mixing_fraction,   &
      gpt,               &
      temperature_Planck
   
   integer :: i, j, k
   integer :: did, vid
   integer :: ierr, varid

   character(len=128) :: error_msg

   character(len=*), parameter :: sub = 'coefs_init'
   !----------------------------------------------------------------------------

   ! Open file
   call getfil(filename, locfn, 0)
   call cam_pio_openfile(file_handle, locfn, PIO_NOWRITE)

   call pio_seterrorhandling(file_handle, PIO_BCAST_ERROR)

   ! Get dimensions
   absorber = get_dimension_length(file_handle, 'absorber')
   atmos_layer = get_dimension_length(file_handle, 'atmos_layer')
   bnd = get_dimension_length(file_handle, 'bnd')
   pressure = get_dimension_length(file_handle, 'pressure')
   temperature = get_dimension_length(file_handle, 'temperature')
   major_absorber = get_dimension_length(file_handle, 'major_absorber')
   pressure_interp = get_dimension_length(file_handle, 'pressure_interp')
   mixing_fraction = get_dimension_length(file_handle, 'mixing_fraction')
   gpt = get_dimension_length(file_handle, 'gpt')
   temperature_Planck = get_dimension_length(file_handle, 'temperature_Planck', stop_on_error=.false.)

   ! Get variables

   ! names of absorbing gases
   allocate(this%gas_names(absorber))
   call handle_pio_error(pio_inq_varid(file_handle, 'gas_names', varid))
   call handle_pio_error(pio_get_var(file_handle, varid, this%gas_names))

   ! key species pair for each band
   allocate(this%key_species(2,atmos_layer,bnd))
   call handle_pio_error(pio_inq_varid(file_handle, 'key_species', varid))
   call handle_pio_error(pio_get_var(file_handle, varid, this%key_species))

   ! beginning and ending gpoint for each band
   allocate(this%bnd_limits_gpt(2,bnd))
   call handle_pio_error(pio_inq_varid(file_handle, 'bnd_limits_gpt', varid))
   call handle_pio_error(pio_get_var(file_handle, varid, this%bnd_limits_gpt))

   ! beginning and ending wavenumber for each band
   allocate(this%bnd_limits_wavenumber(2,bnd))
   call handle_pio_error(pio_inq_varid(file_handle, 'bnd_limits_wavenumber', varid))
   call handle_pio_error(pio_get_var(file_handle, varid, this%bnd_limits_wavenumber))

   ! pressures [hPa] for reference atmosphere; press_ref(# reference layers)
   allocate(this%press_ref(pressure))
   call handle_pio_error(pio_inq_varid(file_handle, 'press_ref', varid))
   call handle_pio_error(pio_get_var(file_handle, varid, this%press_ref))

   ! reference pressure separating the lower and upper atmosphere
   call handle_pio_error(pio_inq_varid(file_handle, 'press_ref_trop', varid))
   call handle_pio_error(pio_get_var(file_handle, varid, this%press_ref_trop))

   ! temperatures [K] for reference atmosphere; temp_ref(# reference layers)
   allocate(this%temp_ref(temperature))
   call handle_pio_error(pio_inq_varid(file_handle, 'temp_ref', varid))
   call handle_pio_error(pio_get_var(file_handle, varid, this%temp_ref))

   ! standard spectroscopic reference temperature [K]
   call handle_pio_error(pio_inq_varid(file_handle, 'absorption_coefficient_ref_T', varid))
   call handle_pio_error(pio_get_var(file_handle, varid, this%absorption_coefficient_ref_T))

   ! standard spectroscopic reference pressure [hPa]
   call handle_pio_error(pio_inq_varid(file_handle, 'absorption_coefficient_ref_P', varid))
   call handle_pio_error(pio_get_var(file_handle, varid, this%absorption_coefficient_ref_P))

   ! volume mixing ratios for reference atmosphere
   allocate(this%vmr_ref(atmos_layer,major_absorber,temperature))
   call handle_pio_error(pio_inq_varid(file_handle, 'vmr_ref', varid))
   call handle_pio_error(pio_get_var(file_handle, varid, this%vmr_ref))

   ! absorption coefficients due to major absorbing gases
   allocate(this%kmajor(gpt,mixing_fraction,pressure_interp,temperature))
   call handle_pio_error(pio_inq_varid(file_handle, 'kmajor', varid))
   call handle_pio_error(pio_get_var(file_handle, varid, this%kmajor))

   ! absorption coefficients due to water vapor self continuum
   allocate(this%wv_self(gpt,mixing_fraction,temperature))
   call handle_pio_error(pio_inq_varid(file_handle, 'wv_self', varid))
   call handle_pio_error(pio_get_var(file_handle, varid, this%wv_self))

   ! absorption coefficients due to water vapor self continuum
   allocate(this%wv_for(gpt,mixing_fraction,temperature))
   call handle_pio_error(pio_inq_varid(file_handle, 'wv_for', varid))
   call handle_pio_error(pio_get_var(file_handle, varid, this%wv_for))

   ! absorption coefficients due to minor absorbing gases in lower part of atmosphere
   allocate(this%kminor_lower(absorber,gpt,mixing_fraction,temperature))
   call handle_pio_error(pio_inq_varid(file_handle, 'kminor_lower', varid))
   call handle_pio_error(pio_get_var(file_handle, varid, this%kminor_lower))

   ! absorption coefficients due to minor absorbing gases in upper part of atmosphere
   allocate(this%kminor_upper(absorber,gpt,mixing_fraction,temperature))
   call handle_pio_error(pio_inq_varid(file_handle, 'kminor_upper', varid))
   call handle_pio_error(pio_get_var(file_handle, varid, this%kminor_upper))

   ! integrated Planck function by band
   if (var_exists(file_handle, 'totplnk')) then
      allocate(this%totplnk(temperature_Planck,bnd))
      call handle_pio_error(pio_inq_varid(file_handle, 'totplnk', varid))
      call handle_pio_error(pio_get_var(file_handle, varid, this%totplnk))
   end if

   ! Planck fractions
   if (var_exists(file_handle, 'planck_fraction')) then
      allocate(this%planck_fraction(gpt,mixing_fraction,pressure_interp,temperature))
      call handle_pio_error(pio_inq_varid(file_handle, 'planck_fraction', varid))
      call handle_pio_error(pio_get_var(file_handle, varid, this%planck_fraction))
   end if

   ! solar_src
   if (var_exists(file_handle, 'solar_src')) then
      allocate(this%solar_src(gpt))
      call handle_pio_error(pio_inq_varid(file_handle, 'solar_src', varid))
      call handle_pio_error(pio_get_var(file_handle, varid, this%solar_src))
   end if

   ! rayleigh scattering contribution in lower part of atmosphere
   if (var_exists(file_handle, 'rayl_lower')) then
      allocate(this%rayl_lower(gpt,mixing_fraction,temperature))
      call handle_pio_error(pio_inq_varid(file_handle, 'rayl_lower', varid))
      call handle_pio_error(pio_get_var(file_handle, varid, this%rayl_lower))
   end if

   ! rayleigh scattering contribution in upper part of atmosphere
   if (var_exists(file_handle, 'rayl_upper')) then
      allocate(this%rayl_upper(gpt,mixing_fraction,temperature))
      call handle_pio_error(pio_inq_varid(file_handle, 'rayl_upper', varid))
      call handle_pio_error(pio_get_var(file_handle, varid, this%rayl_upper))
   end if

   ! Close file
   call pio_closefile(file_handle)

contains
   logical function var_exists(file_handle, var_name)

      type(file_desc_t), intent(in) :: file_handle
      character(len=*), intent(in) :: var_name
      integer :: var_id, error_status

      ! try to find variable in file
      error_status = pio_inq_varid(file_handle, var_name, var_id)

      ! If variable ID query was successful, error_status will be PIO_NOERR
      ! and we set the existence flags to .true., otherwise we set the
      ! existence flags to .false.
      if (error_status == PIO_NOERR) then
         var_exists = .true.
      else
         var_exists = .false.
      end if

   end function var_exists

   subroutine handle_pio_error(error_status, stop_on_error)
      integer, intent(in) :: error_status
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
      if (error_status /= PIO_NOERR) then
         if (stop_on_error_local) call endrun(module_name // ' PIO failure.')
      end if
   end subroutine handle_pio_error

end subroutine coefficients_load_from_file


end module radiation
