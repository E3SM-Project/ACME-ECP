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
use mo_rte_kind, only: wp

! Use my assertion routines to perform sanity checks
use assertions, only: assert, assert_valid, assert_range

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
! TODO: move these to COSP interface instead
integer, public, allocatable :: cosp_cnt(:)       ! counter for cosp
integer, public              :: cosp_cnt_init = 0 ! initial value for cosp counter

! Private module data
! TODO: remove these?
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

! Class to store fluxes on the CAM vertical grid and handle functions specific
! to CAM fluxes, like sending values to the physics buffer, exporting surface
! fluxes to the surface exchange fields, and sending outputs to the history
! buffer.
type cam_flux_type

   ! RRTMGP may operate on a subset of vertical levels in the case that the CAM
   ! vertical grid contains levels with very small or very large values of
   ! pressure that are outside of those assumed reasonable by RRTMGP (this might
   ! happen in the case of a model that extends into the stratosphere maybe?).
   ! To handle these cases, we allow RRTMGP to operate on its own vertical grid.
   ! The variables ktop and kbot are the indices to the top-most and bottom-most
   ! CAM level used in the RRTMGP grid. That is, the RRTMGP grid will consist of
   ! all CAM levels between ktop and kbot, and we will have a mapping between
   ! the CAM and RRTMGP grid where cam_level(ktop) = rrtmgp_level(1),
   ! cam_level(kbot) = rrtmgp_level(nlevels_rrtmgp).
   integer :: ktop  ! Index of top-most level fluxes were calculated at
   integer :: kbot  ! Index of bottom-most level fluxes were calculated at

   ! Fluxes mapped back to CAM vertical grid. Note that those with a level
   ! dimension are defined at the layer *interfaces*, not at the model level
   ! midpoints.
   real(r8), allocatable :: flux_dn(:,:)
   real(r8), allocatable :: flux_up(:,:)
   real(r8), allocatable :: flux_net(:,:)
   real(r8), allocatable :: bnd_flux_dn(:,:,:)
   real(r8), allocatable :: bnd_flux_up(:,:,:)
   real(r8), allocatable :: bnd_flux_net(:,:,:)
   real(r8), allocatable :: flux_net_bot(:)
   real(r8), allocatable :: flux_net_top(:)

   ! Heating rate, to be calculated from the fluxes at model level interfaces.
   real(r8), allocatable :: heating_rate(:,:)

   ! Variable to hold the "type" of flux a particlar instance of this object
   ! contains, being either "shortwave" or "longwave". This will be used to
   ! decide how to calculate net fluxes (because net is defined as down minus up
   ! for shortwave but up minus down for longwave) and how to map fluxes in this
   ! object to outputs.
   character(len=32) :: flux_type

contains

   ! Type-bound procedures for this class.
   procedure :: initialize=>cam_fluxes_initialize
   procedure :: finalize=>cam_fluxes_finalize
   procedure :: set_fluxes=>cam_fluxes_set_fluxes
   procedure :: calculate_heating_rate=>cam_fluxes_calculate_heating_rate
   procedure :: to_pbuf=>cam_fluxes_to_pbuf
   procedure :: export_surface_fluxes=>cam_fluxes_export_surface_fluxes

end type cam_flux_type

! give this module a name
character(len=*), parameter :: module_name = 'radiation'

interface clip_values
   module procedure clip_values_1d, clip_values_2d
end interface clip_values

!===============================================================================

contains

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

   ! Variables defined in namelist
   namelist /radiation_nl/ rrtmgp_coefficients_file_lw, &
                           rrtmgp_coefficients_file_sw, &
                           iradsw, iradlw, irad_always, &
                           use_rad_dt_cosz, spectralflux

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

#ifdef SPMD
   ! Broadcast namelist variables
   ! TODO: do we need to broadcast? Only if masterproc?
   call mpibcast(rrtmgp_coefficients_file_lw, cl, mpi_character, mstrid, mpicom, ierr)
   call mpibcast(rrtmgp_coefficients_file_sw, cl, mpi_character, mstrid, mpicom, ierr)
   call mpibcast(iradsw, 1, mpi_integer, mstrid, mpicom, ierr)
   call mpibcast(iradlw, 1, mpi_integer, mstrid, mpicom, ierr)
   call mpibcast(irad_always, 1, mpi_integer, mstrid, mpicom, ierr)
   call mpibcast(use_rad_dt_cosz, 1, mpi_logical, mstrid, mpicom, ierr)
   call mpibcast(spectralflux, 1, mpi_logical, mstrid, mpicom, ierr)
#endif

   ! Set module data
   coefficients_file_lw = rrtmgp_coefficients_file_lw
   coefficients_file_sw = rrtmgp_coefficients_file_sw

   ! Convert iradsw, iradlw and irad_always from hours to timesteps if necessary
   if (present(dtime_in)) then
      dtime = dtime_in
   else
      dtime = get_step_size()
   end if
   if (iradsw      < 0) iradsw      = nint((-iradsw     *3600._r8)/dtime)
   if (iradlw      < 0) iradlw      = nint((-iradlw     *3600._r8)/dtime)
   if (irad_always < 0) irad_always = nint((-irad_always*3600._r8)/dtime)

   ! Print runtime options to log.
   if (masterproc) then
      write(iulog,*) 'RRTMGP radiation scheme parameters:'
      write(iulog,10) trim(coefficients_file_lw), trim(coefficients_file_sw), &
                      iradsw, iradlw, irad_always, &
                      use_rad_dt_cosz, spectralflux
   end if
10 format('  LW coefficents file: ',                                a/, &
          '  SW coefficents file: ',                                a/, &
          '  Frequency (timesteps) of Shortwave Radiation calc:  ',i5/, &
          '  Frequency (timesteps) of Longwave Radiation calc:   ',i5/, &
          '  SW/LW calc done every timestep for first N steps. N=',i5/, &
          '  Use average zenith angle:                           ',l5/, &
          '  Output spectrally resolved fluxes:                  ',l5/)

end subroutine radiation_readnl

!================================================================================================

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
   use physconst,          only: gravit, epsilo, stebol, &
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
   use mo_rrtmgp_clr_all_sky, only: rrtmgp_sw_init=>rte_sw_init, &
                                    rrtmgp_lw_init=>rte_lw_init
   use mo_load_coefficients, only: rrtmgp_load_coefficients=>load_and_init
   use mo_gas_concentrations, only: ty_gas_concs

   ! For optics
   use cloud_rad_props, only: cloud_rad_props_init

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

   ! Initialize cloud optics
   call cloud_rad_props_init()

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

   ! TODO: what is this?
   call hirsbt_init()

   ! "irad_always" is number of time steps to execute radiation continuously from start of
   ! initial OR restart run
   nstep = get_nstep()
   if (irad_always > 0) then
      irad_always = irad_always + nstep
   end if

   ! Initialize the satellite simulator package (COSP). 
   ! TODO: Should this be moved to a higher level? 
   ! This should probably not be specific to a given radiation
   ! package. Also move most of this to cospsimulator package to handle itself,
   ! rather than relying on radiation driver handling this logic. Too much
   ! duplicate code.
   !if (docosp) call cospsimulator_intr_init()
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

   !
   ! Add fields to history buffer
   !

   ! Shortwave radiation
   call addfld('TOT_CLD_VISTAU', (/ 'lev' /), 'A',   '1', &
               'Total gridbox cloud visible optical depth', &
               sampling_seq='rad_lwsw', flag_xyfill=.true.)
   call addfld('TOT_ICLD_VISTAU', (/ 'lev' /), 'A',  '1', &
               'Total in-cloud visible optical depth', &
               sampling_seq='rad_lwsw', flag_xyfill=.true.)
   call addfld('LIQ_ICLD_VISTAU', (/ 'lev' /), 'A',  '1', &
               'Liquid in-cloud visible optical depth', &
               sampling_seq='rad_lwsw', flag_xyfill=.true.)
   call addfld('ICE_ICLD_VISTAU', (/ 'lev' /), 'A',  '1', &
               'Ice in-cloud visible optical depth', &
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
                     'Solar downward near infrared direct to surface', &
                     sampling_seq='rad_lwsw')
         call addfld('SOLS'//diag(icall),  horiz_only, 'A',   'W/m2', &
                     'Solar downward visible direct to surface', &
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
         call addfld('FNS'//diag(icall),  (/ 'ilev' /),  'I',  'W/m2', &
                     'Shortwave net flux')
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

   ! Cosine of solar zenith angle (primarily for debugging)
   call addfld('COSZRS', horiz_only, 'A', 'None', &
               'Cosine of solar zenith angle', &
               sampling_seq='rad_lwsw')

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
         call addfld('FNL'//diag(icall),     (/'ilev'/), 'I', 'W/m2', &
                     'Longwave net flux')
         call addfld('FULC'//diag(icall),    (/'ilev'/), 'I', 'W/m2', &
                     'Longwave clear-sky upward flux')
         call addfld('FDLC'//diag(icall),    (/'ilev'/), 'I', 'W/m2', &
                     'Longwave clear-sky downward flux')

         call add_default('QRL'//diag(icall),   1, ' ')
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
   ! TODO: why include this here? Can't this be calculated from QRS and QRL,
   ! which are already output?
   call addfld ('HR',(/ 'lev' /), 'A','K/s','Heating rate needed for d(theta)/dt computation')

   if (history_budget .and. history_budget_histfile_num > 1) then
      call add_default ('QRL     ', history_budget_histfile_num, ' ')
      call add_default ('QRS     ', history_budget_histfile_num, ' ')
   end if

   if (history_vdiag) then
      call add_default('FLUT', 2, ' ')
      call add_default('FLUT', 3, ' ')
   end if

   cldfsnow_idx = pbuf_get_index('CLDFSNOW',errcode=err)
   if (cldfsnow_idx > 0) then
      call addfld('CLDFSNOW',(/ 'lev' /),'I','1','CLDFSNOW',flag_xyfill=.true.)
      call addfld('SNOW_ICLD_VISTAU', (/ 'lev' /), 'A', '1', &
                  'Snow in-cloud extinction visible sw optical depth', &
                  sampling_seq='rad_lwsw', flag_xyfill=.true.)
   endif

   ! Fluxes for debugging
   call addfld('FLUX_SW_UP', (/'lev'/), 'I', 'W/m^2', &
               'Upwelling shortwave flux')
   call addfld('FLUX_SW_DN', (/'lev'/), 'I', 'W/m^2', &
               'Downwelling shortwave flux')

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

   ! For calculating radiative heating tendencies
   use radheat,         only: radheat_tend

   ! For getting radiative constituent gases
   use rad_constituents, only: N_DIAG, rad_cnst_get_call_list

   ! RRTMGP radiation drivers and derived types, needed for interacting with the
   ! RRTMGP radiative transfer interface
   use mo_rrtmgp_clr_all_sky, only: rrtmgp_sw=>rte_sw, rrtmgp_lw=>rte_lw
   use mo_gas_concentrations, only: ty_gas_concs
   use mo_optical_props, only: ty_optical_props, &
                               ty_optical_props_1scl, &
                               ty_optical_props_2str
   use mo_fluxes_byband, only:  ty_fluxes_byband

   ! CAM history module provides subroutine to send output data to the history
   ! buffer to be aggregated and written to disk
   use cam_history, only: outfld

   ! DEBUG for COSZRS...this should be moved to subroutine, but copying
   ! implementation from RRTMG to figure out why this is producing different
   ! results right now.
   use phys_grid,    only: get_rlat_all_p, get_rlon_all_p
   use time_manager, only: get_curr_calday, is_first_step
   use orbit,        only: zenith

   ! ---------------------------------------------------------------------------
   ! Arguments
   ! ---------------------------------------------------------------------------

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
   ! ---------------------------------------------------------------------------
   logical :: dosw, dolw
   integer :: nstep  ! current timestep number

   real(r8) :: heating_rate(pcols,pver)  ! Net heating rate

   ! Pointers to heating rates on physics buffer
   real(r8), pointer :: qrs(:,:) => null()  ! shortwave radiative heating rate 
   real(r8), pointer :: qrl(:,:) => null()  ! longwave  radiative heating rate 

   ! Flag to carry (QRS,QRL)*dp across time steps. 
   ! TODO: what does this mean?
   logical :: conserve_energy = .true.

   ! Loop variables; i is used to loop over columns, and k is used to loop over
   ! levels
   integer :: i, k
   integer :: icall

   ! nlay is number of layers used in radiation calculation, nlev is number of
   ! interface levels between those layers
   integer :: ncol, nlay, nlev

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

   ! Cosine solar zenith angle for all columns in chunk
   !real(r8), allocatable :: coszrs(:)
   real(r8) :: coszrs(pcols)

   ! State fields that are passed into RRTMGP. Some of these may need to
   ! modified from what exist in the physics_state object, i.e. to clip
   ! temperatures to make sure they are within the valid range.
   real(r8), allocatable :: coszrs_rad(:)  ! cosine solar zenith angle
   real(r8), allocatable :: ts_rad(:)  ! Surface temperature
   real(r8), allocatable :: t_rad(:,:)  ! Temperature
   real(r8), allocatable :: pmid_rad(:,:)  ! Pressure at level midpoints
   real(r8), allocatable :: pint_rad(:,:)  ! Pressure at level interfaces
   real(r8), allocatable :: surface_emissivity(:,:)

   ! Albedo for shortwave calculations
   real(r8), allocatable :: alb_dir(:,:)
   real(r8), allocatable :: alb_dif(:,:)
   
   ! RRTMGP types
   type(ty_gas_concs) :: gas_concentrations_sw, gas_concentrations_lw
   type(ty_optical_props_1scl) :: aerosol_optics_lw
   type(ty_optical_props_1scl) :: cloud_optics_lw
   type(ty_optical_props_2str) :: aerosol_optics_sw
   type(ty_optical_props_2str) :: cloud_optics_sw
   type(ty_fluxes_byband) :: flux_lw_allsky, flux_lw_clearsky
   type(ty_fluxes_byband) :: flux_sw_allsky, flux_sw_clearsky

   ! Net fluxes needed for radheat_tend
   type(cam_flux_type) :: cam_fluxes_sw, cam_fluxes_lw
   type(cam_flux_type) :: cam_fluxes_sw_clearsky, cam_fluxes_lw_clearsky

   ! DEBUG for COSZRS
   real(r8) :: calday                        ! current calendar day
   real(r8) :: clat(pcols)                   ! current latitudes(radians)
   real(r8) :: clon(pcols)                   ! current longitudes(radians)
   !----------------------------------------------------------------------
   
   ! Number of physics columns in this "chunk"; used in multiple places
   ! throughout this subroutine, so set once for convenience
   ncol = state%ncol

   ! Set pointers to heating rates stored on physics buffer. These will be
   ! modified in this routine.
   call pbuf_get_field(pbuf, pbuf_get_index('QRS'), qrs)
   call pbuf_get_field(pbuf, pbuf_get_index('QRL'), qrl)
  
   ! Figure out which vertical levels to use; only use those between the min and
   ! max of press_ref read from the coefficients file. I think the implementation 
   ! that Brian Eaton at NCAR used was to look at the reference pressure pref. We
   ! could also look at the min and max over all columns in this chunk. That
   ! might be the most transparent thing to do. Note that this code here finds
   ! the indices to the top and bottom interfaces, not model midpoints. Note
   ! also that in the CESM implementation (and also in the RRTMG implementation)
   ! an extra layer is added above the RRTMGP topmost layer. I am not 100% sure
   ! why this is, but I think maybe it is to get "top of atmosphere" vs "top of
   ! model" heating and fluxes? TODO.
   ktop = min_pressure_index(state%pint(:ncol,:pverp), k_dist_lw%get_press_ref_min())
   kbot = max_pressure_index(state%pint(:ncol,:pverp), k_dist_lw%get_press_ref_max())
   
   ! Size of grid for radiation; note kbot > ktop because CAM grid goes from TOA
   ! to surface. RRTMGP grid will go from surface to TOA internally, but this is
   ! handled automatically so we do not worry about this here. Note that nlay is
   ! the number of *layers* used, so the number of *interfaces* should be nlay+1
   nlay = kbot - ktop

   if (ktop < 1) call endrun('ktop < 1')
   if (kbot > pverp) call endrun('kbot > pverp')

   ! Sanity check on nlay (make sure we found some valid pressure levels)
   if (nlay < 1) call endrun(module_name // ': no valid pressure levels?')
   if (nlay > pver) call endrun(module_name // ': too many pressure levels?')

   ! Initialize object to hold output fluxes on the CAM vertical grid, which may
   ! differ from the radiation grid (if CAM grid contains pressures below min
   ! radiation pressure). Note that we are initializing the fluxes to have
   ! vertical dimension pverp, because the fluxes are defined on level
   ! interfaces, not on layers.
   call cam_fluxes_sw%initialize(ncol, pverp, nswbands, ktop, kbot, 'shortwave')
   call cam_fluxes_lw%initialize(ncol, pverp, nlwbands, ktop, kbot, 'longwave')
   call cam_fluxes_sw_clearsky%initialize(ncol, pverp, nswbands, ktop, kbot, 'shortwave')
   call cam_fluxes_lw_clearsky%initialize(ncol, pverp, nlwbands, ktop, kbot, 'longwave')

   ! Get cosine solar zenith angle for current time step. If the swrad_off flag
   ! is set, meaning we should not do SW radiation, then we just set coszrs to
   ! zero everywhere. TODO: why not just set dosw false and skip this
   ! loop?
   !allocate(coszrs(ncol))
   !allocate(coszrs(pcols))
   !call set_cosine_solar_zenith_angle(state, dt_avg, coszrs)
   calday = get_curr_calday()
   call get_rlat_all_p(state%lchnk, ncol, clat)
   call get_rlon_all_p(state%lchnk, ncol, clon)
   call zenith(calday, clat, clon, coszrs, ncol, dt_avg)

   ! Send to history buffer; NOTE: I'm not sure I can do this like
   ! this...outfld in cam_history says that we should use "pcols" for
   ! "phys_decomp"; is this phys_decomp? What does that mean? It would be
   ! great if we did not have to allocate arrays of size pcols just to
   ! only populate ncol of them for output to the history buffer, but I do
   ! no know if this is how outfld works.
   call outfld('COSZRS', coszrs(:ncol), ncol, state%lchnk)

   if (swrad_off) then
      coszrs(:) = 0._r8
   endif

   ! Number of daytime columns in curent chunk.
   nday = count(coszrs(:ncol) > 0)
   allocate(day_indices(nday))

   ! Gather night/day column indices for subsetting SW inputs; we only want to
   ! do the shortwave radiative transfer during the daytime to save
   ! computational cost (and because RRTMGP will fail for cosine solar zenith
   ! angles less than or equal to zero)
   call set_daytime_indices(coszrs(:ncol), day_indices)

   ! Check daytime_indices
   call assert_range(day_indices, 1, ncol, 'radiation_tend: day_indices')
   
   ! Check if we are supposed to do shortwave and longwave stuff at this
   ! timestep, and if we are then we begin setting optical properties and then
   ! doing the radiative transfer separately for shortwave and longwave.
   dosw = radiation_do('sw')
   dolw = radiation_do('lw')
   if (dosw .or. dolw) then

      ! Do shortwave stuff...
      if (dosw) then

         ! DEBUG: Make sure temperatures are within range. This should *not* be
         ! necessary, but the addition of the shortwave cloud optics seems to be
         ! creating shortwave heating rates that are leading to unphysical
         ! temperatures, so we put this code here to fail when that happens.
         call assert_range(state%t(:ncol,:), &
                           k_dist_sw%get_temp_ref_min(), &
                           k_dist_sw%get_temp_ref_max(), 'state%t')

         ! Allocate RRTMGP input variables, only as large as (nday,nlay). We do
         ! this because radiation may be calculated on a reduced vertical grid
         ! (nlay <= pver), and because we only want to pass the daytime columns
         ! to RRTMGP (because RRTMGP does not skip nighttime columns internally
         ! and rather just fails with a fatal error if the cosine solar zenith
         ! angle is less than zero for *any* of the columns).
         allocate(coszrs_rad(nday), ts_rad(nday), t_rad(nday,nlay), &
                  pmid_rad(nday,nlay), pint_rad(nday,nlay+1), &
                  alb_dir(nswbands,nday), alb_dif(nswbands,nday))

         ! Populate RRTMGP input variables. Use the day_indices index array to
         ! map CAM variables on all columns to the daytime-only arrays, and take
         ! only the ktop:kbot vertical levels (mapping CAM vertical grid to
         ! RRTMGP vertical grid).
         do i = 1,nday
            if (day_indices(i) < 1 .or. day_indices(i) > ncol) then
               print *, 'iday = ', i, '; day_indices(i) = ', day_indices(i)
               call endrun('Day_indices are screwy...crashing now.')
            end if
            coszrs_rad(i) = coszrs(day_indices(i))
            ts_rad(i) = cam_in%ts(day_indices(i))
            t_rad(i,:nlay) = state%t(day_indices(i),ktop:kbot-1)
            pmid_rad(i,:nlay) = state%pmid(day_indices(i),ktop:kbot-1)
            pint_rad(i,:nlay+1) = state%pint(day_indices(i),ktop:kbot)
         end do

         ! DEBUG: Make sure temperatures are within range. This should *not* be
         ! necessary, but the addition of the shortwave cloud optics seems to be
         ! creating shortwave heating rates that are leading to unphysical
         ! temperatures, so we put this code here to fail when that happens.
         call assert_range(t_rad, &
                           k_dist_sw%get_temp_ref_min(), &
                           k_dist_sw%get_temp_ref_max(), &
                           't_rad')

         ! Make sure temperatures are within range and clip if they are not.
         ! NOTE: We should *not* be doing this here. Somehow, we encountered
         ! some very low (~150 K) temperatures (in lower levels?) in ne4
         ! simulations, and these are below the allowed range in RRTMGP. We
         ! really need to figure out what is going on here, but for now we
         ! simply clip the temperatures if they are below the limits expected by
         ! RRTMGP to get the code to run. We do however issue a warning when
         ! this happens (warn=.true.).
         call clip_values( &
            t_rad, &
            k_dist_sw%get_temp_ref_min(), &
            k_dist_sw%get_temp_ref_max(), &
            varname='sw t_rad', warn=.true. &
         )
         call clip_values( &
            ts_rad, &
            k_dist_sw%get_temp_ref_min(), &
            k_dist_sw%get_temp_ref_max(), &
            varname='sw ts_rad', warn=.true. &
         )

         ! Get albedo. This uses CAM routines internally and just provides a
         ! nice wrapper to improve readability of the code here.
         call get_albedo(nday, day_indices(:nday), cam_in, alb_dir, alb_dif)

         ! Allocate shortwave outputs; why is this not part of the
         ! ty_fluxes_byband object?
         ! NOTE: fluxes defined at interfaces, so initialize to have vertical
         ! dimension nlay+1, while we initialized the RRTMGP input variables to
         ! have vertical dimension nlay.
         call initialize_rrtmgp_fluxes( &
            nday, nlay+1, nswbands, &
            flux_sw_allsky &
         )
         call initialize_rrtmgp_fluxes( &
            nday, nlay+1, nswbands, &
            flux_sw_clearsky &
         )
 
         ! Initialize cloud optics object. Cloud optics are defined at level
         ! midpoints (to represent optical properties through each layer), so
         ! these are allocated with vertical dimension nlay rather than nlay+1.
         ! TODO: should this go in set_cloud_optics_sw?
         call handle_rrtmgp_error(cloud_optics_sw%init_2str( &
            nday, nlay, k_dist_sw%get_ngpt(), &
            name='shortwave cloud optics' &
         ))

         ! Do shortwave cloud optics calculations. Note that we have to pass
         ! ktop and kbot here. Maybe a better thing to do would be to set cloud
         ! optics on the CAM grid first, and then map to the RRTMGP grid here.
         ! That might improve readability and reduce the potential for mistakes?
         call t_startf('shortwave cloud optics')
         call set_cloud_optics_sw( &
            day_indices(:nday), ktop, kbot, k_dist_sw, state, pbuf, &
            cloud_optics_sw &
         )
         call t_stopf('shortwave cloud optics')

         ! Send cloud optics to history buffer
         !call output_cloud_optics_sw(state, pbuf, cloud_optics_sw)

!        ! Get shortwave aerosol optics
!        call t_startf('shortwave aerosol optics')
!        !call set_aerosol_optics_sw(state, pbuf, aerosol_sw)
!        call t_stopf('shortwave aerosol optics')

         ! Loop over diagnostic calls (TODO: more documentation on what this
         ! means)

         ! get list of active radiation calls
         call rad_cnst_get_call_list(active_calls)

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

               ! DEBUG
               cloud_optics_sw%tau = 0
               cloud_optics_sw%ssa = 1
               cloud_optics_sw%g = 0

               ! Do shortwave radiative transfer calculations
               call t_startf('shortwave radiation calculations')
               call handle_rrtmgp_error(rrtmgp_sw( &
                  k_dist_sw, gas_concentrations_sw, &
                  pmid_rad, t_rad, pint_rad, &
                  coszrs_rad, alb_dir, alb_dif, cloud_optics_sw, &
                  flux_sw_allsky, flux_sw_clearsky &
                  !aer_props=aerosol_sw)
               ))
               call t_stopf('shortwave radiation calculations')

               ! Map RRTMGP fluxes to CAM vertical grid and calculate derived fluxes
               ! NOTE: we need to pass day_indices to properly map columns from
               ! the RRTMGP outputs (only nday long to only run over sunlit
               ! columns) to the CAM outputs (ncol-sized to make the rest of the
               ! routines easier to deal with)
               call cam_fluxes_sw%set_fluxes(flux_sw_allsky, cam_indices=day_indices(:nday))
               call cam_fluxes_sw_clearsky%set_fluxes(flux_sw_clearsky, cam_indices=day_indices(:nday))

               ! Calculate shortwave heating rate
               call cam_fluxes_sw%calculate_heating_rate(state%pdel(:ncol,:))
               call cam_fluxes_sw_clearsky%calculate_heating_rate(state%pdel(:ncol,:))

               ! NOTE: when use_spcam == .true. and if we are looping over
               ! column, then add a routine here to aggregate/average fluxes
               ! over all CRM columns:
               !
               !call cam_fluxes_sw_all%aggregate(cam_fluxes_sw)
               !
               ! Other idea: in the McICA sampling, just randomly select CRM
               ! columns for each g-point, and calculate domain mean fluxes as
               ! done here. This would be much faster, but less accurate..

               ! Copy fluxes to pbuf
               call cam_fluxes_sw%to_pbuf(pbuf)

               ! Write outputs to history tapes
               call radiation_output_sw( &
                  icall, state, cam_fluxes_sw, cam_fluxes_sw_clearsky &
               )
            end if
         end do

         ! Clean up
         !deallocate(coszrs, coszrs_rad, ts_rad, t_rad, &
         !           pmid_rad, pint_rad, day_indices, &
         !           alb_dir, alb_dif)
         deallocate(coszrs_rad, ts_rad, t_rad, &
                    pmid_rad, pint_rad, day_indices, &
                    alb_dir, alb_dif)

         qrs(:ncol,:pver) = cam_fluxes_sw%heating_rate(:ncol,:pver)

         call assert_valid(qrs, 'qrs')
      end if  ! dosw

      ! Do longwave stuff...
      if (dolw) then

         ! initialize aerosol optics object
         ! NOTE: we can initialize the aerosol optical properties by ngpt or by
         ! nlwbands. We initialize by nlwbands here, because that is what the
         ! internal CAM aer_rad_props routines expect. There is logic
         ! encapsulated in RRTMGP to figure out how the optical properties are
         ! split up.
         ! TODO: should this have size nlay or nlay+1?
         call handle_rrtmgp_error(aerosol_optics_lw%init_1scl( &
            ncol, nlay, nlwbands, &
            name='longwave aerosol optics' &
         ))

         ! Initialize cloud optics object
         ! TODO: should this have size nlay or nlay+1?
         call handle_rrtmgp_error(cloud_optics_lw%init_1scl( &
            ncol, nlay, k_dist_lw%get_ngpt(), &
            name='longwave cloud optics' &
         ))

         ! Manually set cloud optical depth to zero since we do not have a
         ! function to set appropriately yet. TODO: move this to the
         ! "set_cloud_optics_lw" routine.
         cloud_optics_lw%tau(:,:,:) = 0.0

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
         allocate(ts_rad(ncol), t_rad(ncol,nlay), pmid_rad(ncol,nlay), pint_rad(ncol,nlay+1))

         ! Set surface emissivity to 1 here. There is a note in the RRTMG
         ! implementation that this is treated in the land model, but the old
         ! RRTMG implementation also sets this to 1. This probably does not make
         ! a lot of difference either way, but if a more intelligent value
         ! exists or is assumed in the model we should use it here as well.
         ! TODO: set this more intelligently?
         allocate(surface_emissivity(nlwbands,ncol))
         surface_emissivity(:,:) = 1.0_r8

         ! get list of active radiation calls
         call rad_cnst_get_call_list(active_calls)

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
               ts_rad(:ncol) = cam_in%ts(:ncol)
               t_rad(:ncol,:nlay) = state%t(:ncol,ktop:kbot-1)
               pmid_rad(:ncol,:nlay) = state%pmid(:ncol,ktop:kbot-1)
               pint_rad(:ncol,:nlay+1) = state%pint(:ncol,ktop:kbot)

               ! Clip temperature values in case they are out of range
               call clip_values( &
                  t_rad, &
                  k_dist_lw%get_temp_ref_min(), &
                  k_dist_lw%get_temp_ref_max() &
               )
               call clip_values( &
                  ts_rad, &
                  k_dist_lw%get_temp_ref_min(), &
                  k_dist_lw%get_temp_ref_max() &
               )

               ! Do longwave radiative transfer calculations
               call t_startf('longwave radiation calculations')
               call handle_rrtmgp_error(rrtmgp_lw( &
                  k_dist_lw, gas_concentrations_lw, &
                  pmid_rad(:ncol,:nlay), t_rad(:ncol,:nlay), &
                  pint_rad(:ncol,:nlay+1), ts_rad(:ncol), &
                  surface_emissivity(:nlwbands,:ncol), &
                  cloud_optics_lw, &
                  flux_lw_allsky, flux_lw_clearsky  &
                  !aer_props=aerosol_optics_lw &
               ))
               call t_stopf('longwave radiation calculations')

               ! Map RRTMGP fluxes to CAM vertical grid and calculate derived fluxes
               call cam_fluxes_lw%set_fluxes(flux_lw_allsky)
               call cam_fluxes_lw_clearsky%set_fluxes(flux_lw_clearsky)

               ! Calculate longwave heating rate
               call cam_fluxes_lw%calculate_heating_rate(state%pdel(:ncol,:))
               call cam_fluxes_lw_clearsky%calculate_heating_rate(state%pdel(:ncol,:))

               ! Copy fluxes to pbuf
               call cam_fluxes_lw%to_pbuf(pbuf)

               ! Write outputs to history tapes
               call radiation_output_lw( &
                  icall, state, cam_fluxes_lw, cam_fluxes_lw_clearsky &
               )

            end if  ! active calls
         end do  ! loop over diagnostic calls
               
         ! Free memory allocated to RRTMGP input state variables
         deallocate(ts_rad, t_rad, pmid_rad, pint_rad, surface_emissivity)

         ! copy heating rate to qrl pointer
         qrl(:ncol,:pver) = cam_fluxes_lw%heating_rate(:ncol,:pver)

         call assert_valid(qrl, 'qrl (after assigment)')
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
   call radheat_tend(state, pbuf, ptend, &
                     qrl, qrs, &
                     cam_fluxes_sw%flux_net_bot, cam_fluxes_sw%flux_net_top, &
                     cam_fluxes_lw%flux_net_bot, cam_fluxes_lw%flux_net_top, &
                     cam_in%asdir, net_flx)
   call t_stopf('radheat_tend')

   ! DEBUG: do not apply radiative heating tendency
   ptend%ls = .false.
   ptend%s(:ncol,:) = 0._r8

   ! DEBUG: set net flux to zero
   net_flx = 0._r8

   ! Compute net heating rate for dtheta/dt
   ! TODO: how is this different than above?
   !call t_startf('heating_rate')
   !do k=1,pver
   !   do i = 1,ncol
   !      heating_rate(i,k) = (qrs(i,k) + qrl(i,k))/cpair * (1.e5_r8/state%pmid(i,k))**cappa
   !   end do
   !end do
   !call t_stopf('heating_rate')

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

   ! copy surface sw fluxes to cam_out structure
   call cam_fluxes_sw%export_surface_fluxes(cam_out)

   ! copy surface lw fluxes to cam_out structure
   call cam_fluxes_lw%export_surface_fluxes(cam_out)

   ! deallocate cam fluxes
   call cam_fluxes_sw%finalize()
   call cam_fluxes_lw%finalize()
   call cam_fluxes_sw_clearsky%finalize()
   call cam_fluxes_lw_clearsky%finalize()

end subroutine radiation_tend


subroutine set_daytime_indices(coszrs, day_indices)
   ! Input: cosine of solar zenith angle
   real(r8), intent(in) :: coszrs(:)

   ! Output: array of indices to daytime columns
   integer, intent(inout) :: day_indices(:)

   ! Loop indices; icol is index to physics columns in current chunk and iday is
   ! index to daytime indices
   integer :: icol, iday

   ! Subroutine name for error messages
   character(len=128) :: sub_name = 'set_daytime_indices'

   ! Initialize array of daytime indices to be all zero. If any zeros exist when
   ! we are done, something went wrong.
   day_indices(:) = 0

   ! Loop over columns and identify daytime columns as those where the cosine
   ! solar zenith angle exceeds zero. Note that we wrap the setting of
   ! day_indices in an if-then to make sure we are not accesing day_indices out
   ! of bounds, and stopping with an informative error message if we do for some
   ! reason.
   iday = 0
   do icol = 1,size(coszrs)
      if (coszrs(icol) > 0._r8) then
         iday = iday + 1
         if (iday <= size(day_indices)) then
            day_indices(iday) = icol
         else
            call endrun(trim(sub_name) // ': iday > size(day_indices)')
         end if
      end if
   end do

   ! Check indices
   call assert_range(day_indices, 1, size(coszrs), trim(sub_name) // 'day_indices')

end subroutine set_daytime_indices


#ifdef SEPARATE_SOLIN_CALC
! Subroutine to calculate the solar insolation assumed by RRTMG? Why did this
! exist in the previous version? Why not just get the incident flux from
! solin(:) = fluxes%flux_dn(:,ktop) ?
subroutine get_solar_insolation(coszrs, day_indices)
   use radconstants,  only: get_ref_solar_band_irrad
   use rad_solar_var, only: get_variability

   real(r8), intent(in) :: coszrs(:)
   integer, intent(in) :: day_indices(:)

   real(r8), allocatable :: solar_variability_factor(:)
   real(r8), allocatable :: solar_band_irradiance(:)
   real(r8) :: eccentricity_factor
   integer :: nday, iday, iband

   ! Define solar incident radiation
   call get_ref_solar_band_irrad(solar_band_irrad)
   call get_variability(sfac)

   ! Get eccentricity factor
   call shr_orb_decl(calday, eccen, mvelpp, lambm0, obliqr, &
                     delta, eccf)

   solin_day = 0._r8
   do i = 1, nday
      do ib = 1, nswbands
         bnd_irrad = sfac(ib) * solar_band_irrad(ib) * eccf * coszrs_day(i)
         solin_day(i) = solin_day(i) + bnd_irrad
      end do
   end do

   solin = 0._r8
   do i = 1, nday
      solin(idxday(i)) = solin_day(i)
   end do

end subroutine get_solar_insolation
#endif


!subroutine output_cloud_optics_sw(state, pbuf, cloud_optics, daytime_indices)
!
!   use physics_types, only: physics_state
!   use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
!
!   type(physics_state), intent(in) :: state
!   type(physics_buffer_desc), pointer :: pbuf(:)
!   type(ty_optical_props_2str), intent(in) :: cloud_optics
!
!   real(r8), pointer :: cloud_fraction(:,:)
!   real(r8), allocatable :: gridbox_cloud_optical_depth(:,:)
!
!   ! Allocate gridbox mean cloud optical depth
!   allocate(gridbox_cloud_optical_depth, source=cloud_fraction)
!
!   ! Multiply by cloud fraction to the the gridbox mean value
!   call pbuf_get_field(pbuf, pbuf_get_index('CLD'), cloud_fraction)
!
!   do i = 1,size(daytime_indices)
!      gridbox_cloud_optical_depth(day_indices(i),:) &
!         = cloud_optics%tau(i,:) &
!         * cloud_fraction(day_indices(i),:)
!   end do
!      
!
!end subroutine



function get_cosine_solar_zenith_angle(state, dt) result(coszrs)
   use physics_types,   only: physics_state
   use phys_grid,       only: get_rlat_all_p, get_rlon_all_p
   use time_manager,    only: get_curr_calday
   use orbit,           only: zenith

   ! Inputs
   type(physics_state), intent(in) :: state
   real(r8), intent(in) :: dt

   ! Output cosine solar zenith angle
   real(r8) :: coszrs(state%ncol)

   ! Local variables
   real(r8) :: calday            ! current calendar day
   real(r8) :: clat(state%ncol)  ! current latitudes(radians)
   real(r8) :: clon(state%ncol)  ! current longitudes(radians)

   calday = get_curr_calday()
   call get_rlat_all_p(state%lchnk, state%ncol, clat)
   call get_rlon_all_p(state%lchnk, state%ncol, clon)
   call zenith(calday, clat, clon, coszrs, state%ncol, dt)
end function get_cosine_solar_zenith_angle


! Get and set cosine of the solar zenith angle or all columns in a physics chuck
! based on input physics_state object and timestep. This routine serves mainly
! as a wrapper for the "zenith" subroutine that handles the task of grabbing the
! appropriate calendar day for this timestep and the latitude and longitude 
! values for the columns in the current chunk, and containing addition CAM
! module use to improve readability of the main radiation_tend routine.
subroutine set_cosine_solar_zenith_angle(state, dt, coszrs)
   use physics_types,   only: physics_state
   use phys_grid,       only: get_rlat_all_p, get_rlon_all_p
   use time_manager,    only: get_curr_calday
   use orbit,           only: zenith

   ! Inputs
   type(physics_state), intent(in) :: state
   real(r8), intent(in) :: dt

   ! Output cosine solar zenith angle. Note that this is declared as an
   ! assumed-shape array, but below we will require coszrs to be at least as
   ! large as state%ncol (because latitudes and longitudes are only defined for
   ! the state%ncol columns that exist in the current physics "chunk")
   real(r8), intent(inout) :: coszrs(:)

   ! Local variables
   real(r8) :: calday  ! current calendar day
   real(r8) :: clat(size(coszrs))  ! current latitudes(radians)
   real(r8) :: clon(size(coszrs))  ! current longitudes(radians)

   ! Make sure desired coszrs has the correct shape. The "zenith" subroutine
   ! expects input arrays to have shape state%ncol, although this should
   ! probably be relaxed to take assumed-shape arrays.
   call assert(size(coszrs) >= state%ncol, 'size(coszrs) < ncol')

   ! Get solar zenith angle from CAM utility routine. The "zenith" routine needs
   ! the current calendar day, and the latitude and longitudes for all columns
   ! in the current physics "chunk", so we use CAM routines to retrieve those
   ! values here.
   calday = get_curr_calday()
   call get_rlat_all_p(state%lchnk, state%ncol, clat(:state%ncol))
   call get_rlon_all_p(state%lchnk, state%ncol, clon(:state%ncol))

   ! Call zenith to calculate cosine solar zenith angle. Note we only pass the
   ! first state%ncol columns in case coszrs was allocated to be larger than
   ! ncol in the calling routine (i.e, pcols).
   call zenith(calday, clat(:state%ncol), clon(:state%ncol), &
               coszrs(:state%ncol), state%ncol, dt)
end subroutine set_cosine_solar_zenith_angle


subroutine check_pbuf_fields(ncol, pbuf, fields)
   use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
   integer, intent(in) :: ncol
   type(physics_buffer_desc), pointer :: pbuf(:)
   character(len=*), intent(in) :: fields(:)
   integer :: i
   real(r8), pointer :: x(:)
   character(len=32) :: sub_name = 'check_pbuf_fields'

   do i = 1,size(fields,1)
      call pbuf_get_field(pbuf, pbuf_get_index(fields(i)), x)
      call assert_valid(x(:ncol), trim(sub_name) // ': ' // fields(i))
   end do
end subroutine check_pbuf_fields


subroutine check_surface_fluxes(ncol, cam_out)
   use camsrfexch, only: cam_out_t
   type(cam_out_t), intent(in) :: cam_out
   integer, intent(in) :: ncol

   call assert_valid(cam_out%netsw(:ncol), 'cam_out%netsw')
   call assert_valid(cam_out%soll(:ncol), 'cam_out%soll')
   call assert_valid(cam_out%sols(:ncol), 'cam_out%sols')
   call assert_valid(cam_out%solld(:ncol), 'cam_out%solld')
   call assert_valid(cam_out%solsd(:ncol), 'cam_out%solsd')
   call assert_valid(cam_out%flwds(:ncol), 'cam_out%flwds')
end subroutine check_surface_fluxes


subroutine get_albedo(nday, day_indices, cam_in, alb_dir, alb_dif)
   use camsrfexch, only: cam_in_t

   integer, intent(in) :: nday
   integer, intent(in) :: day_indices(:)
   type(cam_in_t), intent(in) :: cam_in
   real(r8), intent(inout) :: alb_dir(nswbands,nday)   ! surface albedo, direct radiation
   real(r8), intent(inout) :: alb_dif(nswbands,nday)   ! surface albedo, diffuse radiation

   ! Local namespace
   integer :: i

   ! Surface albedo (band mapping is hardcoded for RRTMG(P) code)
   ! This mapping assumes nswbands=14.
   if (nswbands /= 14) then
      call endrun(module_name //': ERROR: albedo band mapping assumes nswbands=14')
   end if

   do i = 1,nday

      ! Near-IR bands (1-9 and 14), 820-16000 cm-1, 0.625-12.195 microns
      alb_dir(1:8,i) = cam_in%aldir(day_indices(i))
      alb_dif(1:8,i) = cam_in%aldif(day_indices(i))
      alb_dir(14,i)  = cam_in%aldir(day_indices(i))
      alb_dif(14,i)  = cam_in%aldif(day_indices(i))

      ! Set band 24 (or, band 9 counting from 1) to use linear average of UV/visible
      ! and near-IR values, since this band straddles 0.7 microns:
      alb_dir(9,i) = 0.5_r8*(cam_in%aldir(day_indices(i)) + cam_in%asdir(day_indices(i)))
      alb_dif(9,i) = 0.5_r8*(cam_in%aldif(day_indices(i)) + cam_in%asdif(day_indices(i)))

      ! UV/visible bands 25-28 (10-13), 16000-50000 cm-1, 0.200-0.625 micron
      alb_dir(10:13,i) = cam_in%asdir(day_indices(i))
      alb_dif(10:13,i) = cam_in%asdif(day_indices(i))
   enddo

   ! Check values
   call clip_values(alb_dir, 0._r8, 1._r8, varname='alb_dir', warn=.true.)
   call clip_values(alb_dif, 0._r8, 1._r8, varname='alb_dif', warn=.true.)

end subroutine get_albedo


!-------------------------------------------------------------------------------
! Type-bound procedures for cam_fluxes_type class
!-------------------------------------------------------------------------------
subroutine cam_fluxes_export_surface_fluxes(this, cam_out)
   use camsrfexch, only: cam_out_t

   class(cam_flux_type), intent(in) :: this
   type(cam_out_t), intent(inout) :: cam_out
   integer :: i

   ! Export surface fluxes
   ! The RRTMGP output isn't broken down into direct and diffuse.  For first cut
   ! put entire flux into direct and leave the diffuse set to zero.  To break the
   ! fluxes into the UV/vis and near-IR bands use the same scheme as for the albedos
   ! which is hardcoded for 14 spectral bands.
   !
   ! sols(pcols)      Direct solar rad on surface (< 0.7)
   ! soll(pcols)      Direct solar rad on surface (>= 0.7)
   !
   ! Near-IR bands (1-9 and 14), 820-16000 cm-1, 0.625-12.195 microns
   !
   ! Put half of band 9 in each of the UV/visible and near-IR values,
   ! since this band straddles 0.7 microns:
   !
   ! UV/visible bands 10-13, 16000-50000 cm-1, 0.200-0.625 micron
   if (trim(this%flux_type) == 'shortwave') then
      cam_out%soll = 0
      cam_out%sols = 0
      cam_out%solld = 0
      cam_out%solsd = 0
      do i = 1,size(this%bnd_flux_dn, 1)
         cam_out%soll(i) &
            = sum(this%bnd_flux_dn(i,this%kbot,1:8)) &
            + 0.5_r8 * this%bnd_flux_dn(i,this%kbot,9) &
            + this%bnd_flux_dn(i,this%kbot,14)

         cam_out%sols(i) &
            = 0.5_r8 * this%bnd_flux_dn(i,this%kbot,9) &
            + sum(this%bnd_flux_dn(i,this%kbot,10:13))

         cam_out%netsw(i) = this%flux_net_bot(i)
      end do
   else if (trim(this%flux_type) == 'longwave') then
      do i = 1,size(this%flux_dn, 1)
         cam_out%flwds(i) = this%flux_dn(i,this%kbot)
      end do
   else
      call endrun('flux_type ' // this%flux_type // ' not known.')
   end if

end subroutine cam_fluxes_export_surface_fluxes


subroutine cam_fluxes_set_fluxes(this, flux_input, cam_indices)
   use mo_fluxes_byband, only:  ty_fluxes_byband
   class(cam_flux_type), intent(inout) :: this
   type(ty_fluxes_byband), intent(in) :: flux_input
   integer, intent(in), optional :: cam_indices(:)

   integer :: k_cam, k_rad, i
   character(len=32) :: sub_name = 'cam_fluxes_set_fluxes'

   ! DEBUG checks on inputs 
   call assert_valid(flux_input%flux_dn, sub_name // ': flux_input%flux_dn')
   call assert_valid(flux_input%flux_up, sub_name // ': flux_input%flux_up')
   call assert_valid(flux_input%flux_net, sub_name // ': flux_input%flux_net')
   call assert_valid(flux_input%bnd_flux_dn, sub_name // ': flux_input%bnd_flux_dn')
   call assert_valid(flux_input%bnd_flux_up, sub_name // ': flux_input%bnd_flux_up')
   call assert_valid(flux_input%bnd_flux_net, sub_name // ': flux_input%bnd_flux_net')

   ! Map fluxes on radiation grid to cam grid
   if (present(cam_indices)) then

      ! Sanity check on cam_indices
      if (size(cam_indices) /= size(flux_input%flux_dn, 1)) then
         call endrun(sub_name // 'indices do not conform.')
      end if

      do i = 1,size(cam_indices)
         k_rad = 1
         do k_cam = this%ktop,this%kbot
            ! Broadband fluxes
            this%flux_dn(cam_indices(i),k_cam) = flux_input%flux_dn(i,k_rad)
            this%flux_up(cam_indices(i),k_cam) = flux_input%flux_up(i,k_rad)
            this%flux_net(cam_indices(i),k_cam) = flux_input%flux_net(i,k_rad)

            ! Band-by-band fluxes
            this%bnd_flux_dn(cam_indices(i),k_cam,:) = flux_input%bnd_flux_dn(i,k_rad,:)
            this%bnd_flux_up(cam_indices(i),k_cam,:) = flux_input%bnd_flux_up(i,k_rad,:)
            this%bnd_flux_net(cam_indices(i),k_cam,:) = flux_input%bnd_flux_net(i,k_rad,:)

            ! Increment rad level index
            k_rad = k_rad + 1
         end do
      end do
   else
      k_rad = 1
      do k_cam = this%ktop,this%kbot
         ! Broadband fluxes
         this%flux_dn(:,k_cam) = flux_input%flux_dn(:,k_rad)
         this%flux_up(:,k_cam) = flux_input%flux_up(:,k_rad)
         this%flux_net(:,k_cam) = flux_input%flux_net(:,k_rad)

         ! Band-by-band fluxes
         this%bnd_flux_dn(:,k_cam,:) = flux_input%bnd_flux_dn(:,k_rad,:)
         this%bnd_flux_up(:,k_cam,:) = flux_input%bnd_flux_up(:,k_rad,:)
         this%bnd_flux_net(:,k_cam,:) = flux_input%bnd_flux_net(:,k_rad,:)

         ! Increment rad level index
         k_rad = k_rad + 1
      end do
   end if

   ! Compute net fluxes at surface and top of model
   if (trim(this%flux_type) == 'longwave') then
      ! Net fluxes are always down minus up in RRTMGP, but CAM expects upward to
      ! be positive for longwave, so we invert the net flux here
      this%flux_net(:,:) = - this%flux_net(:,:)
      this%flux_net_bot(:) = this%flux_up(:,this%kbot) - this%flux_dn(:,this%kbot)
      this%flux_net_top(:) = this%flux_up(:,this%ktop) - this%flux_dn(:,this%ktop)
   else if (trim(this%flux_type) == 'shortwave') then
      this%flux_net_bot(:) = this%flux_dn(:,this%kbot) - this%flux_up(:,this%kbot)
      this%flux_net_top(:) = this%flux_dn(:,this%ktop) - this%flux_up(:,this%ktop)
   else
      call endrun('flux_type ' // this%flux_type // ' not known')
   end if

   ! Check outputs to make sure we set things right
   call assert_valid(this%flux_up, trim(sub_name) // ': flux_up')
   call assert_valid(this%flux_dn, trim(sub_name) // ': flux_dn')
   call assert_valid(this%flux_net, trim(sub_name) // ': flux_net')
   call assert_valid(this%bnd_flux_up, trim(sub_name) // ': bnd_flux_up')
   call assert_valid(this%bnd_flux_dn, trim(sub_name) // ': bnd_flux_dn')
   call assert_valid(this%bnd_flux_net, trim(sub_name) // ': bnd_flux_net')
   call assert_valid(this%flux_net_bot, trim(sub_name) // ': flux_net_bot')
   call assert_valid(this%flux_net_top, trim(sub_name) // ': flux_net_top')

end subroutine cam_fluxes_set_fluxes


subroutine cam_fluxes_calculate_heating_rate(this, pdel)

   use physconst, only: gravit

   ! Inputs
   class(cam_flux_type), intent(inout) :: this 
   real(r8), intent(in) :: pdel(:,:)

   ! Loop indices
   integer :: i, k

   ! Everyone needs a name
   character(len=32) :: sub_name = 'calculate_heating_rate'

   ! Check inputs
   call assert_valid(this%flux_net, 'calculate_heating_rate: flux_net')

   ! Loop over levels and calculate heating rates; note that the fluxes *should*
   ! be defined at interfaces, so the loop ktop,kbot and grabbing the current
   ! and next value of k should be safe. ktop should be the top interface, and
   ! kbot + 1 should be the bottom interface.
   !
   ! NOTE: to get heating rate in K/day, normally we would use:
   !
   !     H = dF / dp * g * (sec/day) * (1e-5) / (cpair)
   !
   ! Here we just use
   !
   !     H = dF / dp * g
   !
   ! Why? Something to do with convenience with applying the fluxes to the
   ! heating tendency?
   if (trim(this%flux_type) == 'longwave') then
      do i = 1,size(this%flux_net,1)
         do k = 1,size(this%flux_net,2)-1
            this%heating_rate(i,k) = ( &
               this%flux_net(i,k+1) - this%flux_net(i,k) &
            ) * gravit / pdel(i,k)
         end do
      end do
   else if (trim(this%flux_type) == 'shortwave') then
      do i = 1,size(this%flux_net,1)
         do k = 1,size(this%flux_net,2)-1
            this%heating_rate(i,k) = ( &
               this%flux_net(i,k) - this%flux_net(i,k+1) &
            ) * gravit / pdel(i,k)
         end do
      end do
   else
      call endrun(sub_name // ': flux_type ' // trim(this%flux_type) // ' not defined.')
   end if

end subroutine cam_fluxes_calculate_heating_rate
      

subroutine cam_fluxes_to_pbuf(this, pbuf)
   use physics_buffer,  only: physics_buffer_desc, pbuf_get_field, &
                              pbuf_get_index

   class(cam_flux_type), intent(in) :: this
   type(physics_buffer_desc), pointer :: pbuf(:)

   ! Pointers to pbuf fields
   real(r8), pointer :: fsds(:)
   real(r8), pointer :: fsns(:)
   real(r8), pointer :: fsnt(:)
   real(r8), pointer :: flns(:)
   real(r8), pointer :: flnt(:)

   integer :: ncol, nlevels
   character(len=32) :: sub_name = 'cam_fluxes_to_pbuf'

   ! Check inputs
   call assert_valid(this%flux_dn, trim(sub_name) // ': this%flux_dn')
   call assert_valid(this%flux_net_bot, trim(sub_name) // ': this%flux_net_bot')
   call assert_valid(this%flux_net_top, trim(sub_name) // ': this%flux_net_top')

   call assert_valid(this%flux_up, trim(sub_name) // ': flux_up')
   call assert_valid(this%flux_dn, trim(sub_name) // ': flux_dn')
   call assert_valid(this%flux_net, trim(sub_name) // ': flux_net')
   call assert_valid(this%bnd_flux_up, trim(sub_name) // ': bnd_flux_up')
   call assert_valid(this%bnd_flux_dn, trim(sub_name) // ': bnd_flux_dn')
   call assert_valid(this%bnd_flux_net, trim(sub_name) // ': bnd_flux_net')
   call assert_valid(this%flux_net_bot, trim(sub_name) // ': flux_net_bot')
   call assert_valid(this%flux_net_top, trim(sub_name) // ': flux_net_top')

   ncol = size(this%flux_up, 1)
   nlevels = size(this%flux_up, 2)
   if (trim(this%flux_type) == 'shortwave') then
      ! Associate pointers
      call pbuf_get_field(pbuf, pbuf_get_index('FSDS'), fsds)
      call pbuf_get_field(pbuf, pbuf_get_index('FSNS'), fsns)
      call pbuf_get_field(pbuf, pbuf_get_index('FSNT'), fsnt)

      ! Copy data
      if (size(this%flux_dn, 2) >= this%kbot) then
         fsds(:ncol) = this%flux_dn(:ncol,this%kbot)
      else
         call endrun(trim(sub_name) // ': number of levels < kbot')
      end if
      fsns(:ncol) = this%flux_net_bot(:ncol)
      fsnt(:ncol) = this%flux_net_top(:ncol)

      ! Check that we set them right
      call assert_valid(fsds(:ncol), trim(sub_name) // ': fsds')
      call assert_valid(fsns(:ncol), trim(sub_name) // ': fsns')
      call assert_valid(fsnt(:ncol), trim(sub_name) // ': fsnt')
   else if (trim(this%flux_type) == 'longwave') then
      ! Associate pointers
      call pbuf_get_field(pbuf, pbuf_get_index('FLNS'), flns)
      call pbuf_get_field(pbuf, pbuf_get_index('FLNT'), flnt)

      ! Copy data
      flns(:ncol) = this%flux_net_bot(:ncol)
      flnt(:ncol) = this%flux_net_top(:ncol)

      ! Check that we set them right
      call assert_valid(flns(:ncol), trim(sub_name) // ': flns')
      call assert_valid(flnt(:ncol), trim(sub_name) // ': flnt')
   else
      call endrun(trim(sub_name) // ': flux_type ' // this%flux_type // ' not known.')
   end if

end subroutine cam_fluxes_to_pbuf


subroutine cam_fluxes_initialize(this, ncol, nlevels, nbands, ktop, kbot, flux_type)
   class(cam_flux_type), intent(inout) :: this
   integer, intent(in) :: ncol, nlevels, nbands, ktop, kbot
   character(len=*), intent(in) :: flux_type

   this%ktop = ktop
   this%kbot = kbot
   this%flux_type = flux_type

   ! allocate CAM fluxes
   allocate(this%flux_up(ncol,nlevels), &
            this%flux_dn(ncol,nlevels), &
            this%flux_net(ncol,nlevels), &
            this%bnd_flux_up(ncol,nlevels,nbands), &
            this%bnd_flux_dn(ncol,nlevels,nbands), &
            this%bnd_flux_net(ncol,nlevels,nbands), &
            this%flux_net_bot(ncol), &
            this%flux_net_top(ncol))

   ! allocate CAM heating rate
   allocate(this%heating_rate(ncol,nlevels-1))

   ! initialize to zero
   this%flux_up = 0.0
   this%flux_dn = 0.0
   this%flux_net = 0.0
   this%bnd_flux_up = 0.0
   this%bnd_flux_dn = 0.0
   this%bnd_flux_net = 0.0
   this%flux_net_bot = 0.0
   this%flux_net_top = 0.0
   this%heating_rate = 0.0

end subroutine cam_fluxes_initialize


subroutine cam_fluxes_finalize(this)
   class(cam_flux_type), intent(inout) :: this
   deallocate(this%flux_up, this%flux_dn, this%flux_net, &
              this%bnd_flux_up, this%bnd_flux_dn, this%bnd_flux_net, &
              this%flux_net_bot, this%flux_net_top, &
              this%heating_rate)
end subroutine cam_fluxes_finalize


!-------------------------------------------------------------------------------

subroutine clip_values_1d(x, min_x, max_x, varname, warn)
   real(r8), intent(inout) :: x(:)
   real(r8), intent(in) :: min_x
   real(r8), intent(in) :: max_x
   character(len=*), intent(in), optional :: varname
   logical, intent(in), optional :: warn

   logical :: warn_local

   warn_local = .false.
   if (present(warn)) then
      warn_local = warn
   end if

   ! Look for values less than threshold
   if (any(x < min_x)) then
      ! Raise warning?
      if (warn_local) then
         if (present(varname)) then
            print *, module_name // ' warning: ', &
                     count(x < min_x), ' values are below threshold for variable ', &
                     trim(varname), '; min = ', minval(x)
         else
            print *, module_name // ' warning: ', &
                     count(x < min_x), ' values are below threshold; min = ', minval(x)
         end if
      end if

      ! Clip values
      where (x < min_x)
         x = min_x
      endwhere
   end if

   ! Look for values greater than threshold 
   if (any(x > max_x)) then 
      ! Raise warning?
      if (warn_local) then
         if (present(varname)) then
            print *, module_name // ' warning: ', &
                     count(x > max_x), ' values are above threshold for variable ', &
                     trim(varname), '; max = ', maxval(x)
         else
            print *, module_name // ' warning: ', &
                     count(x > max_x), ' values are above threshold; max = ', maxval(x)
         end if
      end if

      ! Clip values
      where (x > max_x)
         x = max_x
      end where
   end if
end subroutine clip_values_1d

subroutine clip_values_2d(x, min_x, max_x, varname, warn)
   real(r8), intent(inout) :: x(:,:)
   real(r8), intent(in) :: min_x
   real(r8), intent(in) :: max_x
   character(len=*), intent(in), optional :: varname
   logical, intent(in), optional :: warn

   integer :: i, j, nbad
   logical :: warn_local

   warn_local = .false.
   if (present(warn)) then
      warn_local = warn
   end if

   ! look for values less than threshold
   if (any(x < min_x)) then
      ! Raise warning?
      if (warn_local) then
         if (present(varname)) then
            print *, module_name // ' warning: ', &
                     count(x < min_x), ' values are below threshold for variable ', &
                     trim(varname), '; min = ', minval(x)
         else
            print *, module_name // ' warning: ', &
                     count(x < min_x), ' values are below threshold; min = ', minval(x)
         end if
      end if

      ! Clip values
      where (x < min_x)
         x = min_x
      endwhere
   end if

   ! Look for values greater than threshold
   if (any(x > max_x)) then 
      ! Raise warning?
      if (warn_local) then
         if (present(varname)) then
            print *, module_name // ' warning: ', &
                     count(x > max_x), ' values are above threshold for variable ', &
                     trim(varname), '; max = ', maxval(x)
         else
            print *, module_name // ' warning: ', &
                     count(x > max_x), ' values are above threshold; max = ', maxval(x)
         end if
      end if

      ! Clip values
      where (x > max_x)
         x = max_x
      end where
   end if

!  nbad = 0
!  do i = 1,size(x,1)
!     do j = 1,size(x,2)
!        if (x(i,j) < min_x) then
!           x(i,j) = min_x
!           nbad = nbad + 1
!           if (warn_local) then
!              if (present(varname)) then
!                 print *, trim(varname) // ' < min at i,j = ', i, j
!              end if
!           end if
!        else if (x(i,j) > max_x) then
!           x(i,j) = max_x
!           nbad = nbad + 1
!           if (warn_local) then
!              if (present(varname)) then
!                 print *, trim(varname) // ' > min at i,j = ', i, j
!              end if
!           end if
!        end if
!     end do
!  end do

!  if (nbad > 0) then
!     if (warn_local) then
!        if (present(varname)) then
!           print *, 'Warning: ' // varname // ' values clipped. min/max: ', &
!              minval(x), maxval(x)
!        else
!           print *, 'Warning: values clipped. min/max: ', &
!              minval(x), maxval(x)
!        end if
!     end if
!  end if
!           
end subroutine clip_values_2d


! Send shortwave fluxes and heating rates to history buffer
subroutine radiation_output_sw(icall, state, flux_allsky, flux_clearsky)
   use physconst, only: cpair
   use physics_types, only: physics_state
   use cam_history, only: outfld
   
   integer, intent(in) :: icall
   type(physics_state), intent(in) :: state
   type(cam_flux_type), intent(in) :: flux_allsky
   type(cam_flux_type), intent(in) :: flux_clearsky

   real(r8) :: fus(pcols,pverp), fds(pcols,pverp), fns(pcols,pverp)
   real(r8) :: fsut(pcols), fsnt(pcols), fsns(pcols), fsds(pcols)
   real(r8) :: fsutc(pcols), fsntc(pcols), fsnsc(pcols)
   integer :: ncol
   integer :: ktop, kbot  ! indices to limits of CAM grid

   ncol = state%ncol

   ! Initialize, since these are size(pcol)
   fus = 0.0
   fds = 0.0
   fns = 0.0
   fsut = 0.0
   fsnt = 0.0
   fsns = 0.0
   fsds = 0.0
   fsutc = 0.0
   fsntc = 0.0
   fsnsc = 0.0
   
   ktop = flux_allsky%ktop
   kbot = flux_allsky%kbot

   fus(:ncol,:) = flux_allsky%flux_up(:ncol,:)
   fds(:ncol,:) = flux_allsky%flux_dn(:ncol,:)
   fns(:ncol,:) = flux_allsky%flux_net(:ncol,:)

   fsnt(:ncol) = flux_allsky%flux_dn(:ncol,ktop) - flux_allsky%flux_up(:ncol,ktop)
   fsns(:ncol) = flux_allsky%flux_dn(:ncol,kbot) - flux_allsky%flux_up(:ncol,kbot)
   fsut(:ncol) = flux_allsky%flux_up(:ncol,ktop)
   fsds(:ncol) = flux_allsky%flux_dn(:ncol,kbot)
   
   fsntc(:ncol) = flux_clearsky%flux_dn(:ncol,ktop) - flux_clearsky%flux_up(:ncol,ktop)
   fsnsc(:ncol) = flux_clearsky%flux_dn(:ncol,kbot) - flux_clearsky%flux_up(:ncol,kbot)
   fsutc(:ncol) = flux_clearsky%flux_up(:ncol,ktop)

   ! Heating rates
   call outfld('QRS'//diag(icall), flux_allsky%heating_rate(:ncol,:pver)/cpair, ncol, state%lchnk)

   ! All-sky fluxes at model interfaces
   call outfld('FDS'//diag(icall), fds, pcols, state%lchnk)
   call outfld('FUS'//diag(icall), fus, pcols, state%lchnk)
   call outfld('FNS'//diag(icall), fns, pcols, state%lchnk)

   ! All-sky fluxes
   call outfld('FSNT'//diag(icall), fsnt, pcols, state%lchnk)
   call outfld('FSNS'//diag(icall), fsns, pcols, state%lchnk)
   call outfld('FSDS'//diag(icall), fsds, pcols, state%lchnk)
   call outfld('FSUTOA'//diag(icall), fsut, pcols, state%lchnk)
   call outfld('FSNTOA'//diag(icall), fsnt, pcols, state%lchnk)

   ! Clear-sky fluxes
   call outfld('FSNSC'//diag(icall), fsnsc, pcols, state%lchnk)
   call outfld('FSUTOAC'//diag(icall), fsutc, pcols, state%lchnk)
   call outfld('FSNTOAC'//diag(icall), fsntc, pcols, state%lchnk)

end subroutine radiation_output_sw


! Send longwave fluxes and heating rates to history buffer
subroutine radiation_output_lw(icall, state, flux_allsky, flux_clearsky)
   use physconst, only: cpair
   use physics_types, only: physics_state
   use cam_history, only: outfld
   
   integer, intent(in) :: icall
   type(physics_state), intent(in) :: state
   type(cam_flux_type), intent(in) :: flux_allsky
   type(cam_flux_type), intent(in) :: flux_clearsky

   real(r8) :: ful(pcols,pverp), fdl(pcols,pverp), fnl(pcols,pverp)
   real(r8) :: flut(pcols), flnt(pcols), flns(pcols), flds(pcols)
   real(r8) :: flntc(pcols), flnsc(pcols)
   integer :: ncol
   integer :: ktop, kbot  ! indices to limits of CAM grid INTERFACES

   ncol = state%ncol

   ! Initialize, since these are size(pcol)
   ful = 0.0
   fdl = 0.0
   fnl = 0.0
   flut = 0.0
   flnt = 0.0
   flns = 0.0
   flds = 0.0
   flntc = 0.0
   flnsc = 0.0
   
   ktop = flux_allsky%ktop
   kbot = flux_allsky%kbot

   ful(:ncol,:) = flux_allsky%flux_up(:ncol,:)
   fdl(:ncol,:) = flux_allsky%flux_dn(:ncol,:)
   fnl(:ncol,:) = flux_allsky%flux_net(:ncol,:)

   flnt(:ncol) = flux_allsky%flux_up(:ncol,ktop) - flux_allsky%flux_dn(:ncol,ktop)
   flns(:ncol) = flux_allsky%flux_up(:ncol,kbot) - flux_allsky%flux_dn(:ncol,kbot)
   flut(:ncol) = flux_allsky%flux_up(:ncol,ktop)
   flds(:ncol) = flux_allsky%flux_dn(:ncol,kbot)
   
   ! Heating rates
   call outfld('QRL'//diag(icall), flux_allsky%heating_rate(:ncol,:pver)/cpair, ncol, state%lchnk)

   ! All-sky fluxes on model levels (for debugging)
   call outfld('FDL'//diag(icall), fdl, pcols, state%lchnk)
   call outfld('FUL'//diag(icall), ful, pcols, state%lchnk)
   call outfld('FNL'//diag(icall), fnl, pcols, state%lchnk)

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

   ! Allocate flux arrays
   ! NOTE: fluxes defined at interfaces, so need to either pass nlevels as
   ! number of model levels plus one, or allocate as nlevels+1 if nlevels
   ! represents number of model levels rather than number of interface levels.

   ! Broadband fluxes
   allocate(fluxes%flux_up(ncol, nlevels))
   allocate(fluxes%flux_dn(ncol, nlevels))
   allocate(fluxes%flux_net(ncol, nlevels))
   !allocate(fluxes%flux_dn_dir(ncol, nlevels))

   ! Fluxes by band
   allocate(fluxes%bnd_flux_up(ncol, nlevels, nbands))
   allocate(fluxes%bnd_flux_dn(ncol, nlevels, nbands))
   allocate(fluxes%bnd_flux_net(ncol, nlevels, nbands))
   !allocate(fluxes%bnd_flux_dn_dir(ncol, nlevels, nbands))

   ! Initialize
   fluxes%flux_up = 0
   fluxes%flux_dn = 0
   fluxes%flux_net = 0

   fluxes%bnd_flux_up = 0
   fluxes%bnd_flux_dn = 0
   fluxes%bnd_flux_net = 0
end subroutine initialize_rrtmgp_fluxes


! Populate RRTMGP optical properties object with CAM-computed aerosol
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
end subroutine set_aerosol_optics_lw


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
   use mo_util_string, only: lower_case

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

      ! Populate the RRTMGP gas concentration object with values for this gas.
      ! NOTE: RRTMGP makes some assumptions about gas names internally, so we
      ! need to pass the gas names as LOWER case here!
      call handle_rrtmgp_error( &
         gas_concentrations%set_vmr(trim(lower_case(gas_species(igas))),vol_mix_ratio) & 
      )

   end do

   ! Free memory allocated to the volume mixing ratio array
   deallocate(vol_mix_ratio)
   
end subroutine set_gas_concentrations


subroutine check_gas_concentrations(gas_concentrations)
   use mo_gas_concentrations, only: ty_gas_concs

   type(ty_gas_concs), intent(in) :: gas_concentrations
   integer :: igas
   do igas = 1,size(gas_concentrations%concs,1)
      call assert_valid(gas_concentrations%concs(igas)%conc, gas_concentrations%gas_name(igas))

      print *, 'check_gas_concentrations ' // gas_concentrations%gas_name(igas) // ': ', &
               minval(gas_concentrations%concs(igas)%conc), &
               maxval(gas_concentrations%concs(igas)%conc)
   end do
end subroutine


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


subroutine set_cloud_optics_sw(day_indices, ktop, kbot, kdist, state, pbuf, cloud_sw)
   
   use ppgrid, only: pcols, pver, pverp
   use physics_types, only: physics_state
   use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
   use mo_optical_props, only: ty_optical_props_2str
   use mo_gas_optics, only: ty_gas_optics_specification
   use mcica_subcol_gen, only: mcica_subcol_sw

   type(ty_gas_optics_specification), intent(in) :: kdist
   type(physics_state), intent(in) :: state
   type(physics_buffer_desc), pointer :: pbuf(:)
   type(ty_optical_props_2str), intent(inout) :: cloud_sw
   integer, intent(in) :: day_indices(:)
   integer, intent(in) :: ktop, kbot

   ! For MCICA sampling routine
   integer, parameter :: changeseed = 1

   ! Dimension sizes
   integer :: nbnd, ngpt, ncol, nday, nlay

   ! Temporary arrays to hold gridbox-mean cloud optics (nbnd,ncol,pver)
   real(r8), allocatable :: cloud_optical_depth(:,:,:)
   real(r8), allocatable :: single_scattering_albedo(:,:,:)
   real(r8), allocatable :: assymmetry_parameter(:,:,:)

   ! Temporary array to hold gridbox-mean cloud fraction (nday,nlev)
   ! NOTE: pointer because this comes from pbuf
   real(r8), pointer :: cloud_fraction(:,:)

   ! Temporary arrays to hold gridbox-mean cloud optics subset to the radiation
   ! grid (nbnd,nday,nlay)
   real(r8), allocatable :: cloud_optical_depth_rad(:,:,:)
   real(r8), allocatable :: single_scattering_albedo_rad(:,:,:)
   real(r8), allocatable :: assymmetry_parameter_rad(:,:,:)

   ! Temporary array to hold gridbox-mean pressure and cloud fraction (nday,nlev)
   real(r8), allocatable :: pmid_rad(:,:)
   real(r8), allocatable :: cloud_fraction_rad(:,:)

   ! Loop variables
   integer :: i

   ! Initialize cloud optics object
   ! TODO: put initialization/allocation here
   cloud_sw%tau = 0.0
   cloud_sw%ssa = 1.0
   cloud_sw%g = 0.0

   ! Allocate temporary arrays
   ! NOTE: get_cloud_optics_sw returns optics in BANDS, not G-POINTS...we will
   ! need to map these to g-points!
   ncol = state%ncol
   nbnd = nswbands
   allocate(cloud_optical_depth(nbnd,ncol,pver), &
            single_scattering_albedo(nbnd,ncol,pver), &
            assymmetry_parameter(nbnd,ncol,pver))

   ! Get cloud optics (for all CAM columns in this chunk and all levels)
   call get_cloud_optics_sw(state, pbuf, &
                            cloud_optical_depth, &
                            single_scattering_albedo, &
                            assymmetry_parameter)

   ! Allocate temporary arrays on radiation grid
   nday = size(day_indices)
   nlay = abs(kbot - ktop)
   allocate(pmid_rad(nday,nlay), &
            cloud_fraction_rad(nday,nlay), &
            cloud_optical_depth_rad(nbnd,nday,nlay), &
            single_scattering_albedo_rad(nbnd,nday,nlay), &
            assymmetry_parameter_rad(nbnd,nday,nlay))

   ! Map inputs to radiation grid
   call pbuf_get_field(pbuf, pbuf_get_index('CLD'), cloud_fraction)
   do i = 1,nday
      pmid_rad(i,:nlay) = state%pmid(day_indices(i),ktop:kbot-1)
      cloud_fraction_rad(i,:nlay) = cloud_fraction(day_indices(i),ktop:kbot-1)
      cloud_optical_depth_rad(:nbnd,i,:nlay) = cloud_optical_depth(:nbnd,day_indices(i),ktop:kbot-1)
      single_scattering_albedo_rad(:nbnd,i,:nlay) = single_scattering_albedo(:nbnd,day_indices(i),ktop:kbot-1)
      assymmetry_parameter_rad(:nbnd,i,:nlay) = assymmetry_parameter(:nbnd,day_indices(i),ktop:kbot-1)
   end do

   ! Deallocate temporary arrays we are no longer using
   deallocate(cloud_optical_depth, &
              single_scattering_albedo, &
              assymmetry_parameter)

   ! Do MCICA sampling of optics here (ONLY those on the radiation grid)
   ! TODO: I think this also maps bands to g-points, but need to confirm this.
   ngpt = kdist%get_ngpt()
   call mcica_subcol_sw( &
      kdist, nbnd, ngpt, nday, nlay, nlay, changeseed, &
      pmid_rad, cloud_fraction_rad, &
      cloud_optical_depth_rad, single_scattering_albedo_rad, assymmetry_parameter_rad, &
      cloud_sw%tau, cloud_sw%ssa, cloud_sw%g &
   )

   ! Deallocate temporary arrays on radiation grid
   deallocate(pmid_rad, &
              cloud_fraction_rad, &
              cloud_optical_depth_rad, &
              single_scattering_albedo_rad, &
              assymmetry_parameter_rad)

end subroutine set_cloud_optics_sw


subroutine get_cloud_optics_sw(state, pbuf, &
                               cloud_optical_depth, &
                               single_scattering_albedo, &
                               assymmetry_parameter) 

   use physics_types, only: physics_state
   use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_get_index
   use cloud_rad_props, only: get_ice_optics_sw, &
                              get_liquid_optics_sw, &
                              get_snow_optics_sw

   type(physics_state), intent(in) :: state
   type(physics_buffer_desc), pointer :: pbuf(:)
   real(r8), intent(inout) :: cloud_optical_depth(:,:,:)
   real(r8), intent(inout) :: single_scattering_albedo(:,:,:)
   real(r8), intent(inout) :: assymmetry_parameter(:,:,:)

   ! Temporary variables to hold liquid cloud optical properties before combining
   ! into output arrays
   real(r8), allocatable :: liquid_tau(:,:,:)
   real(r8), allocatable :: liquid_tau_ssa(:,:,:)
   real(r8), allocatable :: liquid_tau_ssa_g(:,:,:)
   real(r8), allocatable :: liquid_tau_ssa_f(:,:,:)

   ! Temporary variables to hold ice cloud optical properties before combining
   ! into output arrays
   real(r8), allocatable :: ice_tau(:,:,:)
   real(r8), allocatable :: ice_tau_ssa(:,:,:)
   real(r8), allocatable :: ice_tau_ssa_g(:,:,:)
   real(r8), allocatable :: ice_tau_ssa_f(:,:,:)

   ! Temporary variables to hold "snow cloud" optical properties before combining
   ! into output arrays
   real(r8), allocatable :: snow_tau(:,:,:)
   real(r8), allocatable :: snow_tau_ssa(:,:,:)
   real(r8), allocatable :: snow_tau_ssa_g(:,:,:)
   real(r8), allocatable :: snow_tau_ssa_f(:,:,:)

   ! Temporary variables to hold combined optical properties
   real(r8), allocatable :: combined_tau(:,:,:)
   real(r8), allocatable :: combined_tau_ssa(:,:,:)
   real(r8), allocatable :: combined_tau_ssa_g(:,:,:)
   real(r8), allocatable :: combined_tau_ssa_f(:,:,:)

   integer :: nbnd, ncol, nlev
   logical :: do_snow_optics = .false.
   integer :: i, j, k

   ! Allocate temporary arrays
   nbnd = nswbands
   ncol = state%ncol
   nlev = pver
   allocate(liquid_tau(nbnd,pcols,nlev), &
            liquid_tau_ssa(nbnd,pcols,nlev), &
            liquid_tau_ssa_g(nbnd,pcols,nlev), &
            liquid_tau_ssa_f(nbnd,pcols,nlev), &
            ice_tau(nbnd,pcols,nlev), &
            ice_tau_ssa(nbnd,pcols,nlev), &
            ice_tau_ssa_g(nbnd,pcols,nlev), &
            ice_tau_ssa_f(nbnd,pcols,nlev), &
            snow_tau(nbnd,pcols,nlev), &
            snow_tau_ssa(nbnd,pcols,nlev), &
            snow_tau_ssa_g(nbnd,pcols,nlev), &
            snow_tau_ssa_f(nbnd,pcols,nlev), &
            combined_tau(nbnd,pcols,nlev), &
            combined_tau_ssa(nbnd,pcols,nlev), &
            combined_tau_ssa_g(nbnd,pcols,nlev), &
            combined_tau_ssa_f(nbnd,pcols,nlev))

   ! Get ice cloud optics
   call get_ice_optics_sw(state, pbuf, &
                          ice_tau, ice_tau_ssa, &
                          ice_tau_ssa_g, ice_tau_ssa_f)
   
   ! Get liquid cloud optics
   call get_liquid_optics_sw(state, pbuf, &
                             liquid_tau, liquid_tau_ssa, &
                             liquid_tau_ssa_g, liquid_tau_ssa_f)

   ! Should we do snow optics? Check for existence of "cldfsnow" variable
   !call pbuf_get_index(pbuf, 'cldfsnow', err=err)

   ! Get snow cloud optics
   if (do_snow_optics) then
      ! Doing snow optics; call procedure to get these from CAM state and pbuf
      call get_snow_optics_sw(state, pbuf, &
                              snow_tau, snow_tau_ssa, &
                              snow_tau_ssa_g, snow_tau_ssa_f)
   else
      ! We are not doing snow optics, so set these to zero so we can still use 
      ! the arrays without additional logic
      snow_tau(:,:,:) = 0.0
      snow_tau_ssa(:,:,:) = 0.0
      snow_tau_ssa_g(:,:,:) = 0.0
      snow_tau_ssa_f(:,:,:) = 0.0
   end if

   ! Combine all cloud optics from CAM routines
   combined_tau = ice_tau + liquid_tau + snow_tau
   combined_tau_ssa = ice_tau_ssa + liquid_tau_ssa + snow_tau_ssa
   combined_tau_ssa_g = ice_tau_ssa_g + liquid_tau_ssa_g + snow_tau_ssa_g

   ! DEBUG check shapes
   call assert(all(shape(cloud_optical_depth) == (/nbnd,ncol,nlev/)), &
               'cloud_optical_depth size mismatch')
   call assert(all(shape(single_scattering_albedo) == (/nbnd,ncol,nlev/)), &
               'single_scattering_albedo size mismatch')
   call assert(all(shape(assymmetry_parameter) == (/nbnd,ncol,nlev/)), &
               'assymmetry_parameter size mismatch')

   ! Copy to output arrays, converting to optical depth, single scattering
   ! albedo, and assymmetry parameter from the products that the CAM routines
   ! return. Make sure we do not try to divide by zero...
   cloud_optical_depth(:nbnd,:ncol,:nlev) = combined_tau(:nbnd,:ncol,:nlev)
   where (combined_tau(:nbnd,:ncol,:nlev) > 0)
      single_scattering_albedo(:nbnd,:ncol,:nlev) &
         = combined_tau_ssa(:nbnd,:ncol,:nlev) &
         / combined_tau(:nbnd,:ncol,:nlev)
   elsewhere
      single_scattering_albedo(:nbnd,:ncol,:nlev) = 1.0
   endwhere
   where (combined_tau_ssa(:nbnd,:ncol,:nlev) > 0)
      assymmetry_parameter(:nbnd,:ncol,:nlev) &
         = combined_tau_ssa_g(:nbnd,:ncol,:nlev) &
         / combined_tau_ssa(:nbnd,:ncol,:nlev)
   elsewhere
      assymmetry_parameter(:nbnd,:ncol,:nlev) = 0.0
   end where

   ! Free memory allocated to temporary arrays
   deallocate(liquid_tau, &
              liquid_tau_ssa, &
              liquid_tau_ssa_g, &
              liquid_tau_ssa_f, &
              ice_tau, &
              ice_tau_ssa, &
              ice_tau_ssa_g, &
              ice_tau_ssa_f, &
              snow_tau, &
              snow_tau_ssa, &
              snow_tau_ssa_g, &
              snow_tau_ssa_f, &
              combined_tau, &
              combined_tau_ssa, &
              combined_tau_ssa_g, &
              combined_tau_ssa_f)

end subroutine get_cloud_optics_sw


integer function min_pressure_index(pressure, min_pressure)
   real(r8), intent(in) :: pressure(:,:)
   real(r8), intent(in) :: min_pressure
   integer :: k

   do k = 1,size(pressure,2)
      if (all(pressure(:,k) > min_pressure)) then
         min_pressure_index = k
         exit
      end if
   end do
end function

integer function max_pressure_index(pressure, max_pressure)
   real(r8), intent(in) :: pressure(:,:)
   real(r8), intent(in) :: max_pressure
   integer :: k

   do k = size(pressure,2),1,-1
      if (all(pressure(:,k) < max_pressure)) then
         max_pressure_index = k
         exit
      end if
   end do
end function


end module radiation
