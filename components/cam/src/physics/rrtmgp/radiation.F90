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

use shr_kind_mod,     only: r8=>shr_kind_r8, cl=>shr_kind_cl
use ppgrid,           only: pcols, pver, pverp, begchunk, endchunk
use cam_abortutils,   only: endrun
use scamMod,          only: scm_crm_mode, single_column, swrad_off
use rad_constituents, only: N_DIAG

! RRTMGP modules
use mo_gas_optics_specification, only: ty_gas_optics_specification

implicit none
private
save

public :: &
   radiation_register,    &! registers radiation physics buffer fields
   radiation_defaultopts, &! set default values of namelist variables in runtime_opts
   radiation_setopts,     &! set namelist values from runtime_opts
   radiation_printopts,   &! print namelist values to log
   radiation_nextsw_cday, &! calendar day of next radiation calculation
   radiation_do,          &! query which radiation calcs are done this timestep
   radiation_init,        &! calls radini
   radiation_tend          ! moved from radctl.F90

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

! Default values for namelist variables
integer :: iradsw = -1     ! freq. of shortwave radiation calc in time steps (positive)
                           ! or hours (negative).
integer :: iradlw = -1     ! frequency of longwave rad. calc. in time steps (positive)
                           ! or hours (negative).
integer :: irad_always = 0 ! Specifies length of time in timesteps (positive)
                           ! or hours (negative) SW/LW radiation will be
                           ! run continuously from the start of an
                           ! initial or restart run
logical :: spectralflux  = .false. ! calculate fluxes (up and down) per band.

! if true, uses the radiation dt for all cosz calculations
logical :: use_rad_dt_cosz  = .false. 

character(len=16) :: microp_scheme  ! microphysics scheme

character(len=4) :: diag(0:N_DIAG) =(/'    ','_d1 ','_d2 ','_d3 ','_d4 ','_d5 ','_d6 ','_d7 ','_d8 ','_d9 ','_d10'/)

logical :: dohirs = .false. ! diagnostic  brightness temperatures at the top of the
                            ! atmosphere for 7 TOVS/HIRS channels (2,4,6,8,10,11,12) and 4 TOVS/MSU 
                            ! channels (1,2,3,4).
integer :: ihirsfq = 1      ! frequency (timesteps) of brightness temperature calcs

! time step to use for the shr_orb_cosz calculation, if use_rad_dt_cosz set to true
real(r8) :: dt_avg = 0.0_r8  

! k-distribution coefficients
class(ty_gas_optics_specification) :: k_dist_sw, k_dist_lw

! k-distribution coefficients files to read from. These are set via namelist
! variables.
character(len=cl) :: coefficients_file_sw, coefficients_file_lw

!===============================================================================

contains

!===============================================================================

subroutine radiation_readnl(nlfile)

   ! Purpose: Read radiation_nl namelist group.

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_integer, mpi_logical, &
                              mpi_character

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
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'rrtmgp_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, rrtmgp_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(sub // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

   ! Broadcast namelist variables
   call mpi_bcast(rrtmgp_coefficients_file_lw, cl, mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subroutine_name//": FATAL: mpi_bcast: rrtmgp_coefficients_file_lw")
   call mpi_bcast(rrtmgp_coefficients_file_sw, cl, mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subroutine_name//": FATAL: mpi_bcast: rrtmgp_coefficients_file_sw")
   call mpi_bcast(rrtmgp_iradsw, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subroutine_name//": FATAL: mpi_bcast: rrtmgp_iradsw")
   call mpi_bcast(rrtmgp_iradlw, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subroutine_name//": FATAL: mpi_bcast: rrtmgp_iradlw")
   call mpi_bcast(rrtmgp_irad_always, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subroutine_name//": FATAL: mpi_bcast: rrtmgp_irad_always")
   call mpi_bcast(rrtmgp_use_rad_dt_cosz, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subroutine_name//": FATAL: mpi_bcast: rrtmgp_use_rad_dt_cosz")
   call mpi_bcast(rrtmgp_spectralflux, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subroutine_name//": FATAL: mpi_bcast: rrtmgp_spectralflux")

   ! Set module data
   coefficients_file_lw   = rrtmgp_coefficients_file_lw
   coefficients_file_sw   = rrtmgp_coefficients_file_sw
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
   use radconstants,   only: nswbands, nlwbands

   integer :: idx  ! dummy index for adding fields to physics buffer

   ! Heating rate profiles; QRS is the shortwave radiative heating rate, and QRL
   ! is the longwave radiative heating rate
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
   integer :: iradae   ! not used by RRTMG

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
      dosw = .true. 
      dolw = .true.
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
   use time_manager,       only: get_step_size, is_first_restart_step
   use radiation_data,     only: init_rad_data

   ! RRTMGP modules
   use mo_rrtmgp_clr_all_sky, only: &
         rrtmgp_sw_init=>rte_sw_init, &
         rrtmgp_lw_init=>rte_lw_init
   use mo_load_coefficients, only: rrtmgp_load_coefficients=>load_and_init
   
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
   character(len=16) :: SPCAM_microp_scheme  ! SPCAM microphysics scheme
   !-----------------------------------------------------------------------
    
   call init_rad_data() ! initialize output fields for offline driver

   call phys_getopts(microp_scheme_out = microp_scheme)

   ! Read k-distribution coefficients from file and initialize k-distribution
   ! objects.
   call rrtmgp_load_coefficients(k_dist_sw, coefficients_file_sw)
   call rrtmgp_load_coefficients(k_dist_lw, coefficients_file_lw

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

   call phys_getopts(use_SPCAM_out = use_SPCAM)
   call phys_getopts(SPCAM_microp_scheme_out = SPCAM_microp_scheme)
   if (use_SPCAM .and. SPCAM_microp_scheme .eq. 'sam1mom') then
      cldfsnow_idx = 0
   end if

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

   use perf_mod,        only: t_startf, t_stopf
   use physics_buffer, only: physics_buffer_desc, pbuf_get_field, &
                             pbuf_old_tim_idx, pbuf_get_index, pbuf_set_field
   use phys_grid,       only: get_rlat_all_p, get_rlon_all_p
   use physics_types,   only: physics_state, physics_ptend
   use cospsimulator_intr, only: docosp, cospsimulator_intr_run,cosp_nradsteps
   use time_manager,    only: get_curr_calday
   use camsrfexch,      only: cam_out_t, cam_in_t
   use cam_history_support, only: fillvalue
   use parrrtm,         only: nbndlw
   use parrrsw,         only: nbndsw
   use hirsbt,          only: hirsrtm
   use hirsbtpar,       only: pnb_hirs, pnf_msu, hirsname, msuname
   use radheat,         only: radheat_tend
   use ppgrid
   use pspect
   use rad_constituents, only: oldcldoptics, liqcldoptics, icecldoptics
   use aer_rad_props,    only: aer_rad_props_sw, aer_rad_props_lw
   use interpolate_data, only: vertinterp
   use radiation_data,   only: output_rad_data

   use crmdims,              only: crm_nx, crm_ny, crm_nz, crm_nx_rad, crm_ny_rad
   use physconst,            only: cpair, cappa
   use constituents,         only: cnst_get_ind
#ifdef CRM
   use crm_physics,          only: m2005_effradius
#endif
#ifdef MODAL_AERO
  use modal_aero_data,       only: ntot_amode
#endif
   use phys_control,         only: phys_getopts
   use orbit,                only: zenith
   use output_aerocom_aie,   only: do_aerocom_ind3
   use pkg_cldoptics,        only: cldefr

   ! TODO: implement RRTMGP counterparts to these
   !use rrtmg_state,      only: rrtmg_state_create, rrtmg_state_update, &
   !                            rrtmg_state_destroy, rrtmg_state_t, num_rrtmg_levs
   use mo_rrtmgp_clr_all_sky, only: &
         rrtmgp_sw=>rte_sw, rrtmgp_lw=>rte_lw
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

   ! Combined cloud radiative parameters (liquid, ice, and possibly snow). Note
   ! that these are "in cloud" not "in cell"
   real(r8) :: c_cld_tau    (nbndsw,pcols,pver) ! cloud extinction optical depth
   real(r8) :: c_cld_tau_w  (nbndsw,pcols,pver) ! cloud single scattering albedo * tau
   real(r8) :: c_cld_tau_w_g(nbndsw,pcols,pver) ! cloud assymetry parameter * w * tau
   real(r8) :: c_cld_tau_w_f(nbndsw,pcols,pver) ! cloud forward scattered fraction * w * tau
   real(r8) :: c_cld_lw_abs (nbndlw,pcols,pver) ! cloud absorption optics depth (LW)

   ! Cloud radiative parameters (liquid and ice cloud). Note that these are are
   ! "in cloud" not "in cell"
   real(r8) :: cld_tau    (nbndsw,pcols,pver) ! cloud extinction optical depth
   real(r8) :: cld_tau_w  (nbndsw,pcols,pver) ! cloud single scattering albedo * tau
   real(r8) :: cld_tau_w_g(nbndsw,pcols,pver) ! cloud assymetry parameter * w * tau
   real(r8) :: cld_tau_w_f(nbndsw,pcols,pver) ! cloud forward scattered fraction * w * tau
   real(r8) :: cld_lw_abs (nbndlw,pcols,pver) ! cloud absorption optics depth (LW)

   ! Ice cloud radiative parameters. Note that these are "in cloud" not "in cell".
   real(r8) :: ice_tau    (nbndsw,pcols,pver) ! ice extinction optical depth
   real(r8) :: ice_tau_w  (nbndsw,pcols,pver) ! ice single scattering albedo * tau
   real(r8) :: ice_tau_w_g(nbndsw,pcols,pver) ! ice assymetry parameter * tau * w
   real(r8) :: ice_tau_w_f(nbndsw,pcols,pver) ! ice forward scattered fraction * tau * w
   real(r8) :: ice_lw_abs (nbndlw,pcols,pver) ! ice absorption optics depth (LW)

   ! "Snow" cloud radiative parameters. Note that these are "in cloud" not "in cell".
   real(r8) :: snow_tau    (nbndsw,pcols,pver) ! snow extinction optical depth
   real(r8) :: snow_tau_w  (nbndsw,pcols,pver) ! snow single scattering albedo * tau
   real(r8) :: snow_tau_w_g(nbndsw,pcols,pver) ! snow assymetry parameter * tau * w
   real(r8) :: snow_tau_w_f(nbndsw,pcols,pver) ! snow forward scattered fraction * tau * w
   real(r8) :: snow_lw_abs (nbndlw,pcols,pver) ! snow absorption optics depth (LW)

   ! Liquid cloud radiative parameters. Note that these are "in cloud" not "in cell".
   real(r8) :: liq_tau    (nbndsw,pcols,pver) ! liquid extinction optical depth
   real(r8) :: liq_tau_w  (nbndsw,pcols,pver) ! liquid single scattering albedo * tau
   real(r8) :: liq_tau_w_g(nbndsw,pcols,pver) ! liquid assymetry parameter * tau * w
   real(r8) :: liq_tau_w_f(nbndsw,pcols,pver) ! liquid forward scattered fraction * tau * w
   real(r8) :: liq_lw_abs (nbndlw,pcols,pver) ! liquid absorption optics depth (LW)

   ! Aerosol radiative properties
   real(r8) :: aer_tau    (pcols,0:pver,nbndsw) ! aerosol extinction optical depth
   real(r8) :: aer_tau_w  (pcols,0:pver,nbndsw) ! aerosol single scattering albedo * tau
   real(r8) :: aer_tau_w_g(pcols,0:pver,nbndsw) ! aerosol assymetry parameter * w * tau
   real(r8) :: aer_tau_w_f(pcols,0:pver,nbndsw) ! aerosol forward scattered fraction * w * tau
   real(r8) :: aer_lw_abs (pcols,pver,nbndlw)   ! aerosol absorption optics depth (LW)

   ! Stuff for output fields on history files
   ! TODO: move these to a subroutine
   real(r8) :: tot_cld_vistau(pcols,pver)  ! tot gbx cloud visible sw optical depth for output on history files
   real(r8) :: tot_icld_vistau(pcols,pver) ! tot in-cloud visible sw optical depth for output on history files
   real(r8) :: liq_icld_vistau(pcols,pver) ! liq in-cloud visible sw optical depth for output on history files
   real(r8) :: ice_icld_vistau(pcols,pver) ! ice in-cloud visible sw optical depth for output on history files
   real(r8) :: snow_icld_vistau(pcols,pver) ! snow in-cloud visible sw optical depth for output on history files

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
   
   ! pointers to heating rates
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
   ! TODO: this seems naive on RRTMGP's part. Should it not just return a
   ! fillvalue for non-day columns?
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

   dosw = radiation_do('sw')
   dolw = radiation_do('lw')
   if (dosw .or. dolw) then

      ! initialize RRTMGP state
      ! TODO

      ! Do shortwave stuff...
      if (dosw) then
         ! Do shortwave cloud optics calculations
         call t_startf('shortwave cloud optics')
         call set_cloud_optics_sw(state, pbuf, cloud_sw)
         call t_stopf('shortwave cloud optics')

         ! Get shortwave gas optics
         call t_startf('shortwave gas concentrations')
         call set_gas_optics_sw(state, pbuf, gas_concentrations_sw)
         call t_stopf('shortwave gas concentrations')

         ! Get shortwave aerosol optics
         call t_startf('shortwave aerosol optics')
         call set_aerosol_optics_sw(state, pbuf, aerosol_sw)
         call t_stopf('shortwave aerosol optics')

         ! Subset optical properties to get only daytime columns
         call subset_daytime_optics()

         ! Do shortwave radiative transfer calculations
         call t_startf('shortwave radiation calculations')
         errmsg = rrtmgp_sw(k_dist_sw, gas_concentrations_sw, &
                            pmid_day, t_day, pint_day, &
                            coszrs_day, alb_dir, alb_dif, cloud_sw, &
                            allsky_fluxes_sw, clrsky_fluxes_sw, &
                            aer_props=aerosol_sw)
         call t_stopf('shortwave radiation calculations')

         ! Map RRTMGP outputs to CAM outputs
         ! TODO
         
      end if  ! dosw

      ! Do longwave stuff...
      if (dolw) then

         ! Do longwave cloud optics calculations
         call t_startf('longwave cloud optics')
         call set_cloud_optics_lw(state, pbuf, cloud_lw)
         call t_stopf('longwave cloud optics')

         ! Get longwave gas optics
         call t_startf('longwave gas concentrations')
         call set_gas_optics_lw(state, pbuf, gas_concentrations_lw)
         call t_stopf('longwave gas concentrations')

         ! Get longwave aerosol optics
         call t_startf('longwave aerosol optics')
         call set_aerosol_optics_lw(state, pbuf, aerosol_lw)
         call t_stopf('longwave aerosol optics')

         ! Do longwave radiative transfer calculations
         call t_startf('longwave radiation calculations')
         errmsg = rrtmgp_lw(kdist_lw, gas_concentrations_lw, &
                            pmid_rad, t_day, pint_rad, &
                            t_sfc, emis_sfc, cloud_lw, flw, flwc, &
                            aer_props=aerosol_lw)
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
   call radheat_tend(state, pbuf,  ptend, qrl, qrs, fsns, &
                     fsnt, flns, flnt, cam_in%asdir, net_flx)
   call t_stopf('radheat_tend')

   ! Compute net heating rate for dtheta/dt
   call t_startf('heating_rate')
   do k=1,pver
      do i = 1,state%ncol
         heating_rate(i,k) = (qrs(i,k) + qrl(i,k))/cpair * (1.e5_r8/state%pmid(i,k))**cappa
      end do
   end do
   call t_stopf('heating_rate')

   ! convert radiative heating rates to Q*dp for energy conservation
   if (conserve_energy) then
      do k = 1,pver
         do i = 1,state%ncol
            qrs(i,k) = qrs(i,k)*state%pdel(i,k)
            qrl(i,k) = qrl(i,k)*state%pdel(i,k)
         end do
      end do
   end if

   ! copy net sw flux to cam_out structure
   cam_out%netsw(:state%ncol) = fsns(:state%ncol)

end subroutine radiation_tend


subroutine set_cloud_optics_sw(state, pbuf, cloud_sw)
   
   use physics_types, only: physics_state
   use physics_buffer, only: physics_buffer_desc
   use mo_optical_props, only: ty_optical_props_2str

   type(physics_state), intent(in) :: state
   type(physics_buffer_desc), intent(in) :: pbuf
   type(ty_optical_props_2str), intent(out) :: cloud_sw

   if (use_SPCAM) then
      ! get cloud optics from resolved fields
      call get_spcam_optics(state, pbuf, cloud_sw)
   else
      call get_optics_sw(state, pbuf, cld_tau, cld_tau_w, &
                         cld_tau_w_g, cld_tau_w_f)

      ! Do MCICA sampling
      call do_mcica_cloud_sampling(cld_tau, cld_tau_w, cld_tau_w_g, &
                                   cld_tau_w_f, cloud_sw)
   end if

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
   real(r8) cld_crm    (pcols, crm_nx_rad, crm_ny_rad, crm_nz)
   real(r8) cliqwp_crm (pcols, crm_nx_rad, crm_ny_rad, crm_nz)
   real(r8) cicewp_crm (pcols, crm_nx_rad, crm_ny_rad, crm_nz)
   real(r8) rel_crm    (pcols, crm_nx_rad, crm_ny_rad, crm_nz)
   real(r8) rei_crm    (pcols, crm_nx_rad, crm_ny_rad, crm_nz)
   real(r8) cld_tau_crm(pcols, crm_nx_rad, crm_ny_rad, crm_nz)
   real(r8) emis_crm   (pcols, crm_nx_rad, crm_ny_rad, crm_nz)
   real(r8) qrl_crm    (pcols, crm_nx_rad, crm_ny_rad, crm_nz)
   real(r8) qrs_crm    (pcols, crm_nx_rad, crm_ny_rad, crm_nz)
   real(r8) crm_fsnt   (pcols, crm_nx_rad, crm_ny_rad)   ! net shortwave fluxes at TOA at CRM grids
   real(r8) crm_fsntc  (pcols, crm_nx_rad, crm_ny_rad)   ! net clear-sky shortwave fluxes at TOA at CRM grids
   real(r8) crm_fsns   (pcols, crm_nx_rad, crm_ny_rad)   ! net shortwave fluxes at surface at CRM grids
   real(r8) crm_fsnsc  (pcols, crm_nx_rad, crm_ny_rad)   ! net clear-sky shortwave fluxes at surface at CRM grids
   real(r8) crm_flnt   (pcols, crm_nx_rad, crm_ny_rad)   ! net longwave fluxes at TOA at CRM grids
   real(r8) crm_flntc  (pcols, crm_nx_rad, crm_ny_rad)   ! net clear-sky longwave fluxes at TOA at CRM grids
   real(r8) crm_flns   (pcols, crm_nx_rad, crm_ny_rad)   ! net longwave fluxes at surface at CRM grids
   real(r8) crm_flnsc  (pcols, crm_nx_rad, crm_ny_rad)   ! net clear-sky longwave fluxes at surface at CRM grids

   real(r8) crm_aodvisz(pcols, crm_nx_rad, crm_ny_rad, crm_nz)   ! layer aerosol optical depth at 550nm at CRM grids
   real(r8) crm_aodvis (pcols, crm_nx_rad, crm_ny_rad)   ! AOD at 550nm at CRM grids
   real(r8) crm_aod400 (pcols, crm_nx_rad, crm_ny_rad)   ! AOD at 400nm at CRM grids
   real(r8) crm_aod700 (pcols, crm_nx_rad, crm_ny_rad)   ! AOD at 700nm at CRM grids
   real(r8) aod400     (pcols)   ! AOD at 400nm at CRM grids
   real(r8) aod700     (pcols)   ! AOD at 700nm at CRM grids

   integer :: nct_tot_icld_vistau(pcols,pver) ! the number of CRM columns that has in-cloud visible sw optical depth 
   integer :: nct_liq_icld_vistau(pcols,pver) ! the number of CRM column that has liq in-cloud visible sw optical depth 
   integer :: nct_ice_icld_vistau(pcols,pver) ! the number of CRM column that has ice in-cloud visible sw optical depth 
   integer :: nct_snow_icld_vistau(pcols,pver) ! the number of CRM column that has snow in-cloud visible sw optical depth 
 
end subroutine output_radiation_diagnostics()

!===============================================================================

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


end module radiation

