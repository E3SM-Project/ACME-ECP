module tphysbc_sp_mod

   use shr_kind_mod,       only: r8 => shr_kind_r8
   use shr_sys_mod,        only: shr_sys_flush
   use spmd_utils,         only: masterproc, iam
   use physconst,          only: latvap, latice, rh2o
   use physics_types,      only: physics_state, physics_tend, physics_state_set_grid,  &
                                 physics_ptend, physics_tend_init,                     &
                                 physics_type_alloc, physics_ptend_dealloc,            &
                                 physics_state_alloc, physics_state_dealloc,           &
                                 physics_tend_alloc, physics_tend_dealloc
   use physics_update_mod, only: physics_update, physics_update_init, hist_vars,       &
                                 nvars_prtrb_hist, get_var
   use phys_grid,          only: get_ncols_p
   use phys_gmean,         only: gmean_mass
   use ppgrid,             only: begchunk, endchunk, pcols, pver, pverp, psubcols
   use constituents,       only: pcnst, cnst_name, cnst_get_ind
   use camsrfexch,         only: cam_out_t, cam_in_t
   use cam_control_mod,    only: ideal_phys, adiabatic
   use phys_control,       only: phys_do_flux_avg, phys_getopts, waccmx_is
   use zm_conv,            only: trigmem
   use scamMod,            only: single_column, scm_crm_mode
   use flux_avg,           only: flux_avg_init
#ifdef SPMD
   use mpishorthand
#endif
   use perf_mod
   use cam_logfile,        only: iulog
   use camsrfexch,         only: cam_export
   
   use modal_aero_calcsize,   only: modal_aero_calcsize_init, modal_aero_calcsize_diag, &
                                    modal_aero_calcsize_reg, modal_aero_calcsize_sub
   use modal_aero_wateruptake,only: modal_aero_wateruptake_init, modal_aero_wateruptake_dr, &
                                    modal_aero_wateruptake_reg

   implicit none
   private 

   !!!  Physics buffer indices
   integer :: teout_idx            = 0  
   integer :: tini_idx             = 0 
   integer :: qini_idx             = 0 
   integer :: cldliqini_idx        = 0 
   integer :: cldiceini_idx        = 0 
   integer :: static_ener_ac_idx   = 0
   integer :: water_vap_ac_idx     = 0
   ! integer :: species_class(pcnst) = -1 

   save

   public :: tphysbc_sp

   !!! Physics package options
   character(len=16) :: shallow_scheme
   character(len=16) :: macrop_scheme
   character(len=16) :: microp_scheme 
   integer           :: cld_macmic_num_steps    ! Number of macro/micro substeps
   logical           :: do_clubb_sgs
   logical           :: use_subcol_microp   ! if true, use subcolumns in microphysics
   logical           :: state_debug_checks  ! Debug physics_state.
   logical           :: clim_modal_aero     ! climate controled by prognostic or prescribed modal aerosols
   logical           :: prog_modal_aero     ! Prognostic modal aerosols present
   logical           :: micro_do_icesupersat
   logical           :: pergro_test_active= .false.
   logical           :: pergro_mods = .false.
   logical           :: is_cmip6_volc !true if cmip6 style volcanic file is read otherwise false
   
   contains

#if defined( SP_ALT_TPHYSBC )

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
subroutine tphysbc_sp(ztodt, pbuf2d, landm_in, &
                      fsns_in, fsnt_in, flns_in, flnt_in, fsds_in,   & 
                      state_in, tend_in, cam_in_in, cam_out_in,      &
                      species_class )
   !---------------------------------------------------------------------------
   ! Purpose: 
   ! Evaluate and apply physical processes that occur 
   ! BEFORE coupling to land, sea, and ice models.  
   !---------------------------------------------------------------------------
   use physics_buffer,     only: physics_buffer_desc, pbuf_get_field, pbuf_get_chunk
   use physics_buffer,     only: pbuf_get_index, pbuf_old_tim_idx
   use physics_buffer,     only: col_type_subcol, dyn_time_lvls
   use shr_kind_mod,       only: r8 => shr_kind_r8
   use stratiform,         only: stratiform_tend
   use microp_driver,      only: microp_driver_tend
   use microp_aero,        only: microp_aero_run
   use macrop_driver,      only: macrop_driver_tend
   use physics_types,      only: physics_state, physics_tend, physics_ptend, physics_ptend_init, &
                                 physics_ptend_sum, physics_state_check, physics_ptend_scale
   use cam_diagnostics,    only: diag_conv_tend_ini, diag_phys_writeout, diag_conv, diag_export, diag_state_b4_phys_write
   use cam_history,        only: outfld, fieldname_len
   use physconst,          only: cpair, latvap, gravit, rga
   use constituents,       only: pcnst, qmin, cnst_get_ind
   use convect_deep,       only: convect_deep_tend, convect_deep_tend_2, deep_scheme_does_scav_trans
   use time_manager,       only: is_first_step, get_nstep
   use convect_shallow,    only: convect_shallow_tend
   use check_energy,       only: check_energy_chng, check_energy_fix, check_qflx,  & 
                                 check_water, check_prect, check_energy_timestep_init
   use check_energy,       only: check_tracers_data, check_tracers_init, check_tracers_chng
   use dycore,             only: dycore_is
   use aero_model,         only: aero_model_wetdep
   use radiation,          only: radiation_tend
   use cloud_diagnostics,  only: cloud_diagnostics_calc
   use perf_mod
   use mo_gas_phase_chemdr,only: map2chm
   use clybry_fam,         only: clybry_fam_adj
   use clubb_intr,         only: clubb_tend_cam
   use sslt_rebin,         only: sslt_rebin_adv
   use tropopause,         only: tropopause_output
   use output_aerocom_aie, only: do_aerocom_ind3, cloud_top_aerocom
   use cam_abortutils,     only: endrun
   use subcol,             only: subcol_gen, subcol_ptend_avg
   use subcol_utils,       only: subcol_ptend_copy, is_subcol_on
   use phys_control,       only: use_qqflx_fixer, use_mass_borrower

#ifdef CRM
   !!! CRM modules
   use crmdims,               only: crm_nz, crm_nx, crm_ny, crm_dx, crm_dy, crm_dt
   use crm_physics,           only: crm_physics_tend, crm_surface_flux_bypass_tend, &
                                    crm_save_state_tend, crm_recall_state_tend
   use crm_ecpp_output_module,only: crm_ecpp_output_type

#if defined( ECPP )
   use module_ecpp_ppdriver2, only: parampollu_driver2
   use module_data_ecpp1,     only: dtstep_pp_input
   use crmclouds_camaerosols, only: crmclouds_mixnuc_tend
#endif

#if defined( DIFFUSE_PHYS_TEND )
   use phys_hyperviscosity_mod
#endif

#endif /* CRM */

   implicit none
   !---------------------------------------------------------------------------
   !!! Interface arguments
   real(r8),            intent(in   ) :: ztodt                                ! physics time step
   type(physics_buffer_desc), pointer :: pbuf2d   (:,:)                       ! physics buffer
   real(r8),            intent(in   ), target :: landm_in (pcols,begchunk:endchunk)   ! land fraction ramp
   real(r8),            intent(inout), target :: fsns_in  (pcols,begchunk:endchunk)   ! Surface solar absorbed flux
   real(r8),            intent(inout), target :: fsnt_in  (pcols,begchunk:endchunk)   ! Net column abs solar flux at model top
   real(r8),            intent(inout), target :: flns_in  (pcols,begchunk:endchunk)   ! Srf longwave cooling (up-down) flux
   real(r8),            intent(inout), target :: flnt_in  (pcols,begchunk:endchunk)   ! Net outgoing lw flux at model top
   real(r8),            intent(inout), target :: fsds_in  (pcols,begchunk:endchunk)   ! Surface solar down flux
   type(physics_state), intent(inout), target :: state_in       (begchunk:endchunk)   ! physics state
   type(physics_tend ), intent(inout), target :: tend_in        (begchunk:endchunk)   ! physics tend (for dynamics)
   type(cam_in_t),      intent(in   ), target :: cam_in_in      (begchunk:endchunk)   ! atmos comp input structure
   type(cam_out_t),     intent(inout), target :: cam_out_in     (begchunk:endchunk)   ! atmos comp output structure
   integer,             intent(inout), target :: species_class(pcnst)
   !---------------------------------------------------------------------------
   !!! Local variables
   type(physics_ptend) :: ptend              ! indivdual parameterization tendencies
   type(physics_ptend) :: ptend_crm(begchunk:endchunk)    ! CRM tendencies

   !!! temporary variables for given chunk
   type(physics_buffer_desc), pointer :: pbuf  (:) 
   ! real(r8)             :: landm (pcols)
   ! real(r8)             :: fsns  (pcols)
   ! real(r8)             :: fsnt  (pcols)
   ! real(r8)             :: flns  (pcols)
   ! real(r8)             :: flnt  (pcols)
   ! real(r8)             :: fsds  (pcols)
   ! type(physics_state)  :: state
   ! type(physics_tend)   :: tend
   ! type(cam_in_t)       :: cam_in
   ! type(cam_out_t)      :: cam_out
   real(r8)            ,pointer :: landm (:)
   real(r8)            ,pointer :: fsns  (:)
   real(r8)            ,pointer :: fsnt  (:)
   real(r8)            ,pointer :: flns  (:)
   real(r8)            ,pointer :: flnt  (:)
   real(r8)            ,pointer :: fsds  (:)
   type(physics_state) ,pointer :: state
   type(physics_tend)  ,pointer :: tend
   type(cam_in_t)      ,pointer :: cam_in
   type(cam_out_t)     ,pointer :: cam_out

   integer lchnk                             ! chunk identifier
   integer ncol                              ! number of atmospheric columns
   integer  :: nstep                         ! current timestep number
   real(r8) :: net_flx(pcols)                ! radiation flux
   real(r8) :: rtdt                          ! 1./ztodt
   integer  :: c,i,k,m                       ! loop variables
   integer  :: ixcldice, ixcldliq            ! constituent indices for cloud liquid and ice water.

   !!! physics buffer fields to compute tendencies for stratiform package
   integer itim_old, ifld
   real(r8), pointer, dimension(:,:) :: cld        ! cloud fraction
   real(r8), pointer, dimension(:,:) :: cldo       ! old cloud fraction

   character(len=16)  :: deep_scheme      ! Default set in phys_control.F90

   !!! physics buffer fields for total energy and mass adjustment
   real(r8), pointer, dimension(:  )   :: teout
   real(r8), pointer, dimension(:,:)   :: tini
   real(r8), pointer, dimension(:,:)   :: qini
   real(r8), pointer, dimension(:,:)   :: cldliqini
   real(r8), pointer, dimension(:,:)   :: cldiceini
   real(r8), pointer, dimension(:,:)   :: dtcore

   real(r8), pointer, dimension(:,:,:) :: fracis   ! fraction of transported species that are insoluble

   !!! convective precipitation variables
   real(r8),pointer :: prec_dp(:)                  ! total precipitation from ZM convection
   real(r8),pointer :: snow_dp(:)                  ! snow from ZM convection
   real(r8),pointer :: prec_sh(:)                  ! total precipitation from Hack convection
   real(r8),pointer :: snow_sh(:)                  ! snow from Hack convection

   !!! stratiform precipitation variables
   real(r8),pointer :: prec_str(:)                 ! sfc flux of precip from stratiform (m/s)
   real(r8),pointer :: snow_str(:)                 ! sfc flux of snow from stratiform   (m/s)
   real(r8),pointer :: prec_str_sc(:)              ! sfc flux of precip from stratiform (m/s) -- for subcolumns
   real(r8),pointer :: snow_str_sc(:)              ! sfc flux of snow from stratiform   (m/s) -- for subcolumns
   real(r8),pointer :: prec_pcw(:)                 ! total precip from prognostic cloud scheme
   real(r8),pointer :: snow_pcw(:)                 ! snow from prognostic cloud scheme
   real(r8),pointer :: prec_sed(:)                 ! total precip from cloud sedimentation
   real(r8),pointer :: snow_sed(:)                 ! snow from cloud ice sedimentation
   real(r8)         :: sh_e_ed_ratio(pcols,pver)   ! shallow conv [ent/(ent+det)] ratio  

   !!! energy checking variables
   real(r8) :: zero(pcols)                         ! array of zeros
   real(r8) :: zero_tracers(pcols,pcnst)           ! array of zeros
   real(r8) :: flx_heat(pcols)
   type(check_tracers_data) :: tracerint           ! energy integrals and cummulative boundary fluxes
   logical  :: lq(pcnst)
   real(r8) :: ftem(pcols,pver)                    ! tmp space
   real(r8), pointer, dimension(:) :: static_ener_ac_2d ! Vertically integrated static energy
   real(r8), pointer, dimension(:) :: water_vap_ac_2d   ! Vertically integrated water vapor
   real(r8) :: CIDiff(pcols)                            ! Difference in vertically integrated static energy

   logical :: l_bc_energy_fix
   logical :: l_dry_adj
   logical :: l_tracer_aero
   logical :: l_st_mac
   logical :: l_st_mic
   logical :: l_rad

   real(r8) :: qexcess(pcols)

   logical                    :: use_SPCAM
   logical                    :: use_ECPP
   character(len=16)          :: SPCAM_microp_scheme
#ifdef CRM
   integer                    :: phys_stage       ! physics stage indicator (tphysbc => 1)
   real(r8)                   :: crm_run_time     ! length of CRM integration
   real(r8), dimension(pcols) :: sp_qchk_prec_dp  ! CRM precipitation diagostic (liq+ice)  used for check_energy_chng
   real(r8), dimension(pcols) :: sp_qchk_snow_dp  ! CRM precipitation diagostic (ice only) used for check_energy_chng
   real(r8), dimension(pcols) :: sp_rad_flux      ! CRM radiative flux diagnostic used for check_energy_chng
   type(crm_ecpp_output_type) :: crm_ecpp_output  ! CRM output data for ECPP calculations
#if defined( ECPP )
   !!! ECPP variables
   real(r8),pointer,dimension(:)   :: pblh              ! PBL height (for ECPP)
   real(r8),pointer,dimension(:,:) :: acldy_cen_tbeg    ! cloud fraction
   real(r8)                        :: dtstep_pp         ! ECPP time step (seconds)
   integer                         :: necpp             ! number of GCM time steps in which ECPP is called once
#endif /* ECPP */
#endif /* CRM */

   !---------------------------------------------------------------------------
   !---------------------------------------------------------------------------
   call phys_getopts( use_SPCAM_out           = use_SPCAM )
   call phys_getopts( use_ECPP_out            = use_ECPP)
   call phys_getopts( SPCAM_microp_scheme_out = SPCAM_microp_scheme)
   call phys_getopts( microp_scheme_out       = microp_scheme      &
                     ,macrop_scheme_out       = macrop_scheme      &
                     ,use_subcol_microp_out   = use_subcol_microp  &
                     ,deep_scheme_out         = deep_scheme        &
                     ,state_debug_checks_out  = state_debug_checks &
                     ,l_bc_energy_fix_out     = l_bc_energy_fix    &
                     ,l_dry_adj_out           = l_dry_adj          &
                     ,l_tracer_aero_out       = l_tracer_aero      &
                     ,l_st_mac_out            = l_st_mac           &
                     ,l_st_mic_out            = l_st_mic           &
                     ,l_rad_out               = l_rad              &
                     )
   
   teout_idx     = pbuf_get_index('TEOUT')
   tini_idx      = pbuf_get_index('TINI')
   qini_idx      = pbuf_get_index('QINI')
   cldliqini_idx = pbuf_get_index('CLDLIQINI')
   cldiceini_idx = pbuf_get_index('CLDICEINI')

   !---------------------------------------------------------------------------
   ! Initialize ptend,tend,state for CRM
   !---------------------------------------------------------------------------
   lq(:) = .true.
   do c = begchunk, endchunk
      call physics_ptend_init(ptend_crm(c), state_in(c)%psetcols, 'crm', lu=.true., lv=.true., ls=.true., lq=lq, fromcrm=.true.) 
   end do 

   ! call physics_tend_alloc(tend, pcols)

   !---------------------------------------------------------------------------
   ! start loop over chunks
   !---------------------------------------------------------------------------
   do c = begchunk, endchunk
      !---------------------------------------------------------------------------
      ! Set up temporary variables for given chunk 
      !---------------------------------------------------------------------------

      lchnk = state_in(c)%lchnk
      ncol  = state_in(c)%ncol

      !!! need to allocate state here because we need the correct lchnk
      ! call physics_state_alloc(state, lchnk, pcols)

      ! landm   = landm_in(:,c) 
      ! fsns    = fsns_in (:,c) 
      ! fsnt    = fsnt_in (:,c) 
      ! flns    = flns_in (:,c) 
      ! flnt    = flnt_in (:,c) 
      ! fsds    = fsds_in (:,c) 
      ! state   = state_in  (c)
      ! tend    = tend_in   (c)
      ! cam_in  = cam_in_in (c)
      ! cam_out = cam_out_in(c)

      landm   => landm_in(:,c) 
      fsns    => fsns_in (:,c) 
      fsnt    => fsnt_in (:,c) 
      flns    => flns_in (:,c) 
      flnt    => flnt_in (:,c) 
      fsds    => fsds_in (:,c) 
      state   => state_in  (c)
      tend    => tend_in   (c)
      cam_in  => cam_in_in (c)
      cam_out => cam_out_in(c)

      pbuf => pbuf_get_chunk(pbuf2d,c)

      !---------------------------------------------------------------------------
      !---------------------------------------------------------------------------
      call t_startf('tphysbc_sp_init')

      zero = 0._r8
      zero_tracers(:,:) = 0._r8

      rtdt = 1._r8/ztodt

      nstep = get_nstep()

      !---------------------------------------------------------------------------
      !---------------------------------------------------------------------------

      static_ener_ac_idx = pbuf_get_index('static_ener_ac')
      call pbuf_get_field(pbuf, static_ener_ac_idx, static_ener_ac_2d )
      water_vap_ac_idx   = pbuf_get_index('water_vap_ac')
      call pbuf_get_field(pbuf, water_vap_ac_idx, water_vap_ac_2d )

      ! Integrate and compute the difference
      ! CIDiff = difference of column integrated values
      if( nstep == 0 ) then
         CIDiff(:ncol) = 0.0_r8
         call outfld('DTENDTH', CIDiff, pcols, lchnk )
         call outfld('DTENDTQ', CIDiff, pcols, lchnk )
      else
         ! MSE first
         ftem(:ncol,:) = (state%s(:ncol,:) + latvap*state%q(:ncol,:,1)) * state%pdel(:ncol,:)*rga
         do k=2,pver
            ftem(:ncol,1) = ftem(:ncol,1) + ftem(:ncol,k)
         end do
         CIDiff(:ncol) = (ftem(:ncol,1) - static_ener_ac_2d(:ncol))*rtdt

         call outfld('DTENDTH', CIDiff, pcols, lchnk )
         ! Water vapor second
         ftem(:ncol,:) = state%q(:ncol,:,1)*state%pdel(:ncol,:)*rga
         do k=2,pver
            ftem(:ncol,1) = ftem(:ncol,1) + ftem(:ncol,k)
         end do
         CIDiff(:ncol) = (ftem(:ncol,1) - water_vap_ac_2d(:ncol))*rtdt

         call outfld('DTENDTQ', CIDiff, pcols, lchnk )
      end if

      !---------------------------------------------------------------------------
      ! Associate pointers with physics buffer fields
      !---------------------------------------------------------------------------
      itim_old = pbuf_old_tim_idx()
      ifld = pbuf_get_index('CLD')
      call pbuf_get_field(pbuf, ifld, cld, (/1,1,itim_old/),(/pcols,pver,1/))

      call pbuf_get_field(pbuf, teout_idx, teout, (/1,itim_old/), (/pcols,1/))

      call pbuf_get_field(pbuf, tini_idx, tini)
      call pbuf_get_field(pbuf, qini_idx, qini)
      call pbuf_get_field(pbuf, cldliqini_idx, cldliqini)
      call pbuf_get_field(pbuf, cldiceini_idx, cldiceini)

      ifld =  pbuf_get_index('DTCORE')
      call pbuf_get_field(pbuf, ifld, dtcore, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

      ifld = pbuf_get_index('FRACIS')
      call pbuf_get_field(pbuf, ifld, fracis, start=(/1,1,1/), kount=(/pcols, pver, pcnst/)  )
      fracis (:ncol,:,1:pcnst) = 1._r8

      ! Set physics tendencies to 0
      tend %dTdt(:ncol,:pver) = 0._r8
      tend %dudt(:ncol,:pver) = 0._r8
      tend %dvdt(:ncol,:pver) = 0._r8

      call check_qflx (state, tend, "PHYBC01", nstep, ztodt, cam_in%cflx(:,1))
      call check_water(state, tend, "PHYBC01", nstep, ztodt)

      !---------------------------------------------------------------------------
      !---------------------------------------------------------------------------
#if defined(SP_FLUX_BYPASS)
      if(.not.use_qqflx_fixer) then 
         ! Check if latent heat flux exceeds the total moisture content of the
         ! lowest model layer, thereby creating negative moisture.
         call qneg4('TPHYSBC '       ,lchnk               ,ncol  ,ztodt ,        &
                     state%q(1,pver,1),state%rpdel(1,pver) ,cam_in%shf ,         &
                     cam_in%lhf , cam_in%cflx ,qexcess)
      end if 
      call outfld('QEXCESS',qexcess,pcols,lchnk)
#endif
      !---------------------------------------------------------------------------
      !---------------------------------------------------------------------------

      if(use_mass_borrower) then 
         !!! printout diagnostic information
          call qneg3('TPHYSBCb',lchnk  ,ncol    ,pcols   ,pver    , &
               1, pcnst, qmin  ,state%q, .False.)
         !!! tracer borrower for mass conservation 
         do m = 1, pcnst 
            call massborrow("PHYBC01",lchnk,ncol,state%psetcols,m,m,qmin(m),state%q(1,1,m),state%pdel)
         end do
      else
         !!! original fixer to make sure tracers are all positive
         call qneg3('TPHYSBCb',lchnk  ,ncol    ,pcols   ,pver    , &
                    1, pcnst, qmin  ,state%q, .True. )
      end if ! use_mass_borrower

      !!! Validate state coming from the dynamics.
      if (state_debug_checks) &
         call physics_state_check(state, name="before tphysbc (dycore?)")

      call clybry_fam_adj( ncol, lchnk, map2chm, state%q, pbuf )

      if(use_mass_borrower) then
         !!! if use_mass_borrower = True, only printout diagnostic information
         call qneg3('TPHYSBCc',lchnk  ,ncol    ,pcols   ,pver    , &
                    1, pcnst, qmin  ,state%q, .False. )
         !!! tracer borrower for mass conservation 
         do m = 1, pcnst
            call massborrow("PHYBC02",lchnk,ncol,state%psetcols,m,m,qmin(m),state%q(1,1,m),state%pdel)
         end do
      else
         !!! original fixer to make sure tracers are all positive
         call qneg3('TPHYSBCc',lchnk  ,ncol    ,pcols   ,pver    , &
                    1, pcnst, qmin  ,state%q, .True. )
      end if ! use_mass_borrower

      call check_water(state, tend, "PHYBC02", nstep, ztodt)

      !!! Validate output of clybry_fam_adj.
      if (state_debug_checks) &
         call physics_state_check(state, name="clybry_fam_adj")

      !!! Dump out "before physics" state
      call diag_state_b4_phys_write(state)

      !!! compute mass integrals of input tracers state
      call check_tracers_init(state, tracerint)

      call t_stopf('tphysbc_sp_init')

      !---------------------------------------------------------------------------
      ! Global mean total energy fixer
      !---------------------------------------------------------------------------
      if (l_bc_energy_fix) then

         call t_startf('energy_fixer')

         !!!*** BAB's FV heating kludge *** save the initial temperature
         tini(:ncol,:pver) = state%t(:ncol,:pver)
         if (dycore_is('LR') .or. dycore_is('SE'))  then
            call check_energy_fix(state, ptend, nstep, flx_heat)
            call physics_update(state, ptend, ztodt, tend)
            call check_energy_chng(state, tend, "chkengyfix", nstep, ztodt, zero, zero, zero, flx_heat)
         end if
         !!! Save state for convective tendency calculations.
         call diag_conv_tend_ini(state, pbuf)

         call cnst_get_ind('CLDLIQ', ixcldliq)
         call cnst_get_ind('CLDICE', ixcldice)
         qini     (:ncol,:pver) = state%q(:ncol,:pver,       1)
         cldliqini(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq)
         cldiceini(:ncol,:pver) = state%q(:ncol,:pver,ixcldice)

         call outfld('TEOUT', teout       , pcols, lchnk   )
         call outfld('TEINP', state%te_ini, pcols, lchnk   )
         call outfld('TEFIX', state%te_cur, pcols, lchnk   )

         !!! set and output the dse change due to dynpkg
         if( nstep > dyn_time_lvls-1 ) then
            do k = 1,pver
               dtcore(:ncol,k) = (tini(:ncol,k) - dtcore(:ncol,k))/(ztodt) + tend%dTdt(:ncol,k)
            end do
            call outfld( 'DTCORE', dtcore, pcols, lchnk )
         end if

         call t_stopf('energy_fixer')

      end if
      
      !---------------------------------------------------------------------------
      ! Dry adjustment
      !---------------------------------------------------------------------------
      if (l_dry_adj) then

         call t_startf('dry_adjustment')

         ! Copy state info for input to dadadj
         ! This is a kludge, so that dadadj does not have to be correctly reformulated in dry static energy

         lq(:) = .FALSE.
         lq(1) = .TRUE.
         call physics_ptend_init(ptend, state%psetcols, 'dadadj', ls=.true., lq=lq)
         ptend%s(:ncol,:pver)   = state%t(:ncol,:pver)
         ptend%q(:ncol,:pver,1) = state%q(:ncol,:pver,1)

         call dadadj (lchnk, ncol, state%pmid,  &
                      state%pint,  state%pdel,  &
                      ptend%s, ptend%q(1,1,1) )
         ptend%s(:ncol,:)   = (ptend%s(:ncol,:)   - state%t(:ncol,:)  )/ztodt * cpair
         ptend%q(:ncol,:,1) = (ptend%q(:ncol,:,1) - state%q(:ncol,:,1))/ztodt
         call physics_update(state, ptend, ztodt, tend)

         call t_stopf('dry_adjustment')

      end if

      !======================================================================================
      !--------------------------------------------------------------------------------------
      ! CRM Physics
      !--------------------------------------------------------------------------------------
      !======================================================================================
#ifdef CRM

      crm_run_time = ztodt
      !---------------------------------------------------------------------------
      ! Apply surface fluxes if using SP_FLUX_BYPASS
      !---------------------------------------------------------------------------
#if defined( SP_FLUX_BYPASS )
      call crm_surface_flux_bypass_tend(state, cam_in, ptend)
      call physics_update(state, ptend, ztodt, tend)  
      call check_energy_chng(state, tend, "crm_tend", nstep, crm_run_time,  &
                             cam_in%shf(:), zero, zero, cam_in%cflx(:,1)) 
#endif
      !---------------------------------------------------------------------------
      ! Initialize variabale for ECPP data
      !---------------------------------------------------------------------------
#if defined( ECPP )
      if (use_ECPP) then
         call crm_ecpp_output%initialize(pcols,pver)
      end if ! use_ECPP
#endif
      !---------------------------------------------------------------------------
      ! Run the CRM 
      !---------------------------------------------------------------------------
      phys_stage = 1  ! for tphysbc() => phys_stage = 1
      call crm_physics_tend(ztodt, state, tend, ptend, pbuf, cam_in, cam_out,    &
                            species_class, phys_stage, crm_ecpp_output,         &
                            sp_qchk_prec_dp, sp_qchk_snow_dp, sp_rad_flux)

      ptend_crm(c) = ptend
      !---------------------------------------------------------------------------
      !---------------------------------------------------------------------------
      !!! this deallocation is necessary
      ! call physics_state_dealloc(state)
      !---------------------------------------------------------------------------
      !---------------------------------------------------------------------------
   end do ! c = begchunk, endchunk


#if defined( DIFFUSE_PHYS_TEND )
   call phys_hyperviscosity(ptend_crm)
#endif
      

   !!! restart loop over chunks
   do c = begchunk, endchunk
      !---------------------------------------------------------------------------
      ! Set up temporary variables for given chunk 
      !---------------------------------------------------------------------------
      lchnk = state_in(c)%lchnk
      ncol  = state_in(c)%ncol

      !!! need to re-allocate state here because we need the correct lchnk
      ! call physics_state_alloc(state, lchnk, pcols)

      landm   => landm_in(:,c) 
      fsns    => fsns_in (:,c) 
      fsnt    => fsnt_in (:,c) 
      flns    => flns_in (:,c) 
      flnt    => flnt_in (:,c) 
      fsds    => fsds_in (:,c) 
      state   => state_in  (c)
      tend    => tend_in   (c)
      cam_in  => cam_in_in (c)
      cam_out => cam_out_in(c)

      pbuf => pbuf_get_chunk(pbuf2d,c)

      ptend = ptend_crm(c)

      !---------------------------------------------------------------------------
      ! update state with smoothed tendencies
      !---------------------------------------------------------------------------

      call physics_update(state, ptend, crm_run_time, tend)

      call check_energy_chng(state, tend, "crm_tend", nstep, crm_run_time,  &
                             zero, sp_qchk_prec_dp, sp_qchk_snow_dp, sp_rad_flux)

      !---------------------------------------------------------------------------
      ! Modal aerosol wet radius for radiative calculation
      !---------------------------------------------------------------------------
#if defined(MODAL_AERO)  
      !!! temporarily turn on all lq, so it is allocated
      lq(:) = .true.
      call physics_ptend_init(ptend, state%psetcols, 'crm - modal_aero_wateruptake_dr', lq=lq)

      !!! set all ptend%lq to false as they will be set in modal_aero_calcsize_sub
      ptend%lq(:) = .false.
      call modal_aero_calcsize_sub (state, ptend, ztodt, pbuf)
      call modal_aero_wateruptake_dr(state, pbuf)

      ! When ECPP is included, wet deposition is done ECPP,
      ! So tendency from wet depostion is not updated 
      ! in mz_aero_wet_intr (mz_aerosols_intr.F90)
      ! tendency from other parts of crmclouds_aerosol_wet_intr() are still updated here.
      call physics_update (state, ptend, crm_run_time, tend)
         call check_energy_chng(state, tend, "crm_tend", nstep, crm_run_time, zero, zero, zero, zero)
#endif /* MODAL_AERO */

      !---------------------------------------------------------------------------
      ! ECPP - Explicit-Cloud Parameterized-Pollutant
      !---------------------------------------------------------------------------
#if defined( ECPP )
      if (use_ECPP) then

         call pbuf_get_field(pbuf, pbuf_get_index('pblh'), pblh)
         call pbuf_get_field(pbuf, pbuf_get_index('ACLDY_CEN'), acldy_cen_tbeg)
       
         dtstep_pp = dtstep_pp_input
         necpp = dtstep_pp/crm_run_time

         if (nstep.ne.0 .and. mod(nstep, necpp).eq.0) then

            !!! aerosol tendency from droplet activation and mixing
            !!! cldo and cldn are set to be the same in crmclouds_mixnuc_tend,
            !!! So only turbulence mixing is done here.
            call t_startf('crmclouds_mixnuc')
            call crmclouds_mixnuc_tend(state, ptend, dtstep_pp,           &
                                       cam_in%cflx, pblh, pbuf,           &
                                       crm_ecpp_output%wwqui_cen,         &
                                       crm_ecpp_output%wwqui_cloudy_cen,  &
                                       crm_ecpp_output%wwqui_bnd,         &
                                       crm_ecpp_output%wwqui_cloudy_bnd,  &
                                       species_class)
            call physics_update(state, ptend, dtstep_pp, tend)
            call t_stopf('crmclouds_mixnuc')

            !!! ECPP interface
            call t_startf('ecpp')
            call parampollu_driver2(state, ptend, pbuf, dtstep_pp, dtstep_pp,   &
                                    crm_ecpp_output%acen,       crm_ecpp_output%abnd,         &
                                    crm_ecpp_output%acen_tf,    crm_ecpp_output%abnd_tf,      &
                                    crm_ecpp_output%massflxbnd, crm_ecpp_output%rhcen,        &
                                    crm_ecpp_output%qcloudcen,  crm_ecpp_output%qlsinkcen,    &
                                    crm_ecpp_output%precrcen,   crm_ecpp_output%precsolidcen, &
                                    acldy_cen_tbeg)
            call physics_update(state, ptend, dtstep_pp, tend)
            call t_stopf ('ecpp')

         end if ! nstep.ne.0 .and. mod(nstep, necpp).eq.0

         call crm_ecpp_output%finalize()

      end if ! use_ECPP
#endif /* ECPP */

      !---------------------------------------------------------------------------
      ! save old CRM cloud fraction - w/o CRM, this is done in cldwat2m.F90
      !---------------------------------------------------------------------------

      ifld = pbuf_get_index('CLDO')
      call pbuf_get_field(pbuf, ifld, cldo, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

      ifld = pbuf_get_index('CLD')
      call pbuf_get_field(pbuf, ifld, cld , start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

      cldo(1:ncol,1:pver) = cld(1:ncol,1:pver)

      !---------------------------------------------------------------------------
      !---------------------------------------------------------------------------
#endif /* CRM */
      !======================================================================================
      !--------------------------------------------------------------------------------------
      ! end CRM Physics
      !--------------------------------------------------------------------------------------
      !======================================================================================


      !---------------------------------------------------------------------------
      ! Moist physical parameteriztions complete, send data to history file
      !---------------------------------------------------------------------------
      call t_startf('bc_history_write')
      call diag_phys_writeout(state, cam_out%psl)
      call diag_conv(state, ztodt, pbuf)
      call t_stopf('bc_history_write')

      !---------------------------------------------------------------------------
      ! Write cloud diagnostics on history file
      !---------------------------------------------------------------------------
      call t_startf('bc_cld_diag_history_write')
      call cloud_diagnostics_calc(state, pbuf)
      call t_stopf('bc_cld_diag_history_write')

      !---------------------------------------------------------------------------
      ! Radiation computations
      !---------------------------------------------------------------------------
      if (l_rad) then
         
         call t_startf('radiation')
         call radiation_tend( state,ptend, pbuf,                  &
                              cam_out, cam_in,                    &
                              cam_in%landfrac,landm,              &
                              cam_in%icefrac, cam_in%snowhland,   &
                              fsns, fsnt, flns, flnt,             &
                              fsds, net_flx,is_cmip6_volc)

         ! Set net flux used by spectral dycores
         do i=1,ncol
            tend%flx_net(i) = net_flx(i)
         end do

         !!! for super-parameterization, don't add radiative tendency 
         !!! because it was added above as part of crm tendency.
         if (use_SPCAM) ptend%s = 0.

         call physics_update(state, ptend, ztodt, tend)

         if (use_SPCAM) then
            call check_energy_chng(state, tend, "spradheat", nstep, ztodt, zero, zero, zero, zero) 
         else 
            call check_energy_chng(state, tend, "radheat", nstep, ztodt, zero, zero, zero, net_flx)
         endif

         call t_stopf('radiation')

      end if ! l_rad
      !---------------------------------------------------------------------------
      !---------------------------------------------------------------------------

      if(do_aerocom_ind3) then
         call cloud_top_aerocom(state, pbuf) 
      end if

      !!! Diagnose the location of the tropopause and its location to the history file(s).
      call t_startf('tropopause')
      call tropopause_output(state)
      call t_stopf('tropopause')

      !!! Save atmospheric fields to force surface models
      call t_startf('cam_export')
      call cam_export (state,cam_out,pbuf)
      call t_stopf('cam_export')

      !!! Write export state to history file
      call t_startf('diag_export')
      call diag_export(cam_out)
      call t_stopf('diag_export')

      !---------------------------------------------------------------------------
      !---------------------------------------------------------------------------
   end do ! c = begchunk, endchunk

   ! call physics_state_dealloc(state)
   ! call physics_tend_dealloc(tend)
   call physics_ptend_dealloc(ptend)

end subroutine tphysbc_sp
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
#endif /* SP_ALT_TPHYSBC */

end module tphysbc_sp_mod