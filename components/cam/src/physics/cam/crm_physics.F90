#define SP_DIR_NS

module crm_physics
!-----------------------------------------------------------------------
! Purpose: 
! 
!    Provides the CAM interface to the crm code.  
!
! Revision history: 
! June, 2009, Minghuai Wang:  
!          crm_physics_tend 
! July, 2009, Minghuai Wang: m2005_effradius
!
!---------------------------------------------------------------------------
   use shr_kind_mod,    only: r8 => shr_kind_r8
   use cam_abortutils,  only: endrun
   use cam_logfile,     only: iulog
   use physics_types,   only: physics_state, physics_tend
   use ppgrid,          only: pcols, pver, pverp
   use constituents,    only: pcnst
#ifdef MODAL_AERO
   use modal_aero_data, only: ntot_amode
#endif

   implicit none 
   private
   save

   public :: crm_physics_tend, crm_physics_register, crm_physics_init
   public :: crm_save_state_tend
   public :: m2005_effradius

   integer :: crm_u_idx, crm_v_idx, crm_w_idx, crm_t_idx
   integer :: crm_qt_idx, crm_nc_idx, crm_qr_idx, crm_nr_idx, crm_qi_idx, crm_ni_idx 
   integer :: crm_qs_idx, crm_ns_idx, crm_qg_idx, crm_ng_idx, crm_qc_idx, crm_qp_idx, crm_qn_idx
   integer :: crm_t_rad_idx, crm_qv_rad_idx, crm_qc_rad_idx, crm_qi_rad_idx, crm_cld_rad_idx
   integer :: crm_nc_rad_idx, crm_ni_rad_idx, crm_qs_rad_idx, crm_ns_rad_idx, crm_qrad_idx
   integer :: crm_qaerwat_idx, crm_dgnumwet_idx
   integer :: prec_dp_idx, snow_dp_idx, prec_sh_idx, snow_sh_idx
   integer :: prec_sed_idx, snow_sed_idx, prec_pcw_idx, snow_pcw_idx
   integer :: cldo_idx

   integer :: clubb_buffer_idx

   type(physics_state)  :: state_save
   type(physics_tend)   :: tend_save
   real(r8) cldo_save(pcols, pver)  ! save for the old cloud fraction
   real(r8),pointer     :: acldy_cen_tbeg(:,:)        ! cloud fraction
   real(r8), pointer, dimension(:,:) :: cldo

#ifdef MODAL_AERO
   real(r8)          :: qqcw_save(pcols,pver,pcnst)
   real(r8)          :: qqcw_all(pcols,pver,pcnst)
   real(r8)          :: dgnumwet_save(pcols, pver, ntot_amode)
   real(r8),pointer  :: dgnumwet(:,:,:)
#endif

!========================================================================================================
contains
!========================================================================================================

!---------------------------------------------------------------------------------------------------------
subroutine crm_physics_register()
!-------------------------------------------------------------------------------------------------------
! 
! Purpose:  add necessary fileds into physics buffer
!
!--------------------------------------------------------------------------------------------------------
  use spmd_utils,      only: masterproc
  use physconst,       only: mwdry, cpair
  use ppgrid,          only: pcols, pver, pverp
  use physics_buffer,  only: dyn_time_lvls, pbuf_add_field, dtype_r8, pbuf_get_index
  use phys_control,    only: phys_getopts
  use crmdims,         only: crm_nx, crm_ny, crm_nz, crm_dx, crm_dy, crm_dt, nclubbvars
#ifdef CRM
  use microphysics,    only: nmicro_fields
  use setparm_mod,     only: setparm
#endif



! local variables
  integer idx
  logical           :: use_ECPP, use_SPCAM
  character(len=16) :: SPCAM_microp_scheme

  call phys_getopts( use_ECPP_out            = use_ECPP)
  call phys_getopts( SPCAM_microp_scheme_out = SPCAM_microp_scheme)
  call phys_getopts( use_SPCAM_out           = use_SPCAM)

  if(masterproc) then
      print*,'_________________________________________'
      print*,'_ Super-parameterization run ____________'
      print*,'crm_nx=',crm_nx,'   crm_ny=',crm_ny,'   crm_nz=',crm_nz
      print*,'crm_dx=',crm_dx,'   crm_dy=',crm_dy,'   crm_dt=',crm_dt
      if (SPCAM_microp_scheme .eq. 'sam1mom') print*,'Microphysics: SAM1MOM'
      if (SPCAM_microp_scheme .eq. 'm2005') print*,'Microphysics: M2005'
      print*,'_________________________________________'
  end if

#ifdef CLUBB_CRM
  call pbuf_add_field('CLUBB_BUFFER','global', dtype_r8, (/pcols,crm_nx,crm_ny,crm_nz+1,nclubbvars/), clubb_buffer_idx)
#endif

#ifdef CRM 
  if (use_SPCAM) then
     call setparm()
  end if

  call pbuf_add_field('CRM_U',     'global', dtype_r8, (/pcols,crm_nx, crm_ny, crm_nz/),                crm_u_idx)
  call pbuf_add_field('CRM_V',     'global', dtype_r8, (/pcols,crm_nx, crm_ny, crm_nz/),                crm_v_idx)
  call pbuf_add_field('CRM_W',     'global', dtype_r8, (/pcols,crm_nx, crm_ny, crm_nz/),                crm_w_idx)
  call pbuf_add_field('CRM_T',     'global', dtype_r8, (/pcols,crm_nx, crm_ny, crm_nz/),                crm_t_idx)

  call pbuf_add_field('CRM_T_RAD',   'physpkg', dtype_r8, (/pcols,crm_nx, crm_ny, crm_nz/),              crm_t_rad_idx)
  call pbuf_add_field('CRM_QV_RAD',  'physpkg', dtype_r8, (/pcols,crm_nx, crm_ny, crm_nz/),              crm_qv_rad_idx)
  call pbuf_add_field('CRM_QC_RAD',  'physpkg', dtype_r8, (/pcols,crm_nx, crm_ny, crm_nz/),              crm_qc_rad_idx)
  call pbuf_add_field('CRM_QI_RAD',  'physpkg', dtype_r8, (/pcols,crm_nx, crm_ny, crm_nz/),              crm_qi_rad_idx)
  call pbuf_add_field('CRM_CLD_RAD', 'physpkg', dtype_r8, (/pcols,crm_nx, crm_ny, crm_nz/),              crm_cld_rad_idx)
  call pbuf_add_field('CRM_QRAD',    'global',  dtype_r8, (/pcols,crm_nx, crm_ny, crm_nz/),              crm_qrad_idx)


#ifdef MODAL_AERO
  call pbuf_add_field('CRM_QAERWAT', 'physpkg', dtype_r8, (/pcols,crm_nx, crm_ny, crm_nz, ntot_amode/),  crm_qaerwat_idx)
  call pbuf_add_field('CRM_DGNUMWET','physpkg', dtype_r8, (/pcols,crm_nx, crm_ny, crm_nz, ntot_amode/),  crm_dgnumwet_idx)
#endif


  if (SPCAM_microp_scheme .eq. 'm2005') then
    call pbuf_add_field('CRM_NC_RAD','physpkg', dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),               crm_nc_rad_idx)
    call pbuf_add_field('CRM_NI_RAD','physpkg', dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),               crm_ni_rad_idx)
    call pbuf_add_field('CRM_QS_RAD','physpkg', dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),               crm_qs_rad_idx)
    call pbuf_add_field('CRM_NS_RAD','physpkg', dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),               crm_ns_rad_idx)
    call pbuf_add_field('CRM_QT',    'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),               crm_qt_idx)
    call pbuf_add_field('CRM_NC',    'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),               crm_nc_idx)
    call pbuf_add_field('CRM_QR',    'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),               crm_qr_idx)
    call pbuf_add_field('CRM_NR',    'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),               crm_nr_idx)
    call pbuf_add_field('CRM_QI',    'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),               crm_qi_idx)
    call pbuf_add_field('CRM_NI',    'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),               crm_ni_idx)
    call pbuf_add_field('CRM_QS',    'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),               crm_qs_idx)
    call pbuf_add_field('CRM_NS',    'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),               crm_ns_idx)
    call pbuf_add_field('CRM_QG',    'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),               crm_qg_idx)
    call pbuf_add_field('CRM_NG',    'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),               crm_ng_idx)
    call pbuf_add_field('CRM_QC',    'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),               crm_qc_idx)
    cldo_idx = pbuf_get_index('CLDO')
  else
    call pbuf_add_field('CRM_QT',    'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),               crm_qt_idx)
    call pbuf_add_field('CRM_QP',    'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),               crm_qp_idx)
    call pbuf_add_field('CRM_QN',    'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),               crm_qn_idx)
    call pbuf_add_field('CLDO',      'global',  dtype_r8, (/pcols ,pver, dyn_time_lvls/),                      cldo_idx)
  endif
#endif

   
  if (use_ECPP) then
     call pbuf_add_field('MU_CRM', 'physpkg', dtype_r8, (/pcols,pver/), idx)  ! mass flux up
     call pbuf_add_field('MD_CRM', 'physpkg', dtype_r8, (/pcols,pver/), idx)  ! mass flux down
     call pbuf_add_field('DU_CRM', 'physpkg', dtype_r8, (/pcols,pver/), idx)  ! mass detrainment from updraft
     call pbuf_add_field('EU_CRM', 'physpkg', dtype_r8, (/pcols,pver/), idx)  ! mass detrainment from updraft
     call pbuf_add_field('ED_CRM', 'physpkg', dtype_r8, (/pcols,pver/), idx)  ! mass detrainment from downdraft
     call pbuf_add_field('JT_CRM', 'physpkg', dtype_r8, (/pcols/), idx)  ! index of cloud (convection) top for each column
     call pbuf_add_field('MX_CRM', 'physpkg', dtype_r8, (/pcols/), idx)  !  index of cloud (convection) bottom for each column
     call pbuf_add_field('IDEEP_CRM', 'physpkg', dtype_r8, (/pcols/), idx)  ! Gathering array for convective columns

     call pbuf_add_field('TKE_CRM','physpkg', dtype_r8, (/pcols,pver/), idx) ! TKE from CRM  (m2/s2)
     call pbuf_add_field('TK_CRM', 'physpkg', dtype_r8, (/pcols,pver/), idx) ! TK from CRM (m2/s)
!
! ACLDY_CEN has to be gloabl in the physcal buffer to be saved in the restart file??
  call pbuf_add_field('ACLDY_CEN','global', dtype_r8, (/pcols,pver/), idx) ! total (all sub-classes) cloudy fractional area in previous time step 
  endif 


end subroutine crm_physics_register
!=========================================================================================================

subroutine crm_physics_init()
!-------------------------------------------------------------------------------------------------------
! 
! Purpose: initialize some varialbes, and add necessary fileds into output fields 
!
!--------------------------------------------------------------------------------------------------------
  use physics_buffer,  only: pbuf_get_index
  use physics_types,   only: physics_tend_alloc
  use physconst,       only: mwdry, cpair
  use ppgrid,          only: pcols, pver, pverp
  use constituents,    only: pcnst, cnst_name
  use cam_history,     only: addfld, add_default, horiz_only
  use crmdims,         only: crm_nx, crm_ny, crm_nz
  use phys_control,    only: phys_getopts
#ifdef ECPP
  use module_ecpp_ppdriver2,   only: papampollu_init
  use ecppvars,   only: NCLASS_CL,ncls_ecpp_in,NCLASS_PR
#endif
#ifdef MODAL_AERO
    use cam_history,   only: fieldname_len
    use spmd_utils,    only: masterproc
    use modal_aero_data, only:  cnst_name_cw, &
!                               lmassptr_amode, lmassptrcw_amode, lwaterptr_amode, &
                                lmassptr_amode, lmassptrcw_amode, &
                                nspec_amode, ntot_amode, numptr_amode, numptrcw_amode, ntot_amode, &
                                species_class, spec_class_gas
       
    integer :: l, lphase, lspec
    character(len=fieldname_len)   :: tmpname
    character(len=fieldname_len+3) :: fieldname
    character(128)                 :: long_name
    character(8)                   :: unit
#endif


! local variables
  integer :: m
  logical :: use_ECPP
  character(len=16) :: SPCAM_microp_scheme

  call phys_getopts(use_ECPP_out = use_ECPP)
  call phys_getopts(SPCAM_microp_scheme_out = SPCAM_microp_scheme)

#ifdef ECPP
  if (use_ECPP) then
     call papampollu_init ()
  end if
#endif

  if (SPCAM_microp_scheme .eq. 'm2005') then
     call addfld ('SPNC    ',(/ 'lev' /), 'A', '/kg   ','Cloud water dropet number from CRM')
     call addfld ('SPNI    ',(/ 'lev' /), 'A', '/kg   ','Cloud ice crystal number from CRM')
     call addfld ('SPNS    ',(/ 'lev' /), 'A', '/kg   ','Snow particle number from CRM')
     call addfld ('SPNG    ',(/ 'lev' /), 'A', '/kg   ','Graupel particle number from CRM')
     call addfld ('SPNR    ',(/ 'lev' /), 'A', '/kg   ','Rain particle number from CRM')
     call add_default ('SPNC    ', 1, ' ')
     call add_default ('SPNI    ', 1, ' ')
     call add_default ('SPNS    ', 1, ' ')
     call add_default ('SPNG    ', 1, ' ')
     call add_default ('SPNR    ', 1, ' ')

     call addfld ('CRM_FLIQ ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '1      ','Frequency of Occurrence of Liquid'      )
     call addfld ('CRM_FICE ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '1      ','Frequency of Occurrence of Ice'         )
     call addfld ('CRM_FRAIN',(/'crm_nx','crm_ny','crm_nz'/), 'A', '1      ','Frequency of Occurrence of Rain'        )
     call addfld ('CRM_FSNOW',(/'crm_nx','crm_ny','crm_nz'/), 'A', '1      ','Frequency of Occurrence of Snow'        )
     call addfld ('CRM_FGRAP',(/'crm_nx','crm_ny','crm_nz'/), 'A', '1      ','Frequency of Occurrence of Graupel'     )

     call addfld ('CRM_QS  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', 'kg/kg   ','Snow mixing ratio from CRM'             )
     call addfld ('CRM_QG  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', 'kg/kg   ','Graupel mixing ratio from CRM'          )
     call addfld ('CRM_QR  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', 'kg/kg   ','Rain mixing ratio from CRM'             )

     call addfld ('CRM_NC  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/kg     ','Cloud water dropet number from CRM'     )
     call addfld ('CRM_NI  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/kg     ','Cloud ice crystal number from CRM'      )
     call addfld ('CRM_NS  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/kg     ','Snow particle number from CRM'          )
     call addfld ('CRM_NG  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/kg     ','Graupel particle number from CRM'       )
     call addfld ('CRM_NR  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/kg     ','Rain particle number from CRM'          )

! hm 7/26/11, add new output
! below is for *instantaneous* crm output
     call addfld ('CRM_AUT  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/s     ','Autoconversion cloud waterfrom CRM'     )
     call addfld ('CRM_ACC  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/s     ','Accretion cloud water from CRM'         )
     call addfld ('CRM_EVPC ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/s     ','Evaporation cloud water from CRM'       )
     call addfld ('CRM_EVPR ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/s     ','Evaporation rain from CRM'              )
     call addfld ('CRM_MLT  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/s     ','Melting ice snow graupel from CRM'      )
     call addfld ('CRM_SUB  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/s     ','Sublimation ice snow graupel from CRM'  )
     call addfld ('CRM_DEP  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/s     ','Deposition ice snow graupel from CRM'   )
     call addfld ('CRM_CON  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/s     ','Condensation cloud water from CRM'      )

! hm 8/31/11
! below is for *gcm-grid and time-step-avg* process output
     call addfld ('A_AUT  ',(/ 'lev' /), 'A', '/s   ','Avg autoconversion cloud water from CRM'            )
     call addfld ('A_ACC  ',(/ 'lev' /), 'A', '/s   ','Avg accretion cloud water from CRM'                 )
     call addfld ('A_EVPC ',(/ 'lev' /), 'A', '/s   ','Avg evaporation cloud water from CRM'               )
     call addfld ('A_EVPR ',(/ 'lev' /), 'A', '/s   ','Avg evaporation rain from CRM'                      )
     call addfld ('A_MLT  ',(/ 'lev' /), 'A', '/s   ','Avg melting ice snow graupel from CRM'              )
     call addfld ('A_SUB  ',(/ 'lev' /), 'A', '/s   ','Avg sublimation ice snow graupel from CRM'          )
     call addfld ('A_DEP  ',(/ 'lev' /), 'A', '/s   ','Avg deposition ice snow graupel from CRM'           )
     call addfld ('A_CON  ',(/ 'lev' /), 'A', '/s   ','Avg condensation cloud water from CRM'              )

     call addfld ('CRM_REL  ', (/'crm_nx','crm_ny','crm_nz'/), 'A', 'micrometers', 'cloud scale droplet effective radius')
     call addfld ('CRM_REI  ', (/'crm_nx','crm_ny','crm_nz'/), 'A', 'micrometers', 'cloud scale ice crystal effective radius')
     call addfld ('CRM_DEI  ', (/'crm_nx','crm_ny','crm_nz'/), 'A', 'micrometers', 'cloud scale Mitchell ice effective diameter')
     call addfld ('CRM_DES  ', (/'crm_nx','crm_ny','crm_nz'/), 'A', 'micrometers', 'cloud scale snow effective diameter')
     call addfld ('CRM_MU   ', (/'crm_nx','crm_ny','crm_nz'/), 'A', 'micrometers', &
                                                 'cloud scale droplet size distribution shape parameter for radiation')
     call addfld ('CRM_LAMBDA',(/'crm_nx','crm_ny','crm_nz'/), 'A', 'micrometers',  &
                                                 'cloud scale slope of droplet distribution for radiation')
     call addfld ('CRM_TAU  ', (/'crm_nx','crm_ny','crm_nz'/), 'A', '1',           'cloud scale cloud optical depth'  )
     call addfld ('CRM_WVAR' , (/'crm_nx','crm_ny','crm_nz'/), 'A', 'm/s',         'vertical velocity variance from CRM')

     call addfld ('CRM_FSNT',  (/'crm_nx','crm_ny'/),          'A',  'unitless', 'net TOA shortwave fluxes at CRM grids')
     call addfld ('CRM_FSNTC', (/'crm_nx','crm_ny'/),          'A',  'unitless', 'net TOA clear-sky shortwave fluxes at CRM grids')
     call addfld ('CRM_FSNS',  (/'crm_nx','crm_ny'/),          'A',  'unitless', 'net surface shortwave fluxes at CRM grids')
     call addfld ('CRM_FSNSC', (/'crm_nx','crm_ny'/),          'A',  'unitless',  &
                                                 'net surface clear-sky shortwave fluxes at CRM grids')
     call addfld ('CRM_FLNT',  (/'crm_nx','crm_ny'/),          'A',  'unitless', 'net TOA longwave fluxes at CRM grids')
     call addfld ('CRM_FLNTC', (/'crm_nx','crm_ny'/),          'A',  'unitless', 'net TOA clear-sky longwave fluxes at CRM grids')
     call addfld ('CRM_FLNS',  (/'crm_nx','crm_ny'/),          'A',  'unitless', 'net surface longwave fluxes at CRM grids')
     call addfld ('CRM_FLNSC', (/'crm_nx','crm_ny'/),          'A',  'unitless',  &
                                                 'net surface clear-sky longwave fluxes at CRM grids')

     call addfld ('CRM_AODVIS', (/'crm_nx','crm_ny'/),          'A', 'unitless', 'Aerosol optical depth at 550nm in CRM grids',&
                                                                                                               flag_xyfill=.true.)
     call addfld ('CRM_AOD400', (/'crm_nx','crm_ny'/),          'A', 'unitless', 'Aerosol optical depth at 400nm in CRM grids',&
                                                                                                               flag_xyfill=.true.)
     call addfld ('CRM_AOD700', (/'crm_nx','crm_ny'/),          'A', 'unitless', 'Aerosol optical depth at 700nm in CRM grids', &
                                                                                                               flag_xyfill=.true.)
     call addfld ('CRM_AODVISZ',(/'crm_nx','crm_ny','crm_nz'/), 'A', 'unitless',  &
                                                  'Aerosol optical depth at each layer at 500nm in CRM grids', flag_xyfill=.true.)
     call addfld ('AOD400',      horiz_only,                    'A', 'unitless', 'Aerosol optical depth at 400nm', &
                                                                                                               flag_xyfill=.true.)
     call addfld ('AOD700',      horiz_only,                    'A', 'unitless', 'Aerosol optical depth at 700nm', &
                                                                                                               flag_xyfill=.true.)
     call add_default ('AOD400',  1, ' ')
     call add_default ('AOD700',  1, ' ')
   endif
#ifdef CLUBB_CRM
     call addfld ('UP2     ',    (/'crm_nx','crm_ny','crm_z1'/), 'A', 'm^2/s^2',    'u prime ^2 from clubb')
     call addfld ('VP2     ',    (/'crm_nx','crm_ny','crm_z1'/), 'A', 'm^2/s^2',    'v prime ^2 from clubb')
     call addfld ('WPRTP   ',    (/'crm_nx','crm_ny','crm_z1'/), 'A', 'mkg/skg',    'w prime * rt prime from clubb')
     call addfld ('WPTHLP  ',    (/'crm_nx','crm_ny','crm_z1'/), 'A', 'mK/s',       'w prime * th_l prime from clubb')
     call addfld ('WP2     ',    (/'crm_nx','crm_ny','crm_z1'/), 'A', 'm^2/s^2',    'w prime ^2 from clubb')
     call addfld ('WP3     ',    (/'crm_nx','crm_ny','crm_z1'/), 'A', 'm^3/s^3',    'w prime ^3 from clubb')
     call addfld ('RTP2    ',    (/'crm_nx','crm_ny','crm_z1'/), 'A', '(kg/kg)2',   'r_t prime ^2 from clubb')
     call addfld ('THLP2   ',    (/'crm_nx','crm_ny','crm_z1'/), 'A', 'K^2',        'th_l_prime ^2 from clubb')
     call addfld ('RTPTHLP ',    (/'crm_nx','crm_ny','crm_z1'/), 'A', 'kgK/kg',     'r_t prime * th_l prime  from clubb')
     call addfld ('UPWP    ',    (/'crm_nx','crm_ny','crm_z1'/), 'A', 'm^2/s^2',    'u prime * w prime from clubb')
     call addfld ('VPWP    ',    (/'crm_nx','crm_ny','crm_z1'/), 'A', 'm^2/s^2',    'v prime * w prime from clubb')
     call addfld ('CRM_CLD ',    (/'crm_nx','crm_ny','crm_z1'/), 'A', 'fraction',   'cloud fraction from clubb')
#endif

#ifdef ECPP
  if (use_ECPP) then
     call addfld ('ABND    ', (/'ilev        ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'fraction', &
                  'cloud fraction for each sub-sub class for full time period at layer boundary')
     call addfld ('ABND_TF ', (/'ilev        ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'fraction', &
                  'cloud fraction for each sub-sub class for end-portion of time period at layer boundary')
     call addfld ('MASFBND ', (/'ilev        ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'kg/m2/s',  &
                  'sub-class vertical mass flux (kg/m2/s) at layer boundary')
     call addfld ('ACEN    ', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'fraction', &
                  'cloud fraction for each sub-sub class for full time period at layer center')
     call addfld ('ACEN_TF ', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'fraction', &
                  'cloud fraction for each sub-sub class for end-portion of time period at layer center')
     call addfld ('RHCEN   ', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'fraction', &
                  'relative humidity for each sub-sub calss at layer center')
     call addfld ('QCCEN   ', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'kg/kg',    &
                  'cloud water for each sub-sub class at layer center')
     call addfld ('QICEN   ', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'kg/kg',    &
                  'cloud ice for each sub-sub class at layer center')
     call addfld ('QSINK_AFCEN', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', '/s',    &
                  'cloud water loss rate from precip. using cloud water after precip. for each sub-sub class at layer center')
     call addfld ('QSINK_BFCEN', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', '/s',    &
                 'cloud water loss rate from precip. using cloud water before precip. for each sub-sub class at layer center')
     call addfld ('QSINK_AVGCEN', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', '/s',   &
      'cloud water loss rate from precip. using averaged cloud water and precip. rate for each sub-sub class at layer center')
     call addfld ('PRAINCEN', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'kg/kg/s',  &
                  ' cloud water loss rate from precipitation (kg/kg/s) for each sub-sub class at layer center')
     call addfld ('PRECRCEN', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'kg/m2/s',  &
                  'liquid (rain) precipitation rate for each sub-sub class at layer center')
     call addfld ('PRECSCEN', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'kg/m2/s',  &
                  'solid (snow, graupel,...) precipitation rate for each sub-sub class at layer center')
     call addfld ('WUPTHRES',      (/ 'ilev' /), 'A', 'm/s',    'vertical velocity threshold for updraft')
     call addfld ('WDNTHRES',      (/ 'ilev' /), 'A', 'm/s',    'vertical velocity threshold for dndraft')
     call addfld ('WWQUI_CEN',     (/ 'lev' /),  'A', 'm2/s2',  'vertical velocity variance in the quiescent class, layer center')
     call addfld ('WWQUI_CLD_CEN', (/ 'lev' /),  'A', 'm2/s2',  &
                                                      'vertical velocity variance in the cloudy quiescent class, layer center')
     call addfld ('WWQUI_BND',     (/ 'ilev' /), 'A', 'm2/s2',  &
                                                      'vertical velocity variance in the quiescent class, layer boundary')
     call addfld ('WWQUI_CLD_BND', (/ 'ilev' /), 'A', 'm2/s2',  &
                                                      'vertical velocity variance in the cloudy quiescent class, layer boundary')
  endif
#endif

  call addfld ('MU_CRM   ', (/ 'lev' /), 'A', 'Pa/s',     'mass flux up from CRM')
  call addfld ('MD_CRM   ', (/ 'lev' /), 'A', 'Pa/s',     'mass flux down from CRM')
  call addfld ('DU_CRM   ', (/ 'lev' /), 'A', '/s',       'detrainment from updraft from CRM')
  call addfld ('EU_CRM   ', (/ 'lev' /), 'A', '/s',       'entraiment rate from updraft')
  call addfld ('ED_CRM   ', (/ 'lev' /), 'A', '/s',       'entraiment rate from downdraft')

  do m=1, pcnst 
    if(cnst_name(m) == 'DMS') then 
       call addfld('DMSCONV',   (/ 'lev' /), 'A', 'kg/kg/s',  'DMS tendency from ZM convection')
    end if
    if(cnst_name(m) == 'SO2') then 
       call addfld('SO2CONV',   (/ 'lev' /), 'A', 'kg/kg/s',  'SO2 tendency from ZM convection')
     end if
  end do

  call addfld ('SPQRL    ', (/ 'lev' /), 'A', 'K/s',      'long-wave heating rate')
  call addfld ('SPQRS    ', (/ 'lev' /), 'A', 'K/s',      'short-wave heating rate')
  call addfld ('LENGC    ', (/ 'ilev' /), 'A', 'm  ',      'Mixing length scale for the calcuation of vertical difusivity')

  call addfld ('SPKVH     ',(/ 'ilev' /), 'A', 'm2/s    ', 'Vertical diffusivity used in dropmixnuc in the MMF call')
  call addfld ('SPWTKE   ', (/ 'lev' /), 'A', 'm/s',      'Standard deviation of updraft velocity')
  call addfld ('SPLCLOUD  ',(/ 'lev' /), 'A', '        ', 'Liquid cloud fraction')

    !call addfld ('SPNDROPMIX','#/kg/s  ',pver,  'A','Droplet number mixing',phys_decomp)
    !call addfld ('SPNDROPSRC','#/kg/s  ',pver,  'A','Droplet number source',phys_decomp)
    !call addfld ('SPNDROPCOL','#/m2    ',1,     'A','Column droplet number',phys_decomp)

    call add_default ('SPKVH     ', 1, ' ')
    call add_default ('SPWTKE    ', 1, ' ')
    call add_default ('SPLCLOUD  ', 1, ' ')
    !call add_default ('SPNDROPSRC', 1, ' ')
    !call add_default ('SPNDROPMIX', 1, ' ')
    !call add_default ('SPNDROPCOL', 1, ' ')

#ifdef MODAL_AERO
! add dropmixnuc tendencies for all modal aerosol species
    do m = 1, ntot_amode
    do lphase = 1, 2
    do lspec = 0, nspec_amode(m)+1   ! loop over number + chem constituents + water
       unit = 'kg/m2/s'
       if (lspec == 0) then   ! number
          unit = '#/m2/s'
          if (lphase == 1) then
             l = numptr_amode(m)
          else
             l = numptrcw_amode(m)
          endif
       else if (lspec <= nspec_amode(m)) then   ! non-water mass
          if (lphase == 1) then
             l = lmassptr_amode(lspec,m)
          else
             l = lmassptrcw_amode(lspec,m)
          endif
       else   ! water mass
!         if (lphase == 1) then
!            l = lwaterptr_amode(m)
!         else
             cycle
!         end if
       end if
       if (lphase == 1) then
          tmpname = cnst_name(l)
       else
          tmpname = cnst_name_cw(l)
       end if

       fieldname = trim(tmpname) // '_mixnuc1sp'
       long_name = trim(tmpname) // ' dropmixnuc mixnuc column tendency in the mmf one '
       call addfld( fieldname,  horiz_only, 'A', unit, long_name)
       call add_default( fieldname, 1, ' ' )

    end do   ! lspec
    end do   ! lphase
    end do   ! m

    do m=1, pcnst
       if(species_class(m).eq.spec_class_gas) then
          fieldname = trim(cnst_name(m)) // '_mixnuc1sp'
          long_name = trim(cnst_name(m)) // ' dropmixnuc mixnuc column tendency in the mmf one '
          call addfld( fieldname,  horiz_only, 'A', unit, long_name)
          call add_default( fieldname, 1, ' ' )
       end if
    end do
#endif

    prec_dp_idx  =  pbuf_get_index('PREC_DP')
    snow_dp_idx  =  pbuf_get_index('SNOW_DP')
    prec_sh_idx  =  pbuf_get_index('PREC_SH')
    snow_sh_idx  =  pbuf_get_index('SNOW_SH')
    prec_sed_idx =  pbuf_get_index('PREC_SED')
    snow_sed_idx =  pbuf_get_index('SNOW_SED')
    prec_pcw_idx =  pbuf_get_index('PREC_PCW')
    snow_pcw_idx =  pbuf_get_index('SNOW_PCW')
    
!-- mdb spcam:  try putting this here
    call physics_tend_alloc(tend_save, pcols)


end subroutine crm_physics_init
!=========================================================================================================

!---------------------------------------------------------------------------------------------------------
   subroutine crm_physics_tend(ztodt, state, tend, ptend, pbuf, dlf, cam_in, cam_out)

!------------------------------------------------------------------------------------------
!  Purpose: to update state from CRM physics. 
! 
! Revision history: 
!
! June, 2009, Minghuai Wang: 
!          These codes are taken out from tphysbc.F90 
!       in the spcam3.5, developed by Marat Khairoutdinov 
!       (mkhairoutdin@ms.cc.sunysb.edu). Here we try to follow the procedure 
!      in 'Interface to Column Physics and Chemistry packages' to implement 
!      the CRM physics.
! July, 13, 2009, Minghuai Wang: 
!      Hydrometer numbers are outputed from SAM when Morrison's microphysics is used, 
!      and will be used in the radiative transfer code to calculate radius. 
! July, 15, 2009, Minghuai Wang: 
!      Get modal aerosol, and use it in the SAM. 
! 
!-------------------------------------------------------------------------------------------
   use ppgrid
   use physics_buffer,  only: physics_buffer_desc, pbuf_old_tim_idx, pbuf_get_index, dyn_time_lvls, pbuf_get_field, pbuf_set_field
   use physics_types,   only: physics_state, physics_tend, physics_ptend, physics_update, physics_ptend_init
   use camsrfexch,      only: cam_in_t, cam_out_t
   use time_manager,    only: is_first_step, get_nstep
   use cam_history,     only: outfld
   use perf_mod
   use crmdims,       only: crm_nx, crm_ny, crm_nz
   use physconst,       only: cpair, latvap, latice, gravit, cappa
   use constituents,    only: pcnst, qmin, cnst_get_ind, cnst_cam_outfld, bpcnst, cnst_name
#ifdef CRM
   use crm_module,    only: crm
   use microphysics,    only: nmicro_fields
#endif
   use physconst,    only: latvap
   use phys_control,  only: phys_getopts
   use check_energy,    only: check_energy_chng

! modal_aero_data only exists if MODAL_AERO
#if (defined  m2005 && defined MODAL_AERO)  
   use crmclouds_camaerosols, only: crmclouds_mixnuc_tend
   use modal_aero_data, only: ntot_amode, ntot_amode
   use ndrop,  only: loadaer
   use microphysics,  only: iqv, iqci, iqr, iqs, iqg, incl, inci, inr, ing, ins   !!!!!!BE CAUIOUS, these indices can only                                                                                   ! defined before call to crm.
#endif
#if (defined MODAL_AERO)  
   use modal_aero_data, only: qqcw_get_field
#endif
#ifdef ECPP
   use module_ecpp_ppdriver2,  only: parampollu_driver2
   use ecppvars,  only: NCLASS_CL, ncls_ecpp_in, NCLASS_PR
   use module_data_ecpp1, only: dtstep_pp_input
#endif
   use phys_grid,       only: get_rlat_all_p, get_rlon_all_p, get_lon_all_p, get_lat_all_p
   use convect_deep,    only: convect_deep_tend_2, deep_scheme_does_scav_trans

   !!!use aerosol_intr,    only: aerosol_wet_intr


   real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)
   type(physics_state), intent(inout) :: state   !  state should be intent(in), but it is changed in this subroutine,
                                                ! and physics_update is also called in this subroutine. It should be changed 
                                                ! in the future. -Minghuai Wang (Minghuai.Wang@pnl.gov). 
   type(physics_tend), intent(inout) :: tend
   type(physics_ptend ), intent(out) :: ptend
   type(physics_buffer_desc),pointer :: pbuf(:)
   real(r8), intent(in) :: dlf(pcols,pver)           ! shallow+deep convective detrainment [kg/kg/s]
   type (cam_in_t), intent(in)      :: cam_in
   type (cam_out_t), intent(inout)  :: cam_out

! convective precipitation variables
   real(r8), pointer :: prec_dp(:)                ! total precipitation from ZM convection [m/s]
   real(r8), pointer :: snow_dp(:)                ! snow from ZM convection                [m/s]
   real(r8), pointer :: prec_sh(:)                ! total precipitation from Hack convection [m/s]
   real(r8), pointer :: snow_sh(:)                ! snow from Hack convection              [m/s]

! stratiform precipitation variables
   real(r8), pointer :: prec_pcw(:)     ! total precip from prognostic cloud scheme   [m/s]
   real(r8), pointer :: snow_pcw(:)     ! snow from prognostic cloud scheme           [m/s]
   real(r8), pointer :: prec_sed(:)     ! total precip from cloud sedimentation       [m/s]
   real(r8), pointer :: snow_sed(:)     ! snow from cloud ice sedimentation           [m/s]


   real(r8), pointer ::  nc_rad(:,:,:,:) ! rad cloud water droplet number [#/kg]
   real(r8), pointer ::  ni_rad(:,:,:,:) ! rad cloud ice crystal number [#/kg]
   real(r8), pointer ::  qs_rad(:,:,:,:) ! rad cloud snow mass [kg/kg]
   real(r8), pointer ::  ns_rad(:,:,:,:) ! rad cloud snow crystal number [#/kg]
   real(r8), pointer ::  cld_rad(:,:,:,:)  ! cloud fraction

   real(r8), pointer ::  t_rad (:,:,:,:) ! rad temperuture
   real(r8), pointer ::  qv_rad(:,:,:,:) ! rad vapor
   real(r8), pointer ::  qc_rad(:,:,:,:) ! rad cloud water
   real(r8), pointer ::  qi_rad(:,:,:,:) ! rad cloud ice
   real(r8), pointer ::  crm_qrad(:,:,:,:)
   real(r8), pointer ::  clubb_buffer  (:,:,:,:,:)
! at layer center

!
!--------------------------- Local variables -----------------------------------------------------------
!
   integer lchnk                              ! chunk identifier
   integer ncol                               ! number of atmospheric columns

   integer nstep     ! time steps

   real(r8) qc_crm (pcols,crm_nx, crm_ny, crm_nz)
   real(r8) qi_crm (pcols,crm_nx, crm_ny, crm_nz)
   real(r8) qpc_crm(pcols,crm_nx, crm_ny, crm_nz)
   real(r8) qpi_crm(pcols,crm_nx, crm_ny, crm_nz)
#ifdef CLUBB_CRM
   real(r8) crm_cld(pcols,crm_nx, crm_ny, crm_nz+1)
#endif
   real(r8) prec_crm(pcols,crm_nx, crm_ny)
   real(r8) mctot(pcols,pver)        ! total cloud mass flux
   real(r8) mcup(pcols,pver)         ! cloud updraft mass flux
   real(r8) mcdn(pcols,pver)         ! cloud downdraft mass flux
   real(r8) mcuup(pcols,pver)        ! unsaturated updraft mass flux
   real(r8) mcudn(pcols,pver)        ! unsaturated downdraft mass flux
   real(r8) spqc(pcols,pver)         ! cloud water
   real(r8) spqi(pcols,pver)         ! cloud ice
   real(r8) spqs(pcols,pver)         ! snow
   real(r8) spqg(pcols,pver)         ! graupel
   real(r8) spqr(pcols,pver)         ! rain
   real(r8) spnc(pcols,pver)         ! cloud water droplet (#/kg)
   real(r8) spni(pcols,pver)         ! cloud ice crystal number (#/kg)
   real(r8) spns(pcols,pver)         ! snow particle number (#/kg)
   real(r8) spng(pcols,pver)         ! graupel particle number (#/kg)
   real(r8) spnr(pcols,pver)         ! rain particle number (#/kg)
   real(r8) wvar_crm (pcols,crm_nx, crm_ny, crm_nz)   ! vertical velocity variance (m/s)

! hm 7/26/11, add new output
   real(r8) aut_crm (pcols,crm_nx, crm_ny, crm_nz)   ! Cloud water autoconversion (1/s)
   real(r8) acc_crm (pcols,crm_nx, crm_ny, crm_nz)   ! Cloud water accretion by rain (1/s)
   real(r8) evpc_crm (pcols,crm_nx, crm_ny, crm_nz)  ! Cloud water evaporation (1/s)
   real(r8) evpr_crm (pcols,crm_nx, crm_ny, crm_nz)  ! Rain evaporation (1/s)
   real(r8) mlt_crm (pcols,crm_nx, crm_ny, crm_nz)   ! Ice, snow, graupel melting (1/s)
   real(r8) sub_crm (pcols,crm_nx, crm_ny, crm_nz)   ! Ice, snow, graupel sublimation (1/s)
   real(r8) dep_crm (pcols,crm_nx, crm_ny, crm_nz)   ! Ice, snow, graupel deposition (1/s)
   real(r8) con_crm (pcols,crm_nx, crm_ny, crm_nz)   ! Cloud water condensation (1/s)
! hm 8/31/11, add new output
   real(r8) aut_crm_a (pcols,pver)   ! Cloud water autoconversion (1/s)
   real(r8) acc_crm_a (pcols,pver)   ! Cloud water accretion by rain (1/s)
   real(r8) evpc_crm_a (pcols,pver)  ! Cloud water evaporation (1/s)
   real(r8) evpr_crm_a (pcols,pver)  ! Rain evaporation (1/s)
   real(r8) mlt_crm_a (pcols,pver)   ! Ice, snow, graupel melting (1/s)
   real(r8) sub_crm_a (pcols,pver)   ! Ice, snow, graupel sublimation (1/s)
   real(r8) dep_crm_a (pcols,pver)   ! Ice, snow, graupel deposition (1/s)
   real(r8) con_crm_a (pcols,pver)   ! Cloud water condensation (1/s)
 
   character(len=16) :: microp_scheme  ! microphysics scheme
   real(r8) flux_qt(pcols,pver)        ! nonprecipitating water flux
   real(r8) flux_u(pcols,pver)        ! x-momentum flux
   real(r8) flux_v(pcols,pver)        ! y-momentum flux
   real(r8) fluxsgs_qt(pcols,pver)     ! sgs nonprecipitating water flux
   real(r8) tkez(pcols,pver)     ! tke profile [kg/m/s2]
   real(r8) tkesgsz(pcols,pver)     ! sgs tke profile  [kg/m/s2]
   real(r8) tkz(pcols, pver)     ! tk profile [m2/s]
   real(r8) flux_qp(pcols,pver)        ! precipitating water flux
   real(r8) precflux(pcols,pver)       ! precipitation flux
   real(r8) qt_ls(pcols,pver)        ! water tendency due to large-scale
   real(r8) qt_trans(pcols,pver)     ! nonprecip water tendency due to transport
   real(r8) qp_trans(pcols,pver)     ! precip water tendency due to transport
   real(r8) qp_fall(pcols,pver)      ! precip water tendency due to fall-out
   real(r8) qp_evp(pcols,pver)       ! precip water tendency due to evap
   real(r8) qp_src(pcols,pver)       ! precip water tendency due to conversion
   real(r8) t_ls(pcols,pver)        ! tendency of crm's liwse due to large-scale
   real(r8) cldtop(pcols,pver)
   real(r8) cwp   (pcols,pver)      ! in-cloud cloud (total) water path (kg/m2)
   real(r8) gicewp(pcols,pver)      ! grid-box cloud ice water path  (g/m2)
   real(r8) gliqwp(pcols,pver)      ! grid-box cloud liquid water path (g/m2)
   real(r8) gwp   (pcols,pver)      ! grid-box cloud (total) water path (kg/m2)
   real(r8) tgicewp(pcols)          ! Vertically integrated ice water path (kg/m2
   real(r8) tgliqwp(pcols)          ! Vertically integrated liquid water path (kg/m2)
   real(r8) cicewp(pcols,pver)             ! in-cloud cloud ice water path (kg/m2)
   real(r8) cliqwp(pcols,pver)             ! in-cloud cloud liquid water path (kg/m2)
   real(r8) tgwp   (pcols)          ! Vertically integrated (total) cloud water path  (kg/m2)
   real(r8) precc(pcols)            ! convective precipitation [m/s]
   real(r8) precl(pcols)            ! large scale precipitation [m/s]
   real(r8) precsc(pcols)           ! convecitve snow   [m/s]
   real(r8) precsl(pcols)           ! convective snow   [m/s]
   real(r8) cltot(pcols)                      ! Diagnostic total cloud cover
   real(r8) cllow(pcols)                       !     "  low cloud cover
   real(r8) clmed(pcols)                      !       "     mid  cloud cover
   real(r8) clhgh(pcols)                      !       "     hgh  cloud cover
   real(r8) :: ftem(pcols,pver)              ! Temporary workspace for outfld variables
   real(r8) ul(pver)
   real(r8) vl(pver)

   real(r8) :: mu_crm(pcols, pver)   
   real(r8) :: md_crm(pcols, pver) 
   real(r8) :: du_crm(pcols, pver)
   real(r8) :: eu_crm(pcols, pver)
   real(r8) :: ed_crm(pcols, pver)
   real(r8) :: jt_crm(pcols)
   real(r8) :: mx_crm(pcols)
   real(r8) :: ideep_crm(pcols)

! physics buffer fields to compute tendencies for stratiform package
   integer itim 
   real(r8), pointer, dimension(:,:) :: cld        ! cloud fraction

#if (defined m2005 && defined MODAL_AERO)
   real(r8) na(pcols)              ! aerosol number concentration [/m3]
   real(r8) va(pcols)              ! aerosol voume concentration [m3/m3]
   real(r8) hy(pcols)              ! aerosol bulk hygroscopicity
   real(r8) naermod(pver, ntot_amode)     ! Aerosol number concentration [/m3]
   real(r8) vaerosol(pver, ntot_amode)    ! aerosol volume concentration [m3/m3]
   real(r8) hygro(pver, ntot_amode)       ! hygroscopicity of aerosol mode 
   !real(r8) cs(pcols, pver)                     ! air density  [kg/m3]
   integer   phase                  ! phase to determine whether it is interstitial, cloud-borne, or the sum. 
#endif

!-- mark branson fix for cs compiler error using 1-moment microphysics
   real(r8) cs(pcols, pver)                     ! air density  [kg/m3]

#ifdef ECPP
   real(r8) :: qicecen(pcols,pver,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)   ! cloud ice (kg/kg)
   real(r8) :: qlsink_afcen(pcols,pver,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! cloud water loss rate from precipitation calculated
                                                                                      ! cloud water before precipitatinog (/s)
   real(r8) :: qlsink_bfcen(pcols,pver,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! cloud water loss rate from precipitation calculated
                                                                                      ! cloud water before precipitatinog (/s)
   real(r8) :: qlsink_avgcen(pcols,pver,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! cloud water loss rate from precipitation calculated
                                                                                      ! from praincen and qlcoudcen averaged over
                                                                                      ! ntavg1_ss time step (/s)
   real(r8) :: praincen(pcols,pver,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! cloud water loss rate from precipitation (kg/kg/s)
   real(r8) :: wupthresh_bnd(pcols, pverp) 
   real(r8) :: wdownthresh_bnd(pcols, pverp)
#endif

! CRM column radiation stuff:
   real(r8) prectend(pcols) ! tendency in precipitating water and ice
   real(r8) precstend(pcols) ! tendency in precipitating ice
   real(r8) wtricesink(pcols) ! sink of water vapor + cloud water + cloud ice
   real(r8) icesink(pcols) ! sink of
   real(r8) tau00  ! surface stress
   real(r8) wnd  ! surface wnd
   real(r8) bflx   ! surface buoyancy flux (Km/s)
   real(r8) taux_crm(pcols)  ! zonal CRM surface stress perturbation
   real(r8) tauy_crm(pcols)  ! merid CRM surface stress perturbation
   real(r8) z0m(pcols)  ! surface momentum roughness length
   real(r8), pointer, dimension(:,:) :: qrs, qrl        ! rad heating rates
   real(r8), pointer, dimension(:,:) :: tempPtr
   real(r8), pointer, dimension(:,:,:,:)   :: crm_u
   real(r8), pointer, dimension(:,:,:,:)   :: crm_v
   real(r8), pointer, dimension(:,:,:,:)   :: crm_w
   real(r8), pointer, dimension(:,:,:,:)   :: crm_t
   real(r8), pointer, dimension(:,:,:,:)   :: crm_qt
   real(r8), pointer, dimension(:,:,:,:)   :: crm_qp
   real(r8), pointer, dimension(:,:,:,:)   :: crm_qn
   real(r8), pointer, dimension(:,:,:,:)   :: crm_nc
   real(r8), pointer, dimension(:,:,:,:)   :: crm_qr
   real(r8), pointer, dimension(:,:,:,:)   :: crm_nr
   real(r8), pointer, dimension(:,:,:,:)   :: crm_qi
   real(r8), pointer, dimension(:,:,:,:)   :: crm_ni
   real(r8), pointer, dimension(:,:,:,:)   :: crm_qs
   real(r8), pointer, dimension(:,:,:,:)   :: crm_ns
   real(r8), pointer, dimension(:,:,:,:)   :: crm_qg
   real(r8), pointer, dimension(:,:,:,:)   :: crm_ng
   real(r8), pointer, dimension(:,:,:,:)   :: crm_qc

#ifdef CRM
   real(r8), dimension(pcols,crm_nx,crm_ny,crm_nz,nmicro_fields+1) :: crm_micro
#endif

   integer                         :: pblh_idx
   real(r8), pointer, dimension(:) :: pblh
   real(r8), pointer, dimension(:,:) :: qqcw

#ifdef ECPP
   real(r8),allocatable :: acen(:,:,:,:,:)   ! cloud fraction for each sub-sub class for full time period
   real(r8),allocatable :: acen_tf(:,:,:,:,:) ! cloud fraction for end-portion of time period
   real(r8),allocatable :: rhcen(:,:,:,:,:)  ! relative humidity (0-1)
   real(r8),allocatable :: qcloudcen(:,:,:,:,:)  ! cloud water (kg/kg)
   real(r8),allocatable :: qlsinkcen(:,:,:,:,:)  ! cloud water loss rate from precipitation (/s??)
   real(r8),allocatable :: precrcen(:,:,:,:,:)   ! liquid (rain) precipitation rate (kg/m2/s)
   real(r8),allocatable :: precsolidcen(:,:,:,:,:)   ! solid (rain) precipitation rate (kg/m2/s)
   real(r8),allocatable :: wwqui_cen(:,:)                                ! vertical velocity variance in quiescent class (m2/s2)
   real(r8),allocatable :: wwqui_cloudy_cen(:,:)                         ! vertical velocity variance in quiescent, and cloudy class (m2/s2)
! at layer boundary
   real(r8),allocatable :: abnd(:,:,:,:,:)       ! cloud fraction for each sub-sub class for full time period
   real(r8),allocatable :: abnd_tf(:,:,:,:,:)    ! cloud fraction for end-portion of time period
   real(r8),allocatable :: massflxbnd(:,:,:,:,:) ! sub-class vertical mass flux (kg/m2/s) at layer bottom boundary.
   real(r8),allocatable :: wwqui_bnd(:,:)        ! vertical velocity variance in quiescent class (m2/s2)
   real(r8),allocatable :: wwqui_cloudy_bnd(:,:) ! vertical velocity variance in quiescent, and cloudy class (m2/s2)
#endif

   real(r8)             :: qaer(pcols, pver, pcnst)
#ifdef MODAL_AERO
   real(r8)             :: qaerwat_rad(pcols, crm_nx, crm_ny, crm_nz, ntot_amode)  ! aerosol water
   real(r8)             :: dgnumwet_rad(pcols, crm_nx, crm_ny, crm_nz, ntot_amode) ! wet mode dimaeter
#endif

! Surface fluxes +++mhwang
   real(r8) ::  fluxu0           ! surface momenment fluxes
   real(r8) ::  fluxv0           ! surface momenment fluxes
   real(r8) ::  fluxt0           ! surface sensible heat fluxes
   real(r8) ::  fluxq0           ! surface latent heat fluxes
   real(r8) ::  dtstep_pp        ! time step for the ECPP (seconds)
   integer  ::  necpp            ! the number of GCM time step in which ECPP is called once.


   real(r8) radflux(pcols)   ! radiative fluxes from radiation calculation (qrs + qrl)

!+++mhwangtest
   real(r8) qtot(pcols, 3)   ! total water
   real(r8) qt_hydro(pcols, 2)  ! total hydrometer 
   real(r8) qt_cloud(pcols, 3)  ! total cloud water
   real(r8) qtv(pcols, 3)     ! total water vapor
   real(r8) ql_cloud(pcols, 2) ! liquid cloud water
   real(r8) qli_hydro(pcols, 2)  ! column-integraetd rain + snow + graupel 
   real(r8) qi_hydro(pcols, 2)  ! column-integrated snow water + graupel water
   real(r8) sfactor 
!---mhwangtest

   real(r8) zero(pcols)      ! zero
   real(r8) timing_factor(pcols) ! factor for crm cpu-usage: 1 means no subcycling

   real(r8) qtotcrm(pcols, 20)   ! the toal water calculated in crm.F90

   integer nlon(pcols)
   integer nlat(pcols)

   integer ii, jj, mm
   integer iii,lll
   integer ixcldliq, ixcldice, ixnumliq, ixnumice
   integer i, k, m
   integer ifld
   logical :: use_ECPP, use_SPCAM
   character(len=16) :: SPCAM_microp_scheme

!-- MDB 8/2013
   real(r8) tvwle(pcols,pver),buoy(pcols,pver),buoysd(pcols,pver)
   real(r8) msef(pcols,pver),qvw(pcols,pver)

   logical :: ls, lu, lv, lq(pcnst), fromcrm

   zero = 0.0_r8
!========================================================
!========================================================
!  CRM (Superparameterization).
! Author: Marat Khairoutdinov (mkhairoutdin@ms.cc.sunysb.edu)
!========================================================

   call t_startf ('crm')

!-- mdb spcam: ACME version of physics_ptend_init has additional arguments
   !!!call physics_ptend_init(ptend)                ! Initialize output physics_ptend object
   lu = .true. 
   lv = .true.
   ls = .true.
   lq(:) = .true.
   fromcrm = .true.
   call physics_ptend_init(ptend,     state%psetcols, 'crm', lu=lu, lv=lv, ls=ls, lq=lq, fromcrm=fromcrm)  ! Initialize output physics_ptend object
   fromcrm = .false.

   call phys_getopts(use_ECPP_out            = use_ECPP)
   call phys_getopts(use_SPCAM_out           = use_SPCAM)
   call phys_getopts(SPCAM_microp_scheme_out = SPCAM_microp_scheme)

   nstep = get_nstep()

   lchnk = state%lchnk
   ncol  = state%ncol

   if (SPCAM_microp_scheme .eq. 'm2005') then
     call pbuf_get_field(pbuf, crm_nc_rad_idx, nc_rad, start=(/1,1,1,1/), kount=(/pcols,crm_nx, crm_ny, crm_nz/))
     call pbuf_get_field(pbuf, crm_ni_rad_idx, ni_rad, start=(/1,1,1,1/), kount=(/pcols,crm_nx, crm_ny, crm_nz/))
     call pbuf_get_field(pbuf, crm_qs_rad_idx, qs_rad, start=(/1,1,1,1/), kount=(/pcols,crm_nx, crm_ny, crm_nz/))
     call pbuf_get_field(pbuf, crm_ns_rad_idx, ns_rad, start=(/1,1,1,1/), kount=(/pcols,crm_nx, crm_ny, crm_nz/))
     call pbuf_get_field(pbuf, crm_cld_rad_idx, cld_rad, start=(/1,1,1,1/), kount=(/pcols,crm_nx, crm_ny, crm_nz/))
   endif

#ifdef ECPP
    if (use_ECPP) then
       allocate( acen(pcols,pver,NCLASS_CL,ncls_ecpp_in,NCLASS_PR))
       allocate( acen_tf(pcols,pver,NCLASS_CL,ncls_ecpp_in,NCLASS_PR))
       allocate( rhcen(pcols,pver,NCLASS_CL,ncls_ecpp_in,NCLASS_PR))
       allocate( qcloudcen(pcols,pver,NCLASS_CL,ncls_ecpp_in,NCLASS_PR))
       allocate( qlsinkcen(pcols,pver,NCLASS_CL,ncls_ecpp_in,NCLASS_PR))
       allocate( precrcen(pcols,pver,NCLASS_CL,ncls_ecpp_in,NCLASS_PR))
       allocate( precsolidcen(pcols,pver,NCLASS_CL,ncls_ecpp_in,NCLASS_PR))
       allocate( wwqui_cen(pcols, pver))
       allocate( wwqui_cloudy_cen(pcols, pver))
! at layer boundary
       allocate( abnd(pcols,pver+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR))
       allocate( abnd_tf(pcols,pver+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR))
       allocate( massflxbnd(pcols,pver+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR))
       allocate( wwqui_bnd(pcols, pver+1))
       allocate( wwqui_cloudy_bnd(pcols, pver+1))
    end if
#endif

   if(is_first_step())then
      do m=1,dyn_time_lvls
         call pbuf_get_field(pbuf, cldo_idx, cldo, start=(/1,1,m/), kount=(/pcols,pver,1/) )
         cldo(:ncol,:) = 0
      enddo
   endif


! Forget all changes to the state due to conventional physics above:

! gas and aerosol species are updated in by shallow convective transport in CAM.
! Aerosol changes from dropmixnuc in cldwat2m.F90 is discarded. (When dropmixnuc
! is called in cldwat2m.F90, the tenendcy is set to zero. See dropmixnuc for details).
! Minghuai Wang, 2010-01 (Minghuai.Wang@pnl.gov)
!
   qaer = state%q

   state = state_save
   tend  = tend_save
   cldo(:ncol, :) = cldo_save(:ncol, :)

#if ( defined MODAL_AERO )
   do i=1,pcnst
      qqcw   =>  qqcw_get_field(pbuf, i,lchnk, .true.)
      if (associated(qqcw)) qqcw(:, :) = qqcw_save(:,:,i)
   end do
   dgnumwet = dgnumwet_save
#endif

!
! tracer species other than water vapor and cloud water are updated in convetional CAM.
! When ECPP is used, dropmixnuc and all transport(deep and shallow) are done in ECPP.
! So all change in aerosol and gas species in the conventinal CAM are discarded.
! Minghuai Wang, 2010-01 (Minghuai.Wang@pnl.gov)
!
   if (.not. use_ECPP) then
       call phys_getopts(microp_scheme_out=microp_scheme)
       if ( microp_scheme .eq. 'MG' ) then
         state%q(:ncol, :pver, 6:pcnst) = qaer(:ncol, :pver, 6:pcnst)
       else if ( microp_scheme .eq. 'RK' ) then
         state%q(:ncol, :pver, 4:pcnst) = qaer(:ncol, :pver, 4:pcnst)
       end if
   endif

! Associate pointers with physics buffer fields
   itim = pbuf_old_tim_idx()
   ifld = pbuf_get_index('CLD')
   call pbuf_get_field(pbuf, ifld, cld, start=(/1,1,itim/), kount=(/pcols,pver,1/) )

   call pbuf_get_field (pbuf, crm_u_idx,         crm_u)
   call pbuf_get_field (pbuf, crm_v_idx,         crm_v)
   call pbuf_get_field (pbuf, crm_w_idx,         crm_w)
   call pbuf_get_field (pbuf, crm_t_idx,         crm_t)

   if (SPCAM_microp_scheme .eq. 'sam1mom') then
      call pbuf_get_field (pbuf, crm_qt_idx,        crm_qt)
      call pbuf_get_field (pbuf, crm_qp_idx,        crm_qp)
      call pbuf_get_field (pbuf, crm_qn_idx,        crm_qn)
   else if (SPCAM_microp_scheme .eq. 'm2005') then
      call pbuf_get_field (pbuf, crm_qt_idx,        crm_qt)
      call pbuf_get_field (pbuf, crm_nc_idx,        crm_nc)
      call pbuf_get_field (pbuf, crm_qr_idx,        crm_qr)
      call pbuf_get_field (pbuf, crm_nr_idx,        crm_nr)
      call pbuf_get_field (pbuf, crm_qi_idx,        crm_qi)
      call pbuf_get_field (pbuf, crm_ni_idx,        crm_ni)
      call pbuf_get_field (pbuf, crm_qs_idx,        crm_qs)
      call pbuf_get_field (pbuf, crm_ns_idx,        crm_ns)
      call pbuf_get_field (pbuf, crm_qg_idx,        crm_qg)
      call pbuf_get_field (pbuf, crm_ng_idx,        crm_ng)
      call pbuf_get_field (pbuf, crm_qc_idx,        crm_qc)
   end if
   call pbuf_get_field (pbuf, crm_qrad_idx,      crm_qrad)
#ifdef CLUBB_CRM
   call pbuf_get_field (pbuf, clubb_buffer_idx,  clubb_buffer)
#endif

   call pbuf_get_field (pbuf, crm_t_rad_idx,  t_rad)
   call pbuf_get_field (pbuf, crm_qv_rad_idx, qv_rad)
   call pbuf_get_field (pbuf, crm_qc_rad_idx, qc_rad)
   call pbuf_get_field (pbuf, crm_qi_rad_idx, qi_rad)


   
! Initialize stuff:
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)

   if(is_first_step()) then
!      call check_energy_timestep_init(state, tend, pbuf)
       do i=1,ncol
         do k=1,crm_nz
               m = pver-k+1

#ifdef SP_DIR_NS
               if(crm_ny.eq.1) then ! change domain orientation only for 2D CRM
                 crm_u(i,:,:,k) = state%v(i,m)
                 crm_v(i,:,:,k) = state%u(i,m)
               else
                 crm_u(i,:,:,k) = state%u(i,m)
                 crm_v(i,:,:,k) = state%v(i,m)
               end if
#else
               crm_u(i,:,:,k) = state%u(i,m)
               crm_v(i,:,:,k) = state%v(i,m)
#endif
               crm_w(i,:,:,k) = 0.
               crm_t(i,:,:,k) = state%t(i,m)

! Load the crm_micro array for use in the call to crm

#ifdef CRM
               if (SPCAM_microp_scheme .eq. 'sam1mom') then

                  crm_qt(i,:,:,k) = state%q(i,m,1)+state%q(i,m,ixcldliq)+state%q(i,m,ixcldice)
                  crm_qp(i,:,:,k) = 0.0_r8
                  crm_qn(i,:,:,k) = state%q(i,m,ixcldliq)+state%q(i,m,ixcldice)
                  crm_micro(i,:,:,k,1) = crm_qt(i,:,:,k)
                  crm_micro(i,:,:,k,2) = crm_qp(i,:,:,k)
                  crm_micro(i,:,:,k,3) = crm_qn(i,:,:,k)

               else if (SPCAM_microp_scheme .eq. 'm2005') then
                  crm_qt(i,:,:,k) = state%q(i,m,1)+state%q(i,m,ixcldliq)
                  crm_nc(i,:,:,k) = 0.0_r8
                  crm_qr(i,:,:,k) = 0.0_r8
                  crm_nr(i,:,:,k) = 0.0_r8
                  crm_qi(i,:,:,k) = state%q(i,m,ixcldice) 
                  crm_ni(i,:,:,k) = 0.0_r8
                  crm_qs(i,:,:,k) = 0.0_r8
                  crm_ns(i,:,:,k) = 0.0_r8
                  crm_qg(i,:,:,k) = 0.0_r8
                  crm_ng(i,:,:,k) = 0.0_r8
                  crm_qc(i,:,:,k) = state%q(i,m,ixcldliq)

                  crm_micro(i,:,:,k,1)  = crm_qt(i,:,:,k)
                  crm_micro(i,:,:,k,2)  = crm_nc(i,:,:,k)
                  crm_micro(i,:,:,k,3)  = crm_qr(i,:,:,k)
                  crm_micro(i,:,:,k,4)  = crm_nr(i,:,:,k)
                  crm_micro(i,:,:,k,5)  = crm_qi(i,:,:,k)
                  crm_micro(i,:,:,k,6)  = crm_ni(i,:,:,k)
                  crm_micro(i,:,:,k,7)  = crm_qs(i,:,:,k)
                  crm_micro(i,:,:,k,8)  = crm_ns(i,:,:,k)
                  crm_micro(i,:,:,k,9)  = crm_qg(i,:,:,k)
                  crm_micro(i,:,:,k,10) = crm_ng(i,:,:,k)
                  crm_micro(i,:,:,k,11) = crm_qc(i,:,:,k)
               endif
#endif

#ifdef CLUBB_CRM
              clubb_buffer(i,:,:,k,:) = 0.0  ! In the inital run, variables are set in clubb_sgs_setup at the first time step. 
#endif
         end do
       end do

       do k=1,crm_nz
         m = pver-k+1
         do i=1,ncol
             crm_qrad (i,:,:,k)    = 0.
             qc_crm (i,:,:,k)      = 0.
             qi_crm (i,:,:,k)      = 0.
             qpc_crm(i,:,:,k)      = 0.
             qpi_crm(i,:,:,k)      = 0.
             t_rad  (i,:,:,k)      = state%t(i,m)
             qv_rad (i,:,:,k)      = state%q(i,m,1)
             qc_rad (i,:,:,k)      = 0.
             qi_rad (i,:,:,k)      = 0.
#ifdef m2005
             if (SPCAM_microp_scheme .eq. 'm2005') then
                nc_rad(i,:,:,k) = 0.
                ni_rad(i,:,:,k) = 0.       
                qs_rad(i,:,:,k) = 0.0
                ns_rad(i,:,:,k) = 0.0
                wvar_crm(i,:,:,k) = 0.0
! hm 7/26/11, add new output
                aut_crm(i,:,:,k) = 0.0
                acc_crm(i,:,:,k) = 0.0
                evpc_crm(i,:,:,k) = 0.0
                evpr_crm(i,:,:,k) = 0.0
                mlt_crm(i,:,:,k) = 0.0
                sub_crm(i,:,:,k) = 0.0
                dep_crm(i,:,:,k) = 0.0
                con_crm(i,:,:,k) = 0.0
             endif
#endif

         end do
       end do

! use radiation from grid-cell mean radctl on first time step
       prec_crm (:,:,:) = 0.
       ptend%q(:,:,1) = 0.
       ptend%q(:,:,ixcldliq) = 0.
       ptend%q(:,:,ixcldice) = 0.
       ptend%s(:,:) = 0.
       precc(:) = 0.
       precl(:) = 0.
       precsc(:) = 0.
       precsl(:) = 0.
       cltot(:) = 0.
       clhgh(:) = 0.
       clmed(:) = 0.
       cllow(:) = 0.
       cld(:,:) = 0.
       cldtop(:,:) = 0.
       gicewp(:,:)=0
       gliqwp(:,:)=0
       mctot(:,:) = 0.
       mcup(:,:) = 0.
       mcdn(:,:) = 0.
       mcuup(:,:) = 0.
       mcudn(:,:) = 0.
       spqc(:,:) = 0.
       spqi(:,:) = 0.
       spqs(:,:) = 0.
       spqg(:,:) = 0.
       spqr(:,:) = 0.
       if (SPCAM_microp_scheme .eq. 'm2005') then
          spnc(:,:) = 0.
          spni(:,:) = 0.
          spns(:,:) = 0.
          spng(:,:) = 0.
          spnr(:,:) = 0.
! hm 8/31/11, add new output
          aut_crm_a(:,:) = 0.
          acc_crm_a(:,:) = 0.
          evpc_crm_a(:,:) = 0.
          evpr_crm_a(:,:) = 0.
          mlt_crm_a(:,:) = 0.
          sub_crm_a(:,:) = 0.
          dep_crm_a(:,:) = 0.
          con_crm_a(:,:) = 0.
       endif 
       flux_qt(:,:) = 0.
       flux_u(:,:) = 0.
       flux_v(:,:) = 0.
       fluxsgs_qt(:,:) = 0.
       tkez(:,:) = 0.
       tkesgsz(:,:) = 0.
       tkz(:,:) = 0.
       flux_qp(:,:) = 0.
       precflux(:,:) = 0.
       qt_ls(:,:) = 0.
       qt_trans(:,:) = 0.
       qp_trans(:,:) = 0.
       qp_fall(:,:) = 0.
       qp_evp(:,:) = 0.
       qp_src(:,:) = 0.
       z0m(:) = 0.
       taux_crm(:) = 0.
       tauy_crm(:) = 0.
       t_ls(:,:) = 0.
       tvwle(:,:) = 0.    ! MDB 8/2013
       buoy(:,:) = 0.     ! MDB 8/2013
       buoysd(:,:) = 0.   ! MDB 8/2013
       msef(:,:) = 0.     ! MDB 8/2013
       qvw(:,:) = 0.      ! MDB 8/2013


#ifdef ECPP
       if (use_ECPP) then
         abnd=0.0
         abnd_tf=0.0
         massflxbnd=0.0
         acen=0.0
         acen_tf=0.0
         rhcen=0.0
         qcloudcen=0.0
         qicecen=0.0
         qlsinkcen=0.0
         precrcen=0.0
         precsolidcen=0.0
         wupthresh_bnd = 0.0
         wdownthresh_bnd = 0.0
         qlsink_afcen = 0.0
         qlsink_bfcen = 0.0
         qlsink_avgcen = 0.0
         praincen = 0.0
! default is clear, non-precipitating, and quiescent class
         abnd(:,:,1,1,1)=1.0 
         abnd_tf(:,:,1,1,1)=1.0 
         acen(:,:,1,1,1)=1.0 
         acen_tf(:,:,1,1,1)=1.0 
         wwqui_cen = 0.0
         wwqui_bnd = 0.0
         wwqui_cloudy_cen = 0.0
         wwqui_cloudy_bnd = 0.0
! turbulence
         ifld = pbuf_get_index('TKE_CRM')
         cs(:ncol, 1:pver) = state%pmid(:ncol, 1:pver)/(287.15*state%t(:ncol, 1:pver))
         call pbuf_set_field(pbuf, ifld, 0.0_r8, start=(/1,1/), kount=(/pcols, pver/) )

         ifld = pbuf_get_index('TK_CRM')
         call pbuf_set_field(pbuf, ifld, 0.0_r8, start=(/1,1/), kount=(/pcols, pver/) )
       endif 
#endif

   else  ! not is_first_step

#ifdef CRM
       if (SPCAM_microp_scheme .eq. 'sam1mom') then
          crm_micro(:,:,:,:,1) = crm_qt(:,:,:,:)
          crm_micro(:,:,:,:,2) = crm_qp(:,:,:,:)
          crm_micro(:,:,:,:,3) = crm_qn(:,:,:,:)
       else if (SPCAM_microp_scheme .eq. 'm2005') then
          crm_micro(:,:,:,:,1)  = crm_qt(:,:,:,:)
          crm_micro(:,:,:,:,2)  = crm_nc(:,:,:,:)
          crm_micro(:,:,:,:,3)  = crm_qr(:,:,:,:)
          crm_micro(:,:,:,:,4)  = crm_nr(:,:,:,:)
          crm_micro(:,:,:,:,5)  = crm_qi(:,:,:,:)
          crm_micro(:,:,:,:,6)  = crm_ni(:,:,:,:)
          crm_micro(:,:,:,:,7)  = crm_qs(:,:,:,:)
          crm_micro(:,:,:,:,8)  = crm_ns(:,:,:,:)
          crm_micro(:,:,:,:,9)  = crm_qg(:,:,:,:)
          crm_micro(:,:,:,:,10) = crm_ng(:,:,:,:)
          crm_micro(:,:,:,:,11) = crm_qc(:,:,:,:)
       endif
#endif

       ptend%q(:,:,1) = 0.  ! necessary?
       ptend%q(:,:,ixcldliq) = 0.
       ptend%q(:,:,ixcldice) = 0.
       ptend%s(:,:) = 0. ! necessary?
       cwp   = 0.
       gicewp = 0.
       gliqwp = 0.
       cltot = 0.
       clhgh = 0.
       clmed = 0.
       cllow = 0.

       qc_crm  = 0.
       qi_crm  = 0.
       qpc_crm = 0.
       qpi_crm = 0
       prec_crm = 0.

       if (SPCAM_microp_scheme .eq. 'm2005') then
! hm 8/31/11, initialize gcm-time-step-avg output at start of each time step
! IS THIS NECESSARY??? 
         aut_crm_a = 0.0
         acc_crm_a = 0.0
         evpc_crm_a = 0.0
         evpr_crm_a = 0.0
         mlt_crm_a = 0.0
         sub_crm_a = 0.0
         dep_crm_a = 0.0
         con_crm_a = 0.0
       end if

!===================================================================================
!!!!!! should other varialbes also be set to be zero (suc as precc)? !!!!!!!!!!
!  -Minghuai Wang
!===================================================================================

       call t_startf ('crm_call')

       do m=1,crm_nz
         k = pver-m+1
         do i = 1,ncol
             crm_qrad(i,:,:,m) = crm_qrad(i,:,:,m) / state%pdel(i,k) ! for energy conservation
         end do
       end do

#if (defined m2005 && defined MODAL_AERO)
       cs(1:ncol, 1:pver) = state%pmid(1:ncol, 1:pver)/(287.15*state%t(1:ncol, 1:pver))
#endif

       call get_lat_all_p(lchnk, ncol, nlat)
       call get_lon_all_p(lchnk, ncol, nlon)


       do i = 1,ncol

          tau00 = sqrt(cam_in%wsx(i)**2 + cam_in%wsy(i)**2)
          wnd = sqrt(state%u(i,pver)**2 + state%v(i,pver)**2)
          bflx = cam_in%shf(i)/cpair + 0.61*state%t(i,pver)*cam_in%lhf(i)/latvap

!+++mhwang
         fluxu0 = cam_in%wsx(i)     !N/m2
         fluxv0 = cam_in%wsy(i)     !N/m2
         fluxt0 = cam_in%shf(i)/cpair  ! K Kg/ (m2 s)
         fluxq0 = cam_in%lhf(i)/latvap ! Kg/(m2 s)
!---mwwang

!+++mhwangtest
!
! calculate total water before calling crm
! total hydrometer water (rain, snow, and graupel)
!-- mark branson: insert ifdef m2005 around this test so 1-moment micro compile will work
         if (SPCAM_microp_scheme .eq. 'm2005') then
           qt_hydro(i, 1) = 0.0_r8
           qli_hydro(i, 1) = 0.0_r8
           qi_hydro(i, 1) = 0.0_r8
           do m=1, crm_nz
             k=pver-m+1
             do ii=1, crm_nx
             do jj=1, crm_ny
               qt_hydro(i,1) = qt_hydro(i,1)+(crm_qr(i,ii,jj,m)+crm_qs(i,ii,jj,m)+crm_qg(i,ii,jj,m)) * state%pdel(i,k)/gravit 
               qli_hydro(i,1) = qli_hydro(i,1)+(crm_qr(i,ii,jj,m)+crm_qs(i,ii,jj,m)+crm_qg(i,ii,jj,m)) * state%pdel(i,k)/gravit
               qi_hydro(i,1) = qi_hydro(i,1)+(crm_qs(i,ii,jj,m)+crm_qg(i,ii,jj,m)) * state%pdel(i,k)/gravit
             end do
             end do 
           end do
           qt_hydro(i,1) = qt_hydro(i,1)/(crm_nx*crm_ny)
           qli_hydro(i,1) = qli_hydro(i,1)/(crm_nx*crm_ny)
           qi_hydro(i,1) = qi_hydro(i,1)/(crm_nx*crm_ny)

! total cloud water and total water vapor
            qt_cloud(i,1) = 0._r8
            qtv(i,1) = 0._r8
            do k=1, pver
              qt_cloud(i,1) = qt_cloud(i,1) + (state%q(i,k,ixcldliq)+state%q(i,k,ixcldice)) * state%pdel(i,k)/gravit
              qtv(i,1) = qtv(i,1) + state%q(i,k,1) * state%pdel(i,k)/gravit
            end do 
! total water 
            qtot(i,1) = qt_hydro(i,1) + qt_cloud(i,1) + qtv(i,1)
!---mhwangtest
         else if (SPCAM_microp_scheme .eq. 'sam1mom') then
           qli_hydro(i, 1) = 0.0_r8
           qi_hydro(i, 1) = 0.0_r8
           do m=1, crm_nz
             k=pver-m+1
             do ii=1, crm_nx
             do jj=1, crm_ny
               sfactor = max(0._r8,min(1._r8,(crm_t(i,ii,jj,m)-268.16)*1./(283.16-268.16)))
               qli_hydro(i,1) = qli_hydro(i,1)+crm_qp(i,ii,jj,m) * state%pdel(i,k)/gravit
               qi_hydro(i,1) = qi_hydro(i,1)+crm_qp(i,ii,jj,m) * (1-sfactor) * state%pdel(i,k)/gravit
             end do
             end do
           end do
           qli_hydro(i,1) = qli_hydro(i,1)/(crm_nx*crm_ny)
           qi_hydro(i,1) = qi_hydro(i,1)/(crm_nx*crm_ny)    

! total cloud water and total water vapor, and energy
            qt_cloud(i,1) = 0._r8
            ql_cloud(i,1) = 0._r8
            qtv(i,1) = 0._r8
            do k=1, pver
              qt_cloud(i,1) = qt_cloud(i,1) + (state%q(i,k,ixcldliq)+state%q(i,k,ixcldice)) * state%pdel(i,k)/gravit
              ql_cloud(i,1) = ql_cloud(i,1) + state%q(i,k,ixcldliq) * state%pdel(i,k)/gravit
              qtv(i,1) = qtv(i,1) + state%q(i,k,1) * state%pdel(i,k)/gravit
            end do
         endif

#if (defined m2005 && defined MODAL_AERO)
         do k=1, pver
           phase = 1  ! interstital aerosols only
           do m=1, ntot_amode
            call loadaer( &
                  state, pbuf, i, i, k, &
                  m, cs, phase, na, va, &
                  hy)
             naermod(k, m) = na(i)
             vaerosol(k, m) = va(i)
             hygro(k, m) = hy(i)
           end do    
         end do
#endif
#ifdef SP_DIR_NS
         if(crm_ny.eq.1) then
           ul(:) = state%v(i,:)  ! change orientation only if 2D CRM
           vl(:) = state%u(i,:)
         else
           ul(:) = state%u(i,:)
           vl(:) = state%v(i,:)
         end if
#else
           ul(:) = state%u(i,:)
           vl(:) = state%v(i,:)
#endif
#ifdef CRM
           if (.not.allocated(ptend%q)) write(*,*) '=== ptend%q not allocated ===',i
           if (.not.allocated(ptend%s)) write(*,*) '=== ptend%s not allocated ===',i
           call crm (lchnk,           i,                                                                                                   &
             state%t(i,:),            state%q(i,:,1),           state%q(i,:,ixcldliq), state%q(i,:,ixcldice),                              &
             ul(:),                   vl(:),                                                                                               &
             state%ps(i),             state%pmid(i,:),          state%pdel(i,:),       state%phis(i),                                      &
             state%zm(i,:),           state%zi(i,:),            ztodt,                 pver,                                               &
             !ptend%u(i,:),            ptend%v(i,:),                                                                                        &
             ptend%q(i,:,1),          ptend%q(i,:,ixcldliq),    ptend%q(i,:,ixcldice), ptend%s(i,:),                                       &
             crm_u(i,:,:,:),          crm_v(i,:,:,:),           crm_w(i,:,:,:),        crm_t(i,:,:,:),          crm_micro(i,:,:,:,:),      &
             crm_qrad(i,:,:,:),                                                                                                            &
             qc_crm(i,:,:,:),         qi_crm(i,:,:,:),          qpc_crm(i,:,:,:),      qpi_crm(i,:,:,:),                                   &
             prec_crm(i,:,:),         t_rad(i,:,:,:),           qv_rad(i,:,:,:),                                                           &
             qc_rad(i,:,:,:),         qi_rad(i,:,:,:),                                                                                     &
#ifdef m2005
             nc_rad(i,:,:,:),         ni_rad(i,:,:,:),          qs_rad(i,:,:,:),       ns_rad(i,:,:,:),         wvar_crm(i,:,:,:),         &
             aut_crm(i,:,:,:),        acc_crm(i,:,:,:),         evpc_crm(i,:,:,:),     evpr_crm(i,:,:,:),       mlt_crm(i,:,:,:),          &
             sub_crm(i,:,:,:),        dep_crm(i,:,:,:),         con_crm(i,:,:,:),                                                          &
             aut_crm_a(i,:),          acc_crm_a(i,:),           evpc_crm_a(i,:),       evpr_crm_a(i,:),         mlt_crm_a(i,:),            &
             sub_crm_a(i,:),          dep_crm_a(i,:),           con_crm_a(i,:),                                                            &
#endif
             precc(i),                precl(i),                 precsc(i),             precsl(i),                                          &
             cltot(i),                clhgh(i),                 clmed(i),              cllow(i),                cld(i,:),    cldtop(i,:) , &
             gicewp(i,:),             gliqwp(i,:),     &
             mctot(i,:),              mcup(i,:),                mcdn(i,:),             mcuup(i,:),              mcudn(i,:),                &
             spqc(i,:),               spqi(i,:),                spqs(i,:),             spqg(i,:),               spqr(i,:),               &
#ifdef m2005
             spnc(i,:),               spni(i,:),                spns(i,:),             spng(i,:),               spnr(i,:),               &
#ifdef MODAL_AERO
             naermod,                 vaerosol,                 hygro,                                                                     &
#endif 
#endif
#ifdef CLUBB_CRM
             clubb_buffer(i,:,:,:,:),                                                                                                      &
             crm_cld(i,:, :, :),                                                                                                           &
#endif
             mu_crm(i,:),             md_crm(i,:),              du_crm(i,:),           eu_crm(i,:),                                        & 
             ed_crm(i,:),             jt_crm(i),                mx_crm(i),                                                                 &
#ifdef ECPP
             abnd(i,:,:,:,:),         abnd_tf(i,:,:,:,:),       massflxbnd(i,:,:,:,:), acen(i,:,:,:,:),         acen_tf(i,:,:,:,:),        &
             rhcen(i,:,:,:,:),        qcloudcen(i,:,:,:,:),     qicecen(i,:,:,:,:),    qlsink_afcen(i,:,:,:,:),                            &
             precrcen(i,:,:,:,:),     precsolidcen(i,:,:,:,:),                                                                             &
             qlsink_bfcen(i,:,:,:,:), qlsink_avgcen(i,:,:,:,:), praincen(i,:,:,:,:),                                                       &
             wupthresh_bnd(i,:),      wdownthresh_bnd(i,:),                                                                                &
             wwqui_cen(i,:),          wwqui_bnd(i,:),           wwqui_cloudy_cen(i,:), wwqui_cloudy_bnd(i,:),                              &
#endif
             tkez(i,:),               tkesgsz(i,:),             tkz(i, :),                                                                 &
             flux_u(i,:),             flux_v(i,:),              flux_qt(i,:),          fluxsgs_qt(i,:),         flux_qp(i,:),              &
             precflux(i,:),           qt_ls(i,:),               qt_trans(i,:),         qp_trans(i,:),           qp_fall(i,:),              &
             qp_evp(i,:),             qp_src(i,:),              t_ls(i,:),             prectend(i),             precstend(i),              &
             cam_in%ocnfrac(i),       wnd,                      tau00,                 bflx,                                               & 
             fluxu0,                  fluxv0,                   fluxt0,                fluxq0,                                             & 
             taux_crm(i),             tauy_crm(i),              z0m(i),                timing_factor(i),        qtotcrm(i, :),             &
             tvwle(i,:),buoy(i,:),buoysd(i,:),msef(i,:),qvw(i,:) )   ! MDB 8/2013  
           if (SPCAM_microp_scheme .eq. 'sam1mom') then
              crm_qt(i,:,:,:) = crm_micro(i,:,:,:,1)
              crm_qp(i,:,:,:) = crm_micro(i,:,:,:,2)
              crm_qn(i,:,:,:) = crm_micro(i,:,:,:,3)
           else if (SPCAM_microp_scheme .eq. 'm2005') then 
              crm_qt(i,:,:,:) = crm_micro(i,:,:,:,1)
              crm_nc(i,:,:,:) = crm_micro(i,:,:,:,2)
              crm_qr(i,:,:,:) = crm_micro(i,:,:,:,3)
              crm_nr(i,:,:,:) = crm_micro(i,:,:,:,4)
              crm_qi(i,:,:,:) = crm_micro(i,:,:,:,5)
              crm_ni(i,:,:,:) = crm_micro(i,:,:,:,6)
              crm_qs(i,:,:,:) = crm_micro(i,:,:,:,7)
              crm_ns(i,:,:,:) = crm_micro(i,:,:,:,8)
              crm_qg(i,:,:,:) = crm_micro(i,:,:,:,9)
              crm_ng(i,:,:,:) = crm_micro(i,:,:,:,10)
              crm_qc(i,:,:,:) = crm_micro(i,:,:,:,11)
          endif
#endif
       end do ! i (loop over ncol)

       call t_stopf('crm_call')

! There is no separate convective and stratiform precip for CRM:
       precc(:ncol) = precc(:ncol) + precl(:ncol)
       precl(:ncol) = 0.
       precsc(:ncol) = precsc(:ncol) + precsl(:ncol)
       precsl(:ncol) = 0.

       call pbuf_get_field(pbuf, prec_dp_idx, prec_dp )
       call pbuf_get_field(pbuf, snow_dp_idx, snow_dp )
       call pbuf_get_field(pbuf, prec_sh_idx, prec_sh )
       call pbuf_get_field(pbuf, snow_sh_idx, snow_sh )
       call pbuf_get_field(pbuf, prec_sed_idx, prec_sed )
       call pbuf_get_field(pbuf, snow_sed_idx, snow_sed )
       call pbuf_get_field(pbuf, prec_pcw_idx, prec_pcw )
       call pbuf_get_field(pbuf, snow_pcw_idx, snow_pcw )

       prec_dp = precc
       snow_dp = precsc
       prec_sh = 0.
       snow_sh = 0. 
       prec_sed = 0.
       snow_sed = 0.
       prec_pcw = 0
       snow_pcw = 0.

!      print*,'ptend%s:',minval(ptend%s(:ncol,:))/cpair*86400.,maxval(ptend%s(:ncol,:))/cpair*86400.
!      print*,'ptend%q:',minval(ptend%q(:ncol,:,1))*86400000.,maxval(ptend%q(:ncol,:,1))*86400000.
!      print*,'ptend%qc:',minval(ptend%q(:ncol,:,ixcldliq))*86400000.,maxval(ptend%q(:ncol,:,ixcldliq))*86400000.
!      print*,'ptend%qi:',minval(ptend%q(:ncol,:,ixcldice))*86400000.,maxval(ptend%q(:ncol,:,ixcldice))*86400000.
!      print*,'t_rad',minval(t_rad(:ncol,:,:,:)),maxval(t_rad(:ncol,:,:,:))
!      print*,'qv_rad',minval(qv_rad(:ncol,:,:,:))*1000.,maxval(qv_rad(:ncol,:,:,:))*1000.
!      print*,'qc_rad',minval(qc_rad(:ncol,:,:,:))*1000.,maxval(qc_rad(:ncol,:,:,:))*1000.
!      print*,'qi_rad',minval(qi_rad(:ncol,:,:,:))*1000.,maxval(qi_rad(:ncol,:,:,:))*1000.
!      print*,'crm_qrad',minval(crm_qrad(:ncol,:,:,:))*86400.,maxval(crm_qrad(:ncol,:,:,:))*86400.
!      print*,'precc:',minval(precc(:ncol))*86400000.,maxval(precc(:ncol))*86400000.

       do m=1,crm_nz
         k = pver-m+1
         do i = 1,ncol
           crm_qrad(i,:,:,m) = crm_qrad(i,:,:,m) * state%pdel(i,k) ! for energy conservation
         end do
       end do

       call outfld('PRES    ',state%pmid ,pcols   ,lchnk   )
       call outfld('DPRES   ',state%pdel ,pcols   ,lchnk   )
!       call outfld('HEIGHT  ',state%zm   ,pcols   ,lchnk   )

       call outfld('CRM_U   ',crm_u, pcols   ,lchnk   )
       call outfld('CRM_V   ',crm_v, pcols   ,lchnk   )
       call outfld('CRM_W   ',crm_w, pcols   ,lchnk   )
       call outfld('CRM_T   ',crm_t, pcols   ,lchnk   )

       if (SPCAM_microp_scheme .eq. 'sam1mom') then
          call outfld('CRM_QV  ',(crm_qt(:,:,:,:)-qc_crm-qi_crm),pcols   ,lchnk   )
       else if (SPCAM_microp_scheme .eq. 'm2005') then 
          call outfld('CRM_QV  ',crm_qt(:,:,:,:)-qc_crm, pcols   ,lchnk   )
       endif
       call outfld('CRM_QC  ',qc_crm   ,pcols   ,lchnk   )
       call outfld('CRM_QI  ',qi_crm   ,pcols   ,lchnk   )
       call outfld('CRM_QPC ',qpc_crm  ,pcols   ,lchnk   )
       call outfld('CRM_QPI ',qpi_crm  ,pcols   ,lchnk   )
       call outfld('CRM_PREC',prec_crm       ,pcols   ,lchnk   )

#ifdef m2005
       if (SPCAM_microp_scheme .eq. 'm2005') then
! index is defined in ./crm/MICRO_M2005/microphysics.F90
! Be cautious to use them here. They are defined in crm codes, and these codes are called only 
! after the subroutine of crm is called. So they can only be used after the 'crm' subroutine. 
! some initializaiton part of crm codes should be called in the initializaation part of cam 
! in the future.
! incl, inci, ... can not be used here, for they are defined before we call them???
! +++mhwang
          call outfld('CRM_NC ',crm_nc(:, :, :, :)   ,pcols   ,lchnk   )
          call outfld('CRM_NI ',crm_ni(:, :, :, :)   ,pcols   ,lchnk   )
          call outfld('CRM_NR ',crm_nr(:, :, :, :)   ,pcols   ,lchnk   )
          call outfld('CRM_NS ',crm_ns(:, :, :, :)   ,pcols   ,lchnk   )
          call outfld('CRM_NG ',crm_ng(:, :, :, :)   ,pcols   ,lchnk   )

          call outfld('CRM_WVAR', wvar_crm, pcols, lchnk)

          call outfld('CRM_QR ',crm_qr(:, :, :, :)   ,pcols   ,lchnk   )
          call outfld('CRM_QS ',crm_qs(:, :, :, :)   ,pcols   ,lchnk   )
          call outfld('CRM_QG ',crm_qg(:, :, :, :)   ,pcols   ,lchnk   )

! hm 7/26/11, add new output
          call outfld('CRM_AUT', aut_crm, pcols, lchnk)
          call outfld('CRM_ACC', acc_crm, pcols, lchnk)
          call outfld('CRM_EVPC', evpc_crm, pcols, lchnk)
          call outfld('CRM_EVPR', evpr_crm, pcols, lchnk)
          call outfld('CRM_MLT', mlt_crm, pcols, lchnk)
          call outfld('CRM_SUB', sub_crm, pcols, lchnk)
          call outfld('CRM_DEP', dep_crm, pcols, lchnk)
          call outfld('CRM_CON', con_crm, pcols, lchnk)
! hm 8/31/11, add new output for time-mean-avg
          call outfld('A_AUT', aut_crm_a, pcols, lchnk)
          call outfld('A_ACC', acc_crm_a, pcols, lchnk)
          call outfld('A_EVPC', evpc_crm_a, pcols, lchnk)
          call outfld('A_EVPR', evpr_crm_a, pcols, lchnk)
          call outfld('A_MLT', mlt_crm_a, pcols, lchnk)
          call outfld('A_SUB', sub_crm_a, pcols, lchnk)
          call outfld('A_DEP', dep_crm_a, pcols, lchnk)
          call outfld('A_CON', con_crm_a, pcols, lchnk)
       endif ! else of if_first_step
#endif

#ifdef CLUBB_CRM
       call outfld('UP2     ', clubb_buffer(:, :, :, :, 1)   ,pcols   ,lchnk   )
       call outfld('VP2     ', clubb_buffer(:, :, :, :, 2)   ,pcols   ,lchnk   )
       call outfld('WPRTP   ', clubb_buffer(:, :, :, :, 3)   ,pcols   ,lchnk   )
       call outfld('WPTHLP  ', clubb_buffer(:, :, :, :, 4)   ,pcols   ,lchnk   )
       call outfld('WP2     ', clubb_buffer(:, :, :, :, 5)   ,pcols   ,lchnk   )
       call outfld('WP3     ', clubb_buffer(:, :, :, :, 6)   ,pcols   ,lchnk   )
       call outfld('RTP2    ', clubb_buffer(:, :, :, :, 7)   ,pcols   ,lchnk   )
       call outfld('THLP2   ', clubb_buffer(:, :, :, :, 8)   ,pcols   ,lchnk   )
       call outfld('RTPTHLP ', clubb_buffer(:, :, :, :, 9)   ,pcols   ,lchnk   )
       call outfld('UPWP    ', clubb_buffer(:, :, :, :, 10)  ,pcols   ,lchnk   )
       call outfld('VPWP    ', clubb_buffer(:, :, :, :, 11)  ,pcols   ,lchnk   )
       call outfld('CRM_CLD ', clubb_buffer(:, :, :, :, 12)  ,pcols   ,lchnk   )
#endif

       ifld = pbuf_get_index('QRL')
       call pbuf_get_field(pbuf, ifld, qrl)
       ifld = pbuf_get_index('QRS')
       call pbuf_get_field(pbuf, ifld, qrs)
       do k =1 , pver
          do i = 1, ncol
             qrs(i,k) = qrs(i,k)/state%pdel(i,k)
             qrl(i,k) = qrl(i,k)/state%pdel(i,k)
          end do
       end do

!
! add radiation tendencies to levels above CRM domain and 2 top CRM levels
! The radition tendencies in the top 4 GCM levels are set to be zero in the CRM
       ptend%s(:ncol, :pver-crm_nz+2) = qrs(:ncol,:pver-crm_nz+2)+qrl(:ncol,:pver-crm_nz+2)
   

! calculate the radiative fluxes from the radiation calculation
! This will be used to check energe conservations
!+++mhwang, 2012-02-07 (Minghuai.Wang@pnnl.gov)
       radflux(:) = 0.0_r8
       do k=1, pver
          do i=1, ncol
            radflux(i) = radflux(i) + (qrs(i,k)+qrl(i,k)) * state%pdel(i,k)/gravit
          end do
       end do
!        print*,'qrs',minval(qrs(:ncol,:))/cpair*86400.,maxval(qrs(:ncol,:))/cpair*86400.
!        print*,'qrl',minval(qrl(:ncol,:))/cpair*86400.,maxval(qrl(:ncol,:))/cpair*86400.

       ftem(:ncol,:pver) = (ptend%s(:ncol,:pver)-qrs(:ncol,:pver)-qrl(:ncol,:pver))/cpair
!        print*,'ptend%s-qrad:',minval(ftem(:ncol,:))*86400.,maxval(ftem(:ncol,:))*86400.

       call outfld('SPQRL   ',qrl/cpair      ,pcols   ,lchnk   )
       call outfld('SPQRS   ',qrs/cpair      ,pcols   ,lchnk   )

       call outfld('SPDT    ',ftem           ,pcols   ,lchnk   )
       call outfld('SPDQ    ',ptend%q(1,1,1) ,pcols   ,lchnk   )
       call outfld('SPDQC   ',ptend%q(1,1,ixcldliq) ,pcols   ,lchnk   )
       call outfld('SPDQI   ',ptend%q(1,1,ixcldice) ,pcols   ,lchnk   )
       call outfld('SPMC    ',mctot          ,pcols   ,lchnk   )
       call outfld('SPMCUP  ',mcup           ,pcols   ,lchnk   )
       call outfld('SPMCDN  ',mcdn           ,pcols   ,lchnk   )
       call outfld('SPMCUUP ',mcuup          ,pcols   ,lchnk   )
       call outfld('SPMCUDN ',mcudn          ,pcols   ,lchnk   )
       call outfld('SPQC    ',spqc         ,pcols   ,lchnk   )
       call outfld('SPQI    ',spqi         ,pcols   ,lchnk   )
       call outfld('SPQS    ',spqs         ,pcols   ,lchnk   )
       call outfld('SPQG    ',spqg         ,pcols   ,lchnk   )
       call outfld('SPQR    ',spqr         ,pcols   ,lchnk   )

       if (SPCAM_microp_scheme .eq. 'm2005') then
          call outfld('SPNC    ',spnc         ,pcols   ,lchnk   )
          call outfld('SPNI    ',spni         ,pcols   ,lchnk   )
          call outfld('SPNS    ',spns         ,pcols   ,lchnk   )
          call outfld('SPNG    ',spng         ,pcols   ,lchnk   )
          call outfld('SPNR    ',spnr         ,pcols   ,lchnk   )
       endif

       call outfld('SPQTFLX ',flux_qt        ,pcols   ,lchnk   )
       call outfld('SPUFLX  ',flux_u         ,pcols   ,lchnk   )
       call outfld('SPVFLX  ',flux_v         ,pcols   ,lchnk   )
       call outfld('SPTKE   ',tkez           ,pcols   ,lchnk   )
       call outfld('SPTKES  ',tkesgsz        ,pcols   ,lchnk   )
       call outfld('SPTK    ',tkz            ,pcols   ,lchnk   )
       call outfld('SPQTFLXS',fluxsgs_qt     ,pcols   ,lchnk   )
       call outfld('SPQPFLX ',flux_qp        ,pcols   ,lchnk   )
       call outfld('SPPFLX  ',precflux       ,pcols   ,lchnk   )
       call outfld('SPQTLS  ',qt_ls          ,pcols   ,lchnk   )
       call outfld('SPQTTR  ',qt_trans       ,pcols   ,lchnk   )
       call outfld('SPQPTR  ',qp_trans       ,pcols   ,lchnk   )
       call outfld('SPQPEVP ',qp_evp         ,pcols   ,lchnk   )
       call outfld('SPQPFALL',qp_fall        ,pcols   ,lchnk   )
       call outfld('SPQPSRC ',qp_src         ,pcols   ,lchnk   )
       call outfld('SPTLS   ',t_ls           ,pcols   ,lchnk   )

!-- MDB 8/2013
       call outfld('SPTVFLUX',tvwle          ,pcols   ,lchnk   )
       call outfld('SPBUOY  ',buoy           ,pcols   ,lchnk   )
       call outfld('SPBUOYSD',buoysd         ,pcols   ,lchnk   )
       call outfld('SPMSEF  ',msef           ,pcols   ,lchnk   )
       call outfld('SPQVFLUX',qvw            ,pcols   ,lchnk   )

       call outfld('CLOUD   ',cld,  pcols,lchnk)
       call outfld('CLDTOT  ',cltot  ,pcols,lchnk)
       call outfld('CLDHGH  ',clhgh  ,pcols,lchnk)
       call outfld('CLDMED  ',clmed  ,pcols,lchnk)
       call outfld('CLDLOW  ',cllow  ,pcols,lchnk)
       call outfld('CLOUDTOP',cldtop, pcols,lchnk)

!       call outfld('Z0M     ',z0m  ,pcols,lchnk)
!       call outfld('TAUX_CRM',taux_crm  ,pcols,lchnk)
!       call outfld('TAUY_CRM',tauy_crm  ,pcols,lchnk)

       call outfld('TIMINGF ',timing_factor  ,pcols,lchnk)
!
! Compute liquid water paths (for diagnostics only)

       tgicewp(:ncol) = 0.
       tgliqwp(:ncol) = 0.
       do k=1,pver
          do i = 1,ncol
             cicewp(i,k) = gicewp(i,k) * 1.0e-3 / max(0.01_r8,cld(i,k)) ! In-cloud ice water path.  g/m2 --> kg/m2
             cliqwp(i,k) = gliqwp(i,k) * 1.0e-3 / max(0.01_r8,cld(i,k)) ! In-cloud liquid water path. g/m2 --> kg/m2
             tgicewp(i)  = tgicewp(i) + gicewp(i,k) *1.0e-3 ! grid cell mean ice water path.  g/m2 --> kg/m2
             tgliqwp(i)  = tgliqwp(i) + gliqwp(i,k) *1.0e-3 ! grid cell mean ice water path.  g/m2 --> kg/m2
          end do
       end do
       tgwp(:ncol) = tgicewp(:ncol) + tgliqwp(:ncol)
       gwp(:ncol,:pver) = gicewp(:ncol,:pver) + gliqwp(:ncol,:pver)
       cwp(:ncol,:pver) = cicewp(:ncol,:pver) + cliqwp(:ncol,:pver)

       call outfld('GCLDLWP' ,gwp    , pcols,lchnk)
       call outfld('TGCLDCWP',tgwp   , pcols,lchnk)
       call outfld('TGCLDLWP',tgliqwp, pcols,lchnk)
       call outfld('TGCLDIWP',tgicewp, pcols,lchnk)
       call outfld('ICLDTWP' ,cwp    , pcols,lchnk)
       call outfld('ICLDIWP' ,cicewp , pcols,lchnk)

       if (use_ECPP) then
! total clouds and precipiation 
!         ifld = pbuf_get_fld_idx('AST')
!         pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim) = cld(:ncol, :pver)
!         ifld = pbuf_get_fld_idx('QME')
!         pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1) = 0.0_r8
!         ifld = pbuf_get_fld_idx('PRAIN')
!         pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1) = qp_src(:ncol, :pver)
!         ifld = pbuf_get_fld_idx('NEVAPR')
!         pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1) = qp_evp(:ncol, :pver)

! set convective rain to be zero for PRAIN already includes precipitation from convection. 
!         ifld = pbuf_get_fld_idx( 'RPRDTOT' )
!         pbuf(ifld)%fld_ptr(1,1:ncol,1:pver,lchnk,1) = 0.0_r8  
!         ifld = pbuf_get_fld_idx( 'RPRDDP' )
!         pbuf(ifld)%fld_ptr(1,1:ncol,1:pver,lchnk,1) = 0.0_r8
!         ifld = pbuf_get_fld_idx( 'RPRDSH' )
!         pbuf(ifld)%fld_ptr(1,1:ncol,1:pver,lchnk,1) = 0.0_r8
!         ifld = pbuf_get_fld_idx( 'ICWMRDP' )
!         pbuf(ifld)%fld_ptr(1,1:ncol,1:pver,lchnk,1) = 0.0_r8
!         ifld = pbuf_get_fld_idx( 'ICWMRSH' )
!         pbuf(ifld)%fld_ptr(1,1:ncol,1:pver,lchnk,1) = 0.0_r8

! turbulence
          allocate(tempPtr(pcols,pver))
          ifld = pbuf_get_index('TKE_CRM')
          cs(:ncol, 1:pver) = state%pmid(:ncol, 1:pver)/(287.15*state%t(:ncol, 1:pver))
          tempPtr(:ncol, 1:pver) = tkez(:ncol, 1:pver)/cs(:ncol, 1:pver)
          call pbuf_set_field(pbuf, ifld, tempPtr, start=(/1,1/), kount=(/pcols, pver/) )
          deallocate(tempPtr)

          ifld = pbuf_get_index('TK_CRM')
          call pbuf_set_field(pbuf, ifld, tkz, start=(/1,1/), kount=(/pcols, pver/) )

! For convective transport
          do i=1, ncol
           ideep_crm(i) = i*1.0 
          end do
          ifld = pbuf_get_index( 'MU_CRM' )
          call pbuf_set_field(pbuf, ifld, mu_crm, start=(/1,1/), kount=(/pcols, pver/) )
          ifld = pbuf_get_index( 'MD_CRM' )
          call pbuf_set_field(pbuf, ifld, md_crm, start=(/1,1/), kount=(/pcols, pver/) )
          ifld = pbuf_get_index( 'EU_CRM' )
          call pbuf_set_field(pbuf, ifld, eu_crm, start=(/1,1/), kount=(/pcols, pver/) )
          ifld = pbuf_get_index( 'DU_CRM' )
          call pbuf_set_field(pbuf, ifld, du_crm, start=(/1,1/), kount=(/pcols, pver/) )
          ifld = pbuf_get_index( 'ED_CRM' )
          call pbuf_set_field(pbuf, ifld, eu_crm, start=(/1,1/), kount=(/pcols, pver/) )
          ifld = pbuf_get_index( 'JT_CRM' )
          call pbuf_set_field(pbuf, ifld, jt_crm, start=(/1/), kount=(/pcols/) )
          ifld = pbuf_get_index( 'MX_CRM' )
          call pbuf_set_field(pbuf, ifld, mx_crm, start=(/1/), kount=(/pcols/) )
          ifld = pbuf_get_index( 'IDEEP_CRM' )
          call pbuf_set_field(pbuf, ifld, ideep_crm, start=(/1/), kount=(/pcols/) )
       endif
       call outfld('MU_CRM  ', mu_crm, pcols, lchnk)
       call outfld('MD_CRM  ', md_crm, pcols, lchnk)
       call outfld('EU_CRM  ', eu_crm, pcols, lchnk)
       call outfld('DU_CRM  ', du_crm, pcols, lchnk)
       call outfld('ED_CRM  ', ed_crm, pcols, lchnk)

#ifdef ECPP
       if (use_ECPP) then

         qlsinkcen = qlsink_avgcen

         call outfld('ACEN    ', acen, pcols, lchnk)
         call outfld('ABND    ', abnd, pcols, lchnk)
         call outfld('ACEN_TF ', acen_tf, pcols, lchnk)
         call outfld('ABND_TF ', abnd_tf, pcols, lchnk)
         call outfld('MASFBND ', massflxbnd, pcols, lchnk)
         call outfld('RHCEN   ', rhcen, pcols, lchnk)
         call outfld('QCCEN   ', qcloudcen, pcols, lchnk)
         call outfld('QICEN   ', qicecen, pcols, lchnk)
         call outfld('QSINK_AFCEN', qlsink_afcen, pcols, lchnk)
         call outfld('PRECRCEN', precrcen, pcols, lchnk)
         call outfld('PRECSCEN', precsolidcen, pcols, lchnk)
         call outfld('WUPTHRES', wupthresh_bnd, pcols, lchnk)
         call outfld('WDNTHRES', wdownthresh_bnd, pcols, lchnk)
         call outfld('WWQUI_CEN', wwqui_cen, pcols, lchnk)
         call outfld('WWQUI_CLD_CEN', wwqui_cloudy_cen, pcols, lchnk)
         call outfld('WWQUI_BND', wwqui_cen, pcols, lchnk)
         call outfld('WWQUI_CLD_BND', wwqui_cloudy_cen, pcols, lchnk)
         call outfld('QSINK_BFCEN', qlsink_bfcen, pcols, lchnk)
         call outfld('QSINK_AVGCEN', qlsink_avgcen, pcols, lchnk)
         call outfld('PRAINCEN', praincen, pcols, lchnk)
       endif !/*ECPP*/
#endif

       ptend%name  = 'crm'
       ptend%ls    = .TRUE.
       ptend%lq(1) = .TRUE.
       ptend%lq(ixcldliq) = .TRUE.
       ptend%lq(ixcldice) = .TRUE.
       ptend%lu    = .FALSE.
       ptend%lv    = .FALSE.
#ifdef CRM3D
       ptend%lu    = .TRUE.
       ptend%lv    = .TRUE.
       call outfld('UCONVMOM',ptend%u,pcols   ,lchnk   )
       call outfld('VCONVMOM',ptend%v,pcols   ,lchnk   )
#endif

       call phys_getopts(microp_scheme_out=microp_scheme)
       if(microp_scheme .eq. 'MG' ) then
#ifdef m2005
         if (SPCAM_microp_scheme .eq. 'm2005') then
            call cnst_get_ind('NUMLIQ', ixnumliq)
            call cnst_get_ind('NUMICE', ixnumice)
            ptend%lq(ixnumliq) = .TRUE.
            ptend%lq(ixnumice) = .TRUE.
            ptend%q(:, :, ixnumliq) = 0._r8
            ptend%q(:, :, ixnumice) = 0._r8

            do i = 1, ncol
             do k=1, crm_nz 
               m= pver-k+1
               do ii=1, crm_nx
               do jj=1, crm_ny
                 ptend%q(i,m,ixnumliq) = ptend%q(i,m,ixnumliq) + crm_nc(i,ii,jj,k) 
                 ptend%q(i,m,ixnumice) = ptend%q(i,m,ixnumice) + crm_ni(i,ii,jj,k)
               end do
               end do
               ptend%q(i,m,ixnumliq) = (ptend%q(i,m,ixnumliq)/(crm_nx*crm_ny) - state%q(i,m,ixnumliq))/ztodt
               ptend%q(i,m,ixnumice) = (ptend%q(i,m,ixnumice)/(crm_nx*crm_ny) - state%q(i,m,ixnumice))/ztodt
             end do
            end do
         endif
#endif
       end if

! physics_update should be moved to tphysbc. Also, it seems, state is changed in the subroutine of crm. 
! so doese this mean state is updated twice??? -Minghuai Wang (minghuai@pnl.gov).
       !!!call physics_update(state, tend, ptend, ztodt)
       call physics_update(state, ptend, ztodt, tend)

       if(SPCAM_microp_scheme .eq. 'm2005') then
! calculate column water of rain, snow and graupel 
         do i=1, ncol
            qt_hydro(i, 2) = 0.0_r8
            qli_hydro(i, 2) = 0.0_r8
            qi_hydro(i, 2) = 0.0_r8
            qtot(i, 3) = 0.0_r8
            qt_cloud(i, 3) = 0.0_r8
            qtv(i, 3) = 0.0_r8
            do m=1, crm_nz
              k=pver-m+1
              do ii=1, crm_nx
              do jj=1, crm_ny
                 qt_hydro(i,2) = qt_hydro(i,2)+(crm_qr(i,ii,jj,m)+crm_qs(i,ii,jj,m)+crm_qg(i,ii,jj,m)) * state%pdel(i,k)/gravit 
                 qli_hydro(i,2) = qli_hydro(i,2)+(crm_qr(i,ii,jj,m)+crm_qs(i,ii,jj,m)+crm_qg(i,ii,jj,m)) * state%pdel(i,k)/gravit
                 qi_hydro(i,2) = qi_hydro(i,2)+(crm_qs(i,ii,jj,m)+crm_qg(i,ii,jj,m)) * state%pdel(i,k)/gravit
                 qtot(i, 3) = qtot(i,3)+(crm_qr(i,ii,jj,m)+crm_qs(i,ii,jj,m)+crm_qg(i,ii,jj,m)) * state%pdel(i,k)/gravit     &
                      +(crm_qt(i,ii,jj,m)+crm_qi(i,ii,jj,m)) * state%pdel(i,k)/gravit
                 qt_cloud(i, 3) = qt_cloud(i, 3) + (crm_qt(i,ii,jj,m)+crm_qi(i,ii,jj,m)) * state%pdel(i,k)/gravit 
              end do
              end do
            end do
            qt_hydro(i,2) = qt_hydro(i,2)/(crm_nx*crm_ny)
            qli_hydro(i,2) = qli_hydro(i,2)/(crm_nx*crm_ny)
            qi_hydro(i,2) = qi_hydro(i,2)/(crm_nx*crm_ny)
            qtot(i, 3) = qtot(i, 3)/(crm_nx*crm_ny)
            qt_cloud(i, 3) = qt_cloud(i, 3)/(crm_nx*crm_ny)
         end do 
       else if(SPCAM_microp_scheme .eq. 'sam1mom') then 
         do i=1, ncol
            qli_hydro(i, 2) = 0.0_r8
            qi_hydro(i, 2) = 0.0_r8
            do m=1, crm_nz
              k=pver-m+1
              do ii=1, crm_nx
              do jj=1, crm_ny
                sfactor = max(0._r8,min(1._r8,(crm_t(i,ii,jj,m)-268.16)*1./(283.16-268.16)))
                qli_hydro(i,2) = qli_hydro(i,2)+crm_qp(i,ii,jj,m) * state%pdel(i,k)/gravit
                qi_hydro(i,2) = qi_hydro(i,2)+crm_qp(i,ii,jj,m) * (1-sfactor) * state%pdel(i,k)/gravit
              end do
              end do
            end do
            qli_hydro(i,2) = qli_hydro(i,2)/(crm_nx*crm_ny)
            qi_hydro(i,2) = qi_hydro(i,2)/(crm_nx*crm_ny)

! total cloud water and total water vapor, and energy
            qt_cloud(i,2) = 0._r8
            ql_cloud(i,2) = 0._r8
            qtv(i,2) = 0._r8
            do k=1, pver
              qt_cloud(i,2) = qt_cloud(i,2) + (state%q(i,k,ixcldliq)+state%q(i,k,ixcldice)) * state%pdel(i,k)/gravit
              ql_cloud(i,2) = ql_cloud(i,2) + state%q(i,k,ixcldliq) * state%pdel(i,k)/gravit
              qtv(i,2) = qtv(i,2) + state%q(i,k,1) * state%pdel(i,k)/gravit
            end do
         end do
       end if

! check water and energy conservation
       call check_energy_chng(state, tend, "crm_tend", nstep, ztodt, zero, prec_dp+(qli_hydro(:,2)-qli_hydro(:,1))/ztodt/1000._r8,  &
                snow_dp+(qi_hydro(:,2)-qi_hydro(:,1))/ztodt/1000., radflux)
!         call check_energy_chng(state, tend, "crm_tend", nstep, ztodt, zero, prec_dp,  &
!                snow_dp, zero)
!         call check_energy_chng(state, tend, "crm_tend", nstep, ztodt, zero, zero, zero, zero)

!+++mhwangtest
!
! calculate total water after crm update
! total hydrometer water (rain, snow, and graupel)
!-- mark branson: insert ifdef m2005 around this test so 1-moment micro compile will work
#ifdef m2005
       if (SPCAM_microp_scheme .eq. 'm2005') then
          do i=1, ncol

! total cloud water and total water vapor
             qt_cloud(i,2) = 0._r8
             qtv(i,2) = 0._r8
             do k=1, pver
               qt_cloud(i,2) = qt_cloud(i,2) + (state%q(i,k,ixcldliq)+state%q(i,k,ixcldice)) * state%pdel(i,k)/gravit
               qtv(i,2) = qtv(i,2) + state%q(i,k,1) * state%pdel(i,k)/gravit
             end do
! total water
             qtot(i,2) = qt_hydro(i,2) + qt_cloud(i,2) + qtv(i,2)

! to check water conservations
             if(abs((qtot(i,2)+(precc(i)+precl(i))*1000*ztodt)-qtot(i,1))/qtot(i,1).gt.1.0e-5) then 
                write(0, *) 'water before crm call', i, lchnk, qtot(i,1), qtv(i,1), qt_cloud(i,1), qt_hydro(i,1)
                write(0, *) 'water after crm call', i, lchnk, qtot(i,2)+(precc(i)+precl(i))*1000*ztodt, qtv(i,2), qt_cloud(i,2), qt_hydro(i,2), (precc(i)+precl(i))*1000*ztodt
                write(0, *) 'water, nstep, crm call2', nstep, i, lchnk, ((qtot(i,2)+(precc(i)+precl(i))*1000*ztodt)-qtot(i,1))/qtot(i,1)
                write(0, *) 'water, calcualted in crm.F90', i, lchnk, qtotcrm(i, 1), qtotcrm(i, 9), qtot(i, 3)+(precc(i)+precl(i))*1000*ztodt, qt_cloud(i, 3), qtv(i,2)+qt_cloud(i,2)
                write(0, *) 'water, temperature', i, lchnk, state%t(i,pver)
!               call endrun('water conservation error in crm_physics')
             end if
          end do ! end i
!---mhwangtest
       endif
#endif

   end if ! (is_first_step())

! ***  OLD COMMENT ******
! physics_update is called in crm_physics. This is not a good practice. But it seems state is changed in
! crm_physics. The physics_udpate should be moved to here in the future. ---Minghuai Wang (Minghuai.Wang@pnl.gov)
!   call physics_update(state, tend, ptend, ztodt)

! should the energy check are applied???? -Minghuai Wang
!   call check_energy_chng(state, tend, "crm_physics", nstep, ztodt, zero, zero, zero, zero)

   call t_stopf('crm')

!-- mark branson: insert ifdef m2005 block so that 1-moment microphysics will compile
   if (SPCAM_microp_scheme .eq. 'm2005') then
       call t_startf('bc_aerosols_mmf')
!      calculate aerosol water at CRM domain using water vapor at CRM domain +++mhwang
!
       do i=1, ncol
       do ii=1, crm_nx
        do jj=1, crm_ny
         do m=1, crm_nz
          if(qc_rad(i,ii,jj,m)+qi_rad(i,ii,jj,m).gt.1.0e-10) then
            cld_rad(i,ii,jj,m) = 1.0_r8
          else
            cld_rad(i,ii,jj,m) = 0.0_r8
          endif
         enddo
        enddo
       enddo
       enddo
      !!!call aerosol_wet_intr (state, ptend, ztodt, pbuf, cam_out, dlf)

! When ECPP is included, wet deposition is done ECPP,
! So tendency from wet depostion is not updated in mz_aero_wet_intr (mz_aerosols_intr.F90)
! tendency from other parts of crmclouds_aerosol_wet_intr are still updated here.
      !!!call physics_update (state, tend, ptend, ztodt)
      call physics_update (state, ptend, ztodt, tend)


#ifdef ECPP
      if (use_ECPP) then

          pblh_idx  = pbuf_get_index('pblh')
          call pbuf_get_field(pbuf, pblh_idx, pblh)

!
! cldo and cldn are set to be the same in crmclouds_mixnuc_tend,
! So only turbulence mixing is done here.
!
!           call t_startf('crmclouds_mixnuc')
!           call crmclouds_mixnuc_tend (state, ptend, ztodt, cam_in%cflx, pblh, pbuf,  &
!                 wwqui_cen, wwqui_cloudy_cen, wwqui_bnd, wwqui_cloudy_bnd)
!           call physics_update(state, tend, ptend, ztodt)
!           call t_stopf('crmclouds_mixnuc')

!
!   ECPP is called at every 3rd GCM time step.
!   GCM time step is 10 minutes, and ECPP time step is 30 minutes.
!
          dtstep_pp = dtstep_pp_input
          necpp = dtstep_pp/ztodt
          if(nstep.ne.0 .and. mod(nstep, necpp).eq.0) then
            call t_startf('crmclouds_mixnuc')
            call crmclouds_mixnuc_tend (state, ptend, dtstep_pp, cam_in%cflx, pblh, pbuf,  &
                        wwqui_cen, wwqui_cloudy_cen, wwqui_bnd, wwqui_cloudy_bnd)
            !!!call physics_update(state, tend, ptend, dtstep_pp)
            call physics_update(state, ptend, dtstep_pp, tend)
            call t_stopf('crmclouds_mixnuc')

            call t_startf('ecpp')
            call parampollu_driver2(state, ptend, pbuf, dtstep_pp, dtstep_pp,  &
               acen, abnd, acen_tf, abnd_tf, massflxbnd,   &
               rhcen, qcloudcen, qlsinkcen, precrcen, precsolidcen, acldy_cen_tbeg )
            !!!call physics_update(state, tend, ptend, dtstep_pp)
            call physics_update(state, ptend, dtstep_pp, tend)
            call t_stopf ('ecpp')
          end if
      endif ! /*ECPP*/
#endif

      call t_stopf('bc_aerosols_mmf')
   endif ! /*m2005*/

!### mdb spcam:  commenting this out for now just to get sam1mom working
!   if (use_SPCAM .and. use_ECPP .and. SPCAM_microp_scheme .eq.  'm2005') then
!!
!!  call deep convective transport subroutine to get the tendency for testing purpose.
!     if ( .not. deep_scheme_does_scav_trans() ) then
!
!       call convect_deep_tend_2( state,   ptend,  ztodt,  pbuf )
!       do m=1, pcnst
!         if(cnst_name(m) == 'DMS') call outfld(trim(cnst_name(m))//'CONV', ptend%q(:,:,m), pcols, lchnk )
!         if(cnst_name(m) == 'SO2') call outfld(trim(cnst_name(m))//'CONV', ptend%q(:,:,m), pcols, lchnk )
!       end do
!!
!! for testing purpose, q is not updated.
!       ptend%lq(:) = .false.
!       call physics_update (state, tend, ptend, ztodt)
!
!     endif ! deep_scheme
!   endif ! /*use_SPCAM && use_ECPP* && *m2005*/

! save for old cloud fraction in the MMF simulations
! In the CAM model, this is done in cldwat2m.F90.
!
   cldo(:ncol, :) = cld(:ncol, :)


end subroutine crm_physics_tend
!--------------------------------------------------------------------------------------------

!=====================================================================================================
subroutine crm_save_state_tend(state,tend,pbuf)
!-----------------------------------------------------------------------------------------------------
! This subroutine is used to save state, tend, cldo, qqcw and dgnumwet so that they can be used
! after they have been changed by conventional physics.
!-----------------------------------------------------------------------------------------------------

  use infnan,          only: inf
  use physics_types,   only: physics_state, physics_tend
  use time_manager,    only: is_first_step
  use physics_buffer,  only: pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field, physics_buffer_desc
  use phys_control,    only: phys_getopts
#ifdef MODAL_AERO
  use modal_aero_data, only: ntot_amode, qqcw_get_field
#endif


  implicit none

  type(physics_state),    intent(in) :: state
  type(physics_tend),     intent(in) :: tend
  type(physics_buffer_desc), pointer :: pbuf(:)

  integer itim, ifld, ncol, i, lchnk
  real(r8), pointer, dimension(:,:) :: cld        ! cloud fraction
  real(r8), pointer, dimension(:,:) :: qqcw
  logical                           :: use_SPCAM, use_ECPP


   call phys_getopts( use_SPCAM_out           = use_SPCAM )
   call phys_getopts( use_ECPP_out            = use_ECPP)

   ncol  = state%ncol

   if (use_SPCAM) then
      itim = pbuf_old_tim_idx()
      call pbuf_get_field(pbuf, cldo_idx, cldo, start=(/1,1,itim/), kount=(/pcols,pver,1/) )
      if (use_ECPP) then
         ifld = pbuf_get_index('CLD')
         call pbuf_get_field(pbuf, ifld, cld, (/1,1,itim/),(/pcols,pver,1/))
         ifld = pbuf_get_index('ACLDY_CEN')
         call pbuf_get_field(pbuf, ifld, acldy_cen_tbeg)
         if(is_first_step())then
           acldy_cen_tbeg(:ncol,:) = cld(:ncol, :)
         end if
      endif
   endif ! /*use_SPCAM*/

   
   if (use_SPCAM) then
! Save the state and tend variables to overwrite conventional physics effects
! leter before calling the superparameterization. Conventional moist
! physics is allowed to compute tendencies due to conventional
! moist physics for diagnostics purposes. -Marat

       state_save = state
       !tend_save = tend            ! mdb: try uncommenting this line
       
!-- mdb spcam: replace these lines with call to physics_tend_alloc
       !tend_save%dtdt    = inf
       !tend_save%dudt    = inf
       !tend_save%dvdt    = inf
       !tend_save%te_tnd  = inf
       !tend_save%tw_tnd  = inf
       !tend_save%flx_net = inf
       !call physics_tend_alloc(tend_save, pcols)


       tend_save%dtdt(:ncol,:pver) = tend%dtdt(:ncol,:pver)
       tend_save%dudt(:ncol,:pver) = tend%dudt(:ncol,:pver)
       tend_save%dvdt(:ncol,:pver) = tend%dvdt(:ncol,:pver)
       tend_save%te_tnd(:ncol)     = tend%te_tnd(:ncol)
       tend_save%tw_tnd(:ncol)     = tend%tw_tnd(:ncol)
       tend_save%flx_net(:ncol)    = tend_save%flx_net(:ncol)

       cldo_save(:ncol, :) = cldo(:ncol, :)

#if ( defined MODAL_AERO )
       qqcw_all=0_r8
       do i=1,pcnst
          qqcw   =>  qqcw_get_field(pbuf, i,lchnk, .true.)
          if (associated(qqcw)) qqcw_all(:,:,i) = qqcw(:,:)
       end do

       ifld = pbuf_get_index('DGNUMWET')
       call pbuf_get_field(pbuf, ifld, dgnumwet, start=(/1,1,1/), kount=(/pcols,pver,ntot_amode/))

       if(is_first_step())then
         do i=1,pcnst
            qqcw_all(:,:,i) = 1.e-38_r8
         end do
         dgnumwet(1:pcols,1:pver,1:ntot_amode) = 0.0_r8
       endif
   
       qqcw_save = qqcw_all
       dgnumwet_save = dgnumwet
#endif

   endif ! use_SPCAM

end subroutine crm_save_state_tend
!-----------------------------------------------------------------------------------------------------

!=====================================================================================================
subroutine m2005_effradius(ql, nl,qi,ni,qs, ns, cld, pres, tk, effl, effi, effl_fn, deffi, lamcrad, pgamrad, des)
!-----------------------------------------------------------------------------------------------------
! 
! This subroutine is used to calculate droplet and ice crystal effective radius, which will be used 
! in the CAM radiation code. The method to calcualte effective radius is taken out of the Morrision's
! two momenent scheme from M2005MICRO_GRAUPEL. It is also very similar with the subroutine of effradius in 
! the module of cldwat2m in the CAM source codes. 
!
! Adopted by Minghuai Wang (Minghuai.Wang@pnl.gov). 
!
!-----------------------------------------------------------------------------------------------------
   ! ----------------------------------------------------------- !
   ! Calculate effective radius for pass to radiation code       !
   ! If no cloud water, default value is 10 micron for droplets, !
   ! 25 micron for cloud ice.                                    !
   ! Be careful of the unit of effective radius : [micro meter]  !
   ! ----------------------------------------------------------- !
  !!!use gamma_function, only: gamma
  use shr_spfn_mod,    only: gamma => shr_spfn_gamma
  implicit none

     real(r8), intent(in)    :: ql                ! Mean LWC of pixels [ kg/kg ]
     real(r8), intent(in)    :: nl                ! Grid-mean number concentration of cloud liquid droplet [#/kg]
     real(r8), intent(in)    :: qi                ! Mean IWC of pixels [ kg/kg ]
     real(r8), intent(in)    :: ni                ! Grid-mean number concentration of cloud ice    droplet [#/kg]
     real(r8), intent(in)    :: qs                ! mean snow water content [kg/kg]
     real(r8), intent(in)    :: ns                ! Mean snow crystal number concnetration [#/kg]
     real(r8), intent(in)    :: cld               ! Physical stratus fraction
     real(r8), intent(in)    :: pres               ! Air pressure [Pa] 
     real(r8), intent(in)    :: tk                ! air temperature [K]

     real(r8), intent(out)   :: effl              ! Effective radius of cloud liquid droplet [micro-meter]
     real(r8), intent(out)   :: effi              ! Effective radius of cloud ice    droplet [micro-meter]
     real(r8), intent(out)   :: effl_fn           ! effl for fixed number concentration of nlic = 1.e8
     real(r8), intent(out)   :: deffi             ! ice effective diameter for optics (radiation)
     real(r8), intent(out)   :: pgamrad           ! gamma parameter for optics (radiation)
     real(r8), intent(out)   :: lamcrad           ! slope of droplet distribution for optics (radiation)
     real(r8), intent(out)   :: des               ! snow effective diameter for optics (radiation) [micro-meter]

     real(r8)  qlic                               ! In-cloud LWC [kg/m3]
     real(r8)  qiic                               ! In-cloud IWC [kg/m3]
     real(r8)  nlic                               ! In-cloud liquid number concentration [#/kg]
     real(r8)  niic                               ! In-cloud ice    number concentration [#/kg]

     real(r8)  mtime                              ! Factor to account for droplet activation timescale [no]
     real(r8)  cldm                               ! Constrained stratus fraction [no]
     real(r8)  mincld                             ! Minimum stratus fraction [no]

     real(r8)  lami, laml, lammax, lammin, pgam, lams, lammaxs, lammins

     real(r8)  dcs    !autoconversion size threshold   [meter]
     real(r8)  di, ci     ! cloud ice mass-diameter relationship
     real(r8)  ds, cs     ! snow crystal mass-diameter relationship 
     real(r8)  qsmall
     real(r8)  rho       ! air density [kg/m3]
     real(r8)  rhow      ! liquid water density [kg/m3]
     real(r8)  rhoi      ! ice density [kg/m3]
     real(r8)  rhos      ! snow density [kg/m3]
     real(r8)  res       ! effective snow diameters
     real(r8)  pi
     real(r8)  tempnc

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

     pi = 3.1415926535897932384626434
!     qsmall = 1.0e-18  ! in the CAM source code (cldwat2m)
     qsmall = 1.0e-14  ! in the SAM source code (module_mp_graupel)
!     rhow = 1000.      ! in cldwat2m, CAM 
     rhow = 997.       ! in module_mp_graupel, SAM
     rhoi = 500.       ! in both CAM and SAM

!     dcs = 70.e-6_r8    ! in cldwat2m, CAM 
     dcs = 125.e-6_r8   ! in module_mp_graupel, SAM 
     ci = rhoi * pi/6.
     di = 3.

!   for snow water
     rhos = 100.      ! in both SAM and CAM5 
     cs = rhos*pi/6.
     ds = 3.


     rho = pres / (287.15*tk)    ! air density [kg/m3]

     mincld  = 0.0001_r8
     cldm    = max(cld,mincld)
     qlic    = min(5.e-3_r8,max(0._r8,ql/cldm))
     qiic    = min(5.e-3_r8,max(0._r8,qi/cldm))
     nlic    = max(nl,0._r8)/cldm
     niic    = max(ni,0._r8)/cldm

!------------------------------------------------------
!     Effective diameters of snow crystals
!------------------------------------------------------
     if(qs.gt.1.0e-7) then 
       lammaxs=1._r8/10.e-6_r8
       lammins=1._r8/2000.e-6_r8
       lams = (gamma(1._r8+ds)*cs * ns/qs)**(1._r8/ds)
       lams = min(lammaxs,max(lams,lammins))
       res = 1.5/lams*1.0e6_r8
     else
       res = 500._r8 
     end if 
!
! from Hugh Morrision: rhos/917 accouts for assumptions about 
! ice density in the Mitchell optics. 
!
     des = res * rhos/917._r8 *2._r8

   ! ------------------------------------- !
   ! Effective radius of cloud ice droplet !
   ! ------------------------------------- !

     if( qiic.ge.qsmall ) then
         niic   = min(niic,qiic*1.e20_r8)
!         lammax = 1._r8/10.e-6_r8      ! in cldwat2m, CAM
!         lammin = 1._r8/(2._r8*dcs)    ! in cldwat2m, CAM
         lammax = 1._r8/1.e-6_r8      ! in module_mp_graupel, SAM 
         lammin = 1._r8/(2._r8*dcs+100.e-6_r8)    ! in module_mp_graupel, SAM 
         lami   = (gamma(1._r8+di)*ci*niic/qiic)**(1._r8/di)
         lami   = min(lammax,max(lami,lammin))
         effi   = 1.5_r8/lami*1.e6_r8
     else
         effi   = 25._r8
     endif

!--hm ice effective radius for david mitchell's optics
!--ac morrison indicates that this is effective diameter
!--ac morrison indicates 917 (for the density of pure ice..)
     deffi  = effi *rhoi/917._r8*2._r8

   ! ---------------------------------------- !
   ! Effective radius of cloud liquid droplet !
   ! ---------------------------------------- !

     if( qlic.ge.qsmall ) then
!        Matin et al., 1994 (JAS) formula for pgam (the same is used in both CAM and SAM).
!        See also Morrison and Grabowski (2007, JAS, Eq. (2))
         nlic   = min(nlic,qlic*1.e20_r8)

! set the minimum droplet number as 20/cm3.
!         nlic   = max(nlic,20.e6_r8/rho) ! sghan minimum in #/cm3
         tempnc = nlic/rho/1.0e6    ! #/kg --> #/cm3
!         if (tempnc.gt.100._r8) then 
!           write(0, *) 'nc larger than 100  ', tempnc, rho
!         end if

!!!!!! ????? Should be the in-cloud dropelt number calculated as nlic*rho/1.0e6_r8 ????!!!! +++mhwang
!         pgam   = 0.0005714_r8*(nlic/1.e6_r8/rho) + 0.2714_r8  !wrong, confirmed with Hugh Morrison. fixed in the latest SAM. 
         pgam   = 0.0005714_r8*(nlic*rho/1.e6_r8) + 0.2714_r8
         pgam   = 1._r8/(pgam**2)-1._r8
!         pgam   = min(15._r8,max(pgam,2._r8))   ! in cldwat2m, CAM
         pgam   = min(10._r8,max(pgam,2._r8))   ! in module_mp_graupel, SAM
!         if(pgam.gt.2.01_r8 .and.pgam.lt.9.99_r8) then
!           write(0, *) 'pgam', pgam
!         end if
         laml   = (pi/6._r8*rhow*nlic*gamma(pgam+4._r8)/(qlic*gamma(pgam+1._r8)))**(1._r8/3._r8)
         lammin = (pgam+1._r8)/50.e-6_r8    ! in cldwat2m, CAM
         lammax = (pgam+1._r8)/2.e-6_r8     ! in cldwat2m, CAM   ! cldwat2m should be used, 
                                                                 ! if lammax is too large, this will lead to crash in 
                                                                 ! src/physics/rrtmg/cloud_rad_props.F90 because 
                                                                 ! klambda-1 can be zero in gam_liquid_lw and gam_liquid_sw
                                                                 !  and g_lambda(kmu,klambda-1) will not be defined. 
!         lammin = (pgam+1._r8)/60.e-6_r8    ! in module_mp_graupel, SAM
!         lammax = (pgam+1._r8)/1.e-6_r8     ! in module_mp_graupel, SAM

         laml   = min(max(laml,lammin),lammax)
!         effl   = gamma(qcvar+1._r8/3._r8)/(gamma(qcvar)*qcvar**(1._r8/3._r8))* &
!                  gamma(pgam+4._r8)/gamma(pgam+3._r8)/laml/2._r8*1.e6_r8      ! in cldwat2m, CAM
         effl   =  gamma(pgam+4._r8)/gamma(pgam+3._r8)/laml/2._r8*1.e6_r8  ! in module_mp_graupel, SAM
         lamcrad  = laml 
         pgamrad  = pgam
     else
! we chose 10. over 25, since 10 is a more reasonable value for liquid droplet. +++mhwang
         effl   = 10._r8     ! in cldwat2m, CAM
!         effl   = 25._r8     ! in module_mp_graupel, SAM
         lamcrad  = 0.0_r8
         pgamrad  = 0.0_r8
     endif

   ! ---------------------------------------------------------------------- !
   ! Recalculate effective radius for constant number, in order to separate !
   ! first and second indirect effects. Assume constant number of 10^8 kg-1 !
   ! ---------------------------------------------------------------------- !

     nlic = 1.e8
     if( qlic.ge.qsmall ) then
!        Matin et al., 1994 (JAS) formula for pgam (the same is used in both CAM and SAM). 
!        See also Morrison and Grabowski (2007, JAS, Eq. (2))  
         nlic   = min(nlic,qlic*1.e20_r8)
         pgam   = 0.0005714_r8*(nlic/1.e6_r8/rho) + 0.2714_r8
         pgam   = 1._r8/(pgam**2)-1._r8
!         pgam   = min(15._r8,max(pgam,2._r8))   ! in cldwat2m, CAM
         pgam   = min(10._r8,max(pgam,2._r8))   ! in module_mp_graupel, SAM
         laml   = (pi/6._r8*rhow*nlic*gamma(pgam+4._r8)/(qlic*gamma(pgam+1._r8)))**(1._r8/3._r8)
!         lammin = (pgam+1._r8)/50.e-6_r8    ! in cldwat2m, CAM
!         lammax = (pgam+1._r8)/2.e-6_r8     ! in cldwat2m, CAM
         lammin = (pgam+1._r8)/60.e-6_r8    ! in module_mp_graupel, SAM
         lammax = (pgam+1._r8)/1.e-6_r8     ! in module_mp_graupel, SAM

         laml   = min(max(laml,lammin),lammax)
!         effl_fn   = gamma(qcvar+1._r8/3._r8)/(gamma(qcvar)*qcvar**(1._r8/3._r8))* &
!                  gamma(pgam+4._r8)/gamma(pgam+3._r8)/laml/2._r8*1.e6_r8      ! in cldwat2m, CAM
         effl_fn   =  gamma(pgam+4._r8)/gamma(pgam+3._r8)/laml/2._r8*1.e6_r8  ! in module_mp_graupel, SAM
     else
! we chose 10. over 25, since 10 is a more reasonable value for liquid droplet. +++mhwang
         effl_fn   = 10._r8     ! in cldwat2m, CAM
!         effl_fn   = 25._r8     ! in module_mp_graupel, SAM
     endif

   return
end subroutine m2005_effradius

end module crm_physics
