
!================================================================================================
! output data necessary to drive radiation offline
! Francis Vitt -- Created 15 Dec 2009
!================================================================================================
module radiation_data

  use shr_kind_mod,     only: r8=>shr_kind_r8
  use ppgrid,           only: pcols, pver, pverp
  use cam_history,      only: addfld, horiz_only, add_default, outfld
  use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_gas, rad_cnst_get_aer_mmr
  use radconstants,     only: nradgas, gaslist
  use cam_history_support, only: fieldname_len, fillvalue
  use spmd_utils,       only: masterproc
  use cam_abortutils,       only: endrun

  implicit none
  private

  public :: output_rad_data
  public :: init_rad_data
  public :: rad_data_readnl

  integer :: cld_ifld,concld_ifld,rel_ifld,rei_ifld
  integer :: dei_ifld,mu_ifld,lambdac_ifld,iciwp_ifld,iclwp_ifld,rel_fn_ifld
  integer :: des_ifld,icswp_ifld,cldfsnow_ifld

  character(len=fieldname_len), public, parameter :: &
       lndfrc_fldn    = 'rad_lndfrc      ' , &
       icefrc_fldn    = 'rad_icefrc      ' , &
       snowh_fldn     = 'rad_snowh       ' , &
       landm_fldn     = 'rad_landm       ' , &
       asdir_fldn     = 'rad_asdir       ' , &
       asdif_fldn     = 'rad_asdif       ' , &
       aldir_fldn     = 'rad_aldir       ' , &
       aldif_fldn     = 'rad_aldif       ' , &
       coszen_fldn    = 'rad_coszen      ' , &
       asdir_pos_fldn = 'rad_asdir_pos   ' , &
       asdif_pos_fldn = 'rad_asdif_pos   ' , &
       aldir_pos_fldn = 'rad_aldir_pos   ' , &
       aldif_pos_fldn = 'rad_aldif_pos   ' , &
       lwup_fldn      = 'rad_lwup        ' , &
       ts_fldn        = 'rad_ts          ' , &
       temp_fldn      = 'rad_temp        ' , &
       pdel_fldn      = 'rad_pdel        ' , &
       pdeldry_fldn   = 'rad_pdeldry     ' , &
       pmid_fldn      = 'rad_pmid        ' , &
       watice_fldn    = 'rad_watice      ' , &
       watliq_fldn    = 'rad_watliq      ' , &
       watvap_fldn    = 'rad_watvap      ' , &
       zint_fldn      = 'rad_zint        ' , &
       pint_fldn      = 'rad_pint        ' , &
       cld_fldn       = 'rad_cld         ' , &
       cldfsnow_fldn  = 'rad_cldfsnow    ' , &
       concld_fldn    = 'rad_concld      ' , &
       rel_fldn       = 'rad_rel         ' , &
       rei_fldn       = 'rad_rei         ' , &
       dei_fldn       = 'rad_dei         ' , &
       des_fldn       = 'rad_des         ' , &
       mu_fldn        = 'rad_mu          ' , &
       lambdac_fldn   = 'rad_lambdac     ' , &
       iciwp_fldn     = 'rad_iciwp       ' , &
       iclwp_fldn     = 'rad_iclwp       ' , &
       icswp_fldn     = 'rad_icswp       '

  ! rad constituents mixing ratios
  integer :: ngas, naer
  character(len=64), allocatable :: gasnames(:)
  character(len=64), allocatable :: aernames(:)
 
  ! control options  
  logical          :: rad_data_output = .false.
  integer          :: rad_data_histfile_num = 2
  character(len=1) :: rad_data_avgflag = 'A'

  ! MG microphys check
  logical, public :: mg_microphys

contains

!================================================================================================
!================================================================================================
  subroutine rad_data_readnl(nlfile)

    ! Read rad_data_nl namelist group.  Parse input.

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr, i
    character(len=*), parameter :: subname = 'rad_data_readnl'

    namelist /rad_data_nl/ rad_data_output, rad_data_histfile_num, rad_data_avgflag

    !-----------------------------------------------------------------------------

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'rad_data_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, rad_data_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast (rad_data_output,       1,   mpilog ,  0, mpicom)
    call mpibcast (rad_data_histfile_num, 1,   mpiint ,  0, mpicom)
    call mpibcast (rad_data_avgflag,      1,   mpichar , 0, mpicom)
#endif
    
  end subroutine rad_data_readnl

  !================================================================================================
  !================================================================================================
  subroutine init_rad_data
    use phys_control,     only: phys_getopts
    use physics_buffer, only: pbuf_get_index
#ifdef MAML1
        use crmdims,       only: crm_nx,crm_ny
#endif
    
    implicit none
    
    integer :: i
    character(len=64) :: name
    character(len=128):: long_name
    character(len=64) :: long_name_description
    character(len=16) :: microp_scheme  ! microphysics scheme

    if (.not.rad_data_output) return
   
    call phys_getopts(microp_scheme_out=microp_scheme)
    mg_microphys =  (trim(microp_scheme) == 'MG')

    cld_ifld    = pbuf_get_index('CLD')
    concld_ifld = pbuf_get_index('CONCLD')
    rel_ifld    = pbuf_get_index('REL')
    rei_ifld    = pbuf_get_index('REI')
    if (mg_microphys) then
       dei_ifld      = pbuf_get_index('DEI')
       des_ifld      = pbuf_get_index('DES')
       mu_ifld       = pbuf_get_index('MU')
       lambdac_ifld  = pbuf_get_index('LAMBDAC')
       iciwp_ifld    = pbuf_get_index('ICIWP')
       iclwp_ifld    = pbuf_get_index('ICLWP')
       icswp_ifld    = pbuf_get_index('ICSWP')
       cldfsnow_ifld = pbuf_get_index('CLDFSNOW')
    endif

    call addfld (lndfrc_fldn, horiz_only,    rad_data_avgflag, 'fraction',&
         'radiation input: land fraction')
    call addfld (icefrc_fldn, horiz_only,    rad_data_avgflag, 'fraction',&
         'radiation input: ice fraction')
    call addfld (landm_fldn,     horiz_only,    rad_data_avgflag,  'none',&
         'radiation input: land mask: ocean(0), continent(1), transition(0-1)')

#ifdef MAML1    
    call addfld (snowh_fldn,    (/'crm_nx','crm_ny'/),    rad_data_avgflag, 'm', &
         'radiation input: water equivalent snow depth',  flag_xyfill=.true.)
    call addfld (asdir_fldn,    (/'crm_nx','crm_ny'/),    rad_data_avgflag, '1', &
         'radiation input: short wave direct albedo',  flag_xyfill=.true.)
    call addfld (asdif_fldn,    (/'crm_nx','crm_ny'/),    rad_data_avgflag, '1', &
         'radiation input: short wave diffuse albedo',  flag_xyfill=.true.)
    call addfld (aldir_fldn,    (/'crm_nx','crm_ny'/),    rad_data_avgflag, '1', &
         'radiation input: long wave direct albedo',  flag_xyfill=.true.)
    call addfld (aldif_fldn,    (/'crm_nx','crm_ny'/),    rad_data_avgflag, '1', &
         'radiation input: long wave diffuse albedo',  flag_xyfill=.true.)
    
    call addfld (asdir_pos_fldn,    (/'crm_nx','crm_ny'/),    rad_data_avgflag, '1', &
         'radiation input: short wave direct albedo weighted by coszen',  flag_xyfill=.true.)
    call addfld (asdif_pos_fldn,    (/'crm_nx','crm_ny'/),    rad_data_avgflag, '1', &
         'radiation input: short wave diffuse albedo weighted by coszen',  flag_xyfill=.true.)
    call addfld (aldir_pos_fldn,    (/'crm_nx','crm_ny'/),    rad_data_avgflag, '1', &
         'radiation input: long wave direct albedo weighted by coszen',  flag_xyfill=.true.)
    call addfld (aldif_pos_fldn,    (/'crm_nx','crm_ny'/),    rad_data_avgflag, '1', &
         'radiation input: long wave diffuse albedo weighted by coszen',  flag_xyfill=.true.)
    
    call addfld (lwup_fldn,    (/'crm_nx','crm_ny'/),    rad_data_avgflag, 'W/m2', &
         'radiation input: long wave up radiation flux',  flag_xyfill=.true.)
    ![lee1046] error. cam_in%ts is not CRM level field (camsrfexch.F90)
    !call addfld (ts_fldn,    (/'crm_nx','crm_ny'/),    rad_data_avgflag, 'K', &
    !     'radiation input: surface temperature',  flag_xyfill=.true.)
#else
    call addfld (snowh_fldn,        horiz_only,    rad_data_avgflag,  'm',&
         'radiation input: water equivalent snow depth')
    call addfld (asdir_fldn,        horiz_only,    rad_data_avgflag,  '1',&
         'radiation input: short wave direct albedo', flag_xyfill=.true.)
    call addfld (asdif_fldn,        horiz_only,    rad_data_avgflag,  '1',&
         'radiation input: short wave difuse albedo', flag_xyfill=.true.)
    call addfld (aldir_fldn,        horiz_only,    rad_data_avgflag,  '1',&
         'radiation input: long wave direct albedo', flag_xyfill=.true.)
    call addfld (aldif_fldn,        horiz_only,    rad_data_avgflag,  '1',&
         'radiation input: long wave difuse albedo', flag_xyfill=.true.)

    call addfld (asdir_pos_fldn, horiz_only,    rad_data_avgflag,  '1',&
         'radiation input: short wave direct albedo weighted by coszen', flag_xyfill=.true.)
    call addfld (asdif_pos_fldn, horiz_only,    rad_data_avgflag,  '1',&
         'radiation input: short wave difuse albedo weighted by coszen', flag_xyfill=.true.)
    call addfld (aldir_pos_fldn, horiz_only,    rad_data_avgflag,  '1',&
         'radiation input: long wave direct albedo weighted by coszen', flag_xyfill=.true.)
    call addfld (aldif_pos_fldn, horiz_only,    rad_data_avgflag,  '1',&
         'radiation input: long wave difuse albedo weighted by coszen', flag_xyfill=.true.)
    
    call addfld (lwup_fldn,     horiz_only,    rad_data_avgflag,   'W/m2',&
         'radiation input: long wave up radiation flux ')
#endif

    call addfld (ts_fldn,        horiz_only,    rad_data_avgflag,     'K',&
         'radiation input: surface temperature')
    call addfld (coszen_fldn, horiz_only,    rad_data_avgflag,     '1',&
         'radiation input: cosine solar zenith when positive', flag_xyfill=.true.)
    call addfld (temp_fldn,        (/ 'lev' /), rad_data_avgflag,   'K',&
         'radiation input: midpoint temperature')
    call addfld (pdel_fldn,       (/ 'lev' /), rad_data_avgflag,   'Pa',&
         'radiation input: pressure layer thickness')
    call addfld (pdeldry_fldn,       (/ 'lev' /), rad_data_avgflag,'Pa',&
         'radiation input: dry pressure layer thickness')
    call addfld (pmid_fldn,       (/ 'lev' /), rad_data_avgflag,   'Pa',&
         'radiation input: midpoint pressure')
    call addfld (watice_fldn,    (/ 'lev' /), rad_data_avgflag, 'kg/kg',&
         'radiation input: cloud ice')
    call addfld (watliq_fldn,    (/ 'lev' /), rad_data_avgflag, 'kg/kg',&
         'radiation input: cloud liquid water')
    call addfld (watvap_fldn,    (/ 'lev' /), rad_data_avgflag, 'kg/kg',&
         'radiation input: water vapor')

    call addfld (zint_fldn,       (/ 'ilev' /),rad_data_avgflag,   'km',&
         'radiation input: interface height')
    call addfld (pint_fldn,       (/ 'ilev' /),rad_data_avgflag,   'Pa',&
         'radiation input: interface pressure')

    call addfld (cld_fldn, (/ 'lev' /), rad_data_avgflag,    'fraction',&
         'radiation input: cloud fraction')
    call addfld (concld_fldn, (/ 'lev' /), rad_data_avgflag, 'fraction',&
         'radiation input: convective cloud fraction')
    call addfld (rel_fldn,   (/ 'lev' /), rad_data_avgflag,    'micron',&
         'radiation input: effective liquid drop radius')
    call addfld (rei_fldn,   (/ 'lev' /), rad_data_avgflag,    'micron',&
         'radiation input: effective ice partical radius')
    
    if (mg_microphys) then
       call addfld (dei_fldn,   (/ 'lev' /), rad_data_avgflag,    'micron',&
            'radiation input: effective ice partical diameter')
       call addfld (des_fldn,   (/ 'lev' /), rad_data_avgflag,    'micron',&
            'radiation input: effective snow partical diameter')
       call addfld (mu_fldn,        (/ 'lev' /), rad_data_avgflag,     ' ',&
            'radiation input: ice gamma parameter for optics (radiation)')
       call addfld (lambdac_fldn,        (/ 'lev' /), rad_data_avgflag,' ',&
            'radiation input: slope of droplet distribution for optics (radiation)')
       call addfld (iciwp_fldn,    (/ 'lev' /), rad_data_avgflag,  'kg/m2',&
            'radiation input: In-cloud ice water path')
       call addfld (iclwp_fldn,    (/ 'lev' /), rad_data_avgflag,  'kg/m2',&
            'radiation input: In-cloud liquid water path')
       call addfld (icswp_fldn,    (/ 'lev' /), rad_data_avgflag,  'kg/m2',&
            'radiation input: In-cloud snow water path')
       call addfld (cldfsnow_fldn, (/ 'lev' /), rad_data_avgflag, 'fraction',&
            'radiation input: cloud liquid drops + snow')
    endif

    call add_default (lndfrc_fldn,    rad_data_histfile_num, ' ')
    call add_default (icefrc_fldn,    rad_data_histfile_num, ' ')
    call add_default (snowh_fldn,     rad_data_histfile_num, ' ')
    call add_default (landm_fldn,     rad_data_histfile_num, ' ')
    call add_default (asdir_fldn,     rad_data_histfile_num, ' ')
    call add_default (asdif_fldn,     rad_data_histfile_num, ' ')
    call add_default (aldir_fldn,     rad_data_histfile_num, ' ')
    call add_default (aldif_fldn,     rad_data_histfile_num, ' ')

    call add_default (coszen_fldn,    rad_data_histfile_num, ' ')
    call add_default (asdir_pos_fldn, rad_data_histfile_num, ' ')
    call add_default (asdif_pos_fldn, rad_data_histfile_num, ' ')
    call add_default (aldir_pos_fldn, rad_data_histfile_num, ' ')
    call add_default (aldif_pos_fldn, rad_data_histfile_num, ' ')

    call add_default (lwup_fldn,      rad_data_histfile_num, ' ')
    call add_default (ts_fldn,        rad_data_histfile_num, ' ')
    call add_default (temp_fldn,      rad_data_histfile_num, ' ')
    call add_default (pdel_fldn,      rad_data_histfile_num, ' ')
    call add_default (pdeldry_fldn,   rad_data_histfile_num, ' ')
    call add_default (pmid_fldn,      rad_data_histfile_num, ' ')
    call add_default (watice_fldn,    rad_data_histfile_num, ' ')
    call add_default (watliq_fldn,    rad_data_histfile_num, ' ')
    call add_default (watvap_fldn,    rad_data_histfile_num, ' ')
    call add_default (zint_fldn,      rad_data_histfile_num, ' ')
    call add_default (pint_fldn,      rad_data_histfile_num, ' ')

    call add_default (cld_fldn,       rad_data_histfile_num, ' ')
    call add_default (concld_fldn,    rad_data_histfile_num, ' ')
    call add_default (rel_fldn,       rad_data_histfile_num, ' ')
    call add_default (rei_fldn,       rad_data_histfile_num, ' ')
    
    if (mg_microphys) then
       call add_default (dei_fldn,       rad_data_histfile_num, ' ')
       call add_default (des_fldn,       rad_data_histfile_num, ' ')
       call add_default (mu_fldn,        rad_data_histfile_num, ' ')
       call add_default (lambdac_fldn,   rad_data_histfile_num, ' ')
       call add_default (iciwp_fldn,     rad_data_histfile_num, ' ')
       call add_default (iclwp_fldn,     rad_data_histfile_num, ' ')
       call add_default (icswp_fldn,     rad_data_histfile_num, ' ')
       call add_default (cldfsnow_fldn,  rad_data_histfile_num, ' ')
    endif

    ! rad constituents

    call rad_cnst_get_info(0, ngas=ngas, naero=naer)
    long_name_description = ' mass mixing ratio used in rad climate calculation'

    ! The code to output the gases assumes that the rad_constituents module has
    ! ordered them in the same way that they are ordered in the "gaslist" array
    ! in module radconstants, and that there are nradgas of them.  This ordering 
    ! is performed in the internal init_lists routine in rad_constituents.
    if (ngas /= nradgas) then
       call endrun('init_rad_data: ERROR: ngas /= nradgas')
    end if

    allocate( gasnames(ngas) )
    call rad_cnst_get_info(0, gasnames=gasnames)

    do i = 1, ngas
       long_name = trim(gasnames(i))//trim(long_name_description)
       name = 'rad_'//gasnames(i)
       call addfld(trim(name), (/ 'lev' /), rad_data_avgflag, 'kg/kg', trim(long_name))
       call add_default (trim(name), rad_data_histfile_num, ' ')
    end do

    if (naer > 0) then
       allocate( aernames(naer) )
       call rad_cnst_get_info(0, aernames=aernames)

       do i = 1, naer
          long_name = trim(aernames(i))//trim(long_name_description)
          name = 'rad_'//aernames(i)
          call addfld(trim(name), (/ 'lev' /), rad_data_avgflag, 'kg/kg', trim(long_name))
          call add_default (trim(name), rad_data_histfile_num, ' ')
       end do
    end if

  end subroutine init_rad_data

  !================================================================================================
  !================================================================================================
  subroutine output_rad_data(  pbuf, state, cam_in, landm, coszen )

    use physics_types,    only: physics_state
    use camsrfexch,       only: cam_in_t     
    
    use constituents,     only: cnst_get_ind
    use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx
    use seq_comm_mct, only : num_inst_atm
    implicit none
    type(physics_buffer_desc), pointer :: pbuf(:)
    
    type(physics_state), intent(in), target :: state
    type(cam_in_t),      intent(in) :: cam_in
    real(r8),            intent(in) :: landm(pcols)
    real(r8),            intent(in) :: coszen(pcols)

    ! Local variables
    integer :: i, j
    character(len=32) :: name
    real(r8), pointer :: mmr(:,:)

    integer :: lchnk, itim_old, ifld
    integer :: ixcldice              ! cloud ice water index
    integer :: ixcldliq              ! cloud liquid water index
    integer :: icol
    integer :: ncol

    ! surface albedoes weighted by (positive cosine zenith angle)
    real(r8):: coszrs_pos(pcols)    ! = max(coszrs,0)
#ifdef MAML1
    real(r8):: asdir_pos (pcols,num_inst_atm)    !
    real(r8):: asdif_pos (pcols,num_inst_atm)    !
    real(r8):: aldir_pos (pcols,num_inst_atm)    !
    real(r8):: aldif_pos (pcols,num_inst_atm)    !
    integer :: j
#else
    real(r8):: snowh_avg (pcols)    ! avg of cam_in%snowhland   
    real(r8):: lwup_avg  (pcols)    ! avg of cam_in%lwup 
    real(r8):: asdir_avg (pcols)    ! cam_in%albedoes (pcols, num_inst_atm) averaged across num_inst_atm
    real(r8):: asdif_avg (pcols)    !
    real(r8):: aldir_avg (pcols)    !
    real(r8):: aldif_avg (pcols)    !
    real(r8):: asdir_pos (pcols)    ! averaged albedo (above 4) weighed by cosine of zenith angle
    real(r8):: asdif_pos (pcols)    !
    real(r8):: aldir_pos (pcols)    !
    real(r8):: aldif_pos (pcols)    !
    real(r8):: avgfac
#endif
    real(r8), pointer, dimension(:,:)  :: ptr
    
    logical           :: use_MAML ! switch to enable MAML

    call phys_getopts(use_MAML_out=use_MAML)

    if (.not.rad_data_output) return

    ! get index of (liquid+ice) cloud water
    call cnst_get_ind('CLDICE', ixcldice)
    call cnst_get_ind('CLDLIQ', ixcldliq)

    lchnk = state%lchnk
    ncol = state%ncol

#ifdef MAML1
    do icol = 1, ncol
       coszrs_pos(icol)  = max(coszen(icol),0._r8)
       do j= 1, num_inst_atm
          ![lee1046] coszrs_cos, index error. index i changed to index icol
          asdir_pos(icol,j) = cam_in%asdir(icol,j) * coszrs_pos(icol)
          asdif_pos(icol,j) = cam_in%asdif(icol,j) * coszrs_pos(icol)
          aldir_pos(icol,j) = cam_in%aldir(icol,j) * coszrs_pos(icol)
          aldif_pos(icol,j) = cam_in%aldif(icol,j) * coszrs_pos(icol)
       enddo
    enddo
#else 
    do icol = 1, ncol
       coszrs_pos(icol)  = max(coszen(icol),0._r8)
    enddo
    ! cam_in albedoes have 2 dimensions (pcols, num_inst_atm)
    ! for output purpose, these albedoes are averaged across num_inst_atm (CRM
    ! domain mean)
    avgfac = 1._r8/real(num_inst_atm,r8)
    snowh_avg(1:ncol)  = sum(cam_in%snowhland(1:ncol,1:num_inst_atm),dim=2)*avgfac
    lwup_avg (1:ncol)   = sum(cam_in%lwup (1:ncol,1:num_inst_atm),dim=2)*avgfac
    asdir_avg(1:ncol)  = sum(cam_in%asdir(1:ncol,1:num_inst_atm),dim=2)*avgfac
    asdif_avg(1:ncol)  = sum(cam_in%asdif(1:ncol,1:num_inst_atm),dim=2)*avgfac
    aldir_avg(1:ncol)  = sum(cam_in%aldir(1:ncol,1:num_inst_atm),dim=2)*avgfac
    aldif_avg(1:ncol)  = sum(cam_in%aldif(1:ncol,1:num_inst_atm),dim=2)*avgfac
    
    asdir_pos(1:ncol)  = asdir_avg(1:ncol) * coszrs_pos(1:ncol)
    asdif_pos(1:ncol)  = asdif_avg(1:ncol) * coszrs_pos(1:ncol)
    aldir_pos(1:ncol)  = aldir_avg(1:ncol) * coszrs_pos(1:ncol)
    aldif_pos(1:ncol)  = aldif_avg(1:ncol) * coszrs_pos(1:ncol)
#endif

    call outfld(lndfrc_fldn, cam_in%landfrac,  pcols, lchnk)
    call outfld(icefrc_fldn, cam_in%icefrac,   pcols, lchnk)
    call outfld(landm_fldn,  landm,            pcols, lchnk)
    call outfld(temp_fldn,   state%t,               pcols, lchnk   )
    call outfld(pdel_fldn,   state%pdel,            pcols, lchnk   )
    call outfld(pdeldry_fldn,state%pdeldry,         pcols, lchnk   )
    call outfld(watice_fldn, state%q(:,:,ixcldice), pcols, lchnk   )
    call outfld(watliq_fldn, state%q(:,:,ixcldliq), pcols, lchnk   )
    call outfld(watvap_fldn, state%q(:,:,1),        pcols, lchnk   )
    call outfld(zint_fldn,   state%zi,              pcols, lchnk   )
    call outfld(pint_fldn,   state%pint,            pcols, lchnk   )
    call outfld(pmid_fldn,   state%pmid,            pcols, lchnk   )
    call outfld(coszen_fldn, coszrs_pos, pcols, lchnk   )
    call outfld(ts_fldn,    cam_in%ts,    pcols, lchnk   )

#ifdef MAML1
    ![lee1046] in MAML1, output as CRM-level
    call outfld(snowh_fldn,  cam_in%snowhland, pcols, lchnk)
    call outfld(asdir_fldn, cam_in%asdir, pcols, lchnk   )
    call outfld(asdif_fldn, cam_in%asdif, pcols, lchnk   )
    call outfld(aldir_fldn, cam_in%aldir, pcols, lchnk   )
    call outfld(aldif_fldn, cam_in%aldif, pcols, lchnk   )
    call outfld(asdir_pos_fldn, asdir_pos, pcols, lchnk   )
    call outfld(asdif_pos_fldn, asdif_pos, pcols, lchnk   )
    call outfld(aldir_pos_fldn, aldir_pos, pcols, lchnk   )
    call outfld(aldif_pos_fldn, aldif_pos, pcols, lchnk   )
    call outfld(lwup_fldn,  cam_in%lwup,  pcols, lchnk   )
#else
    call outfld(snowh_fldn,     snowh_avg, pcols, lchnk)
    call outfld(asdir_fldn,     asdir_avg, pcols, lchnk   )
    call outfld(asdif_fldn,     asdif_avg, pcols, lchnk   )
    call outfld(aldir_fldn,     aldir_avg, pcols, lchnk   )
    call outfld(aldif_fldn,     aldif_avg, pcols, lchnk   )
    call outfld(asdir_pos_fldn, asdir_pos, pcols, lchnk   )
    call outfld(asdif_pos_fldn, asdif_pos, pcols, lchnk   )
    call outfld(aldir_pos_fldn, aldir_pos, pcols, lchnk   )
    call outfld(aldif_pos_fldn, aldif_pos, pcols, lchnk   )
    call outfld(lwup_fldn,      lwup_avg,  pcols, lchnk   )
#endif

    itim_old = pbuf_old_tim_idx()

    call pbuf_get_field(pbuf, cld_ifld,    ptr, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
    call outfld(cld_fldn,    ptr, pcols, lchnk )

    call pbuf_get_field(pbuf, concld_ifld, ptr, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
    call outfld(concld_fldn, ptr, pcols, lchnk )

    call pbuf_get_field(pbuf, rel_ifld,    ptr )
    call outfld(rel_fldn,    ptr, pcols, lchnk )

    call pbuf_get_field(pbuf, rei_ifld,    ptr )
    call outfld(rei_fldn,    ptr, pcols, lchnk )

    if (mg_microphys) then

       call pbuf_get_field(pbuf,  dei_ifld, ptr     )
       call outfld(dei_fldn,      ptr, pcols, lchnk   )       

       call pbuf_get_field(pbuf,  des_ifld, ptr     )
       call outfld(des_fldn,      ptr, pcols, lchnk   )       

       call pbuf_get_field(pbuf,  mu_ifld, ptr      )
       call outfld(mu_fldn,       ptr, pcols, lchnk   ) 

       call pbuf_get_field(pbuf,  lambdac_ifld, ptr )
       call outfld(lambdac_fldn,  ptr, pcols, lchnk   )       

       call pbuf_get_field(pbuf,  iciwp_ifld, ptr   )
       call outfld(iciwp_fldn,    ptr, pcols, lchnk   )       

       call pbuf_get_field(pbuf,  iclwp_ifld, ptr   )
       call outfld(iclwp_fldn,    ptr, pcols, lchnk   )       

       call pbuf_get_field(pbuf,  icswp_ifld, ptr   )
       call outfld(icswp_fldn,    ptr, pcols, lchnk   )       

       call pbuf_get_field(pbuf,  cldfsnow_ifld, ptr, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
       call outfld(cldfsnow_fldn, ptr, pcols, lchnk   )

    endif

    ! output mixing ratio of rad constituents 

    do i = 1, ngas
       name = 'rad_'//gasnames(i)
       call rad_cnst_get_gas(0, gaslist(i), state, pbuf, mmr)
       call outfld(trim(name), mmr, pcols, lchnk)
    end do

    do i = 1, naer
       name = 'rad_'//aernames(i)
       call rad_cnst_get_aer_mmr(0, i, state, pbuf, mmr)
       call outfld(trim(name), mmr, pcols, lchnk)
    end do

  end subroutine output_rad_data


end module radiation_data
