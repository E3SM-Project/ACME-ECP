module camsrfexch
!-----------------------------------------------------------------------
!
! Module to handle data that is exchanged between the CAM atmosphere
! model and the surface models (land, sea-ice, and ocean).
!
!-----------------------------------------------------------------------
!
! USES:
!
  use shr_kind_mod,  only: r8 => shr_kind_r8, r4 => shr_kind_r4
  use constituents,  only: pcnst
  use ppgrid,        only: pcols, begchunk, endchunk
  use phys_grid,     only: get_ncols_p, phys_grid_initialized
  use infnan,        only: posinf, assignment(=)
  use cam_abortutils,    only: endrun
  use cam_logfile,   only: iulog
  use seq_comm_mct,  only : num_inst_atm
  
  implicit none

!----------------------------------------------------------------------- 
! PRIVATE: Make default data and interfaces private
!----------------------------------------------------------------------- 
  private     ! By default all data is private to this module
!
! Public interfaces
!
  public atm2hub_alloc              ! Atmosphere to surface data allocation method
  public hub2atm_alloc              ! Merged hub surface to atmosphere data allocation method
  public atm2hub_deallocate
  public hub2atm_deallocate
  public cam_export
!
! Public data types
!
  public cam_out_t                  ! Data from atmosphere
  public cam_in_t                   ! Merged surface data

!---------------------------------------------------------------------------
! This is the data that is sent from the atmosphere to the surface models
!---------------------------------------------------------------------------

  type cam_out_t 
     integer  :: lchnk               ! chunk index
     integer  :: ncol                ! number of columns in chunk
      
     real(r8) :: tbot(pcols,num_inst_atm)         ! bottom level temperature
     real(r8) :: thbot(pcols,num_inst_atm)        ! bottom level potential temperature 
     real(r8) :: zbot(pcols,num_inst_atm)         ! bottom level height above surface
     real(r8) :: ubot(pcols,num_inst_atm)         ! bottom level u wind
     real(r8) :: vbot(pcols,num_inst_atm)         ! bottom level v wind
     real(r8) :: qbot(pcols,pcnst,num_inst_atm)   ! bottom level specific humidity
     real(r8) :: pbot(pcols,num_inst_atm)         ! bottom level pressure
     real(r8) :: rho(pcols,num_inst_atm)          ! bottom level density  
     real(r8) :: precsc(pcols,num_inst_atm)       ! convective snow rate (m/sec)
     real(r8) :: precsl(pcols,num_inst_atm)       ! large-scale (stable) snow rate (m/sec)
     real(r8) :: precc(pcols,num_inst_atm)        ! convective precipitationrate (m/sec)
     real(r8) :: precl(pcols,num_inst_atm)        ! large-scale (stable) precipitation rate(m/sec)
     real(r8) :: soll(pcols,num_inst_atm)         ! solar downward NIR beam (W/m2) 
     real(r8) :: sols(pcols,num_inst_atm)         ! solar downward VIS beam (W/m2)
     real(r8) :: solld(pcols,num_inst_atm)        ! solar downward NIR diffuse (W/m2)
     real(r8) :: solsd(pcols,num_inst_atm)        ! solar downward VIS diffuse (W/m2)
     real(r8) :: netsw(pcols,num_inst_atm)        ! Surface solar, net absorbed [W/m2]
     real(r8) :: flwds(pcols,num_inst_atm)        ! longwave down [W/m2]
     
     real(r8) :: co2prog(pcols)      ! prognostic co2
     real(r8) :: co2diag(pcols)      ! diagnostic co2
     real(r8) :: psl(pcols)
     real(r8) :: bcphiwet(pcols)     ! wet deposition of hydrophilic black carbon
     real(r8) :: bcphidry(pcols)     ! dry deposition of hydrophilic black carbon
     real(r8) :: bcphodry(pcols)     ! dry deposition of hydrophobic black carbon
     real(r8) :: ocphiwet(pcols)     ! wet deposition of hydrophilic organic carbon
     real(r8) :: ocphidry(pcols)     ! dry deposition of hydrophilic organic carbon
     real(r8) :: ocphodry(pcols)     ! dry deposition of hydrophobic organic carbon
     real(r8) :: dstwet1(pcols)      ! wet deposition of dust (bin1)
     real(r8) :: dstdry1(pcols)      ! dry deposition of dust (bin1)
     real(r8) :: dstwet2(pcols)      ! wet deposition of dust (bin2)
     real(r8) :: dstdry2(pcols)      ! dry deposition of dust (bin2)
     real(r8) :: dstwet3(pcols)      ! wet deposition of dust (bin3)
     real(r8) :: dstdry3(pcols)      ! dry deposition of dust (bin3)
     real(r8) :: dstwet4(pcols)      ! wet deposition of dust (bin4)
     real(r8) :: dstdry4(pcols)      ! dry deposition of dust (bin4)
  end type cam_out_t 

!---------------------------------------------------------------------------
! This is the merged state of sea-ice, land and ocean surface parameterizations
!---------------------------------------------------------------------------

  type cam_in_t    
     integer  :: lchnk                   ! chunk index
     integer  :: ncol                    ! number of active columns
     
     !modifying the CLM-input vars to reflect the added CRM columns
     real(r8) :: asdir(pcols,num_inst_atm)            ! albedo: shortwave, direct
     real(r8) :: asdif(pcols,num_inst_atm)            ! albedo: shortwave, diffuse
     real(r8) :: aldir(pcols,num_inst_atm)            ! albedo: longwave, direct
     real(r8) :: aldif(pcols,num_inst_atm)            ! albedo: longwave, diffuse
     real(r8) :: lwup(pcols,num_inst_atm)             ! longwave up radiative flux
     real(r8) :: lhf(pcols,num_inst_atm)              ! latent heat flux
     real(r8) :: shf(pcols,num_inst_atm)              ! sensible heat flux
     real(r8) :: wsx(pcols,num_inst_atm)              ! surface u-stress (N)
     real(r8) :: wsy(pcols,num_inst_atm)              ! surface v-stress (N)
     real(r8) :: snowhland(pcols,num_inst_atm)        ! snow depth (liquid water equivalent) over land 
     !modifying the CLM-input vars to reflect the added CRM columns
     
     real(r8) :: tref(pcols)             ! ref height surface air temp
     real(r8) :: qref(pcols)             ! ref height specific humidity 
     real(r8) :: u10(pcols)              ! 10m wind speed
     real(r8) :: ts(pcols)               ! merged surface temp 
     real(r8) :: sst(pcols)              ! sea surface temp
     real(r8) :: snowhice(pcols)         ! snow depth over ice
     real(r8) :: fco2_lnd(pcols)         ! co2 flux from lnd
     real(r8) :: fco2_ocn(pcols)         ! co2 flux from ocn
     real(r8) :: fdms(pcols)             ! dms flux
     real(r8) :: landfrac(pcols)         ! land area fraction
     real(r8) :: icefrac(pcols)          ! sea-ice areal fraction
     real(r8) :: ocnfrac(pcols)          ! ocean areal fraction
     real(r8), pointer, dimension(:) :: ram1  !aerodynamical resistance (s/m) (pcols)
     real(r8), pointer, dimension(:) :: fv    !friction velocity (m/s) (pcols)
     real(r8), pointer, dimension(:) :: soilw !volumetric soil water (m3/m3)
     real(r8) :: cflx(pcols,pcnst)       ! constituent flux (emissions)
     real(r8) :: ustar(pcols)            ! atm/ocn saved version of ustar
     real(r8) :: re(pcols)               ! atm/ocn saved version of re
     real(r8) :: ssq(pcols)              ! atm/ocn saved version of ssq
     real(r8), pointer, dimension(:,:) :: depvel ! deposition velocities
     real(r8), pointer, dimension(:,:) :: dstflx ! dust fluxes
     real(r8), pointer, dimension(:,:) :: meganflx ! MEGAN fluxes

  end type cam_in_t    

!===============================================================================
CONTAINS
!===============================================================================

!----------------------------------------------------------------------- 
! 
! BOP
!
! !IROUTINE: hub2atm_alloc
!
! !DESCRIPTION:
!
!   Allocate space for the surface to atmosphere data type. And initialize
!   the values.
! 
!-----------------------------------------------------------------------
!
! !INTERFACE
!
  subroutine hub2atm_alloc( cam_in )
    use seq_drydep_mod,  only: lnd_drydep, n_drydep
    use cam_cpl_indices, only: index_x2a_Sl_ram1, index_x2a_Sl_fv, index_x2a_Sl_soilw, index_x2a_Fall_flxdst1
    use cam_cpl_indices, only: index_x2a_Fall_flxvoc
    use shr_megan_mod,   only: shr_megan_mechcomps_n

!
!!ARGUMENTS:
!
   type(cam_in_t), pointer ::  cam_in(:)     ! Merged surface state
!
!!LOCAL VARIABLES:
!
    integer :: c        ! chunk index
    integer :: ierror   ! Error code
!----------------------------------------------------------------------- 
! 
! EOP
!
    if ( .not. phys_grid_initialized() ) call endrun( "HUB2ATM_ALLOC error: phys_grid not called yet" )
    allocate (cam_in(begchunk:endchunk), stat=ierror)
    if ( ierror /= 0 )then
      write(iulog,*) 'Allocation error: ', ierror
      call endrun('HUB2ATM_ALLOC error: allocation error')
    end if

    do c = begchunk,endchunk
       nullify(cam_in(c)%ram1)
       nullify(cam_in(c)%fv)
       nullify(cam_in(c)%soilw)
       nullify(cam_in(c)%depvel)
       nullify(cam_in(c)%dstflx)
       nullify(cam_in(c)%meganflx)
    enddo  
    do c = begchunk,endchunk 
       if (index_x2a_Sl_ram1>0) then
          allocate (cam_in(c)%ram1(pcols), stat=ierror)
          if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error ram1')
       endif
       if (index_x2a_Sl_fv>0) then
          allocate (cam_in(c)%fv(pcols), stat=ierror)
          if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error fv')
       endif
       if (index_x2a_Sl_soilw /= 0) then
          allocate (cam_in(c)%soilw(pcols), stat=ierror)
          if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error soilw')
       end if
       if (index_x2a_Fall_flxdst1>0) then
          ! Assume 4 bins from surface model ....
          allocate (cam_in(c)%dstflx(pcols,4), stat=ierror)
          if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error dstflx')
       endif
       if ( index_x2a_Fall_flxvoc>0 .and. shr_megan_mechcomps_n>0 ) then
          allocate (cam_in(c)%meganflx(pcols,shr_megan_mechcomps_n), stat=ierror)
          if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error meganflx')
       endif
    end do

    if (lnd_drydep .and. n_drydep>0) then
       do c = begchunk,endchunk 
          allocate (cam_in(c)%depvel(pcols,n_drydep), stat=ierror)
          if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error depvel')
       end do
    endif

    do c = begchunk,endchunk
       cam_in(c)%lchnk = c
       cam_in(c)%ncol  = get_ncols_p(c)
       cam_in(c)%asdir    (:,:) = 0._r8
       cam_in(c)%asdif    (:,:) = 0._r8
       cam_in(c)%aldir    (:,:) = 0._r8
       cam_in(c)%aldif    (:,:) = 0._r8
       cam_in(c)%lwup     (:,:) = 0._r8
       cam_in(c)%lhf      (:,:) = 0._r8
       cam_in(c)%shf      (:,:) = 0._r8
       cam_in(c)%wsx      (:,:) = 0._r8
       cam_in(c)%wsy      (:,:) = 0._r8
       cam_in(c)%snowhland(:,:) = 0._r8
       cam_in(c)%tref     (:) = 0._r8
       cam_in(c)%qref     (:) = 0._r8
       cam_in(c)%u10      (:) = 0._r8
       cam_in(c)%ts       (:) = 0._r8
       cam_in(c)%sst      (:) = 0._r8
       cam_in(c)%snowhice (:) = 0._r8
       cam_in(c)%fco2_lnd (:) = 0._r8
       cam_in(c)%fco2_ocn (:) = 0._r8
       cam_in(c)%fdms     (:) = 0._r8
       cam_in(c)%landfrac (:) = posinf
       cam_in(c)%icefrac  (:) = posinf
       cam_in(c)%ocnfrac  (:) = posinf

       if (associated(cam_in(c)%ram1)) &
            cam_in(c)%ram1  (:) = 0.1_r8
       if (associated(cam_in(c)%fv)) &
            cam_in(c)%fv    (:) = 0.1_r8
       if (associated(cam_in(c)%soilw)) &
            cam_in(c)%soilw (:) = 0.0_r8
       if (associated(cam_in(c)%dstflx)) &
            cam_in(c)%dstflx(:,:) = 0.0_r8
       if (associated(cam_in(c)%meganflx)) &
            cam_in(c)%meganflx(:,:) = 0.0_r8

       cam_in(c)%cflx   (:,:) = 0._r8
       cam_in(c)%ustar    (:) = 0._r8
       cam_in(c)%re       (:) = 0._r8
       cam_in(c)%ssq      (:) = 0._r8
       if (lnd_drydep .and. n_drydep>0) then
          cam_in(c)%depvel (:,:) = 0._r8
       endif

    end do

  end subroutine hub2atm_alloc

!
!===============================================================================
!

!----------------------------------------------------------------------- 
! 
! BOP
!
! !IROUTINE: atm2hub_alloc
!
! !DESCRIPTION:
!
!   Allocate space for the atmosphere to surface data type. And initialize
!   the values.
! 
!-----------------------------------------------------------------------
!
! !INTERFACE
!
  subroutine atm2hub_alloc( cam_out )
!
!!USES:
!
!
!!ARGUMENTS:
!
   type(cam_out_t), pointer :: cam_out(:)    ! Atmosphere to surface input
!
!!LOCAL VARIABLES:
!
    integer :: c            ! chunk index
    integer :: ierror       ! Error code
    !----------------------------------------------------------------------- 

    if ( .not. phys_grid_initialized() ) call endrun( "ATM2HUB_ALLOC error: phys_grid not called yet" )
    allocate (cam_out(begchunk:endchunk), stat=ierror)
    if ( ierror /= 0 )then
      write(iulog,*) 'Allocation error: ', ierror
      call endrun('ATM2HUB_ALLOC error: allocation error')
    end if

    do c = begchunk,endchunk
       cam_out(c)%lchnk       = c
       cam_out(c)%ncol        = get_ncols_p(c)
       cam_out(c)%thbot(:,:)    = 0._r8
       cam_out(c)%tbot(:,:)     = 0._r8
       cam_out(c)%zbot(:,:)     = 0._r8
       cam_out(c)%ubot(:,:)     = 0._r8
       cam_out(c)%vbot(:,:)     = 0._r8
       cam_out(c)%qbot(:,:,:)   = 0._r8
       cam_out(c)%pbot(:,:)     = 0._r8
       cam_out(c)%rho(:,:)      = 0._r8
       cam_out(c)%precsc(:,:)   = 0._r8
       cam_out(c)%precsl(:,:)   = 0._r8
       cam_out(c)%precc(:,:)    = 0._r8
       cam_out(c)%precl(:,:)    = 0._r8
       cam_out(c)%soll(:,:)     = 0._r8
       cam_out(c)%sols(:,:)     = 0._r8
       cam_out(c)%solld(:,:)    = 0._r8
       cam_out(c)%solsd(:,:)    = 0._r8
       cam_out(c)%flwds(:,:)    = 0._r8
       cam_out(c)%netsw(:,:)    = 0._r8
       cam_out(c)%co2prog(:)  = 0._r8
       cam_out(c)%co2diag(:)  = 0._r8
       cam_out(c)%psl(:)      = 0._r8
       cam_out(c)%bcphidry(:) = 0._r8
       cam_out(c)%bcphodry(:) = 0._r8
       cam_out(c)%bcphiwet(:) = 0._r8
       cam_out(c)%ocphidry(:) = 0._r8
       cam_out(c)%ocphodry(:) = 0._r8
       cam_out(c)%ocphiwet(:) = 0._r8
       cam_out(c)%dstdry1(:)  = 0._r8
       cam_out(c)%dstwet1(:)  = 0._r8
       cam_out(c)%dstdry2(:)  = 0._r8
       cam_out(c)%dstwet2(:)  = 0._r8
       cam_out(c)%dstdry3(:)  = 0._r8
       cam_out(c)%dstwet3(:)  = 0._r8
       cam_out(c)%dstdry4(:)  = 0._r8
       cam_out(c)%dstwet4(:)  = 0._r8
    end do

  end subroutine atm2hub_alloc

  subroutine atm2hub_deallocate(cam_out)
    type(cam_out_t), pointer :: cam_out(:)    ! Atmosphere to surface input
    if(associated(cam_out)) then
       deallocate(cam_out)
    end if
    nullify(cam_out)

  end subroutine atm2hub_deallocate
  subroutine hub2atm_deallocate(cam_in)
    type(cam_in_t), pointer :: cam_in(:)    ! Atmosphere to surface input
    integer :: c

    if(associated(cam_in)) then
       do c=begchunk,endchunk
          if(associated(cam_in(c)%ram1)) then
             deallocate(cam_in(c)%ram1)
             nullify(cam_in(c)%ram1)
          end if
          if(associated(cam_in(c)%fv)) then
             deallocate(cam_in(c)%fv)
             nullify(cam_in(c)%fv)
          end if
          if(associated(cam_in(c)%soilw)) then
             deallocate(cam_in(c)%soilw)
             nullify(cam_in(c)%soilw)
          end if
          if(associated(cam_in(c)%dstflx)) then
             deallocate(cam_in(c)%dstflx)
             nullify(cam_in(c)%dstflx)
          end if
          if(associated(cam_in(c)%meganflx)) then
             deallocate(cam_in(c)%meganflx)
             nullify(cam_in(c)%meganflx)
          end if
          if(associated(cam_in(c)%depvel)) then
             deallocate(cam_in(c)%depvel)
             nullify(cam_in(c)%depvel)
          end if
          
       enddo

       deallocate(cam_in)
    end if
    nullify(cam_in)

  end subroutine hub2atm_deallocate


!======================================================================

subroutine cam_export(state,cam_out,pbuf)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Transfer atmospheric fields into necessary surface data structures
! 
! Author: L. Bath  CMS Contact: M. Vertenstein
! 
!-----------------------------------------------------------------------
   use physics_types,    only: physics_state
   use ppgrid,           only: pver
   use cam_history,      only: outfld
   use comsrf,           only: psm1, srfrpdel, prcsnw
   use chem_surfvals,    only: chem_surfvals_get
   use co2_cycle,        only: co2_transport, c_i
   use physconst,        only: mwdry, mwco2
   use constituents,     only: pcnst
   use cam_control_mod,  only: rair
   use physics_buffer,   only: pbuf_get_index, pbuf_get_field, physics_buffer_desc
   use phys_control,     only: phys_getopts
   implicit none

   !------------------------------Arguments--------------------------------
   !
   ! Input arguments
   !
   type(physics_state),  intent(in)    :: state
   type (cam_out_t),     intent(inout) :: cam_out
   type(physics_buffer_desc), pointer  :: pbuf(:)

   !
   !---------------------------Local variables-----------------------------
   !
   integer :: i              ! Longitude index
   integer :: m              ! constituent index
   integer :: lchnk          ! Chunk index
   integer :: ncol
   integer :: prec_dp_idx, snow_dp_idx, prec_sh_idx, snow_sh_idx
   integer :: prec_sed_idx,snow_sed_idx,prec_pcw_idx,snow_pcw_idx
   integer :: vmag_gust_idx
   real(r8) :: umb(pcols), vmb(pcols),vmag(pcols)
   logical :: use_SPCAM, use_MAML ! flag for MMF and MAML, respectively
! [lee1046] TODO: decide whether to keep #ifdef MAML or not
#ifdef MAML
   !CRM-level variables 
   integer :: j
   integer :: crm_t_idx
   integer :: crm_qv_idx
   integer :: crm_u_idx
   integer :: crm_v_idx
   integer :: crm_pcp_idx
   integer :: crm_snw_idx
   integer :: crm_angle_idx
   real(r8), pointer :: crm_t(:,:,:,:)
   real(r8), pointer :: crm_qv(:,:,:,:)
   real(r8), pointer :: crm_u(:,:,:,:)
   real(r8), pointer :: crm_v(:,:,:,:)
   real(r8), pointer :: crm_pcp(:,:,:)
   real(r8), pointer :: crm_snw(:,:,:)
   real(r8), pointer :: crm_angle(:)
#endif


   real(r8), pointer :: prec_dp(:)                 ! total precipitation   from ZM convection
   real(r8), pointer :: snow_dp(:)                 ! snow from ZM   convection
   real(r8), pointer :: prec_sh(:)                 ! total precipitation   from Hack convection
   real(r8), pointer :: snow_sh(:)                 ! snow from   Hack   convection
   real(r8), pointer :: prec_sed(:)                ! total precipitation   from ZM convection
   real(r8), pointer :: snow_sed(:)                ! snow from ZM   convection
   real(r8), pointer :: prec_pcw(:)                ! total precipitation   from Hack convection
   real(r8), pointer :: snow_pcw(:)                ! snow from Hack   convection
   real(r8), pointer :: vmag_gust(:)

   !-----------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol

   prec_dp_idx = pbuf_get_index('PREC_DP')
   snow_dp_idx = pbuf_get_index('SNOW_DP')
   prec_sh_idx = pbuf_get_index('PREC_SH')
   snow_sh_idx = pbuf_get_index('SNOW_SH')
   prec_sed_idx = pbuf_get_index('PREC_SED')
   snow_sed_idx = pbuf_get_index('SNOW_SED')
   prec_pcw_idx = pbuf_get_index('PREC_PCW')
   snow_pcw_idx = pbuf_get_index('SNOW_PCW')
   vmag_gust_idx = pbuf_get_index('vmag_gust')
   
   call pbuf_get_field(pbuf, prec_dp_idx, prec_dp)
   call pbuf_get_field(pbuf, snow_dp_idx, snow_dp)
   call pbuf_get_field(pbuf, prec_sh_idx, prec_sh)
   call pbuf_get_field(pbuf, snow_sh_idx, snow_sh)
   call pbuf_get_field(pbuf, prec_sed_idx, prec_sed)
   call pbuf_get_field(pbuf, snow_sed_idx, snow_sed)
   call pbuf_get_field(pbuf, prec_pcw_idx, prec_pcw)
   call pbuf_get_field(pbuf, snow_pcw_idx, snow_pcw)
   call pbuf_get_field(pbuf, vmag_gust_idx, vmag_gust)
   
   do j=1,num_inst_atm
      cam_out%pbot(1:ncol,j)         = state%pmid(1:ncol,pver)
      cam_out%zbot(1:ncol,j)         = state%zm  (1:ncol,pver)
      cam_out%qbot(1:ncol,2:pcnst,j) = state%q   (1:ncol,pver,2:pcnst)
   end do
   psm1           (1:ncol,lchnk)  = state%ps(1:ncol)
   srfrpdel       (1:ncol,lchnk)  = state%rpdel(1:ncol,pver)
   cam_out%co2diag(1:ncol)        = chem_surfvals_get('CO2VMR') * 1.0e+6_r8 
   if (co2_transport()) then
      cam_out%co2prog(1:ncol) = state%q(1:ncol,pver,c_i(4)) * 1.0e+6_r8 *mwdry/mwco2
   end if

   !
   ! temperature, qv, u, v and precipitation fields in the surface interface
   ! are assigned from CRM (CAM) in MMF with MAML (non-MMF or MMF without MAML)
   !
   
   call phys_getopts(use_SPCAM_out = use_SPCAM)
   call phys_getopts(use_MAML_out  = use_MAML)

   if(use_MMF .and. use_MAML) then  ! [lee1046] Or if(num_inst_atm.gt.1)
      crm_t_idx     = pbuf_get_index('CRM_T')
      crm_qv_idx    = pbuf_get_index('CRM_QV_RAD')
      crm_u_idx     = pbuf_get_index('CRM_U')
      crm_v_idx     = pbuf_get_index('CRM_V')
      crm_pcp_idx   = pbuf_get_index('CRM_PCP')
      crm_snw_idx   = pbuf_get_index('CRM_SNW')
      crm_angle_idx = pbuf_get_index('CRM_ANGLE')

      call pbuf_get_field(pbuf, crm_t_idx    , crm_t)
      call pbuf_get_field(pbuf, crm_qv_idx   , crm_qv)
      call pbuf_get_field(pbuf, crm_u_idx    , crm_u)
      call pbuf_get_field(pbuf, crm_v_idx    , crm_v)
      call pbuf_get_field(pbuf, crm_pcp_idx  , crm_pcp)
      call pbuf_get_field(pbuf, crm_snw_idx  , crm_snw)
      call pbuf_get_field(pbuf, crm_angle_idx, crm_angle)
      
      cam_out%precl (1:ncol,1:num_inst_atm) = 0._r8     ! large-scale precip set to zero
      cam_out%precsl(1:ncol,1:num_inst_atm) = 0._r8     ! large-scale precip set to zero

      do j=1,num_inst_atm
         cam_out%tbot  (1:ncol,j)   = crm_t(1:ncol,j,1,1)
         cam_out%qbot  (1:ncol,1,j) = crm_qv(1:ncol,j,1,1)
         ! u and v will use CRM value (must transform because of CRM orientation)
         cam_out%ubot  (1:ncol,j)   = crm_u(1:ncol,j,1,1) * cos(crm_angle(1:ncol)) - crm_v(1:ncol,j,1,1) * sin(crm_angle(1:ncol))
         cam_out%vbot  (1:ncol,j)   = crm_v(1:ncol,j,1,1) * cos(crm_angle(1:ncol)) + crm_u(1:ncol,j,1,1) * sin(crm_angle(1:ncol))
         cam_out%precc (1:ncol,j)   = crm_pcp(1:ncol,j,1)  
         cam_out%precsc(1:ncol,j)   = crm_snw(1:ncol,j,1)
      end do 
   else ! if(use_MMF .and. use_MAML)
   !
   ! when use_MAML = .false. => num_inst_atm = 1
   !
      cam_out%tbot(1:ncol,1)   = state%t(1:ncol,pver)
      cam_out%qbot(1:ncol,1,1) = state%q(1:ncol,pver,1) 
      umb         (1:ncol)     = state%u(1:ncol,pver)
      vmb         (1:ncol)     = state%v(1:ncol,pver)
      vmag        (1:ncol)     = max(1.e-5_r8,sqrt( umb(1:ncol)**2._r8 + vmb(1:ncol)**2._r8))            
      !PMA adds gustiness to surface scheme c20181128
      cam_out%ubot(1:ncol,1)   = state%u(1:ncol,pver) * ((vmag_gust(1:ncol)+vmag(1:ncol))/vmag(1:ncol))
      cam_out%vbot(1:ncol,1)   = state%v(1:ncol,pver) * ((vmag_gust(1:ncol)+vmag(1:ncol))/vmag(1:ncol))
      !
      ! Precipation and snow rates from shallow convection, deep convection and stratiform processes.
      ! Compute total convective and stratiform precipitation and snow rates
      !
      cam_out%precc (1:ncol,1) = prec_dp (1:ncol) + prec_sh (1:ncol)
      cam_out%precl (1:ncol,1) = prec_sed(1:ncol) + prec_pcw(1:ncol)
      cam_out%precsc(1:ncol,1) = snow_dp (1:ncol) + snow_sh (1:ncol)
      cam_out%precsl(1:ncol,1) = snow_sed(1:ncol) + snow_pcw(1:ncol)

      ! SPCAM says below checks are unnecessary because it is already done in
      ! the parameterization. Can we delete it?
      do i=1,ncol
         if (cam_out%precc(i) .lt.0._r8) cam_out%precc(i)  = 0._r8
         if (cam_out%precl(i) .lt.0._r8) cam_out%precl(i)  = 0._r8
         if (cam_out%precsc(i).lt.0._r8) cam_out%precsc(i) = 0._r8
         if (cam_out%precsl(i).lt.0._r8) cam_out%precsl(i) = 0._r8
         if (cam_out%precsc(i).gt.cam_out%precc(i)) cam_out%precsc(i) = cam_out%precc(i)
         if (cam_out%precsl(i).gt.cam_out%precl(i)) cam_out%precsl(i) = cam_out%precl(i)
      end do
   end if
   !
   ! total snowfall rate: needed by slab ocean model
   ! for MAML, ocean model is not multi-instance, therefore we average total
   ! snowfall rate across the instances. This averaging doesn't affect for
   ! non-MAML case as num_inst_atm = 1
   !
   prcsnw(1:ncol,lchnk) = 0._r8
   do j = 1, num_inst_atm
      prcsnw       (1:ncol,lchnk) = prcsnw(1:ncol,lchnk) + cam_out%precsc(1:ncol,j) + cam_out%precsl(1:ncol,j)   
      cam_out%thbot(1:ncol,j)     = cam_out%tbot(1:ncol,j,1,1) * state%exner(1:ncol,pver) ! potential temperature
      cam_out%rho  (1:ncol,j)     = cam_out%pbot(1:ncol,j)/(rair*cam_out%tbot(1:ncol,j))  ! air density
   end if
   prcsnw(1:ncol,lchnk) = prcsnw(1:ncol,lchnk)/num_inst_atm

end subroutine cam_export

end module camsrfexch
