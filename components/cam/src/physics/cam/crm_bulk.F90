module crm_bulk_mod

! the chec here just allows the module to be hidden during development
#if defined( SP_CRM_BULK )
#if defined( CRM )

   !-------------------------------------------------------------------------------------
   ! Purpose: 
   !  Provides a method for convective transport of tracers using 
   !  cloud statistics from the explicit convection of the CRM.
   !  This is meant to be an interim alternative to the ECPP 
   !  (explicit-cloud parameterized-pollutant) approach. 
   !
   !  Revision history: 
   !     2018: Walter Hannah - initial version based on Minghuai Wang's 
   !                           "crmclouds_camaerosols" module
   !-------------------------------------------------------------------------------------
   use shr_kind_mod,       only: r8 => shr_kind_r8
   use cam_abortutils,     only: endrun
   use ppgrid
   use params

   implicit none

   ! public :: crm_bulk_aerosol_activation
   public :: crm_bulk_aerosol_wet_removal
   public :: crm_bulk_transport
   public :: crm_bulk_aero_mix_nuc 
   ! public :: crm_bulk_diagnose_cloudy_clear

   private :: crm_bulk_transport_tend 

contains 
!==================================================================================================
!==================================================================================================

!!! the ifdef is used to hide crm_bulk_aerosol_activation() during development

! #if defined( SP_BULK_AERO_ACT )

! subroutine crm_bulk_aerosol_activation(state, pbuf, ptend)
!    !-------------------------------------------------------------------------------------
!    ! Purpose: calculate activation of aerosols using CRM cloud information
!    ! Author: Walter Hannah (LLNL), 2018
!    !-------------------------------------------------------------------------------------
!    use physics_types,   only: physics_state, physics_ptend, physics_ptend_init
!    use physics_buffer,  only: physics_buffer_desc, pbuf_old_tim_idx, &
!                               pbuf_get_index, pbuf_get_field
!    use constituents,    only: pcnst, cnst_get_ind
!    use modal_aero_data,  only: ntot_amode

!    !!! Input Arguments
!    type(physics_state), intent(in )   :: state              ! Physics state variables
!    type(physics_buffer_desc), pointer :: pbuf(:)            ! physics buffer

!    !!! Output Arguments
!    type(physics_ptend), intent(out)   :: ptend              ! indivdual parameterization tendencies

!    !!! Local variables
!    integer  :: i, k, m, lchnk, ncol
!    integer  :: ixcldice, ixcldliq                           ! constituent indices for cloud liquid and ice water.

!    real(r8)                        :: wbar             ! mean updraft velocity (cm/s)
!    real(r8)                        :: sigw             ! standard deviation of updraft velocity (cm/s)
!    real(r8)                        :: wdiab            ! diabatic vertical velocity (cm/s)
!    real(r8)                        :: wminf            ! lower limit for integration over updraft spectrum (cm/s)
!    real(r8)                        :: wmaxf            ! upper limit for integration over updraft spectrum (cm/s)
!    real(r8), dimension(ntot_amode) :: naerosol         ! interstitial+activated number conc (#/m3)
!    real(r8), dimension(ntot_amode) :: vaerosol         ! interstitial + activated aerosol volume mixing ratio (m3/m3)
!    real(r8), dimension(ntot_amode) :: hygro            ! current hygroscopicity for int+act

!    real(r8), dimension(ntot_amode) :: fn               ! number fraction of aerosols activated
!    real(r8), dimension(ntot_amode) :: fm               ! mass fraction of aerosols activated
!    real(r8), dimension(ntot_amode) :: fluxn            ! see activate_modal() in ndrop.F90
!    real(r8), dimension(ntot_amode) :: fluxm            ! see activate_modal() in ndrop.F90
!    real(r8)                        :: flux_fullact     ! see activate_modal() in ndrop.F90
   
   
!    !-------------------------------------------------------------------------------------
!    !-------------------------------------------------------------------------------------
   
!    lchnk = state%lchnk
!    ncol  = state%ncol

!    !!! Initialize ptend
!    call cnst_get_ind('CLDLIQ', ixcldliq)
!    call cnst_get_ind('CLDICE', ixcldice)
!    lq(:)        = .true.
!    lq(1)        = .false.  ! vapor
!    lq(ixcldliq) = .false.  ! liquid
!    lq(ixcldice) = .false.  ! ice
!    call physics_ptend_init(ptend,state%psetcols,'crm_bulk_aerosol_activation',lq=lq)

!    !-------------------------------------------------------------------------------------
!    !-------------------------------------------------------------------------------------

!    !!! subroutine activate_modal Interface:
!    !!!    real(r8), intent(in) :: wbar          ! grid cell mean vertical velocity (m/s)
!    !!!    real(r8), intent(in) :: sigw          ! subgrid standard deviation of vertical vel (m/s)
!    !!!    real(r8), intent(in) :: wdiab         ! diabatic vertical velocity (0 if adiabatic)
!    !!!    real(r8), intent(in) :: wminf         ! minimum updraft velocity for integration (m/s)
!    !!!    real(r8), intent(in) :: wmaxf         ! maximum updraft velocity for integration (m/s)
!    !!!    real(r8), intent(in) :: tair          ! air temperature (K)
!    !!!    real(r8), intent(in) :: rhoair        ! air density (kg/m3)
!    !!!    real(r8), intent(in) :: na(:)         ! aerosol number concentration (/m3)
!    !!!    integer,  intent(in) :: nmode         ! number of aerosol modes
!    !!!    real(r8), intent(in) :: volume(:)     ! aerosol volume concentration (m3/m3)
!    !!!    real(r8), intent(in) :: hygro(:)      ! hygroscopicity of aerosol mode
!    !!!    real(r8), intent(out) :: fn(:)        ! number fraction of aerosols activated
!    !!!    real(r8), intent(out) :: fm(:)        ! mass fraction of aerosols activated
!    !!!    real(r8), intent(out) :: fluxn(:)     ! flux of activated aerosol number fraction into cloud (cm/s)
!    !!!    real(r8), intent(out) :: fluxm(:)     ! flux of activated aerosol mass fraction into cloud (cm/s)
!    !!!    real(r8), intent(out) :: flux_fullact ! flux of activated aerosol fraction assuming 100% activation (cm/s)

!    tair = state%t(:,:)
!    use physconst, only: gravit, rair, rhoh2o
!    rhoair_i(k) = pmid(icol,k)/(rair*t(icol,k))
!    rhoair = 


!    call activate_modal(wbar, sigw, wdiab, wminf, wmaxf, tair, rhoair,   &
!                        naerosol, ntot_amode, vaerosol, hygro,           &
!                        fn, fm, fluxn, fluxm, flux_fullact               )

!    !-------------------------------------------------------------------------------------
!    !-------------------------------------------------------------------------------------

! end subroutine crm_bulk_aerosol_activation

! #endif /* SP_BULK_AERO_ACT */

!==================================================================================================
!==================================================================================================
subroutine crm_bulk_aerosol_wet_removal(state, pbuf, ptend)
   !-------------------------------------------------------------------------------------
   ! Purpose: calculate wet removal (i.e. scavenging) of aerosols
   !          Scheme is based on supplementary material from:
   !              Wang et al., 2013: Sensitivity of remote aerosol distributions to 
   !                 representation of cloud–aerosol interactions in a global climate model, 
   !                 Geosci. Model Dev., 6, 765-782, https://doi.org/10.5194/gmd-6-765-2013.
   ! Author: Walter Hannah (LLNL), 2018
   !-------------------------------------------------------------------------------------
   use physics_types,   only: physics_state, physics_ptend, physics_ptend_init
   use physics_buffer,  only: physics_buffer_desc, pbuf_old_tim_idx, &
                              pbuf_get_index, pbuf_get_field
   use constituents,    only: pcnst, cnst_get_ind

   !!! Input Arguments
   type(physics_state), intent(in )   :: state              ! Physics state variables
   type(physics_buffer_desc), pointer :: pbuf(:)            ! physics buffer

   !!! Output Arguments
   type(physics_ptend), intent(out)   :: ptend              ! indivdual parameterization tendencies

   !!! Local variables
   integer  :: i, k, m, lchnk, ncol
   integer  :: ixcldice, ixcldliq                           ! constituent indices for cloud liquid and ice water.
   real(r8) :: available_water                              ! available cloud water to be rained out
   logical,  dimension(pcnst)            :: lq              ! flag for ptend
   real(r8), dimension(pcols,pver)       :: aero_loss_rate  ! 
   real(r8), dimension(pcols,pver,pcnst) :: cld_borne_aero  ! cloud borne aerosol mass mixing ratio 

   !!! physics buffer fields 
   integer itim, ifld
   real(r8), pointer, dimension(:,:) :: cld_fraction        ! cloud fraction (current time step)
   real(r8), pointer, dimension(:,:) :: rain_production     ! rain production rate 
   
   !-------------------------------------------------------------------------------------
   !-------------------------------------------------------------------------------------
   
   lchnk = state%lchnk
   ncol  = state%ncol

   !!! Initialize ptend
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)
   lq(:)        = .true.
   lq(1)        = .false.  ! vapor
   lq(ixcldliq) = .false.  ! liquid
   lq(ixcldice) = .false.  ! ice
   call physics_ptend_init(ptend,state%psetcols,'crm_bulk_aerosol_wet_removal',lq=lq)

   !!! retreive pbuf fields
   ! itim = pbuf_old_tim_idx()
   ! ifld = pbuf_get_index('CLD')
   ! call pbuf_get_field(pbuf,ifld, cld_fraction    ,start=(/1,1,itim/), kount=(/pcols,pver,1/) )
   ! ifld = pbuf_get_index('PRAIN_CRM' )
   ! call pbuf_get_field(pbuf,ifld, rain_production ,start=(/1,1/)     , kount=(/pcols,pver/) )

   !!! retreive pbuf fields
   call pbuf_get_field(pbuf,pbuf_get_index('CLD'),        cld_fraction     )
   call pbuf_get_field(pbuf,pbuf_get_index('PRAIN_CRM' ), rain_production  )

   !!! TEMPORARY - this needs to change to reflect activated aerosols - use crm_bulk_aerosol_activation()
   cld_borne_aero = 0.5_r8 * state%q(:,:,:)

   do i = 1,ncol
      do k = 1,pver
         
         !!! check if there sufficient cloud water available
         available_water = cld_fraction(i,k) * state%q(i,k,ixcldliq)

         if ( available_water > 0.0_r8 ) then

            !!! Calculate updraft wet removal loss rate
            aero_loss_rate(i,k) = rain_production(i,k) / available_water

            !!! Calculate updraft wet removal tendency
            do m = 1, pcnst
               ptend%q(i,k,m) = -1.0_r8 * aero_loss_rate(i,k) * cld_borne_aero(i,k,m)
            end do ! m=1,pcnst  

         end if ! available_water

      end do ! k=1,pver
   end do ! i=1,ncol

   !-------------------------------------------------------------------------------------
   !-------------------------------------------------------------------------------------

end subroutine crm_bulk_aerosol_wet_removal

!==================================================================================================
!==================================================================================================
subroutine crm_bulk_transport(state, pbuf, ptend)
   !-------------------------------------------------------------------------------------
   ! Purpose:  to do convective transport of tracer species using the 
   !           cloud updraft and downdraft mass fluxes from the CRM.
   !           Adopted from crmclouds_convect_tend() by Minghuai Wang
   !           and convtran() from the Zhang-McFarlane scheme.
   !           Original version used the "gathered" array terminology 
   !           but for the CRM we omit this and operate on all columns.
   ! Walter Hannah (LLNL), 2018
   !-------------------------------------------------------------------------------------
   use physics_types,   only: physics_state, physics_ptend, physics_ptend_init
   use physics_buffer,  only: physics_buffer_desc, pbuf_old_tim_idx, pbuf_get_index, pbuf_get_field
   use constituents,    only: pcnst, cnst_get_ind, cnst_get_type_byind
   use error_messages,  only: alloc_err  

   !!! Input Arguments
   type(physics_state), intent(in )   :: state          ! Physics state variables
   type(physics_buffer_desc), pointer :: pbuf(:)      ! physics buffer

   !!! Output Arguments
   type(physics_ptend), intent(out)   :: ptend        ! indivdual parameterization tendencies

   !!! Local variables
   integer :: i, m, lchnk, ncol
   integer :: ixcldice, ixcldliq                      ! constituent indices for cloud liquid and ice water.
   real(r8), dimension(pcols,pver)  :: dpdry          ! pressure thickness for dry air
   real(r8), dimension(pcols,pver)  :: dp             ! pressure thickness in mbs (between upper/lower interface).
   logical,  dimension(pcols,pcnst) :: dry_const_flag ! flag for dry constituents

   !!! physics buffer fields 
   ! integer itim, ifld
   real(r8), pointer, dimension(:,:,:) :: frac_insol  ! fraction of transported species that are insoluble
   real(r8), pointer, dimension(:,:)   :: mu          ! updraft mass flux       (pcols,pver,begchunk:endchunk)
   real(r8), pointer, dimension(:,:)   :: eu          ! updraft entrainment     (pcols,pver,begchunk:endchunk)
   real(r8), pointer, dimension(:,:)   :: du          ! updraft detrainment     (pcols,pver,begchunk:endchunk)
   real(r8), pointer, dimension(:,:)   :: md          ! downdraft mass flux     (pcols,pver,begchunk:endchunk)
   real(r8), pointer, dimension(:,:)   :: ed          ! downdraft entrainment   (pcols,pver,begchunk:endchunk)

   ! real(r8), pointer, dimension(:) :: cld_top_idx_flt ! index of cloud top      (pcols,begchunk:endchunk)
   ! real(r8), pointer, dimension(:) :: cld_bot_idx_flt ! index of cloud bottom   (pcols,begchunk:endchunk)

   real(r8), pointer, dimension(:) :: cld_top_idx_ptr ! index of cloud top      (pcols,begchunk:endchunk)
   real(r8), pointer, dimension(:) :: cld_bot_idx_ptr ! index of cloud bottom   (pcols,begchunk:endchunk)

   integer, dimension(pcols) :: cld_top_idx           ! index of cloud top
   integer, dimension(pcols) :: cld_bot_idx           ! index of cloud bottom
   logical, dimension(pcnst) :: lq                    ! flag to indicate whether to do transport calculation
   !-------------------------------------------------------------------------------------
   !-------------------------------------------------------------------------------------
   lchnk = state%lchnk
   ncol  = state%ncol

   !!! Initialize ptend
   ! transport of all trace species except 
   ! cloud liquid and cloud ice is done here 
   ! because we need to do the scavenging first 
   ! to determine the interstitial fraction.
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)
   lq(:)        = .true.
   lq(1)        = .false.  ! vapor
   lq(ixcldliq) = .false.  ! liquid
   lq(ixcldice) = .false.  ! ice
   
   call physics_ptend_init(ptend,state%psetcols,'crm_bulk_transport',lq=lq)

   !!! set the flag for dry constituents
   dry_const_flag(:,:) = .false.
   do m = 1,pcnst
      if (cnst_get_type_byind(m).eq.'dry') then
         dry_const_flag(:,m) = .true.
      end if
   end do

   !!! Associate pointers with physics buffer fields
   call pbuf_get_field(pbuf, pbuf_get_index('FRACIS'), frac_insol, start=(/1,1,1/), kount=(/pcols,pver,pcnst/) )
   call pbuf_get_field(pbuf, pbuf_get_index('MU_CRM'), mu )
   call pbuf_get_field(pbuf, pbuf_get_index('MD_CRM'), md )
   call pbuf_get_field(pbuf, pbuf_get_index('DU_CRM'), du )
   call pbuf_get_field(pbuf, pbuf_get_index('EU_CRM'), eu )
   call pbuf_get_field(pbuf, pbuf_get_index('ED_CRM'), ed )
   call pbuf_get_field(pbuf, pbuf_get_index('JT_CRM'), cld_top_idx_ptr )
   call pbuf_get_field(pbuf, pbuf_get_index('MX_CRM'), cld_bot_idx_ptr )

   !!! get integer cloud indices - add small constant since int() rounds down to nearest integer
   ! cld_top_idx = int( cld_top_idx_flt + 0.1_r8)
   ! cld_bot_idx = int( cld_bot_idx_flt + 0.1_r8)
   cld_top_idx = cld_top_idx_ptr
   cld_bot_idx = cld_bot_idx_ptr

   !!! initialize dp and dpdry - used for tracers of dry mixing ratio type
   dpdry = 0._r8

   !!! convert pressure thickness from Pa to hPa
   dpdry(1:ncol,:) = state%pdeldry(1:ncol,:)/100._r8
   dp   (1:ncol,:) = state%pdel   (1:ncol,:)/100._r8

   !-------------------------------------------------------------------------------------
   ! Calculate the transport tendencies
   !-------------------------------------------------------------------------------------
   call crm_bulk_transport_tend (lchnk, ncol, pcnst, lq, state%q,    &
                                 dry_const_flag, frac_insol,         &
                                 mu, md, du, eu, ed, dp, dpdry,      &
                                 cld_top_idx, cld_bot_idx,           &
                                 ptend%q )   
   !-------------------------------------------------------------------------------------
   !-------------------------------------------------------------------------------------

end subroutine crm_bulk_transport

!==================================================================================================
!==================================================================================================

subroutine crm_bulk_transport_tend( lchnk, ncol, ncnst, do_transport, q, &
                                    dry_const_flag, frac_insol,          &
                                    mu, md, du, eu, ed, dp, dpdry,       &
                                    cld_top_idx, cld_bot_idx,            &
                                    q_tend_out )
   !-------------------------------------------------------------------------------------
   ! Purpose: to do convective transport of tracer species using CRM cloud statistics
   ! Walter Hannah (LLNL), 2018: based on convtran() from the Zhang-McFarlane scheme
   !-------------------------------------------------------------------------------------
   use ppgrid
   ! use cam_abortutils,  only: endrun
   use shr_kind_mod,    only: r8 => shr_kind_r8

   implicit none
   
   !!! Input Arguments
   integer,  intent(in) :: lchnk                                        ! chunk identifier
   integer,  intent(in) :: ncol                                         ! number of columns
   integer,  intent(in) :: ncnst                                        ! number of tracers to transport
   logical,  intent(in), dimension(ncnst)            :: do_transport    ! flag for doing convective transport
   real(r8), intent(in), dimension(pcols,pver,ncnst) :: q               ! Tracer array including moisture
   logical,  intent(in), dimension(pcols,ncnst)      :: dry_const_flag  ! flag to indicate dry constituent
   real(r8), intent(in), dimension(pcols,pver,ncnst) :: frac_insol      ! fraction of tracer that is insoluble
   real(r8), intent(in), dimension(pcols,pver)       :: mu              ! updraft mass flux     [mb/s]
   real(r8), intent(in), dimension(pcols,pver)       :: eu              ! updraft entrainment
   real(r8), intent(in), dimension(pcols,pver)       :: du              ! updraft detrainment
   real(r8), intent(in), dimension(pcols,pver)       :: md              ! downdraft mass flux   [mb/s]
   real(r8), intent(in), dimension(pcols,pver)       :: ed              ! downdraft entrainment
   real(r8), intent(in), dimension(pcols,pver)       :: dp              ! Delta pressure between interfaces
   real(r8), intent(in), dimension(pcols,pver)       :: dpdry           ! Delta pressure between interfaces - for dry constituents
   integer,  intent(in), dimension(pcols)            :: cld_top_idx     ! Index of cloud top    for each column
   integer,  intent(in), dimension(pcols)            :: cld_bot_idx     ! Index of cloud bottom for each column
   
   !!! Output Arguments
   real(r8), intent(out), dimension(pcols,pver,ncnst) :: q_tend_out       ! Tracer tendency array

   !!! Local Variables
   integer i,m,k             ! Work indices
   integer km1,kp1           ! k-1, k+1

   real(r8) q_above          ! Mix ratio of constituent above
   real(r8) q_below          ! Mix ratio of constituent below
   real(r8) q_diff_norm      ! Normalized diff between q_above and q_below

   real(r8), parameter :: q_small          = 1.e-36_r8   ! small constituent mixing ratio threshold
   real(r8), parameter :: min_mf           = 1.e-15_r8   ! threshold below which we treat the mass fluxes as zero (in mb/s)
   real(r8), parameter :: min_q_diff       = 1.E-6_r8    ! minimum q_diff_norm
   real(r8), parameter :: max_q_factor     = 1.e-12_r8   ! scaling factor for max_q
   real(r8), parameter :: max_flux_factor  = 1.e-12_r8   ! scaling factor for max flux limiter

   real(r8) min_q            ! temporary variable used to determine q_diff_norm
   real(r8) max_q            ! temporary variable used to determine q_diff_norm
   real(r8) flux_in          ! flux into layer k
   real(r8) flux_out         ! flux out of layer k
   real(r8) flux_net         ! net flux for layer k

   real(r8), dimension(pver) :: mu_p_duxdp   ! updraft mass flux + detrainment * dpdry

   real(r8), dimension(pver) :: q_intrp      ! constituent mixing ratio in env       at interfaces
   real(r8), dimension(pver) :: q_up         ! constituent mixing ratio in updraft   at interfaces
   real(r8), dimension(pver) :: q_dn         ! constituent mixing ratio in downdraft at interfaces

   real(r8), dimension(pver) :: eu_tmp       ! Mass entraining from updraft
   real(r8), dimension(pver) :: du_tmp       ! Mass detraining from updraft
   real(r8), dimension(pver) :: ed_tmp       ! Mass entraining from downdraft
   real(r8), dimension(pver) :: dp_tmp       ! Delta pressure between interfaces
   !-------------------------------------------------------------------------------------
   !-------------------------------------------------------------------------------------
   q_tend_out(:,:,:) = 0._r8   ! Initialize output tendency array   
   !-------------------------------------------------------------------------------------
   ! Loop ever each constituent
   !-------------------------------------------------------------------------------------
   do m = 1, ncnst
      if ( do_transport(m) ) then
         do i = 1,ncol
            !--------------------------------------------------------
            ! prepare data for transport tendency calculation
            !--------------------------------------------------------
            do k = 1,pver

               !!! make adjustment for for "dry" constituents
               if ( dry_const_flag(i,m) ) then
                  dp_tmp(k) = dpdry(i,k)
                  du_tmp(k) = du(i,k)*dp(i,k)/dpdry(i,k)
                  eu_tmp(k) = eu(i,k)*dp(i,k)/dpdry(i,k)
                  ed_tmp(k) = ed(i,k)*dp(i,k)/dpdry(i,k)
               else
                  dp_tmp(k) = dp(i,k)
                  du_tmp(k) = du(i,k)
                  eu_tmp(k) = eu(i,k)
                  ed_tmp(k) = ed(i,k)
               endif

               !!! calculate q_diff_norm to determine averaging prodcedure
               km1 = max(1,k-1)
               min_q = min( q(i,km1,m), q(i,k,m) )
               max_q = max( q(i,km1,m), q(i,k,m) )
               if (min_q < 0) then
                  q_diff_norm = 0._r8
               else
                  q_diff_norm = abs( q(i,k,m) - q(i,km1,m) ) / max( max_q, q_small )
               endif
               
               !!! Interpolate tracer values to interfaces
               if ( q_diff_norm > min_q_diff ) then
                  !!! If the two layers differ significantly use a geometric mean
                  q_above    = max( q(i,km1,m), max_q*max_q_factor )
                  q_below    = max( q(i,k  ,m), max_q*max_q_factor )
                  q_intrp(k) = log(q_above/q_below) /(q_above-q_below) *q_above*q_below
               else             
                  !!! If the two layers are sufficiently small just use arithmetic mean
                  q_intrp(k) = 0.5_r8* ( q(i,k,m) + q(i,km1,m) )
               end if
               !!! provisional up and down draft tracer values
               q_up(k) = q_intrp(k)
               q_dn(k) = q_intrp(k)

               !!! set mu_p_duxdp (updraft mass flux + detrainment * dpdry)
               mu_p_duxdp(k) = mu(i,k) + du_tmp(k)*dp_tmp(k)

            end do ! k
            !--------------------------------------------------------
            ! Determine updraft tracer values
            !--------------------------------------------------------
            !!! set initial updraft tracer values at surface
            k = pver 
            if ( mu_p_duxdp(k) > min_mf ) then
               q_up(k) = ( eu_tmp(k) * frac_insol(i,k,m) * q(i,k,m) * dp_tmp(k) )/mu_p_duxdp(k)
            endif

            !!! Updraft from bottom to top
            do k = pver-1,1,-1
               kp1 = min(pver,k+1)
               if (mu_p_duxdp(k) > min_mf) then
                  q_up(k) = (  mu(i,kp1) * q_up(kp1)          &
                              + eu_tmp(k) * q(i,k,m) * frac_insol(i,k,m) * dp_tmp(k) &
                             )/mu_p_duxdp(k)
               endif
            end do
            !--------------------------------------------------------
            ! Determine downdraft tracer values
            !--------------------------------------------------------
            !!! set initial downdraft tracer values at top of model
            k = 2
            km1 = 1
            if ( md(i,k) < -min_mf ) then
               q_dn(k) =  ( -ed_tmp(km1) * frac_insol(i,km1,m) * q(i,km1,m) * dp_tmp(km1) )/md(i,k)
            endif

            !!! Downdraft from top to bottom
            do k = 3,pver
               km1 = max(1,k-1)
               if (md(i,k) < -min_mf) then
                  q_dn(k) =  (  md(i,km1)   * q_dn(km1)            &
                              - ed_tmp(km1) * q(i,km1,m) * frac_insol(i,km1,m) * dp_tmp(km1) &
                             )/md(i,k)
               endif
            end do ! k
            !--------------------------------------------------------
            ! Calculate fluxes of q by convective transport
            !--------------------------------------------------------
            do k = cld_top_idx(i),pver
               km1 = max(1,k-1)
               kp1 = min(pver,k+1)

               !!! flux into layer k
               flux_in =   mu(i,kp1)*q_up(kp1) + mu(i,k  )*min(q_intrp(k  ),q(i,km1,m)) &
                        -( md(i,k  )*q_dn(k  ) + md(i,kp1)*min(q_intrp(kp1),q(i,kp1,m)) )

               !!! flux out of layer k
               flux_out =   mu(i,k  )*q_up(k  ) + mu(i,kp1)*min(q_intrp(kp1),q(i,k,m)) &
                         -( md(i,kp1)*q_dn(kp1) + md(i,k  )*min(q_intrp(k  ),q(i,k,m)))

               ! limit fluxes outside convection to mass in appropriate layer
               ! limiters are "probably" only safe for positive definite quantitities
               ! it assumes that updrarft (mu) and dowdraft (md) mass fluxes 
               ! already satisfy a courant number limit of 1

               flux_net = flux_in - flux_out
               if ( abs(flux_net) < ( max(flux_in,flux_out)*max_flux_factor ) ) then
                  flux_net = 0._r8
               endif

               q_tend_out(i,k,m) = flux_net/dp_tmp(k)

            end do ! k
            !--------------------------------------------------------
            ! calculate fluxes for lowest cloud layer
            !--------------------------------------------------------
            do k = cld_bot_idx(i),pver
               km1 = max(1,k-1)
               if ( k == cld_bot_idx(i) ) then

                  !!! calculate fluxes in and out of cloud base
                  flux_in  = mu(i,k)*min(q_intrp(k),q(i,km1,m)) - md(i,k)*q_dn(k)
                  flux_out = mu(i,k)*q_up(k) - md(i,k)*min(q_intrp(k),q(i,k,m))

                  !!! limit fluxes similar as above
                  flux_net = flux_in - flux_out
                  if ( abs(flux_net) < max(flux_in,flux_out)*max_flux_factor ) then
                     flux_net = 0._r8
                  endif

                  q_tend_out(i,k,m) = flux_net/dp_tmp(k)

               else if ( k > cld_bot_idx(i) ) then
                  !!! otherwise tendency is zero
                  q_tend_out(i,k,m) = 0._r8
               end if
            end do ! k
            !--------------------------------------------------------
            !--------------------------------------------------------
         end do ! i=1,ncol
      end if ! do_transport
   end do ! m=1,pcnst
   !-------------------------------------------------------------------------------------
   !-------------------------------------------------------------------------------------
end subroutine crm_bulk_transport_tend

!==================================================================================================
!==================================================================================================

subroutine crm_bulk_aero_mix_nuc( state, ptend, pbuf, dtime    &
                                  ! ,cflx, ww_cen, ww_bnd        &
                                  ,pblht, species_class )
   !-------------------------------------------------------------------------------------
   ! Purpose: to calculate aerosol tendency from droplet activation and mixing
   ! Walter Hannah (LLNL), 2018: based on Minghuai Wang's crmclouds_mixnuc_tend()
   !                             Adopted from mmicro_pcond in cldwat2m.F90
   !-------------------------------------------------------------------------------------
   use physics_types,    only: physics_state, physics_ptend, physics_tend, physics_ptend_init
   use physics_buffer,   only: physics_buffer_desc, pbuf_old_tim_idx, pbuf_get_index, pbuf_get_field
   use physconst,        only: gravit, rair, karman, spec_class_gas
   use constituents,     only: cnst_get_ind, pcnst
   use time_manager,     only: is_first_step
   use cam_history,      only: outfld
   use ndrop,            only: dropmixnuc
   use modal_aero_data,  only: ntot_amode, numptr_amode, lmassptr_amode, nspec_amode
   use rad_constituents, only: rad_cnst_get_info

   !!! Input 
   type(physics_state), intent(in)                  :: state         ! state variable
   type(physics_buffer_desc), pointer, dimension(:) :: pbuf          ! physics buffer
   real(r8), intent(in)                             :: dtime         ! timestep for microphysics - passed to dropmixnuc()
   real(r8), intent(in), dimension(pcols)           :: pblht         ! PBL height - used for mixing length calculation [m]
   ! real(r8), intent(in), dimension(pcols,pcnst)     :: cflx          ! constituent flux from surface (not used)
   ! real(r8), intent(in), dimension(pcols, pver)     :: ww_cen        ! vert. velocity variance (at layer center)    [m2/s2]
   ! real(r8), intent(in), dimension(pcols, pver+1)   :: ww_bnd        ! vert. velocity variance (at layer interface) [m2/s2]
   integer,  intent(in), dimension(pcnst)           :: species_class ! 

   !!! output arguments
   type(physics_ptend), intent(out) :: ptend   ! output physics tendencies

   !!! Local variables
   integer i, k, m, k1, k2
   ! integer ifld, itim
   integer ixcldliq, ixcldice, ixnumliq
   integer l,lnum,lnumcw,lmass,lmasscw
   integer :: lchnk                    ! chunk identifier
   integer :: ncol                     ! number of atmospheric columns
   integer :: nmodes
   logical :: do_mmf = .true.          ! value insignificant, if present, means that dropmixnuc is called the mmf part. 

!TODO: calculate the vertical velocites using crm grid data on the pbuf
   real(r8), dimension(pcols,pver ) :: ww_cen         ! vert. velocity variance (at layer center)    [m2/s2]
   real(r8), dimension(pcols,pverp) :: ww_bnd         ! vert. velocity variance (at layer interface) [m2/s2]
   ! real(r8), dimension(pcols,pver ) :: nc             ! droplet number concentration                 [#/kg]
   ! real(r8), dimension(pcols,pver ) :: nctend         ! change in droplet number concentration
   ! real(r8), dimension(pcols,pver ) :: omega          ! grid-averaaged vertical velocity - not used?
   real(r8), dimension(pcols,pver ) :: qc             ! liquid water content                          [kg/kg]
   real(r8), dimension(pcols,pver ) :: qi             ! ice water content                             [kg/kg] 
   real(r8), dimension(pcols,pver ) :: lcldn          ! new liquid cloud fraction (current time step)
   real(r8), dimension(pcols,pver ) :: lcldo          ! old liquid cloud fraction (previous time step)
   real(r8), dimension(pcols,pver ) :: wsub           ! subgrid vertical velocity
   real(r8), dimension(pcols,pver ) :: dz_inv         ! inverse of distance between levels (meter)
   real(r8), dimension(pcols,pver ) :: dz             ! layer depth                                   [m]
   real(r8), dimension(pcols,pver ) :: cs             ! air density
   real(r8), dimension(pcols,pverp) :: ekd_crm        ! diffusivity
   real(r8), dimension(pcols,pverp) :: kkvh_crm       ! eddy diffusivity
   real(r8), dimension(pcols,pverp) :: mix_len        ! mixing length                                 [m]
   real(r8), dimension(pcols,pverp) :: zheight        ! height at lay interface                       [m]
   real(r8), dimension(pcols,pverp) :: asym_length    ! asymptotic length scale                       [m]
   real(r8), dimension(pcols,pver ) :: tend_ndrop     ! tendency of cloud droplet number concentrations (not used in the MMF) 
   
   real(r8), allocatable, dimension(:,:,:) :: aero_act_frac_num ! activation fraction for aerosol number

   real(r8) :: qcld                                   ! temporary cloud condensate - liquid + ice

   !!! physics buffer variables
   real(r8), pointer, dimension(:,:)   :: cldn        ! cloud fraction                  (current time step)
   real(r8), pointer, dimension(:,:)   :: cldo        ! cloud fraction                  (previous time step)
   real(r8), pointer, dimension(:,:)   :: acldy_cen   ! liquid cloud fraction from ECPP (previous time step )
   real(r8), pointer, dimension(:,:)   :: kkvh        ! vertical diffusivity
   real(r8), pointer, dimension(:,:)   :: tke         ! turbulence kenetic energy 
   real(r8), pointer, dimension(:,:)   :: tk_crm      ! m2/s
   logical,           dimension(pcnst) :: lq          ! flag for ptend

   !!! fixed parameters
   real(r8), parameter :: qcld_threshold = 1.e-18_r8  ! lower bound for cloud condesnate
   real(r8), parameter :: wsub_min       =  0.2_r8    ! lower bound for subgrid vertical velocity
   real(r8), parameter :: wsub_max       = 10.0_r8    ! upper bound for subgrid vertical velocity
   !-------------------------------------------------------------------------------------
   !-------------------------------------------------------------------------------------
   lchnk = state%lchnk
   ncol  = state%ncol

   call rad_cnst_get_info(0, nmodes = nmodes )

   allocate( aero_act_frac_num(pcols,pver,nmodes) )

   !--------------------------------------------------------
   ! Initialize ptend
   !--------------------------------------------------------
   lq(:)=.false.

   do m = 1,ntot_amode
      lnum = numptr_amode(m)
      if ( lnum > 0 ) then
         lq(lnum) = .true.
      end if
      do l = 1,nspec_amode(m)
         lmass = lmassptr_amode(l,m)
         lq(lmass) = .true.
      end do
   end do ! m=1,ntot_amode

   call physics_ptend_init(ptend,state%psetcols,'crm_bulk_aero_mix_nuc', lq=lq)

   !!! In the MMF model, turbulent mixing for tracer species are turned off in tphysac.
   !!! So the turbulent for gas species mixing are added here.
   do m = 1,pcnst
      if ( species_class(m) .eq. spec_class_gas ) then
         ptend%lq(m) = .true.
      end if
   end do

   !--------------------------------------------------------
   ! get pbuf variables
   !--------------------------------------------------------
   call pbuf_get_field(pbuf,pbuf_get_index('CLD'),       cldn      )
   call pbuf_get_field(pbuf,pbuf_get_index('CLDO'),      cldo      )
!TODO: need a different version of cloud fraction derived from CRM data directly
   call pbuf_get_field(pbuf,pbuf_get_index('ACLDY_CEN'), acldy_cen )    
   call pbuf_get_field(pbuf,pbuf_get_index('kvh'),       kkvh      )
   call pbuf_get_field(pbuf,pbuf_get_index('tke'),       tke       )
   call pbuf_get_field(pbuf,pbuf_get_index('TK_CRM'),    tk_crm    )

   !!! initialize diffusivity and tke
   if (is_first_step()) then
      kkvh(:,:) = 0.0_r8
      tke(:,:)  = 0.0_r8
   endif

   !--------------------------------------------------------
   ! Mixing length calculation
   !--------------------------------------------------------
   do i=1,ncol
      do k=1,pver-1
         dz_inv(i,k) = 1._r8/(state%zm(i,k)-state%zm(i,k+1))
      end do ! k=1,pver-1
      dz_inv(i,pver) = dz_inv(i,pver-1)

      !!! calculate height at layer interface (simple calculation)
      zheight(i,pverp) = 0.0
      do k=pver,1,-1
         zheight(i,k) = zheight(i,k+1) + state%pdel(i,k)/state%pmid(i,k)*(rair*state%t(i,k)/gravit)
      end do ! k=pver,1,-1

      !!! calculate mixing length - Holtslag and Boville, 1993, J. Climate. 
      do k=1,pverp
         if ( zheight(i,k).le.pblht(i) ) then
            asym_length(i,k) = 300._r8
         else
            asym_length(i,k) = 30._r8 + 270._r8 * exp( 1._r8 - zheight(i,k) / pblht(i) )
         endif
         mix_len(i,k) = asym_length(i,k)*karman*zheight(i,k)/(asym_length(i,k)+karman*zheight(i,k))
      end do ! k=1,pverp
   end do ! i=1,ncol

   !!! output mixing length?
   ! call outfld('LENGC'      , mix_len, pcols, lchnk)
   ! call outfld('CRM_MIXLENC', mix_len, pcols, lchnk)

   !--------------------------------------------------------
   ! Calculate vert. velocity and eddy diffusivity
   !--------------------------------------------------------
   kkvh_crm = 0._r8
   do i=1,ncol

      !!! Calculate subgrid vertical velocity variance
      do k=1,pver
         ! wsub(i,k) = sqrt(tke(i,k)/3.)  

         !!! tke from CRM is located in the middle of model level, should be tke or tke/3. 
         !!! in cldwat2m.F90, it is tke. wsub seems too large from this approach. 

         ! wsub(i,k) = tk_crm(i,k) * dz_inv(i,k)

         !!! vertical variance in the quiescent class excludes contribution from strong updraft and downdraft. 

         ! wsub(i,k) = sqrt(wwqui_cen(i,k))           ! variance in quiescent area
         ! wsub(i,k) = sqrt(wwqui_cloudy_cen(i,k))    ! variance in cloudy quiescent area  (this version is used for ECPP case)

         wsub(i,k) = sqrt( ww_cen(i,k) )

         !!! from tke in CAM
         ! wsub(i,k) = sqrt(0.5_r8*(tke(i,k)+tke(i,k+1)))

         wsub(i,k) = min( wsub(i,k), wsub_max )
         wsub(i,k) = max( wsub(i,k), wsub_min )
      end do ! k=1,pver

      !!! calculate eddy diffusivity (ekd_crm) from subgrid vertical velocity (wsub) 
      do k=1,pver+1
         k1=min(k, pver)
         k2=max(k-1, 1)

         !!! following ndrop.F90
         ! ekd_crm(i,k) = wsub(i,k) / dz_inv(i,k)

         !!! calculate wsub at boundary
         ! wsub_bnd(i,k) = sqrt(wwqui_cloudy_bnd(i,k))
         wsub_bnd(i,k) = sqrt( ww_bnd(i,k) )

         wsub_bnd(i,k) = min( wsub_bnd(i,k), wsub_max )
         wsub_bnd(i,k) = max( wsub_bnd(i,k), wsub_min )

         !!! use cloudy quiescent class - use wsub at layer boundary - large ekd in free troposphere
         ! ekd_crm(i,k) = wsub_bnd(i,k) * 2.0 / (dz_inv(i,k1)+dz_inv(i,k2))  
         ekd_crm(i,k) = wsub_bnd(i,k) * mix_len(i,k) 
         kkvh_crm(i,k) = ekd_crm(i,k)

         !!! kkvh_crm is from kvh in CAM
         ! kkvh_crm(i,k) = kkvh(i,k)

         !!! set kkvh to kkvh_crm so this will be used in dropmixnuc in the mmf call
         kkvh(i,k) = kkvh_crm(i,k)

      end do ! k=1,pver

   end do ! i=1,ncol

   !--------------------------------------------------------
   ! Adjust new and old liquid cloud fraction
   !--------------------------------------------------------
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)
   ! call cnst_get_ind('NUMLIQ', ixnumliq)

   qc(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq)
   qi(:ncol,:pver) = state%q(:ncol,:pver,ixcldice)
   ! nc(:ncol,:pver) = state%q(:ncol,:pver,ixnumliq)

   do k=1,pver
      do i=1,ncol

         qcld = qc(i,k) + qi(i,k)

         !!! check if total cloud condensate is sufficiently large
         if ( qcld.gt.qcld_threshold ) then
            lcldn(i,k) = cldn(i,k)*qc(i,k)/qcld
            lcldo(i,k) = cldo(i,k)*qc(i,k)/qcld
         else
            lcldn(i,k) = 0._r8
            lcldo(i,k) = 0._r8
         end if

      end do ! i=1,ncol
   end do ! k=1,pver

   !!! should we set large-scale omega to be zero ??
   ! omega(:ncol,:) = state%omega(:ncol,:)

   !--------------------------------------------------------
   ! Calculate mixing and nucleation tendencies
   !--------------------------------------------------------
   
   call dropmixnuc(state, ptend, dtime, pbuf,      &
                   wsub, lcldn, lcldo,             &
                   tend_ndrop, aero_act_frac_num,  &
                   species_class, do_mmf )
   
   !--------------------------------------------------------
   ! clean-up
   !--------------------------------------------------------
   deallocate(aero_act_frac_num)

end subroutine crm_bulk_aero_mix_nuc

!==================================================================================================
!==================================================================================================
! subroutine crm_bulk_diagnose_cloudy_clear()
   !-------------------------------------------------------------------------------------
   ! Purpose: diagnose statistics of cloudy and clear classes from CRM 
   !          for input data to crm_bulk_aero_mix_nuc()
   ! Walter Hannah (LLNL), 2018
   !-------------------------------------------------------------------------------------

   !!! initialize running average quantities when CRM starts? (or use separate subroutine) 

   !!! determine clear/cloudy staet of each CRM gridcell

   !!! determine current vertical velocity variance

   !!! add stats to running average variables

! end subroutine 

!==================================================================================================
!==================================================================================================

#endif /* CRM */
#endif /* SP_CRM_BULK */

end module crm_bulk_mod
