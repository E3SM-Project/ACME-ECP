module crm_bulk_transport_mod

! the check for "CBT" just allows the module to be hidden
#ifdef CBT
#if (defined CRM)

   !------------------------------------------------------------------------
   ! Purpose: 
   !  Provides a method for convective transport of tracers using 
   !  cloud statistics from the explicit convection of the CRM.
   !  This is meant to be an interim alternative to the ECPP 
   !  (explicit-cloud parameterized-pollutant) approach. 
   !
   !  Revision history: 
   !     2018: Walter Hannah - initial version based on Minghuai Wang's "crmclouds_camaerosols" module
   ! 
   !-------------------------------------------------------------------------------------------- 

   use shr_kind_mod,       only: r8 => shr_kind_r8
   use cam_abortutils,     only: endrun
   use ppgrid
   use params

   implicit none

   public :: crm_bulk_transport
   ! public :: crm_bulk_transport_tend 
   ! public :: crm_bulk_cloud_diag 

   contains 

!======================================================================================================
!======================================================================================================

subroutine crm_bulk_transport(state,  ptend,  ztodt,  pbuf)
   !-----------------------------------------------------------------------------------------
   ! Purpose:  to do convective transport of tracer species using the 
   !           cloud updraft and downdraft mass fluxes from the CRM.
   !           Adopted from crmclouds_convect_tend() by Minghuai Wang
   !           and convtran() from the Zhang-McFarlane scheme.
   !           Original version used the "gathered" array terminology 
   !           but for the CRM we omit this and operate on all columns.
   ! Walter Hannah (LLNL), 2018
   !-----------------------------------------------------------------------------------------
   use physics_types, only: physics_state, physics_ptend, physics_ptend_init
   use time_manager,  only: get_nstep
   use physics_buffer, only: physics_buffer_desc, pbuf_old_tim_idx, pbuf_get_index, pbuf_get_field
   use constituents,  only: pcnst, cnst_get_ind
   use zm_conv,       only: convtran
   use error_messages, only: alloc_err  

   !!! Input Arguments
   type(physics_state), intent(in ) :: state       ! Physics state variables
   real(r8),            intent(in ) :: ztodt       ! GCM time step

   !!! Output Arguments
   type(physics_ptend), intent(out)   :: ptend     ! indivdual parameterization tendencies
   type(physics_buffer_desc), pointer :: pbuf(:)   ! physics buffer

   !!! Local variables
   integer :: i, lchnk, ncol
   integer :: ixcldice, ixcldliq                   ! constituent indices for cloud liquid and ice water.
   real(r8), dimension(pcols,pver) :: dpdry        ! ?
   real(r8), dimension(pcols,pver) :: dp           ! layer thickness in mbs (between upper/lower interface).

   !!! physics buffer fields 
   integer itim, ifld
   real(r8), pointer, dimension(:,:,:) :: fracis   ! fraction of transported species that are insoluble
   real(r8), pointer, dimension(:,:)   :: mu       ! updraft mass flux       (pcols,pver,begchunk:endchunk)
   real(r8), pointer, dimension(:,:)   :: eu       ! updraft entrainment     (pcols,pver,begchunk:endchunk)
   real(r8), pointer, dimension(:,:)   :: du       ! updraft detrainment     (pcols,pver,begchunk:endchunk)
   real(r8), pointer, dimension(:,:)   :: md       ! downdraft mass flux     (pcols,pver,begchunk:endchunk)
   real(r8), pointer, dimension(:,:)   :: ed       ! downdraft entrainment   (pcols,pver,begchunk:endchunk)

   real(r8), pointer, dimension(:) :: cld_top_idx_flt    ! index of cloud top      (pcols,begchunk:endchunk)
   real(r8), pointer, dimension(:) :: cld_bot_idx_flt    ! index of cloud bottom   (pcols,begchunk:endchunk)

   integer, dimension(pcols) :: cld_top_idx              ! index of cloud top
   integer, dimension(pcols) :: cld_bot_idx              ! index of cloud bottom
   logical, dimension(pcnst) :: lq  

   !-----------------------------------------------------------------------------------------
   !-----------------------------------------------------------------------------------------
   
   lchnk = state%lchnk
   ncol  = state%ncol

   !!! Initialize ptend
   ! transport of all trace species except 
   ! cloud liquid and cloud ice is done here 
   ! because we need to do the scavenging first 
   ! to determine the interstitial fraction.
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)
   lq(:) = .true.
   lq(1)        = .false.
   lq(ixcldice) = .false.
   lq(ixcldliq) = .false.

   call physics_ptend_init(ptend,state%psetcols,'crm_bulk_transport',lq=lq)

   !!! Associate pointers with physics buffer fields
   call pbuf_get_field(pbuf, pbuf_get_index('FRACIS'), fracis, start=(/1,1,1/), kount=(/pcols,pver,pcnst/) )
   call pbuf_get_field(pbuf, pbuf_get_index('MU_CRM'), mu )
   call pbuf_get_field(pbuf, pbuf_get_index('MD_CRM'), md )
   call pbuf_get_field(pbuf, pbuf_get_index('DU_CRM'), du )
   call pbuf_get_field(pbuf, pbuf_get_index('EU_CRM'), eu )
   call pbuf_get_field(pbuf, pbuf_get_index('ED_CRM'), ed )
   call pbuf_get_field(pbuf, pbuf_get_index('JT_CRM'), cld_top_idx_flt )
   call pbuf_get_field(pbuf, pbuf_get_index('MX_CRM'), cld_bot_idx_flt )


   !!! get integer cloud indices - add small constant since int() rounds down to nearest integer
   cld_top_idx = int( cld_top_idx_flt + 0.1_r8)
   cld_bot_idx = int( cld_bot_idx_flt + 0.1_r8)

   !!! initialize dpdry for call to convtran - it is used for tracers of dry smixing ratio type
   dpdry = 0._r8
   do i = 1,ncol
      dpdry(i,:) = state%pdeldry(i,:)/100._r8
      dp(i,:)    = state%pdel(i,:)/100._r8
   end do

   !-----------------------------------------------------------------------------------------
   !!! Calculate the transport tendencies
   !-----------------------------------------------------------------------------------------
   

   call crm_bulk_transport_tend (lchnk, pcnst, lq,                    &
                                 state%q, fracis, mu, md, du, eu, ed, &
                                 dp, dpdry, cld_top_idx, cld_bot_idx, &
                                 ptend%q )   
   
   !-----------------------------------------------------------------------------------------
   !-----------------------------------------------------------------------------------------

end subroutine crm_bulk_transport

!=====================================================================================================
!======================================================================================================

subroutine crm_bulk_transport_tend( lchnk, ncnst, do_transport,            &
                                    q, fracis, mu, md, du, eu, ed,       &
                                    dp, dpdry, cld_top_idx, cld_bot_idx, &
                                    dqdt )
   !-----------------------------------------------------------------------------------------
   ! Purpose:  to do convective transport of tracer species using CRM cloud statistics
   ! Walter Hannah (LLNL), 2018: based on convtran() from the Zhang-McFarlane scheme
   !-----------------------------------------------------------------------------------------
   use ppgrid
   use shr_kind_mod,    only: r8 => shr_kind_r8
   use constituents,    only: cnst_get_type_byind
   use cam_abortutils,  only: endrun

   implicit none
   
   !!! Input Arguments
   integer,  intent(in) :: lchnk                                     ! chunk identifier
   integer,  intent(in) :: ncnst                                     ! number of tracers to transport
   logical,  intent(in), dimension(ncnst)            :: do_transport ! flag for doing convective transport
   real(r8), intent(in), dimension(pcols,pver,ncnst) :: q            ! Tracer array including moisture
   real(r8), intent(in), dimension(pcols,pver,ncnst) :: frac_insol   ! fraction of tracer that is insoluble

   real(r8), intent(in), dimension(pcols,pver)       :: mu           ! updraft mass flux     [mb/s]
   real(r8), intent(in), dimension(pcols,pver)       :: eu           ! updraft entrainment
   real(r8), intent(in), dimension(pcols,pver)       :: du           ! updraft detrainment
   real(r8), intent(in), dimension(pcols,pver)       :: md           ! downdraft mass flux   [mb/s]
   real(r8), intent(in), dimension(pcols,pver)       :: ed           ! downdraft entrainment
   
   real(r8), intent(in), dimension(pcols,pver)       :: dp           ! Delta pressure between interfaces
   real(r8), intent(in), dimension(pcols,pver)       :: dpdry        ! Delta pressure between interfaces

   integer,  intent(in), dimension(pcols)            :: cld_top_idx  ! Index of cloud top for each column
   integer,  intent(in), dimension(pcols)            :: cld_bot_idx  ! Index of cloud top for each column
   
   !!! Output Arguments
   real(r8), intent(out), dimension(pcols,pver,ncnst) :: tend_out    ! Tracer tendency array

   !!! Local Variables

   integer i,,m,k            ! Work indices
   integer kbm               ! Highest altitude index of cloud base
   integer ktm               ! Highest altitude index of cloud top
   integer kk,kkp1,km1,kp1   ! Work index

   real(r8) q_above          ! Mix ratio of constituent above
   real(r8) q_below          ! Mix ratio of constituent below
   real(r8) q_diff_norm      ! Normalized diff between q_above and q_below
   real(r8) small            ! A small number
   real(r8) mbsth            ! Threshold for mass fluxes

   real(r8) mupdudp          ! A work variable
   real(r8) minc             ! A work variable
   real(r8) maxc             ! A work variable
   real(r8) fluxin           ! A work variable
   real(r8) fluxout          ! A work variable
   real(r8) netflux          ! A work variable


   real(r8), dimension(pver) :: q_int        ! constituent mixing ratio in env       at interfaces
   real(r8), dimension(pver) :: q_up         ! constituent mixing ratio in updraft   at interfaces
   real(r8), dimension(pver) :: q_dn         ! constituent mixing ratio in downdraft at interfaces
   real(r8), dimension(pver) :: q_tend       ! provisional tendencies

   real(r8), dimension(pver) :: eu_tmp       ! Mass entraining from updraft
   real(r8), dimension(pver) :: du_tmp       ! Mass detraining from updraft
   real(r8), dimension(pver) :: ed_tmp       ! Mass entraining from downdraft
   real(r8), dimension(pver) :: dp_tmp       ! Delta pressure between interfaces

   !-----------------------------------------------------------------------------------------
   !-----------------------------------------------------------------------------------------
   
   !!! Initialize output tendencies
   tend_out(:,:,:) = 0._r8

   small = 1.e-36_r8    ! threshold for ???????
   mbsth = 1.e-15_r8    ! threshold below which we treat the mass fluxes as zero (in mb/s)

   !!! Find the highest level top and bottom levels of convection
   ktm = pver
   kbm = pver
   do i = 1,ncol
      ktm = min( ktm, cld_top_idx(i) )
      kbm = min( kbm, cld_bot_idx(i) )
   end do

   !-----------------------------------------------------------------------------------------
   !!! Loop ever each constituent
   !-----------------------------------------------------------------------------------------
   do m = 1, pcnst
      if ( do_transport(m) ) then
         do i = 1,ncol

            !!! check for "dry" constituents
            if (cnst_get_type_byind(m).eq.'dry') then
               do k = 1,pver
                  dp_tmp(i,k) = dpdry(i,k)
                  du_tmp(i,k) = du(i,k)*dp(i,k)/dpdry(i,k)
                  eu_tmp(i,k) = eu(i,k)*dp(i,k)/dpdry(i,k)
                  ed_tmp(i,k) = ed(i,k)*dp(i,k)/dpdry(i,k)
               end do
            else
               do k = 1,pver
                  dp_tmp(i,k) = dp(i,k)
                  du_tmp(i,k) = du(i,k)
                  eu_tmp(i,k) = eu(i,k)
                  ed_tmp(i,k) = ed(i,k)
               end do
            endif

            !!! Interpolate environment tracer values to interfaces
            do k = 1,pver
               km1 = max(1,k-1)
               
               minc = min( q(i,km1,m), q(i,k,m) )
               maxc = max( q(i,km1,m), q(i,k,m) )
               if (minc < 0) then
                  q_diff_norm = 0._r8
               else
                  q_diff_norm = abs( q(i,k,m) - q(i,km1,m) ) / max( maxc, small )
               endif

               !!! If the two layers differ significantly use a geometric averaging procedure
               if (q_diff_norm > 1.E-6_r8) then
                  q_above = max( q(i,km1,m), maxc*1.e-12_r8 )
                  q_below = max( q(i,k  ,m), maxc*1.e-12_r8 )
                  q_int(i,k) = log(q_above/q_below) /(q_above-q_below) *q_above*q_below
               else             
                  !!! otherwise just use arithmetic mean
                  q_int(i,k) = 0.5_r8* ( q(i,k,m) + q(i,km1,m) )
               end if

               !!! Provisional up and down draft values
               q_up(i,k) = q_int(i,k)
               q_dn(i,k) = q_int(i,k)

               !!! provisional tendencies
               q_tend(k) = 0._r8

            end do

            !!! Do levels adjacent to top and bottom
            k = 2
            km1 = 1
            kk = pver

            mupdudp = mu(i,kk) + du_tmp(i,kk)*dptmp(i,kk)
            if (mupdudp > mbsth) then
               q_up(i,kk) = ( +eu_tmp(i,kk) * frac_insol(i,kk) * q(i,kk,m) * dptmp(i,kk) )/mupdudp
            endif
            if (md(i,k) < -mbsth) then
               q_dn(i,k) =  ( -ed_tmp(i,km1) * frac_insol(i,km1) * q(i,km1,m) * dptmp(i,km1) )/md(i,k)
            endif

            !!! Updraft from bottom to top
            do kk = pver-1,1,-1
               kkp1 = min(pver,kk+1)
               mupdudp = mu(i,kk) + du_tmp(i,kk)*dptmp(i,kk)
               if (mupdudp > mbsth) then
                  q_up(i,kk) = (  mu(i,kkp1)*q_up(i,kkp1)   &
                                + eu_tmp(i,kk) * frac_insol(i,kk) * q(i,kk,m) * dptmp(i,kk) &
                               )/mupdudp
               endif
            end do

            !!! Downdraft from top to bottom
            do k = 3,pver
               km1 = max(1,k-1)
               if (md(i,k) < -mbsth) then
                  q_dn(i,k) =  (  md(i,km1)*q_dn(i,km1)-ed_tmp(i,km1)*frac_insol(i,km1)*q(i,km1,m) &
                                  *dptmp(i,km1) )/md(i,k)
               endif
            end do

            !!! WTF does this do?
            do k = ktm,pver
               km1 = max(1,k-1)
               kp1 = min(pver,k+1)

               ! limit fluxes outside convection to mass in appropriate layer
               ! these limiters are probably only safe for positive definite quantitities
               ! it assumes that mu and md already satisfy a courant number limit of 1

               fluxin =  mu(i,kp1)*q_up(i,kp1)+ mu(i,k)*min(q_int(i,k),q(i,km1,m)) &
                         -(md(i,k)  *q_dn(i,k) + md(i,kp1)*min(q_int(i,kp1),q(i,kp1,m)))

               fluxout = mu(i,k)*q_up(i,k) + mu(i,kp1)*min(q_int(i,kp1),q(i,k,m)) &
                         -(md(i,kp1)*q_dn(i,kp1) + md(i,k)*min(q_int(i,k),q(i,k,m)))

               netflux = fluxin - fluxout
               if ( abs(netflux) < ( max(fluxin,fluxout)*1.e-12_r8 ) ) then
                  netflux = 0._r8
               endif
               q_tend(k) = netflux/dptmp(i,k)

            end do

            !!! WTF does this do?
            do k = kbm,pver
               km1 = max(1,k-1)
               if (k == cld_bot_idx(i)) then

                  !!! version 3
                  fluxin =  mu(i,k)*min(q_int(i,k),q(i,km1,m)) - md(i,k)*q_dn(i,k)
                  fluxout = mu(i,k)*q_up(i,k) - md(i,k)*min(q_int(i,k),q(i,k,m))

                  netflux = fluxin - fluxout
                  if (abs(netflux) < max(fluxin,fluxout)*1.e-12_r8) then
                     netflux = 0._r8
                  endif
                  q_tend(k) = netflux/dptmp(i,k)
               else if (k > mx(i)) then
                  q_tend(k) = 0._r8
               end if
            end do

            !!! put provisional tendency in output array
            tend_out(i,:,m) = q_tend(:)
            
         end if ! i=1,ncol
      end if ! do_transport
   end do ! m=1,pcnst
   !-----------------------------------------------------------------------------------------
   !-----------------------------------------------------------------------------------------

   return

end subroutine crm_bulk_transport_tend


!======================================================================================================
!======================================================================================================
! subroutine crm_bulk_cloud_diag
   !-----------------------------------------------------------------------------------------
   ! Purpose:  ?????
   ! Walter Hannah (LLNL), 2018
   !-----------------------------------------------------------------------------------------

! end subroutine crm_bulk_cloud_diag
!======================================================================================================
!======================================================================================================

#endif /*CRM*/
#endif /*CBT*/

end module crm_bulk_transport_mod
