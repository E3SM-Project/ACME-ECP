module crm_output_module
   use params,       only: crm_rknd
   use crmdims,      only: crm_nx, crm_ny, crm_nz
   implicit none
   public crm_output_type
   type crm_output_type
      ! Derived type to encapsulate CRM output fields (things that are output
      ! only, intent(out) in the original implementation)

      ! These are copies of the SAM cloud and precip liquid and ice, previously
      ! passed in and out of crm() via qc_crm, qi_crm, etc.
      real(crm_rknd), allocatable :: qcl(:,:,:,:)
      real(crm_rknd), allocatable :: qci(:,:,:,:)
      real(crm_rknd), allocatable :: qpl(:,:,:,:)
      real(crm_rknd), allocatable :: qpi(:,:,:,:)

      real(crm_rknd), allocatable :: tk (:,:,:,:)
      real(crm_rknd), allocatable :: tkh(:,:,:,:)
      real(crm_rknd), allocatable :: prec_crm(:,:,:) ! CRM precipiation rate (surface)

      ! 2-moment process rates
      real(crm_rknd), allocatable :: wvar(:,:,:,:) ! vertical velocity variance (m/s)
      real(crm_rknd), allocatable :: aut (:,:,:,:) ! cloud water autoconversion (1/s)
      real(crm_rknd), allocatable :: acc (:,:,:,:) ! cloud water accretion (1/s)
      real(crm_rknd), allocatable :: evpc(:,:,:,:) ! cloud water evaporation (1/s)
      real(crm_rknd), allocatable :: evpr(:,:,:,:) ! rain evaporation (1/s)
      real(crm_rknd), allocatable :: mlt (:,:,:,:) ! ice, snow, graupel melting (1/s)
      real(crm_rknd), allocatable :: sub (:,:,:,:) ! ice, snow, graupel sublimation (1/s)
      real(crm_rknd), allocatable :: dep (:,:,:,:) ! ice, snow, graupel deposition (1/s)
      real(crm_rknd), allocatable :: con (:,:,:,:) ! cloud water condensation(1/s)

      ! Cloud area fractions
      real(crm_rknd), allocatable :: cltot(:)  ! shaded cloud fraction
      real(crm_rknd), allocatable :: clhgh(:)  ! shaded cloud fraction
      real(crm_rknd), allocatable :: clmed(:)  ! shaded cloud fraction
      real(crm_rknd), allocatable :: cllow(:)  ! shaded cloud fraction

      real(crm_rknd), allocatable :: cldtop(:,:)  ! cloud top ... pressure???
      real(crm_rknd), allocatable :: precc(:)   ! convective precipitation rate
      real(crm_rknd), allocatable :: precl(:)   ! stratiform precipitation rate
      real(crm_rknd), allocatable :: precsc(:)   ! convective snow precipitation rate
      real(crm_rknd), allocatable :: precsl(:)   ! stratiform snow precipitation rate

      ! TODO: These diagnostics are currently all on the GCM vertical grid. I think this is
      ! misleading though, and overly complicates crm_module. I think the better thing to do would
      ! be to define everything within crm_module on the CRM grid, and then do the
      ! mapping/interpolation at the GCM (crm_physics) level. For now though, I am just copying
      ! these over directly to minimize chances of me making a mistake. Also, some of these probably
      ! do not need to be calculated here, and might make more sense to calculate at the
      ! crm_physics_tend level. For example, I think tendencies should be calculated in
      ! crm_physics_tend, from, for example, something like crm_output%uwind - crm_input%uwind.
      real(crm_rknd), allocatable :: ta_mean(:,:)  ! mean cloud water
      real(crm_rknd), allocatable :: qv_mean(:,:)  ! mean cloud water
      real(crm_rknd), allocatable :: qc_mean(:,:)  ! mean cloud water
      real(crm_rknd), allocatable :: qi_mean(:,:)  ! mean cloud ice
      real(crm_rknd), allocatable :: qs_mean(:,:)  ! mean snow
      real(crm_rknd), allocatable :: qg_mean(:,:)  ! mean graupel
      real(crm_rknd), allocatable :: qr_mean(:,:)  ! mean rain
#ifdef m2005
      real(crm_rknd), allocatable :: nc_mean(:,:)  ! mean cloud water  (#/kg)
      real(crm_rknd), allocatable :: ni_mean(:,:)  ! mean cloud ice    (#/kg)
      real(crm_rknd), allocatable :: ns_mean(:,:)  ! mean snow         (#/kg)
      real(crm_rknd), allocatable :: ng_mean(:,:)  ! mean graupel      (#/kg)
      real(crm_rknd), allocatable :: nr_mean(:,:)  ! mean rain         (#/kg)

      ! Time and domain averaged process rates
      real(crm_rknd), allocatable :: aut_a (:,:)  ! cloud water autoconversion (1/s)
      real(crm_rknd), allocatable :: acc_a (:,:)  ! cloud water accretion (1/s)
      real(crm_rknd), allocatable :: evpc_a(:,:)  ! cloud water evaporation (1/s)
      real(crm_rknd), allocatable :: evpr_a(:,:)  ! rain evaporation (1/s)
      real(crm_rknd), allocatable :: mlt_a (:,:)  ! ice, snow, graupel melting (1/s)
      real(crm_rknd), allocatable :: sub_a (:,:)  ! ice, snow, graupel sublimation (1/s)
      real(crm_rknd), allocatable :: dep_a (:,:)  ! ice, snow, graupel deposition (1/s)
      real(crm_rknd), allocatable :: con_a (:,:)  ! cloud water condensation(1/s)
#endif /* m2005 */

#if defined( SPMOMTRANS )
      real(crm_rknd), allocatable :: ultend(:,:)            ! tendency of ul
      real(crm_rknd), allocatable :: vltend(:,:)            ! tendency of vl
#endif

#if defined( SP_ESMT )
      real(crm_rknd), allocatable :: u_tend_esmt(:,:)       ! CRM scalar u-momentum tendency
      real(crm_rknd), allocatable :: v_tend_esmt(:,:)       ! CRM scalar v-momentum tendency
#endif

      real(crm_rknd), allocatable :: sltend  (:,:)          ! CRM output tendency of static energy
      real(crm_rknd), allocatable :: qltend  (:,:)          ! CRM output tendency of water vapor
      real(crm_rknd), allocatable :: qcltend (:,:)          ! CRM output tendency of cloud liquid water
      real(crm_rknd), allocatable :: qiltend (:,:)          ! CRM output tendency of cloud ice

      ! These are all time and spatial averages, on the GCM grid
      real(crm_rknd), allocatable :: cld   (:,:)      ! cloud fraction
      real(crm_rknd), allocatable :: gicewp(:,:)      ! ice water path
      real(crm_rknd), allocatable :: gliqwp(:,:)      ! ice water path
      real(crm_rknd), allocatable :: mctot (:,:)      ! cloud mass flux
      real(crm_rknd), allocatable :: mcup  (:,:)      ! updraft cloud mass flux
      real(crm_rknd), allocatable :: mcdn  (:,:)      ! downdraft cloud mass flux
      real(crm_rknd), allocatable :: mcuup (:,:)      ! unsat updraft cloud mass flux
      real(crm_rknd), allocatable :: mcudn (:,:)      ! unsat downdraft cloud mass flux

      ! For convective transport
      real(crm_rknd), allocatable :: mu_crm(:,:)      ! mass flux up
      real(crm_rknd), allocatable :: md_crm(:,:)      ! mass flux down
      real(crm_rknd), allocatable :: du_crm(:,:)      ! mass detrainment from updraft
      real(crm_rknd), allocatable :: eu_crm(:,:)      ! mass entrainment from updraft
      real(crm_rknd), allocatable :: ed_crm(:,:)      ! mass detrainment from downdraft
      real(crm_rknd), allocatable :: jt_crm(:)        ! index of cloud (convection) top
      real(crm_rknd), allocatable :: mx_crm(:)        ! index of cloud (convection) bottom

      ! Other stuff...
      real(crm_rknd), allocatable :: flux_qt      (:,:)  ! nonprecip water flux        [kg/m2/s]
      real(crm_rknd), allocatable :: fluxsgs_qt   (:,:)  ! sgs non-precip water flux   [kg/m2/s]
      real(crm_rknd), allocatable :: tkez         (:,:)  ! tke profile                 [kg/m/s2]
      real(crm_rknd), allocatable :: tkesgsz      (:,:)  ! sgs tke profile             [kg/m/s2]
      real(crm_rknd), allocatable :: tkz          (:,:)  ! tk profile                  [m2/s]
      real(crm_rknd), allocatable :: flux_u       (:,:)  ! x-momentum flux             [m2/s2]
      real(crm_rknd), allocatable :: flux_v       (:,:)  ! y-momentum flux             [m2/s2]
      real(crm_rknd), allocatable :: flux_qp      (:,:)  ! precipitating water flux    [kg/m2/s or mm/s]
      real(crm_rknd), allocatable :: precflux     (:,:)  ! precipitation flux          [m/s]
      real(crm_rknd), allocatable :: qt_ls        (:,:)  ! tend of nonprec water due to large-scale   [kg/kg/s]
      real(crm_rknd), allocatable :: qt_trans     (:,:)  ! tend of nonprec water due to transport     [kg/kg/s]
      real(crm_rknd), allocatable :: qp_trans     (:,:)  ! tend of    prec water due to transport     [kg/kg/s]
      real(crm_rknd), allocatable :: qp_fall      (:,:)  ! tend of    prec water due to fall-out      [kg/kg/s]
      real(crm_rknd), allocatable :: qp_src       (:,:)  ! tend of    prec water due to conversion    [kg/kg/s]
      real(crm_rknd), allocatable :: qp_evp       (:,:)  ! tend of    prec water due to evp           [kg/kg/s]
      real(crm_rknd), allocatable :: t_ls         (:,:)  ! tend of lwse  due to large-scale           [kg/kg/s] ???
      real(crm_rknd), allocatable :: prectend     (:)    ! column integrated tend in precip water+ice [kg/m2/s]
      real(crm_rknd), allocatable :: precstend    (:)    ! column integrated tend in precip ice       [kg/m2/s]
      real(crm_rknd), allocatable :: taux     (:)    ! zonal CRM surface stress perturbation      [N/m2]
      real(crm_rknd), allocatable :: tauy     (:)    ! merid CRM surface stress perturbation      [N/m2]
      real(crm_rknd), allocatable :: z0m          (:)    ! surface stress                             [N/m2]
      real(crm_rknd), allocatable :: timing_factor(:)    ! crm cpu efficiency

   contains
      procedure, public :: initialize=>crm_output_initialize
      procedure, public :: finalize=>crm_output_finalize
   end type crm_output_type

contains

   !------------------------------------------------------------------------------------------------
   ! Type-bound procedures for crm_output_type
   subroutine crm_output_initialize(this, ncol, nlev)
      class(crm_output_type), intent(inout) :: this
      integer, intent(in), optional :: ncol, nlev

      ! Allocate arrays if dimensions are passed as input
      if (present(ncol)) then

         ! Allocate instantaneous outputs
         if (.not. allocated(this%qcl)) allocate(this%qcl(ncol,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(this%qci)) allocate(this%qci(ncol,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(this%qpl)) allocate(this%qpl(ncol,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(this%qpi)) allocate(this%qpi(ncol,crm_nx,crm_ny,crm_nz))

         if (.not. allocated(this%tk )) allocate(this%tk (ncol,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(this%tkh)) allocate(this%tkh(ncol,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(this%prec_crm)) allocate(this%prec_crm(ncol,crm_nx,crm_ny))

         if (.not. allocated(this%wvar)) allocate(this%wvar(ncol,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(this%aut))  allocate(this%aut (ncol,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(this%acc))  allocate(this%acc (ncol,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(this%evpc)) allocate(this%evpc(ncol,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(this%evpr)) allocate(this%evpr(ncol,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(this%mlt))  allocate(this%mlt (ncol,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(this%sub))  allocate(this%sub (ncol,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(this%dep))  allocate(this%dep (ncol,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(this%con))  allocate(this%con (ncol,crm_nx,crm_ny,crm_nz))


         ! Allocate domain and time-averaged fields
         if (.not. allocated(this%cltot)) allocate(this%cltot(ncol))
         if (.not. allocated(this%cllow)) allocate(this%cllow(ncol))
         if (.not. allocated(this%clmed)) allocate(this%clmed(ncol))
         if (.not. allocated(this%clhgh)) allocate(this%clhgh(ncol))

         if (.not. allocated(this%precc))  allocate(this%precc(ncol))
         if (.not. allocated(this%precl))  allocate(this%precl(ncol))
         if (.not. allocated(this%precsc)) allocate(this%precsc(ncol))
         if (.not. allocated(this%precsl)) allocate(this%precsl(ncol))

         ! NOTE: this output had a bug in the previous implementation
         if (.not. allocated(this%cldtop)) allocate(this%cldtop(ncol,nlev))

         if (.not. allocated(this%ta_mean)) allocate(this%ta_mean(ncol,nlev))
         if (.not. allocated(this%qv_mean)) allocate(this%qv_mean(ncol,nlev))
         if (.not. allocated(this%qc_mean)) allocate(this%qc_mean(ncol,nlev))
         if (.not. allocated(this%qi_mean)) allocate(this%qi_mean(ncol,nlev))
         if (.not. allocated(this%qs_mean)) allocate(this%qs_mean(ncol,nlev))
         if (.not. allocated(this%qg_mean)) allocate(this%qg_mean(ncol,nlev))
         if (.not. allocated(this%qr_mean)) allocate(this%qr_mean(ncol,nlev))
#ifdef m2005
         if (.not. allocated(this%nc_mean)) allocate(this%nc_mean(ncol,nlev))
         if (.not. allocated(this%ni_mean)) allocate(this%ni_mean(ncol,nlev))
         if (.not. allocated(this%ns_mean)) allocate(this%ns_mean(ncol,nlev))
         if (.not. allocated(this%ng_mean)) allocate(this%ng_mean(ncol,nlev))
         if (.not. allocated(this%nr_mean)) allocate(this%nr_mean(ncol,nlev))

         if (.not. allocated(this%aut_a )) allocate(this%aut_a (ncol,nlev))
         if (.not. allocated(this%acc_a )) allocate(this%acc_a (ncol,nlev))
         if (.not. allocated(this%evpc_a)) allocate(this%evpc_a(ncol,nlev))
         if (.not. allocated(this%evpr_a)) allocate(this%evpr_a(ncol,nlev))
         if (.not. allocated(this%mlt_a )) allocate(this%mlt_a (ncol,nlev))
         if (.not. allocated(this%sub_a )) allocate(this%sub_a (ncol,nlev))
         if (.not. allocated(this%dep_a )) allocate(this%dep_a (ncol,nlev))
         if (.not. allocated(this%con_a )) allocate(this%con_a (ncol,nlev))
#endif /* m2005 */

#if defined( SPMOMTRANS )
         if (.not. allocated(this%ultend )) allocate(this%ultend (ncol,nlev))
         if (.not. allocated(this%vltend )) allocate(this%vltend (ncol,nlev))
#endif

#if defined( SP_ESMT )
         if (.not. allocated(this%u_tend_esmt )) allocate(this%u_tend_esmt (ncol,nlev))
         if (.not. allocated(this%v_tend_esmt )) allocate(this%v_tend_esmt (ncol,nlev))
#endif
         
         if (.not. allocated(this%sltend ))  allocate(this%sltend (ncol,nlev))
         if (.not. allocated(this%qltend ))  allocate(this%qltend (ncol,nlev))
         if (.not. allocated(this%qcltend))  allocate(this%qcltend(ncol,nlev))
         if (.not. allocated(this%qiltend))  allocate(this%qiltend(ncol,nlev))

         if (.not. allocated(this%cld   )) allocate(this%cld   (ncol,nlev))  ! cloud fraction
         if (.not. allocated(this%gicewp)) allocate(this%gicewp(ncol,nlev))  ! ice water path
         if (.not. allocated(this%gliqwp)) allocate(this%gliqwp(ncol,nlev))  ! ice water path
         if (.not. allocated(this%mctot )) allocate(this%mctot (ncol,nlev))  ! cloud mass flux
         if (.not. allocated(this%mcup  )) allocate(this%mcup  (ncol,nlev))  ! updraft cloud mass flux
         if (.not. allocated(this%mcdn  )) allocate(this%mcdn  (ncol,nlev))  ! downdraft cloud mass flux
         if (.not. allocated(this%mcuup )) allocate(this%mcuup (ncol,nlev))  ! unsat updraft cloud mass flux
         if (.not. allocated(this%mcudn )) allocate(this%mcudn (ncol,nlev))  ! unsat downdraft cloud mass flux

         if (.not. allocated(this%mu_crm)) allocate(this%mu_crm(ncol,nlev))  ! mass flux up
         if (.not. allocated(this%md_crm)) allocate(this%md_crm(ncol,nlev))  ! mass flux down
         if (.not. allocated(this%du_crm)) allocate(this%du_crm(ncol,nlev))  ! mass detrainment from updraft
         if (.not. allocated(this%eu_crm)) allocate(this%eu_crm(ncol,nlev))  ! mass entrainment from updraft
         if (.not. allocated(this%ed_crm)) allocate(this%ed_crm(ncol,nlev))  ! mass detrainment from downdraft
         if (.not. allocated(this%jt_crm)) allocate(this%jt_crm(ncol))       ! index of cloud (convection) top
         if (.not. allocated(this%mx_crm)) allocate(this%mx_crm(ncol))       ! index of cloud (convection) bottom

         if (.not. allocated(this%flux_qt      )) allocate(this%flux_qt      (ncol,nlev))
         if (.not. allocated(this%fluxsgs_qt   )) allocate(this%fluxsgs_qt   (ncol,nlev))
         if (.not. allocated(this%tkez         )) allocate(this%tkez         (ncol,nlev))
         if (.not. allocated(this%tkesgsz      )) allocate(this%tkesgsz      (ncol,nlev))
         if (.not. allocated(this%tkz          )) allocate(this%tkz          (ncol,nlev))
         if (.not. allocated(this%flux_u       )) allocate(this%flux_u       (ncol,nlev))
         if (.not. allocated(this%flux_v       )) allocate(this%flux_v       (ncol,nlev))
         if (.not. allocated(this%flux_qp      )) allocate(this%flux_qp      (ncol,nlev))
         if (.not. allocated(this%precflux     )) allocate(this%precflux     (ncol,nlev))
         if (.not. allocated(this%qt_ls        )) allocate(this%qt_ls        (ncol,nlev))
         if (.not. allocated(this%qt_trans     )) allocate(this%qt_trans     (ncol,nlev))
         if (.not. allocated(this%qp_trans     )) allocate(this%qp_trans     (ncol,nlev))
         if (.not. allocated(this%qp_fall      )) allocate(this%qp_fall      (ncol,nlev))
         if (.not. allocated(this%qp_src       )) allocate(this%qp_src       (ncol,nlev))
         if (.not. allocated(this%qp_evp       )) allocate(this%qp_evp       (ncol,nlev))
         if (.not. allocated(this%t_ls         )) allocate(this%t_ls         (ncol,nlev))
         if (.not. allocated(this%prectend     )) allocate(this%prectend     (ncol))
         if (.not. allocated(this%precstend    )) allocate(this%precstend    (ncol))
         if (.not. allocated(this%taux         )) allocate(this%taux         (ncol))
         if (.not. allocated(this%tauy         )) allocate(this%tauy         (ncol))
         if (.not. allocated(this%z0m          )) allocate(this%z0m          (ncol))
         if (.not. allocated(this%timing_factor)) allocate(this%timing_factor(ncol))

      end if ! present(ncol)

      ! Initialize 
      this%qcl = 0
      this%qci = 0
      this%qpl = 0
      this%qpi = 0

      this%tk = 0
      this%tkh = 0
      this%prec_crm = 0

      ! 2-moment process rates
      this%wvar = 0
      this%aut  = 0 
      this%acc  = 0
      this%evpc = 0
      this%evpr = 0
      this%mlt  = 0
      this%sub  = 0
      this%dep  = 0
      this%con  = 0

      this%cltot = 0
      this%cllow = 0
      this%clmed = 0
      this%clhgh = 0

      this%cldtop = 0
      this%precc = 0
      this%precl = 0
      this%precsc = 0
      this%precsl = 0

      this%ta_mean = 0
      this%qv_mean = 0
      this%qc_mean = 0
      this%qi_mean = 0
      this%qs_mean = 0
      this%qg_mean = 0
      this%qr_mean = 0
#ifdef m2005
      this%nc_mean = 0
      this%ni_mean = 0
      this%ns_mean = 0
      this%ng_mean = 0
      this%nr_mean = 0

      this%aut_a = 0
      this%acc_a = 0
      this%evpc_a = 0
      this%evpr_a = 0
      this%mlt_a = 0
      this%sub_a = 0
      this%dep_a = 0
      this%con_a = 0
#endif

#if defined( SPMOMTRANS )
      this%ultend = 0
      this%vltend = 0
#endif

#if defined( SP_ESMT )
      this%u_tend_esmt = 0
      this%v_tend_esmt = 0
#endif

      this%sltend  = 0
      this%qltend  = 0
      this%qcltend = 0
      this%qiltend = 0

      this%cld    = 0
      this%gicewp = 0
      this%gliqwp = 0
      this%mctot  = 0
      this%mcup   = 0
      this%mcdn   = 0
      this%mcuup  = 0
      this%mcudn  = 0

      ! Convective transport
      this%mu_crm = 0
      this%md_crm = 0
      this%eu_crm = 0
      this%du_crm = 0
      this%ed_crm = 0
      this%jt_crm = 0
      this%mx_crm = 0

      ! Other stuff...
      this%flux_qt       = 0
      this%fluxsgs_qt    = 0
      this%tkez          = 0
      this%tkesgsz       = 0
      this%tkz           = 0
      this%flux_u        = 0
      this%flux_v        = 0
      this%flux_qp       = 0
      this%precflux      = 0
      this%qt_ls         = 0
      this%qt_trans      = 0
      this%qp_trans      = 0
      this%qp_fall       = 0
      this%qp_src        = 0
      this%qp_evp        = 0
      this%t_ls          = 0
      this%prectend      = 0
      this%precstend     = 0
      this%taux      = 0
      this%tauy      = 0
      this%z0m           = 0
      this%timing_factor = 0

   end subroutine crm_output_initialize
   !------------------------------------------------------------------------------------------------
   subroutine crm_output_finalize(this)
      class(crm_output_type), intent(inout) :: this
      if (allocated(this%qcl)) deallocate(this%qcl)
      if (allocated(this%qci)) deallocate(this%qci)
      if (allocated(this%qpl)) deallocate(this%qpl)
      if (allocated(this%qpi)) deallocate(this%qpi)
      if (allocated(this%tk )) deallocate(this%tk )
      if (allocated(this%tkh)) deallocate(this%tkh)
      if (allocated(this%prec_crm)) deallocate(this%prec_crm)

      if (allocated(this%wvar)) deallocate(this%wvar)
      if (allocated(this%aut)) deallocate(this%aut)
      if (allocated(this%acc)) deallocate(this%acc)
      if (allocated(this%evpc)) deallocate(this%evpc)
      if (allocated(this%evpr)) deallocate(this%evpr)
      if (allocated(this%mlt)) deallocate(this%mlt)
      if (allocated(this%sub)) deallocate(this%sub)
      if (allocated(this%dep)) deallocate(this%dep)
      if (allocated(this%con)) deallocate(this%con)

      if (allocated(this%cltot)) deallocate(this%cltot)
      if (allocated(this%cllow)) deallocate(this%cllow)
      if (allocated(this%clmed)) deallocate(this%clmed)
      if (allocated(this%clhgh)) deallocate(this%clhgh)
      if (allocated(this%cldtop)) deallocate(this%cldtop)
      if (allocated(this%precc)) deallocate(this%precc)
      if (allocated(this%precl)) deallocate(this%precl)
      if (allocated(this%precsc)) deallocate(this%precsc)
      if (allocated(this%precsl)) deallocate(this%precsl)

      if (allocated(this%ta_mean)) deallocate(this%ta_mean)
      if (allocated(this%qv_mean)) deallocate(this%qv_mean)
      if (allocated(this%qc_mean)) deallocate(this%qc_mean)
      if (allocated(this%qi_mean)) deallocate(this%qi_mean)
      if (allocated(this%qs_mean)) deallocate(this%qs_mean)
      if (allocated(this%qg_mean)) deallocate(this%qg_mean)
      if (allocated(this%qr_mean)) deallocate(this%qr_mean)
#ifdef m2005
      if (allocated(this%nc_mean)) deallocate(this%nc_mean)
      if (allocated(this%ni_mean)) deallocate(this%ni_mean)
      if (allocated(this%ns_mean)) deallocate(this%ns_mean)
      if (allocated(this%ng_mean)) deallocate(this%ng_mean)
      if (allocated(this%nr_mean)) deallocate(this%nr_mean)

      ! Time and domain-averaged process rates
      if (allocated(this%aut_a)) deallocate(this%aut_a)
      if (allocated(this%acc_a)) deallocate(this%acc_a)
      if (allocated(this%evpc_a)) deallocate(this%evpc_a)
      if (allocated(this%evpr_a)) deallocate(this%evpr_a)
      if (allocated(this%mlt_a)) deallocate(this%mlt_a)
      if (allocated(this%sub_a)) deallocate(this%sub_a)
      if (allocated(this%dep_a)) deallocate(this%dep_a)
      if (allocated(this%con_a)) deallocate(this%con_a)
#endif

#if defined( SPMOMTRANS )
      if (allocated(this%ultend)) deallocate(this%ultend)
      if (allocated(this%vltend)) deallocate(this%vltend)
#endif

#if defined( SP_ESMT )
      if (allocated(this%u_tend_esmt)) deallocate(this%u_tend_esmt)
      if (allocated(this%v_tend_esmt)) deallocate(this%v_tend_esmt)
#endif

      if (allocated(this%sltend)) deallocate(this%sltend)
      if (allocated(this%qltend)) deallocate(this%qltend)
      if (allocated(this%qcltend)) deallocate(this%qcltend)
      if (allocated(this%qiltend)) deallocate(this%qiltend)

      if (allocated(this%cld)) deallocate(this%cld)
      if (allocated(this%gicewp)) deallocate(this%gicewp)
      if (allocated(this%gliqwp)) deallocate(this%gliqwp)
      if (allocated(this%mctot)) deallocate(this%mctot)
      if (allocated(this%mcup)) deallocate(this%mcup)
      if (allocated(this%mcdn)) deallocate(this%mcdn)
      if (allocated(this%mcuup)) deallocate(this%mcuup)
      if (allocated(this%mcudn)) deallocate(this%mcudn)

      if (allocated(this%mu_crm)) deallocate(this%mu_crm)
      if (allocated(this%md_crm)) deallocate(this%md_crm)
      if (allocated(this%du_crm)) deallocate(this%du_crm)
      if (allocated(this%eu_crm)) deallocate(this%eu_crm)
      if (allocated(this%ed_crm)) deallocate(this%ed_crm)
      if (allocated(this%jt_crm)) deallocate(this%jt_crm)
      if (allocated(this%mx_crm)) deallocate(this%mx_crm)

      if (allocated(this%flux_qt)) deallocate(this%flux_qt)
      if (allocated(this%fluxsgs_qt)) deallocate(this%fluxsgs_qt)
      if (allocated(this%tkez)) deallocate(this%tkez)
      if (allocated(this%tkesgsz)) deallocate(this%tkesgsz)
      if (allocated(this%tkz)) deallocate(this%tkz)
      if (allocated(this%flux_u)) deallocate(this%flux_u)
      if (allocated(this%flux_v)) deallocate(this%flux_v)
      if (allocated(this%flux_qp)) deallocate(this%flux_qp)
      if (allocated(this%precflux)) deallocate(this%precflux)
      if (allocated(this%qt_ls)) deallocate(this%qt_ls)
      if (allocated(this%qt_trans)) deallocate(this%qt_trans)
      if (allocated(this%qp_trans)) deallocate(this%qp_trans)
      if (allocated(this%qp_fall)) deallocate(this%qp_fall)
      if (allocated(this%qp_src)) deallocate(this%qp_src)
      if (allocated(this%qp_evp)) deallocate(this%qp_evp)
      if (allocated(this%t_ls)) deallocate(this%t_ls)
      if (allocated(this%prectend)) deallocate(this%prectend)
      if (allocated(this%precstend)) deallocate(this%precstend)
      if (allocated(this%taux)) deallocate(this%taux)
      if (allocated(this%tauy)) deallocate(this%tauy)
      if (allocated(this%z0m)) deallocate(this%z0m)
      if (allocated(this%timing_factor)) deallocate(this%timing_factor)

   end subroutine crm_output_finalize
   !------------------------------------------------------------------------------------------------

end module crm_output_module
