module crm_types

   use params, only: crm_rknd
   use crmdims, only: nclubbvars, crm_nx_rad, crm_ny_rad, crm_nx, crm_ny, crm_nz

#if defined( m2005 ) && defined( MODAL_AERO )
   use modal_aero_data, only: ntot_amode
#endif

   implicit none
   private

   public crm_state_type
   public crm_input_type
   public crm_output_type

   !------------------------------------------------------------------------------------------------
   type crm_state_type
      ! Purpose: Define a class that will encapulate the CRM-level data that needs to be passed
      ! between the GCM and the CRM, and to encapulate operations on that data (i.e.,
      ! initialization, writing fields to netCDF files, displaying info, et.c).

      ! TODO: should these be allocatable, and initialized once per node (at, maybe,
      ! crm_physics_init)?
      !---------------------------------------------------------------------------------------------
      ! CRM-scale fields

      ! NOTE: these were intent(inout) before, so these need to persist across calls; pointers so
      ! they can be used without making a bunch of temporary arrays. Dimensions should be
      ! (pcols,crm_nx,crm_ny,crm_nz)
      real(crm_rknd), pointer :: u_wind(:,:,:,:)       ! CRM u-wind component
      real(crm_rknd), pointer :: v_wind(:,:,:,:)       ! CRM v-wind component
      real(crm_rknd), pointer :: w_wind(:,:,:,:)       ! CRM w-wind component
      real(crm_rknd), pointer :: temperature(:,:,:,:)  ! CRM temperuture

      ! Microphysics
      ! NOTE: These are terrible variable names...replace with more descriptive names.
      ! for m2005...
      real(crm_rknd), pointer :: qt(:,:,:,:) 
      real(crm_rknd), pointer :: nc(:,:,:,:)
      real(crm_rknd), pointer :: qr(:,:,:,:)
      real(crm_rknd), pointer :: nr(:,:,:,:)
      real(crm_rknd), pointer :: qi(:,:,:,:)
      real(crm_rknd), pointer :: ni(:,:,:,:)
      real(crm_rknd), pointer :: qs(:,:,:,:)
      real(crm_rknd), pointer :: ns(:,:,:,:)
      real(crm_rknd), pointer :: qg(:,:,:,:)
      real(crm_rknd), pointer :: ng(:,:,:,:)
      real(crm_rknd), pointer :: qc(:,:,:,:)

      ! for sam1mom...
      real(crm_rknd), pointer :: qp(:,:,:,:)
      real(crm_rknd), pointer :: qn(:,:,:,:)

      ! Quantities used by the radiation code. Note that these are strange in that they are 
      ! time-averages, but spatially-resolved.
      ! TODO: can these be instantaneous fields from the end of the CRM run instead? Or would it be
      ! better to leave them as time averages, but apply an overlap assumption to the individual
      ! columns, since the clouds could be "smeared out" by the time averaging? This may be especially
      ! true when using a reduced grid for the radiation (crm_nx_rad < crm_nx), since in this case
      ! spatial averaging is done as well.
      !
      ! TODO: should these exist in crm_state or crm_output? Or something else entirely?
      real(crm_rknd), pointer :: t_rad  (:,:,:,:) ! rad temperuture
      real(crm_rknd), pointer :: qv_rad (:,:,:,:) ! rad vapor
      real(crm_rknd), pointer :: qc_rad (:,:,:,:) ! rad cloud water
      real(crm_rknd), pointer :: qi_rad (:,:,:,:) ! rad cloud ice
      real(crm_rknd), pointer :: cld_rad(:,:,:,:) ! rad cloud fraction
#ifdef m2005
      real(crm_rknd), pointer :: nc_rad (:,:,:,:) ! rad cloud droplet number (#/kg)
      real(crm_rknd), pointer :: ni_rad (:,:,:,:) ! rad cloud ice crystal number (#/kg)
      real(crm_rknd), pointer :: qs_rad (:,:,:,:) ! rad cloud snow (kg/kg)
      real(crm_rknd), pointer :: ns_rad (:,:,:,:) ! rad cloud snow crystal number (#/kg)
#endif

      ! These are copies of the SAM cloud and precip liquid and ice, previously
      ! passed in and out of crm() via qc_crm, qi_crm, etc. How are these
      ! different from the above microphysics variables? It looks like these are
      ! derived from the above, so maybe we don't need to pass these in and out?
      real(crm_rknd), allocatable :: qcl(:,:,:,:)
      real(crm_rknd), allocatable :: qci(:,:,:,:)
      real(crm_rknd), allocatable :: qpl(:,:,:,:)
      real(crm_rknd), allocatable :: qpi(:,:,:,:)

      real(crm_rknd), allocatable :: crm_tk (:,:,:,:)
      real(crm_rknd), allocatable :: crm_tkh(:,:,:,:)

      real(crm_rknd), allocatable :: prec_crm(:,:,:) ! CRM precipiation rate (surface)

      ! Stuff for 2-moment
      real(crm_rknd), allocatable :: wvar(:,:,:,:) ! vertical velocity variance (m/s)
      real(crm_rknd), allocatable :: aut (:,:,:,:) ! cloud water autoconversion (1/s)
      real(crm_rknd), allocatable :: acc (:,:,:,:) ! cloud water accretion (1/s)
      real(crm_rknd), allocatable :: evpc(:,:,:,:) ! cloud water evaporation (1/s)
      real(crm_rknd), allocatable :: evpr(:,:,:,:) ! rain evaporation (1/s)
      real(crm_rknd), allocatable :: mlt (:,:,:,:) ! ice, snow, graupel melting (1/s)
      real(crm_rknd), allocatable :: sub (:,:,:,:) ! ice, snow, graupel sublimation (1/s)
      real(crm_rknd), allocatable :: dep (:,:,:,:) ! ice, snow, graupel deposition (1/s)
      real(crm_rknd), allocatable :: con (:,:,:,:) ! cloud water condensation(1/s)
   contains
      ! Type-bound procedures. Initialization should nullify fields
      procedure, public :: initialize=>crm_state_initialize
      procedure, public :: finalize=>crm_state_finalize
      !procedure, public :: dump=>crm_state_dump

   end type crm_state_type
   !------------------------------------------------------------------------------------------------
   type crm_input_type

      real(crm_rknd), allocatable :: zmid(:,:)           ! Global grid height (m)
      real(crm_rknd), allocatable :: zint(:,:)           ! Global grid interface height (m)
      real(crm_rknd), allocatable :: tl(:,:)             ! Global grid temperature (K)
      real(crm_rknd), allocatable :: ql(:,:)             ! Global grid water vapor (g/g)
      real(crm_rknd), allocatable :: qccl(:,:)           ! Global grid cloud liquid water (g/g)
      real(crm_rknd), allocatable :: qiil(:,:)           ! Global grid cloud ice (g/g)
      real(crm_rknd), allocatable :: ps(:)               ! Global grid surface pressure (Pa)
      real(crm_rknd), allocatable :: pmid(:,:)           ! Global grid pressure (Pa)
      real(crm_rknd), allocatable :: pdel(:,:)           ! Layer's pressure thickness (Pa)
      real(crm_rknd), allocatable :: phis(:)             ! Global grid surface geopotential (m2/s2)

      real(crm_rknd), allocatable :: ul(:,:)             ! Global grid u (m/s)
      real(crm_rknd), allocatable :: vl(:,:)             ! Global grid v (m/s)

      real(crm_rknd), pointer     :: qrad(:,:,:,:)       ! CRM rad. heating

      real(crm_rknd), allocatable :: ocnfrac(:)          ! area fraction of the ocean
      real(crm_rknd), allocatable :: tau00  (:)          ! large-scale surface stress (N/m2)
      real(crm_rknd), allocatable :: wndls  (:)          ! large-scale surface wind (m/s)
      real(crm_rknd), allocatable :: bflxls (:)          ! large-scale surface buoyancy flux (K m/s)
      real(crm_rknd), allocatable :: fluxu00(:)          ! surface momenent fluxes [N/m2]
      real(crm_rknd), allocatable :: fluxv00(:)          ! surface momenent fluxes [N/m2]
      real(crm_rknd), allocatable :: fluxt00(:)          ! surface sensible heat fluxes [K Kg/ (m2 s)]
      real(crm_rknd), allocatable :: fluxq00(:)          ! surface latent heat fluxes [ kg/(m2 s)]

#if defined( m2005 ) && defined( MODAL_AERO )
      real(crm_rknd), allocatable :: naermod (:,:,:)     ! Aerosol number concentration [/m3]
      real(crm_rknd), allocatable :: vaerosol(:,:,:)     ! aerosol volume concentration [m3/m3]
      real(crm_rknd), allocatable :: hygro   (:,:,:)     ! hygroscopicity of aerosol mode 
#endif

#if defined(SP_ESMT)
      real(crm_rknd), allocatable :: ul_esmt(:,:)        ! input u for ESMT
      real(crm_rknd), allocatable :: vl_esmt(:,:)        ! input v for ESMT
#endif

   contains
      procedure, public :: initialize=>crm_input_initialize
      procedure, public :: finalize=>crm_input_finalize
   end type crm_input_type
   !------------------------------------------------------------------------------------------------
   type crm_output_type
      ! Derived type to encapsulate CRM output fields (things that are either
      ! time-averaged, have reduced spatial dimensions, or both)

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
#endif

#ifdef m2005
      ! Time and domain averaged process rates
      real(crm_rknd), allocatable :: aut_crm_a (:,:)  ! cloud water autoconversion (1/s)
      real(crm_rknd), allocatable :: acc_crm_a (:,:)  ! cloud water accretion (1/s)
      real(crm_rknd), allocatable :: evpc_crm_a(:,:)  ! cloud water evaporation (1/s)
      real(crm_rknd), allocatable :: evpr_crm_a(:,:)  ! rain evaporation (1/s)
      real(crm_rknd), allocatable :: mlt_crm_a (:,:)  ! ice, snow, graupel melting (1/s)
      real(crm_rknd), allocatable :: sub_crm_a (:,:)  ! ice, snow, graupel sublimation (1/s)
      real(crm_rknd), allocatable :: dep_crm_a (:,:)  ! ice, snow, graupel deposition (1/s)
      real(crm_rknd), allocatable :: con_crm_a (:,:)  ! cloud water condensation(1/s)
#endif /* m2005 */

      ! These are all time and spatial averages, on the GCM grid
      real(crm_rknd), allocatable :: cld   (:,:)  ! cloud fraction
      real(crm_rknd), allocatable :: gicewp(:,:)  ! ice water path
      real(crm_rknd), allocatable :: gliqwp(:,:)  ! ice water path
      real(crm_rknd), allocatable :: mctot (:,:)  ! cloud mass flux
      real(crm_rknd), allocatable :: mcup  (:,:)  ! updraft cloud mass flux
      real(crm_rknd), allocatable :: mcdn  (:,:)  ! downdraft cloud mass flux
      real(crm_rknd), allocatable :: mcuup (:,:)  ! unsat updraft cloud mass flux
      real(crm_rknd), allocatable :: mcudn (:,:)  ! unsat downdraft cloud mass flux

      ! For convective transport
      real(crm_rknd), allocatable :: mu_crm(:,:)  ! mass flux up
      real(crm_rknd), allocatable :: md_crm(:,:)  ! mass flux down
      real(crm_rknd), allocatable :: du_crm(:,:)  ! mass detrainment from updraft
      real(crm_rknd), allocatable :: eu_crm(:,:)  ! mass entrainment from updraft
      real(crm_rknd), allocatable :: ed_crm(:,:)  ! mass detrainment from downdraft
      real(crm_rknd), allocatable :: jt_crm(:)       ! index of cloud (convection) top
      real(crm_rknd), allocatable :: mx_crm(:)       ! index of cloud (convection) bottom

      ! Other stuff...
      real(crm_rknd), allocatable :: flux_qt      (:,:)       ! nonprecipitating water flux           [kg/m2/s]
      real(crm_rknd), allocatable :: fluxsgs_qt   (:,:)       ! sgs nonprecipitating water flux    [kg/m2/s]
      real(crm_rknd), allocatable :: tkez         (:,:)       ! tke profile               [kg/m/s2]
      real(crm_rknd), allocatable :: tkesgsz      (:,:)       ! sgs tke profile        [kg/m/s2]
      real(crm_rknd), allocatable :: tkz          (:,:)       ! tk profile                [m2/s]
      real(crm_rknd), allocatable :: flux_u       (:,:)       ! x-momentum flux          [m2/s2]
      real(crm_rknd), allocatable :: flux_v       (:,:)       ! y-momentum flux          [m2/s2]
      real(crm_rknd), allocatable :: flux_qp      (:,:)       ! precipitating water flux [kg/m2/s or mm/s]
      real(crm_rknd), allocatable :: precflux     (:,:)       ! precipitation flux      [m/s]
      real(crm_rknd), allocatable :: qt_ls        (:,:)       ! tendency of nonprec water due to large-scale  [kg/kg/s]
      real(crm_rknd), allocatable :: qt_trans     (:,:)       ! tendency of nonprec water due to transport  [kg/kg/s]
      real(crm_rknd), allocatable :: qp_trans     (:,:)       ! tendency of prec water due to transport [kg/kg/s]
      real(crm_rknd), allocatable :: qp_fall      (:,:)       ! tendency of prec water due to fall-out   [kg/kg/s]
      real(crm_rknd), allocatable :: qp_src       (:,:)       ! tendency of prec water due to conversion  [kg/kg/s]
      real(crm_rknd), allocatable :: qp_evp       (:,:)       ! tendency of prec water due to evp         [kg/kg/s]
      real(crm_rknd), allocatable :: t_ls         (:,:)       ! tendency of lwse  due to large-scale        [kg/kg/s] ???
      real(crm_rknd), allocatable :: prectend     (:)            ! column integrated tendency in precipitating water+ice (kg/m2/s)
      real(crm_rknd), allocatable :: precstend    (:)            ! column integrated tendency in precipitating ice (kg/m2/s)
      real(crm_rknd), allocatable :: taux_crm     (:)            ! zonal CRM surface stress perturbation (N/m2)
      real(crm_rknd), allocatable :: tauy_crm     (:)            ! merid CRM surface stress perturbation (N/m2)
      real(crm_rknd), allocatable :: z0m          (:)            ! surface stress (N/m2)
      real(crm_rknd), allocatable :: timing_factor(:)            ! crm cpu efficiency


   contains
      procedure, public :: initialize=>crm_output_initialize
      procedure, public :: finalize=>crm_output_finalize
   end type crm_output_type

contains

   !------------------------------------------------------------------------------------------------
   ! Type-bound procedures for crm_state_type
   subroutine crm_state_initialize(this, ncrms)
      class(crm_state_type), intent(inout) :: this
      integer, intent(in) :: ncrms

      ! Nullify pointers
      this%u_wind => null()
      this%v_wind => null()
      this%w_wind => null()
      this%temperature => null()

      this%qt => null()
      this%qc => null()
      this%qi => null()
      this%qr => null()
      this%qs => null()
      this%qg => null()
      this%nc => null()
      this%ni => null()
      this%nr => null()
      this%ns => null()
      this%ng => null()

      this%qp => null()
      this%qn => null()

      allocate(this%qcl(ncrms,crm_nx,crm_ny,crm_nz))
      allocate(this%qci(ncrms,crm_nx,crm_ny,crm_nz))
      allocate(this%qpl(ncrms,crm_nx,crm_ny,crm_nz))
      allocate(this%qpi(ncrms,crm_nx,crm_ny,crm_nz))

      if (.not. allocated(this%crm_tk )) allocate(this%crm_tk (ncrms,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(this%crm_tkh)) allocate(this%crm_tkh(ncrms,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(this%prec_crm)) allocate(this%prec_crm(ncrms,crm_nx,crm_ny))

      if (.not. allocated(this%wvar)) allocate(this%wvar(ncrms,crm_nx, crm_ny, crm_nz))
      if (.not. allocated(this%aut)) allocate(this%aut (ncrms,crm_nx, crm_ny, crm_nz))
      if (.not. allocated(this%acc)) allocate(this%acc (ncrms,crm_nx, crm_ny, crm_nz))
      if (.not. allocated(this%evpc)) allocate(this%evpc(ncrms,crm_nx, crm_ny, crm_nz))
      if (.not. allocated(this%evpr)) allocate(this%evpr(ncrms,crm_nx, crm_ny, crm_nz))
      if (.not. allocated(this%mlt)) allocate(this%mlt (ncrms,crm_nx, crm_ny, crm_nz))
      if (.not. allocated(this%sub)) allocate(this%sub (ncrms,crm_nx, crm_ny, crm_nz))
      if (.not. allocated(this%dep)) allocate(this%dep (ncrms,crm_nx, crm_ny, crm_nz))
      if (.not. allocated(this%con)) allocate(this%con (ncrms,crm_nx, crm_ny, crm_nz))

      ! Initialize 
      this%crm_tk = 0.0
      this%crm_tkh = 0.0
      this%prec_crm = 0.0
      
      this%wvar = 0.0
      this%aut  = 0.0
      this%acc  = 0.0
      this%evpc = 0.0
      this%evpr = 0.0
      this%mlt  = 0.0
      this%sub  = 0.0
      this%dep  = 0.0
      this%con  = 0.0

   end subroutine crm_state_initialize
   !------------------------------------------------------------------------------------------------
   subroutine crm_state_finalize(this)
      class(crm_state_type), intent(inout) :: this

      ! Nullify pointers
      this%u_wind => null()
      this%v_wind => null()
      this%w_wind => null()
      this%temperature => null()

#ifdef m2005
      this%qt => null()
      this%qc => null()
      this%qi => null()
      this%qr => null()
      this%qs => null()
      this%qg => null()
      this%nc => null()
      this%ni => null()
      this%nr => null()
      this%ns => null()
      this%ng => null()
#else
      this%qt => null()
      this%qp => null()
      this%qn => null()
#endif

      deallocate(this%qcl)
      deallocate(this%qci)
      deallocate(this%qpl)
      deallocate(this%qpi)
      deallocate(this%crm_tk )
      deallocate(this%crm_tkh)
      deallocate(this%prec_crm)
   end subroutine crm_state_finalize
   !------------------------------------------------------------------------------------------------
   ! Type-bound procedures for crm_state_type
   subroutine crm_input_initialize(this, ncrms, nlev)
      class(crm_input_type), intent(inout) :: this
      integer, intent(in) :: ncrms, nlev
      
      this%qrad => null()

      if (.not. allocated(this%zmid))     allocate(this%zmid(ncrms,nlev))
      if (.not. allocated(this%zint))     allocate(this%zint(ncrms,nlev+1))
      if (.not. allocated(this%tl))       allocate(this%tl(ncrms,nlev))
      if (.not. allocated(this%ql))       allocate(this%ql(ncrms,nlev))
      if (.not. allocated(this%qccl))     allocate(this%qccl(ncrms,nlev))
      if (.not. allocated(this%qiil))     allocate(this%qiil(ncrms,nlev))
      if (.not. allocated(this%ps))       allocate(this%ps(ncrms))
      if (.not. allocated(this%pmid))     allocate(this%pmid(ncrms,nlev))
      if (.not. allocated(this%pdel))     allocate(this%pdel(ncrms,nlev))
      if (.not. allocated(this%phis))     allocate(this%phis(ncrms))
      if (.not. allocated(this%ul))       allocate(this%ul(ncrms,nlev))
      if (.not. allocated(this%vl))       allocate(this%vl(ncrms,nlev))

      if (.not. allocated(this%ocnfrac))  allocate(this%ocnfrac(ncrms))
      if (.not. allocated(this%tau00))    allocate(this%tau00(ncrms))
      if (.not. allocated(this%wndls))    allocate(this%wndls(ncrms))
      if (.not. allocated(this%bflxls))   allocate(this%bflxls(ncrms))
      if (.not. allocated(this%fluxu00))  allocate(this%fluxu00(ncrms))
      if (.not. allocated(this%fluxv00))  allocate(this%fluxv00(ncrms))
      if (.not. allocated(this%fluxt00))  allocate(this%fluxt00(ncrms))
      if (.not. allocated(this%fluxq00))  allocate(this%fluxq00(ncrms))

#if defined( m2005 ) && defined( MODAL_AERO )
      if (.not. allocated(this%naermod))  allocate(this%naermod(ncrms,nlev,ntot_amode))
      if (.not. allocated(this%vaerosol)) allocate(this%vaerosol(ncrms,nlev,ntot_amode))
      if (.not. allocated(this%hygro))    allocate(this%hygro(ncrms,nlev,ntot_amode))
#endif

#if defined(SP_ESMT)
      if (.not. allocated(this%ul_esmt))  allocate(this%ul_esmt(ncrms,nlev))
      if (.not. allocated(this%vl_esmt))  allocate(this%vl_esmt(ncrms,nlev))
#endif

   end subroutine crm_input_initialize
   !------------------------------------------------------------------------------------------------
   subroutine crm_input_finalize(this)
      class(crm_input_type), intent(inout) :: this

      this%qrad => null()

      deallocate(this%zmid)
      deallocate(this%zint)
      deallocate(this%tl)
      deallocate(this%ql)
      deallocate(this%qccl)
      deallocate(this%qiil)
      deallocate(this%ps)
      deallocate(this%pmid)
      deallocate(this%pint)
      deallocate(this%pdel)
      deallocate(this%phis)
      deallocate(this%ul)
      deallocate(this%vl)

      deallocate(this%ocnfrac)
      deallocate(this%tau00)
      deallocate(this%wndls)
      deallocate(this%bflxls)
      deallocate(this%fluxu00)
      deallocate(this%fluxv00)
      deallocate(this%fluxt00)
      deallocate(this%fluxq00)

#if defined( m2005 ) && defined( MODAL_AERO )
      deallocate(this%naermod)
      deallocate(this%vaerosol)
      deallocate(this%hygro)
#endif

#if defined(SP_ESMT)
      deallocate(this%ul_esmt)
      deallocate(this%vl_esmt)
#endif

   end subroutine crm_input_finalize 
   !------------------------------------------------------------------------------------------------
   ! Type-bound procedures for crm_output_type
   subroutine crm_output_initialize(this, ncrms, nlev)
      class(crm_output_type), intent(inout) :: this
      integer, intent(in), optional :: ncrms, nlev

      ! Allocate arrays if dimensions are passed as input
      if (present(ncrms)) then
         ! Allocate (time-averaged?) fields
         
         ! Allocate domain and time-averaged fields
         if (.not. allocated(this%cltot)) allocate(this%cltot(ncrms))
         if (.not. allocated(this%cllow)) allocate(this%cllow(ncrms))
         if (.not. allocated(this%clmed)) allocate(this%clmed(ncrms))
         if (.not. allocated(this%clhgh)) allocate(this%clhgh(ncrms))

         if (.not. allocated(this%precc)) allocate(this%precc(ncrms))
         if (.not. allocated(this%precl)) allocate(this%precl(ncrms))
         if (.not. allocated(this%precsc)) allocate(this%precsc(ncrms))
         if (.not. allocated(this%precsl)) allocate(this%precsl(ncrms))

         if (.not. allocated(this%cldtop)) allocate(this%cldtop(ncrms,crm_nz))

         if (.not. allocated(this%cldtop)) allocate(this%cldtop(ncrms,crm_nz))

         if (.not. allocated(this%qc_mean)) allocate(this%qc_mean(ncrms,nlev))
         if (.not. allocated(this%qi_mean)) allocate(this%qi_mean(ncrms,nlev))
         if (.not. allocated(this%qs_mean)) allocate(this%qs_mean(ncrms,nlev))
         if (.not. allocated(this%qg_mean)) allocate(this%qg_mean(ncrms,nlev))
         if (.not. allocated(this%qr_mean)) allocate(this%qr_mean(ncrms,nlev))
#ifdef m2005
         if (.not. allocated(this%nc_mean)) allocate(this%nc_mean(ncrms,nlev))
         if (.not. allocated(this%ni_mean)) allocate(this%ni_mean(ncrms,nlev))
         if (.not. allocated(this%ns_mean)) allocate(this%ns_mean(ncrms,nlev))
         if (.not. allocated(this%ng_mean)) allocate(this%ng_mean(ncrms,nlev))
         if (.not. allocated(this%nr_mean)) allocate(this%nr_mean(ncrms,nlev))
#endif

#ifdef m2005
         if (.not. allocated(this%aut_crm_a )) allocate(this%aut_crm_a (ncrms,nlev))
         if (.not. allocated(this%acc_crm_a )) allocate(this%acc_crm_a (ncrms,nlev))
         if (.not. allocated(this%evpc_crm_a)) allocate(this%evpc_crm_a(ncrms,nlev))
         if (.not. allocated(this%evpr_crm_a)) allocate(this%evpr_crm_a(ncrms,nlev))
         if (.not. allocated(this%mlt_crm_a )) allocate(this%mlt_crm_a (ncrms,nlev))
         if (.not. allocated(this%sub_crm_a )) allocate(this%sub_crm_a (ncrms,nlev))
         if (.not. allocated(this%dep_crm_a )) allocate(this%dep_crm_a (ncrms,nlev))
         if (.not. allocated(this%con_crm_a )) allocate(this%con_crm_a (ncrms,nlev))
#endif /* m2005 */

         if (.not. allocated(this%cld   )) allocate(this%cld   (ncrms,nlev))  ! cloud fraction
         if (.not. allocated(this%gicewp)) allocate(this%gicewp(ncrms,nlev))  ! ice water path
         if (.not. allocated(this%gliqwp)) allocate(this%gliqwp(ncrms,nlev))  ! ice water path
         if (.not. allocated(this%mctot )) allocate(this%mctot (ncrms,nlev))  ! cloud mass flux
         if (.not. allocated(this%mcup  )) allocate(this%mcup  (ncrms,nlev))  ! updraft cloud mass flux
         if (.not. allocated(this%mcdn  )) allocate(this%mcdn  (ncrms,nlev))  ! downdraft cloud mass flux
         if (.not. allocated(this%mcuup )) allocate(this%mcuup (ncrms,nlev))  ! unsat updraft cloud mass flux
         if (.not. allocated(this%mcudn )) allocate(this%mcudn (ncrms,nlev))  ! unsat downdraft cloud mass flux

         if (.not. allocated(this%mu_crm)) allocate(this%mu_crm(ncrms,nlev))  ! mass flux up
         if (.not. allocated(this%md_crm)) allocate(this%md_crm(ncrms,nlev))  ! mass flux down
         if (.not. allocated(this%du_crm)) allocate(this%du_crm(ncrms,nlev))  ! mass detrainment from updraft
         if (.not. allocated(this%eu_crm)) allocate(this%eu_crm(ncrms,nlev))  ! mass entrainment from updraft
         if (.not. allocated(this%ed_crm)) allocate(this%ed_crm(ncrms,nlev))  ! mass detrainment from downdraft
         if (.not. allocated(this%jt_crm)) allocate(this%jt_crm(ncrms))       ! index of cloud (convection) top
         if (.not. allocated(this%mx_crm)) allocate(this%mx_crm(ncrms))       ! index of cloud (convection) bottom

         if (.not. allocated(this%flux_qt      )) allocate(this%flux_qt             (ncrms,nlev))
         if (.not. allocated(this%fluxsgs_qt   )) allocate(this%fluxsgs_qt          (ncrms,nlev))
         if (.not. allocated(this%tkez         )) allocate(this%tkez                (ncrms,nlev))
         if (.not. allocated(this%tkesgsz      )) allocate(this%tkesgsz             (ncrms,nlev))
         if (.not. allocated(this%tkz          )) allocate(this%tkz                 (ncrms,nlev))
         if (.not. allocated(this%flux_u       )) allocate(this%flux_u              (ncrms,nlev))
         if (.not. allocated(this%flux_v       )) allocate(this%flux_v              (ncrms,nlev))
         if (.not. allocated(this%flux_qp      )) allocate(this%flux_qp             (ncrms,nlev))
         if (.not. allocated(this%precflux     )) allocate(this%precflux            (ncrms,nlev))
         if (.not. allocated(this%qt_ls        )) allocate(this%qt_ls               (ncrms,nlev))
         if (.not. allocated(this%qt_trans     )) allocate(this%qt_trans            (ncrms,nlev))
         if (.not. allocated(this%qp_trans     )) allocate(this%qp_trans            (ncrms,nlev))
         if (.not. allocated(this%qp_fall      )) allocate(this%qp_fall             (ncrms,nlev))
         if (.not. allocated(this%qp_src       )) allocate(this%qp_src              (ncrms,nlev))
         if (.not. allocated(this%qp_evp       )) allocate(this%qp_evp              (ncrms,nlev))
         if (.not. allocated(this%t_ls         )) allocate(this%t_ls                (ncrms,nlev))
         if (.not. allocated(this%prectend     )) allocate(this%prectend            (ncrms))
         if (.not. allocated(this%precstend    )) allocate(this%precstend           (ncrms))
         if (.not. allocated(this%taux_crm     )) allocate(this%taux_crm            (ncrms))
         if (.not. allocated(this%tauy_crm     )) allocate(this%tauy_crm            (ncrms))
         if (.not. allocated(this%z0m          )) allocate(this%z0m                 (ncrms))
         if (.not. allocated(this%timing_factor)) allocate(this%timing_factor       (ncrms))


      end if

      this%cltot = 0.0
      this%cllow = 0.0
      this%clmed = 0.0
      this%clhgh = 0.0

      this%cldtop = 0.0
      this%precc = 0.0
      this%precl = 0.0
      this%precsc = 0.0
      this%precsl = 0.0

      this%qc_mean = 0.0
      this%qi_mean = 0.0
      this%qs_mean = 0.0
      this%qg_mean = 0.0
      this%qr_mean = 0.0
#ifdef m2005
      this%nc_mean = 0.0
      this%ni_mean = 0.0
      this%ns_mean = 0.0
      this%ng_mean = 0.0
      this%nr_mean = 0.0
#endif

#ifdef m2005
      this%aut_crm_a = 0.0
      this%acc_crm_a = 0.0
      this%evpc_crm_a = 0.0
      this%evpr_crm_a = 0.0
      this%mlt_crm_a = 0.0
      this%sub_crm_a = 0.0
      this%dep_crm_a = 0.0
      this%con_crm_a = 0.0
#endif

      this%cld    = 0.
      this%gicewp = 0
      this%gliqwp = 0
      this%mctot  = 0.
      this%mcup   = 0.
      this%mcdn   = 0.
      this%mcuup  = 0.
      this%mcudn  = 0.

      ! Convective transport
      this%mu_crm = 0.
      this%md_crm = 0.
      this%eu_crm = 0.
      this%du_crm = 0.
      this%ed_crm = 0.
      this%jt_crm = 0.
      this%mx_crm = 0.

   end subroutine crm_output_initialize
   !------------------------------------------------------------------------------------------------
   subroutine crm_output_finalize(this)
      class(crm_output_type), intent(inout) :: this

      deallocate(this%cltot)
      deallocate(this%cllow)
      deallocate(this%clmed)
      deallocate(this%clhgh)
      deallocate(this%cldtop)
      deallocate(this%precc)
      deallocate(this%precl)
      deallocate(this%precsc)
      deallocate(this%precsl)
   end subroutine crm_output_finalize
   !------------------------------------------------------------------------------------------------

end module crm_types
