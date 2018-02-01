module vars
  use grid
  use params, only: crm_rknd
#ifdef CRM
#ifdef MODAL_AERO
  use modal_aero_data,   only: ntot_amode
#endif
#endif

  implicit none
  !--------------------------------------------------------------------
  ! prognostic variables:

  real(crm_rknd), allocatable :: u  (:,:,:) ! x-wind
  real(crm_rknd), allocatable :: v  (:,:,:) ! y-wind
  real(crm_rknd), allocatable :: w  (:,:,:) ! z-wind
  real(crm_rknd), allocatable :: t  (:,:,:) ! liquid/ice water static energy

  !--------------------------------------------------------------------
  ! diagnostic variables:

  real(crm_rknd), allocatable :: p   (:,:,:)         ! perturbation pressure (from Poison eq)
  real(crm_rknd), allocatable :: tabs(:,:,:)         ! temperature
  real(crm_rknd), allocatable :: qv  (:,:,:)        ! water vapor
  real(crm_rknd), allocatable :: qcl (:,:,:)        ! liquid water  (condensate)
  real(crm_rknd), allocatable :: qpl (:,:,:)        ! liquid water  (precipitation)
  real(crm_rknd), allocatable :: qci (:,:,:)        ! ice water  (condensate)
  real(crm_rknd), allocatable :: qpi (:,:,:)        ! ice water  (precipitation)

  real(crm_rknd), allocatable :: tke2(:,:,:)       ! SGS TKE
  real(crm_rknd), allocatable :: tk2 (:,:,:)     ! SGS eddyviscosity

  !--------------------------------------------------------------------
  ! time-tendencies for prognostic variables

  real(crm_rknd), allocatable :: dudt(:,:,:,:)
  real(crm_rknd), allocatable :: dvdt(:,:,:,:)
  real(crm_rknd), allocatable :: dwdt(:,:,:,:)

  !----------------------------------------------------------------
  ! Temporary storage array:

  real(crm_rknd), allocatable :: misc(:,:,:)
  !------------------------------------------------------------------
  ! fluxes at the top and bottom of the domain:

  real(crm_rknd), allocatable :: fluxbu  (:,:)
  real(crm_rknd), allocatable :: fluxbv  (:,:)
  real(crm_rknd), allocatable :: fluxbt  (:,:)
  real(crm_rknd), allocatable :: fluxbq  (:,:)
  real(crm_rknd), allocatable :: fluxtu  (:,:)
  real(crm_rknd), allocatable :: fluxtv  (:,:)
  real(crm_rknd), allocatable :: fluxtt  (:,:)
  real(crm_rknd), allocatable :: fluxtq  (:,:)
  real(crm_rknd), allocatable :: fzero   (:,:)
  real(crm_rknd), allocatable :: precsfc (:,:) ! surface precip. rate
  real(crm_rknd), allocatable :: precssfc(:,:) ! surface ice precip. rate

  !-----------------------------------------------------------------
  ! profiles

  real(crm_rknd), allocatable :: t0   (:)
  real(crm_rknd), allocatable :: q0   (:)
  real(crm_rknd), allocatable :: qv0  (:)
  real(crm_rknd), allocatable :: tabs0(:)
  real(crm_rknd), allocatable :: tl0  (:)
  real(crm_rknd), allocatable :: tv0  (:)
  real(crm_rknd), allocatable :: u0   (:)
  real(crm_rknd), allocatable :: v0   (:)
  real(crm_rknd), allocatable :: tg0  (:)
  real(crm_rknd), allocatable :: qg0  (:)
  real(crm_rknd), allocatable :: ug0  (:)
  real(crm_rknd), allocatable :: vg0  (:)
  real(crm_rknd), allocatable :: p0   (:)
  real(crm_rknd), allocatable :: tke0 (:)
  real(crm_rknd), allocatable :: t01  (:)
  real(crm_rknd), allocatable :: q01  (:)
  real(crm_rknd), allocatable :: qp0  (:)
  real(crm_rknd), allocatable :: qn0  (:)

  !-----------------------------------------------------------------
  ! reference vertical profiles:

  real(crm_rknd), allocatable :: prespot(:)  ! (1000./pres)**R/cp
  real(crm_rknd), allocatable :: rho    (:)  ! air density at pressure levels,kg/m3
  real(crm_rknd), allocatable :: rhow   (:)  ! air density at vertical velocity levels,kg/m3
  real(crm_rknd), allocatable :: bet    (:)  ! = ggr/tv0
  real(crm_rknd), allocatable :: gamaz  (:)  ! ggr/cp*z
  real(crm_rknd), allocatable :: wsub   (:)  ! Large-scale subsidence velocity,m/s
  real(crm_rknd), allocatable :: qtend  (:)  ! Large-scale tendency for total water
  real(crm_rknd), allocatable :: ttend  (:)  ! Large-scale tendency for temp.
  real(crm_rknd), allocatable :: utend  (:)  ! Large-scale tendency for u
  real(crm_rknd), allocatable :: vtend  (:)  ! Large-scale tendency for v

  !---------------------------------------------------------------------
  !  Horizontally varying stuff (as a function of xy)
  !
  real(crm_rknd), allocatable :: sstxy    (:,:)  !  surface temperature xy-distribution
  real(crm_rknd), allocatable :: fcory    (:)    !  Coriolis parameter xy-distribution
  real(crm_rknd), allocatable :: fcorzy   (:)    !  z-Coriolis parameter xy-distribution
  real(crm_rknd), allocatable :: latitude (:,:)  ! latitude (degrees)
  real(crm_rknd), allocatable :: longitude(:,:)  ! longitude(degrees)
  real(crm_rknd), allocatable :: prec_xy  (:,:)  ! mean precip. rate for outout
  real(crm_rknd), allocatable :: pw_xy    (:,:)  ! precipitable water
  real(crm_rknd), allocatable :: cw_xy    (:,:)  ! cloud water path
  real(crm_rknd), allocatable :: iw_xy    (:,:)  ! ice water path
  real(crm_rknd), allocatable :: cld_xy   (:,:)  ! cloud frequency
  real(crm_rknd), allocatable :: u200_xy  (:,:)  ! u-wind at 200 mb
  real(crm_rknd), allocatable :: usfc_xy  (:,:)  ! u-wind at at the surface
  real(crm_rknd), allocatable :: v200_xy  (:,:)  ! v-wind at 200 mb
  real(crm_rknd), allocatable :: vsfc_xy  (:,:)  ! v-wind at the surface
  real(crm_rknd), allocatable :: w500_xy  (:,:)  ! w at 500 mb

  !----------------------------------------------------------------------
  !  Vertical profiles of quantities sampled for statitistics purposes:

  real(crm_rknd), allocatable :: twle    (:)
  real(crm_rknd), allocatable :: twsb    (:)
  real(crm_rknd), allocatable :: precflux(:)
  real(crm_rknd), allocatable :: uwle    (:)
  real(crm_rknd), allocatable :: uwsb    (:)
  real(crm_rknd), allocatable :: vwle    (:)
  real(crm_rknd), allocatable :: vwsb    (:)
  real(crm_rknd)              :: w_max   
  real(crm_rknd)              :: u_max   
  real(crm_rknd), allocatable :: tkelediss(:)
  real(crm_rknd), allocatable :: t2leadv  (:)
  real(crm_rknd), allocatable :: t2legrad (:)
  real(crm_rknd), allocatable :: t2lediff (:)
  real(crm_rknd), allocatable :: t2lediss (:)
  real(crm_rknd), allocatable :: twleadv  (:)
  real(crm_rknd), allocatable :: twlediff (:)
  real(crm_rknd), allocatable :: tadv     (:)
  real(crm_rknd), allocatable :: tdiff    (:)
  real(crm_rknd), allocatable :: tlat     (:)
  real(crm_rknd), allocatable :: tlatqi   (:)
  real(crm_rknd), allocatable :: qifall   (:)
  real(crm_rknd), allocatable :: qpfall   (:)


  ! register functions:


  !real(crm_rknd), external :: esatw_crm,esati_crm,dtesatw_crm,dtesati_crm
  !real(crm_rknd), external :: qsatw_crm,qsati_crm,dtqsatw_crm,dtqsati_crm
  !integer, external :: lenstr, bytes_in_rec

  ! energy conservation diagnostics:

  real(8) :: total_water_before 
  real(8) :: total_water_after  
  real(8) :: total_water_evap   
  real(8) :: total_water_prec   
  real(8) :: total_water_ls     
  real(8) :: total_water_clubb  
  real(8) :: total_energy_before
  real(8) :: total_energy_after 
  real(8) :: total_energy_evap  
  real(8) :: total_energy_prec  
  real(8) :: total_energy_ls    
  real(8) :: total_energy_clubb 
  real(8) :: total_energy_rad   
  real(8), allocatable :: qtotmicro(:)  ! total water for water conservation test in microphysics +++mhwang


  !===========================================================================
  ! UW ADDITIONS

  ! conditional average statistics, subsumes cloud_factor, core_factor, coredn_factor
  real(crm_rknd), allocatable :: CF3D(:,:,:)  ! Cloud fraction
  ! =1.0 when there is no fractional cloudiness scheme
  ! = cloud fraction produced by fractioal cloudiness scheme when avaiable

  ! 850 mbar horizontal winds
  real(crm_rknd), allocatable :: u850_xy(:,:) ! zonal velocity at 850 mb
  real(crm_rknd), allocatable :: v850_xy(:,:) ! meridional velocity at 850 mb

  ! Surface pressure
  real(crm_rknd), allocatable :: psfc_xy(:,:) ! pressure (in millibar) at lowest grid point

  ! Saturated water vapor path, useful for computing column relative humidity
  real(crm_rknd), allocatable :: swvp_xy(:,:)  ! saturated water vapor path (wrt water)

  ! Cloud and echo top heights, and cloud top temperature (instantaneous)
  real(crm_rknd), allocatable :: cloudtopheight(:,:)
  real(crm_rknd), allocatable :: echotopheight (:,:)
  real(crm_rknd), allocatable :: cloudtoptemp  (:,:)

  ! END UW ADDITIONS
#if (defined CRM && defined MODAL_AERO)
  real(crm_rknd), allocatable :: naer (:,:)    ! Aerosol number concentration [/m3]
  real(crm_rknd), allocatable :: vaer (:,:)    ! aerosol volume concentration [m3/m3]
  real(crm_rknd), allocatable :: hgaer(:,:)    ! hygroscopicity of aerosol mode
#endif

  integer :: ncondavg, icondavg_cld, icondavg_cor, icondavg_cordn, &
  icondavg_satdn, icondavg_satup, icondavg_env
  real(crm_rknd), allocatable :: condavg_factor(:,:) ! replaces cloud_factor, core_factor
  real(crm_rknd), allocatable :: condavg_mask(:,:,:,:) ! indicator array for various conditional averages
  character(LEN=8), allocatable :: condavgname(:) ! array of short names
  character(LEN=25), allocatable :: condavglongname(:) ! array of long names

contains

  subroutine allocate_vars(ncrms)
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) :: zero
    allocate( u (dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm) ) ! x-wind
    allocate( v (dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm) ) ! y-wind
    allocate( w (dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz ) ) ! z-wind
    allocate( t (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ) ! liquid/ice water static energy
    allocate( p   (0:nx, (1-YES3D):ny, nzm) )     ! perturbation pressure (from Poison eq)
    allocate( tabs(nx, ny, nzm) )                 ! temperature
    allocate( qv  (nx, ny, nzm) )                ! water vapor
    allocate( qcl (nx, ny, nzm) )                ! liquid water  (condensate)
    allocate( qpl (nx, ny, nzm) )                ! liquid water  (precipitation)
    allocate( qci (nx, ny, nzm) )                ! ice water  (condensate)
    allocate( qpi (nx, ny, nzm) )                ! ice water  (precipitation)
    allocate( tke2(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) )   ! SGS TKE
    allocate( tk2 (0:nxp1, (1-YES3D):nyp1, nzm) ) ! SGS eddyviscosity
    allocate( dudt   (nxp1, ny, nzm, 3) )
    allocate( dvdt   (nx, nyp1, nzm, 3) )
    allocate( dwdt   (nx, ny, nz,  3) )
    allocate( misc   (nx, ny, nz) )
    allocate( fluxbu (nx, ny) )
    allocate( fluxbv (nx, ny) )
    allocate( fluxbt (nx, ny) )
    allocate( fluxbq (nx, ny) )
    allocate( fluxtu (nx, ny) )
    allocate( fluxtv (nx, ny) )
    allocate( fluxtt (nx, ny) )
    allocate( fluxtq (nx, ny) )
    allocate( fzero  (nx, ny) )
    allocate( precsfc(nx,ny) ) ! surface precip. rate
    allocate( precssfc(nx,ny) ) ! surface ice precip. rate
    allocate( t0   (nzm) )
    allocate( q0   (nzm) )
    allocate( qv0  (nzm) )
    allocate( tabs0(nzm) )
    allocate( tl0  (nzm) )
    allocate( tv0  (nzm) )
    allocate( u0   (nzm) )
    allocate( v0   (nzm) )
    allocate( tg0  (nzm) )
    allocate( qg0  (nzm) )
    allocate( ug0  (nzm) )
    allocate( vg0  (nzm) )
    allocate( p0   (nzm) )
    allocate( tke0 (nzm) )
    allocate( t01  (nzm) )
    allocate( q01  (nzm) )
    allocate( qp0  (nzm) )
    allocate( qn0  (nzm) )
    allocate( prespot  (nzm) ) ! (1000./pres)**R/cp
    allocate( rho      (nzm) ) ! air density at pressure levels,kg/m3
    allocate( rhow     (nz ) ) ! air density at vertical velocity levels,kg/m3
    allocate( bet      (nzm) ) ! = ggr/tv0
    allocate( gamaz    (nzm) ) ! ggr/cp*z
    allocate( wsub     (nz ) ) ! Large-scale subsidence velocity,m/s
    allocate( qtend    (nzm) ) ! Large-scale tendency for total water
    allocate( ttend    (nzm) ) ! Large-scale tendency for temp.
    allocate( utend    (nzm) ) ! Large-scale tendency for u
    allocate( vtend    (nzm) ) ! Large-scale tendency for v
    allocate( sstxy    (0:nx,(1-YES3D):ny) ) !  surface temperature xy-distribution
    allocate( fcory    (0:ny) )             !  Coriolis parameter xy-distribution
    allocate( fcorzy   (ny) )               !  z-Coriolis parameter xy-distribution
    allocate( latitude (nx,ny) )            ! latitude (degrees)
    allocate( longitude(nx,ny) )            ! longitude(degrees)
    allocate( prec_xy  (nx,ny) )            ! mean precip. rate for outout
    allocate( pw_xy    (nx,ny) )            ! precipitable water
    allocate( cw_xy    (nx,ny) )            ! cloud water path
    allocate( iw_xy    (nx,ny) )            ! ice water path
    allocate( cld_xy   (nx,ny) )            ! cloud frequency
    allocate( u200_xy  (nx,ny) )            ! u-wind at 200 mb
    allocate( usfc_xy  (nx,ny) )            ! u-wind at at the surface
    allocate( v200_xy  (nx,ny) )            ! v-wind at 200 mb
    allocate( vsfc_xy  (nx,ny) )            ! v-wind at the surface
    allocate( w500_xy  (nx,ny) )            ! w at 500 mb
    allocate( twle     (nz) )
    allocate( twsb     (nz) )
    allocate( precflux (nz) )
    allocate( uwle     (nz) )
    allocate( uwsb     (nz) )
    allocate( vwle     (nz) )
    allocate( vwsb     (nz) )
    allocate( tkelediss(nz) )
    allocate( t2leadv  (nz) )
    allocate( t2legrad (nz) )
    allocate( t2lediff (nz) )
    allocate( t2lediss (nz) )
    allocate( twleadv  (nz) )
    allocate( twlediff (nz) )
    allocate( tadv     (nz) )
    allocate( tdiff    (nz) )
    allocate( tlat     (nz) )
    allocate( tlatqi   (nz) )
    allocate( qifall   (nz) )
    allocate( qpfall   (nz) )
    allocate( qtotmicro(5) )
    allocate( CF3D          (1:nx, 1:ny, 1:nzm) )
    allocate( u850_xy       (nx,ny) )
    allocate( v850_xy       (nx,ny) )
    allocate( psfc_xy       (nx,ny) )
    allocate( swvp_xy       (nx,ny) )
    allocate( cloudtopheight(nx,ny) )
    allocate( echotopheight (nx,ny) )
    allocate( cloudtoptemp  (nx,ny) )
#if (defined CRM && defined MODAL_AERO)
    allocate( naer          (nzm, ntot_amode) )
    allocate( vaer          (nzm, ntot_amode) )
    allocate( hgaer         (nzm, ntot_amode) )
#endif

    zero = 0.

    u              = zero
    v              = zero
    w              = zero
    t              = zero
    p              = zero
    tabs           = zero
    qv             = zero
    qcl            = zero
    qpl            = zero
    qci            = zero
    qpi            = zero
    tke2           = zero
    tk2            = zero
    dudt           = zero
    dvdt           = zero
    dwdt           = zero
    misc           = zero
    fluxbu         = zero
    fluxbv         = zero
    fluxbt         = zero
    fluxbq         = zero
    fluxtu         = zero
    fluxtv         = zero
    fluxtt         = zero
    fluxtq         = zero
    fzero          = zero
    precsfc        = zero
    precssfc       = zero
    t0             = zero
    q0             = zero
    qv0            = zero
    tabs0          = zero
    tl0            = zero
    tv0            = zero
    u0             = zero
    v0             = zero
    tg0            = zero
    qg0            = zero
    ug0            = zero
    vg0            = zero
    p0             = zero
    tke0           = zero
    t01            = zero
    q01            = zero
    qp0            = zero
    qn0            = zero
    prespot        = zero
    rho            = zero
    rhow           = zero
    bet            = zero
    gamaz          = zero
    wsub           = zero
    qtend          = zero
    ttend          = zero
    utend          = zero
    vtend          = zero
    sstxy          = zero
    fcory          = zero
    fcorzy         = zero
    latitude       = zero
    longitude      = zero
    prec_xy        = zero
    pw_xy          = zero
    cw_xy          = zero
    iw_xy          = zero
    cld_xy         = zero
    u200_xy        = zero
    usfc_xy        = zero
    v200_xy        = zero
    vsfc_xy        = zero
    w500_xy        = zero
    twle           = zero
    twsb           = zero
    precflux       = zero
    uwle           = zero
    uwsb           = zero
    vwle           = zero
    vwsb           = zero
    tkelediss      = zero
    t2leadv        = zero
    t2legrad       = zero
    t2lediff       = zero
    t2lediss       = zero
    twleadv        = zero
    twlediff       = zero
    tadv           = zero
    tdiff          = zero
    tlat           = zero
    tlatqi         = zero
    qifall         = zero
    qpfall         = zero
    qtotmicro      = zero
    CF3D           = zero
    u850_xy        = zero
    v850_xy        = zero
    psfc_xy        = zero
    swvp_xy        = zero
    cloudtopheight = zero
    echotopheight  = zero
    cloudtoptemp   = zero
#if (defined CRM && defined MODAL_AERO)
    naer           = zero
    vaer           = zero
    hgaer          = zero
#endif
  end subroutine allocate_vars


  subroutine deallocate_vars
    implicit none
    deallocate( u ) ! x-wind
    deallocate( v ) ! y-wind
    deallocate( w ) ! z-wind
    deallocate( t ) ! liquid/ice water static energy
    deallocate( p    )   ! perturbation pressure (from Poison eq)
    deallocate( tabs )   ! temperature
    deallocate( qv   )   ! water vapor
    deallocate( qcl  )   ! liquid water  (condensate)
    deallocate( qpl  )   ! liquid water  (precipitation)
    deallocate( qci  )   ! ice water  (condensate)
    deallocate( qpi  )   ! ice water  (precipitation)
    deallocate( tke2 )   ! SGS TKE
    deallocate( tk2  )   ! SGS eddyviscosity
    deallocate( dudt    )
    deallocate( dvdt    )
    deallocate( dwdt    )
    deallocate( misc    )
    deallocate( fluxbu  )
    deallocate( fluxbv  )
    deallocate( fluxbt  )
    deallocate( fluxbq  )
    deallocate( fluxtu  )
    deallocate( fluxtv  )
    deallocate( fluxtt  )
    deallocate( fluxtq  )
    deallocate( fzero   )
    deallocate( precsfc ) ! surface precip. rate
    deallocate( precssfc ) ! surface ice precip. rate
    deallocate( t0    )
    deallocate( q0    )
    deallocate( qv0   )
    deallocate( tabs0 )
    deallocate( tl0   )
    deallocate( tv0   )
    deallocate( u0    )
    deallocate( v0    )
    deallocate( tg0   )
    deallocate( qg0   )
    deallocate( ug0   )
    deallocate( vg0   )
    deallocate( p0    )
    deallocate( tke0  )
    deallocate( t01   )
    deallocate( q01   )
    deallocate( qp0   )
    deallocate( qn0   )
    deallocate( prespot   ) ! (1000./pres)**R/cp
    deallocate( rho       ) ! air density at pressure levels,kg/m3
    deallocate( rhow      ) ! air density at vertical velocity levels,kg/m3
    deallocate( bet       ) ! = ggr/tv0
    deallocate( gamaz     ) ! ggr/cp*z
    deallocate( wsub      ) ! Large-scale subsidence velocity,m/s
    deallocate( qtend     ) ! Large-scale tendency for total water
    deallocate( ttend     ) ! Large-scale tendency for temp.
    deallocate( utend     ) ! Large-scale tendency for u
    deallocate( vtend     ) ! Large-scale tendency for v
    deallocate( sstxy     ) !  surface temperature xy-distribution
    deallocate( fcory     )            !  Coriolis parameter xy-distribution
    deallocate( fcorzy    )            !  z-Coriolis parameter xy-distribution
    deallocate( latitude  )            ! latitude (degrees)
    deallocate( longitude )            ! longitude(degrees)
    deallocate( prec_xy   )            ! mean precip. rate for outout
    deallocate( pw_xy     )            ! precipitable water
    deallocate( cw_xy     )            ! cloud water path
    deallocate( iw_xy     )            ! ice water path
    deallocate( cld_xy    )            ! cloud frequency
    deallocate( u200_xy   )            ! u-wind at 200 mb
    deallocate( usfc_xy   )            ! u-wind at at the surface
    deallocate( v200_xy   )            ! v-wind at 200 mb
    deallocate( vsfc_xy   )            ! v-wind at the surface
    deallocate( w500_xy   )            ! w at 500 mb
    deallocate( twle      )
    deallocate( twsb      )
    deallocate( precflux  )
    deallocate( uwle      )
    deallocate( uwsb      )
    deallocate( vwle      )
    deallocate( vwsb      )
    deallocate( tkelediss )
    deallocate( t2leadv   )
    deallocate( t2legrad  )
    deallocate( t2lediff  )
    deallocate( t2lediss  )
    deallocate( twleadv   )
    deallocate( twlediff  )
    deallocate( tadv      )
    deallocate( tdiff     )
    deallocate( tlat      )
    deallocate( tlatqi    )
    deallocate( qifall    )
    deallocate( qpfall    )
    deallocate( qtotmicro )
    deallocate( CF3D           )
    deallocate( u850_xy        )
    deallocate( v850_xy        )
    deallocate( psfc_xy        )
    deallocate( swvp_xy        )
    deallocate( cloudtopheight )
    deallocate( echotopheight  )
    deallocate( cloudtoptemp   )
#if (defined CRM && defined MODAL_AERO)
    deallocate( naer           )
    deallocate( vaer           )
    deallocate( hgaer          )
#endif
  end subroutine deallocate_vars

end module vars
