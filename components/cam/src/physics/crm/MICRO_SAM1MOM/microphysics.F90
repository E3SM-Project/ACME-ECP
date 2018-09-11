module microphysics
  use cloud_mod
  use precip_init_mod
  use precip_proc_mod

  ! module for original SAM bulk microphysics
  ! Marat Khairoutdinov, 2006

  use grid, only: nx,ny,nzm,nz, dimx1_s,dimx2_s,dimy1_s,dimy2_s ! subdomain grid information
  use params, only: doprecip, docloud, doclubb, crm_rknd
  use micro_params
  implicit none

  !----------------------------------------------------------------------
  !!! required definitions:

  integer, parameter :: nmicro_fields = 2   ! total number of prognostic water vars

  !!! microphysics prognostic variables are storred in this array:


  integer, parameter :: flag_wmass(nmicro_fields) = (/1,1/)
  integer, parameter :: index_water_vapor = 1 ! index for variable that has water vapor
  integer, parameter :: index_cloud_ice = 1   ! index for cloud ice (sedimentation)
  integer, parameter :: flag_precip(nmicro_fields) = (/0,1/)

  ! both variables correspond to mass, not number
  integer, parameter :: flag_number(nmicro_fields) = (/0,0/)

  ! SAM1MOM 3D microphysical fields are output by default.
  integer, parameter :: flag_micro3Dout(nmicro_fields) = (/0,0/)


  !!! these arrays are needed for output statistics:


  !======================================================================
  ! UW ADDITIONS

  !bloss: arrays with names/units for microphysical outputs in statistics.

  ! END UW ADDITIONS
  !======================================================================

  !------------------------------------------------------------------
  ! Optional (internal) definitions)

  ! make aliases for prognostic variables:
  ! note that the aliases should be local to microphysics

  real(crm_rknd) vrain, vsnow, vgrau, crain, csnow, cgrau  ! precomputed coefs for precip terminal velocity

  real(crm_rknd), allocatable, target :: micro_field(:,:,:,:,:)
  real(crm_rknd), allocatable :: fluxbmk (:,:,:,:) ! surface flux of tracers
  real(crm_rknd), allocatable :: fluxtmk (:,:,:,:) ! top boundary flux of tracers
  real(crm_rknd), allocatable :: mkwle  (:,:,:)  ! resolved vertical flux
  real(crm_rknd), allocatable :: mkwsb  (:,:,:)  ! SGS vertical flux
  real(crm_rknd), allocatable :: mkadv  (:,:,:)  ! tendency due to vertical advection
  real(crm_rknd), allocatable :: mklsadv(:,:,:)  ! tendency due to large-scale vertical advection
  real(crm_rknd), allocatable :: mkdiff (:,:,:)  ! tendency due to vertical diffusion
  real(crm_rknd), allocatable :: mstor  (:,:,:)  ! storage terms of microphysical variables
  character*3   , allocatable :: mkname       (:)
  character*80  , allocatable :: mklongname   (:)
  character*10  , allocatable :: mkunits      (:)
  real(crm_rknd), allocatable :: mkoutputscale(:)
  real(crm_rknd), allocatable :: qn(:,:,:,:)  ! cloud condensate (liquid + ice)
  real(crm_rknd), allocatable :: qpsrc(:,:)  ! source of precipitation microphysical processes
  real(crm_rknd), allocatable :: qpevp(:,:)  ! sink of precipitating water due to evaporation
  real(crm_rknd), pointer :: q (:,:,:,:)   ! total nonprecipitating water
  real(crm_rknd), pointer :: qp(:,:,:,:)  ! total precipitating water


CONTAINS


  subroutine allocate_micro(ncrms)
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) :: zero
    allocate( micro_field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, nmicro_fields,ncrms))
    allocate( fluxbmk (nx, ny, 1:nmicro_fields,ncrms) )
    allocate( fluxtmk (nx, ny, 1:nmicro_fields,ncrms) )
    allocate( mkwle  (nz,1:nmicro_fields,ncrms)  )
    allocate( mkwsb  (nz,1:nmicro_fields,ncrms)  )
    allocate( mkadv  (nz,1:nmicro_fields,ncrms)  )
    allocate( mklsadv(nz,1:nmicro_fields,ncrms)  )
    allocate( mkdiff (nz,1:nmicro_fields,ncrms)  )
    allocate( mstor  (nz,1:nmicro_fields,ncrms)  )
    allocate( mkname       (nmicro_fields))
    allocate( mklongname   (nmicro_fields))
    allocate( mkunits      (nmicro_fields))
    allocate( mkoutputscale(nmicro_fields))
    allocate( qn(nx,ny,nzm,ncrms)  )
    allocate( qpsrc(nz,ncrms)  )
    allocate( qpevp(nz,ncrms)  )

    q (dimx1_s:,dimy1_s:,1:,1:) => micro_field(:,:,:,1,:)
    qp(dimx1_s:,dimy1_s:,1:,1:) => micro_field(:,:,:,2,:)

    zero = 0

    micro_field = zero
    fluxbmk  = zero
    fluxtmk  = zero
    mkwle   = zero
    mkwsb   = zero
    mkadv   = zero
    mklsadv = zero
    mkdiff  = zero
    mstor   = zero
    mkname        = ''
    mklongname    = ''
    mkunits       = ''
    mkoutputscale = zero
    qn = zero
    qpsrc = zero
    qpevp = zero
  end subroutine allocate_micro


  subroutine deallocate_micro()
    implicit none
    deallocate(micro_field  )
    deallocate(fluxbmk   )
    deallocate(fluxtmk   )
    deallocate(mkwle    )
    deallocate(mkwsb    )
    deallocate(mkadv    )
    deallocate(mklsadv  )
    deallocate(mkdiff   )
    deallocate(mstor    )
    deallocate(mkname         )
    deallocate(mklongname     )
    deallocate(mkunits        )
    deallocate(mkoutputscale  )
    deallocate(qn  )
    deallocate(qpsrc  )
    deallocate(qpevp  )
    nullify(q )
    nullify(qp)
  end subroutine deallocate_micro


  ! required microphysics subroutines and function:
  !----------------------------------------------------------------------
  !!! Read microphysics options from prm file

  subroutine micro_setparm()
    ! no user-definable options in SAM1MOM microphysics.
  end subroutine micro_setparm

  !----------------------------------------------------------------------
  !!! Initialize microphysics:


  subroutine micro_init(ncrms,icrm)

#ifdef CLUBB_CRM
    use params, only: doclubb, doclubbnoninter ! dschanen UWM 21 May 2008
    use params, only: nclubb
#endif
    use grid, only: nrestart
    use vars, only: q0
    use params, only: dosmoke
    implicit none
    integer, intent(in) :: ncrms,icrm
    integer k, n
#ifdef CLUBB_CRM
    !  if ( nclubb /= 1 ) then
    !    write(0,*) "The namelist parameter nclubb is not equal to 1,",  &
    !      " but SAM single moment microphysics is enabled."
    !    write(0,*) "This will create unrealistic results in subsaturated grid boxes. ", &
    !      "Exiting..."
    !    call task_abort()
    !  end if
#endif

    a_bg = 1./(tbgmax-tbgmin)
    a_pr = 1./(tprmax-tprmin)
    a_gr = 1./(tgrmax-tgrmin)

    ! if(doprecip) call precip_init()

    if(nrestart.eq.0) then

#ifndef CRM
      micro_field(:,:,:,:,icrm) = 0.
      do k=1,nzm
        q(:,:,k,icrm) = q0(k,icrm)
      end do
      qn(:,:,:,icrm) = 0.
#endif

      fluxbmk(:,:,:,icrm) = 0.
      fluxtmk(:,:,:,icrm) = 0.

#ifdef CLUBB_CRM
      if ( docloud .or. doclubb ) then
#else
      if(docloud) then
#endif
#ifndef CRM
        call cloud(ncrms,icrm,micro_field,qn)
#endif
        call micro_diagnose(ncrms,icrm)
      end if
      if(dosmoke) then
        call micro_diagnose(ncrms,icrm)
      end if

    end if

    mkwle  (:,:,icrm) = 0.
    mkwsb  (:,:,icrm) = 0.
    mkadv  (:,:,icrm) = 0.
    mkdiff (:,:,icrm) = 0.
    mklsadv(:,:,icrm) = 0.
    mstor  (:,:,icrm) = 0.

    qpsrc(:,icrm) = 0.
    qpevp(:,icrm) = 0.

    mkname(1) = 'QT'
    mklongname(1) = 'TOTAL WATER (VAPOR + CONDENSATE)'
    mkunits(1) = 'g/kg'
    mkoutputscale(1) = 1.e3

    mkname(2) = 'QP'
    mklongname(2) = 'PRECIPITATING WATER'
    mkunits(2) = 'g/kg'
    mkoutputscale(2) = 1.e3

    ! set mstor to be the inital microphysical mixing ratios
    do n=1, nmicro_fields
      do k=1, nzm
        mstor(k, n,icrm) = SUM(micro_field(1:nx,1:ny,k,n,icrm))
      end do
    end do

  end subroutine micro_init

  !----------------------------------------------------------------------
  !!! fill-in surface and top boundary fluxes:
  !
  subroutine micro_flux(ncrms,icrm)

    use vars, only: fluxbq, fluxtq
    implicit none
    integer, intent(in) :: ncrms,icrm

#ifdef CLUBB_CRM
    ! Added by dschanen UWM
    use params, only: doclubb, doclubb_sfc_fluxes, docam_sfc_fluxes
    if ( doclubb .and. (doclubb_sfc_fluxes .or. docam_sfc_fluxes) ) then
      ! Add this in later
      fluxbmk(:,:,index_water_vapor,icrm) = 0.0
    else
      fluxbmk(:,:,index_water_vapor,icrm) = fluxbq(:,:,icrm)
    end if
#else
    fluxbmk(:,:,index_water_vapor,icrm) = fluxbq(:,:,icrm)
#endif /*CLUBB_CRM*/
    fluxtmk(:,:,index_water_vapor,icrm) = fluxtq(:,:,icrm)

  end subroutine micro_flux

  !----------------------------------------------------------------------
  !!! compute local microphysics processes (bayond advection and SGS diffusion):
  !
  subroutine micro_proc(ncrms,icrm)

    use grid, only: nstep,dt,icycle
    use params, only: dosmoke
    use cloud_mod
    use precip_init_mod
    use precip_proc_mod
#ifdef CLUBB_CRM
    use params, only: doclubb, doclubbnoninter ! dschanen UWM 21 May 2008
    use clubbvars, only: cloud_frac
    use vars, only:  CF3D
    use grid, only: nzm
#endif
    implicit none
    integer, intent(in) :: ncrms,icrm

    ! Update bulk coefficient
    if(doprecip.and.icycle.eq.1) call precip_init(ncrms,icrm)

    if(docloud) then
      call cloud(ncrms,icrm,micro_field,qn)
      if(doprecip) call precip_proc(ncrms,icrm,qpsrc,qpevp,micro_field,qn)
      call micro_diagnose(ncrms,icrm)
    end if
    if(dosmoke) then
      call micro_diagnose(ncrms,icrm)
    end if
#ifdef CLUBB_CRM
    if ( doclubb ) then ! -dschanen UWM 21 May 2008
      CF3D(:,:, 1:nzm,icrm) = cloud_frac(:,:,2:nzm+1) ! CF3D is used in precip_proc_clubb,
      ! so it is set here first  +++mhwang
      if(doprecip) call precip_proc_clubb(ncrms,icrm)
      call micro_diagnose(ncrms,icrm)
    end if
#endif /*CLUBB_CRM*/

  end subroutine micro_proc

  !----------------------------------------------------------------------
  !!! Diagnose arrays nessesary for dynamical core and statistics:
  !
  subroutine micro_diagnose(ncrms,icrm)

    use vars
    implicit none
    integer, intent(in) :: ncrms,icrm

    real(crm_rknd) omn, omp
    integer i,j,k

    do k=1,nzm
      do j=1,ny
        do i=1,nx
          qv(i,j,k,icrm) = q(i,j,k,icrm) - qn(i,j,k,icrm)
          omn = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(i,j,k,icrm)-tbgmin)*a_bg))
          qcl(i,j,k,icrm) = qn(i,j,k,icrm)*omn
          qci(i,j,k,icrm) = qn(i,j,k,icrm)*(1.-omn)
          omp = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(i,j,k,icrm)-tprmin)*a_pr))
          qpl(i,j,k,icrm) = qp(i,j,k,icrm)*omp
          qpi(i,j,k,icrm) = qp(i,j,k,icrm)*(1.-omp)
        end do
      end do
    end do



  end subroutine micro_diagnose

#ifdef CLUBB_CRM
  !---------------------------------------------------------------------
  subroutine micro_update()

    ! Description:
    ! This subroutine essentially does what micro_proc does but does not
    ! call any microphysics subroutines.  We need this so that CLUBB gets a
    ! properly updated value of ice fed in.
    !
    ! dschanen UWM 7 Jul 2008
    !---------------------------------------------------------------------

    !   call micro_diagnose()

    call micro_diagnose_clubb()

  end subroutine micro_update

  !---------------------------------------------------------------------
  subroutine micro_adjust( new_qv, new_qc )
    ! Description:
    ! Adjust vapor and liquid water.
    ! Microphysical variables are stored separately in
    !    SAM's dynamics + CLUBB ( e.g. qv, qcl, qci) and
    !    SAM's microphysics. (e.g. q and qn).
    ! This subroutine stores values of qv, qcl updated by CLUBB
    !   in the single-moment microphysical variables q and qn.
    !
    ! dschanen UWM 20 May 2008
    !---------------------------------------------------------------------

    use vars, only: qci

    implicit none

    real(crm_rknd), dimension(nx,ny,nzm), intent(in) :: &
    new_qv, & ! Water vapor mixing ratio that has been adjusted by CLUBB [kg/kg]
    new_qc    ! Cloud water mixing ratio that has been adjusted by CLUBB [kg/kg].
    ! For the single moment microphysics, it is liquid + ice

    q(1:nx,1:ny,1:nzm,icrm) = new_qv + new_qc ! Vapor + Liquid + Ice
    qn(1:nx,1:ny,1:nzm,icrm) = new_qc ! Liquid + Ice

    return
  end subroutine micro_adjust

  subroutine micro_diagnose_clubb()

    use vars
    use constants_clubb, only: fstderr, zero_threshold
    use error_code, only: clubb_at_least_debug_level ! Procedur

    real(crm_rknd) omn, omp
    integer i,j,k

    do k=1,nzm
      do j=1,ny
        do i=1,nx
          ! For CLUBB,  water vapor and liquid water is used
          ! so set qcl to qn while qci to zero. This also allows us to call CLUBB
          ! every nclubb th time step  (see sgs_proc in sgs.F90)

          qv(i,j,k,icrm) = q(i,j,k,icrm) - qn(i,j,k,icrm)
          ! Apply local hole-filling to vapor by converting liquid to vapor. Moist
          ! static energy should be conserved, so updating temperature is not
          ! needed here. -dschanen 31 August 2011
          if ( qv(i,j,k,icrm) < zero_threshold ) then
            qn(i,j,k,icrm) = qn(i,j,k,icrm) + qv(i,j,k,icrm)
            qv(i,j,k,icrm) = zero_threshold
            if ( qn(i,j,k,icrm) < zero_threshold ) then
              if ( clubb_at_least_debug_level( 1 ) ) then
                write(fstderr,*) "Total water at", "i =", i, "j =", j, "k =", k, "is negative.", &
                "Applying non-conservative hard clipping."
              end if
              qn(i,j,k,icrm) = zero_threshold
            end if ! cloud_liq < 0
          end if ! qv < 0

          qcl(i,j,k,icrm) = qn(i,j,k,icrm)
          qci(i,j,k,icrm) = 0.0
          omp = max(0.,min(1.,(tabs(i,j,k,icrm)-tprmin)*a_pr))
          qpl(i,j,k,icrm) = qp(i,j,k,icrm)*omp
          qpi(i,j,k,icrm) = qp(i,j,k,icrm)*(1.-omp)
        end do
      end do
    end do

  end subroutine micro_diagnose_clubb

#endif /*CLUBB_CRM*/
  !----------------------------------------------------------------------
  !!! function to compute terminal velocity for precipitating variables:
  ! In this particular case there is only one precipitating variable.

  real(crm_rknd) function term_vel_qp(ncrms,icrm,i,j,k,ind)

    use vars
    integer, intent(in) :: ncrms,icrm
    integer, intent(in) :: i,j,k,ind
    real(crm_rknd) wmax, omp, omg, qrr, qss, qgg

    term_vel_qp = 0.
    if(qp(i,j,k,icrm).gt.qp_threshold) then
      omp = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(i,j,k,icrm)-tprmin)*a_pr))
      if(omp.eq.1.) then
        term_vel_qp = vrain*(rho(k,icrm)*qp(i,j,k,icrm))**crain
      elseif(omp.eq.0.) then
        omg = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(i,j,k,icrm)-tgrmin)*a_gr))
        qgg=omg*qp(i,j,k,icrm)
        qss=qp(i,j,k,icrm)-qgg
        term_vel_qp = (omg*vgrau*(rho(k,icrm)*qgg)**cgrau &
        +(1.-omg)*vsnow*(rho(k,icrm)*qss)**csnow)
      else
        omg = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(i,j,k,icrm)-tgrmin)*a_gr))
        qrr=omp*qp(i,j,k,icrm)
        qss=qp(i,j,k,icrm)-qrr
        qgg=omg*qss
        qss=qss-qgg
        term_vel_qp = (omp*vrain*(rho(k,icrm)*qrr)**crain &
        +(1.-omp)*(omg*vgrau*(rho(k,icrm)*qgg)**cgrau &
        +(1.-omg)*vsnow*(rho(k,icrm)*qss)**csnow))
      endif
    end if
  end function term_vel_qp

  !----------------------------------------------------------------------
  !!! compute sedimentation
  !
  subroutine micro_precip_fall(ncrms,icrm)

    use vars
    use params, only : pi
    use precip_fall_mod
    implicit none
    integer, intent(in) :: ncrms,icrm

    real(crm_rknd) omega(nx,ny,nzm)
    integer ind
    integer i,j,k

    crain = b_rain / 4.
    csnow = b_snow / 4.
    cgrau = b_grau / 4.
    vrain = a_rain * gamr3 / 6. / (pi * rhor * nzeror) ** crain
    vsnow = a_snow * gams3 / 6. / (pi * rhos * nzeros) ** csnow
    vgrau = a_grau * gamg3 / 6. / (pi * rhog * nzerog) ** cgrau

    do k=1,nzm
      do j=1,ny
        do i=1,nx
          omega(i,j,k) = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(i,j,k,icrm)-tprmin)*a_pr))
        end do
      end do
    end do

    call precip_fall(ncrms,icrm, micro_field, term_vel_qp, 2, omega, ind)


  end subroutine micro_precip_fall

  !----------------------------------------------------------------------
  ! called when stepout() called

  !-----------------------------------------------------------------------
  ! Supply function that computes total water in a domain:
  !
  real(8) function total_water(ncrms,icrm)
    use vars, only : nstep,nprint,adz,dz,rho
    implicit none
    integer, intent(in) :: ncrms,icrm
    real(8) tmp
    integer i,j,k,m

    total_water = 0.
    do m=1,nmicro_fields
      if(flag_wmass(m).eq.1) then
        do k=1,nzm
          tmp = 0.
          do j=1,ny
            do i=1,nx
              tmp = tmp + micro_field(i,j,k,m,icrm)
            end do
          end do
          total_water = total_water + tmp*adz(k,icrm)*dz(icrm)*rho(k,icrm)
        end do
      end if
    end do

  end function total_water

end module microphysics
