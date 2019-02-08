module microphysics
  use cloud_mod
  use precip_init_mod
  use precip_proc_mod

  ! module for original SAM bulk microphysics
  ! Marat Khairoutdinov, 2006

  use grid, only: nx,ny,nzm,nz, dimx1_s,dimx2_s,dimy1_s,dimy2_s ! subdomain grid information
  use params, only: doprecip, docloud, doclubb, crm_rknd, asyncid
  use micro_params
  implicit none

  !----------------------------------------------------------------------
  !!! required definitions:

  integer, parameter :: nmicro_fields = 2   ! total number of prognostic water vars

  !!! microphysics prognostic variables are storred in this array:


  integer, parameter :: index_water_vapor = 1 ! index for variable that has water vapor
  integer, parameter :: index_cloud_ice = 1   ! index for cloud ice (sedimentation)

  ! both variables correspond to mass, not number
  ! SAM1MOM 3D microphysical fields are output by default.
  integer, allocatable :: flag_micro3Dout(:,:)
  integer, allocatable :: flag_precip    (:)
  integer, allocatable :: flag_wmass     (:,:)
  integer, allocatable :: flag_number    (:,:)


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
    integer :: icrm
    real(crm_rknd) :: zero
    allocate( micro_field(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, nmicro_fields))
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
    allocate( flag_micro3Dout(nmicro_fields,ncrms) )
    allocate( flag_precip    (nmicro_fields) )
    allocate( flag_wmass     (nmicro_fields,ncrms) )
    allocate( flag_number    (nmicro_fields,ncrms) )

    q (1:,dimx1_s:,dimy1_s:,1:) => micro_field(:,:,:,:,1)
    qp(1:,dimx1_s:,dimy1_s:,1:) => micro_field(:,:,:,:,2)

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
    flag_precip    (:)  = (/0,1/)
    do icrm = 1 , ncrms
      flag_micro3Dout(:,icrm)  = (/0,0/)
      flag_wmass     (:,icrm)  = (/1,1/)
      flag_number    (:,icrm)  = (/0,0/)
    enddo
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
    deallocate(flag_micro3Dout)
    deallocate(flag_precip    )
    deallocate(flag_wmass     )
    deallocate(flag_number    )
  end subroutine deallocate_micro


  ! required microphysics subroutines and function:
  !----------------------------------------------------------------------
  !!! Read microphysics options from prm file

  subroutine micro_setparm()
    ! no user-definable options in SAM1MOM microphysics.
  end subroutine micro_setparm

  !----------------------------------------------------------------------
  !!! Initialize microphysics:


  subroutine micro_init(ncrms)

#ifdef CLUBB_CRM
    use params, only: doclubb, doclubbnoninter ! dschanen UWM 21 May 2008
    use params, only: nclubb
#endif
    use grid, only: nrestart
    use vars, only: q0
    use params, only: dosmoke
    implicit none
    integer, intent(in) :: ncrms
    integer k, n,icrm
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

    if(nrestart.eq.0) then

#ifndef CRM
      micro_field(icrm,:,:,:,:) = 0.
      do k=1,nzm
        q(icrm,:,:,k) = q0(k,icrm)
      end do
      qn(:,:,:,icrm) = 0.
#endif

      do icrm = 1 , ncrms
        fluxbmk(:,:,:,icrm) = 0.
        fluxtmk(:,:,:,icrm) = 0.
      enddo

#ifdef CLUBB_CRM
      if ( docloud .or. doclubb ) then
#else
      if(docloud) then
#endif
#ifndef CRM
        call cloud(ncrms,q,qp,qn)
        !$acc wait(asyncid)
#endif
        call micro_diagnose(ncrms)
        !$acc wait(asyncid)
      end if
      if(dosmoke) then
        call micro_diagnose(ncrms)
        !$acc wait(asyncid)
      end if
    end if

    do icrm = 1 , ncrms
      mkwle  (:,:,icrm) = 0.
      mkwsb  (:,:,icrm) = 0.
      mkadv  (:,:,icrm) = 0.
      mkdiff (:,:,icrm) = 0.
      mklsadv(:,:,icrm) = 0.
      mstor  (:,:,icrm) = 0.

      qpsrc(:,icrm) = 0.
      qpevp(:,icrm) = 0.

      ! set mstor to be the inital microphysical mixing ratios
      do n=1, nmicro_fields
        do k=1, nzm
          mstor(k, n,icrm) = SUM(micro_field(icrm,1:nx,1:ny,k,n))
        end do
      end do
    enddo

    mkname(1) = 'QT'
    mklongname(1) = 'TOTAL WATER (VAPOR + CONDENSATE)'
    mkunits(1) = 'g/kg'
    mkoutputscale(1) = 1.e3

    mkname(2) = 'QP'
    mklongname(2) = 'PRECIPITATING WATER'
    mkunits(2) = 'g/kg'
    mkoutputscale(2) = 1.e3

  end subroutine micro_init

  !----------------------------------------------------------------------
  !!! fill-in surface and top boundary fluxes:
  subroutine micro_flux(ncrms)
    use vars, only: fluxbq, fluxtq
    implicit none
    integer, intent(in) :: ncrms
    integer :: icrm, i, j

#ifdef CLUBB_CRM
    ! Added by dschanen UWM
    use params, only: doclubb, doclubb_sfc_fluxes, docam_sfc_fluxes
    do icrm = 1 , ncrms
      if ( doclubb .and. (doclubb_sfc_fluxes .or. docam_sfc_fluxes) ) then
        ! Add this in later
        fluxbmk(:,:,index_water_vapor,icrm) = 0.0
      else
        fluxbmk(:,:,index_water_vapor,icrm) = fluxbq(:,:,icrm)
      end if
    enddo
#else
    !$acc parallel loop collapse(3) copyin(fluxbq) copy(fluxbmk) async(asyncid)
    do icrm = 1 , ncrms
      do j = 1 , ny
        do i = 1 , nx
          fluxbmk(i,j,index_water_vapor,icrm) = fluxbq(i,j,icrm)
        enddo
      enddo
    enddo
#endif
    !$acc parallel loop collapse(3) copyin(fluxtq) copy(fluxtmk) async(asyncid)
    do icrm = 1 , ncrms
      do j = 1 , ny
        do i = 1 , nx
          fluxtmk(i,j,index_water_vapor,icrm) = fluxtq(i,j,icrm)
        enddo
      enddo
    enddo
  end subroutine micro_flux

  !----------------------------------------------------------------------
  !!! compute local microphysics processes (bayond advection and SGS diffusion):
  !
  subroutine micro_proc(ncrms)
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
    integer, intent(in) :: ncrms
    integer :: icrm

    ! Update bulk coefficient
    if(doprecip.and.icycle.eq.1) call precip_init(ncrms)

    if(docloud) then
#ifdef __PGI
      !Passing q and qp via the first element because PGI has a bug with pointers here
      call cloud(ncrms,q(1,dimx1_s,dimy1_s,1),qp(1,dimx1_s,dimy1_s,1),qn)
      if(doprecip) call precip_proc(ncrms,qpsrc,qpevp,q(1,dimx1_s,dimy1_s,1),qp(1,dimx1_s,dimy1_s,1),qn)
#else
      call cloud(ncrms, q, qp, qn)
      if(doprecip) call precip_proc(ncrms, qpsrc, qpevp, q, qp, qn)
#endif
      call micro_diagnose(ncrms)
    end if
    if(dosmoke) then
      call micro_diagnose(ncrms)
    end if
#ifdef CLUBB_CRM
    if ( doclubb ) then ! -dschanen UWM 21 May 2008
      do icrm = 1 , ncrms
        CF3D(:,:, 1:nzm,icrm) = cloud_frac(:,:,2:nzm+1) ! CF3D is used in precip_proc_clubb,
        ! so it is set here first  +++mhwang
        if(doprecip) call precip_proc_clubb(ncrms,icrm)
      enddo
      call micro_diagnose(ncrms)
    end if
#endif /*CLUBB_CRM*/
  end subroutine micro_proc

  !----------------------------------------------------------------------
  !!! Diagnose arrays nessesary for dynamical core and statistics:
  !
  subroutine micro_diagnose(ncrms)
    use vars
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) omn, omp
    integer i,j,k,icrm

    !$acc parallel loop collapse(4) copy(qv,q,qn,tabs,qp,qpl,qpi,qcl,qci) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            qv(i,j,k,icrm) = q(icrm,i,j,k) - qn(i,j,k,icrm)
            omn = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(i,j,k,icrm)-tbgmin)*a_bg))
            qcl(i,j,k,icrm) = qn(i,j,k,icrm)*omn
            qci(i,j,k,icrm) = qn(i,j,k,icrm)*(1.-omn)
            omp = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(i,j,k,icrm)-tprmin)*a_pr))
            qpl(i,j,k,icrm) = qp(icrm,i,j,k)*omp
            qpi(i,j,k,icrm) = qp(icrm,i,j,k)*(1.-omp)
          end do
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

    q(icrm,1:nx,1:ny,1:nzm) = new_qv + new_qc ! Vapor + Liquid + Ice
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

          qv(i,j,k,icrm) = q(icrm,i,j,k) - qn(i,j,k,icrm)
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
          qpl(i,j,k,icrm) = qp(icrm,i,j,k)*omp
          qpi(i,j,k,icrm) = qp(icrm,i,j,k)*(1.-omp)
        end do
      end do
    end do

  end subroutine micro_diagnose_clubb

#endif /*CLUBB_CRM*/
  !----------------------------------------------------------------------
  !!! function to compute terminal velocity for precipitating variables:
  ! In this particular case there is only one precipitating variable.

  real(crm_rknd) function term_vel_qp(ncrms,icrm,i,j,k,ind,qploc,rho,tabs,qp_threshold,tprmin,&
                                      a_pr,vrain,crain,tgrmin,a_gr,vgrau,cgrau,vsnow,csnow)
    !$acc routine seq
    implicit none
    integer, intent(in) :: ncrms,icrm
    integer, intent(in) :: i,j,k,ind
    real(crm_rknd), intent(in) :: qploc
    real(crm_rknd), intent(in) :: rho(nzm,ncrms), tabs(nx, ny, nzm, ncrms)
    real(crm_rknd), intent(in) :: qp_threshold,tprmin,a_pr,vrain,crain,tgrmin,a_gr,vgrau,cgrau,vsnow,csnow
    real(crm_rknd) wmax, omp, omg, qrr, qss, qgg

    term_vel_qp = 0.
    if(qploc.gt.qp_threshold) then
      omp = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(i,j,k,icrm)-tprmin)*a_pr))
      if(omp.eq.1.) then
        term_vel_qp = vrain*(rho(k,icrm)*qploc)**crain
      elseif(omp.eq.0.) then
        omg = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(i,j,k,icrm)-tgrmin)*a_gr))
        qgg=omg*qploc
        qss=qploc-qgg
        term_vel_qp = (omg*vgrau*(rho(k,icrm)*qgg)**cgrau &
        +(1.-omg)*vsnow*(rho(k,icrm)*qss)**csnow)
      else
        omg = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(i,j,k,icrm)-tgrmin)*a_gr))
        qrr=omp*qploc
        qss=qploc-qrr
        qgg=omg*qss
        qss=qss-qgg
        term_vel_qp = (omp*vrain*(rho(k,icrm)*qrr)**crain + (1.-omp)*(omg*vgrau*(rho(k,icrm)*qgg)**cgrau + &
                      (1.-omg)*vsnow*(rho(k,icrm)*qss)**csnow))
      endif
    end if
  end function term_vel_qp

  !----------------------------------------------------------------------
  !!! compute sedimentation
  !
  subroutine micro_precip_fall(ncrms)
    use vars
    use params, only : pi
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) omega(nx,ny,nzm,ncrms)
    integer ind
    integer i,j,k,icrm

    !$acc enter data create(omega) async(asyncid)

    crain = b_rain / 4.
    csnow = b_snow / 4.
    cgrau = b_grau / 4.
    vrain = a_rain * gamr3 / 6. / (pi * rhor * nzeror) ** crain
    vsnow = a_snow * gams3 / 6. / (pi * rhos * nzeros) ** csnow
    vgrau = a_grau * gamg3 / 6. / (pi * rhog * nzerog) ** cgrau

    !$acc parallel loop collapse(4) copyin(tabs) copy(omega) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            omega(i,j,k,icrm) = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(i,j,k,icrm)-tprmin)*a_pr))
          end do
        end do
      end do
    end do

    call precip_fall(ncrms, 2, omega, ind)

    !$acc exit data delete(omega) async(asyncid)

  end subroutine micro_precip_fall


  subroutine precip_fall(ncrms,hydro_type, omega, ind)
    !     positively definite monotonic advection with non-oscillatory option
    !     and gravitational sedimentation
    use vars
    use params
    implicit none
    integer, intent(in) :: ncrms
    integer :: hydro_type   ! 0 - all liquid, 1 - all ice, 2 - mixed
    real(crm_rknd) :: omega(nx,ny,nzm,ncrms)   !  = 1: liquid, = 0: ice;  = 0-1: mixed : used only when hydro_type=2
    integer :: ind
    ! Terminal velocity fnction
    ! Local:
    real(crm_rknd) :: mx(nx,ny,nzm,ncrms),mn(nx,ny,nzm,ncrms), lfac(nx,ny,nz,ncrms)
    real(crm_rknd) :: www(nx,ny,nz,ncrms),fz(nx,ny,nz,ncrms)
    real(crm_rknd) :: eps
    integer :: i,j,k,kc,kb,icrm
    logical :: nonos
    real(crm_rknd) :: y,pp,pn
    real(crm_rknd) :: lat_heat, wmax
    real(crm_rknd) :: wp(nx,ny,nzm,ncrms), tmp_qp(nx,ny,nzm,ncrms), irhoadz(nzm,ncrms), iwmax(nzm,ncrms), &
                      rhofac(nzm,ncrms), prec_cfl
    integer nprec, iprec
    real(crm_rknd) :: flagstat, tmp
    real(crm_rknd), pointer :: qp(:,:,:,:)  ! total precipitating water

    !Statement functions
    pp(y)= max(real(0.,crm_rknd),y)
    pn(y)=-min(real(0.,crm_rknd),y)

    qp(1:,dimx1_s:,dimy1_s:,1:) => micro_field(:,:,:,:,2)

    eps = 1.e-10
    nonos = .true.

    !$acc enter data create(mx,mn,lfac,www,fz,wp,tmp_qp,irhoadz,iwmax,rhofac) async(asyncid)

    !$acc parallel loop gang vector collapse(2) copyin(rho,adz,dz) copy(rhofac,irhoadz,iwmax) async(asyncid)
    do icrm = 1 , ncrms
      do k = 1,nzm
        rhofac(k,icrm) = sqrt(1.29/rho(k,icrm))
        irhoadz(k,icrm) = 1./(rho(k,icrm)*adz(k,icrm)) ! Useful factor
        kb = max(1,k-1)
        wmax       = dz(icrm)*adz(kb,icrm)/dtn   ! Velocity equivalent to a cfl of 1.0.
        iwmax(k,icrm)   = 1./wmax
      enddo
    enddo

    ! 	Add sedimentation of precipitation field to the vert. vel.
    prec_cfl = 0.
    !$acc parallel loop gang vector collapse(4) copyin(omega,rhofac,micro_field,rho,tabs,iwmax,rhow,dz) &
    !$acc&                                      copy(prec_cfl,wp,fz,www,lfac,flagstat) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            select case (hydro_type)
            case(0)
              lfac(i,j,k,icrm) = fac_cond
              flagstat = 1.
            case(1)
              lfac(i,j,k,icrm) = fac_sub
              flagstat = 1.
            case(2)
              lfac(i,j,k,icrm) = fac_cond + (1-omega(i,j,k,icrm))*fac_fus
              flagstat = 1.
            case(3)
              lfac(i,j,k,icrm) = 0.
              flagstat = 0.
            case default
              if(masterproc) then
                !print*, 'unknown hydro_type in precip_fall. exitting ...'
                !call task_abort
              endif
            end select
            wp(i,j,k,icrm)=rhofac(k,icrm)*term_vel_qp(ncrms,icrm,i,j,k,ind,micro_field(icrm,i,j,k,2),rho(:,:),&
                                                      tabs(:,:,:,:),qp_threshold,tprmin,a_pr,vrain,crain,tgrmin,&
                                                      a_gr,vgrau,cgrau,vsnow,csnow)
            tmp = wp(i,j,k,icrm)*iwmax(k,icrm)
            !$acc atomic update
            prec_cfl = max(prec_cfl,tmp) ! Keep column maximum CFL
            wp(i,j,k,icrm) = -wp(i,j,k,icrm)*rhow(k,icrm)*dtn/dz(icrm)
            if (k == 1) then
              fz(i,j,nz,icrm)=0.
              www(i,j,nz,icrm)=0.
              lfac(i,j,nz,icrm)=0
            endif
          enddo  ! k
        enddo
      enddo
    enddo

    ! If maximum CFL due to precipitation velocity is greater than 0.9,
    ! take more than one advection step to maintain stability.
    if (prec_cfl.gt.0.9) then
      nprec = CEILING(prec_cfl/0.9)
      !$acc parallel loop gang vector collapse(4) copy(wp) async(asyncid)
      do icrm = 1 , ncrms
        do k = 1,nzm
          do j=1,ny
            do i=1,nx
              ! wp already includes factor of dt, so reduce it by a
              ! factor equal to the number of precipitation steps.
              wp(i,j,k,icrm) = wp(i,j,k,icrm)/real(nprec,crm_rknd)
            enddo
          enddo
        enddo
      enddo
    else
      nprec = 1
    endif

    !  loop over iterations
    do iprec = 1,nprec
      !$acc parallel loop gang vector collapse(4) copyin(qp) copy(tmp_qp) async(asyncid)
      do icrm = 1 , ncrms
        do k = 1,nzm
          do j=1,ny
            do i=1,nx
              tmp_qp(i,j,k,icrm) = qp(icrm,i,j,k) ! Temporary array for qp in this column
            enddo
          enddo
        enddo
      enddo

      !$acc parallel loop gang vector collapse(4) copyin(tmp_qp,wp) copy(mn,mx,fz) async(asyncid)
      do icrm = 1 , ncrms
        do k=1,nzm
          do j=1,ny
            do i=1,nx
              if(nonos) then
                kc=min(nzm,k+1)
                kb=max(1,k-1)
                mx(i,j,k,icrm)=max(tmp_qp(i,j,kb,icrm),tmp_qp(i,j,kc,icrm),tmp_qp(i,j,k,icrm))
                mn(i,j,k,icrm)=min(tmp_qp(i,j,kb,icrm),tmp_qp(i,j,kc,icrm),tmp_qp(i,j,k,icrm))
              endif  ! nonos
              ! Define upwind precipitation flux
              fz(i,j,k,icrm)=tmp_qp(i,j,k,icrm)*wp(i,j,k,icrm)
            enddo
          enddo
        enddo
      enddo

      !$acc parallel loop gang vector collapse(4) copyin(fz,irhoadz) copy(tmp_qp) async(asyncid)
      do icrm = 1 , ncrms
        do k=1,nzm
          do j=1,ny
            do i=1,nx
              kc=k+1
              tmp_qp(i,j,k,icrm)=tmp_qp(i,j,k,icrm)-(fz(i,j,kc,icrm)-fz(i,j,k,icrm))*irhoadz(k,icrm) !Update temporary qp
            enddo
          enddo
        enddo
      enddo

      !$acc parallel loop gang vector collapse(4) copyin(wp,irhoadz,tmp_qp,wp) copy(www) async(asyncid)
      do icrm = 1 , ncrms
        do k=1,nzm
          do j=1,ny
            do i=1,nx
              ! Also, compute anti-diffusive correction to previous
              ! (upwind) approximation to the flux
              kb=max(1,k-1)
              ! The precipitation velocity is a cell-centered quantity,
              ! since it is computed from the cell-centered
              ! precipitation mass fraction.  Therefore, a reformulated
              ! anti-diffusive flux is used here which accounts for
              ! this and results in reduced numerical diffusion.
              www(i,j,k,icrm) = 0.5*(1.+wp(i,j,k,icrm)*irhoadz(k,icrm))*(tmp_qp(i,j,kb,icrm)*wp(i,j,kb,icrm) - &
                                     tmp_qp(i,j,k,icrm)*wp(i,j,k,icrm)) ! works for wp(k)<0
            enddo
          enddo
        enddo
      enddo

      !---------- non-osscilatory option ---------------
      if(nonos) then
        !$acc parallel loop gang vector collapse(4) copyin(tmp_qp,rho,adz,www) copy(mn,mx) async(asyncid)
        do icrm = 1 , ncrms
          do k=1,nzm
            do j=1,ny
              do i=1,nx
                kc=min(nzm,k+1)
                kb=max(1,k-1)
                mx(i,j,k,icrm)=max(tmp_qp(i,j,kb,icrm),tmp_qp(i,j,kc,icrm),tmp_qp(i,j,k,icrm),mx(i,j,k,icrm))
                mn(i,j,k,icrm)=min(tmp_qp(i,j,kb,icrm),tmp_qp(i,j,kc,icrm),tmp_qp(i,j,k,icrm),mn(i,j,k,icrm))
                kc=min(nzm,k+1)
                mx(i,j,k,icrm)=rho(k,icrm)*adz(k,icrm)*(mx(i,j,k,icrm)-tmp_qp(i,j,k,icrm))/(pn(www(i,j,kc,icrm)) + pp(www(i,j,k,icrm))+eps)
                mn(i,j,k,icrm)=rho(k,icrm)*adz(k,icrm)*(tmp_qp(i,j,k,icrm)-mn(i,j,k,icrm))/(pp(www(i,j,kc,icrm)) + pn(www(i,j,k,icrm))+eps)
              enddo
            enddo
          enddo
        enddo
        !$acc parallel loop gang vector collapse(4) copyin(www,mn,mx) copy(fz) async(asyncid)
        do icrm = 1 , ncrms
          do k=1,nzm
            do j=1,ny
              do i=1,nx
                kb=max(1,k-1)
                ! Add limited flux correction to fz(k).
                fz(i,j,k,icrm) = fz(i,j,k,icrm) + pp(www(i,j,k,icrm))*min(real(1.,crm_rknd),mx(i,j,k,icrm), mn(i,j,kb,icrm)) - &
                                                  pn(www(i,j,k,icrm))*min(real(1.,crm_rknd),mx(i,j,kb,icrm),mn(i,j,k,icrm)) ! Anti-diffusive flux
              enddo
            enddo
          enddo
        enddo
      endif ! nonos

      ! Update precipitation mass fraction and liquid-ice static
      ! energy using precipitation fluxes computed in this column.
      !$acc parallel loop gang vector collapse(4) copyin(fz,irhoadz,lfac,flagstat,omega) &
      !$acc&                                      copy(qp,qpfall,t,tlat,precflux,precsfc,precssfc,prec_xy) async(asyncid)
      do icrm = 1 , ncrms
        do j=1,ny
          do i=1,nx
            do k=1,nzm
              kc=k+1
              ! Update precipitation mass fraction.
              ! Note that fz is the total flux, including both the
              ! upwind flux and the anti-diffusive correction.
              qp(icrm,i,j,k)=qp(icrm,i,j,k)-(fz(i,j,kc,icrm)-fz(i,j,k,icrm))*irhoadz(k,icrm)
              qpfall(k,icrm)=qpfall(k,icrm)-(fz(i,j,kc,icrm)-fz(i,j,k,icrm))*irhoadz(k,icrm)*flagstat  ! For qp budget
              lat_heat = -(lfac(i,j,kc,icrm)*fz(i,j,kc,icrm)-lfac(i,j,k,icrm)*fz(i,j,k,icrm))*irhoadz(k,icrm)
              t(icrm,i,j,k)=t(icrm,i,j,k)-lat_heat
              !$acc atomic update
              tlat(k,icrm)=tlat(k,icrm)-lat_heat            ! For energy budget
              tmp = fz(i,j,k,icrm)*flagstat
              !$acc atomic update
              precflux(k,icrm) = precflux(k,icrm) - tmp   ! For statistics
              if (k == 1) then
                precsfc(i,j,icrm) = precsfc(i,j,icrm) - fz(i,j,1,icrm)*flagstat ! For statistics
                precssfc(i,j,icrm) = precssfc(i,j,icrm) - fz(i,j,1,icrm)*(1.-omega(i,j,1,icrm))*flagstat ! For statistics
                prec_xy(i,j,icrm) = prec_xy(i,j,icrm) - fz(i,j,1,icrm)*flagstat ! For 2D output
              endif
            enddo
          enddo
        enddo
      enddo

      if (iprec.lt.nprec) then
        ! Re-compute precipitation velocity using new value of qp.
        !$acc parallel loop gang vector collapse(4) copyin(rhofac,micro_field,rho,tabs,rhow,dz) copy(wp,fz,www,lfac) async(asyncid)
        do icrm = 1 , ncrms
          do j=1,ny
            do i=1,nx
              do k=1,nzm
                !Passing variables via first index because of PGI bug with pointers
                wp(i,j,k,icrm) = rhofac(k,icrm)*term_vel_qp(ncrms,icrm,i,j,k,ind,micro_field(icrm,i,j,k,2),rho(1,1),&
                                 tabs(1,1,1,1),qp_threshold,tprmin,a_pr,vrain,crain,tgrmin,a_gr,vgrau,cgrau,vsnow,csnow)
                ! Decrease precipitation velocity by factor of nprec
                wp(i,j,k,icrm) = -wp(i,j,k,icrm)*rhow(k,icrm)*dtn/dz(icrm)/real(nprec,crm_rknd)
                ! Note: Don't bother checking CFL condition at each
                ! substep since it's unlikely that the CFL will
                ! increase very much between substeps when using
                ! monotonic advection schemes.
                if (k == 1) then
                  fz(i,j,nz,icrm)=0.
                  www(i,j,nz,icrm)=0.
                  lfac(i,j,nz,icrm)=0.
                endif
              enddo
            enddo
          enddo
        enddo
      endif

    enddo
    
    !$acc exit data delete(mx,mn,lfac,www,fz,wp,tmp_qp,irhoadz,iwmax,rhofac) async(asyncid)

  end subroutine precip_fall

  !----------------------------------------------------------------------
  ! called when stepout() called

  !-----------------------------------------------------------------------
  ! Supply function that computes total water in a domain:
  !
  real(8) function total_water(ncrms,icrm)
    use vars, only : nstep,adz,dz,rho
    implicit none
    integer, intent(in) :: ncrms,icrm
    real(8) tmp
    integer i,j,k,m

    total_water = 0.
    do m=1,nmicro_fields
      if(flag_wmass(m,icrm).eq.1) then
        do k=1,nzm
          tmp = 0.
          do j=1,ny
            do i=1,nx
              tmp = tmp + micro_field(icrm,i,j,k,m)
            end do
          end do
          total_water = total_water + tmp*adz(k,icrm)*dz(icrm)*rho(k,icrm)
        end do
      end if
    end do

  end function total_water

end module microphysics
