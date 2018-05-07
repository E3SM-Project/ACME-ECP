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

  real(crm_rknd), allocatable, target :: micro_field(:,:,:,:,:) !REDIM

  integer, parameter :: flag_wmass(nmicro_fields) = (/1,1/)
  integer, parameter :: index_water_vapor = 1 ! index for variable that has water vapor
  integer, parameter :: index_cloud_ice = 1   ! index for cloud ice (sedimentation)
  integer, parameter :: flag_precip(nmicro_fields) = (/0,1/)

  ! both variables correspond to mass, not number
  integer, parameter :: flag_number(nmicro_fields) = (/0,0/)

  ! SAM1MOM 3D microphysical fields are output by default.
  integer, parameter :: flag_micro3Dout(nmicro_fields) = (/0,0/)

  real(crm_rknd), allocatable :: fluxbmk (:,:,:,:) !REDIM ! surface flux of tracers
  real(crm_rknd), allocatable :: fluxtmk (:,:,:,:) !REDIM ! top boundary flux of tracers

  !!! these arrays are needed for output statistics:

  real(crm_rknd), allocatable :: mkwle  (:,:,:) !REDIM  ! resolved vertical flux
  real(crm_rknd), allocatable :: mkwsb  (:,:,:) !REDIM  ! SGS vertical flux
  real(crm_rknd), allocatable :: mkadv  (:,:,:) !REDIM  ! tendency due to vertical advection
  real(crm_rknd), allocatable :: mklsadv(:,:,:) !REDIM  ! tendency due to large-scale vertical advection
  real(crm_rknd), allocatable :: mkdiff (:,:,:) !REDIM  ! tendency due to vertical diffusion
  real(crm_rknd), allocatable :: mstor  (:,:,:) !REDIM  ! storage terms of microphysical variables

  !======================================================================
  ! UW ADDITIONS

  !bloss: arrays with names/units for microphysical outputs in statistics.
  character*3, dimension(nmicro_fields) :: mkname
  character*80, dimension(nmicro_fields) :: mklongname
  character*10, dimension(nmicro_fields) :: mkunits
  real(crm_rknd), dimension(nmicro_fields) :: mkoutputscale

  ! END UW ADDITIONS
  !======================================================================

  !------------------------------------------------------------------
  ! Optional (internal) definitions)

  ! make aliases for prognostic variables:
  ! note that the aliases should be local to microphysics

  real(crm_rknd), pointer :: q (:,:,:,:) !REDIM   ! total nonprecipitating water
  real(crm_rknd), pointer :: qp(:,:,:,:) !REDIM  ! total precipitating water
  real(crm_rknd), allocatable :: qn(:,:,:,:) !REDIM  ! cloud condensate (liquid + ice)

  real(crm_rknd), allocatable :: qpsrc(:,:) !REDIM  ! source of precipitation microphysical processes
  real(crm_rknd), allocatable :: qpevp(:,:) !REDIM  ! sink of precipitating water due to evaporation

  real(crm_rknd) vrain, vsnow, vgrau, crain, csnow, cgrau  ! precomputed coefs for precip terminal velocity

CONTAINS

  subroutine allocate_microphysics(ncrms)
    use openacc_utils
    implicit none
    integer, intent(in) :: ncrms
    allocate( micro_field(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, nmicro_fields) )
    allocate( fluxbmk (ncrms,nx, ny, 1:nmicro_fields) )
    allocate( fluxtmk (ncrms,nx, ny, 1:nmicro_fields) )
    allocate( mkwle   (ncrms,nz,1:nmicro_fields) )
    allocate( mkwsb   (ncrms,nz,1:nmicro_fields) )
    allocate( mkadv   (ncrms,nz,1:nmicro_fields) )
    allocate( mklsadv (ncrms,nz,1:nmicro_fields) )
    allocate( mkdiff  (ncrms,nz,1:nmicro_fields) )
    allocate( mstor   (ncrms,nz,1:nmicro_fields) )
    allocate( qn      (ncrms,nx,ny,nzm) )
    allocate( qpsrc   (ncrms,nz) )
    allocate( qpevp   (ncrms,nz) )
    q (1:,dimx1_s:,dimy1_s:,1:) => micro_field(1:ncrms,dimx1_s:dimx2_s,dimy1_s:dimy2_s,1:nzm,1)
    qp(1:,dimx1_s:,dimy1_s:,1:) => micro_field(1:ncrms,dimx1_s:dimx2_s,dimy1_s:dimy2_s,1:nzm,2)

    micro_field  = 0
    fluxbmk      = 0
    fluxtmk      = 0
    mkwle        = 0
    mkwsb        = 0
    mkadv        = 0
    mklsadv      = 0
    mkdiff       = 0
    mstor        = 0
    q            = 0
    qp           = 0
    qn           = 0
    qpsrc        = 0
    qpevp        = 0
  end subroutine allocate_microphysics

  subroutine deallocate_microphysics
    implicit none
    deallocate( micro_field )
    deallocate( fluxbmk  )
    deallocate( fluxtmk  )
    deallocate( mkwle )
    deallocate( mkwsb )
    deallocate( mkadv )
    deallocate( mklsadv )
    deallocate( mkdiff )
    deallocate( mstor )
    deallocate( qn )
    deallocate( qpsrc )
    deallocate( qpevp )
  end subroutine deallocate_microphysics

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
    use vars
    use params, only: dosmoke
    implicit none
    integer, intent(in) :: ncrms
    integer k, n, icrm
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

    ! if(doprecip) call precip_init(ncrms,icrm)

    if(nrestart.eq.0) then

#ifndef CRM
    !$acc parallel loop gang vector collapse(4) default(present) async(1)
    do k = 1 , nzm
      do j = dimy1_s,dimy2_s
        do i = dimx1_x,dimx2_s
          do icrm = 1 , ncrms
            micro_field(icrm,i,j,k,:) = 0.
          enddo
        enddo
      enddo
    enddo
    !$acc parallel loop gang vector collapse(4) default(present) async(1)
    do k = 1 , nzm
      do j = dimy1_s,dimy2_s
        do i = dimx1_x,dimx2_s
          do icrm = 1 , ncrms
            q(icrm,i,j,k) = q0(icrm,k)
          enddo
        enddo
      enddo
    enddo
    !$acc parallel loop gang vector collapse(4) default(present) async(1)
    do k = 1 , nzm
      do j = 1 , nx
        do i = 1 , ny
          do icrm = 1 , ncrms
            qn(icrm,i,j,k) = 0.
          enddo
        enddo
      enddo
    enddo
#endif

#ifdef CLUBB_CRM
      if ( docloud .or. doclubb ) then
#else
      if(docloud) then
#endif
#ifndef CRM
        call cloud(q,qn,qp,ncrms)
#endif
        call micro_diagnose(ncrms)
      end if
      if(dosmoke) then
        call micro_diagnose(ncrms)
      end if

    end if

    mkname(1) = 'QT'
    mklongname(1) = 'TOTAL WATER (VAPOR + CONDENSATE)'
    mkunits(1) = 'g/kg'
    mkoutputscale(1) = 1.e3

    mkname(2) = 'QP'
    mklongname(2) = 'PRECIPITATING WATER'
    mkunits(2) = 'g/kg'
    mkoutputscale(2) = 1.e3

    ! set mstor to be the inital microphysical mixing ratios
    !$acc parallel loop gang vector collapse(3) default(present) async(1)
    do n=1, nmicro_fields
      do k=1, nzm
        do icrm = 1 , ncrms
          mstor(icrm,k, n) = SUM(micro_field(icrm,1:nx,1:ny,k,n))
        end do
      end do
    end do

  end subroutine micro_init

  !----------------------------------------------------------------------
  !!! fill-in surface and top boundary fluxes:
  !
  subroutine micro_flux(ncrms)

    use vars, only: fluxbq, fluxtq

#ifdef CLUBB_CRM
    ! Added by dschanen UWM
    use params, only: doclubb, doclubb_sfc_fluxes, docam_sfc_fluxes
#endif
    implicit none
    integer, intent(in) :: ncrms
    integer :: i,j,icrm
    !$acc parallel loop gang vector collapse(3) default(present) async(1)
    do j = 1 , ny
      do i = 1 , nx
        do icrm = 1 , ncrms
#ifdef CLUBB_CRM
          if ( doclubb .and. (doclubb_sfc_fluxes .or. docam_sfc_fluxes) ) then
            ! Add this in later
            fluxbmk(icrm,i,j,index_water_vapor) = 0.0
          else
            fluxbmk(icrm,i,j,index_water_vapor) = fluxbq(icrm,i,j)
          end if
#else
          fluxbmk(icrm,i,j,index_water_vapor) = fluxbq(icrm,i,j)
#endif /*CLUBB_CRM*/
          fluxtmk(icrm,i,j,index_water_vapor) = fluxtq(icrm,i,j)
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

    ! Update bulk coefficient
    if(doprecip.and.icycle.eq.1) call precip_init(ncrms)

    if(docloud) then
      call cloud(q,qn,qp,ncrms)
      if(doprecip) call precip_proc(qpsrc,qpevp,qp,q,qn,ncrms)
      call micro_diagnose(ncrms)
    end if
    if(dosmoke) then
      call micro_diagnose(ncrms)
    end if
#ifdef CLUBB_CRM
    if ( doclubb ) then ! -dschanen UWM 21 May 2008
      CF3D(:,:,:, 1:nzm) = cloud_frac(:,:,2:nzm+1) ! CF3D is used in precip_proc_clubb,
      ! so it is set here first  +++mhwang
      !     if(doprecip) call precip_proc()
      do icrm = 1 , ncrms
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

    !$acc parallel loop gang vector collapse(4) default(present) async(1)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            qv(icrm,i,j,k) = q(icrm,i,j,k) - qn(icrm,i,j,k)
            omn = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(icrm,i,j,k)-tbgmin)*a_bg))
            qcl(icrm,i,j,k) = qn(icrm,i,j,k)*omn
            qci(icrm,i,j,k) = qn(icrm,i,j,k)*(1.-omn)
            omp = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(icrm,i,j,k)-tprmin)*a_pr))
            qpl(icrm,i,j,k) = qp(icrm,i,j,k)*omp
            qpi(icrm,i,j,k) = qp(icrm,i,j,k)*(1.-omp)
          end do
        end do
      end do
    end do
  end subroutine micro_diagnose

#ifdef CLUBB_CRM
  !---------------------------------------------------------------------
  subroutine micro_update(ncrms,icrm)

    ! Description:
    ! This subroutine essentially does what micro_proc does but does not
    ! call any microphysics subroutines.  We need this so that CLUBB gets a
    ! properly updated value of ice fed in.
    !
    ! dschanen UWM 7 Jul 2008
    !---------------------------------------------------------------------

    !   call cloud()
    !   call micro_diagnose()

    call micro_diagnose_clubb(ncrms,icrm)

  end subroutine micro_update

  !---------------------------------------------------------------------
  subroutine micro_adjust( new_qv, new_qc ,ncrms,icrm)
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
    integer, intent(in) :: ncrms,icrm

    real(crm_rknd), dimension(ncrms,nx,ny,nzm), intent(in) :: &
    new_qv, & ! Water vapor mixing ratio that has been adjusted by CLUBB [kg/kg]
    new_qc    ! Cloud water mixing ratio that has been adjusted by CLUBB [kg/kg].
    ! For the single moment microphysics, it is liquid + ice

    q(icrm,1:nx,1:ny,1:nzm) = new_qv + new_qc ! Vapor + Liquid + Ice
    qn(icrm,1:nx,1:ny,1:nzm) = new_qc ! Liquid + Ice

    return
  end subroutine micro_adjust

  subroutine micro_diagnose_clubb(ncrms,icrm)

    use vars
    use constants_clubb, only: fstderr, zero_threshold
    use error_code, only: clubb_at_least_debug_level ! Procedur
    implicit none
    integer, intent(in) :: ncrms,icrm

    real(crm_rknd) omn, omp
    integer i,j,k

    do k=1,nzm
      do j=1,ny
        do i=1,nx
          ! For CLUBB,  water vapor and liquid water is used
          ! so set qcl to qn while qci to zero. This also allows us to call CLUBB
          ! every nclubb th time step  (see sgs_proc in sgs.F90)

          qv(icrm,i,j,k) = q(icrm,i,j,k) - qn(icrm,i,j,k)
          ! Apply local hole-filling to vapor by converting liquid to vapor. Moist
          ! static energy should be conserved, so updating temperature is not
          ! needed here. -dschanen 31 August 2011
          if ( qv(icrm,i,j,k) < zero_threshold ) then
            qn(icrm,i,j,k) = qn(icrm,i,j,k) + qv(icrm,i,j,k)
            qv(icrm,i,j,k) = zero_threshold
            if ( qn(icrm,i,j,k) < zero_threshold ) then
              if ( clubb_at_least_debug_level( 1 ) ) then
                write(fstderr,*) "Total water at", "i =", i, "j =", j, "k =", k, "is negative.", &
                "Applying non-conservative hard clipping."
              end if
              qn(icrm,i,j,k) = zero_threshold
            end if ! cloud_liq < 0
          end if ! qv < 0

          qcl(icrm,i,j,k) = qn(icrm,i,j,k)
          qci(icrm,i,j,k) = 0.0
          omp = max(0.,min(1.,(tabs(icrm,i,j,k)-tprmin)*a_pr))
          qpl(icrm,i,j,k) = qp(icrm,i,j,k)*omp
          qpi(icrm,i,j,k) = qp(icrm,i,j,k)*(1.-omp)
        end do
      end do
    end do

  end subroutine micro_diagnose_clubb

#endif /*CLUBB_CRM*/
  !----------------------------------------------------------------------
  !!! function to compute terminal velocity for precipitating variables:
  ! In this particular case there is only one precipitating variable.

  real(crm_rknd) function term_vel_qp(i,j,k,ind,ncrms,icrm)

    use vars
    implicit none
    integer, intent(in) :: ncrms,icrm
    integer, intent(in) :: i,j,k,ind
    real(crm_rknd) wmax, omp, omg, qrr, qss, qgg

    term_vel_qp = 0.
    if(qp(icrm,i,j,k).gt.qp_threshold) then
      omp = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(icrm,i,j,k)-tprmin)*a_pr))
      if(omp.eq.1.) then
        term_vel_qp = vrain*(rho(icrm,k)*qp(icrm,i,j,k))**crain
      elseif(omp.eq.0.) then
        omg = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(icrm,i,j,k)-tgrmin)*a_gr))
        qgg=omg*qp(icrm,i,j,k)
        qss=qp(icrm,i,j,k)-qgg
        term_vel_qp = (omg*vgrau*(rho(icrm,k)*qgg)**cgrau &
        +(1.-omg)*vsnow*(rho(icrm,k)*qss)**csnow)
      else
        omg = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(icrm,i,j,k)-tgrmin)*a_gr))
        qrr=omp*qp(icrm,i,j,k)
        qss=qp(icrm,i,j,k)-qrr
        qgg=omg*qss
        qss=qss-qgg
        term_vel_qp = (omp*vrain*(rho(icrm,k)*qrr)**crain &
        +(1.-omp)*(omg*vgrau*(rho(icrm,k)*qgg)**cgrau &
        +(1.-omg)*vsnow*(rho(icrm,k)*qss)**csnow))
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

    real(crm_rknd), allocatable :: omega(:,:,:,:)
    integer ind
    integer i,j,k,icrm

    allocate(omega(ncrms,nx,ny,nzm))

    !$acc enter data create(omega) async(1)

    crain = b_rain / 4.
    csnow = b_snow / 4.
    cgrau = b_grau / 4.
    vrain = a_rain * gamr3 / 6. / (pi * rhor * nzeror) ** crain
    vsnow = a_snow * gams3 / 6. / (pi * rhos * nzeros) ** csnow
    vgrau = a_grau * gamg3 / 6. / (pi * rhog * nzerog) ** cgrau

    !$acc parallel loop gang vector collapse(4) default(present) async(1)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            omega(icrm,i,j,k) = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(icrm,i,j,k)-tprmin)*a_pr))
          end do
        end do
      end do
    end do

    call precip_fall(qp, 2, omega(:,:,:,:), ind, ncrms)

    !$acc exit data delete(omega) async(1)

    deallocate(omega)
  end subroutine micro_precip_fall

  !----------------------------------------------------------------------
  ! called when stepout() called

  subroutine micro_print()
  end subroutine micro_print

  !-----------------------------------------------------------------------
  ! Supply function that computes total water in a domain:
  !
  subroutine total_water(ncrms,tw)
    use vars, only : nstep,nprint,adz,dz,rho
    implicit none
    integer, intent(in) :: ncrms
    real(8), intent(out) :: tw(ncrms)
    real(8) tmp
    integer i,j,k,m,icrm

    !$acc parallel loop gang vector default(present) async(1)
    do icrm = 1 , ncrms
      tw(icrm) = 0.
    enddo
    do m=1,nmicro_fields
      if(flag_wmass(m).eq.1) then
        !$acc parallel loop gang vector collapse(4) default(present) async(1)
        do k=1,nzm
          do j=1,ny
            do i=1,nx
              do icrm = 1 , ncrms
                tmp = micro_field(icrm,i,j,k,m)*adz(icrm,k)*dz(icrm)*rho(icrm,k)
                !$acc atomic update
                tw(icrm) = tw(icrm) + tmp
              enddo
            end do
          end do
        end do
      end if
    end do

  end subroutine total_water

  ! -------------------------------------------------------------------------------
  ! dummy effective radius functions:

  function Get_reffc() ! liquid water
    real(crm_rknd), pointer, dimension(:,:,:) :: Get_reffc
  end function Get_reffc

  function Get_reffi() ! ice
    real(crm_rknd), pointer, dimension(:,:,:) :: Get_reffi
  end function Get_reffi



  subroutine precip_fall(qp, hydro_type, omega, ind, ncrms)
    !     positively definite monotonic advection with non-oscillatory option
    !     and gravitational sedimentation
    use vars
    use params
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd), pointer :: qp(:,:,:,:) ! falling hydrometeor
    integer hydro_type   ! 0 - all liquid, 1 - all ice, 2 - mixed
    real(crm_rknd) omega(ncrms,nx,ny,nzm)   !  = 1: liquid, = 0: ice;  = 0-1: mixed : used only when hydro_type=2
    integer ind

    ! Local:
    real(crm_rknd) wmax, omp, omg, qrr, qss, qgg
    real(crm_rknd), allocatable :: mx(:,:,:,:),mn(:,:,:,:), lfac(:,:,:,:)
    real(crm_rknd), allocatable :: www(:,:,:,:),fz(:,:,:,:)
    real(crm_rknd) eps
    integer i,j,k,kc,kb,icrm
    logical nonos

    real(crm_rknd) y,pp,pn
    pp(y)= max(real(0.,crm_rknd),y)
    pn(y)=-min(real(0.,crm_rknd),y)

    real(crm_rknd) lat_heat, term_vel_qp, flagstat, tmp
    real(crm_rknd), allocatable :: wp(:,:,:,:), tmp_qp(:,:,:,:)
    integer iprec, nprec, inttmp

    !--------------------------------------------------------
    !call t_startf ('precip_fall')

    allocate(mx    (ncrms,nx,ny,nzm))
    allocate(mn    (ncrms,nx,ny,nzm))
    allocate(lfac  (ncrms,nx,ny,nz ))
    allocate(www   (ncrms,nx,ny,nz ))
    allocate(fz    (ncrms,nx,ny,nz ))
    allocate(wp    (ncrms,nx,ny,nzm))
    allocate(tmp_qp(ncrms,nx,ny,nzm))

    !$acc enter data create(mx,mn,lfac,www,fz,wp,tmp_qp) async(1)

    eps = 1.e-10
    nonos = .true.
    nprec = 1

    ! 	Add sedimentation of precipitation field to the vert. vel.
    !$acc parallel loop gang vector collapse(4) present(lfac,wp,omega,tabs,rho,qp,dz,adz,rhow) async(1)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            if (k <= nzm) then

              kb = max(1,k-1)

              select case (hydro_type)
              case(0)
                lfac(icrm,i,j,k) = fac_cond
                flagstat = 1.
              case(1)
                lfac(icrm,i,j,k) = fac_sub
                flagstat = 1.
              case(2)
                lfac(icrm,i,j,k) = fac_cond + (1-omega(icrm,i,j,k))*fac_fus
                flagstat = 1.
              case(3)
                lfac(icrm,i,j,k) = 0.
                flagstat = 0.
              case default
                if(masterproc) then
                  print*, 'unknown hydro_type in precip_fall. exitting ...'
                  ! call task_abort
                end if
              end select

              term_vel_qp = 0.
              if(qp(icrm,i,j,k).gt.qp_threshold) then
                omp = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(icrm,i,j,k)-tprmin)*a_pr))
                if(omp.eq.1.) then
                  term_vel_qp = vrain*(rho(icrm,k)*qp(icrm,i,j,k))**crain
                elseif(omp.eq.0.) then
                  omg = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(icrm,i,j,k)-tgrmin)*a_gr))
                  qgg=omg*qp(icrm,i,j,k)
                  qss=qp(icrm,i,j,k)-qgg
                  term_vel_qp = (omg*vgrau*(rho(icrm,k)*qgg)**cgrau &
                  +(1.-omg)*vsnow*(rho(icrm,k)*qss)**csnow)
                else
                  omg = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(icrm,i,j,k)-tgrmin)*a_gr))
                  qrr=omp*qp(icrm,i,j,k)
                  qss=qp(icrm,i,j,k)-qrr
                  qgg=omg*qss
                  qss=qss-qgg
                  term_vel_qp = (omp*vrain*(rho(icrm,k)*qrr)**crain &
                  +(1.-omp)*(omg*vgrau*(rho(icrm,k)*qgg)**cgrau &
                  +(1.-omg)*vsnow*(rho(icrm,k)*qss)**csnow))
                endif
              end if

              wp(icrm,i,j,k)=sqrt(1.29/rho(icrm,k))*term_vel_qp
              tmp = wp(icrm,i,j,k)/(dz(icrm)*adz(icrm,kb)/dtn)
              if (tmp > 0.9) then
                inttmp = int(CEILING(tmp/0.9))
                !$acc atomic update
                nprec = max(nprec,inttmp)
              endif
              wp(icrm,i,j,k) = -wp(icrm,i,j,k)*rhow(icrm,k)*dtn/dz(icrm)

            elseif (k == nz) then
              lfac(icrm,i,j,nz)=0
            endif
          enddo
        enddo
      enddo
    end do  ! k

    if (nprec.gt.1) then
      !$acc parallel loop gang vector collapse(4) default(present) async(1)
      do k = 1,nzm
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              ! If maximum CFL due to precipitation velocity is greater than 0.9,
              ! take more than one advection step to maintain stability.
              ! wp already includes factor of dt, so reduce it by a
              ! factor equal to the number of precipitation steps.
              wp(icrm,i,j,k) = wp(icrm,i,j,k)/real(nprec,crm_rknd)
            end do
          enddo
        enddo
      enddo
    endif

    do iprec = 1,nprec
      !$acc parallel loop gang vector collapse(4) present(mx,mn,fz,qp,wp) async(1)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              if (k <= nzm) then
                if(nonos) then
                  kc=min(nzm,k+1)
                  kb=max(1,k-1)
                  mx(icrm,i,j,k)=max(qp(icrm,i,j,kb),qp(icrm,i,j,kc),qp(icrm,i,j,k))
                  mn(icrm,i,j,k)=min(qp(icrm,i,j,kb),qp(icrm,i,j,kc),qp(icrm,i,j,k))
                end if  ! nonos
                ! Define upwind precipitation flux
                fz(icrm,i,j,k)=qp(icrm,i,j,k)*wp(icrm,i,j,k)
              elseif (k == nz) then
                fz(icrm,i,j,nz)=0.
              endif
            enddo
          enddo
        enddo
      enddo

      !$acc parallel loop gang vector collapse(4) default(present) async(1)
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              kc=k+1
              tmp_qp(icrm,i,j,k)=qp(icrm,i,j,k)-(fz(icrm,i,j,kc)-fz(icrm,i,j,k))/(rho(icrm,k)*adz(icrm,k)) !Update temporary qp
            end do
          enddo
        enddo
      enddo

      !$acc parallel loop gang vector collapse(4) present(www,wp,rho,adz,tmp_qp) async(1)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              if (k <= nzm) then
                ! Also, compute anti-diffusive correction to previous
                ! (upwind) approximation to the flux
                kb=max(1,k-1)
                ! The precipitation velocity is a cell-centered quantity,
                ! since it is computed from the cell-centered
                ! precipitation mass fraction.  Therefore, a reformulated
                ! anti-diffusive flux is used here which accounts for
                ! this and results in reduced numerical diffusion.
                www(icrm,i,j,k) = 0.5*(1.+wp(icrm,i,j,k)/(rho(icrm,k)*adz(icrm,k))) *(tmp_qp(icrm,i,j,kb)*wp(icrm,i,j,kb) - tmp_qp(icrm,i,j,k)*wp(icrm,i,j,k)) ! works for wp(icrm,i,j,k)<0
              elseif (k == nz) then
                www(icrm,i,j,nz)=0.
              endif
            enddo
          enddo
        enddo
      end do

      !---------- non-osscilatory option ---------------

      if(nonos) then

        !$acc parallel loop gang vector collapse(4) default(present) async(1)
        do k=1,nzm
          do j=1,ny
            do i=1,nx
              do icrm = 1 , ncrms
                kc=min(nzm,k+1)
                kb=max(1,k-1)
                mx(icrm,i,j,k)=max(tmp_qp(icrm,i,j,kb),tmp_qp(icrm,i,j,kc),tmp_qp(icrm,i,j,k),mx(icrm,i,j,k))
                mn(icrm,i,j,k)=min(tmp_qp(icrm,i,j,kb),tmp_qp(icrm,i,j,kc),tmp_qp(icrm,i,j,k),mn(icrm,i,j,k))
                mx(icrm,i,j,k)=rho(icrm,k)*adz(icrm,k)*(mx(icrm,i,j,k)-tmp_qp(icrm,i,j,k))/(pn(www(icrm,i,j,kc)) + pp(www(icrm,i,j,k))+eps)
                mn(icrm,i,j,k)=rho(icrm,k)*adz(icrm,k)*(tmp_qp(icrm,i,j,k)-mn(icrm,i,j,k))/(pp(www(icrm,i,j,kc)) + pn(www(icrm,i,j,k))+eps)
              end do
            end do
          end do
        end do

        !$acc parallel loop gang vector collapse(4) default(present) async(1)
        do k=1,nzm
          do j=1,ny
            do i=1,nx
              do icrm = 1 , ncrms
                kb=max(1,k-1)
                ! Add limited flux correction to fz(icrm,i,j,k).
                fz(icrm,i,j,k) = fz(icrm,i,j,k) &                        ! Upwind flux
                + pp(www(icrm,i,j,k))*min(real(1.,crm_rknd),mx(icrm,i,j,k), mn(icrm,i,j,kb)) &
                - pn(www(icrm,i,j,k))*min(real(1.,crm_rknd),mx(icrm,i,j,kb),mn(icrm,i,j,k)) ! Anti-diffusive flux
              end do
            end do
          end do
        end do

      endif ! nonos

      !$acc parallel loop gang vector collapse(4) default(present) async(1)
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              ! Update precipitation mass fraction and liquid-ice static
              ! energy using precipitation fluxes computed in this column.
              kc=k+1
              ! Update precipitation mass fraction.
              ! Note that fz is the total flux, including both the
              ! upwind flux and the anti-diffusive correction.
              qp(icrm,i,j,k)=qp(icrm,i,j,k)-(fz(icrm,i,j,kc)-fz(icrm,i,j,k))/(rho(icrm,k)*adz(icrm,k))
              tmp = (fz(icrm,i,j,kc)-fz(icrm,i,j,k))/(rho(icrm,k)*adz(icrm,k))*flagstat
              !$acc atomic update
              qpfall(icrm,k)=qpfall(icrm,k)-tmp  ! For qp budget
              lat_heat = -(lfac(icrm,i,j,kc)*fz(icrm,i,j,kc)-lfac(icrm,i,j,k)*fz(icrm,i,j,k))/(rho(icrm,k)*adz(icrm,k))
              t(icrm,i,j,k)=t(icrm,i,j,k)-lat_heat
              !$acc atomic update
              tlat(icrm,k)=tlat(icrm,k)-lat_heat            ! For energy budget
              tmp = fz(icrm,i,j,k)*flagstat
              !$acc atomic update
              precflux(icrm,k) = precflux(icrm,k) - tmp   ! For statistics
              if (k == 1) then
                tmp = fz(icrm,i,j,1)*flagstat
                !$acc atomic update
                precsfc(icrm,i,j) = precsfc(icrm,i,j) - tmp ! For statistics
              endif
              if (k == 2) then
                tmp = fz(icrm,i,j,1)*(1.-omega(icrm,i,j,1))*flagstat
                !$acc atomic update
                precssfc(icrm,i,j) = precssfc(icrm,i,j) - tmp ! For statistics
              endif
              if (k == 3) then
                tmp = fz(icrm,i,j,1)*flagstat
                !$acc atomic update
                prec_xy(icrm,i,j) = prec_xy(icrm,i,j) - tmp ! For 2D output
              endif
            enddo
          enddo
        enddo
      end do


      !$acc parallel loop gang vector collapse(4) present(wp,fz,www,lfac,tabs,rho,qp,rhow,dz) async(1)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              if (k <= nzm) then
                ! Re-compute precipitation velocity using new value of qp.
                if (iprec.lt.nprec) then
                  term_vel_qp = 0.
                  if(qp(icrm,i,j,k).gt.qp_threshold) then
                    omp = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(icrm,i,j,k)-tprmin)*a_pr))
                    if(omp.eq.1.) then
                      term_vel_qp = vrain*(rho(icrm,k)*qp(icrm,i,j,k))**crain
                    elseif(omp.eq.0.) then
                      omg = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(icrm,i,j,k)-tgrmin)*a_gr))
                      qgg=omg*qp(icrm,i,j,k)
                      qss=qp(icrm,i,j,k)-qgg
                      term_vel_qp = (omg*vgrau*(rho(icrm,k)*qgg)**cgrau &
                      +(1.-omg)*vsnow*(rho(icrm,k)*qss)**csnow)
                    else
                      omg = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(icrm,i,j,k)-tgrmin)*a_gr))
                      qrr=omp*qp(icrm,i,j,k)
                      qss=qp(icrm,i,j,k)-qrr
                      qgg=omg*qss
                      qss=qss-qgg
                      term_vel_qp = (omp*vrain*(rho(icrm,k)*qrr)**crain &
                      +(1.-omp)*(omg*vgrau*(rho(icrm,k)*qgg)**cgrau &
                      +(1.-omg)*vsnow*(rho(icrm,k)*qss)**csnow))
                    endif
                  end if
                  wp(icrm,i,j,k) = sqrt(1.29/rho(icrm,k))*term_vel_qp
                  ! Decrease precipitation velocity by factor of nprec
                  wp(icrm,i,j,k) = -wp(icrm,i,j,k)*rhow(icrm,k)*dtn/dz(icrm)/real(nprec,crm_rknd)
                  ! Note: Don't bother checking CFL condition at each
                  ! substep since it's unlikely that the CFL will
                  ! increase very much between substeps when using
                  ! monotonic advection schemes.
                endif
              elseif (k == nz) then
                if (iprec.lt.nprec) then
                  fz(icrm,i,j,nz)=0.
                  www(icrm,i,j,nz)=0.
                  lfac(icrm,i,j,nz)=0.
                endif
              endif
            end do
          end do
        end do
      end do
    end do !iprec

    !$acc exit data delete(mx,mn,lfac,www,fz,wp,tmp_qp) async(1)

    deallocate(mx    )
    deallocate(mn    )
    deallocate(lfac  )
    deallocate(www   )
    deallocate(fz    )
    deallocate(wp    )
    deallocate(tmp_qp)

    !call t_stopf ('precip_fall')

  end subroutine precip_fall


end module microphysics
