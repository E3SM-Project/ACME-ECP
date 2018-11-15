! -----------------------------------------------------------------------------
! MODULE  accelerate_crm_mod
!    This module provides functionality to apply mean-state acceleration
!    (Jones et al., 2015) to CRM. ... (to be continued)
!
! author: Christopher R Jones
! email: christopher.jones@pnnl.gov
! date: 1/2018
!
!  Contains subroutines:
!    crm_accel_nstop() - adjusts nstop based on crm_accel_factor
!    accelerate_crm()  - applies mean state acceleration tendency to crm fields
!
! updated: 11/2018
!    need to account for additional "ncrms" dimension to most variables
! -----------------------------------------------------------------------------
module accelerate_crm_mod
#ifdef UNITTEST
  use unittest_mod
#else
  use grid, only: nx, ny
  use shr_kind_mod, only: r8=>shr_kind_r8
  use params, only: rc=>crm_rknd
#endif

  implicit none

  ! variables
  real(r8) :: coef = 1._r8 / dble(nx * ny)  ! coefficient for horizontal avg
  ! specify this stuff I need ...
  logical   :: use_crm_accel
  logical   :: crm_accel_uv
  logical   :: distribute_qneg
  real(r8)  :: crm_accel_factor

#ifdef UNITTEST
  public :: coef
  ! subroutines
  public :: accelerate_scalars
  public :: accelerate_momentum
  public :: accelerate_t
  public :: accelerate_micro
  public :: crm_horiz_mean
  public :: apply_accel_tend_micro
  public :: partition_micro
#else
  private :: coef, distribute_qneg, crm_accel_uv
  ! subroutines
  private :: accelerate_scalars
  private :: accelerate_momentum
  private :: accelerate_t
  private :: accelerate_micro
  private :: crm_horiz_mean
  private :: apply_accel_tend_micro
  private :: partition_micro
#endif
  public :: use_crm_accel, crm_accel_factor
  public :: accelerate_crm
  public :: crm_accel_nstop
  public :: crm_accel_reset_nstop
  public :: crm_accel_init
  public :: crm_accel_verbose_debug
  public :: accelerate_crm_orig
contains

  subroutine crm_accel_init()
    use phys_control, only: phys_getopts
    use cam_logfile, only: iulog
    use spmd_utils,  only: masterproc

    implicit none
    integer :: crm_accel_micro_opt = 0

    ! initialize defaults
    use_crm_accel = .false.
    crm_accel_factor = 0.
    crm_accel_uv = .false.

#ifdef CRMACCEL
    call phys_getopts(use_crm_accel_out = use_crm_accel, &
                      crm_accel_factor_out = crm_accel_factor, &
                      crm_accel_uv_out = crm_accel_uv, &
                      crm_accel_micro_opt_out = crm_accel_micro_opt)
#endif
    if (crm_accel_micro_opt .eq. 1) then
      distribute_qneg = .true.
    else
      distribute_qneg = .false.
    endif

    if (masterproc) then
       write(iulog, *) 'USING CRM MEAN STATE ACCELERATION'
       write(iulog, *) 'crm_accel: use_crm_accel = ', use_crm_accel
       write(iulog, *) 'crm_accel: crm_accel_factor = ', crm_accel_factor
       write(iulog, *) 'crm_accel: crm_accel_uv = ', crm_accel_uv
       write(iulog, *) 'crm_accel: crm_accel_micro_opt = ', crm_accel_micro_opt
       write(iulog, *) 'crm_accel: setting distribute_qneg = ', distribute_qneg
    end if
  end subroutine crm_accel_init

  subroutine crm_accel_nstop(nstop)
  ! implementation option (1): use grid, only: nstop
  ! option (2): pass nstop in as inout variable (that's my preference)
  ! note: should write a unit test to make sure it (a) evenly divides when it
  !       should and (b) quits when it can't
#ifdef UNITTEST
    use unittest_mod
#else
    use cam_abortutils, only: endrun
    use cam_logfile,  only: iulog
#endif

    implicit none

    integer, intent(inout) :: nstop

    if (mod(real(nstop), (1. + crm_accel_factor)) .ne. 0) then
      write(iulog,*) "CRM acceleration unexpected exception:"
      write(iulog,*) "(1+crm_accel_factor) does not divide equally into nstop"
      write(iulog,*) "nstop = ", nstop
      write(iulog,*) "crm_accel_factor = ", crm_accel_factor
      call endrun('crm main: bad crm_accel_factor and nstop pair')
      nstop = -1  ! hack for testing
    else
      nstop = nstop / (1.D0 + crm_accel_factor)
    end if
  end subroutine crm_accel_nstop


  subroutine crm_accel_reset_nstop(nstop, nstep)
    ! If acceleration was turned off for remainder of steps at step "nstep",
    ! need to adjust nstop to account for remaining steps advancing at
    ! original (unaccelerated) pace.
#ifdef UNITTEST
    use unittest_mod, only: iulog
#else
    use cam_abortutils, only: endrun
    use cam_logfile,  only: iulog
#endif
    implicit none

    integer, intent(inout) :: nstop
    integer, intent(in) :: nstep
    write (iulog,*) 'crm: nstop increased from ', nstop, ' to ', &
            int(nstop+(nstop-nstep+1)*crm_accel_factor)
    nstop = nstop + (nstop - nstep + 1)*crm_accel_factor ! only can happen once
  end subroutine crm_accel_reset_nstop


  subroutine accelerate_crm(crm_accel_ceaseflag, icrm)
    implicit none

    logical, intent(out) :: crm_accel_ceaseflag
    integer, intent(in) :: icrm

    ! namelist variable: crm_accel_uv, crm_accel_micro_opt

    ! accelerate scalars t and q (and micro_field(:,:,:, index_water_vapor, icrm))
    ! raise crm_accel_ceaseflag and cancel mean-state acceleration
    !   if magnitude of t-tend is too great
    call accelerate_scalars(crm_accel_ceaseflag, icrm)

    ! accelerate velocity u, v
    if (crm_accel_uv .and. .not. crm_accel_ceaseflag) then
      call accelerate_momentum(icrm)
    end if
  end subroutine accelerate_crm

  subroutine accelerate_scalars(ceaseflag, icrm)
    ! accelerates the scalar fields (t and q)
    implicit none

    logical, intent(out) :: ceaseflag
    integer, intent(in) :: icrm

    ! initializations
    ceaseflag = .false.    ! flag to stop applying accelerat_crm

    call accelerate_t(ceaseflag, icrm)
    if (.not. ceaseflag) then
      call accelerate_micro(icrm)
    end if
  end subroutine accelerate_scalars

  subroutine accelerate_t(ceaseflag, icrm)
    ! accelerates liquid static energy (t)
#ifdef UNITTEST
    use unittest_mod
#else
    use grid, only: nx, ny, nzm
    use vars, only: t, t0
    use cam_logfile,  only: iulog
#endif

    implicit none
    logical, intent (out) :: ceaseflag
    integer, intent(in) :: icrm

    ! local variables
    integer i, j, k
    real(rc) :: tbaccel(nzm)   ! t before acceleration
    real(rc) :: ttend_acc(nzm) ! mean-state acceleration tendency

    ! initializations
    ceaseflag = .false.    ! flag to stop applying accelerate_crm

    ! calculate tendency * dtn
    call crm_horiz_mean(tbaccel, t(1:nx, 1:ny, :, icrm))
    ttend_acc(1:nzm) = tbaccel(1:nzm) - t0(1:nzm, icrm)

    ! stop accelerating if acceleration tendency is too large anywhere
    if(maxval(abs(ttend_acc)) .gt. 5.) then ! special clause for dT/dt too large
      ceaseflag = .true.
      write (iulog, *) 'accelerate_crm: |dT|>5K; dT = ', ttend_acc  ! write full array
      write (iulog, *) 'mean-state acceleration not applied this step'
    else
      ttend_acc = crm_accel_factor * ttend_acc
      do k = 1, nzm
        do j = 1, ny
          do i = 1, nx
            ! don't let T go negative!
            t(i, j, k, icrm) = max(50._rc, t(i, j, k, icrm) + ttend_acc(k))
          end do
        end do
      end do
    end if
  end subroutine accelerate_t

  subroutine accelerate_micro(icrm)
#ifdef UNITTEST
    use unittest_mod
#else
    use grid, only: nx, ny, nzm
    use vars, only: qcl, qci, qv, q0
    use microphysics, only: micro_field, ixw=>index_water_vapor
    use cam_logfile,  only: iulog
#endif
    implicit none
    integer, intent(in) :: icrm

    ! local variables
    integer i, j, k
    real(rc) :: qtbaccel(nzm) ! t and q before acceleration
    real(rc) :: qtend_acc(nzm), neg_qacc(nzm) ! accel tendencies

    do k = 1, nzm
      ! crjones: better to work with micro_field direcly?
      qtbaccel(k) = sum(qcl(1:nx, 1:ny, k, icrm) + qci(1:nx, 1:ny, k, icrm) + qv(1:nx, 1:ny, k, icrm)) * coef
      ! neg_qacc(k) = 0.
    end do

    ! tendency * dtn
    qtend_acc(1:nzm) = qtbaccel(1:nzm) - q0(1:nzm, icrm)
    qtend_acc = crm_accel_factor * qtend_acc

    if(distribute_qneg) then
      ! redistribute moistre in level by removing from cells with positive q
      call apply_accel_tend_micro(qtend_acc, icrm)
    else
      ! original version: simply remove negative moisture
      ! crjones concern: since qcl and qci are not touched, this could lead to
      !   inconsistent acceleration tendency, since q0 at start of next step is
      !   diagnosed from qv + qcl + qci after micro_init, but saturation
      !   adjustment is not applied until micro_proc (so it's possible that
      !   micro_field(:,:,k,index_water_vapor) = 0, but q0(k) > 0, leading to
      !   incorrect acceleration tendency that may even have incorrect sign).
      do k = 1, nzm
        do j = 1, ny
          do i = 1, nx
            micro_field(i,j,k,ixw, icrm) = &
              micro_field(i,j,k,ixw, icrm)+qtend_acc(k)
            ! enforce positivity and accumulate (negative) excess
            if(micro_field(i,j,k,ixw, icrm) .lt. 0.) then
              ! neg_qacc(k)=neg_qacc(k)+micro_field(i,j,k)
              micro_field(i,j,k,ixw, icrm)=0.
            end if
            qv(i, j, k, icrm) = max(0._rc, qv(i, j, k, icrm) + qtend_acc(k))
          end do
        end do
      end do
    end if
  end subroutine accelerate_micro

  subroutine crm_horiz_mean(out1D, in3D)
#ifdef UNITTEST
    use unittest_mod
#else
    use grid, only: nx, ny, nzm
#endif

    implicit none

    real(rc), intent(in) :: in3D(:, :, :)
    real(rc), intent(out) :: out1D(nzm)

    integer k
    do k = 1, nzm
      out1D(k) = sum(in3D(1:nx, 1:ny, k)) * coef
    end do
  end subroutine crm_horiz_mean

  subroutine apply_accel_tend_micro(deltaq, icrm)
#ifdef UNITTEST
    use unittest_mod
#else
    use vars, only: qv, qcl, qci
    use microphysics, only: micro_field, ixw=>index_water_vapor
    use grid, only: nx, ny, nzm
#endif
    implicit none
    real(rc), intent(in) :: deltaq(nzm)
    integer, intent(in) :: icrm

    integer i, j, k, nneg
    real(r8) :: qpoz, qneg, factor

    do k = 1, nzm
      qpoz = 0.
      qneg = 0.

      ! update micro_field
      micro_field(1:nx, 1:ny, k, ixw, icrm) = &
        micro_field(1:nx, 1:ny, k, ixw, icrm) + deltaq(k)

      if (deltaq(k) .ge. 0.) then
        ! skip the hole-filling logic, dump all deltaq(k) all into qv
        qv(1:nx, 1:ny, k, icrm) = qv(1:nx, 1:ny, k, icrm) + deltaq(k)
        cycle
      end if

      ! accumulate negative excess (if any)
      do j = 1, ny
        do i = 1, nx
          if (micro_field(i, j, k, ixw, icrm) .lt. 0.) then
            qneg = qneg + micro_field(i, j, k, ixw, icrm)
          else
            qpoz = qpoz + micro_field(i, j, k, ixw, icrm)
          end if
        end do
      end do

      ! warn if all q is removed and set appropriately.
      ! note to self: can possibly improve performance by adding loops in future
      if (qpoz + qneg .lt. 0.) then
        write(*, *) "crm_accel_warning: all moisture depleted in layer ",k
        micro_field(1:nx, 1:ny, k, ixw, icrm) = 0._rc
        qv(1:nx, 1:ny, k, icrm) = 0._rc
        qcl(1:nx, 1:ny, k, icrm) = 0._rc
        qci(1:nx, 1:ny, k, icrm) = 0._rc
      else
        factor = 1._r8 + qneg / qpoz
        ! apply to micro_field and partition qv, qcl, qci appropriately
        micro_field(1:nx, 1:ny, k, ixw, icrm) = &
          max(0._rc, micro_field(1:nx, 1:ny, k, ixw, icrm) * factor)
        call partition_micro(qv(1:nx, 1:ny, k, icrm), qcl(1:nx, 1:ny, k, icrm), qci(1:nx, 1:ny, k, icrm), &
                             micro_field(1:nx, 1:ny, k, ixw, icrm))
      end if
    end do  ! k = 1, nzm
  end subroutine apply_accel_tend_micro

  elemental subroutine partition_micro(qvx, qclx, qcix, micro_fieldx)
    ! partition micro_field into qv, qcl, and qci
    ! such that micro_field = qv + qcl + qci, following these rules:
    !  (1) adjust qv first
    !  (2) adjust qcl and qci only if needed to ensure positivity
#ifdef UNITTEST
    use unittest_mod, only: rc
#endif
    implicit none

    real(rc), intent(out) :: qvx
    real(rc), intent(inout) :: qclx, qcix
    real(rc), intent(in) :: micro_fieldx

    real(r8) :: qfactor

    if (micro_fieldx .le. 0._rc) then
      qvx = 0.
      qclx = 0.
      qcix = 0.
    else
      qfactor = micro_fieldx - qclx - qcix
      qvx = max(0._rc, qfactor)
      if (qfactor .lt. 0._r8) then
        ! note: this theoretically is guaranteed to work
        qfactor = 1._r8 + qfactor / (qclx + qcix)
        qclx = qclx * qfactor
        qcix = qcix * qfactor
      end if
    end if
  end subroutine partition_micro

  subroutine accelerate_momentum(icrm)
  ! accelerates the velocity fields (u and v)
#ifdef UNITTEST
    use unittest_mod
#else
    use domain, only: yes3d
    use grid, only: nzm
    use vars, only: u, v, u0, v0
#endif

  implicit none
  integer, intent(in) :: icrm

  ! local variables
  integer k
  real(rc) :: ubaccel(nzm), vbaccel(nzm)      ! u and v before acceleration
  real(rc) :: utend_acc(nzm), vtend_acc(nzm)  ! accel tendencies

  ! always accelerate u
  call crm_horiz_mean(ubaccel, u(1:nx, 1:ny, :, icrm))
  utend_acc(1:nzm) = crm_accel_factor * (ubaccel(1:nzm) - u0(1:nzm, icrm))
  do k = 1, nzm
    u(1:nx, 1:ny, k, icrm) = u(1:nx, 1:ny, k, icrm) + utend_acc(k)
  end do

  ! only mess with v if 3D:
  if (yes3d .gt. 0) then
    call crm_horiz_mean(vbaccel, v(1:nx, 1:ny, :, icrm))
    vtend_acc(1:nzm) = crm_accel_factor * (vbaccel(1:nzm) - v0(1:nzm, icrm))
    do k = 1, nzm
      v(1:nx, 1:ny, k, icrm) = v(1:nx, 1:ny, k, icrm) + vtend_acc(k)
    end do
  endif
end subroutine accelerate_momentum

subroutine crm_accel_verbose_debug(icrm)
  ! verbose log output (kinda crazy)
  use vars
  use microphysics, only: micro_field, ixw=>index_water_vapor
  use grid, only: nzm

  implicit none
  integer, intent(in) :: icrm

  integer k
  do k=1,nzm
    write (0,*) '(debug) accel_skipping: icrm = ', icrm
    write (0,*) '(debug) accel_skipping: lev k =', k
    write (0,*) '(debug) accel_skipping: t = ', t(1:nx, 1:ny, k, icrm)
    write (0,*) '(debug) accel_skipping: t0 = ', t0(k, icrm)
    write (0,*) '(debug) accel_skipping: qv = ', qv(1:nx, 1:ny, k, icrm)
    write (0,*) '(debug) accel_skipping: qcl = ', qcl(1:nx, 1:ny, k, icrm)
    write (0,*) '(debug) accel_skipping: qci = ', qci(1:nx, 1:ny, k, icrm)
    write (0,*) '(debug) accel_skipping: q0 = ', q0(k, icrm)
    write (0,*) '(debug) accel_skipping: micro_field(:,:,k,1) = ', micro_field(1:nx, 1:ny, k, ixw, icrm)
    write (0,*) '(debug) accel_skipping: u = ', u(1:nx, 1:ny, k, icrm)
    write (0,*) '(debug) accel_skipping: u0 = ', u0(k, icrm)
  end do
end subroutine crm_accel_verbose_debug

subroutine accelerate_crm_orig(ceaseflag, icrm)
! original version of accelerate_crm from UPCAM port ...
  use shr_kind_mod, only: r8=>shr_kind_r8
  use vars
  use params
  use microphysics, only: micro_field, index_water_vapor
  implicit none

  logical, intent (out):: ceaseflag
  integer, intent(in) :: icrm

  integer i,j,k
  real(r8) :: coefxy
  real(r8) :: dq_accel,tbaccel(nzm),qtbaccel(nzm)
  real(r8) :: ttend_acc(nzm), qtend_acc(nzm), neg_qacc(nzm)
  real(r8) :: accel_factor

  accel_factor = crm_accel_factor
  coefxy = 1./dble(nx*ny)

  ! NOTE:
  ! neg_qacc(k) now equals horizontal mean
  ! qt that was not removed from the system because acceleration
  ! tendency would have driven qt negative.

  ! calculate horizontal means
  do k=1,nzm
    tbaccel(k)=0.    ! lse after physics (before applying accel)
    qtbaccel(k)=0.   ! qt after physics
    neg_qacc(k) = 0. ! excess qt that cannot be depleted
    do j=1,ny
      do i=1,nx
       tbaccel(k) = tbaccel(k)+t(i,j,k,icrm)
       qtbaccel(k) = qtbaccel(k) + qcl(i,j,k,icrm)+qci(i,j,k,icrm)+qv(i,j,k,icrm)
      end do
    end do
    tbaccel(k)=tbaccel(k)*coefxy
    qtbaccel(k)=qtbaccel(k)*coefxy
  end do ! k

  ceaseflag = .false.
  do k=1,nzm
    do j=1,ny
        do i=1,nx
           if (abs(tbaccel(k)-t0(k,icrm)) .gt. 5.) then ! special clause for cases when dTdt is too large
              ceaseflag = .true.
  ! Note host crm_module receives this and adjusts number of integration steps accordingly
              write (6,*) 'MDEBUG: |dT|>5K; dT,i,k,icrm=',tbaccel(k)-t0(k,icrm),i,k,icrm
           endif
        end do
    end do
  end do

  if (.not. ceaseflag) then
  ! apply acceleration tendency
  do k=1,nzm
       ! pritch notes t0 and q0 are profiles of horizontal average field
          ! available in common.inc

             ! pritch asks  - what is dtn?
  !              dtn = dt/ncycle (from crm.F)
  !              dynamically adjusted timestep, modified based on
  !              convergence issues

             ! pritch asks - what is t0 and when is it updated?
  !               t0,q0 = mean domain profiles prior to CRM time
  !               integration loop.

     ttend_acc(k) = accel_factor*(tbaccel(k)-t0(k,icrm))/dtn
     dq_accel = accel_factor*(qtbaccel(k) - q0(k,icrm))
     qtend_acc(k) = dq_accel/dtn
     do j=1,ny
        do i=1,nx
           t(i,j,k,icrm) = max(50.,t(i,j,k,icrm) +accel_factor*(tbaccel(k)-t0(k,icrm))) 
           ! pritch, avoid abs T going negative in cases of extreme horizontal mean temperature change
           micro_field(i,j,k,index_water_vapor,icrm) = &
                micro_field(i,j,k,index_water_vapor,icrm)+dq_accel

           ! enforce positivity and accumulate (negative) excess
           if(micro_field(i,j,k,index_water_vapor,icrm) .lt. 0.) then
              neg_qacc(k)=neg_qacc(k)+micro_field(i,j,k,index_water_vapor,icrm)
              micro_field(i,j,k,index_water_vapor,icrm)=0.
           end if

           ! add qt tendency to qv
           qv(i,j,k,icrm) = max(0.,qv(i,j,k,icrm)+dq_accel)   !hparish: the current calcultion of qv is not consistent with line 1185 of
  !                                                 microphysics.f90. because dq_qccel represents the total water and not the water vapor.
  !                                                 the possible reason is unknown. The introduction of parallel variables is questionble.
  !                                                 like micro_field and qv, qcl, etc. which requires clean up. The variable doubling
  !                                                 will also jeopardize the performance. date: 8/31/2015.
        end do
     end do
     neg_qacc(k) = neg_qacc(k)*coefxy
  end do
  endif
end subroutine accelerate_crm_orig

end module accelerate_crm_mod
