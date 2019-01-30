! -----------------------------------------------------------------------------
! MODULE  accelerate_crm_mod
!    This module provides functionality to apply mean-state acceleration
!    (Jones et al., 2015) to CRM. ... (to be continued)
!
! author: Christopher R Jones
! email: christopher.jones@pnnl.gov
! date: 1/2019
!
!  Contains subroutines:
!    crm_accel_init() - initialize quantities needed to apply acceleration
!    crm_accel_nstop() - adjusts nstop based on crm_accel_factor
!    accelerate_crm()  - applies mean state acceleration tendency to crm fields
!
! updated: 11/2018
!    need to account for additional "ncrms" dimension to most variables
! -----------------------------------------------------------------------------
module accelerate_crm_mod
    use grid, only: nx, ny
    use shr_kind_mod, only: r8=>shr_kind_r8
    use params, only: rc=>crm_rknd
  
    implicit none
  
    ! variables
    real(r8) :: coef = 1._r8 / dble(nx * ny)  ! coefficient for horizontal avg
    ! specify this stuff I need ...
    logical   :: use_crm_accel
    logical   :: crm_accel_uv
    logical   :: distribute_qneg
    real(r8)  :: crm_accel_factor
  
    ! private constants
    private :: coef, distribute_qneg, crm_accel_uv
    
    ! private subroutines
    private :: accelerate_scalars
    private :: accelerate_momentum
    private :: accelerate_t
    private :: accelerate_micro
    private :: crm_horiz_mean
    private :: apply_accel_tend_micro
    private :: partition_micro

    ! public subroutines
    public :: use_crm_accel, crm_accel_factor
    public :: accelerate_crm
    public :: crm_accel_nstop
    public :: crm_accel_reset_nstop
    public :: crm_accel_init
  contains

    subroutine crm_accel_init()
    ! initialize namelist options for CRM mean-state acceleration
      use phys_control, only: phys_getopts
      use cam_logfile, only: iulog
      use spmd_utils,  only: masterproc
      use cam_abortutils, only: endrun
  
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
#if defined(CRMACCEL) && !defined(sam1mom)
      ! ensure CRMACCEL runs with sam1mom only
      if (masterproc) then
        write(0,*) "CRMACCEL is only compatible with sam1mom microphysics"
        call endrun('crm main')
      endif
#endif
    
    end subroutine crm_accel_init
  
    subroutine crm_accel_nstop(nstop)
    ! reduce nstop to appropriate value give crm_accel_factor
      use cam_abortutils, only: endrun
      use cam_logfile,  only: iulog
  
      implicit none
  
      integer, intent(inout) :: nstop
  
      if (mod(real(nstop), (1. + crm_accel_factor)) .ne. 0) then
        write(iulog,*) "CRM acceleration unexpected exception:"
        write(iulog,*) "(1+crm_accel_factor) does not divide equally into nstop"
        write(iulog,*) "nstop = ", nstop
        write(iulog,*) "crm_accel_factor = ", crm_accel_factor
        call endrun('crm main: bad crm_accel_factor and nstop pair')
      else
        nstop = nstop / (1.D0 + crm_accel_factor)
      end if
    end subroutine crm_accel_nstop


    subroutine crm_accel_reset_nstop(nstop, nstep)
      ! If acceleration was turned off for remainder of steps at step "nstep",
      ! need to adjust nstop to account for remaining steps advancing at
      ! original (unaccelerated) pace.
      use cam_logfile,  only: iulog
  
      implicit none
  
      integer, intent(inout) :: nstop
      integer, intent(in) :: nstep
      write (iulog,*) 'crm: nstop increased from ', nstop, ' to ', &
              int(nstop+(nstop-nstep+1)*crm_accel_factor)
      nstop = nstop + (nstop - nstep + 1)*crm_accel_factor ! only can happen once
    end subroutine crm_accel_reset_nstop


    subroutine accelerate_crm(ncrms, crm_accel_ceaseflag)
      ! accelerate scalars t and q (and micro_field(:,:,:, index_water_vapor, icrm))
      ! raise crm_accel_ceaseflag and cancel mean-state acceleration
      !       if magnitude of t-tendency is too great
      implicit none
  
      integer, intent(in) :: ncrms
      logical, intent(out) :: crm_accel_ceaseflag    
      call accelerate_scalars(ncrms, crm_accel_ceaseflag)
  
      ! accelerate velocity u, v
      if (crm_accel_uv .and. .not. crm_accel_ceaseflag) then
        call accelerate_momentum(ncrms)
      end if
    end subroutine accelerate_crm
  
    subroutine accelerate_scalars(ncrms, ceaseflag)
      ! accelerates the scalar fields (t and q)
      ! aborts and returns ceaseflag = .true. if t-tendency exceeds threshhold anywhere
      implicit none
  
      integer, intent(in) :: ncrms
      logical, intent(out) :: ceaseflag
  
      ! initializations
      ceaseflag = .false.    ! flag to stop applying accelerate_crm
  
      call accelerate_t(ncrms, ceaseflag)
      if (.not. ceaseflag) then
        call accelerate_micro(ncrms)
      end if
    end subroutine accelerate_scalars
  
    subroutine accelerate_t(ncrms, ceaseflag)
      ! accelerates liquid static energy (t)
      use grid, only: nx, ny, nzm
      use vars, only: t, t0
      use cam_logfile,  only: iulog
  
      implicit none
      integer, intent(in) :: ncrms
      logical, intent(out) :: ceaseflag
  
      ! local variables
      integer i, j, k, icrm
      real(rc) :: tbaccel(nzm, ncrms)   ! t before acceleration
      real(rc) :: ttend_acc(nzm, ncrms) ! mean-state acceleration tendency
  
      ! initializations
      ceaseflag = .false.    ! flag to stop applying accelerate_crm
  
      do icrm = 1, ncrms
        ! calculate tendency * dtn
        call crm_horiz_mean(tbaccel(:, icrm), t(1:nx, 1:ny, :, icrm))
        ttend_acc(1:nzm, icrm) = tbaccel(1:nzm, icrm) - t0(1:nzm, icrm)
      end do

      ! stop accelerating if acceleration tendency is too large anywhere
      if(maxval(abs(ttend_acc)) .gt. 5.) then ! special clause for dT/dt too large
        ceaseflag = .true.
        write (iulog, *) 'accelerate_crm: |dT|>5K; dT = ', ttend_acc  ! write full array
        write (iulog, *) 'mean-state acceleration not applied this step'
      else
        do icrm = 1, ncrms
          ttend_acc(1:nzm, icrm) = crm_accel_factor * ttend_acc(1:nzm, icrm)
          do k = 1, nzm
            do j = 1, ny
              do i = 1, nx
                ! don't let T go negative!
                t(i, j, k, icrm) = max(50._rc, t(i, j, k, icrm) + ttend_acc(k, icrm))
              end do
            end do
          end do
        end do
      end if
    end subroutine accelerate_t
  
    subroutine accelerate_micro(ncrms)
      use grid, only: nx, ny, nzm
      use vars, only: qcl, qci, qv, q0
      use microphysics, only: micro_field, ixw=>index_water_vapor
      use cam_logfile,  only: iulog

      implicit none
      integer, intent(in) :: ncrms
  
      ! local variables
      integer i, j, k, icrm
      real(rc) :: qtbaccel(nzm, ncrms)  ! t and q before acceleration
      real(rc) :: qtend_acc(nzm, ncrms) ! accel tendencies

      do icrm = 1, ncrms
        do k = 1, nzm
          ! crjones: better to work with micro_field direcly?
          qtbaccel(k, icrm) = sum(qcl(1:nx, 1:ny, k, icrm) + qci(1:nx, 1:ny, k, icrm) + qv(1:nx, 1:ny, k, icrm)) * coef
        end do
      end do
  
      do icrm = 1, ncrms
        ! tendency * dtn
        qtend_acc(1:nzm, icrm) = qtbaccel(1:nzm, icrm) - q0(1:nzm, icrm)
        qtend_acc(1:nzm, icrm) = crm_accel_factor * qtend_acc(1:nzm, icrm)
      end do
  
      if(distribute_qneg) then
        ! redistribute moistre in level by removing from cells with positive q
        call apply_accel_tend_micro(ncrms, qtend_acc)
      else
        ! original version: simply remove negative moisture
        ! crjones concern: since qcl and qci are not touched, this could lead to
        !   inconsistent acceleration tendency, since q0 at start of next step is
        !   diagnosed from qv + qcl + qci after micro_init, but saturation
        !   adjustment is not applied until micro_proc (so it's possible that
        !   micro_field(:,:,k,index_water_vapor) = 0, but q0(k) > 0, leading to
        !   incorrect acceleration tendency that may even have incorrect sign).
        do icrm = 1, ncrms
          do k = 1, nzm
            do j = 1, ny
              do i = 1, nx
                micro_field(i,j,k,ixw, icrm) = &
                  micro_field(i,j,k,ixw, icrm)+qtend_acc(k, icrm)
                ! enforce positivity and accumulate (negative) excess
                if(micro_field(i,j,k,ixw, icrm) .lt. 0.) then
                  micro_field(i,j,k,ixw, icrm)=0.
                end if
                qv(i, j, k, icrm) = max(0._rc, qv(i, j, k, icrm) + qtend_acc(k, icrm))
              end do
            end do
          end do
        end do
      end if
    end subroutine accelerate_micro
  
    subroutine crm_horiz_mean(out1D, in3D)
      ! returns horizontal mean of in3D(1:nx, 1:ny, 1:nzm) as out1D
      use grid, only: nx, ny, nzm
  
      implicit none
  
      real(rc), intent(in) :: in3D(:, :, :)
      real(rc), intent(out) :: out1D(nzm)
  
      integer k
      do k = 1, nzm
        out1D(k) = sum(in3D(1:nx, 1:ny, k)) * coef
      end do
    end subroutine crm_horiz_mean
  

    subroutine apply_accel_tend_micro(ncrms, deltaq)
      use vars, only: qv, qcl, qci
      use microphysics, only: micro_field, ixw=>index_water_vapor
      use grid, only: nx, ny, nzm
      implicit none
      integer, intent(in) :: ncrms
      real(rc), intent(in) :: deltaq(nzm, ncrms)
  
      integer i, j, k, nneg, icrm
      real(r8) :: qpoz, qneg, factor

      do icrm = 1, ncrms
        do k = 1, nzm
          qpoz = 0.
          qneg = 0.
  
          ! update micro_field
          micro_field(1:nx, 1:ny, k, ixw, icrm) = &
            micro_field(1:nx, 1:ny, k, ixw, icrm) + deltaq(k, icrm)
  
          ! crjones note: this should probably change with gpu acceleration
          if (deltaq(k, icrm) .ge. 0.) then
            ! skip the hole-filling logic, dump all deltaq(k) all into qv
            qv(1:nx, 1:ny, k, icrm) = qv(1:nx, 1:ny, k, icrm) + deltaq(k, icrm)
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
  
          if (qpoz + qneg .lt. 0.) then
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
      end do ! icrm = 1, ncrms
    end subroutine apply_accel_tend_micro
  
    elemental subroutine partition_micro(qvx, qclx, qcix, micro_fieldx)
      ! partition micro_field into qv, qcl, and qci
      ! such that micro_field = qv + qcl + qci, following these rules:
      !  (1) adjust qv first
      !  (2) adjust qcl and qci only if needed to ensure positivity
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
  
    subroutine accelerate_momentum(ncrms)
    ! accelerates the velocity fields (u and v)
      use domain, only: yes3d
      use grid, only: nzm
      use vars, only: u, v, u0, v0
  
    implicit none
    integer, intent(in) :: ncrms
  
    ! local variables
    integer k, icrm
    real(rc) :: ubaccel(nzm, ncrms), vbaccel(nzm, ncrms)      ! u and v before acceleration
    real(rc) :: utend_acc(nzm, ncrms), vtend_acc(nzm, ncrms)  ! accel tendencies
  
    do icrm = 1, ncrms
      ! always accelerate u
      call crm_horiz_mean(ubaccel(1:nzm, icrm), u(1:nx, 1:ny, :, icrm))
      utend_acc(1:nzm, icrm) = crm_accel_factor * (ubaccel(1:nzm, icrm) - u0(1:nzm, icrm))
      do k = 1, nzm
        u(1:nx, 1:ny, k, icrm) = u(1:nx, 1:ny, k, icrm) + utend_acc(k, icrm)
      end do
    
      ! only mess with v if 3D:
      if (yes3d .gt. 0) then
        call crm_horiz_mean(vbaccel(1:nzm, icrm), v(1:nx, 1:ny, :, icrm))
        vtend_acc(1:nzm, icrm) = crm_accel_factor * (vbaccel(1:nzm, icrm) - v0(1:nzm, icrm))
        do k = 1, nzm
          v(1:nx, 1:ny, k, icrm) = v(1:nx, 1:ny, k, icrm) + vtend_acc(k, icrm)
        end do
      endif
    end do
  end subroutine accelerate_momentum
    
end module accelerate_crm_mod
  
