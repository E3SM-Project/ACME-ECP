! -----------------------------------------------------------------------------
! MODULE  accelerate_crm_mod
!   This module provides functionality to apply mean-state acceleration (MSA)
!   (Jones et al., 2015, doi:10.1002/2015MS000488) to the CRMs when using
!   superparameterization.
!
! PUBLIC SUBROUTINES:
!   crm_accel_init: initialize quantities needed to apply MSA
!   crm_accel_nstop: adjusts 'nstop' in crm_module based on crm_accel_factor
!   accelerate_crm: calculates and applies MSA tendency to CRM
!
! PUBLIC MODULE VARIABLES:
!   logical  :: use_crm_accel - apply MSA if true (cam namelist variable)
!   real(r8) :: crm_accel_factor - MSA factor to use (cam namelist variable)
!
! REVISION HISTORY:
!   2018-Nov-01: Initial implementation
!   2019-Jan-30: Initial subroutine port to GPU using openacc directives
!
! CONTACT: Christopher Jones (christopher.jones@pnnl.gov)
! -----------------------------------------------------------------------------
module accelerate_crm_mod
    use grid, only: nx, ny
    use shr_kind_mod, only: r8=>shr_kind_r8
    use params, only: asyncid, rc=>crm_rknd

    implicit none

    ! private module variables
    real(r8), parameter :: coef = 1._r8 / dble(nx * ny)  ! coefficient for horizontal averaging
    logical :: crm_accel_uv  ! (false) apply MSA only to scalar fields (T and QT)
                             ! (true) apply MSA to winds (U/V) and scalar fields

    ! public module variables
    logical :: use_crm_accel  ! use MSA if true
    real(r8) :: crm_accel_factor  ! 1 + crm_accel_factor = 'a' in Jones etal (2015)

    private :: coef, crm_accel_uv
    public :: use_crm_accel, crm_accel_factor
    public :: accelerate_crm
    public :: crm_accel_nstop
    public :: crm_accel_init

  contains
    subroutine crm_accel_init()
    ! initialize namelist options for CRM mean-state acceleration
      use phys_control, only: phys_getopts
      use cam_logfile, only: iulog
      use spmd_utils,  only: masterproc
      use cam_abortutils, only: endrun
  
      implicit none
  
      ! initialize defaults
      use_crm_accel = .false.
      crm_accel_factor = 0.
      crm_accel_uv = .false.
  
      call phys_getopts(use_crm_accel_out = use_crm_accel, &
                        crm_accel_factor_out = crm_accel_factor, &
                        crm_accel_uv_out = crm_accel_uv)
  
      if (masterproc) then
         write(iulog, *) 'USING CRM MEAN STATE ACCELERATION'
         write(iulog, *) 'crm_accel: use_crm_accel = ', use_crm_accel
         write(iulog, *) 'crm_accel: crm_accel_factor = ', crm_accel_factor
         write(iulog, *) 'crm_accel: crm_accel_uv = ', crm_accel_uv
      endif
#if !defined(sam1mom)
      if (masterproc) then
        write(0,*) "CRM time step relaxation is only compatible with sam1mom microphysics"
        call endrun('crm main')
      endif
#endif
    end subroutine crm_accel_init


    subroutine crm_accel_nstop(nstop)
      ! Reduces nstop to appropriate value give crm_accel_factor.
      ! 
      ! To correctly apply mean-state acceleration in the crm_module/crm
      ! subroutine, nstop must be reduced to nstop / (1 + crm_accel_factor)
      !
      ! Argument(s):
      !  nstop (inout) - number of crm iterations to apply MSA
      ! -----------------------------------------------------------------------
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
      endif
    end subroutine crm_accel_nstop


    subroutine accelerate_crm(ncrms, nstep, nstop, ceaseflag)
      ! Applies mean-state acceleration (MSA) to CRM
      !
      ! Applies MSA to the following crm fields:
      !   t, qv, qcl, qci, micro_field(:,:,:,index_water_vapor,:),
      !   u (optional), v (optional)
      ! Raises ceaseflag and aborts MSA if the magnitude of 
      ! the change in "t" exceeds 5K at any point.
      ! 
      ! Arguments:
      !   ncrms (in) - number of crm columns in this group
      !   nstep (in) - current crm iteration, needed only if 
      !                ceaseflag is triggered
      !   nstop (inout) - number of crm iterations, adjusted only
      !                   if ceaseflag is triggered
      !   ceaseflag (inout) - returns true if accelerate_crm aborted
      !                       before MSA applied; otherwise false
      ! Notes:
      !   Intended to be called from crm subroutine in crm_module
      ! -----------------------------------------------------------------------
      use grid, only: nzm
      use vars, only: u, v, u0, v0, t0,q0, t,qcl,qci,qv
      use microphysics, only: micro_field, ixw=>index_water_vapor
      use cam_logfile,  only: iulog
      implicit none
      integer, intent(in   ) :: ncrms
      integer, intent(in   ) :: nstep
      integer, intent(inout) :: nstop
      logical, intent(inout) :: ceaseflag
      real(rc) :: ubaccel(nzm,ncrms)   ! u before applying MSA tendency
      real(rc) :: vbaccel(nzm,ncrms)   ! v before applying MSA tendency
      real(rc) :: tbaccel(nzm,ncrms)   ! t before applying MSA tendency
      real(rc) :: qtbaccel(nzm,ncrms)  ! Non-precipitating qt before applying MSA tendency
      real(rc) :: ttend_acc(nzm,ncrms) ! MSA adjustment of t
      real(rc) :: qtend_acc(nzm,ncrms) ! MSA adjustment of qt
      real(rc) :: utend_acc(nzm,ncrms) ! MSA adjustment of u
      real(rc) :: vtend_acc(nzm,ncrms) ! MSA adjustment of v
      real(rc) :: tmp  ! temporary variable for atomic updates
      integer i, j, k, icrm  ! iteration variables
      real(r8) :: qpoz(nzm,ncrms) ! total positive micro_field(:,:,k,ixw,:) in level k
      real(r8) :: qneg(nzm,ncrms) ! total negative micro_field(:,:,k,ixw,:) in level k
      real(r8) :: factor, qfactor ! local variables for redistributing moisture
      real(rc) :: ttend_threshold ! threshold for ttend_acc at which MSA aborts
      real(rc) :: tmin  ! mininum value of t allowed (sanity factor)

      ttend_threshold = 5.  ! 5K, following UP-CAM implementation
      tmin = 50.  ! should never get below 50K in crm, following UP-CAM implementation

      !$acc enter data create(qpoz,qneg,ubaccel,vbaccel,tbaccel,qtbaccel,ttend_acc,qtend_acc,utend_acc,vtend_acc) async(asyncid)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Compute the average among horizontal columns for each variable
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !$acc parallel loop collapse(2) async(asyncid)
      do icrm = 1, ncrms
        do k = 1, nzm
          tbaccel (k,icrm) = 0
          qtbaccel(k,icrm) = 0
          if (crm_accel_uv) then
            ubaccel (k,icrm) = 0
            vbaccel(k,icrm) = 0
          endif
        enddo
      enddo
      !$acc parallel loop collapse(4) copyin(t,qcl,qci,qv,u,v) async(asyncid)
      do icrm = 1, ncrms
        do k = 1, nzm
          do j = 1 , ny
            do i = 1 , nx
              ! calculate tendency * dtn
              tmp = t(i,j,k,icrm) * coef
              !$acc atomic update
              tbaccel (k,icrm) = tbaccel (k,icrm) + tmp
              tmp = (qcl(i,j, k, icrm) + qci(i,j, k, icrm) + qv(i,j, k, icrm)) * coef
              !$acc atomic update
              qtbaccel(k,icrm) = qtbaccel(k,icrm) + tmp
              if (crm_accel_uv) then
                tmp = u(i,j,k,icrm) * coef
                !$acc atomic update
                ubaccel(k,icrm) = ubaccel(k,icrm) + tmp
                tmp = v(i,j,k,icrm) * coef
                !$acc atomic update
                vbaccel(k,icrm) = vbaccel(k,icrm) + tmp
              endif
            enddo
          enddo
        enddo
      enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Compute the accelerated tendencies
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !$acc parallel loop collapse(2) copyin(t0,q0,u0,v0) copy(ceaseflag) async(asyncid)
      do icrm = 1, ncrms
        do k = 1, nzm
          ttend_acc(k,icrm) = tbaccel (k,icrm) - t0(k,icrm)
          qtend_acc(k,icrm) = qtbaccel(k,icrm) - q0(k,icrm)
          if (crm_accel_uv) then
            utend_acc(k,icrm) = ubaccel(k,icrm) - u0(k,icrm)
            vtend_acc(k,icrm) = vbaccel(k,icrm) - v0(k,icrm)
          endif
          if (ttend_acc(k,icrm) > ttend_threshold) then
            ceaseflag = .true.
          endif
        enddo
      enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Make sure it isn't insane
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !$acc wait(asyncid)
      if (ceaseflag) then ! special case for dT/dt too large
        write (iulog, *) 'accelerate_crm: mean-state acceleration not applied this step'
        ! reset nstop so remainder of this crm integration is carried out without MSA
        write (iulog,*) 'crm: nstop increased from ', nstop, ' to ', int(nstop+(nstop-nstep+1)*crm_accel_factor)
        nstop = nstop + (nstop - nstep + 1)*crm_accel_factor ! only can happen once
        !$acc exit data delete(qpoz,qneg,ubaccel,vbaccel,tbaccel,qtbaccel,ttend_acc,qtend_acc,utend_acc,vtend_acc) async(asyncid)
        return
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Apply the accelerated tendencies
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !$acc parallel loop collapse(4) copy(t,u,v,micro_field) async(asyncid)
      do icrm = 1, ncrms
        do k = 1, nzm
          do j = 1, ny
            do i = 1, nx
              ! don't let T go negative!
              t(i,j,k,icrm) = max(tmin, t(i,j,k,icrm) + crm_accel_factor * ttend_acc(k,icrm))
              if (crm_accel_uv) then
                u(i,j,k,icrm) = u(i,j,k,icrm) + crm_accel_factor * utend_acc(k,icrm) 
                v(i,j,k,icrm) = v(i,j,k,icrm) + crm_accel_factor * vtend_acc(k,icrm) 
              endif
              micro_field(i,j,k,icrm,ixw) = micro_field(i,j,k,icrm,ixw) + crm_accel_factor * qtend_acc(k,icrm)
            enddo
          enddo
        enddo
      enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Fix negative micro and readjust among micro components
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !$acc parallel loop collapse(2) async(asyncid)
      do icrm = 1, ncrms
        do k = 1, nzm
          qpoz(k,icrm) = 0.
          qneg(k,icrm) = 0.
        enddo
      enddo
      !$acc parallel loop collapse(4) copyin(micro_field) async(asyncid)
      do icrm = 1, ncrms
        do k = 1, nzm
          do j = 1, ny
            do i = 1, nx
              if (micro_field(i,j,k,icrm,ixw) < 0.) then
                !$acc atomic update
                qneg(k,icrm) = qneg(k,icrm) + micro_field(i,j,k,icrm,ixw)
              else
                !$acc atomic update
                qpoz(k,icrm) = qpoz(k,icrm) + micro_field(i,j,k,icrm,ixw)
              endif
            enddo
          enddo
        enddo
      enddo
      !$acc parallel loop collapse(4) copy(micro_field,qv,qcl,qci) async(asyncid)
      do icrm = 1, ncrms
        do k = 1, nzm
          do j = 1 , ny
            do i = 1 , nx
              if (qpoz(k,icrm) + qneg(k,icrm) <= 0.) then
                micro_field(i,j,k,icrm,ixw) = 0.
                qv         (i,j,k,icrm    ) = 0.
                qcl        (i,j,k,icrm    ) = 0.
                qci        (i,j,k,icrm    ) = 0.
              else
                factor = 1._r8 + qneg(k,icrm) / qpoz(k,icrm)
                ! apply to micro_field
                micro_field(i,j,k,icrm,ixw) = max(0._rc, micro_field(i,j,k,icrm,ixw) * factor)
                ! partition micro_field into qv, qcl, and qci
                ! such that micro_field = qv + qcl + qci, following these rules:
                !  (1) adjust qv first
                !  (2) adjust qcl and qci only if needed to ensure positivity
                if (micro_field(i,j,k,icrm,ixw) <= 0._rc) then
                  qv (i,j,k,icrm) = 0.
                  qcl(i,j,k,icrm) = 0.
                  qci(i,j,k,icrm) = 0.
                else
                  qfactor = micro_field(i,j,k,icrm,ixw) - qcl(i,j,k,icrm) - qci(i,j,k,icrm)
                  qv(i,j,k,icrm) = max(0._rc, qfactor)
                  if (qfactor < 0._r8) then
                    qfactor = 1._r8 + qfactor / (qcl(i,j,k,icrm) + qci(i,j,k,icrm))
                    qcl(i,j,k,icrm) = qcl(i,j,k,icrm) * qfactor
                    qci(i,j,k,icrm) = qci(i,j,k,icrm) * qfactor
                  endif
                endif
              endif ! qpoz + qneg < 0.
            enddo
          enddo
        enddo  ! k = 1, nzm
      enddo ! icrm = 1, ncrms

      !$acc exit data delete(qpoz,qneg,ubaccel,vbaccel,tbaccel,qtbaccel,ttend_acc,qtend_acc,utend_acc,vtend_acc) async(asyncid)
    end subroutine accelerate_crm
    
end module accelerate_crm_mod
