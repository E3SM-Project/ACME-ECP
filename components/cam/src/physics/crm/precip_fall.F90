module precip_fall_mod
  use task_util_mod
  use bound_duvdt_mod
  implicit none

contains

  subroutine precip_fall(ncrms,micro_field, term_vel, hydro_type, omega, ind)
    !     positively definite monotonic advection with non-oscillatory option
    !     and gravitational sedimentation
    use vars
    use params
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd), target :: micro_field(dimx1_s:,dimy1_s:, :, :,:)
    integer :: hydro_type   ! 0 - all liquid, 1 - all ice, 2 - mixed
    real(crm_rknd) :: omega(nx,ny,nzm,ncrms)   !  = 1: liquid, = 0: ice;  = 0-1: mixed : used only when hydro_type=2
    integer :: ind
    ! Terminal velocity fnction
    real(crm_rknd), external :: term_vel  ! terminal velocity function
    ! Local:
    real(crm_rknd) :: mx(nx,ny,nzm,ncrms),mn(nx,ny,nzm,ncrms), lfac(nx,ny,nz,ncrms)
    real(crm_rknd) :: www(nx,ny,nz,ncrms),fz(nx,ny,nz,ncrms)
    real(crm_rknd) :: eps
    integer :: i,j,k,kc,kb,icrm
    logical :: nonos
    real(crm_rknd) :: y,pp,pn
    real(crm_rknd) :: lat_heat, wmax
    real(crm_rknd) :: wp(nx,ny,nzm,ncrms), tmp_qp(nx,ny,nzm,ncrms), irhoadz(nzm,ncrms), iwmax(nzm,ncrms), rhofac(nzm,ncrms), prec_cfl
    integer nprec, iprec
    real(crm_rknd) :: flagstat
    real(crm_rknd), pointer :: qp(:,:,:,:)  ! total precipitating water

    !Statement functions
    pp(y)= max(real(0.,crm_rknd),y)
    pn(y)=-min(real(0.,crm_rknd),y)

    qp(dimx1_s:,dimy1_s:,1:,1:) => micro_field(:,:,:,2,:)

    eps = 1.e-10
    nonos = .true.

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
                print*, 'unknown hydro_type in precip_fall. exitting ...'
                call task_abort
              endif
            end select
            wp(i,j,k,icrm)=rhofac(k,icrm)*term_vel(ncrms,icrm,i,j,k,ind)
            prec_cfl = max(prec_cfl,wp(i,j,k,icrm)*iwmax(k,icrm)) ! Keep column maximum CFL
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
      do icrm = 1 , ncrms
        do k = 1,nzm
          do j=1,ny
            do i=1,nx
              tmp_qp(i,j,k,icrm) = qp(i,j,k,icrm) ! Temporary array for qp in this column
            enddo
          enddo
        enddo
      enddo

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
              www(i,j,k,icrm) = 0.5*(1.+wp(i,j,k,icrm)*irhoadz(k,icrm))*(tmp_qp(i,j,kb,icrm)*wp(i,j,kb,icrm) - tmp_qp(i,j,k,icrm)*wp(i,j,k,icrm)) ! works for wp(k)<0
            enddo
          enddo
        enddo
      enddo

      !---------- non-osscilatory option ---------------
      if(nonos) then
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
        do icrm = 1 , ncrms
          do k=1,nzm
            do j=1,ny
              do i=1,nx
                kb=max(1,k-1)
                ! Add limited flux correction to fz(k).
                fz(i,j,k,icrm) = fz(i,j,k,icrm) + pp(www(i,j,k,icrm))*min(real(1.,crm_rknd),mx(i,j,k,icrm), mn(i,j,kb,icrm)) - pn(www(i,j,k,icrm))*min(real(1.,crm_rknd),mx(i,j,kb,icrm),mn(i,j,k,icrm)) ! Anti-diffusive flux
              enddo
            enddo
          enddo
        enddo
      endif ! nonos

      ! Update precipitation mass fraction and liquid-ice static
      ! energy using precipitation fluxes computed in this column.
      do icrm = 1 , ncrms
        do j=1,ny
          do i=1,nx
            do k=1,nzm
              kc=k+1
              ! Update precipitation mass fraction.
              ! Note that fz is the total flux, including both the
              ! upwind flux and the anti-diffusive correction.
              qp(i,j,k,icrm)=qp(i,j,k,icrm)-(fz(i,j,kc,icrm)-fz(i,j,k,icrm))*irhoadz(k,icrm)
              qpfall(k,icrm)=qpfall(k,icrm)-(fz(i,j,kc,icrm)-fz(i,j,k,icrm))*irhoadz(k,icrm)*flagstat  ! For qp budget
              lat_heat = -(lfac(i,j,kc,icrm)*fz(i,j,kc,icrm)-lfac(i,j,k,icrm)*fz(i,j,k,icrm))*irhoadz(k,icrm)
              t(i,j,k,icrm)=t(i,j,k,icrm)-lat_heat
              tlat(k,icrm)=tlat(k,icrm)-lat_heat            ! For energy budget
              precflux(k,icrm) = precflux(k,icrm) - fz(i,j,k,icrm)*flagstat   ! For statistics
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
        do icrm = 1 , ncrms
          do j=1,ny
            do i=1,nx
              do k=1,nzm
                wp(i,j,k,icrm) = rhofac(k,icrm)*term_vel(ncrms,icrm,i,j,k,ind)
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

  end subroutine precip_fall

end module precip_fall_mod
