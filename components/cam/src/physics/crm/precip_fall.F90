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
    real(crm_rknd) :: mx(nzm),mn(nzm), lfac(nz,ncrms)
    real(crm_rknd) :: www(nz),fz(nz)
    real(crm_rknd) :: df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real(crm_rknd) :: f0(nzm),df0(nzm)
    real(crm_rknd) :: eps
    integer :: i,j,k,kc,kb,icrm
    logical :: nonos
    real(crm_rknd) :: y,pp,pn
    real(crm_rknd) :: lat_heat, wmax
    real(crm_rknd) :: wp(nzm), tmp_qp(nzm), irhoadz(nzm), iwmax(nzm), rhofac(nzm), prec_cfl
    integer nprec, iprec
    real(crm_rknd) :: flagstat
    real(crm_rknd), pointer :: qp(:,:,:,:)  ! total precipitating water

    !Statement functions
    pp(y)= max(real(0.,crm_rknd),y)
    pn(y)=-min(real(0.,crm_rknd),y)

    do icrm = 1 , ncrms

    qp(dimx1_s:,dimy1_s:,1:,1:) => micro_field(:,:,:,2,:)

    eps = 1.e-10
    nonos = .true.

    do k = 1,nzm
      rhofac(k) = sqrt(1.29/rho(k,icrm))
      irhoadz(k) = 1./(rho(k,icrm)*adz(k,icrm)) ! Useful factor
      kb = max(1,k-1)
      wmax       = dz(icrm)*adz(kb,icrm)/dtn   ! Velocity equivalent to a cfl of 1.0.
      iwmax(k)   = 1./wmax
    enddo

    ! 	Add sedimentation of precipitation field to the vert. vel.
    do j=1,ny
      do i=1,nx
        ! Compute precipitation velocity and flux column-by-column
        prec_cfl = 0.
        do k=1,nzm

          select case (hydro_type)
          case(0)
            lfac(k,icrm) = fac_cond
            flagstat = 1.
          case(1)
            lfac(k,icrm) = fac_sub
            flagstat = 1.
          case(2)
            lfac(k,icrm) = fac_cond + (1-omega(i,j,k,icrm))*fac_fus
            flagstat = 1.
          case(3)
            lfac(k,icrm) = 0.
            flagstat = 0.
          case default
            if(masterproc) then
              print*, 'unknown hydro_type in precip_fall. exitting ...'
              call task_abort
            endif
          end select

          wp(k)=rhofac(k)*term_vel(ncrms,icrm,i,j,k,ind)
          prec_cfl = max(prec_cfl,wp(k)*iwmax(k)) ! Keep column maximum CFL
          wp(k) = -wp(k)*rhow(k,icrm)*dtn/dz(icrm)

        enddo  ! k

        fz(nz)=0.
        www(nz)=0.
        lfac(nz,icrm)=0

        ! If maximum CFL due to precipitation velocity is greater than 0.9,
        ! take more than one advection step to maintain stability.
        if (prec_cfl.gt.0.9) then
          nprec = CEILING(prec_cfl/0.9)
          do k = 1,nzm
            ! wp already includes factor of dt, so reduce it by a
            ! factor equal to the number of precipitation steps.
            wp(k) = wp(k)/real(nprec,crm_rknd)
          enddo
        else
          nprec = 1
        endif

        !  loop over iterations
        do iprec = 1,nprec
          do k = 1,nzm
            tmp_qp(k) = qp(i,j,k,icrm) ! Temporary array for qp in this column
          enddo

          if(nonos) then
            do k=1,nzm
              kc=min(nzm,k+1)
              kb=max(1,k-1)
              mx(k)=max(tmp_qp(kb),tmp_qp(kc),tmp_qp(k))
              mn(k)=min(tmp_qp(kb),tmp_qp(kc),tmp_qp(k))
            enddo
          endif  ! nonos

          do k=1,nzm
            ! Define upwind precipitation flux
            fz(k)=tmp_qp(k)*wp(k)
          enddo

          do k=1,nzm
            kc=k+1
            tmp_qp(k)=tmp_qp(k)-(fz(kc)-fz(k))*irhoadz(k) !Update temporary qp
          enddo

          do k=1,nzm
            ! Also, compute anti-diffusive correction to previous
            ! (upwind) approximation to the flux
            kb=max(1,k-1)
            ! The precipitation velocity is a cell-centered quantity,
            ! since it is computed from the cell-centered
            ! precipitation mass fraction.  Therefore, a reformulated
            ! anti-diffusive flux is used here which accounts for
            ! this and results in reduced numerical diffusion.
            www(k) = 0.5*(1.+wp(k)*irhoadz(k))*(tmp_qp(kb)*wp(kb) - tmp_qp(k)*wp(k)) ! works for wp(k)<0
          enddo

          !---------- non-osscilatory option ---------------
          if(nonos) then
            do k=1,nzm
              kc=min(nzm,k+1)
              kb=max(1,k-1)
              mx(k)=max(tmp_qp(kb),tmp_qp(kc),tmp_qp(k),mx(k))
              mn(k)=min(tmp_qp(kb),tmp_qp(kc),tmp_qp(k),mn(k))
            enddo
            do k=1,nzm
              kc=min(nzm,k+1)
              mx(k)=rho(k,icrm)*adz(k,icrm)*(mx(k)-tmp_qp(k))/(pn(www(kc)) + pp(www(k))+eps)
              mn(k)=rho(k,icrm)*adz(k,icrm)*(tmp_qp(k)-mn(k))/(pp(www(kc)) + pn(www(k))+eps)
            enddo
            do k=1,nzm
              kb=max(1,k-1)
              ! Add limited flux correction to fz(k).
              fz(k) = fz(k) + pp(www(k))*min(real(1.,crm_rknd),mx(k), mn(kb)) - pn(www(k))*min(real(1.,crm_rknd),mx(kb),mn(k)) ! Anti-diffusive flux
            enddo
          endif ! nonos

          ! Update precipitation mass fraction and liquid-ice static
          ! energy using precipitation fluxes computed in this column.
          do k=1,nzm
            kc=k+1
            ! Update precipitation mass fraction.
            ! Note that fz is the total flux, including both the
            ! upwind flux and the anti-diffusive correction.
            qp(i,j,k,icrm)=qp(i,j,k,icrm)-(fz(kc)-fz(k))*irhoadz(k)
            qpfall(k,icrm)=qpfall(k,icrm)-(fz(kc)-fz(k))*irhoadz(k)*flagstat  ! For qp budget
            lat_heat = -(lfac(kc,icrm)*fz(kc)-lfac(k,icrm)*fz(k))*irhoadz(k)
            t(i,j,k,icrm)=t(i,j,k,icrm)-lat_heat
            tlat(k,icrm)=tlat(k,icrm)-lat_heat            ! For energy budget
            precflux(k,icrm) = precflux(k,icrm) - fz(k)*flagstat   ! For statistics
          enddo
          precsfc(i,j,icrm) = precsfc(i,j,icrm) - fz(1)*flagstat ! For statistics
          precssfc(i,j,icrm) = precssfc(i,j,icrm) - fz(1)*(1.-omega(i,j,1,icrm))*flagstat ! For statistics
          prec_xy(i,j,icrm) = prec_xy(i,j,icrm) - fz(1)*flagstat ! For 2D output

          if (iprec.lt.nprec) then
            ! Re-compute precipitation velocity using new value of qp.
            do k=1,nzm
              wp(k) = rhofac(k)*term_vel(ncrms,icrm,i,j,k,ind)
              ! Decrease precipitation velocity by factor of nprec
              wp(k) = -wp(k)*rhow(k,icrm)*dtn/dz(icrm)/real(nprec,crm_rknd)
              ! Note: Don't bother checking CFL condition at each
              ! substep since it's unlikely that the CFL will
              ! increase very much between substeps when using
              ! monotonic advection schemes.
            enddo
            fz(nz)=0.
            www(nz)=0.
            lfac(nz,icrm)=0.
          endif

        enddo !iprec
      enddo
    enddo

    enddo

  end subroutine precip_fall

end module precip_fall_mod
