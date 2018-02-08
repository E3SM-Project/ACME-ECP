module precip_fall_mod
  use task_util_mod
  use bound_duvdt_mod
  implicit none

contains

  subroutine precip_fall(qp, term_vel, hydro_type, omega, ind, ncrms)
    !     positively definite monotonic advection with non-oscillatory option
    !     and gravitational sedimentation
    use vars
    use params
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) qp(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! falling hydrometeor
    integer hydro_type   ! 0 - all liquid, 1 - all ice, 2 - mixed
    real(crm_rknd) omega(ncrms,nx,ny,nzm)   !  = 1: liquid, = 0: ice;  = 0-1: mixed : used only when hydro_type=2
    integer ind

    ! Terminal velocity fnction
    real(crm_rknd), external :: term_vel  ! terminal velocity function

    ! Local:
    real(crm_rknd), allocatable :: mx(:,:),mn(:,:), lfac(:,:)
    real(crm_rknd), allocatable :: www(:,:),fz(:,:)
    real(crm_rknd) eps
    integer i,j,k,kc,kb,icrm
    logical nonos

    real(crm_rknd) y,pp,pn
    pp(y)= max(real(0.,crm_rknd),y)
    pn(y)=-min(real(0.,crm_rknd),y)

    real(crm_rknd) lat_heat, wmax
    real(crm_rknd), allocatable :: wp(:,:), tmp_qp(:,:), irhoadz(:,:), iwmax(:,:), rhofac(:,:), prec_cfl(:)
    integer, allocatable :: nprec(:)
    integer iprec, nprec_max
    real(crm_rknd), allocatable :: flagstat(:)

    allocate(mx(ncrms,nzm))
    allocate(mn(ncrms,nzm))
    allocate(lfac(ncrms,nz))
    allocate(www(ncrms,nz))
    allocate(fz(ncrms,nz))
    allocate(wp(ncrms,nzm))
    allocate(tmp_qp(ncrms,nzm))
    allocate(irhoadz(ncrms,nzm))
    allocate(iwmax(ncrms,nzm))
    allocate(rhofac(ncrms,nzm))
    allocate(prec_cfl(ncrms))
    allocate(nprec(ncrms))
    allocate(flagstat(ncrms))

    !--------------------------------------------------------
    !call t_startf ('precip_fall')

    eps = 1.e-10
    nonos = .true.

    do k = 1,nzm
      do icrm = 1 , ncrms
        rhofac(icrm,k) = sqrt(1.29/rho(icrm,k))
        irhoadz(icrm,k) = 1./(rho(icrm,k)*adz(icrm,k)) ! Useful factor
        kb = max(1,k-1)
        wmax       = dz(icrm)*adz(icrm,kb)/dtn   ! Velocity equivalent to a cfl of 1.0.
        iwmax(icrm,k)   = 1./wmax
      end do
    end do

    ! 	Add sedimentation of precipitation field to the vert. vel.
    do j=1,ny
      do i=1,nx
        ! Compute precipitation velocity and flux column-by-column
        prec_cfl(:) = 0.

        do k=1,nzm
          do icrm = 1 , ncrms

            select case (hydro_type)
            case(0)
              lfac(icrm,k) = fac_cond
              flagstat(icrm) = 1.
            case(1)
              lfac(icrm,k) = fac_sub
              flagstat(icrm) = 1.
            case(2)
              lfac(icrm,k) = fac_cond + (1-omega(icrm,i,j,k))*fac_fus
              flagstat(icrm) = 1.
            case(3)
              lfac(icrm,k) = 0.
              flagstat(icrm) = 0.
            case default
              if(masterproc) then
                print*, 'unknown hydro_type in precip_fall. exitting ...'
                call task_abort
              end if
            end select

            wp(icrm,k)=rhofac(icrm,k)*term_vel(i,j,k,ind,ncrms,icrm)
            prec_cfl(icrm) = max(prec_cfl(icrm),wp(icrm,k)*iwmax(icrm,k)) ! Keep column maximum CFL
            wp(icrm,k) = -wp(icrm,k)*rhow(icrm,k)*dtn/dz(icrm)

          end do  ! k
        end do  ! k

        do icrm = 1 , ncrms
          fz(icrm,nz)=0.
          www(icrm,nz)=0.
          lfac(icrm,nz)=0
        end do  ! k

        nprec_max = 0
        do icrm = 1 , ncrms
          ! If maximum CFL due to precipitation velocity is greater than 0.9,
          ! take more than one advection step to maintain stability.
          if (prec_cfl(icrm).gt.0.9) then
            nprec(icrm) = CEILING(prec_cfl(icrm)/0.9)
            do k = 1,nzm
              ! wp already includes factor of dt, so reduce it by a
              ! factor equal to the number of precipitation steps.
              wp(icrm,k) = wp(icrm,k)/real(nprec(icrm),crm_rknd)
            end do
          else
            nprec(icrm) = 1
          end if
          nprec_max = max(nprec(icrm),nprec_max)
        enddo

        do iprec = 1,nprec_max

          do k = 1,nzm
            do icrm = 1 , ncrms
              if (iprec <= nprec(icrm)) then
                tmp_qp(icrm,k) = qp(icrm,i,j,k) ! Temporary array for qp in this column
              endif
            end do
          end do

          !-----------------------------------------

          if(nonos) then

            do k=1,nzm
              do icrm = 1 , ncrms
                if (iprec <= nprec(icrm)) then
                  kc=min(nzm,k+1)
                  kb=max(1,k-1)
                  mx(icrm,k)=max(tmp_qp(icrm,kb),tmp_qp(icrm,kc),tmp_qp(icrm,k))
                  mn(icrm,k)=min(tmp_qp(icrm,kb),tmp_qp(icrm,kc),tmp_qp(icrm,k))
                endif
              end do
            end do

          end if  ! nonos

          !  loop over iterations

          do k=1,nzm
            do icrm = 1 , ncrms
              if (iprec <= nprec(icrm)) then
                ! Define upwind precipitation flux
                fz(icrm,k)=tmp_qp(icrm,k)*wp(icrm,k)
              endif
            end do
          end do

          do k=1,nzm
            do icrm = 1 , ncrms
              if (iprec <= nprec(icrm)) then
                kc=k+1
                tmp_qp(icrm,k)=tmp_qp(icrm,k)-(fz(icrm,kc)-fz(icrm,k))*irhoadz(icrm,k) !Update temporary qp
              endif
            end do
          end do

          do k=1,nzm
            do icrm = 1 , ncrms
              if (iprec <= nprec(icrm)) then
                ! Also, compute anti-diffusive correction to previous
                ! (upwind) approximation to the flux
                kb=max(1,k-1)
                ! The precipitation velocity is a cell-centered quantity,
                ! since it is computed from the cell-centered
                ! precipitation mass fraction.  Therefore, a reformulated
                ! anti-diffusive flux is used here which accounts for
                ! this and results in reduced numerical diffusion.
                www(icrm,k) = 0.5*(1.+wp(icrm,k)*irhoadz(icrm,k)) &
                *(tmp_qp(icrm,kb)*wp(icrm,kb) - tmp_qp(icrm,k)*wp(icrm,k)) ! works for wp(k)<0
              endif
            end do
          end do

          !---------- non-osscilatory option ---------------

          if(nonos) then

            do k=1,nzm
              do icrm = 1 , ncrms
                if (iprec <= nprec(icrm)) then
                  kc=min(nzm,k+1)
                  kb=max(1,k-1)
                  mx(icrm,k)=max(tmp_qp(icrm,kb),tmp_qp(icrm,kc),tmp_qp(icrm,k),mx(icrm,k))
                  mn(icrm,k)=min(tmp_qp(icrm,kb),tmp_qp(icrm,kc),tmp_qp(icrm,k),mn(icrm,k))
                endif
              end do
            end do

            do k=1,nzm
              do icrm = 1 , ncrms
                if (iprec <= nprec(icrm)) then
                  kc=min(nzm,k+1)
                  mx(icrm,k)=rho(icrm,k)*adz(icrm,k)*(mx(icrm,k)-tmp_qp(icrm,k))/(pn(www(icrm,kc)) + pp(www(icrm,k))+eps)
                  mn(icrm,k)=rho(icrm,k)*adz(icrm,k)*(tmp_qp(icrm,k)-mn(icrm,k))/(pp(www(icrm,kc)) + pn(www(icrm,k))+eps)
                endif
              end do
            end do

            do k=1,nzm
              do icrm = 1 , ncrms
                if (iprec <= nprec(icrm)) then
                  kb=max(1,k-1)
                  ! Add limited flux correction to fz(icrm,k).
                  fz(icrm,k) = fz(icrm,k) &                        ! Upwind flux
                  + pp(www(icrm,k))*min(real(1.,crm_rknd),mx(icrm,k), mn(icrm,kb)) &
                  - pn(www(icrm,k))*min(real(1.,crm_rknd),mx(icrm,kb),mn(icrm,k)) ! Anti-diffusive flux
                endif
              end do
            end do

          endif ! nonos

          ! Update precipitation mass fraction and liquid-ice static
          ! energy using precipitation fluxes computed in this column.
          do k=1,nzm
            do icrm = 1 , ncrms
              if (iprec <= nprec(icrm)) then
                kc=k+1
                ! Update precipitation mass fraction.
                ! Note that fz is the total flux, including both the
                ! upwind flux and the anti-diffusive correction.
                qp(icrm,i,j,k)=qp(icrm,i,j,k)-(fz(icrm,kc)-fz(icrm,k))*irhoadz(icrm,k)
                qpfall(icrm,k)=qpfall(icrm,k)-(fz(icrm,kc)-fz(icrm,k))*irhoadz(icrm,k)*flagstat(icrm)  ! For qp budget
                lat_heat = -(lfac(icrm,kc)*fz(icrm,kc)-lfac(icrm,k)*fz(icrm,k))*irhoadz(icrm,k)
                t(icrm,i,j,k)=t(icrm,i,j,k)-lat_heat
                tlat(icrm,k)=tlat(icrm,k)-lat_heat            ! For energy budget
                precflux(icrm,k) = precflux(icrm,k) - fz(icrm,k)*flagstat(icrm)   ! For statistics
              endif
            end do
          end do
          do icrm = 1 , ncrms
            if (iprec <= nprec(icrm)) then
              precsfc(icrm,i,j) = precsfc(icrm,i,j) - fz(icrm,1)*flagstat(icrm) ! For statistics
              precssfc(icrm,i,j) = precssfc(icrm,i,j) - fz(icrm,1)*(1.-omega(icrm,i,j,1))*flagstat(icrm) ! For statistics
              prec_xy(icrm,i,j) = prec_xy(icrm,i,j) - fz(icrm,1)*flagstat(icrm) ! For 2D output
            endif
          end do


          ! Re-compute precipitation velocity using new value of qp.
          do k=1,nzm
            do icrm = 1 , ncrms
              if (iprec.lt.nprec(icrm)) then
                wp(icrm,k) = rhofac(icrm,k)*term_vel(i,j,k,ind,ncrms,icrm)
                ! Decrease precipitation velocity by factor of nprec(icrm)
                wp(icrm,k) = -wp(icrm,k)*rhow(icrm,k)*dtn/dz(icrm)/real(nprec(icrm),crm_rknd)
                ! Note: Don't bother checking CFL condition at each
                ! substep since it's unlikely that the CFL will
                ! increase very much between substeps when using
                ! monotonic advection schemes.
              endif
            end do
          end do

          do icrm = 1 , ncrms
            if (iprec.lt.nprec(icrm)) then
              fz(icrm,nz)=0.
              www(icrm,nz)=0.
              lfac(icrm,nz)=0.
            endif
          end do

        end do !iprec

      end do
    end do

    deallocate(mx)
    deallocate(mn)
    deallocate(lfac)
    deallocate(www)
    deallocate(fz)
    deallocate(wp)
    deallocate(tmp_qp)
    deallocate(irhoadz)
    deallocate(iwmax)
    deallocate(rhofac)
    deallocate(prec_cfl)
    deallocate(nprec)
    deallocate(flagstat)


    !call t_stopf ('precip_fall')

  end subroutine precip_fall


end module precip_fall_mod
