module ice_fall_mod
	implicit none

contains

  subroutine ice_fall(ncrms)
    ! Sedimentation of ice:
    use vars
    use microphysics, only: micro_field, index_cloud_ice
    !use micro_params
    use params
    implicit none
    integer, intent(in) :: ncrms
    integer i,j,k, kb, kc, kmax(ncrms), kmin(ncrms), ici,icrm
    real(crm_rknd) coef,dqi,lat_heat,vt_ice
    real(crm_rknd) omnu, omnc, omnd, qiu, qic, qid, tmp_theta, tmp_phi
    real(crm_rknd) fz(nx,ny,nz,ncrms)

    do icrm = 1 , ncrms
      kmax(icrm)=0
      kmin(icrm)=nzm+1
    enddo
    do icrm = 1 , ncrms
      !$acc loop seq
      do k = 1,nzm
        !$acc loop seq
        do j = 1, ny
          !$acc loop seq
          do i = 1, nx
            if(qcl(i,j,k,icrm)+qci(i,j,k,icrm).gt.0..and. tabs(i,j,k,icrm).lt.273.15) then
              kmin(icrm) = min(kmin(icrm),k)
              kmax(icrm) = max(kmax(icrm),k)
            end if
          end do
        end do
      end do
    end do

    do icrm = 1 , ncrms
      do k = 1,nzm
        qifall(k,icrm) = 0.
        tlatqi(k,icrm) = 0.
      end do
    end do

    if(index_cloud_ice.eq.-1) return

    do icrm = 1 , ncrms
      do k = 1,nz
        do j = 1, ny
          do i = 1, nx
            fz(i,j,k,icrm) = 0.
          end do
        end do
      end do
    end do

    ! Compute cloud ice flux (using flux limited advection scheme, as in
    ! chapter 6 of Finite Volume Methods for Hyperbolic Problems by R.J.
    ! LeVeque, Cambridge University Press, 2002).
    do icrm = 1 , ncrms
      do k = max(1,kmin(icrm)-1),kmax(icrm)
        do j = 1,ny
          do i = 1,nx
            ! Set up indices for x-y planes above and below current plane.
            kc = min(nzm,k+1)
            kb = max(1,k-1)
            ! CFL number based on grid spacing interpolated to interface i,j,k-1/2
            coef = dtn/(0.5*(adz(kb,icrm)+adz(k,icrm))*dz(icrm))

            ! Compute cloud ice density in this cell and the ones above/below.
            ! Since cloud ice is falling, the above cell is u (upwind,icrm),
            ! this cell is c (center) and the one below is d (downwind).
            qiu = rho(kc,icrm)*qci(i,j,kc,icrm)
            qic = rho(k,icrm) *qci(i,j,k,icrm)
            qid = rho(kb,icrm)*qci(i,j,kb,icrm)

            ! Ice sedimentation velocity depends on ice content. The fiting is
            ! based on the data by Heymsfield (JAS,2003). -Marat
            vt_ice = min(real(0.4,crm_rknd),8.66*(max(real(0.,crm_rknd),qic)+1.e-10)**0.24)   ! Heymsfield (JAS, 2003, p.2607)

            ! Use MC flux limiter in computation of flux correction.
            ! (MC = monotonized centered difference).
            !         if (qic.eq.qid) then
            if (abs(qic-qid).lt.1.0e-25) then  ! when qic, and qid is very small, qic_qid can still be zero
              ! even if qic is not equal to qid. so add a fix here +++mhwang
              tmp_phi = 0.
            else
              tmp_theta = (qiu-qic)/(qic-qid)
              tmp_phi = max(real(0.,crm_rknd),min(0.5*(1.+tmp_theta),real(2.,crm_rknd),2.*tmp_theta))
            end if

            ! Compute limited flux.
            ! Since falling cloud ice is a 1D advection problem, this
            ! flux-limited advection scheme is monotonic.
            fz(i,j,k,icrm) = -vt_ice*(qic - 0.5*(1.-coef*vt_ice)*tmp_phi*(qic-qid))
          end do
        end do
      end do
    enddo
    do icrm = 1 , ncrms
      do j = 1, ny
        do i = 1, nx
          fz(i,j,nz,icrm) = 0.
        end do
      end do
    end do

    ici = index_cloud_ice

    do icrm = 1 , ncrms
      do k=max(1,kmin(icrm)-2),kmax(icrm)
        do j=1,ny
          do i=1,nx
            coef=dtn/(dz(icrm)*adz(k,icrm)*rho(k,icrm))
            ! The cloud ice increment is the difference of the fluxes.
            dqi=coef*(fz(i,j,k,icrm)-fz(i,j,k+1,icrm))
            ! Add this increment to both non-precipitating and total water.
            micro_field(i,j,k,ici,icrm)  = micro_field(i,j,k,ici,icrm)  + dqi
            ! Include this effect in the total moisture budget.
            qifall(k,icrm) = qifall(k,icrm) + dqi

            ! The latent heat flux induced by the falling cloud ice enters
            ! the liquid-ice static energy budget in the same way as the
            ! precipitation.  Note: use latent heat of sublimation.
            lat_heat  = (fac_cond+fac_fus)*dqi
            ! Add divergence of latent heat flux to liquid-ice static energy.
            t(i,j,k,icrm)  = t(i,j,k,icrm)  - lat_heat
            ! Add divergence to liquid-ice static energy budget.
            tlatqi(k,icrm) = tlatqi(k,icrm) - lat_heat
          end do
        end do
      end do
    end do

    do icrm = 1 , ncrms
      do j=1,ny
        do i=1,nx
          coef=dtn/dz(icrm)
          dqi=-coef*fz(i,j,1,icrm)
          precsfc(i,j,icrm) = precsfc(i,j,icrm)+dqi
          precssfc(i,j,icrm) = precssfc(i,j,icrm)+dqi
        end do
      end do
    end do

  end subroutine ice_fall

end module ice_fall_mod
