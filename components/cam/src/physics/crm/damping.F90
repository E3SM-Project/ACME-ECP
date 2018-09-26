module damping_mod
	use task_util_mod
	implicit none

contains

  subroutine damping(ncrms,icrm)

    !  "Spange"-layer damping at the domain top region

    use vars
    use microphysics, only: micro_field, index_water_vapor
    use params, only: crm_rknd
    implicit none
    integer, intent(in) :: ncrms,icrm
    real(crm_rknd) tau_min	! minimum damping time-scale (at the top)
    real(crm_rknd) tau_max    ! maxim damping time-scale (base of damping layer)
    real(crm_rknd) damp_depth ! damping depth as a fraction of the domain height
    parameter(tau_min=60., tau_max=450., damp_depth=0.4)
    real(crm_rknd) tau(nzm)
    integer i, j, k, n_damp

    if(tau_min.lt.2*dt) then
      print*,'Error: in damping() tau_min is too small!'
      call task_abort()
    end if

    do k=nzm,1,-1
      if(z(nzm,icrm)-z(k,icrm).lt.damp_depth*z(nzm,icrm)) then
        n_damp=nzm-k+1
      endif
    end do

    do k=nzm,nzm-n_damp,-1
      tau(k) = tau_min *(tau_max/tau_min)**((z(nzm,icrm)-z(k,icrm))/(z(nzm,icrm)-z(nzm-n_damp,icrm)))
      tau(k)=1./tau(k)
    end do

    !+++mhwang recalculate grid-mean u0, v0, t0 first,
    ! as t have been updated. No need for qv0, as
    ! qv has not been updated yet the calculation of qv0.
    do k=1, nzm
      u0(k,icrm)=0.0
      v0(k,icrm)=0.0
      t0(k,icrm)=0.0
      do j=1, ny
        do i=1, nx
          u0(k,icrm) = u0(k,icrm) + u(i,j,k,icrm)/(nx*ny)
          v0(k,icrm) = v0(k,icrm) + v(i,j,k,icrm)/(nx*ny)
          t0(k,icrm) = t0(k,icrm) + t(i,j,k,icrm)/(nx*ny)
        end do
      end do
    end do
    !---mhwang

    do k = nzm, nzm-n_damp, -1
      do j=1,ny
        do i=1,nx
          dudt(i,j,k,na(icrm),icrm)= dudt(i,j,k,na(icrm),icrm)-(u(i,j,k,icrm)-u0(k,icrm)) * tau(k)
          dvdt(i,j,k,na(icrm),icrm)= dvdt(i,j,k,na(icrm),icrm)-(v(i,j,k,icrm)-v0(k,icrm)) * tau(k)
          dwdt(i,j,k,na(icrm),icrm)= dwdt(i,j,k,na(icrm),icrm)-w(i,j,k,icrm) * tau(k)
          t(i,j,k,icrm)= t(i,j,k,icrm)-dtn*(t(i,j,k,icrm)-t0(k,icrm)) * tau(k)
          ! In the old version (SAM7.5?) of SAM, water vapor is the prognostic variable for the two-moment microphyscs.
          ! So the following damping approach can lead to the negative water vapor.
          !      micro_field(i,j,k,index_water_vapor,icrm)= micro_field(i,j,k,index_water_vapor,icrm)- &
          !                                    dtn*(qv(i,j,k,icrm)+qcl(i,j,k,icrm)+qci(i,j,k,icrm)-q0(k,icrm)) * tau(k)
          ! a simple fix (Minghuai Wang, 2011-08):
          micro_field(i,j,k,index_water_vapor,icrm)= micro_field(i,j,k,index_water_vapor,icrm)- &
          dtn*(qv(i,j,k,icrm)-qv0(k,icrm)) * tau(k)
        end do! i
      end do! j
    end do ! k

  end subroutine damping

end module damping_mod
