module damping_mod
  use task_util_mod
  implicit none

contains

  subroutine damping(ncrms)

    !  "Spange"-layer damping at the domain top region

    use vars
    use microphysics, only: micro_field, index_water_vapor
    use params, only: crm_rknd
    use openacc_pool
    implicit none
    integer, intent(in) :: ncrms

    real(crm_rknd) tau_min  ! minimum damping time-scale (at the top)
    real(crm_rknd) tau_max    ! maxim damping time-scale (base of damping layer)
    real(crm_rknd) damp_depth ! damping depth as a fraction of the domain height
    parameter(tau_min=60., tau_max=450., damp_depth=0.4)
    real(crm_rknd), pointer :: tau(:,:)
    integer i, j, k, max_depth, icrm
    integer, pointer :: n_damp(:)

    ! allocate(n_damp(ncrms))
    ! allocate(tau(ncrms,nzm))
    call pool_push(n_damp,(/ncrms/))
    call pool_push(tau,(/ncrms,nzm/))

    if(tau_min.lt.2*dt) then
      print*,'Error: in damping() tau_min is too small!'
      call task_abort()
    end if

    max_depth = 0
    !$acc parallel loop gang vector reduction(max:max_depth)
    do icrm = 1 , ncrms
      !$acc loop seq
      do k=nzm,1,-1
        if(z(icrm,nzm)-z(icrm,k).lt.damp_depth*z(icrm,nzm)) then
          n_damp(icrm)=nzm-k+1
        endif
      end do
      max_depth = max(max_depth,n_damp(icrm))
    end do

    !$acc parallel loop gang vector collapse(2)
    do k=nzm,nzm-max_depth,-1
      do icrm = 1 , ncrms
        if (k >= (nzm-n_damp(icrm))) then
          tau(icrm,k) = tau_min *(tau_max/tau_min)**((z(icrm,nzm)-z(icrm,k))/(z(icrm,nzm)-z(icrm,nzm-n_damp(icrm))))
          tau(icrm,k)=1./tau(icrm,k)
        endif
      enddo
    end do

    !+++mhwang recalculate grid-mean u0, v0, t0 first,
    ! as t have been updated. No need for qv0, as
    ! qv has not been updated yet the calculation of qv0.
    !$acc parallel loop gang vector collapse(2)
    do k=1, nzm
      do icrm = 1 , ncrms
        u0(icrm,k)=0.0
        v0(icrm,k)=0.0
        t0(icrm,k)=0.0
        !$acc loop seq
        do j=1, ny
          !$acc loop seq
          do i=1, nx
            u0(icrm,k) = u0(icrm,k) + u(icrm,i,j,k)/(nx*ny)
            v0(icrm,k) = v0(icrm,k) + v(icrm,i,j,k)/(nx*ny)
            t0(icrm,k) = t0(icrm,k) + t(icrm,i,j,k)/(nx*ny)
          end do
        end do
      end do
    end do
    !---mhwang

    !$acc parallel loop gang vector collapse(4)
    do k = nzm, nzm-max_depth, -1
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            if (k >= (nzm-n_damp(icrm))) then
              dudt(icrm,i,j,k,na)= dudt(icrm,i,j,k,na)-(u(icrm,i,j,k)-u0(icrm,k)) * tau(icrm,k)
              dvdt(icrm,i,j,k,na)= dvdt(icrm,i,j,k,na)-(v(icrm,i,j,k)-v0(icrm,k)) * tau(icrm,k)
              dwdt(icrm,i,j,k,na)= dwdt(icrm,i,j,k,na)-w(icrm,i,j,k) * tau(icrm,k)
              t(icrm,i,j,k)= t(icrm,i,j,k)-dtn*(t(icrm,i,j,k)-t0(icrm,k)) * tau(icrm,k)
              ! In the old version (SAM7.5?) of SAM, water vapor is the prognostic variable for the two-moment microphyscs.
              ! So the following damping approach can lead to the negative water vapor.
              !      micro_field(icrm,i,j,k,index_water_vapor)= micro_field(icrm,i,j,k,index_water_vapor)- &
              !                                    dtn*(qv(icrm,i,j,k)+qcl(icrm,i,j,k)+qci(icrm,i,j,k)-q0(icrm,k)) * tau(icrm,k)
              ! a simple fix (Minghuai Wang, 2011-08)icrm
              micro_field(icrm,i,j,k,index_water_vapor)= micro_field(icrm,i,j,k,index_water_vapor)- &
              dtn*(qv(icrm,i,j,k)-qv0(icrm,k)) * tau(icrm,k)
            endif
          enddo
        end do! i
      end do! j
    end do ! k


    ! deallocate(n_damp)
    ! deallocate(tau)
    call pool_pop_multiple(2)

  end subroutine damping

end module damping_mod
