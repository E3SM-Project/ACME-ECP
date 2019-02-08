
module damping_mod
  use params, only: asyncid
  use task_util_mod
  implicit none

contains

  subroutine damping(ncrms)
    !  "Spange"-layer damping at the domain top region
    use vars
    use microphysics, only: micro_field, index_water_vapor
    use params, only: crm_rknd
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) tau_min    ! minimum damping time-scale (at the top)
    real(crm_rknd) tau_max    ! maxim damping time-scale (base of damping layer)
    real(crm_rknd) damp_depth ! damping depth as a fraction of the domain height
    parameter(tau_min=60., tau_max=450., damp_depth=0.4)
    real(crm_rknd) tau(nzm,ncrms), tmp
    integer i, j, k, n_damp(ncrms), icrm
    integer :: numgangs  !For working around PGI OpenACC bug where it didn't create enough gangs
    ! crjones tests: make changes to u0, v0, t0 local instead of shared with vars
    real(crm_rknd) :: t0loc(nzm, ncrms), u0loc(nzm, ncrms), v0loc(nzm, ncrms)
   
    !$acc enter data create(tau,n_damp,t0loc, u0loc, v0loc) async(asyncid)

    if(tau_min.lt.2*dt) then
      print*,'Error: in damping() tau_min is too small!'
      call task_abort()
    end if

    !$acc parallel loop copyin(z) copyout(n_damp) async(asyncid)
    do icrm = 1 , ncrms
      do k=nzm,1,-1
        if(z(nzm,icrm)-z(k,icrm).lt.damp_depth*z(nzm,icrm)) then
          n_damp(icrm)=nzm-k+1
        endif
      end do
    end do

    !$acc parallel loop copyin(z,n_damp) copy(tau) async(asyncid)
    do icrm = 1 , ncrms
      do k=nzm,nzm-n_damp(icrm),-1
        tau(k,icrm) = tau_min *(tau_max/tau_min)**((z(nzm,icrm)-z(k,icrm))/(z(nzm,icrm)-z(nzm-n_damp(icrm),icrm)))
        tau(k,icrm)=1./tau(k,icrm)
      end do
    end do

    ! recalculate grid-mean u0, v0, t0 first,
    ! as t has been updated. No need for qv0, as
    ! qv has not been updated yet the calculation of qv0.
    !$acc parallel loop collapse(2) async(asyncid)
    do icrm = 1 , ncrms
      do k=1, nzm
        u0loc(k,icrm)=0.0
        v0loc(k,icrm)=0.0
        t0loc(k,icrm)=0.0
      end do
    end do
    !$acc parallel loop collapse(4) copyin(u,v,t) async(asyncid)
    do icrm = 1 , ncrms
      do k=1, nzm
        do j=1, ny
          do i=1, nx
            tmp = u(i,j,k,icrm)/(nx*ny)
            !$acc atomic update
            u0loc(k,icrm) = u0loc(k,icrm) + tmp
            tmp = v(i,j,k,icrm)/(nx*ny)
            !$acc atomic update
            v0loc(k,icrm) = v0loc(k,icrm) + tmp
            tmp = t(i,j,k,icrm)/(nx*ny)
            !$acc atomic update
            t0loc(k,icrm) = t0loc(k,icrm) + tmp
          end do
        end do
      end do
    end do

   !For working around PGI OpenACC bug where it didn't create enough gangs 
    numgangs = ceiling(ncrms*ny*nx/128.)
    !$acc parallel loop collapse(3) vector_length(128) num_gangs(numgangs) copy(dudt,dvdt,dwdt,t,micro_field) copyin(n_damp,u,v,tau,w,qv,qv0) async(asyncid)
    do icrm = 1 , ncrms
      do j=1,ny
        do i=1,nx
          do k = nzm, nzm-n_damp(icrm), -1
            dudt(i,j,k,na,icrm)= dudt(i,j,k,na,icrm)-(u(i,j,k,icrm)-u0loc(k,icrm)) * tau(k,icrm)
            dvdt(i,j,k,na,icrm)= dvdt(i,j,k,na,icrm)-(v(i,j,k,icrm)-v0loc(k,icrm)) * tau(k,icrm)
            dwdt(i,j,k,na,icrm)= dwdt(i,j,k,na,icrm)-w(i,j,k,icrm) * tau(k,icrm)
            t(i,j,k,icrm)= t(i,j,k,icrm)-dtn*(t(i,j,k,icrm)-t0loc(k,icrm)) * tau(k,icrm)
            micro_field(i,j,k,icrm,index_water_vapor)= micro_field(i,j,k,icrm,index_water_vapor)-dtn*(qv(i,j,k,icrm)-qv0(k,icrm)) * tau(k,icrm)
          end do! i
        end do! j
      end do ! k
    end do

    !$acc exit data delete(tau,n_damp,t0loc,u0loc,v0loc) async(asyncid)

  end subroutine damping

end module damping_mod
