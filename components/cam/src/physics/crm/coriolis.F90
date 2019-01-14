module coriolis_mod
  use params, only: asyncid
  implicit none

contains

  subroutine coriolis(ncrms)
    use vars
    use params, only: crm_rknd
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) u_av, v_av, w_av
    integer i,j,k,ib,ic,jb,jc,kc,icrm

    if(RUN3D) then
      !$acc parallel loop collapse(4) copyin(w,vg0,ug0,v,fcory,fcorzy,u) copy(dvdt,dudt) async(asyncid)
      do icrm = 1 , ncrms
        do k=1,nzm
          do j=1,ny
            do i=1,nx
              kc=k+1
              jb=j-1
              jc=j+1
              ib=i-1
              ic=i+1
              v_av=0.25*(v(i,j,k,icrm)+v(i,jc,k,icrm)+v(ib,j,k,icrm)+v(ib,jc,k,icrm))
              w_av=0.25*(w(i,j,kc,icrm)+w(ib,j,kc,icrm)+w(i,j,k,icrm)+w(ib,j,k,icrm))
              dudt(i,j,k,na,icrm)=dudt(i,j,k,na,icrm)+fcory(j,icrm)*(v_av-vg0(k,icrm))-fcorzy(j,icrm)*w_av
              u_av=0.25*(u(i,j,k,icrm)+u(ic,j,k,icrm)+u(i,jb,k,icrm)+u(ic,jb,k,icrm))
              dvdt(i,j,k,na,icrm)=dvdt(i,j,k,na,icrm)-0.5*(fcory(j,icrm)+fcory(jb,icrm))*(u_av-ug0(k,icrm))
            end do ! i
          end do ! j
        end do ! k
      end do
    else
      !$acc parallel loop collapse(4) copyin(u,v,w,ug0,vg0,fcorzy,fcory) copy(dudt,dvdt) async(asyncid)
      do icrm = 1 , ncrms
        do k=1,nzm
          do j=1,ny
            do i=1,nx
              kc=k+1
              ib=i-1
              ic=i+1
              w_av=0.25*(w(i,j,kc,icrm)+w(ib,j,kc,icrm)+w(i,j,k,icrm)+w(ib,j,k,icrm))
              dudt(i,j,k,na,icrm)=dudt(i,j,k,na,icrm)+fcory(j,icrm)*(v(i,j,k,icrm)-vg0(k,icrm))-fcorzy(j,icrm)*w_av
              dvdt(i,j,k,na,icrm)=dvdt(i,j,k,na,icrm)-fcory(j,icrm)*(u(i,j,k,icrm)-ug0(k,icrm))
            end do ! i
          end do ! i
        end do ! k
      end do
    endif

  end subroutine coriolis

end module coriolis_mod
