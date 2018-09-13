module coriolis_mod
  implicit none

contains

  subroutine coriolis(ncrms,icrm)

    use vars
    use params, only: crm_rknd

    implicit none
    integer, intent(in) :: ncrms,icrm
    real(crm_rknd) u_av, v_av, w_av
    integer i,j,k,ib,ic,jb,jc,kc

    if(RUN3D) then

      do k=1,nzm
        kc=k+1
        do j=1,ny
          jb=j-1
          jc=j+1
          do i=1,nx
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

    else

      do k=1,nzm
        kc=k+1
        do j=1,ny
          do i=1,nx
            ib=i-1
            ic=i+1
            w_av=0.25*(w(i,j,kc,icrm)+w(ib,j,kc,icrm)+w(i,j,k,icrm)+w(ib,j,k,icrm))
            dudt(i,j,k,na,icrm)=dudt(i,j,k,na,icrm)+fcory(j,icrm)*(v(i,j,k,icrm)-vg0(k,icrm))-fcorzy(j,icrm)*w_av
            dvdt(i,j,k,na,icrm)=dvdt(i,j,k,na,icrm)-fcory(j,icrm)*(u(i,j,k,icrm)-ug0(k,icrm))
          end do ! i
        end do ! i
      end do ! k

    endif

  end subroutine coriolis

end module coriolis_mod
