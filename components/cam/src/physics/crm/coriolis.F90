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
            v_av=0.25*(v(icrm,i,j,k)+v(icrm,i,jc,k)+v(icrm,ib,j,k)+v(icrm,ib,jc,k))
            w_av=0.25*(w(i,j,kc)+w(ib,j,kc)+w(i,j,k)+w(ib,j,k))
            dudt(i,j,k,na)=dudt(i,j,k,na)+fcory(j)*(v_av-vg0(k))-fcorzy(j)*w_av
            u_av=0.25*(u(icrm,i,j,k)+u(icrm,ic,j,k)+u(icrm,i,jb,k)+u(icrm,ic,jb,k))
            dvdt(i,j,k,na)=dvdt(i,j,k,na)-0.5*(fcory(j)+fcory(jb))*(u_av-ug0(k))
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
            w_av=0.25*(w(i,j,kc)+w(ib,j,kc)+w(i,j,k)+w(ib,j,k))
            dudt(i,j,k,na)=dudt(i,j,k,na)+fcory(j)*(v(icrm,i,j,k)-vg0(k))-fcorzy(j)*w_av
            dvdt(i,j,k,na)=dvdt(i,j,k,na)-fcory(j)*(u(icrm,i,j,k)-ug0(k))
          end do ! i
        end do ! i
      end do ! k

    endif

  end subroutine coriolis

end module coriolis_mod
