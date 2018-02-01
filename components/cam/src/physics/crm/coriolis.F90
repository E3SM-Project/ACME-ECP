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
            w_av=0.25*(w(icrm,i,j,kc)+w(icrm,ib,j,kc)+w(icrm,i,j,k)+w(icrm,ib,j,k))
            dudt(icrm,i,j,k,na(icrm))=dudt(icrm,i,j,k,na(icrm))+fcory(icrm,j)*(v_av-vg0(icrm,k))-fcorzy(icrm,j)*w_av
            u_av=0.25*(u(icrm,i,j,k)+u(icrm,ic,j,k)+u(icrm,i,jb,k)+u(icrm,ic,jb,k))
            dvdt(icrm,i,j,k,na(icrm))=dvdt(icrm,i,j,k,na(icrm))-0.5*(fcory(icrm,j)+fcory(icrm,jb))*(u_av-ug0(icrm,k))
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
            w_av=0.25*(w(icrm,i,j,kc)+w(icrm,ib,j,kc)+w(icrm,i,j,k)+w(icrm,ib,j,k))
            dudt(icrm,i,j,k,na(icrm))=dudt(icrm,i,j,k,na(icrm))+fcory(icrm,j)*(v(icrm,i,j,k)-vg0(icrm,k))-fcorzy(icrm,j)*w_av
            dvdt(icrm,i,j,k,na(icrm))=dvdt(icrm,i,j,k,na(icrm))-fcory(icrm,j)*(u(icrm,i,j,k)-ug0(icrm,k))
          end do ! i
        end do ! i
      end do ! k

    endif

  end subroutine coriolis

end module coriolis_mod
