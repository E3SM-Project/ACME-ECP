module buoyancy_mod
  implicit none

contains

  subroutine buoyancy(ncrms,icrm)
    use vars
    use params
    implicit none
    integer, intent(in) :: ncrms,icrm
    integer i,j,k,kb
    real(crm_rknd) betu, betd

    if(docolumn) return

    do k=2,nzm
      kb=k-1
      betu=adz(kb)/(adz(k)+adz(kb))
      betd=adz(k)/(adz(k)+adz(kb))
      do j=1,ny
        do i=1,nx

          dwdt(i,j,k,na,icrm)=dwdt(i,j,k,na,icrm) +  &
          bet(k)*betu* &
          ( tabs0(k)*(epsv*(qv(i,j,k,icrm)-qv0(k))-(qcl(i,j,k,icrm)+qci(i,j,k,icrm)-qn0(k)+qpl(i,j,k,icrm)+qpi(i,j,k,icrm)-qp0(k))) &
          +(tabs(i,j,k,icrm)-tabs0(k))*(1.+epsv*qv0(k)-qn0(k)-qp0(k)) ) &
          + bet(kb)*betd* &
          ( tabs0(kb)*(epsv*(qv(i,j,kb,icrm)-qv0(kb))-(qcl(i,j,kb,icrm)+qci(i,j,kb,icrm)-qn0(kb)+qpl(i,j,kb,icrm)+qpi(i,j,kb,icrm)-qp0(kb))) &
          +(tabs(i,j,kb,icrm)-tabs0(kb))*(1.+epsv*qv0(kb)-qn0(kb)-qp0(kb)) )

        end do ! i
      end do ! j
    end do ! k

  end subroutine buoyancy

end module buoyancy_mod
