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

          dwdt(icrm,i,j,k,na)=dwdt(icrm,i,j,k,na) +  &
          bet(k)*betu* &
          ( tabs0(k)*(epsv*(qv(icrm,i,j,k)-qv0(k))-(qcl(icrm,i,j,k)+qci(icrm,i,j,k)-qn0(k)+qpl(icrm,i,j,k)+qpi(icrm,i,j,k)-qp0(k))) &
          +(tabs(icrm,i,j,k)-tabs0(k))*(1.+epsv*qv0(k)-qn0(k)-qp0(k)) ) &
          + bet(kb)*betd* &
          ( tabs0(kb)*(epsv*(qv(icrm,i,j,kb)-qv0(kb))-(qcl(icrm,i,j,kb)+qci(icrm,i,j,kb)-qn0(kb)+qpl(icrm,i,j,kb)+qpi(icrm,i,j,kb)-qp0(kb))) &
          +(tabs(icrm,i,j,kb)-tabs0(kb))*(1.+epsv*qv0(kb)-qn0(kb)-qp0(kb)) )

        end do ! i
      end do ! j
    end do ! k

  end subroutine buoyancy

end module buoyancy_mod
