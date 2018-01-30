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
          bet(icrm,k)*betu* &
          ( tabs0(icrm,k)*(epsv*(qv(icrm,i,j,k)-qv0(icrm,k))-(qcl(icrm,i,j,k)+qci(icrm,i,j,k)-qn0(icrm,k)+qpl(icrm,i,j,k)+qpi(icrm,i,j,k)-qp0(icrm,k))) &
          +(tabs(icrm,i,j,k)-tabs0(icrm,k))*(1.+epsv*qv0(icrm,k)-qn0(icrm,k)-qp0(icrm,k)) ) &
          + bet(icrm,kb)*betd* &
          ( tabs0(icrm,kb)*(epsv*(qv(icrm,i,j,kb)-qv0(icrm,kb))-(qcl(icrm,i,j,kb)+qci(icrm,i,j,kb)-qn0(icrm,kb)+qpl(icrm,i,j,kb)+qpi(icrm,i,j,kb)-qp0(icrm,kb))) &
          +(tabs(icrm,i,j,kb)-tabs0(icrm,kb))*(1.+epsv*qv0(icrm,kb)-qn0(icrm,kb)-qp0(icrm,kb)) )

        end do ! i
      end do ! j
    end do ! k

  end subroutine buoyancy

end module buoyancy_mod
