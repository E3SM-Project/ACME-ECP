module shear_prod2D_mod
  use params, only: asyncid
  implicit none

contains

  subroutine shear_prod2D(ncrms,def2)
    use vars
    use params, only: crm_rknd
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) def2(nx,ny,nzm,ncrms)
    real(crm_rknd) rdx0,rdx,rdx_up,rdx_dn
    real(crm_rknd) rdz,rdzw_up,rdzw_dn
    integer i,j,k,ib,ic,kb,kc,icrm

    rdx0=1./dx
    j=1

    !$acc parallel loop collapse(3) copyin(adz,w,v0,dz,u,v,u0,adzw) copy(def2) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        do i=1,nx
          if ( k >= 2 .and. k <= nzm-1) then
            kb=k-1
            kc=k+1
            rdz = 1./(dz(icrm)*adz(k,icrm))
            rdzw_up = 1./(dz(icrm)*adzw(icrm,kc))
            rdzw_dn = 1./(dz(icrm)*adzw(icrm,k))
            rdx=rdx0 * sqrt(dx*rdz) ! take into account grid anisotropy
            rdx_up=rdx0 * sqrt(dx*rdzw_up)
            rdx_dn=rdx0 * sqrt(dx*rdzw_dn)
            ib=i-1
            ic=i+1
            def2(i,j,k,icrm)=2.* ( &
            ( (u(icrm,ic,j,k)-u(icrm,i,j,k))*rdx)**2+ &
            ( (w(icrm,i,j,kc)-w(icrm,i,j,k))*rdz)**2 ) &
            + 0.5 * ( &
            ( (v(icrm,ic,j ,k)-v(icrm,i ,j ,k))*rdx )**2 +  &
            ( (v(icrm,i ,j ,k)-v(icrm,ib,j ,k))*rdx )**2 +   &
            ( (u(icrm,ic,j,kc)-u0(kc,icrm)-u(icrm,ic,j, k)+u0(k,icrm))*rdzw_up+ &
            (w(icrm,ic,j,kc)-w(icrm,i ,j,kc))*rdx_up )**2 + &
            ( (u(icrm,i ,j,kc)-u0(kc,icrm)-u(icrm,i ,j, k)+u0(k,icrm))*rdzw_up+ &
            (w(icrm,i ,j,kc)-w(icrm,ib,j,kc))*rdx_up )**2 + &
            ( (u(icrm,ic,j,k )-u0(k,icrm)-u(icrm,ic,j,kb)+u0(kb,icrm))*rdzw_dn+ &
            (w(icrm,ic,j,k )-w(icrm,i ,j,k ))*rdx_dn )**2 + &
            ( (u(icrm,i ,j,k )-u0(k,icrm)-u(icrm,i ,j,kb)+u0(kb,icrm))*rdzw_dn+ &
            (w(icrm,i ,j,k )-w(icrm,ib,j,k ))*rdx_dn )**2 + &
            ( (v(icrm,i,j ,kc)-v0(kc,icrm)-v(icrm,i,j , k)+v0(k,icrm))*rdzw_up )**2 + &
            ( (v(icrm,i,j ,k )-v0(k,icrm)-v(icrm,i,j ,kb)+v0(kb,icrm))*rdzw_dn )**2 )
          elseif (k == 1) then
            kc=k+1
            rdz = 1./(dz(icrm)*adz(k,icrm))
            rdzw_up = 1./(dz(icrm)*adzw(icrm,kc))
            rdx=rdx0 * sqrt(dx*rdz) ! take into account grid anisotropy
            rdx_up=rdx0 * sqrt(dx*rdzw_up)
            ib=i-1
            ic=i+1
            def2(i,j,k,icrm)=2.* ( &
            ( (u(icrm,ic,j,k)-u(icrm,i,j,k))*rdx)**2+ &
            ( (w(icrm,i,j,kc)-w(icrm,i,j,k))*rdz)**2 ) &
            + 0.5 * ( &
            ( (v(icrm,ic,j ,k)-v(icrm,i ,j ,k))*rdx )**2 + &
            ( (v(icrm,i ,j ,k)-v(icrm,ib,j ,k))*rdx )**2 ) &
            +( (v(icrm,i,j ,kc)-v0(kc,icrm)-v(icrm,i,j,k)+v0(k,icrm))*rdzw_up )**2 &
            + 0.5 * ( &
            ( (u(icrm,ic,j,kc)-u0(kc,icrm)-u(icrm,ic,j, k)+u0(k,icrm))*rdzw_up+ &
            (w(icrm,ic,j,kc)-w(icrm,i ,j,kc))*rdx_up )**2 + &
            ( (u(icrm,i ,j,kc)-u0(kc,icrm)-u(icrm,i ,j, k)+u0(k,icrm))*rdzw_up+ &
            (w(icrm,i ,j,kc)-w(icrm,ib,j,kc))*rdx_up )**2 )
          elseif (k == nzm) then
            kc=k+1
            kb=k-1
            rdz = 1./(dz(icrm)*adz(k,icrm))
            rdzw_dn = 1./(dz(icrm)*adzw(icrm,k))
            rdx=rdx0 * sqrt(dx*rdz) ! take into account grid anisotropy
            rdx_dn=rdx0 * sqrt(dx*rdzw_dn)
            ib=i-1
            ic=i+1
            def2(i,j,k,icrm)=2.* ( &
            ( (u(icrm,ic,j,k)-u(icrm,i,j,k))*rdx)**2+ &
            ( (w(icrm,i,j,kc)-w(icrm,i,j,k))*rdz)**2 ) &
            + 0.5 * ( &
            ( (v(icrm,ic,j ,k)-v(icrm,i ,j ,k))*rdx )**2 +  &
            ( (v(icrm,i ,j ,k)-v(icrm,ib,j ,k))*rdx )**2 )   &
            + ( (v(icrm,i,j ,k )-v0(k,icrm)-v(icrm,i,j ,kb)+v0(kb,icrm))*rdzw_dn )**2 &
            + 0.5 * ( &
            ( (u(icrm,ic,j,k )-u0(k,icrm)-u(icrm,ic,j,kb)+u0(kb,icrm))*rdzw_dn+ &
            (w(icrm,ic,j,k )-w(icrm,i ,j,k ))*rdx_dn )**2 + &
            ( (u(icrm,i ,j,k )-u0(k,icrm)-u(icrm,i ,j,kb)+u0(kb,icrm))*rdzw_dn+ &
            (w(icrm,i ,j,k )-w(icrm,ib,j,k ))*rdx_dn )**2 )
          endif
        end do
      end do ! k
    end do
  end subroutine shear_prod2D

end module shear_prod2D_mod
