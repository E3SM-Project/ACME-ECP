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
            rdzw_up = 1./(dz(icrm)*adzw(kc,icrm))
            rdzw_dn = 1./(dz(icrm)*adzw(k,icrm))
            rdx=rdx0 * sqrt(dx*rdz) ! take into account grid anisotropy
            rdx_up=rdx0 * sqrt(dx*rdzw_up)
            rdx_dn=rdx0 * sqrt(dx*rdzw_dn)
            ib=i-1
            ic=i+1
            def2(i,j,k,icrm)=2.* ( &
            ( (u(ic,j,k,icrm)-u(i,j,k,icrm))*rdx)**2+ &
            ( (w(i,j,kc,icrm)-w(i,j,k,icrm))*rdz)**2 ) &
            + 0.5 * ( &
            ( (v(ic,j ,k,icrm)-v(i ,j ,k,icrm))*rdx )**2 +  &
            ( (v(i ,j ,k,icrm)-v(ib,j ,k,icrm))*rdx )**2 +   &
            ( (u(ic,j,kc,icrm)-u0(kc,icrm)-u(ic,j, k,icrm)+u0(k,icrm))*rdzw_up+ &
            (w(ic,j,kc,icrm)-w(i ,j,kc,icrm))*rdx_up )**2 + &
            ( (u(i ,j,kc,icrm)-u0(kc,icrm)-u(i ,j, k,icrm)+u0(k,icrm))*rdzw_up+ &
            (w(i ,j,kc,icrm)-w(ib,j,kc,icrm))*rdx_up )**2 + &
            ( (u(ic,j,k ,icrm)-u0(k,icrm)-u(ic,j,kb,icrm)+u0(kb,icrm))*rdzw_dn+ &
            (w(ic,j,k ,icrm)-w(i ,j,k ,icrm))*rdx_dn )**2 + &
            ( (u(i ,j,k ,icrm)-u0(k,icrm)-u(i ,j,kb,icrm)+u0(kb,icrm))*rdzw_dn+ &
            (w(i ,j,k ,icrm)-w(ib,j,k ,icrm))*rdx_dn )**2 + &
            ( (v(i,j ,kc,icrm)-v0(kc,icrm)-v(i,j , k,icrm)+v0(k,icrm))*rdzw_up )**2 + &
            ( (v(i,j ,k ,icrm)-v0(k,icrm)-v(i,j ,kb,icrm)+v0(kb,icrm))*rdzw_dn )**2 )
          elseif (k == 1) then
            kc=k+1
            rdz = 1./(dz(icrm)*adz(k,icrm))
            rdzw_up = 1./(dz(icrm)*adzw(kc,icrm))
            rdx=rdx0 * sqrt(dx*rdz) ! take into account grid anisotropy
            rdx_up=rdx0 * sqrt(dx*rdzw_up)
            ib=i-1
            ic=i+1
            def2(i,j,k,icrm)=2.* ( &
            ( (u(ic,j,k,icrm)-u(i,j,k,icrm))*rdx)**2+ &
            ( (w(i,j,kc,icrm)-w(i,j,k,icrm))*rdz)**2 ) &
            + 0.5 * ( &
            ( (v(ic,j ,k,icrm)-v(i ,j ,k,icrm))*rdx )**2 + &
            ( (v(i ,j ,k,icrm)-v(ib,j ,k,icrm))*rdx )**2 ) &
            +( (v(i,j ,kc,icrm)-v0(kc,icrm)-v(i,j,k,icrm)+v0(k,icrm))*rdzw_up )**2 &
            + 0.5 * ( &
            ( (u(ic,j,kc,icrm)-u0(kc,icrm)-u(ic,j, k,icrm)+u0(k,icrm))*rdzw_up+ &
            (w(ic,j,kc,icrm)-w(i ,j,kc,icrm))*rdx_up )**2 + &
            ( (u(i ,j,kc,icrm)-u0(kc,icrm)-u(i ,j, k,icrm)+u0(k,icrm))*rdzw_up+ &
            (w(i ,j,kc,icrm)-w(ib,j,kc,icrm))*rdx_up )**2 )
          elseif (k == nzm) then
            kc=k+1
            kb=k-1
            rdz = 1./(dz(icrm)*adz(k,icrm))
            rdzw_dn = 1./(dz(icrm)*adzw(k,icrm))
            rdx=rdx0 * sqrt(dx*rdz) ! take into account grid anisotropy
            rdx_dn=rdx0 * sqrt(dx*rdzw_dn)
            ib=i-1
            ic=i+1
            def2(i,j,k,icrm)=2.* ( &
            ( (u(ic,j,k,icrm)-u(i,j,k,icrm))*rdx)**2+ &
            ( (w(i,j,kc,icrm)-w(i,j,k,icrm))*rdz)**2 ) &
            + 0.5 * ( &
            ( (v(ic,j ,k,icrm)-v(i ,j ,k,icrm))*rdx )**2 +  &
            ( (v(i ,j ,k,icrm)-v(ib,j ,k,icrm))*rdx )**2 )   &
            + ( (v(i,j ,k ,icrm)-v0(k,icrm)-v(i,j ,kb,icrm)+v0(kb,icrm))*rdzw_dn )**2 &
            + 0.5 * ( &
            ( (u(ic,j,k ,icrm)-u0(k,icrm)-u(ic,j,kb,icrm)+u0(kb,icrm))*rdzw_dn+ &
            (w(ic,j,k ,icrm)-w(i ,j,k ,icrm))*rdx_dn )**2 + &
            ( (u(i ,j,k ,icrm)-u0(k,icrm)-u(i ,j,kb,icrm)+u0(kb,icrm))*rdzw_dn+ &
            (w(i ,j,k ,icrm)-w(ib,j,k ,icrm))*rdx_dn )**2 )
          endif
        end do
      end do ! k
    end do
  end subroutine shear_prod2D

end module shear_prod2D_mod
