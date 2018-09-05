module shear_prod2D_mod
  implicit none

contains

  subroutine shear_prod2D(ncrms,icrm,def2)

    use vars
    use params, only: crm_rknd
    implicit none
    integer, intent(in) :: ncrms,icrm
    real(crm_rknd) def2(nx,ny,nzm)

    real(crm_rknd) rdx0,rdx,rdx_up,rdx_dn
    real(crm_rknd) rdz,rdzw_up,rdzw_dn
    integer i,j,k,ib,ic,kb,kc

    rdx0=1./dx
    j=1


    do k=2,nzm-1

      kb=k-1
      kc=k+1
      rdz = 1./(dz(icrm)*adz(k,icrm))
      rdzw_up = 1./(dz(icrm)*adzw(kc,icrm))
      rdzw_dn = 1./(dz(icrm)*adzw(k,icrm))
      rdx=rdx0 * sqrt(dx*rdz) ! take into account grid anisotropy
      rdx_up=rdx0 * sqrt(dx*rdzw_up)
      rdx_dn=rdx0 * sqrt(dx*rdzw_dn)

      do i=1,nx
        ib=i-1
        ic=i+1

        def2(i,j,k)=2.* ( &
        ( (u(ic,j,k)-u(i,j,k))*rdx)**2+ &
        ( (w(i,j,kc)-w(i,j,k))*rdz)**2 ) &
        + 0.5 * ( &
        ( (v(ic,j ,k)-v(i ,j ,k))*rdx )**2 +  &
        ( (v(i ,j ,k)-v(ib,j ,k))*rdx )**2 +   &
        ( (u(ic,j,kc)-u0(kc,icrm)-u(ic,j, k)+u0(k,icrm))*rdzw_up+ &
        (w(ic,j,kc)-w(i ,j,kc))*rdx_up )**2 + &
        ( (u(i ,j,kc)-u0(kc,icrm)-u(i ,j, k)+u0(k,icrm))*rdzw_up+ &
        (w(i ,j,kc)-w(ib,j,kc))*rdx_up )**2 + &
        ( (u(ic,j,k )-u0(k,icrm)-u(ic,j,kb)+u0(kb,icrm))*rdzw_dn+ &
        (w(ic,j,k )-w(i ,j,k ))*rdx_dn )**2 + &
        ( (u(i ,j,k )-u0(k,icrm)-u(i ,j,kb)+u0(kb,icrm))*rdzw_dn+ &
        (w(i ,j,k )-w(ib,j,k ))*rdx_dn )**2 + &
        ( (v(i,j ,kc)-v0(kc,icrm)-v(i,j , k)+v0(k,icrm))*rdzw_up )**2 + &
        ( (v(i,j ,k )-v0(k,icrm)-v(i,j ,kb)+v0(kb,icrm))*rdzw_dn )**2 )

      end do
    end do ! k


    k=1
    kc=k+1

    rdz = 1./(dz(icrm)*adz(k,icrm))
    rdzw_up = 1./(dz(icrm)*adzw(kc,icrm))
    rdx=rdx0 * sqrt(dx*rdz) ! take into account grid anisotropy
    rdx_up=rdx0 * sqrt(dx*rdzw_up)

    do i=1,nx
      ib=i-1
      ic=i+1

      def2(i,j,k)=2.* ( &
      ( (u(ic,j,k)-u(i,j,k))*rdx)**2+ &
      ( (w(i,j,kc)-w(i,j,k))*rdz)**2 ) &
      + 0.5 * ( &
      ( (v(ic,j ,k)-v(i ,j ,k))*rdx )**2 + &
      ( (v(i ,j ,k)-v(ib,j ,k))*rdx )**2 ) &
      +( (v(i,j ,kc)-v0(kc,icrm)-v(i,j,k)+v0(k,icrm))*rdzw_up )**2 &
      + 0.5 * ( &
      ( (u(ic,j,kc)-u0(kc,icrm)-u(ic,j, k)+u0(k,icrm))*rdzw_up+ &
      (w(ic,j,kc)-w(i ,j,kc))*rdx_up )**2 + &
      ( (u(i ,j,kc)-u0(kc,icrm)-u(i ,j, k)+u0(k,icrm))*rdzw_up+ &
      (w(i ,j,kc)-w(ib,j,kc))*rdx_up )**2 )
    end do

    k=nzm
    kc=k+1
    kb=k-1

    rdz = 1./(dz(icrm)*adz(k,icrm))
    rdzw_dn = 1./(dz(icrm)*adzw(k,icrm))
    rdx=rdx0 * sqrt(dx*rdz) ! take into account grid anisotropy
    rdx_dn=rdx0 * sqrt(dx*rdzw_dn)


    do i=1,nx
      ib=i-1
      ic=i+1

      def2(i,j,k)=2.* ( &
      ( (u(ic,j,k)-u(i,j,k))*rdx)**2+ &
      ( (w(i,j,kc)-w(i,j,k))*rdz)**2 ) &
      + 0.5 * ( &
      ( (v(ic,j ,k)-v(i ,j ,k))*rdx )**2 +  &
      ( (v(i ,j ,k)-v(ib,j ,k))*rdx )**2 )   &
      + ( (v(i,j ,k )-v0(k,icrm)-v(i,j ,kb)+v0(kb,icrm))*rdzw_dn )**2 &
      + 0.5 * ( &
      ( (u(ic,j,k )-u0(k,icrm)-u(ic,j,kb)+u0(kb,icrm))*rdzw_dn+ &
      (w(ic,j,k )-w(i ,j,k ))*rdx_dn )**2 + &
      ( (u(i ,j,k )-u0(k,icrm)-u(i ,j,kb)+u0(kb,icrm))*rdzw_dn+ &
      (w(i ,j,k )-w(ib,j,k ))*rdx_dn )**2 )

    end do

  end

end module shear_prod2D_mod
