module shear_prod2D_mod
  implicit none

contains

  subroutine shear_prod2D(def2,ncrms,icrm)

    use vars
    use params, only: crm_rknd
    implicit none
    integer, intent(in) :: ncrms, icrm

    real(crm_rknd) def2(nx,ny,nzm)

    real(crm_rknd) rdx0,rdx,rdx_up,rdx_dn
    real(crm_rknd) rdz,rdzw_up,rdzw_dn
    integer i,j,k,ib,ic,kb,kc

    rdx0=1./dx
    j=1


    do k=2,nzm-1

      kb=k-1
      kc=k+1
      rdz = 1./(dz*adz(k))
      rdzw_up = 1./(dz*adzw(kc))
      rdzw_dn = 1./(dz*adzw(k))
      rdx=rdx0 * sqrt(dx*rdz) ! take into account grid anisotropy
      rdx_up=rdx0 * sqrt(dx*rdzw_up)
      rdx_dn=rdx0 * sqrt(dx*rdzw_dn)

      do i=1,nx
        ib=i-1
        ic=i+1

        def2(i,j,k)=2.* ( &
        ( (u(icrm,ic,j,k)-u(icrm,i,j,k))*rdx)**2+ &
        ( (w(icrm,i,j,kc)-w(icrm,i,j,k))*rdz)**2 ) &
        + 0.5 * ( &
        ( (v(icrm,ic,j ,k)-v(icrm,i ,j ,k))*rdx )**2 +  &
        ( (v(icrm,i ,j ,k)-v(icrm,ib,j ,k))*rdx )**2 +   &
        ( (u(icrm,ic,j,kc)-u0(icrm,kc)-u(icrm,ic,j, k)+u0(icrm,k))*rdzw_up+ &
        (w(icrm,ic,j,kc)-w(icrm,i ,j,kc))*rdx_up )**2 + &
        ( (u(icrm,i ,j,kc)-u0(icrm,kc)-u(icrm,i ,j, k)+u0(icrm,k))*rdzw_up+ &
        (w(icrm,i ,j,kc)-w(icrm,ib,j,kc))*rdx_up )**2 + &
        ( (u(icrm,ic,j,k )-u0(icrm,k)-u(icrm,ic,j,kb)+u0(icrm,kb))*rdzw_dn+ &
        (w(icrm,ic,j,k )-w(icrm,i ,j,k ))*rdx_dn )**2 + &
        ( (u(icrm,i ,j,k )-u0(icrm,k)-u(icrm,i ,j,kb)+u0(icrm,kb))*rdzw_dn+ &
        (w(icrm,i ,j,k )-w(icrm,ib,j,k ))*rdx_dn )**2 + &
        ( (v(icrm,i,j ,kc)-v0(icrm,kc)-v(icrm,i,j , k)+v0(icrm,k))*rdzw_up )**2 + &
        ( (v(icrm,i,j ,k )-v0(icrm,k)-v(icrm,i,j ,kb)+v0(icrm,kb))*rdzw_dn )**2 )

      end do
    end do ! k


    k=1
    kc=k+1

    rdz = 1./(dz*adz(k))
    rdzw_up = 1./(dz*adzw(kc))
    rdx=rdx0 * sqrt(dx*rdz) ! take into account grid anisotropy
    rdx_up=rdx0 * sqrt(dx*rdzw_up)

    do i=1,nx
      ib=i-1
      ic=i+1

      def2(i,j,k)=2.* ( &
      ( (u(icrm,ic,j,k)-u(icrm,i,j,k))*rdx)**2+ &
      ( (w(icrm,i,j,kc)-w(icrm,i,j,k))*rdz)**2 ) &
      + 0.5 * ( &
      ( (v(icrm,ic,j ,k)-v(icrm,i ,j ,k))*rdx )**2 + &
      ( (v(icrm,i ,j ,k)-v(icrm,ib,j ,k))*rdx )**2 ) &
      +( (v(icrm,i,j ,kc)-v0(icrm,kc)-v(icrm,i,j,k)+v0(icrm,k))*rdzw_up )**2 &
      + 0.5 * ( &
      ( (u(icrm,ic,j,kc)-u0(icrm,kc)-u(icrm,ic,j, k)+u0(icrm,k))*rdzw_up+ &
      (w(icrm,ic,j,kc)-w(icrm,i ,j,kc))*rdx_up )**2 + &
      ( (u(icrm,i ,j,kc)-u0(icrm,kc)-u(icrm,i ,j, k)+u0(icrm,k))*rdzw_up+ &
      (w(icrm,i ,j,kc)-w(icrm,ib,j,kc))*rdx_up )**2 )
    end do

    k=nzm
    kc=k+1
    kb=k-1

    rdz = 1./(dz*adz(k))
    rdzw_dn = 1./(dz*adzw(k))
    rdx=rdx0 * sqrt(dx*rdz) ! take into account grid anisotropy
    rdx_dn=rdx0 * sqrt(dx*rdzw_dn)


    do i=1,nx
      ib=i-1
      ic=i+1

      def2(i,j,k)=2.* ( &
      ( (u(icrm,ic,j,k)-u(icrm,i,j,k))*rdx)**2+ &
      ( (w(icrm,i,j,kc)-w(icrm,i,j,k))*rdz)**2 ) &
      + 0.5 * ( &
      ( (v(icrm,ic,j ,k)-v(icrm,i ,j ,k))*rdx )**2 +  &
      ( (v(icrm,i ,j ,k)-v(icrm,ib,j ,k))*rdx )**2 )   &
      + ( (v(icrm,i,j ,k )-v0(icrm,k)-v(icrm,i,j ,kb)+v0(icrm,kb))*rdzw_dn )**2 &
      + 0.5 * ( &
      ( (u(icrm,ic,j,k )-u0(icrm,k)-u(icrm,ic,j,kb)+u0(icrm,kb))*rdzw_dn+ &
      (w(icrm,ic,j,k )-w(icrm,i ,j,k ))*rdx_dn )**2 + &
      ( (u(icrm,i ,j,k )-u0(icrm,k)-u(icrm,i ,j,kb)+u0(icrm,kb))*rdzw_dn+ &
      (w(icrm,i ,j,k )-w(icrm,ib,j,k ))*rdx_dn )**2 )

    end do

  end

end module shear_prod2D_mod
