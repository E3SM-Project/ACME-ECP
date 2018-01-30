module shear_prod3D_mod
  implicit none

contains

  subroutine shear_prod3D(def2,ncrms,icrm)

    use vars
    use params, only: crm_rknd
    implicit none
    integer, intent(in) :: ncrms, icrm

    real(crm_rknd) def2(nx,ny,nzm)

    real(crm_rknd) rdx0,rdx,rdx_up,rdx_dn
    real(crm_rknd) rdy0,rdy,rdy_up,rdy_dn
    real(crm_rknd) rdz,rdzw_up,rdzw_dn
    integer i,j,k,ib,ic,jb,jc,kb,kc

    rdx0=1./dx
    rdy0=1./dy

    do k=2,nzm-1

      kb=k-1
      kc=k+1
      rdz = 1./(dz*adz(k))
      rdzw_up = 1./(dz*adzw(kc))
      rdzw_dn = 1./(dz*adzw(k))
      rdx=rdx0 * sqrt(dx*rdz) ! take into account grid anisotropy
      rdy=rdy0 * sqrt(dy*rdz)
      rdx_up=rdx0 * sqrt(dx*rdzw_up)
      rdy_up=rdy0 * sqrt(dy*rdzw_up)
      rdx_dn=rdx0 * sqrt(dx*rdzw_dn)
      rdy_dn=rdy0 * sqrt(dy*rdzw_dn)

      do j=1,ny
        jb=j-YES3D
        jc=j+YES3D
        do i=1,nx
          ib=i-1
          ic=i+1

          def2(i,j,k)=2.* ( &
          ( (u(icrm,ic,j,k)-u(icrm,i,j,k))*rdx)**2+ &
          ( (v(icrm,i,jc,k)-v(icrm,i,j,k))*rdy)**2+ &
          ( (w(i,j,kc)-w(i,j,k))*rdz)**2 ) &
          + 0.25 * ( &
          ( (u(icrm,ic,jc,k)-u(icrm,ic,j ,k))*rdy+(v(icrm,ic,jc,k)-v(icrm,i ,jc,k))*rdx )**2 +  &
          ( (u(icrm,i ,jc,k)-u(icrm,i ,j ,k))*rdy+(v(icrm,i ,jc,k)-v(icrm,ib,jc,k))*rdx )**2 +  &
          ( (u(icrm,ic,j ,k)-u(icrm,ic,jb,k))*rdy+(v(icrm,ic,j ,k)-v(icrm,i ,j ,k))*rdx )**2 +  &
          ( (u(icrm,i ,j ,k)-u(icrm,i ,jb,k))*rdy+(v(icrm,i ,j ,k)-v(icrm,ib,j ,k))*rdx )**2 )
          def2(i,j,k)=def2(i,j,k) &
          + 0.25 * ( &
          ( (u(icrm,ic,j,kc)-u0(kc)-u(icrm,ic,j, k)+u0(k))*rdzw_up+ &
          (w(ic,j,kc)-w(i ,j,kc))*rdx_up )**2 + &
          ( (u(icrm,i ,j,kc)-u0(kc)-u(icrm,i ,j, k)+u0(k))*rdzw_up+ &
          (w(i ,j,kc)-w(ib,j,kc))*rdx_up )**2 + &
          ( (u(icrm,ic,j,k )-u0(k)-u(icrm,ic,j,kb)+u0(kb))*rdzw_dn+ &
          (w(ic,j,k )-w(i ,j,k ))*rdx_dn )**2 + &
          ( (u(icrm,i ,j,k )-u0(k)-u(icrm,i ,j,kb)+u0(kb))*rdzw_dn+ &
          (w(i ,j,k )-w(ib,j,k ))*rdx_dn )**2 )
          def2(i,j,k)=def2(i,j,k) &
          + 0.25 * ( &
          ( (v(icrm,i,jc,kc)-v0(kc)-v(icrm,i,jc, k)+v0(k))*rdzw_up+ &
          (w(i,jc,kc)-w(i,j ,kc))*rdy_up )**2 + &
          ( (v(icrm,i,j ,kc)-v0(kc)-v(icrm,i,j , k)+v0(k))*rdzw_up+ &
          (w(i,j ,kc)-w(i,jb,kc))*rdy_up )**2 + &
          ( (v(icrm,i,jc,k )-v0(k)-v(icrm,i,jc,kb)+v0(kb))*rdzw_dn+ &
          (w(i,jc,k )-w(i,j ,k ))*rdy_dn )**2 + &
          ( (v(icrm,i,j ,k )-v0(k)-v(icrm,i,j ,kb)+v0(kb))*rdzw_dn+ &
          (w(i,j ,k )-w(i,jb,k ))*rdy_dn )**2 )

        end do
      end do
    end do ! k


    k=1
    kc=k+1

    rdz = 1./(dz*adz(k))
    rdzw_up = 1./(dz*adzw(kc))
    rdx=rdx0 * sqrt(dx*rdz) ! take into account grid anisotropy
    rdy=rdy0 * sqrt(dy*rdz)
    rdx_up=rdx0 * sqrt(dx*rdzw_up)
    rdy_up=rdy0 * sqrt(dy*rdzw_up)

    do j=1,ny
      jb=j-YES3D
      jc=j+YES3D
      do i=1,nx
        ib=i-1
        ic=i+1

        def2(i,j,k)=2.* ( &
        ( (u(icrm,ic,j,k)-u(icrm,i,j,k))*rdx)**2+ &
        ( (v(icrm,i,jc,k)-v(icrm,i,j,k))*rdy)**2+ &
        ( (w(i,j,kc)-w(i,j,k))*rdz)**2 ) &
        + 0.25 * ( &
        ( (u(icrm,ic,jc,k)-u(icrm,ic,j ,k))*rdy+(v(icrm,ic,jc,k)-v(icrm,i ,jc,k))*rdx )**2 +  &
        ( (u(icrm,i ,jc,k)-u(icrm,i ,j ,k))*rdy+(v(icrm,i ,jc,k)-v(icrm,ib,jc,k))*rdx )**2 +  &
        ( (u(icrm,ic,j ,k)-u(icrm,ic,jb,k))*rdy+(v(icrm,ic,j ,k)-v(icrm,i ,j ,k))*rdx )**2 +  &
        ( (u(icrm,i ,j ,k)-u(icrm,i ,jb,k))*rdy+(v(icrm,i ,j ,k)-v(icrm,ib,j ,k))*rdx )**2 )   &
        + 0.5 * ( &
        ( (v(icrm,i,jc,kc)-v0(kc)-v(icrm,i,jc, k)+v0(k))*rdzw_up+ &
        (w(i,jc,kc)-w(i,j ,kc))*rdy_up )**2 + &
        ( (v(icrm,i,j ,kc)-v0(kc)-v(icrm,i,j , k)+v0(k))*rdzw_up+ &
        (w(i,j ,kc)-w(i,jb,kc))*rdy_up )**2 ) &
        + 0.5 * ( &
        ( (u(icrm,ic,j,kc)-u0(kc)-u(icrm,ic,j, k)+u0(k))*rdzw_up+ &
        (w(ic,j,kc)-w(i ,j,kc))*rdx_up )**2 + &
        ( (u(icrm,i ,j,kc)-u0(kc)-u(icrm,i ,j, k)+u0(k))*rdzw_up+ &
        (w(i ,j,kc)-w(ib,j,kc))*rdx_up )**2 )


      end do
    end do


    k=nzm
    kc=k+1
    kb=k-1

    rdz = 1./(dz*adz(k))
    rdzw_dn = 1./(dz*adzw(k))
    rdx=rdx0 * sqrt(dx*rdz) ! take into account grid anisotropy
    rdy=rdy0 * sqrt(dy*rdz)
    rdx_dn=rdx0 * sqrt(dx*rdzw_dn)
    rdy_dn=rdy0 * sqrt(dy*rdzw_dn)

    do j=1,ny
      jb=j-1*YES3D
      jc=j+1*YES3D
      do i=1,nx
        ib=i-1
        ic=i+1
        def2(i,j,k)=2.* ( &
        ( (u(icrm,ic,j,k)-u(icrm,i,j,k))*rdx)**2+ &
        ( (v(icrm,i,jc,k)-v(icrm,i,j,k))*rdy)**2+ &
        ( (w(i,j,kc)-w(i,j,k))*rdz)**2 ) &
        + 0.25 * ( &
        ( (u(icrm,ic,jc,k)-u(icrm,ic,j ,k))*rdy+(v(icrm,ic,jc,k)-v(icrm,i ,jc,k))*rdx )**2 +  &
        ( (u(icrm,i ,jc,k)-u(icrm,i ,j ,k))*rdy+(v(icrm,i ,jc,k)-v(icrm,ib,jc,k))*rdx )**2 +  &
        ( (u(icrm,ic,j ,k)-u(icrm,ic,jb,k))*rdy+(v(icrm,ic,j ,k)-v(icrm,i ,j ,k))*rdx )**2 +  &
        ( (u(icrm,i ,j ,k)-u(icrm,i ,jb,k))*rdy+(v(icrm,i ,j ,k)-v(icrm,ib,j ,k))*rdx )**2 )   &
        + 0.5 * ( &
        ( (v(icrm,i,jc,k )-v0(k)-v(icrm,i,jc,kb)+v0(kb))*rdzw_dn+ &
        (w(i,jc,k )-w(i,j ,k ))*rdy_dn )**2 + &
        ( (v(icrm,i,j ,k )-v0(k)-v(icrm,i,j ,kb)+v0(kb))*rdzw_dn+ &
        (w(i,j ,k )-w(i,jb,k ))*rdy_dn )**2 ) &
        + 0.5 * ( &
        ( (u(icrm,ic,j,k )-u0(k)-u(icrm,ic,j,kb)+u0(kb))*rdzw_dn+ &
        (w(ic,j,k )-w(i ,j,k ))*rdx_dn )**2 + &
        ( (u(icrm,i ,j,k )-u0(k)-u(icrm,i ,j,kb)+u0(kb))*rdzw_dn+ &
        (w(i ,j,k )-w(ib,j,k ))*rdx_dn )**2 )
      end do
    end do

  end

end module shear_prod3D_mod
