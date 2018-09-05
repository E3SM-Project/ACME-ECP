module shear_prod3D_mod
  implicit none

contains

  subroutine shear_prod3D(ncrms,icrm,def2)

    use vars
    use params, only: crm_rknd
    implicit none
    integer, intent(in) :: ncrms,icrm
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
      rdz = 1./(dz(icrm)*adz(k,icrm))
      rdzw_up = 1./(dz(icrm)*adzw(kc,icrm))
      rdzw_dn = 1./(dz(icrm)*adzw(k,icrm))
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
          ( (u(ic,j,k,icrm)-u(i,j,k,icrm))*rdx)**2+ &
          ( (v(i,jc,k,icrm)-v(i,j,k,icrm))*rdy)**2+ &
          ( (w(i,j,kc,icrm)-w(i,j,k,icrm))*rdz)**2 ) &
          + 0.25 * ( &
          ( (u(ic,jc,k,icrm)-u(ic,j ,k,icrm))*rdy+(v(ic,jc,k,icrm)-v(i ,jc,k,icrm))*rdx )**2 +  &
          ( (u(i ,jc,k,icrm)-u(i ,j ,k,icrm))*rdy+(v(i ,jc,k,icrm)-v(ib,jc,k,icrm))*rdx )**2 +  &
          ( (u(ic,j ,k,icrm)-u(ic,jb,k,icrm))*rdy+(v(ic,j ,k,icrm)-v(i ,j ,k,icrm))*rdx )**2 +  &
          ( (u(i ,j ,k,icrm)-u(i ,jb,k,icrm))*rdy+(v(i ,j ,k,icrm)-v(ib,j ,k,icrm))*rdx )**2 )
          def2(i,j,k)=def2(i,j,k) &
          + 0.25 * ( &
          ( (u(ic,j,kc,icrm)-u0(kc,icrm)-u(ic,j, k,icrm)+u0(k,icrm))*rdzw_up+ &
          (w(ic,j,kc,icrm)-w(i ,j,kc,icrm))*rdx_up )**2 + &
          ( (u(i ,j,kc,icrm)-u0(kc,icrm)-u(i ,j, k,icrm)+u0(k,icrm))*rdzw_up+ &
          (w(i ,j,kc,icrm)-w(ib,j,kc,icrm))*rdx_up )**2 + &
          ( (u(ic,j,k ,icrm)-u0(k,icrm)-u(ic,j,kb,icrm)+u0(kb,icrm))*rdzw_dn+ &
          (w(ic,j,k ,icrm)-w(i ,j,k ,icrm))*rdx_dn )**2 + &
          ( (u(i ,j,k ,icrm)-u0(k,icrm)-u(i ,j,kb,icrm)+u0(kb,icrm))*rdzw_dn+ &
          (w(i ,j,k ,icrm)-w(ib,j,k ,icrm))*rdx_dn )**2 )
          def2(i,j,k)=def2(i,j,k) &
          + 0.25 * ( &
          ( (v(i,jc,kc,icrm)-v0(kc,icrm)-v(i,jc, k,icrm)+v0(k,icrm))*rdzw_up+ &
          (w(i,jc,kc,icrm)-w(i,j ,kc,icrm))*rdy_up )**2 + &
          ( (v(i,j ,kc,icrm)-v0(kc,icrm)-v(i,j , k,icrm)+v0(k,icrm))*rdzw_up+ &
          (w(i,j ,kc,icrm)-w(i,jb,kc,icrm))*rdy_up )**2 + &
          ( (v(i,jc,k ,icrm)-v0(k,icrm)-v(i,jc,kb,icrm)+v0(kb,icrm))*rdzw_dn+ &
          (w(i,jc,k ,icrm)-w(i,j ,k ,icrm))*rdy_dn )**2 + &
          ( (v(i,j ,k ,icrm)-v0(k,icrm)-v(i,j ,kb,icrm)+v0(kb,icrm))*rdzw_dn+ &
          (w(i,j ,k ,icrm)-w(i,jb,k ,icrm))*rdy_dn )**2 )

        end do
      end do
    end do ! k


    k=1
    kc=k+1

    rdz = 1./(dz(icrm)*adz(k,icrm))
    rdzw_up = 1./(dz(icrm)*adzw(kc,icrm))
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
        ( (u(ic,j,k,icrm)-u(i,j,k,icrm))*rdx)**2+ &
        ( (v(i,jc,k,icrm)-v(i,j,k,icrm))*rdy)**2+ &
        ( (w(i,j,kc,icrm)-w(i,j,k,icrm))*rdz)**2 ) &
        + 0.25 * ( &
        ( (u(ic,jc,k,icrm)-u(ic,j ,k,icrm))*rdy+(v(ic,jc,k,icrm)-v(i ,jc,k,icrm))*rdx )**2 +  &
        ( (u(i ,jc,k,icrm)-u(i ,j ,k,icrm))*rdy+(v(i ,jc,k,icrm)-v(ib,jc,k,icrm))*rdx )**2 +  &
        ( (u(ic,j ,k,icrm)-u(ic,jb,k,icrm))*rdy+(v(ic,j ,k,icrm)-v(i ,j ,k,icrm))*rdx )**2 +  &
        ( (u(i ,j ,k,icrm)-u(i ,jb,k,icrm))*rdy+(v(i ,j ,k,icrm)-v(ib,j ,k,icrm))*rdx )**2 )   &
        + 0.5 * ( &
        ( (v(i,jc,kc,icrm)-v0(kc,icrm)-v(i,jc, k,icrm)+v0(k,icrm))*rdzw_up+ &
        (w(i,jc,kc,icrm)-w(i,j ,kc,icrm))*rdy_up )**2 + &
        ( (v(i,j ,kc,icrm)-v0(kc,icrm)-v(i,j , k,icrm)+v0(k,icrm))*rdzw_up+ &
        (w(i,j ,kc,icrm)-w(i,jb,kc,icrm))*rdy_up )**2 ) &
        + 0.5 * ( &
        ( (u(ic,j,kc,icrm)-u0(kc,icrm)-u(ic,j, k,icrm)+u0(k,icrm))*rdzw_up+ &
        (w(ic,j,kc,icrm)-w(i ,j,kc,icrm))*rdx_up )**2 + &
        ( (u(i ,j,kc,icrm)-u0(kc,icrm)-u(i ,j, k,icrm)+u0(k,icrm))*rdzw_up+ &
        (w(i ,j,kc,icrm)-w(ib,j,kc,icrm))*rdx_up )**2 )


      end do
    end do


    k=nzm
    kc=k+1
    kb=k-1

    rdz = 1./(dz(icrm)*adz(k,icrm))
    rdzw_dn = 1./(dz(icrm)*adzw(k,icrm))
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
        ( (u(ic,j,k,icrm)-u(i,j,k,icrm))*rdx)**2+ &
        ( (v(i,jc,k,icrm)-v(i,j,k,icrm))*rdy)**2+ &
        ( (w(i,j,kc,icrm)-w(i,j,k,icrm))*rdz)**2 ) &
        + 0.25 * ( &
        ( (u(ic,jc,k,icrm)-u(ic,j ,k,icrm))*rdy+(v(ic,jc,k,icrm)-v(i ,jc,k,icrm))*rdx )**2 +  &
        ( (u(i ,jc,k,icrm)-u(i ,j ,k,icrm))*rdy+(v(i ,jc,k,icrm)-v(ib,jc,k,icrm))*rdx )**2 +  &
        ( (u(ic,j ,k,icrm)-u(ic,jb,k,icrm))*rdy+(v(ic,j ,k,icrm)-v(i ,j ,k,icrm))*rdx )**2 +  &
        ( (u(i ,j ,k,icrm)-u(i ,jb,k,icrm))*rdy+(v(i ,j ,k,icrm)-v(ib,j ,k,icrm))*rdx )**2 )   &
        + 0.5 * ( &
        ( (v(i,jc,k ,icrm)-v0(k,icrm)-v(i,jc,kb,icrm)+v0(kb,icrm))*rdzw_dn+ &
        (w(i,jc,k ,icrm)-w(i,j ,k ,icrm))*rdy_dn )**2 + &
        ( (v(i,j ,k ,icrm)-v0(k,icrm)-v(i,j ,kb,icrm)+v0(kb,icrm))*rdzw_dn+ &
        (w(i,j ,k ,icrm)-w(i,jb,k ,icrm))*rdy_dn )**2 ) &
        + 0.5 * ( &
        ( (u(ic,j,k ,icrm)-u0(k,icrm)-u(ic,j,kb,icrm)+u0(kb,icrm))*rdzw_dn+ &
        (w(ic,j,k ,icrm)-w(i ,j,k ,icrm))*rdx_dn )**2 + &
        ( (u(i ,j,k ,icrm)-u0(k,icrm)-u(i ,j,kb,icrm)+u0(kb,icrm))*rdzw_dn+ &
        (w(i ,j,k ,icrm)-w(ib,j,k ,icrm))*rdx_dn )**2 )
      end do
    end do

  end

end module shear_prod3D_mod
