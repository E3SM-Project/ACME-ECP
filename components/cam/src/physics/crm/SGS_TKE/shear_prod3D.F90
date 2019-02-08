module shear_prod3D_mod
  use params, only: asyncid
  implicit none

contains

  subroutine shear_prod3D(ncrms,def2)
    use vars
    use params, only: crm_rknd
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) def2(nx,ny,nzm,ncrms)
    real(crm_rknd) rdx0,rdx,rdx_up,rdx_dn
    real(crm_rknd) rdy0,rdy,rdy_up,rdy_dn
    real(crm_rknd) rdz,rdzw_up,rdzw_dn
    integer i,j,k,ib,ic,jb,jc,kb,kc,icrm

    rdx0=1./dx
    rdy0=1./dy

    !$acc parallel loop collapse(4) copyin(v,u,adzw,w,v0,u0,dz,adz) copy(def2) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            if ( k >= 2 .and. k <= nzm-1 ) then
              kb=k-1
              kc=k+1
              rdz = 1./(dz(icrm)*adz(k,icrm))
              rdzw_up = 1./(dz(icrm)*adzw(icrm,kc))
              rdzw_dn = 1./(dz(icrm)*adzw(icrm,k))
              rdx=rdx0 * sqrt(dx*rdz) ! take into account grid anisotropy
              rdy=rdy0 * sqrt(dy*rdz)
              rdx_up=rdx0 * sqrt(dx*rdzw_up)
              rdy_up=rdy0 * sqrt(dy*rdzw_up)
              rdx_dn=rdx0 * sqrt(dx*rdzw_dn)
              rdy_dn=rdy0 * sqrt(dy*rdzw_dn)
              jb=j-YES3D
              jc=j+YES3D
              ib=i-1
              ic=i+1
              def2(i,j,k,icrm)=2.* ( &
              ( (u(ic,j,k,icrm)-u(i,j,k,icrm))*rdx)**2+ &
              ( (v(i,jc,k,icrm)-v(i,j,k,icrm))*rdy)**2+ &
              ( (w(i,j,kc,icrm)-w(i,j,k,icrm))*rdz)**2 ) &
              + 0.25 * ( &
              ( (u(ic,jc,k,icrm)-u(ic,j ,k,icrm))*rdy+(v(ic,jc,k,icrm)-v(i ,jc,k,icrm))*rdx )**2 +  &
              ( (u(i ,jc,k,icrm)-u(i ,j ,k,icrm))*rdy+(v(i ,jc,k,icrm)-v(ib,jc,k,icrm))*rdx )**2 +  &
              ( (u(ic,j ,k,icrm)-u(ic,jb,k,icrm))*rdy+(v(ic,j ,k,icrm)-v(i ,j ,k,icrm))*rdx )**2 +  &
              ( (u(i ,j ,k,icrm)-u(i ,jb,k,icrm))*rdy+(v(i ,j ,k,icrm)-v(ib,j ,k,icrm))*rdx )**2 )
              def2(i,j,k,icrm)=def2(i,j,k,icrm) &
              + 0.25 * ( &
              ( (u(ic,j,kc,icrm)-u0(kc,icrm)-u(ic,j, k,icrm)+u0(k,icrm))*rdzw_up+ &
              (w(ic,j,kc,icrm)-w(i ,j,kc,icrm))*rdx_up )**2 + &
              ( (u(i ,j,kc,icrm)-u0(kc,icrm)-u(i ,j, k,icrm)+u0(k,icrm))*rdzw_up+ &
              (w(i ,j,kc,icrm)-w(ib,j,kc,icrm))*rdx_up )**2 + &
              ( (u(ic,j,k ,icrm)-u0(k,icrm)-u(ic,j,kb,icrm)+u0(kb,icrm))*rdzw_dn+ &
              (w(ic,j,k ,icrm)-w(i ,j,k ,icrm))*rdx_dn )**2 + &
              ( (u(i ,j,k ,icrm)-u0(k,icrm)-u(i ,j,kb,icrm)+u0(kb,icrm))*rdzw_dn+ &
              (w(i ,j,k ,icrm)-w(ib,j,k ,icrm))*rdx_dn )**2 )
              def2(i,j,k,icrm)=def2(i,j,k,icrm) &
              + 0.25 * ( &
              ( (v(i,jc,kc,icrm)-v0(kc,icrm)-v(i,jc, k,icrm)+v0(k,icrm))*rdzw_up+ &
              (w(i,jc,kc,icrm)-w(i,j ,kc,icrm))*rdy_up )**2 + &
              ( (v(i,j ,kc,icrm)-v0(kc,icrm)-v(i,j , k,icrm)+v0(k,icrm))*rdzw_up+ &
              (w(i,j ,kc,icrm)-w(i,jb,kc,icrm))*rdy_up )**2 + &
              ( (v(i,jc,k ,icrm)-v0(k,icrm)-v(i,jc,kb,icrm)+v0(kb,icrm))*rdzw_dn+ &
              (w(i,jc,k ,icrm)-w(i,j ,k ,icrm))*rdy_dn )**2 + &
              ( (v(i,j ,k ,icrm)-v0(k,icrm)-v(i,j ,kb,icrm)+v0(kb,icrm))*rdzw_dn+ &
              (w(i,j ,k ,icrm)-w(i,jb,k ,icrm))*rdy_dn )**2 )
            elseif (k == 1) then
              kc=k+1
              rdz = 1./(dz(icrm)*adz(k,icrm))
              rdzw_up = 1./(dz(icrm)*adzw(icrm,kc))
              rdx=rdx0 * sqrt(dx*rdz) ! take into account grid anisotropy
              rdy=rdy0 * sqrt(dy*rdz)
              rdx_up=rdx0 * sqrt(dx*rdzw_up)
              rdy_up=rdy0 * sqrt(dy*rdzw_up)
              jb=j-YES3D
              jc=j+YES3D
              ib=i-1
              ic=i+1
              def2(i,j,k,icrm)=2.* ( &
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
            elseif (k == nzm) then
              kc=k+1
              kb=k-1
              rdz = 1./(dz(icrm)*adz(k,icrm))
              rdzw_dn = 1./(dz(icrm)*adzw(icrm,k))
              rdx=rdx0 * sqrt(dx*rdz) ! take into account grid anisotropy
              rdy=rdy0 * sqrt(dy*rdz)
              rdx_dn=rdx0 * sqrt(dx*rdzw_dn)
              rdy_dn=rdy0 * sqrt(dy*rdzw_dn)
              jb=j-1*YES3D
              jc=j+1*YES3D
              ib=i-1
              ic=i+1
              def2(i,j,k,icrm)=2.* ( &
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
            endif
          end do
        end do
      end do ! k
    end do

  end subroutine shear_prod3D

end module shear_prod3D_mod
