module diffuse_mom3D_mod
  implicit none

contains

  subroutine diffuse_mom3D(ncrms,icrm,grdf_x, grdf_y, grdf_z, dimx1_d, dimx2_d, dimy1_d, dimy2_d, tk)

    !        momentum tendency due to SGS diffusion

    use vars
    use params, only: docolumn, crm_rknd
    implicit none
    integer, intent(in) :: ncrms,icrm
    integer :: dimx1_d, dimx2_d, dimy1_d, dimy2_d
    real(crm_rknd) tk  (dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm,ncrms) ! SGS eddy viscosity
    real(crm_rknd) grdf_x(nzm,ncrms)! grid factor for eddy diffusion in x
    real(crm_rknd) grdf_y(nzm,ncrms)! grid factor for eddy diffusion in y
    real(crm_rknd) grdf_z(nzm,ncrms)! grid factor for eddy diffusion in z

    real(crm_rknd) rdx2,rdy2,rdz2,rdz,rdx25,rdy25
    real(crm_rknd) rdx21,rdy21,rdx251,rdy251,rdz25
    real(crm_rknd) dxy,dxz,dyx,dyz,dzx,dzy

    integer i,j,k,ic,ib,jb,jc,kc,kcu
    real(crm_rknd) tkx, tky, tkz, rhoi, iadzw, iadz
    real(crm_rknd) fu(0:nx,0:ny,nz),fv(0:nx,0:ny,nz),fw(0:nx,0:ny,nz)

    rdx2=1./(dx*dx)
    rdy2=1./(dy*dy)

    rdx25=0.25*rdx2
    rdy25=0.25*rdy2

    dxy=dx/dy
    dxz=dx/dz(icrm)
    dyx=dy/dx
    dyz=dy/dz(icrm)


    do k=1,nzm
      kc=k+1
      kcu=min(kc,nzm)
      dxz=dx/(dz(icrm)*adzw(kc,icrm))
      dyz=dy/(dz(icrm)*adzw(kc,icrm))
      rdx21=rdx2    * grdf_x(k,icrm)
      rdy21=rdy2    * grdf_y(k,icrm)
      rdx251=rdx25  * grdf_x(k,icrm)
      rdy251=rdy25  * grdf_y(k,icrm)
      do j=1,ny
        jb=j-1
        do i=0,nx
          ic=i+1
          tkx=rdx21*tk(i,j,k,icrm)
          fu(i,j,k)=-2.*tkx*(u(ic,j,k)-u(i,j,k))
          tkx=rdx251*(tk(i,j,k,icrm)+tk(i,jb,k,icrm)+tk(ic,j,k,icrm)+tk(ic,jb,k,icrm))
          fv(i,j,k)=-tkx*(v(ic,j,k)-v(i,j,k)+(u(ic,j,k)-u(ic,jb,k))*dxy)
          tkx=rdx251*(tk(i,j,k,icrm)+tk(ic,j,k,icrm)+tk(i,j,kcu,icrm)+tk(ic,j,kcu,icrm))
          fw(i,j,k)=-tkx*(w(ic,j,kc)-w(i,j,kc)+(u(ic,j,kcu)-u(ic,j,k))*dxz)
        end do
        do i=1,nx
          ib=i-1
          dudt(i,j,k,na,icrm)=dudt(i,j,k,na,icrm)-(fu(i,j,k)-fu(ib,j,k))
          dvdt(i,j,k,na,icrm)=dvdt(i,j,k,na,icrm)-(fv(i,j,k)-fv(ib,j,k))
          dwdt(i,j,kc,na,icrm)=dwdt(i,j,kc,na,icrm)-(fw(i,j,k)-fw(ib,j,k))
        end do
      end do

      do j=0,ny
        jc=j+1
        do i=1,nx
          ib=i-1
          tky=rdy21*tk(i,j,k,icrm)
          fv(i,j,k)=-2.*tky*(v(i,jc,k)-v(i,j,k))
          tky=rdy251*(tk(i,j,k,icrm)+tk(ib,j,k,icrm)+tk(i,jc,k,icrm)+tk(ib,jc,k,icrm))
          fu(i,j,k)=-tky*(u(i,jc,k)-u(i,j,k)+(v(i,jc,k)-v(ib,jc,k))*dyx)
          tky=rdy251*(tk(i,j,k,icrm)+tk(i,jc,k,icrm)+tk(i,j,kcu,icrm)+tk(i,jc,kcu,icrm))
          fw(i,j,k)=-tky*(w(i,jc,kc)-w(i,j,kc)+(v(i,jc,kcu)-v(i,jc,k))*dyz)
        end do
      end do
      do j=1,ny
        jb=j-1
        do i=1,nx
          dudt(i,j,k,na,icrm)=dudt(i,j,k,na,icrm)-(fu(i,j,k)-fu(i,jb,k))
          dvdt(i,j,k,na,icrm)=dvdt(i,j,k,na,icrm)-(fv(i,j,k)-fv(i,jb,k))
          dwdt(i,j,kc,na,icrm)=dwdt(i,j,kc,na,icrm)-(fw(i,j,k)-fw(i,jb,k))
        end do
      end do

    end do

    !-------------------------
    rdz=1./dz(icrm)
    dzx=dz(icrm)/dx
    dzy=dz(icrm)/dy

    do k=1,nzm-1
      kc=k+1
      uwsb(kc,icrm)=0.
      vwsb(kc,icrm)=0.
      iadz = 1./adz(k,icrm)
      iadzw= 1./adzw(kc,icrm)
      rdz2 = rdz*rdz * grdf_z(k,icrm)
      rdz25 = 0.25*rdz2
      do j=1,ny
        jb=j-1
        do i=1,nx
          ib=i-1
          tkz=rdz2*tk(i,j,k,icrm)
          fw(i,j,kc)=-2.*tkz*(w(i,j,kc)-w(i,j,k))*rho(k,icrm)*iadz
          tkz=rdz25*(tk(i,j,k,icrm)+tk(ib,j,k,icrm)+tk(i,j,kc,icrm)+tk(ib,j,kc,icrm))
          fu(i,j,kc)=-tkz*( (u(i,j,kc)-u(i,j,k))*iadzw + &
          (w(i,j,kc)-w(ib,j,kc))*dzx)*rhow(kc,icrm)
          tkz=rdz25*(tk(i,j,k,icrm)+tk(i,jb,k,icrm)+tk(i,j,kc,icrm)+tk(i,jb,kc,icrm))
          fv(i,j,kc)=-tkz*( (v(i,j,kc)-v(i,j,k))*iadzw + &
          (w(i,j,kc)-w(i,jb,kc))*dzy)*rhow(kc,icrm)
          uwsb(kc,icrm)=uwsb(kc,icrm)+fu(i,j,kc)
          vwsb(kc,icrm)=vwsb(kc,icrm)+fv(i,j,kc)
        end do
      end do
    end do

    uwsb(1,icrm) = 0.
    vwsb(1,icrm) = 0.

    do j=1,ny
      do i=1,nx
        tkz=rdz2*grdf_z(nzm,icrm)*tk(i,j,nzm,icrm)
        fw(i,j,nz)=-2.*tkz*(w(i,j,nz)-w(i,j,nzm))/adz(nzm,icrm)*rho(nzm,icrm)
        fu(i,j,1)=fluxbu(i,j,icrm) * rdz * rhow(1,icrm)
        fv(i,j,1)=fluxbv(i,j,icrm) * rdz * rhow(1,icrm)
        fu(i,j,nz)=fluxtu(i,j,icrm) * rdz * rhow(nz,icrm)
        fv(i,j,nz)=fluxtv(i,j,icrm) * rdz * rhow(nz,icrm)
        uwsb(1,icrm) = uwsb(1,icrm) + fu(i,j,1)
        vwsb(1,icrm) = vwsb(1,icrm) + fv(i,j,1)
      end do
    end do

    do k=1,nzm
      kc=k+1
      rhoi = 1./(rho(k,icrm)*adz(k,icrm))
      do j=1,ny
        do i=1,nx
          dudt(i,j,k,na,icrm)=dudt(i,j,k,na,icrm)-(fu(i,j,kc)-fu(i,j,k))*rhoi
          dvdt(i,j,k,na,icrm)=dvdt(i,j,k,na,icrm)-(fv(i,j,kc)-fv(i,j,k))*rhoi
        end do
      end do
    end do ! k

    do k=2,nzm
      rhoi = 1./(rhow(k,icrm)*adzw(k,icrm))
      do j=1,ny
        do i=1,nx
          dwdt(i,j,k,na,icrm)=dwdt(i,j,k,na,icrm)-(fw(i,j,k+1)-fw(i,j,k))*rhoi
        end do
      end do
    end do ! k


  end subroutine diffuse_mom3D

end module diffuse_mom3D_mod
