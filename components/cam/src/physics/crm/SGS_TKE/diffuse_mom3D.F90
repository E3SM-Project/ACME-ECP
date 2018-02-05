module diffuse_mom3D_mod
  implicit none

contains

  subroutine diffuse_mom3D(grdf_x, grdf_y, grdf_z, dimx1_d, dimx2_d, dimy1_d, dimy2_d, tk, ncrms)

    !        momentum tendency due to SGS diffusion

    use vars
    use params, only: docolumn, crm_rknd
    implicit none
    integer, intent(in) :: ncrms
    integer :: dimx1_d, dimx2_d, dimy1_d, dimy2_d
    real(crm_rknd) tk  (ncrms,dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm) ! SGS eddy viscosity
    real(crm_rknd) grdf_x(ncrms,nzm)! grid factor for eddy diffusion in x
    real(crm_rknd) grdf_y(ncrms,nzm)! grid factor for eddy diffusion in y
    real(crm_rknd) grdf_z(ncrms,nzm)! grid factor for eddy diffusion in z

    real(crm_rknd) rdx2,rdy2,rdz2(ncrms),rdz(ncrms),rdx25,rdy25
    real(crm_rknd) rdx21(ncrms),rdy21(ncrms),rdx251(ncrms),rdy251(ncrms),rdz25(ncrms)
    real(crm_rknd) dxy,dxz(ncrms),dyx,dyz(ncrms),dzx(ncrms),dzy(ncrms)

    integer i,j,k,ic,ib,jb,jc,kc,kcu, icrm
    real(crm_rknd) tkx, tky, tkz, rhoi(ncrms), iadzw(ncrms), iadz(ncrms)
    real(crm_rknd) fu(ncrms,0:nx,0:ny,nz),fv(ncrms,0:nx,0:ny,nz),fw(ncrms,0:nx,0:ny,nz)


    rdx2=1./(dx*dx)
    rdy2=1./(dy*dy)

    rdx25=0.25*rdx2
    rdy25=0.25*rdy2

    dxy=dx/dy
    dxz(:)=dx/dz(:)
    dyx=dy/dx
    dyz(:)=dy/dz(:)


    do k=1,nzm
      kc=k+1
      kcu=min(kc,nzm)
      dxz(:)=dx/(dz(:)*adzw(:,kc))
      dyz(:)=dy/(dz(:)*adzw(:,kc))
      rdx21(:)=rdx2    * grdf_x(:,k)
      rdy21(:)=rdy2    * grdf_y(:,k)
      rdx251(:)=rdx25  * grdf_x(:,k)
      rdy251(:)=rdy25  * grdf_y(:,k)
      do j=1,ny
        jb=j-1
        do i=0,nx
          do icrm = 1 , ncrms
            ic=i+1
            tkx=rdx21(icrm)*tk(icrm,i,j,k)
            fu(icrm,i,j,k)=-2.*tkx*(u(icrm,ic,j,k)-u(icrm,i,j,k))
            tkx=rdx251(icrm)*(tk(icrm,i,j,k)+tk(icrm,i,jb,k)+tk(icrm,ic,j,k)+tk(icrm,ic,jb,k))
            fv(icrm,i,j,k)=-tkx*(v(icrm,ic,j,k)-v(icrm,i,j,k)+(u(icrm,ic,j,k)-u(icrm,ic,jb,k))*dxy)
            tkx=rdx251(icrm)*(tk(icrm,i,j,k)+tk(icrm,ic,j,k)+tk(icrm,i,j,kcu)+tk(icrm,ic,j,kcu))
            fw(icrm,i,j,k)=-tkx*(w(icrm,ic,j,kc)-w(icrm,i,j,kc)+(u(icrm,ic,j,kcu)-u(icrm,ic,j,k))*dxz(icrm))
          end do
        end do
        do i=1,nx
          do icrm = 1 , ncrms
            ib=i-1
            dudt(icrm,i,j,k,na)=dudt(icrm,i,j,k,na)-(fu(icrm,i,j,k)-fu(icrm,ib,j,k))
            dvdt(icrm,i,j,k,na)=dvdt(icrm,i,j,k,na)-(fv(icrm,i,j,k)-fv(icrm,ib,j,k))
            dwdt(icrm,i,j,kc,na)=dwdt(icrm,i,j,kc,na)-(fw(icrm,i,j,k)-fw(icrm,ib,j,k))
          enddo
        end do
      end do

      do j=0,ny
        jc=j+1
        do i=1,nx
          do icrm = 1 , ncrms
            ib=i-1
            tky=rdy21(icrm)*tk(icrm,i,j,k)
            fv(icrm,i,j,k)=-2.*tky*(v(icrm,i,jc,k)-v(icrm,i,j,k))
            tky=rdy251(icrm)*(tk(icrm,i,j,k)+tk(icrm,ib,j,k)+tk(icrm,i,jc,k)+tk(icrm,ib,jc,k))
            fu(icrm,i,j,k)=-tky*(u(icrm,i,jc,k)-u(icrm,i,j,k)+(v(icrm,i,jc,k)-v(icrm,ib,jc,k))*dyx)
            tky=rdy251(icrm)*(tk(icrm,i,j,k)+tk(icrm,i,jc,k)+tk(icrm,i,j,kcu)+tk(icrm,i,jc,kcu))
            fw(icrm,i,j,k)=-tky*(w(icrm,i,jc,kc)-w(icrm,i,j,kc)+(v(icrm,i,jc,kcu)-v(icrm,i,jc,k))*dyz(icrm))
          enddo
        end do
      end do
      do j=1,ny
        jb=j-1
        do i=1,nx
          do icrm = 1 , ncrms
            dudt(icrm,i,j,k,na)=dudt(icrm,i,j,k,na)-(fu(icrm,i,j,k)-fu(icrm,i,jb,k))
            dvdt(icrm,i,j,k,na)=dvdt(icrm,i,j,k,na)-(fv(icrm,i,j,k)-fv(icrm,i,jb,k))
            dwdt(icrm,i,j,kc,na)=dwdt(icrm,i,j,kc,na)-(fw(icrm,i,j,k)-fw(icrm,i,jb,k))
          enddo
        end do
      end do

    end do

    !-------------------------
    rdz(:)=1./dz(:)
    dzx(:)=dz(:)/dx
    dzy(:)=dz(:)/dy

    do k=1,nzm-1
      kc=k+1
      uwsb(:,kc)=0.
      vwsb(:,kc)=0.
      iadz(:) = 1./adz(:,k)
      iadzw(:)= 1./adzw(:,kc)
      rdz2(:) = rdz(:)*rdz(:) * grdf_z(:,k)
      rdz25(:) = 0.25*rdz2(:)
      do j=1,ny
        jb=j-1
        do i=1,nx
          do icrm = 1 , ncrms
            ib=i-1
            tkz=rdz2(icrm)*tk(icrm,i,j,k)
            fw(icrm,i,j,kc)=-2.*tkz*(w(icrm,i,j,kc)-w(icrm,i,j,k))*rho(icrm,k)*iadz(icrm)
            tkz=rdz25(icrm)*(tk(icrm,i,j,k)+tk(icrm,ib,j,k)+tk(icrm,i,j,kc)+tk(icrm,ib,j,kc))
            fu(icrm,i,j,kc)=-tkz*( (u(icrm,i,j,kc)-u(icrm,i,j,k))*iadzw(icrm) + &
            (w(icrm,i,j,kc)-w(icrm,ib,j,kc))*dzx(icrm))*rhow(icrm,kc)
            tkz=rdz25(icrm)*(tk(icrm,i,j,k)+tk(icrm,i,jb,k)+tk(icrm,i,j,kc)+tk(icrm,i,jb,kc))
            fv(icrm,i,j,kc)=-tkz*( (v(icrm,i,j,kc)-v(icrm,i,j,k))*iadzw(icrm) + &
            (w(icrm,i,j,kc)-w(icrm,i,jb,kc))*dzy(icrm))*rhow(icrm,kc)
            uwsb(icrm,kc)=uwsb(icrm,kc)+fu(icrm,i,j,kc)
            vwsb(icrm,kc)=vwsb(icrm,kc)+fv(icrm,i,j,kc)
          enddo
        end do
      end do
    end do

    uwsb(:,1) = 0.
    vwsb(:,1) = 0.

    do j=1,ny
      do i=1,nx
        do icrm = 1 , ncrms
          tkz=rdz2(icrm)*grdf_z(icrm,nzm)*tk(icrm,i,j,nzm)
          fw(icrm,i,j,nz)=-2.*tkz*(w(icrm,i,j,nz)-w(icrm,i,j,nzm))/adz(icrm,nzm)*rho(icrm,nzm)
          fu(icrm,i,j,1)=fluxbu(icrm,i,j) * rdz(icrm) * rhow(icrm,1)
          fv(icrm,i,j,1)=fluxbv(icrm,i,j) * rdz(icrm) * rhow(icrm,1)
          fu(icrm,i,j,nz)=fluxtu(icrm,i,j) * rdz(icrm) * rhow(icrm,nz)
          fv(icrm,i,j,nz)=fluxtv(icrm,i,j) * rdz(icrm) * rhow(icrm,nz)
          uwsb(icrm,1) = uwsb(icrm,1) + fu(icrm,i,j,1)
          vwsb(icrm,1) = vwsb(icrm,1) + fv(icrm,i,j,1)
        enddo
      end do
    end do

    do k=1,nzm
      kc=k+1
      rhoi(:) = 1./(rho(:,k)*adz(:,k))
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            dudt(icrm,i,j,k,na)=dudt(icrm,i,j,k,na)-(fu(icrm,i,j,kc)-fu(icrm,i,j,k))*rhoi(icrm)
            dvdt(icrm,i,j,k,na)=dvdt(icrm,i,j,k,na)-(fv(icrm,i,j,kc)-fv(icrm,i,j,k))*rhoi(icrm)
          enddo
        end do
      end do
    end do ! k

    do k=2,nzm
      rhoi(:) = 1./(rhow(:,k)*adzw(:,k))
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            dwdt(icrm,i,j,k,na)=dwdt(icrm,i,j,k,na)-(fw(icrm,i,j,k+1)-fw(icrm,i,j,k))*rhoi(icrm)
          enddo
        end do
      end do
    end do ! k


  end subroutine diffuse_mom3D

end module diffuse_mom3D_mod
