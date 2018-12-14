module diffuse_mom3D_mod
  implicit none

contains

  subroutine diffuse_mom3D(ncrms,grdf_x, grdf_y, grdf_z, dimx1_d, dimx2_d, dimy1_d, dimy2_d, tk)

    !        momentum tendency due to SGS diffusion

    use vars
    use params, only: docolumn, crm_rknd
    implicit none
    integer, intent(in) :: ncrms
    integer :: dimx1_d, dimx2_d, dimy1_d, dimy2_d
    real(crm_rknd) tk  (dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm,ncrms) ! SGS eddy viscosity
    real(crm_rknd) grdf_x(nzm,ncrms)! grid factor for eddy diffusion in x
    real(crm_rknd) grdf_y(nzm,ncrms)! grid factor for eddy diffusion in y
    real(crm_rknd) grdf_z(nzm,ncrms)! grid factor for eddy diffusion in z
    real(crm_rknd) rdx2,rdy2,rdz2,rdz,rdx25,rdy25
    real(crm_rknd) rdx21,rdy21,rdx251,rdy251,rdz25
    real(crm_rknd) dxy,dxz,dyx,dyz,dzx,dzy
    integer i,j,k,ic,ib,jb,jc,kc,kcu,icrm
    real(crm_rknd) tkx, tky, tkz, rhoi, iadzw, iadz
    real(crm_rknd) :: fu(0:nx,0:ny,nz,ncrms)
    real(crm_rknd) :: fv(0:nx,0:ny,nz,ncrms)
    real(crm_rknd) :: fw(0:nx,0:ny,nz,ncrms)

    !$acc enter data create(fu,fv,fw) async(1)

    rdx2=1./(dx*dx)
    rdy2=1./(dy*dy)
    rdx25=0.25*rdx2
    rdy25=0.25*rdy2
    dxy=dx/dy
    dyx=dy/dx

    !$acc parallel loop collapse(4) copyin(w,v,u,grdf_x,dz,tk,adzw) copy(fv,fu,fw) async(1)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          do i=0,nx
            jb=j-1
            kc=k+1
            kcu=min(kc,nzm)
            ic=i+1
            dxz=dx/(dz(icrm)*adzw(kc,icrm))
            rdx21=rdx2    * grdf_x(k,icrm)
            rdx251=rdx25  * grdf_x(k,icrm)
            tkx=rdx21*tk(i,j,k,icrm)
            fu(i,j,k,icrm)=-2.*tkx*(u(ic,j,k,icrm)-u(i,j,k,icrm))
            tkx=rdx251*(tk(i,j,k,icrm)+tk(i,jb,k,icrm)+tk(ic,j,k,icrm)+tk(ic,jb,k,icrm))
            fv(i,j,k,icrm)=-tkx*(v(ic,j,k,icrm)-v(i,j,k,icrm)+(u(ic,j,k,icrm)-u(ic,jb,k,icrm))*dxy)
            tkx=rdx251*(tk(i,j,k,icrm)+tk(ic,j,k,icrm)+tk(i,j,kcu,icrm)+tk(ic,j,kcu,icrm))
            fw(i,j,k,icrm)=-tkx*(w(ic,j,kc,icrm)-w(i,j,kc,icrm)+(u(ic,j,kcu,icrm)-u(ic,j,k,icrm))*dxz)
          enddo
        enddo
      enddo
    enddo
    !$acc parallel loop collapse(4) copyin(fv,fw,fu) copy(dvdt,dudt,dwdt) async(1)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            kc=k+1
            ib=i-1
            dudt(i,j,k,na,icrm)=dudt(i,j,k,na,icrm)-(fu(i,j,k,icrm)-fu(ib,j,k,icrm))
            dvdt(i,j,k,na,icrm)=dvdt(i,j,k,na,icrm)-(fv(i,j,k,icrm)-fv(ib,j,k,icrm))
            dwdt(i,j,kc,na,icrm)=dwdt(i,j,kc,na,icrm)-(fw(i,j,k,icrm)-fw(ib,j,k,icrm))
          enddo
        enddo
      enddo
    enddo

    !$acc parallel loop collapse(4) copyin(adzw,tk,dz,grdf_y,u,w,v) copy(fw,fu,fv) async(1)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=0,ny
          do i=1,nx
            jc=j+1
            kc=k+1
            kcu=min(kc,nzm)
            ib=i-1
            dyz=dy/(dz(icrm)*adzw(kc,icrm))
            rdy21=rdy2    * grdf_y(k,icrm)
            rdy251=rdy25  * grdf_y(k,icrm)
            tky=rdy21*tk(i,j,k,icrm)
            fv(i,j,k,icrm)=-2.*tky*(v(i,jc,k,icrm)-v(i,j,k,icrm))
            tky=rdy251*(tk(i,j,k,icrm)+tk(ib,j,k,icrm)+tk(i,jc,k,icrm)+tk(ib,jc,k,icrm))
            fu(i,j,k,icrm)=-tky*(u(i,jc,k,icrm)-u(i,j,k,icrm)+(v(i,jc,k,icrm)-v(ib,jc,k,icrm))*dyx)
            tky=rdy251*(tk(i,j,k,icrm)+tk(i,jc,k,icrm)+tk(i,j,kcu,icrm)+tk(i,jc,kcu,icrm))
            fw(i,j,k,icrm)=-tky*(w(i,jc,kc,icrm)-w(i,j,kc,icrm)+(v(i,jc,kcu,icrm)-v(i,jc,k,icrm))*dyz)
          enddo
        enddo
      enddo
    enddo
    !$acc parallel loop collapse(4) copyin(fu,fw,fv) copy(dwdt,dudt,dvdt) async(1)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            jb=j-1
            kc=k+1
            dudt(i,j,k,na,icrm)=dudt(i,j,k,na,icrm)-(fu(i,j,k,icrm)-fu(i,jb,k,icrm))
            dvdt(i,j,k,na,icrm)=dvdt(i,j,k,na,icrm)-(fv(i,j,k,icrm)-fv(i,jb,k,icrm))
            dwdt(i,j,kc,na,icrm)=dwdt(i,j,kc,na,icrm)-(fw(i,j,k,icrm)-fw(i,jb,k,icrm))
          enddo
        enddo
      enddo
    enddo

    !$acc parallel loop collapse(2) copy(uwsb,vwsb) async(1)
    do icrm = 1 , ncrms
      do k = 1 , nzm
        uwsb(k,icrm)=0.
        vwsb(k,icrm)=0.
      enddo
    enddo

    !-------------------------
    !$acc parallel loop collapse(4) copyin(v,dz,tk,rho,rhow,w,grdf_z,adz,adzw,u) copy(fu,uwsb,fv,vwsb,fw) async(1)
    do icrm = 1 , ncrms
      do k=1,nzm-1
        do j=1,ny
          do i=1,nx
            jb=j-1
            kc=k+1
            ib=i-1
            rdz=1./dz(icrm)
            rdz2 = rdz*rdz * grdf_z(k,icrm)
            rdz25 = 0.25*rdz2
            iadz = 1./adz(k,icrm)
            iadzw= 1./adzw(kc,icrm)
            dzx=dz(icrm)/dx
            dzy=dz(icrm)/dy
            tkz=rdz2*tk(i,j,k,icrm)
            fw(i,j,kc,icrm)=-2.*tkz*(w(i,j,kc,icrm)-w(i,j,k,icrm))*rho(k,icrm)*iadz
            tkz=rdz25*(tk(i,j,k,icrm)+tk(ib,j,k,icrm)+tk(i,j,kc,icrm)+tk(ib,j,kc,icrm))
            fu(i,j,kc,icrm)=-tkz*( (u(i,j,kc,icrm)-u(i,j,k,icrm))*iadzw + (w(i,j,kc,icrm)-w(ib,j,kc,icrm))*dzx)*rhow(kc,icrm)
            tkz=rdz25*(tk(i,j,k,icrm)+tk(i,jb,k,icrm)+tk(i,j,kc,icrm)+tk(i,jb,kc,icrm))
            fv(i,j,kc,icrm)=-tkz*( (v(i,j,kc,icrm)-v(i,j,k,icrm))*iadzw + (w(i,j,kc,icrm)-w(i,jb,kc,icrm))*dzy)*rhow(kc,icrm)
            !$acc atomic update
            uwsb(kc,icrm)=uwsb(kc,icrm)+fu(i,j,kc,icrm)
            !$acc atomic update
            vwsb(kc,icrm)=vwsb(kc,icrm)+fv(i,j,kc,icrm)
          enddo
        enddo
      enddo
    enddo

    !$acc parallel loop collapse(3) copyin(tk,adz,rhow,grdf_z,w,rho,dz,fluxtv,fluxbu,fluxbv,fluxtu) copy(fw,vwsb,fv,fu,uwsb) async(1)
    do icrm = 1 , ncrms
      do j=1,ny
        do i=1,nx
          rdz=1./dz(icrm)
          rdz2 = rdz*rdz * grdf_z(k,icrm)
          tkz=rdz2*grdf_z(nzm,icrm)*tk(i,j,nzm,icrm)
          fw(i,j,nz,icrm)=-2.*tkz*(w(i,j,nz,icrm)-w(i,j,nzm,icrm))/adz(nzm,icrm)*rho(nzm,icrm)
          fu(i,j,1,icrm)=fluxbu(i,j,icrm) * rdz * rhow(1,icrm)
          fv(i,j,1,icrm)=fluxbv(i,j,icrm) * rdz * rhow(1,icrm)
          fu(i,j,nz,icrm)=fluxtu(i,j,icrm) * rdz * rhow(nz,icrm)
          fv(i,j,nz,icrm)=fluxtv(i,j,icrm) * rdz * rhow(nz,icrm)
          !$acc atomic update
          uwsb(1,icrm) = uwsb(1,icrm) + fu(i,j,1,icrm)
          !$acc atomic update
          vwsb(1,icrm) = vwsb(1,icrm) + fv(i,j,1,icrm)
        enddo
      enddo
    enddo

    !$acc parallel loop collapse(4) copyin(rho,fv,adz,fu) copy(dvdt,dudt) async(1)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            kc=k+1
            rhoi = 1./(rho(k,icrm)*adz(k,icrm))
            dudt(i,j,k,na,icrm)=dudt(i,j,k,na,icrm)-(fu(i,j,kc,icrm)-fu(i,j,k,icrm))*rhoi
            dvdt(i,j,k,na,icrm)=dvdt(i,j,k,na,icrm)-(fv(i,j,kc,icrm)-fv(i,j,k,icrm))*rhoi
          enddo
        enddo
      enddo ! k
    enddo

    !$acc parallel loop collapse(4) copyin(adzw,rhow,fw) copy(dwdt) async(1)
    do icrm = 1 , ncrms
      do k=2,nzm
        do j=1,ny
          do i=1,nx
            rhoi = 1./(rhow(k,icrm)*adzw(k,icrm))
            dwdt(i,j,k,na,icrm)=dwdt(i,j,k,na,icrm)-(fw(i,j,k+1,icrm)-fw(i,j,k,icrm))*rhoi
          enddo
        enddo
      enddo ! k
    enddo

    !$acc exit data delete(fu,fv,fw) async(1)

  end subroutine diffuse_mom3D

end module diffuse_mom3D_mod
