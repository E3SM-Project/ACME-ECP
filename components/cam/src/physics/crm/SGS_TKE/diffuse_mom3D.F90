module diffuse_mom3D_mod
  use params, only: asyncid
  implicit none

contains

  subroutine diffuse_mom3D(ncrms,grdf_x, grdf_y, grdf_z, dimx1_d, dimx2_d, dimy1_d, dimy2_d, tk)

    !        momentum tendency due to SGS diffusion

    use vars
    use params, only: docolumn, crm_rknd
    implicit none
    integer, intent(in) :: ncrms
    integer :: dimx1_d, dimx2_d, dimy1_d, dimy2_d
    real(crm_rknd) tk(ncrms,dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm) ! SGS eddy viscosity
    real(crm_rknd) grdf_x(ncrms,nzm)! grid factor for eddy diffusion in x
    real(crm_rknd) grdf_y(ncrms,nzm)! grid factor for eddy diffusion in y
    real(crm_rknd) grdf_z(ncrms,nzm)! grid factor for eddy diffusion in z
    real(crm_rknd) rdx2,rdy2,rdz2,rdz,rdx25,rdy25
    real(crm_rknd) rdx21,rdy21,rdx251,rdy251,rdz25
    real(crm_rknd) dxy,dxz,dyx,dyz,dzx,dzy
    integer i,j,k,ic,ib,jb,jc,kc,kcu,icrm
    real(crm_rknd) tkx, tky, tkz, rhoi, iadzw, iadz
    real(crm_rknd) :: fu(0:nx,0:ny,nz,ncrms)
    real(crm_rknd) :: fv(0:nx,0:ny,nz,ncrms)
    real(crm_rknd) :: fw(0:nx,0:ny,nz,ncrms)

    !$acc enter data create(fu,fv,fw) async(asyncid)

    rdx2=1./(dx*dx)
    rdy2=1./(dy*dy)
    rdx25=0.25*rdx2
    rdy25=0.25*rdy2
    dxy=dx/dy
    dyx=dy/dx

    !$acc parallel loop collapse(4) copyin(w,v,u,grdf_x,dz,tk,adzw) copy(fv,fu,fw) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          do i=0,nx
            jb=j-1
            kc=k+1
            kcu=min(kc,nzm)
            ic=i+1
            dxz=dx/(dz(icrm)*adzw(icrm,kc))
            rdx21=rdx2    * grdf_x(icrm,k)
            rdx251=rdx25  * grdf_x(icrm,k)
            tkx=rdx21*tk(icrm,i,j,k)
            fu(i,j,k,icrm)=-2.*tkx*(u(icrm,ic,j,k)-u(icrm,i,j,k))
            tkx=rdx251*(tk(icrm,i,j,k)+tk(icrm,i,jb,k)+tk(icrm,ic,j,k)+tk(icrm,ic,jb,k))
            fv(i,j,k,icrm)=-tkx*(v(icrm,ic,j,k)-v(icrm,i,j,k)+(u(icrm,ic,j,k)-u(icrm,ic,jb,k))*dxy)
            tkx=rdx251*(tk(icrm,i,j,k)+tk(icrm,ic,j,k)+tk(icrm,i,j,kcu)+tk(icrm,ic,j,kcu))
            fw(i,j,k,icrm)=-tkx*(w(icrm,ic,j,kc)-w(icrm,i,j,kc)+(u(icrm,ic,j,kcu)-u(icrm,ic,j,k))*dxz)
          enddo
        enddo
      enddo
    enddo
    !$acc parallel loop collapse(4) copyin(fv,fw,fu) copy(dvdt,dudt,dwdt) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            kc=k+1
            ib=i-1
            dudt(icrm,i,j,k,na)=dudt(icrm,i,j,k,na)-(fu(i,j,k,icrm)-fu(ib,j,k,icrm))
            dvdt(icrm,i,j,k,na)=dvdt(icrm,i,j,k,na)-(fv(i,j,k,icrm)-fv(ib,j,k,icrm))
            dwdt(icrm,i,j,kc,na)=dwdt(icrm,i,j,kc,na)-(fw(i,j,k,icrm)-fw(ib,j,k,icrm))
          enddo
        enddo
      enddo
    enddo

    !$acc parallel loop collapse(4) copyin(adzw,tk,dz,grdf_y,u,w,v) copy(fw,fu,fv) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=0,ny
          do i=1,nx
            jc=j+1
            kc=k+1
            kcu=min(kc,nzm)
            ib=i-1
            dyz=dy/(dz(icrm)*adzw(icrm,kc))
            rdy21=rdy2    * grdf_y(icrm,k)
            rdy251=rdy25  * grdf_y(icrm,k)
            tky=rdy21*tk(icrm,i,j,k)
            fv(i,j,k,icrm)=-2.*tky*(v(icrm,i,jc,k)-v(icrm,i,j,k))
            tky=rdy251*(tk(icrm,i,j,k)+tk(icrm,ib,j,k)+tk(icrm,i,jc,k)+tk(icrm,ib,jc,k))
            fu(i,j,k,icrm)=-tky*(u(icrm,i,jc,k)-u(icrm,i,j,k)+(v(icrm,i,jc,k)-v(icrm,ib,jc,k))*dyx)
            tky=rdy251*(tk(icrm,i,j,k)+tk(icrm,i,jc,k)+tk(icrm,i,j,kcu)+tk(icrm,i,jc,kcu))
            fw(i,j,k,icrm)=-tky*(w(icrm,i,jc,kc)-w(icrm,i,j,kc)+(v(icrm,i,jc,kcu)-v(icrm,i,jc,k))*dyz)
          enddo
        enddo
      enddo
    enddo
    !$acc parallel loop collapse(4) copyin(fu,fw,fv) copy(dwdt,dudt,dvdt) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            jb=j-1
            kc=k+1
            dudt(icrm,i,j,k,na)=dudt(icrm,i,j,k,na)-(fu(i,j,k,icrm)-fu(i,jb,k,icrm))
            dvdt(icrm,i,j,k,na)=dvdt(icrm,i,j,k,na)-(fv(i,j,k,icrm)-fv(i,jb,k,icrm))
            dwdt(icrm,i,j,kc,na)=dwdt(icrm,i,j,kc,na)-(fw(i,j,k,icrm)-fw(i,jb,k,icrm))
          enddo
        enddo
      enddo
    enddo

    !$acc parallel loop collapse(2) copy(uwsb,vwsb) async(asyncid)
    do icrm = 1 , ncrms
      do k = 1 , nzm
        uwsb(k,icrm)=0.
        vwsb(k,icrm)=0.
      enddo
    enddo

    !-------------------------
    !$acc parallel loop collapse(4) copyin(v,dz,tk,rho,rhow,w,grdf_z,adz,adzw,u) copy(fu,uwsb,fv,vwsb,fw) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm-1
        do j=1,ny
          do i=1,nx
            jb=j-1
            kc=k+1
            ib=i-1
            rdz=1./dz(icrm)
            rdz2 = rdz*rdz * grdf_z(icrm,k)
            rdz25 = 0.25*rdz2
            iadz = 1./adz(icrm,k)
            iadzw= 1./adzw(icrm,kc)
            dzx=dz(icrm)/dx
            dzy=dz(icrm)/dy
            tkz=rdz2*tk(icrm,i,j,k)
            fw(i,j,kc,icrm)=-2.*tkz*(w(icrm,i,j,kc)-w(icrm,i,j,k))*rho(icrm,k)*iadz
            tkz=rdz25*(tk(icrm,i,j,k)+tk(icrm,ib,j,k)+tk(icrm,i,j,kc)+tk(icrm,ib,j,kc))
            fu(i,j,kc,icrm)=-tkz*( (u(icrm,i,j,kc)-u(icrm,i,j,k))*iadzw + (w(icrm,i,j,kc)-w(icrm,ib,j,kc))*dzx)*rhow(icrm,kc)
            tkz=rdz25*(tk(icrm,i,j,k)+tk(icrm,i,jb,k)+tk(icrm,i,j,kc)+tk(icrm,i,jb,kc))
            fv(i,j,kc,icrm)=-tkz*( (v(icrm,i,j,kc)-v(icrm,i,j,k))*iadzw + (w(icrm,i,j,kc)-w(icrm,i,jb,kc))*dzy)*rhow(icrm,kc)
            !$acc atomic update
            uwsb(kc,icrm)=uwsb(kc,icrm)+fu(i,j,kc,icrm)
            !$acc atomic update
            vwsb(kc,icrm)=vwsb(kc,icrm)+fv(i,j,kc,icrm)
          enddo
        enddo
      enddo
    enddo

    !$acc parallel loop collapse(3) copyin(tk,adz,rhow,grdf_z,w,rho,dz,fluxtv,fluxbu,fluxbv,fluxtu) copy(fw,vwsb,fv,fu,uwsb) async(asyncid)
    do icrm = 1 , ncrms
      do j=1,ny
        do i=1,nx
          rdz=1./dz(icrm)
          rdz2 = rdz*rdz * grdf_z(icrm,k)
          tkz=rdz2*grdf_z(icrm,nzm)*tk(icrm,i,j,nzm)
          fw(i,j,nz,icrm)=-2.*tkz*(w(icrm,i,j,nz)-w(icrm,i,j,nzm))/adz(icrm,nzm)*rho(icrm,nzm)
          fu(i,j,1,icrm)=fluxbu(icrm,i,j) * rdz * rhow(icrm,1)
          fv(i,j,1,icrm)=fluxbv(icrm,i,j) * rdz * rhow(icrm,1)
          fu(i,j,nz,icrm)=fluxtu(i,j,icrm) * rdz * rhow(icrm,nz)
          fv(i,j,nz,icrm)=fluxtv(i,j,icrm) * rdz * rhow(icrm,nz)
          !$acc atomic update
          uwsb(1,icrm) = uwsb(1,icrm) + fu(i,j,1,icrm)
          !$acc atomic update
          vwsb(1,icrm) = vwsb(1,icrm) + fv(i,j,1,icrm)
        enddo
      enddo
    enddo

    !$acc parallel loop collapse(4) copyin(rho,fv,adz,fu) copy(dvdt,dudt) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            kc=k+1
            rhoi = 1./(rho(icrm,k)*adz(icrm,k))
            dudt(icrm,i,j,k,na)=dudt(icrm,i,j,k,na)-(fu(i,j,kc,icrm)-fu(i,j,k,icrm))*rhoi
            dvdt(icrm,i,j,k,na)=dvdt(icrm,i,j,k,na)-(fv(i,j,kc,icrm)-fv(i,j,k,icrm))*rhoi
          enddo
        enddo
      enddo ! k
    enddo

    !$acc parallel loop collapse(4) copyin(adzw,rhow,fw) copy(dwdt) async(asyncid)
    do icrm = 1 , ncrms
      do k=2,nzm
        do j=1,ny
          do i=1,nx
            rhoi = 1./(rhow(icrm,k)*adzw(icrm,k))
            dwdt(icrm,i,j,k,na)=dwdt(icrm,i,j,k,na)-(fw(i,j,k+1,icrm)-fw(i,j,k,icrm))*rhoi
          enddo
        enddo
      enddo ! k
    enddo

    !$acc exit data delete(fu,fv,fw) async(asyncid)

  end subroutine diffuse_mom3D

end module diffuse_mom3D_mod
