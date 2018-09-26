module diffuse_mom2D_mod
  implicit none

contains

  subroutine diffuse_mom2D(ncrms,icrm,grdf_x, grdf_z, dimx1_d, dimx2_d, dimy1_d, dimy2_d, tk)

    !        momentum tendency due to SGS diffusion

    use vars
    use params, only: docolumn, crm_rknd
    implicit none
    integer, intent(in) :: ncrms,icrm
    integer :: dimx1_d, dimx2_d, dimy1_d, dimy2_d
    real(crm_rknd) tk  (dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm,ncrms) ! SGS eddy viscosity
    real(crm_rknd) grdf_x(nzm,ncrms)! grid factor for eddy diffusion in x
    real(crm_rknd) grdf_z(nzm,ncrms)! grid factor for eddy diffusion in z

    real(crm_rknd) rdx2,rdz2,rdz,rdx25,rdz25,rdx21,rdx251
    real(crm_rknd) dxz,dzx

    integer i,j,k,ic,ib,kc,kcu
    real(crm_rknd) tkx, tkz, rhoi, iadzw, iadz
    real(crm_rknd) fu(0:nx,1,nz),fv(0:nx,1,nz),fw(0:nx,1,nz)

    rdx2=1./dx/dx
    rdx25=0.25*rdx2

    dxz=dx/dz(icrm)

    j=1

    if(.not.docolumn) then
      do k=1,nzm

        kc=k+1
        kcu=min(kc,nzm)
        dxz=dx/(dz(icrm)*adzw(kc,icrm))
        rdx21=rdx2 * grdf_x(k,icrm)
        rdx251=rdx25 * grdf_x(k,icrm)

        do i=0,nx
          ic=i+1
          tkx=rdx21*tk(i,j,k,icrm)
          fu(i,j,k)=-2.*tkx*(u(ic,j,k,icrm)-u(i,j,k,icrm))
          fv(i,j,k)=-tkx*(v(ic,j,k,icrm)-v(i,j,k,icrm))
          tkx=rdx251*(tk(i,j,k,icrm)+tk(ic,j,k,icrm)+tk(i,j,kcu,icrm)+tk(ic,j,kcu,icrm))
          fw(i,j,k)=-tkx*(w(ic,j,kc,icrm)-w(i,j,kc,icrm)+(u(ic,j,kcu,icrm)-u(ic,j,k,icrm))*dxz)
        end do
        do i=1,nx
          ib=i-1
          dudt(i,j,k,na(icrm),icrm)=dudt(i,j,k,na(icrm),icrm)-(fu(i,j,k)-fu(ib,j,k))
          dvdt(i,j,k,na(icrm),icrm)=dvdt(i,j,k,na(icrm),icrm)-(fv(i,j,k)-fv(ib,j,k))
          dwdt(i,j,kc,na(icrm),icrm)=dwdt(i,j,kc,na(icrm),icrm)-(fw(i,j,k)-fw(ib,j,k))
        end do

      end do

    end if

    !-------------------------
    rdz=1./dz(icrm)
    dzx=dz(icrm)/dx

    do k=1,nzm-1
      kc=k+1
      uwsb(kc,icrm)=0.
      vwsb(kc,icrm)=0.
      iadz = 1./adz(k,icrm)
      iadzw= 1./adzw(kc,icrm)
      rdz2=rdz*rdz *grdf_z(k,icrm)
      rdz25=0.25*rdz2
      do i=1,nx
        ib=i-1
        tkz=rdz2*tk(i,j,k,icrm)
        fw(i,j,kc)=-2.*tkz*(w(i,j,kc,icrm)-w(i,j,k,icrm))*rho(k,icrm)*iadz
        tkz=rdz25*(tk(i,j,k,icrm)+tk(ib,j,k,icrm)+tk(i,j,kc,icrm)+tk(ib,j,kc,icrm))
        fu(i,j,kc)=-tkz*( (u(i,j,kc,icrm)-u(i,j,k,icrm))*iadzw + &
        (w(i,j,kc,icrm)-w(ib,j,kc,icrm))*dzx)*rhow(kc,icrm)
        fv(i,j,kc)=-tkz*(v(i,j,kc,icrm)-v(i,j,k,icrm))*iadzw*rhow(kc,icrm)
        uwsb(kc,icrm)=uwsb(kc,icrm)+fu(i,j,kc)
        vwsb(kc,icrm)=vwsb(kc,icrm)+fv(i,j,kc)
      end do
    end do

    uwsb(1,icrm) = 0.
    vwsb(1,icrm) = 0.

    do i=1,nx
      tkz=rdz2*grdf_z(nzm,icrm)*tk(i,j,nzm,icrm)
      fw(i,j,nz)=-2.*tkz*(w(i,j,nz,icrm)-w(i,j,nzm,icrm))/adz(nzm,icrm)*rho(nzm,icrm)
      fu(i,j,1)=fluxbu(i,j,icrm) * rdz * rhow(1,icrm)
      fv(i,j,1)=fluxbv(i,j,icrm) * rdz * rhow(1,icrm)
      fu(i,j,nz)=fluxtu(i,j,icrm) * rdz * rhow(nz,icrm)
      fv(i,j,nz)=fluxtv(i,j,icrm) * rdz * rhow(nz,icrm)
      uwsb(1,icrm) = uwsb(1,icrm) + fu(i,j,1)
      vwsb(1,icrm) = vwsb(1,icrm) + fv(i,j,1)
    end do


    do k=1,nzm
      kc=k+1
      rhoi = 1./(rho(k,icrm)*adz(k,icrm))
      do i=1,nx
        dudt(i,j,k,na(icrm),icrm)=dudt(i,j,k,na(icrm),icrm)-(fu(i,j,kc)-fu(i,j,k))*rhoi
        dvdt(i,j,k,na(icrm),icrm)=dvdt(i,j,k,na(icrm),icrm)-(fv(i,j,kc)-fv(i,j,k))*rhoi
      end do
    end do ! k

    do k=2,nzm
      rhoi = 1./(rhow(k,icrm)*adzw(k,icrm))
      do i=1,nx
        dwdt(i,j,k,na(icrm),icrm)=dwdt(i,j,k,na(icrm),icrm)-(fw(i,j,k+1)-fw(i,j,k))*rhoi
      end do
    end do ! k


  end subroutine diffuse_mom2D


end module diffuse_mom2D_mod
