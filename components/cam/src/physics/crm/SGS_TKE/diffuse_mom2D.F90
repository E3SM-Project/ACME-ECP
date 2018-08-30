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
    real(crm_rknd) tk  (dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm) ! SGS eddy viscosity
    real(crm_rknd) grdf_x(nzm)! grid factor for eddy diffusion in x
    real(crm_rknd) grdf_z(nzm)! grid factor for eddy diffusion in z

    real(crm_rknd) rdx2,rdz2,rdz,rdx25,rdz25,rdx21,rdx251
    real(crm_rknd) dxz,dzx

    integer i,j,k,ic,ib,kc,kcu
    real(crm_rknd) tkx, tkz, rhoi, iadzw, iadz
    real(crm_rknd) fu(0:nx,1,nz),fv(0:nx,1,nz),fw(0:nx,1,nz)

    rdx2=1./dx/dx
    rdx25=0.25*rdx2

    dxz=dx/dz

    j=1

    if(.not.docolumn) then
      do k=1,nzm

        kc=k+1
        kcu=min(kc,nzm)
        dxz=dx/(dz*adzw(kc))
        rdx21=rdx2 * grdf_x(k)
        rdx251=rdx25 * grdf_x(k)

        do i=0,nx
          ic=i+1
          tkx=rdx21*tk(i,j,k)
          fu(i,j,k)=-2.*tkx*(u(ic,j,k)-u(i,j,k))
          fv(i,j,k)=-tkx*(v(ic,j,k)-v(i,j,k))
          tkx=rdx251*(tk(i,j,k)+tk(ic,j,k)+tk(i,j,kcu)+tk(ic,j,kcu))
          fw(i,j,k)=-tkx*(w(ic,j,kc)-w(i,j,kc)+(u(ic,j,kcu)-u(ic,j,k))*dxz)
        end do
        do i=1,nx
          ib=i-1
          dudt(i,j,k,na,icrm)=dudt(i,j,k,na,icrm)-(fu(i,j,k)-fu(ib,j,k))
          dvdt(i,j,k,na,icrm)=dvdt(i,j,k,na,icrm)-(fv(i,j,k)-fv(ib,j,k))
          dwdt(i,j,kc,na,icrm)=dwdt(i,j,kc,na,icrm)-(fw(i,j,k)-fw(ib,j,k))
        end do

      end do

    end if

    !-------------------------
    rdz=1./dz
    dzx=dz/dx

    do k=1,nzm-1
      kc=k+1
      uwsb(kc,icrm)=0.
      vwsb(kc,icrm)=0.
      iadz = 1./adz(k)
      iadzw= 1./adzw(kc)
      rdz2=rdz*rdz *grdf_z(k)
      rdz25=0.25*rdz2
      do i=1,nx
        ib=i-1
        tkz=rdz2*tk(i,j,k)
        fw(i,j,kc)=-2.*tkz*(w(i,j,kc)-w(i,j,k))*rho(k,icrm)*iadz
        tkz=rdz25*(tk(i,j,k)+tk(ib,j,k)+tk(i,j,kc)+tk(ib,j,kc))
        fu(i,j,kc)=-tkz*( (u(i,j,kc)-u(i,j,k))*iadzw + &
        (w(i,j,kc)-w(ib,j,kc))*dzx)*rhow(kc,icrm)
        fv(i,j,kc)=-tkz*(v(i,j,kc)-v(i,j,k))*iadzw*rhow(kc,icrm)
        uwsb(kc,icrm)=uwsb(kc,icrm)+fu(i,j,kc)
        vwsb(kc,icrm)=vwsb(kc,icrm)+fv(i,j,kc)
      end do
    end do

    uwsb(1,icrm) = 0.
    vwsb(1,icrm) = 0.

    do i=1,nx
      tkz=rdz2*grdf_z(nzm)*tk(i,j,nzm)
      fw(i,j,nz)=-2.*tkz*(w(i,j,nz)-w(i,j,nzm))/adz(nzm)*rho(nzm,icrm)
      fu(i,j,1)=fluxbu(i,j,icrm) * rdz * rhow(1,icrm)
      fv(i,j,1)=fluxbv(i,j,icrm) * rdz * rhow(1,icrm)
      fu(i,j,nz)=fluxtu(i,j,icrm) * rdz * rhow(nz,icrm)
      fv(i,j,nz)=fluxtv(i,j,icrm) * rdz * rhow(nz,icrm)
      uwsb(1,icrm) = uwsb(1,icrm) + fu(i,j,1)
      vwsb(1,icrm) = vwsb(1,icrm) + fv(i,j,1)
    end do


    do k=1,nzm
      kc=k+1
      rhoi = 1./(rho(k,icrm)*adz(k))
      do i=1,nx
        dudt(i,j,k,na,icrm)=dudt(i,j,k,na,icrm)-(fu(i,j,kc)-fu(i,j,k))*rhoi
        dvdt(i,j,k,na,icrm)=dvdt(i,j,k,na,icrm)-(fv(i,j,kc)-fv(i,j,k))*rhoi
      end do
    end do ! k

    do k=2,nzm
      rhoi = 1./(rhow(k,icrm)*adzw(k))
      do i=1,nx
        dwdt(i,j,k,na,icrm)=dwdt(i,j,k,na,icrm)-(fw(i,j,k+1)-fw(i,j,k))*rhoi
      end do
    end do ! k


  end subroutine diffuse_mom2D


end module diffuse_mom2D_mod
