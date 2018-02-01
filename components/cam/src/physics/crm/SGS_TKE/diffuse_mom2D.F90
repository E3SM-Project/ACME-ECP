module diffuse_mom2D_mod
  implicit none

contains

  subroutine diffuse_mom2D(grdf_x, grdf_z, dimx1_d, dimx2_d, dimy1_d, dimy2_d, tk, ncrms, icrm)

    !        momentum tendency due to SGS diffusion

    use vars
    use params, only: docolumn, crm_rknd
    implicit none
    integer, intent(in) :: ncrms, icrm
    integer :: dimx1_d, dimx2_d, dimy1_d, dimy2_d
    real(crm_rknd) tk  (ncrms,dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm) ! SGS eddy viscosity
    real(crm_rknd) grdf_x(ncrms,nzm)! grid factor for eddy diffusion in x
    real(crm_rknd) grdf_z(ncrms,nzm)! grid factor for eddy diffusion in z

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
        dxz=dx/(dz(icrm)*adzw(icrm,kc))
        rdx21=rdx2 * grdf_x(icrm,k)
        rdx251=rdx25 * grdf_x(icrm,k)

        do i=0,nx
          ic=i+1
          tkx=rdx21*tk(icrm,i,j,k)
          fu(i,j,k)=-2.*tkx*(u(icrm,ic,j,k)-u(icrm,i,j,k))
          fv(i,j,k)=-tkx*(v(icrm,ic,j,k)-v(icrm,i,j,k))
          tkx=rdx251*(tk(icrm,i,j,k)+tk(icrm,ic,j,k)+tk(icrm,i,j,kcu)+tk(icrm,ic,j,kcu))
          fw(i,j,k)=-tkx*(w(icrm,ic,j,kc)-w(icrm,i,j,kc)+(u(icrm,ic,j,kcu)-u(icrm,ic,j,k))*dxz)
        end do
        do i=1,nx
          ib=i-1
          dudt(icrm,i,j,k,na)=dudt(icrm,i,j,k,na)-(fu(i,j,k)-fu(ib,j,k))
          dvdt(icrm,i,j,k,na)=dvdt(icrm,i,j,k,na)-(fv(i,j,k)-fv(ib,j,k))
          dwdt(icrm,i,j,kc,na)=dwdt(icrm,i,j,kc,na)-(fw(i,j,k)-fw(ib,j,k))
        end do

      end do

    end if

    !-------------------------
    rdz=1./dz(icrm)
    dzx=dz(icrm)/dx

    do k=1,nzm-1
      kc=k+1
      uwsb(icrm,kc)=0.
      vwsb(icrm,kc)=0.
      iadz = 1./adz(icrm,k)
      iadzw= 1./adzw(icrm,kc)
      rdz2=rdz*rdz *grdf_z(icrm,k)
      rdz25=0.25*rdz2
      do i=1,nx
        ib=i-1
        tkz=rdz2*tk(icrm,i,j,k)
        fw(i,j,kc)=-2.*tkz*(w(icrm,i,j,kc)-w(icrm,i,j,k))*rho(icrm,k)*iadz
        tkz=rdz25*(tk(icrm,i,j,k)+tk(icrm,ib,j,k)+tk(icrm,i,j,kc)+tk(icrm,ib,j,kc))
        fu(i,j,kc)=-tkz*( (u(icrm,i,j,kc)-u(icrm,i,j,k))*iadzw + &
        (w(icrm,i,j,kc)-w(icrm,ib,j,kc))*dzx)*rhow(icrm,kc)
        fv(i,j,kc)=-tkz*(v(icrm,i,j,kc)-v(icrm,i,j,k))*iadzw*rhow(icrm,kc)
        uwsb(icrm,kc)=uwsb(icrm,kc)+fu(i,j,kc)
        vwsb(icrm,kc)=vwsb(icrm,kc)+fv(i,j,kc)
      end do
    end do

    uwsb(icrm,1) = 0.
    vwsb(icrm,1) = 0.

    do i=1,nx
      tkz=rdz2*grdf_z(icrm,nzm)*tk(icrm,i,j,nzm)
      fw(i,j,nz)=-2.*tkz*(w(icrm,i,j,nz)-w(icrm,i,j,nzm))/adz(icrm,nzm)*rho(icrm,nzm)
      fu(i,j,1)=fluxbu(icrm,i,j) * rdz * rhow(icrm,1)
      fv(i,j,1)=fluxbv(icrm,i,j) * rdz * rhow(icrm,1)
      fu(i,j,nz)=fluxtu(icrm,i,j) * rdz * rhow(icrm,nz)
      fv(i,j,nz)=fluxtv(icrm,i,j) * rdz * rhow(icrm,nz)
      uwsb(icrm,1) = uwsb(icrm,1) + fu(i,j,1)
      vwsb(icrm,1) = vwsb(icrm,1) + fv(i,j,1)
    end do


    do k=1,nzm
      kc=k+1
      rhoi = 1./(rho(icrm,k)*adz(icrm,k))
      do i=1,nx
        dudt(icrm,i,j,k,na)=dudt(icrm,i,j,k,na)-(fu(i,j,kc)-fu(i,j,k))*rhoi
        dvdt(icrm,i,j,k,na)=dvdt(icrm,i,j,k,na)-(fv(i,j,kc)-fv(i,j,k))*rhoi
      end do
    end do ! k

    do k=2,nzm
      rhoi = 1./(rhow(icrm,k)*adzw(icrm,k))
      do i=1,nx
        dwdt(icrm,i,j,k,na)=dwdt(icrm,i,j,k,na)-(fw(i,j,k+1)-fw(i,j,k))*rhoi
      end do
    end do ! k


  end subroutine diffuse_mom2D


end module diffuse_mom2D_mod
