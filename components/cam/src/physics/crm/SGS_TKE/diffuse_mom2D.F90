module diffuse_mom2D_mod
  implicit none

contains

  subroutine diffuse_mom2D(ncrms,grdf_x, grdf_z, dimx1_d, dimx2_d, dimy1_d, dimy2_d, tk)

    !        momentum tendency due to SGS diffusion

    use vars
    use params, only: docolumn, crm_rknd
    implicit none
    integer, intent(in) :: ncrms
    integer :: dimx1_d, dimx2_d, dimy1_d, dimy2_d
    real(crm_rknd) tk  (dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm,ncrms) ! SGS eddy viscosity
    real(crm_rknd) grdf_x(nzm,ncrms)! grid factor for eddy diffusion in x
    real(crm_rknd) grdf_z(nzm,ncrms)! grid factor for eddy diffusion in z

    real(crm_rknd) rdx2,rdz2,rdz,rdx25,rdz25,rdx21,rdx251
    real(crm_rknd) dxz,dzx

    integer i,j,k,ic,ib,kc,kcu,icrm
    real(crm_rknd) tkx, tkz, rhoi, iadzw, iadz
    real(crm_rknd) fu(0:nx,1,nz,ncrms)
    real(crm_rknd) fv(0:nx,1,nz,ncrms)
    real(crm_rknd) fw(0:nx,1,nz,ncrms)

    !$acc enter data create(fu,fv,fw) async(1)

    rdx2=1./dx/dx
    rdx25=0.25*rdx2

    j=1

    if( .not. docolumn ) then
      !$acc parallel loop collapse(3) async(1)
      do icrm = 1 , ncrms
        do k=1,nzm
          do i=0,nx
            kc=k+1
            kcu=min(kc,nzm)
            dxz=dx/(dz(icrm)*adzw(kc,icrm))
            rdx21=rdx2 * grdf_x(k,icrm)
            rdx251=rdx25 * grdf_x(k,icrm)
            ic=i+1
            tkx=rdx21*tk(i,j,k,icrm)
            fu(i,j,k,icrm)=-2.*tkx*(u(ic,j,k,icrm)-u(i,j,k,icrm))
            fv(i,j,k,icrm)=-tkx*(v(ic,j,k,icrm)-v(i,j,k,icrm))
            tkx=rdx251*(tk(i,j,k,icrm)+tk(ic,j,k,icrm)+tk(i,j,kcu,icrm)+tk(ic,j,kcu,icrm))
            fw(i,j,k,icrm)=-tkx*(w(ic,j,kc,icrm)-w(i,j,kc,icrm)+(u(ic,j,kcu,icrm)-u(ic,j,k,icrm))*dxz)
          end do
        end do
      end do
      !$acc parallel loop collapse(3) async(1)
      do icrm = 1 , ncrms
        do k=1,nzm
          do i=1,nx
            kc=k+1
            ib=i-1
            dudt(i,j,k,na,icrm)=dudt(i,j,k,na,icrm)-(fu(i,j,k,icrm)-fu(ib,j,k,icrm))
            dvdt(i,j,k,na,icrm)=dvdt(i,j,k,na,icrm)-(fv(i,j,k,icrm)-fv(ib,j,k,icrm))
            dwdt(i,j,kc,na,icrm)=dwdt(i,j,kc,na,icrm)-(fw(i,j,k,icrm)-fw(ib,j,k,icrm))
          end do
        end do
      end do
    end if

    !-------------------------

    !$acc parallel loop collapse(2) async(1)
    do icrm = 1 , ncrms
      do k = 1 , nzm
        uwsb(k,icrm)=0.
        vwsb(k,icrm)=0.
      enddo
    enddo

    !$acc parallel loop collapse(3) async(1)
    do icrm = 1 , ncrms
      do k=1,nzm-1
        do i=1,nx
          kc=k+1
          rdz=1./dz(icrm)
          rdz2=rdz*rdz *grdf_z(k,icrm)
          rdz25=0.25*rdz2
          dzx=dz(icrm)/dx
          iadz = 1./adz(k,icrm)
          iadzw= 1./adzw(kc,icrm)
          ib=i-1
          tkz=rdz2*tk(i,j,k,icrm)
          fw(i,j,kc,icrm)=-2.*tkz*(w(i,j,kc,icrm)-w(i,j,k,icrm))*rho(k,icrm)*iadz
          tkz=rdz25*(tk(i,j,k,icrm)+tk(ib,j,k,icrm)+tk(i,j,kc,icrm)+tk(ib,j,kc,icrm))
          fu(i,j,kc,icrm)=-tkz*( (u(i,j,kc,icrm)-u(i,j,k,icrm))*iadzw + (w(i,j,kc,icrm)-w(ib,j,kc,icrm))*dzx)*rhow(kc,icrm)
          fv(i,j,kc,icrm)=-tkz*(v(i,j,kc,icrm)-v(i,j,k,icrm))*iadzw*rhow(kc,icrm)
          !$acc atomic update
          uwsb(kc,icrm)=uwsb(kc,icrm)+fu(i,j,kc,icrm)
          !$acc atomic update
          vwsb(kc,icrm)=vwsb(kc,icrm)+fv(i,j,kc,icrm)
        end do
      end do
    end do

    !$acc parallel loop collapse(2) async(1)
    do icrm = 1 , ncrms
      do i=1,nx
        rdz=1./dz(icrm)
        rdz2=rdz*rdz *grdf_z(k,icrm)
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
      end do
    end do

    !$acc parallel loop collapse(3) async(1)
    do icrm = 1 , ncrms
      do k=1,nzm
        do i=1,nx
          kc=k+1
          rhoi = 1./(rho(k,icrm)*adz(k,icrm))
          dudt(i,j,k,na,icrm)=dudt(i,j,k,na,icrm)-(fu(i,j,kc,icrm)-fu(i,j,k,icrm))*rhoi
          dvdt(i,j,k,na,icrm)=dvdt(i,j,k,na,icrm)-(fv(i,j,kc,icrm)-fv(i,j,k,icrm))*rhoi
        end do
      end do ! k
    end do ! k

    !$acc parallel loop collapse(3) async(1)
    do icrm = 1 , ncrms
      do k=2,nzm
        do i=1,nx
          rhoi = 1./(rhow(k,icrm)*adzw(k,icrm))
          dwdt(i,j,k,na,icrm)=dwdt(i,j,k,na,icrm)-(fw(i,j,k+1,icrm)-fw(i,j,k,icrm))*rhoi
        end do
      end do ! k
    end do ! k

    !$acc exit data delete(fu,fv,fw) async(1)

  end subroutine diffuse_mom2D


end module diffuse_mom2D_mod
