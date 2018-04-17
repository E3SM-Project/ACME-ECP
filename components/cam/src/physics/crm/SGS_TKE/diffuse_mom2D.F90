module diffuse_mom2D_mod
  implicit none

contains

  subroutine diffuse_mom2D(grdf_x, grdf_z, dimx1_d, dimx2_d, dimy1_d, dimy2_d, tk, ncrms)

    !        momentum tendency due to SGS diffusion

    use vars
    use params, only: docolumn, crm_rknd
    implicit none
    integer, intent(in) :: ncrms
    integer :: dimx1_d, dimx2_d, dimy1_d, dimy2_d
    real(crm_rknd), pointer :: tk  (:,:,:,:) ! SGS eddy viscosity
    real(crm_rknd) grdf_x(ncrms,nzm)! grid factor for eddy diffusion in x
    real(crm_rknd) grdf_z(ncrms,nzm)! grid factor for eddy diffusion in z

    real(crm_rknd) rdx2,rdz,rdx25, tmp1, tmp2

    integer i,j,k,ic,ib,kc,kcu, icrm
    real(crm_rknd) tkz, rhoi
    real(crm_rknd) fu(ncrms,0:nx,1,nz),fv(ncrms,0:nx,1,nz),fw(ncrms,0:nx,1,nz)

    !$acc enter data create(fu,fv,fw) async(1)

    rdx2=1./dx/dx
    rdx25=0.25*rdx2

    j=1

    if(.not.docolumn) then
      !$acc parallel loop gang vector collapse(3) default(present) async(1)
      do k=1,nzm
        do i=1,nx
          do icrm = 1 , ncrms
            kc=k+1
            kcu=min(kc,nzm)
            ib=i-1
            ic=i+1

            tmp1 = -2.*rdx2*grdf_x(icrm,k)*tk(icrm,i ,j,k)*(u(icrm,ic,j,k)-u(icrm,i ,j,k))
            tmp2 = -2.*rdx2*grdf_x(icrm,k)*tk(icrm,ib,j,k)*(u(icrm,i ,j,k)-u(icrm,ib,j,k))
            dudt(icrm,i,j,k,na)=dudt(icrm,i,j,k,na)-(tmp1-tmp2)

            tmp1 = -rdx2*grdf_x(icrm,k)*tk(icrm,i ,j,k)*(v(icrm,ic,j,k)-v(icrm,i ,j,k))
            tmp2 = -rdx2*grdf_x(icrm,k)*tk(icrm,ib,j,k)*(v(icrm,i ,j,k)-v(icrm,ib,j,k))
            dvdt(icrm,i,j,k,na)=dvdt(icrm,i,j,k,na)-(tmp1-tmp2)

            tmp1 = -rdx25*grdf_x(icrm,k)*(tk(icrm,i ,j,k)+tk(icrm,ic,j,k)+tk(icrm,i ,j,kcu)+tk(icrm,ic,j,kcu))*(w(icrm,ic,j,kc)-w(icrm,i ,j,kc)+(u(icrm,ic,j,kcu)-u(icrm,ic,j,k))*dx/(dz(icrm)*adzw(icrm,kc)))
            tmp2 = -rdx25*grdf_x(icrm,k)*(tk(icrm,ib,j,k)+tk(icrm,i ,j,k)+tk(icrm,ib,j,kcu)+tk(icrm,i ,j,kcu))*(w(icrm,i ,j,kc)-w(icrm,ib,j,kc)+(u(icrm,i ,j,kcu)-u(icrm,i ,j,k))*dx/(dz(icrm)*adzw(icrm,kc)))
            dwdt(icrm,i,j,kc,na)=dwdt(icrm,i,j,kc,na)-(tmp1-tmp2)
          enddo
        end do
      end do
    end if

    !-------------------------

    !$acc parallel loop gang vector collapse(2) default(present) async(1)
    do k = 1 , nz
      do icrm = 1 , ncrms
        uwsb(icrm,k) = 0.
        vwsb(icrm,k) = 0.
      enddo
    enddo

    !$acc parallel loop gang vector collapse(3) default(present) async(1)
    do k=1,nzm-1
      do i=1,nx
        do icrm = 1 , ncrms
          kc=k+1
          ib=i-1
          rdz = 1 / dz(icrm)
          tkz=rdz*rdz *grdf_z(icrm,k)*tk(icrm,i,j,k)
          fw(icrm,i,j,kc)=-2.*tkz*(w(icrm,i,j,kc)-w(icrm,i,j,k))*rho(icrm,k)/adz(icrm,k)
          tkz=0.25*rdz*rdz *grdf_z(icrm,k)*(tk(icrm,i,j,k)+tk(icrm,ib,j,k)+tk(icrm,i,j,kc)+tk(icrm,ib,j,kc))
          fu(icrm,i,j,kc)=-tkz*( (u(icrm,i,j,kc)-u(icrm,i,j,k))/adzw(icrm,kc) + (w(icrm,i,j,kc)-w(icrm,ib,j,kc))*dz(icrm) / dx)*rhow(icrm,kc)
          fv(icrm,i,j,kc)=-tkz*(v(icrm,i,j,kc)-v(icrm,i,j,k))/adzw(icrm,kc)*rhow(icrm,kc)
          !$acc atomic update
          uwsb(icrm,kc)=uwsb(icrm,kc)+fu(icrm,i,j,kc)
          !$acc atomic update
          vwsb(icrm,kc)=vwsb(icrm,kc)+fv(icrm,i,j,kc)
        enddo
      end do
    end do


    !$acc parallel loop gang vector collapse(2) default(present) async(1)
    do i=1,nx
      do icrm = 1 , ncrms
        rdz = 1 / dz(icrm)
        tkz=rdz*rdz *grdf_z(icrm,nzm-1)*grdf_z(icrm,nzm)*tk(icrm,i,j,nzm)
        fw(icrm,i,j,nz)=-2.*tkz*(w(icrm,i,j,nz)-w(icrm,i,j,nzm))/adz(icrm,nzm)*rho(icrm,nzm)
        fu(icrm,i,j,1)=fluxbu(icrm,i,j) * rdz * rhow(icrm,1)
        fv(icrm,i,j,1)=fluxbv(icrm,i,j) * rdz * rhow(icrm,1)
        fu(icrm,i,j,nz)=fluxtu(icrm,i,j) * rdz * rhow(icrm,nz)
        fv(icrm,i,j,nz)=fluxtv(icrm,i,j) * rdz * rhow(icrm,nz)
          !$acc atomic update
        uwsb(icrm,1) = uwsb(icrm,1) + fu(icrm,i,j,1)
          !$acc atomic update
        vwsb(icrm,1) = vwsb(icrm,1) + fv(icrm,i,j,1)
      enddo
    end do


    !$acc parallel loop gang vector collapse(3) default(present) async(1)
    do k=1,nzm
      do i=1,nx
        do icrm = 1 , ncrms
          kc=k+1
          rhoi = 1./(rho(icrm,k)*adz(icrm,k))
          dudt(icrm,i,j,k,na)=dudt(icrm,i,j,k,na)-(fu(icrm,i,j,kc)-fu(icrm,i,j,k))*rhoi
          dvdt(icrm,i,j,k,na)=dvdt(icrm,i,j,k,na)-(fv(icrm,i,j,kc)-fv(icrm,i,j,k))*rhoi
          if (k >= 2) then
            rhoi = 1./(rhow(icrm,k)*adzw(icrm,k))
            dwdt(icrm,i,j,k,na)=dwdt(icrm,i,j,k,na)-(fw(icrm,i,j,k+1)-fw(icrm,i,j,k))*rhoi
          endif
        enddo
      end do
    end do ! k

    !$acc exit data delete(fu,fv,fw) async(1)

  end subroutine diffuse_mom2D


end module diffuse_mom2D_mod
