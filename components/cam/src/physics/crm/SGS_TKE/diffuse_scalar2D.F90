module diffuse_scalar2D_mod
  implicit none

contains
  subroutine diffuse_scalar2D (ncrms,icrm,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_z,field,fluxb,fluxt,tkh,rho,rhow,flux)

    use grid
    use params
    implicit none
    integer, intent(in) :: ncrms,icrm
    ! input
    integer :: dimx1_d,dimx2_d,dimy1_d,dimy2_d
    real(crm_rknd) grdf_x(nzm,ncrms)! grid factor for eddy diffusion in x
    real(crm_rknd) grdf_z(nzm,ncrms)! grid factor for eddy diffusion in z
    real(crm_rknd) field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! scalar
    real(crm_rknd) tkh(0:nxp1, 1-YES3D:nyp1, nzm,ncrms) ! eddy conductivity
    real(crm_rknd) fluxb(nx,ny)   ! bottom flux
    real(crm_rknd) fluxt(nx,ny)   ! top flux
    real(crm_rknd) rho(nzm,ncrms)
    real(crm_rknd) rhow(nz,ncrms)
    real(crm_rknd) flux(nz)
    ! local
    real(crm_rknd) flx(0:nx,1,0:nzm)
    real(crm_rknd) dfdt(nx,ny,nzm)
    real(crm_rknd) rdx2,rdz2,rdz,rdx5,rdz5,tmp
    real(crm_rknd) dxz,dzx,tkx,tkz,rhoi
    integer i,j,k,ib,ic,kc,kb

    if(.not.dosgs.and..not.docolumn) return

    rdx2=1./(dx*dx)
    rdz2=1./(dz(icrm)*dz(icrm))
    rdz=1./dz(icrm)
    dxz=dx/dz(icrm)
    dzx=dz(icrm)/dx

    j=1

    dfdt(:,:,:)=0.

    if(dowallx) then

      if(mod(rank,nsubdomains_x).eq.0) then
        do k=1,nzm
          field(0,j,k) = field(1,j,k)
        end do
      end if
      if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
        do k=1,nzm
          field(nx+1,j,k) = field(nx,j,k)
        end do
      end if

    end if


    if(.not.docolumn) then


      do k=1,nzm

        rdx5=0.5*rdx2  *grdf_x(k,icrm)

        do i=0,nx
          ic=i+1
          tkx=rdx5*(tkh(i,j,k,icrm)+tkh(ic,j,k,icrm))
          flx(i,j,k)=-tkx*(field(ic,j,k)-field(i,j,k))
        end do
        do i=1,nx
          ib=i-1
          dfdt(i,j,k)=dfdt(i,j,k)-(flx(i,j,k)-flx(ib,j,k))
        end do

      end do

    end if

    flux(1) = 0.
    tmp=1./adzw(nz,icrm)
    do i=1,nx
      flx(i,j,0)=fluxb(i,j)*rdz*rhow(1,icrm)
      flx(i,j,nzm)=fluxt(i,j)*rdz*tmp*rhow(nz,icrm)
      flux(1) = flux(1) + flx(i,j,0)
    end do


    do k=1,nzm-1
      kc=k+1
      flux(kc)=0.
      rhoi = rhow(kc,icrm)/adzw(kc,icrm)
      rdz5=0.5*rdz2 * grdf_z(k,icrm)
      do i=1,nx
        tkz=rdz5*(tkh(i,j,k,icrm)+tkh(i,j,kc,icrm))
        flx(i,j,k)=-tkz*(field(i,j,kc)-field(i,j,k))*rhoi
        flux(kc) = flux(kc) + flx(i,j,k)
      end do
    end do

    do k=1,nzm
      kb=k-1
      rhoi = 1./(adz(k,icrm)*rho(k,icrm))
      do i=1,nx
        dfdt(i,j,k)=dtn*(dfdt(i,j,k)-(flx(i,j,k)-flx(i,j,kb))*rhoi)
        field(i,j,k)=field(i,j,k) + dfdt(i,j,k)
      end do
    end do

  end subroutine diffuse_scalar2D

end module diffuse_scalar2D_mod
