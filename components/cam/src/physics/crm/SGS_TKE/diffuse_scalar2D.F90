module diffuse_scalar2D_mod
  implicit none

contains
  subroutine diffuse_scalar2D (ncrms,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_z,field,fluxb,fluxt,tkh,rho,rhow,flux)

    use grid
    use params
    implicit none
    integer, intent(in) :: ncrms
    ! input
    integer :: dimx1_d,dimx2_d,dimy1_d,dimy2_d
    real(crm_rknd) grdf_x(nzm,ncrms)! grid factor for eddy diffusion in x
    real(crm_rknd) grdf_z(nzm,ncrms)! grid factor for eddy diffusion in z
    real(crm_rknd) field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm,ncrms) ! scalar
    real(crm_rknd) tkh(0:nxp1, 1-YES3D:nyp1, nzm,ncrms) ! eddy conductivity
    real(crm_rknd) fluxb(nx,ny,ncrms)   ! bottom flux
    real(crm_rknd) fluxt(nx,ny,ncrms)   ! top flux
    real(crm_rknd) rho(nzm,ncrms)
    real(crm_rknd) rhow(nz,ncrms)
    real(crm_rknd) flux(nz,ncrms)
    ! local
    real(crm_rknd) flx(0:nx,1,0:nzm,ncrms)
    real(crm_rknd) dfdt(nx,ny,nzm,ncrms)
    real(crm_rknd) rdx2,rdz2,rdz,rdx5,rdz5,tmp
    real(crm_rknd) dxz,dzx,tkx,tkz,rhoi
    integer i,j,k,ib,ic,kc,kb,icrm

    do icrm = 1 , ncrms

    if(.not.dosgs.and..not.docolumn) return

    rdx2=1./(dx*dx)
    rdz2=1./(dz(icrm)*dz(icrm))
    rdz=1./dz(icrm)
    dxz=dx/dz(icrm)
    dzx=dz(icrm)/dx

    j=1

    dfdt(:,:,:,icrm)=0.

    if(dowallx) then

      if(mod(rank,nsubdomains_x).eq.0) then
        do k=1,nzm
          field(0,j,k,icrm) = field(1,j,k,icrm)
        end do
      end if
      if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
        do k=1,nzm
          field(nx+1,j,k,icrm) = field(nx,j,k,icrm)
        end do
      end if

    end if


    if(.not.docolumn) then


      do k=1,nzm

        rdx5=0.5*rdx2  *grdf_x(k,icrm)

        do i=0,nx
          ic=i+1
          tkx=rdx5*(tkh(i,j,k,icrm)+tkh(ic,j,k,icrm))
          flx(i,j,k,icrm)=-tkx*(field(ic,j,k,icrm)-field(i,j,k,icrm))
        end do
        do i=1,nx
          ib=i-1
          dfdt(i,j,k,icrm)=dfdt(i,j,k,icrm)-(flx(i,j,k,icrm)-flx(ib,j,k,icrm))
        end do

      end do

    end if

    flux(1,icrm) = 0.
    tmp=1./adzw(nz,icrm)
    do i=1,nx
      flx(i,j,0,icrm)=fluxb(i,j,icrm)*rdz*rhow(1,icrm)
      flx(i,j,nzm,icrm)=fluxt(i,j,icrm)*rdz*tmp*rhow(nz,icrm)
      flux(1,icrm) = flux(1,icrm) + flx(i,j,0,icrm)
    end do


    do k=1,nzm-1
      kc=k+1
      flux(kc,icrm)=0.
      rhoi = rhow(kc,icrm)/adzw(kc,icrm)
      rdz5=0.5*rdz2 * grdf_z(k,icrm)
      do i=1,nx
        tkz=rdz5*(tkh(i,j,k,icrm)+tkh(i,j,kc,icrm))
        flx(i,j,k,icrm)=-tkz*(field(i,j,kc,icrm)-field(i,j,k,icrm))*rhoi
        flux(kc,icrm) = flux(kc,icrm) + flx(i,j,k,icrm)
      end do
    end do

    do k=1,nzm
      kb=k-1
      rhoi = 1./(adz(k,icrm)*rho(k,icrm))
      do i=1,nx
        dfdt(i,j,k,icrm)=dtn*(dfdt(i,j,k,icrm)-(flx(i,j,k,icrm)-flx(i,j,kb,icrm))*rhoi)
        field(i,j,k,icrm)=field(i,j,k,icrm) + dfdt(i,j,k,icrm)
      end do
    end do

    enddo

  end subroutine diffuse_scalar2D

end module diffuse_scalar2D_mod
