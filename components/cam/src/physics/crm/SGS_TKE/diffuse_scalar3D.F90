module diffuse_scalar3D_mod
  implicit none

contains

  subroutine diffuse_scalar3D (ncrms,icrm,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_y,grdf_z,field,fluxb,fluxt,tkh,rho,rhow,flux)

    use grid
    use params
    use task_util_mod, only: task_rank_to_index
    implicit none
    integer, intent(in) :: ncrms,icrm
    ! input
    integer :: dimx1_d,dimx2_d,dimy1_d,dimy2_d
    real(crm_rknd) grdf_x(nzm,ncrms)! grid factor for eddy diffusion in x
    real(crm_rknd) grdf_y(nzm,ncrms)! grid factor for eddy diffusion in y
    real(crm_rknd) grdf_z(nzm,ncrms)! grid factor for eddy diffusion in z
    real(crm_rknd) field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! scalar
    real(crm_rknd) tkh(0:nxp1,1-YES3D:nyp1,nzm,ncrms) ! eddy conductivity
    real(crm_rknd) fluxb(nx,ny)   ! bottom flux
    real(crm_rknd) fluxt(nx,ny)   ! top flux
    real(crm_rknd) rho(nzm,ncrms)
    real(crm_rknd) rhow(nz,ncrms)
    real(crm_rknd) flux(nz)
    ! local
    real(crm_rknd) flx(0:nx,0:ny,0:nzm)
    real(crm_rknd) dfdt(nx,ny,nz)
    real(crm_rknd) rdx2,rdy2,rdz2,rdz,rdx5,rdy5,rdz5,tmp
    real(crm_rknd) dxy,dxz,dyx,dyz,dzx,dzy,tkx,tky,tkz,rhoi
    integer i,j,k,ib,ic,jb,jc,kc,kb


    if(.not.dosgs) return

    rdx2=1./(dx*dx)
    rdy2=1./(dy*dy)
    rdz2=1./(dz(icrm)*dz(icrm))
    rdz=1./dz(icrm)
    dxy=dx/dy
    dxz=dx/dz(icrm)
    dyx=dy/dx
    dyz=dy/dz(icrm)
    dzx=dz(icrm)/dx
    dzy=dz(icrm)/dy

    dfdt(:,:,:)=0.

    !-----------------------------------------
    if(dowallx) then

      if(mod(rank,nsubdomains_x).eq.0) then
        do k=1,nzm
          do j=1,ny
            field(0,j,k) = field(1,j,k)
          end do
        end do
      end if
      if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
        do k=1,nzm
          do j=1,ny
            field(nx+1,j,k) = field(nx,j,k)
          end do
        end do
      end if

    end if

    if(dowally) then

      if(rank.lt.nsubdomains_x) then
        do k=1,nzm
          do i=1,nx
            field(i,1-YES3D,k) = field(i,1,k)
          end do
        end do
      end if
      if(rank.gt.nsubdomains-nsubdomains_x-1) then
        do k=1,nzm
          do i=1,ny
            field(i,ny+YES3D,k) = field(i,ny,k)
          end do
        end do
      end if

    end if



    if(dowally) then

      call task_rank_to_index(rank, ib, jb)
      if(jb.eq.0) then
        do k=1,nzm
          do i=1,nx
            field(i,1-YES3D,k) = field(i,1,k)
          end do
        end do
      end if
      if(jb.eq.nsubdomains_y-1) then
        do k=1,nzm
          do i=1,nx
            field(i,ny+YES3D,k) = field(i,ny,k)
          end do
        end do
      end if

    end if

    !-----------------------------------------


    !  Horizontal diffusion:


    do k=1,nzm

      rdx5=0.5*rdx2  * grdf_x(k,icrm)
      rdy5=0.5*rdy2  * grdf_y(k,icrm)

      do j=1,ny
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

      do j=0,ny
        jc=j+1
        do i=1,nx
          tky=rdy5*(tkh(i,j,k,icrm)+tkh(i,jc,k,icrm))
          flx(i,j,k)=-tky*(field(i,jc,k)-field(i,j,k))
        end do
      end do
      do j=1,ny
        jb=j-1
        do i=1,nx
          dfdt(i,j,k)=dfdt(i,j,k)-(flx(i,j,k)-flx(i,jb,k))
        end do
      end do

    end do ! k


    !  Vertical diffusion:

    flux(1) = 0.
    tmp=1./adzw(nz,icrm)
    do j=1,ny
      do i=1,nx
        flx(i,j,0)=fluxb(i,j)*rdz*rhow(1,icrm)
        flx(i,j,nzm)=fluxt(i,j)*rdz*tmp*rhow(nz,icrm)
        flux(1) = flux(1) + flx(i,j,0)
      end do
    end do


    do k=1,nzm-1
      kc=k+1
      flux(kc)=0.
      rhoi = rhow(kc,icrm)/adzw(kc,icrm)
      rdz5=0.5*rdz2 * grdf_z(k,icrm)
      do j=1,ny
        do i=1,nx
          tkz=rdz5*(tkh(i,j,k,icrm)+tkh(i,j,kc,icrm))
          flx(i,j,k)=-tkz*(field(i,j,kc)-field(i,j,k))*rhoi
          flux(kc) = flux(kc) + flx(i,j,k)
        end do
      end do
    end do

    do k=1,nzm
      kb=k-1
      rhoi = 1./(adz(k,icrm)*rho(k,icrm))
      do j=1,ny
        do i=1,nx
          dfdt(i,j,k)=dtn*(dfdt(i,j,k)-(flx(i,j,k)-flx(i,j,kb))*rhoi)
          field(i,j,k)=field(i,j,k)+dfdt(i,j,k)
        end do
      end do
    end do

  end subroutine diffuse_scalar3D

end module diffuse_scalar3D_mod
