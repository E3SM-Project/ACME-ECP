module diffuse_scalar3D_mod
  implicit none

contains

  subroutine diffuse_scalar3D (dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_y,grdf_z,field,fluxb,fluxt,tkh,rho,rhow,flux,ncrms)

    use grid
    use params
    use task_util_mod, only: task_rank_to_index
    implicit none
    integer, intent(in) :: ncrms
    ! input
    integer :: dimx1_d,dimx2_d,dimy1_d,dimy2_d
    real(crm_rknd) grdf_x(ncrms,nzm)! grid factor for eddy diffusion in x
    real(crm_rknd) grdf_y(ncrms,nzm)! grid factor for eddy diffusion in y
    real(crm_rknd) grdf_z(ncrms,nzm)! grid factor for eddy diffusion in z
    real(crm_rknd) field (ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)	! scalar
    real(crm_rknd), pointer :: tkh(:,:,:,:)	! eddy conductivity
    real(crm_rknd) fluxb (ncrms,nx,ny)		! bottom flux
    real(crm_rknd) fluxt (ncrms,nx,ny)		! top flux
    real(crm_rknd) rho   (ncrms,nzm)
    real(crm_rknd) rhow  (ncrms,nz)
    real(crm_rknd) flux  (ncrms,nz)
    ! local
    real(crm_rknd) flx (ncrms,0:nx,0:ny,0:nzm)
    real(crm_rknd) dfdt(ncrms,nx,ny,nz)
    real(crm_rknd) rdx2,rdy2,rdz2(ncrms),rdz(ncrms),rdx5(ncrms),rdy5(ncrms),rdz5(ncrms),tmp(ncrms)
    real(crm_rknd) dxy,dyx,tkx,tky,tkz,rhoi(ncrms)
    integer i,j,k,ib,ic,jb,jc,kc,kb, icrm


    if(.not.dosgs) return

    rdx2=1./(dx*dx)
    rdy2=1./(dy*dy)
    rdz2(:)=1./(dz(:)*dz(:))
    rdz(:)=1./dz(:)
    dxy=dx/dy
    dyx=dy/dx

    dfdt(:,:,:,:)=0.

    !-----------------------------------------
    if(dowallx) then

      if(mod(rank,nsubdomains_x).eq.0) then
        do k=1,nzm
          do j=1,ny
            do icrm = 1 , ncrms
              field(icrm,0,j,k) = field(icrm,1,j,k)
            end do
          end do
        end do
      end if
      if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
        do k=1,nzm
          do j=1,ny
            do icrm = 1 , ncrms
              field(icrm,nx+1,j,k) = field(icrm,nx,j,k)
            end do
          end do
        end do
      end if

    end if

    if(dowally) then

      if(rank.lt.nsubdomains_x) then
        do k=1,nzm
          do i=1,nx
            do icrm = 1 , ncrms
              field(icrm,i,1-YES3D,k) = field(icrm,i,1,k)
            end do
          end do
        end do
      end if
      if(rank.gt.nsubdomains-nsubdomains_x-1) then
        do k=1,nzm
          do i=1,ny
            do icrm = 1 , ncrms
              field(icrm,i,ny+YES3D,k) = field(icrm,i,ny,k)
            end do
          end do
        end do
      end if

    end if



    if(dowally) then

      call task_rank_to_index(rank, ib, jb)
      if(jb.eq.0) then
        do k=1,nzm
          do i=1,nx
            do icrm = 1 , ncrms
              field(icrm,i,1-YES3D,k) = field(icrm,i,1,k)
            end do
          end do
        end do
      end if
      if(jb.eq.nsubdomains_y-1) then
        do k=1,nzm
          do i=1,nx
            do icrm = 1 , ncrms
              field(icrm,i,ny+YES3D,k) = field(icrm,i,ny,k)
            end do
          end do
        end do
      end if

    end if

    !-----------------------------------------


    !  Horizontal diffusion:


    do k=1,nzm

      rdx5(:)=0.5*rdx2  * grdf_x(:,k)
      rdy5(:)=0.5*rdy2  * grdf_y(:,k)

      do j=1,ny
        do i=0,nx
          do icrm = 1 , ncrms
            ic=i+1
            tkx=rdx5(icrm)*(tkh(icrm,i,j,k)+tkh(icrm,ic,j,k))
            flx(icrm,i,j,k)=-tkx*(field(icrm,ic,j,k)-field(icrm,i,j,k))
          end do
        end do
        do i=1,nx
          do icrm = 1 , ncrms
            ib=i-1
            dfdt(icrm,i,j,k)=dfdt(icrm,i,j,k)-(flx(icrm,i,j,k)-flx(icrm,ib,j,k))
          end do
        end do
      end do

      do j=0,ny
        jc=j+1
        do i=1,nx
          do icrm = 1 , ncrms
            tky=rdy5(icrm)*(tkh(icrm,i,j,k)+tkh(icrm,i,jc,k))
            flx(icrm,i,j,k)=-tky*(field(icrm,i,jc,k)-field(icrm,i,j,k))
          end do
        end do
      end do
      do j=1,ny
        jb=j-1
        do i=1,nx
          do icrm = 1 , ncrms
            dfdt(icrm,i,j,k)=dfdt(icrm,i,j,k)-(flx(icrm,i,j,k)-flx(icrm,i,jb,k))
          end do
        end do
      end do

    end do ! k


    !  Vertical diffusion:

    flux(:,1) = 0.
    tmp(:)=1./adzw(:,nz)
    do j=1,ny
      do i=1,nx
        do icrm = 1 , ncrms
          flx(icrm,i,j,0)=fluxb(icrm,i,j)*rdz(icrm)*rhow(icrm,1)
          flx(icrm,i,j,nzm)=fluxt(icrm,i,j)*rdz(icrm)*tmp(icrm)*rhow(icrm,nz)
          flux(icrm,1) = flux(icrm,1) + flx(icrm,i,j,0)
        end do
      end do
    end do


    do k=1,nzm-1
      kc=k+1
      flux(:,kc)=0.
      rhoi(:) = rhow(:,kc)/adzw(:,kc)
      rdz5(:)=0.5*rdz2(:) * grdf_z(:,k)
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            tkz=rdz5(icrm)*(tkh(icrm,i,j,k)+tkh(icrm,i,j,kc))
            flx(icrm,i,j,k)=-tkz*(field(icrm,i,j,kc)-field(icrm,i,j,k))*rhoi(icrm)
            flux(icrm,kc) = flux(icrm,kc) + flx(icrm,i,j,k)
          end do
        end do
      end do
    end do

    do k=1,nzm
      kb=k-1
      rhoi(:) = 1./(adz(:,k)*rho(:,k))
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            dfdt(icrm,i,j,k)=dtn*(dfdt(icrm,i,j,k)-(flx(icrm,i,j,k)-flx(icrm,i,j,kb))*rhoi(icrm))
            field(icrm,i,j,k)=field(icrm,i,j,k)+dfdt(icrm,i,j,k)
          end do
        end do
      end do
    end do

  end subroutine diffuse_scalar3D

end module diffuse_scalar3D_mod
