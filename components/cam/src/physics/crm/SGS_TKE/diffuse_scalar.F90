module diffuse_scalar_mod
  use diffuse_scalar2D_mod
  use diffuse_scalar3D_mod
  implicit none

contains

  subroutine diffuse_scalar (ncrms,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_y,grdf_z,tkh,f,fluxb,fluxt,fdiff,flux)
    use grid
    use vars, only: rho, rhow
    use params
    implicit none
    integer, intent(in) :: ncrms
    ! input:
    integer :: dimx1_d,dimx2_d,dimy1_d,dimy2_d
    real(crm_rknd) grdf_x(nzm,ncrms)! grid factor for eddy diffusion in x
    real(crm_rknd) grdf_y(nzm,ncrms)! grid factor for eddy diffusion in y
    real(crm_rknd) grdf_z(nzm,ncrms)! grid factor for eddy diffusion in z
    real(crm_rknd) tkh (dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm,ncrms) ! SGS eddy conductivity
    real(crm_rknd) f(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, ncrms) ! scalar
    real(crm_rknd) fluxb(nx,ny,ncrms)   ! bottom flux
    real(crm_rknd) fluxt(nx,ny,ncrms)   ! top flux
    real(crm_rknd) fdiff(nz,ncrms)
    real(crm_rknd) flux (nz,ncrms)
    ! Local
    real(crm_rknd) df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm,ncrms)  ! scalar
    integer i,j,k,icrm

    do icrm = 1 , ncrms

    df(:,:,:,icrm) = f(:,:,:,icrm)

    if(RUN3D) then
      call diffuse_scalar3D (ncrms,icrm,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_y,grdf_z,f(:,:,:,icrm),fluxb(:,:,icrm),fluxt(:,:,icrm),tkh,rho,rhow,flux(:,icrm))
    else
      call diffuse_scalar2D (ncrms,icrm,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_y,       f(:,:,:,icrm),fluxb(:,:,icrm),fluxt(:,:,icrm),tkh,rho,rhow,flux(:,icrm))
    endif

    do k=1,nzm
      fdiff(k,icrm)=0.
      do j=1,ny
        do i=1,nx
          fdiff(k,icrm)=fdiff(k,icrm)+f(i,j,k,icrm)-df(i,j,k,icrm)
        end do
      end do
    end do

    enddo

  end subroutine diffuse_scalar

end module diffuse_scalar_mod
