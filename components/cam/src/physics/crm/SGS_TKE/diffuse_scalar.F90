module diffuse_scalar_mod
  use diffuse_scalar2D_mod
  use diffuse_scalar3D_mod
  implicit none

contains

  subroutine diffuse_scalar (ncrms,icrm,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_y,grdf_z,tkh,f,fluxb,fluxt,fdiff,flux,f2lediff,f2lediss,fwlediff,doit)
    use grid
    use vars, only: rho, rhow
    use params
    implicit none
    integer, intent(in) :: ncrms,icrm
    ! input:
    integer :: dimx1_d,dimx2_d,dimy1_d,dimy2_d
    real(crm_rknd) grdf_x(nzm,ncrms)! grid factor for eddy diffusion in x
    real(crm_rknd) grdf_y(nzm,ncrms)! grid factor for eddy diffusion in y
    real(crm_rknd) grdf_z(nzm,ncrms)! grid factor for eddy diffusion in z
    real(crm_rknd) tkh (dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm,ncrms) ! SGS eddy conductivity
    real(crm_rknd) f(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! scalar
    real(crm_rknd) fluxb(nx,ny)   ! bottom flux
    real(crm_rknd) fluxt(nx,ny)   ! top flux
    real(crm_rknd) flux(nz)
    real(crm_rknd) f2lediff(nz),f2lediss(nz),fwlediff(nz)
    real(crm_rknd) fdiff(nz)
    logical doit
    ! Local
    real(crm_rknd) df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)  ! scalar
    integer i,j,k

    !call t_startf ('diffuse_scalars')

    df(:,:,:) = f(:,:,:)

    if(RUN3D) then
      call diffuse_scalar3D (ncrms,icrm,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_y,grdf_z,f,fluxb,fluxt,tkh,rho,rhow,flux)
    else
      call diffuse_scalar2D (ncrms,icrm,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_y,       f,fluxb,fluxt,tkh,rho,rhow,flux)
    endif

    do k=1,nzm
      fdiff(k)=0.
      do j=1,ny
        do i=1,nx
          fdiff(k)=fdiff(k)+f(i,j,k)-df(i,j,k)
        end do
      end do
    end do

    !call t_stopf ('diffuse_scalars')

  end subroutine diffuse_scalar

end module diffuse_scalar_mod
