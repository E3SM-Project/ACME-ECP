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
    real(crm_rknd) tkh(ncrms,dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm) ! SGS eddy conductivity
    real(crm_rknd) f(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! scalar
    real(crm_rknd) fluxb(nx,ny,ncrms)   ! bottom flux
    real(crm_rknd) fluxt(nx,ny,ncrms)   ! top flux
    real(crm_rknd) fdiff(nz,ncrms)
    real(crm_rknd) flux (nz,ncrms)
    ! Local
    real(crm_rknd) df(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)  ! scalar
    real(crm_rknd) :: tmp
    integer i,j,k,icrm

    !$acc enter data create(df) async(asyncid)

    !$acc parallel loop collapse(4) copyin(f) copy(df) async(asyncid)
    do icrm = 1 , ncrms
      do k = 1 , nzm
        do j = dimy1_s , dimy2_s
          do i = dimx1_s , dimx2_s
            df(icrm,i,j,k) = f(icrm,i,j,k)
          enddo
        enddo
      enddo
    enddo

    if(RUN3D) then
      call diffuse_scalar3D (ncrms,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_y,grdf_z,f,fluxb,fluxt,tkh,rho,rhow,flux)
    else
      call diffuse_scalar2D (ncrms,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,       grdf_z,f,fluxb,fluxt,tkh,rho,rhow,flux)
    endif

    !$acc parallel loop collapse(2) copy(fdiff) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        fdiff(k,icrm)=0.
      enddo
    enddo
    !$acc parallel loop collapse(2) copyin(f,df) copy(fdiff) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            tmp = f(icrm,i,j,k)-df(icrm,i,j,k)
            !$acc atomic update
            fdiff(k,icrm)=fdiff(k,icrm)+tmp
          end do
        end do
      end do
    enddo

    !$acc exit data delete(df) async(asyncid)

  end subroutine diffuse_scalar

end module diffuse_scalar_mod
