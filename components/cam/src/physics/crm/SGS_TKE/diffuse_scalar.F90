module diffuse_scalar_mod
  use diffuse_scalar2D_mod
  use diffuse_scalar3D_mod
  implicit none

contains

  subroutine diffuse_scalar (dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_y,grdf_z,tkh,f,fluxb,fluxt,fdiff,flux,f2lediff,f2lediss,fwlediff,doit,ncrms)

    use grid
    use vars, only: rho, rhow
    use params
    implicit none
    integer, intent(in) :: ncrms

    ! input:
    integer :: dimx1_d,dimx2_d,dimy1_d,dimy2_d
    real(crm_rknd) grdf_x  (ncrms,nzm)! grid factor for eddy diffusion in x
    real(crm_rknd) grdf_y  (ncrms,nzm)! grid factor for eddy diffusion in y
    real(crm_rknd) grdf_z  (ncrms,nzm)! grid factor for eddy diffusion in z
    real(crm_rknd), pointer :: tkh(:,:,:,:) ! SGS eddy conductivity
    real(crm_rknd) f       (ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)	! scalar
    real(crm_rknd) fluxb   (ncrms,nx,ny)		! bottom flux
    real(crm_rknd) fluxt   (ncrms,nx,ny)		! top flux
    real(crm_rknd) flux    (ncrms,nz)
    real(crm_rknd) f2lediff(ncrms,nz)
    real(crm_rknd) f2lediss(ncrms,nz)
    real(crm_rknd) fwlediff(ncrms,nz)
    real(crm_rknd) fdiff   (ncrms,nz)
    logical doit
    ! Local
    real(crm_rknd), allocatable :: df(:,:,:,:)	! scalar
    integer i,j,k, icrm

    allocate(df(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm))

    !$acc enter data create(df) async(1)

    !call t_startf ('diffuse_scalars')

    !$acc parallel loop gang vector collapse(4) default(present) async(1)
    do k = 1 , nzm
      do j = dimy1_s,dimy2_s
        do i = dimx1_s,dimx2_s
          do icrm = 1 , ncrms
            df(icrm,i,j,k) = f(icrm,i,j,k)
          enddo
        enddo
      enddo
    enddo

    if(RUN3D) then
      call diffuse_scalar3D (dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_y,grdf_z,f,fluxb,fluxt,tkh,rho,rhow,flux,ncrms)
    else
      call diffuse_scalar2D (dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_y,       f,fluxb,fluxt,tkh,rho,rhow,flux,ncrms)
    endif

    !$acc parallel loop gang vector collapse(2) default(present) async(1)
    do k=1,nzm
      do icrm = 1 , ncrms
        fdiff(icrm,k)=0.
        !$acc loop seq
        do j=1,ny
          !$acc loop seq
          do i=1,nx
            fdiff(icrm,k)=fdiff(icrm,k)+f(icrm,i,j,k)-df(icrm,i,j,k)
          end do
        end do
      end do
    end do

    !call t_stopf ('diffuse_scalars')

    !$acc exit data delete(df) async(1)

    deallocate(df)

  end subroutine diffuse_scalar

end module diffuse_scalar_mod
