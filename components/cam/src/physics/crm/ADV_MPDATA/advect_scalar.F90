module advect_scalar_mod
  use advect_scalar2D_mod
  use advect_scalar3D_mod
  implicit none

contains

  subroutine advect_scalar (f,fadv,flux,f2leadv,f2legrad,fwleadv,doit,ncrms)

    !     positively definite monotonic advection with non-oscillatory option

    use grid
    use vars, only: u, v, w, rho, rhow
    use params, only: docolumn, crm_rknd

    implicit none
    integer, intent(in) :: ncrms

    real(crm_rknd) f(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real(crm_rknd) flux(ncrms,nz), fadv(ncrms,nz)
    real(crm_rknd) f2leadv(ncrms,nz),f2legrad(ncrms,nz),fwleadv(ncrms,nz)
    logical doit

    real(crm_rknd) df(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    integer i,j,k,icrm

    if(docolumn) then
      flux = 0.
      return
    end if

    !call t_startf ('advect_scalars')

    df(:,:,:,:) = f(:,:,:,:)

    if(RUN3D) then
      call advect_scalar3D(f, u, v, w, rho, rhow, flux, ncrms)
    else
      call advect_scalar2D(f, u, w, rho, rhow, flux, ncrms)
    endif

    do k=1,nzm
      fadv(:,k)=0.
      do j=1,ny
        do i=1,nx
          fadv(:,k)=fadv(:,k)+f(:,i,j,k)-df(:,i,j,k)
        end do
      end do
    end do

    !call t_stopf ('advect_scalars')

  end subroutine advect_scalar

end module advect_scalar_mod
