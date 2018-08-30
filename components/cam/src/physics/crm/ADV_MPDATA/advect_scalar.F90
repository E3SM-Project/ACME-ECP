module advect_scalar_mod
  use advect_scalar2D_mod
  use advect_scalar3D_mod
  implicit none

contains

  subroutine advect_scalar (ncrms,icrm,f,fadv,flux,f2leadv,f2legrad,fwleadv,doit)

    !     positively definite monotonic advection with non-oscillatory option

    use grid
    use vars, only: u, v, w, rho, rhow
    use params, only: docolumn, crm_rknd

    implicit none
    integer, intent(in) :: ncrms,icrm
    real(crm_rknd) f(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real(crm_rknd) flux(nz), fadv(nz)
    real(crm_rknd) f2leadv(nz),f2legrad(nz),fwleadv(nz)
    logical doit

    real(crm_rknd) df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    integer i,j,k

    if(docolumn) then
      flux = 0.
      return
    end if

    !call t_startf ('advect_scalars')

    df(:,:,:) = f(:,:,:)

    if(RUN3D) then
      call advect_scalar3D(ncrms, icrm, f, u, v, w, rho, rhow, flux)
    else
      call advect_scalar2D(ncrms, icrm, f, u, w, rho, rhow, flux)
    endif

    do k=1,nzm
      fadv(k)=0.
      do j=1,ny
        do i=1,nx
          fadv(k)=fadv(k)+f(i,j,k)-df(i,j,k)
        end do
      end do
    end do

    !call t_stopf ('advect_scalars')

  end subroutine advect_scalar

end module advect_scalar_mod
