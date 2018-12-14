module advect_scalar_mod
  use advect_scalar2D_mod
  use advect_scalar3D_mod
  implicit none

contains

  subroutine advect_scalar (ncrms,f,fadv,flux)

    !     positively definite monotonic advection with non-oscillatory option

    use grid
    use vars, only: u, v, w, rho, rhow
    use params, only: docolumn, crm_rknd

    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) f(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm,ncrms)
    real(crm_rknd) flux(nz,ncrms), fadv(nz,ncrms)
    real(crm_rknd) f0(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm,ncrms)
    integer i,j,k,icrm

    if(docolumn) then
      do icrm = 1 , ncrms
        flux(:,icrm) = 0.
      enddo
      return
    end if

    do icrm = 1 , ncrms
      f0(:,:,:,icrm) = f(:,:,:,icrm)
    enddo

    if(RUN3D) then
      call advect_scalar3D(ncrms, f, u, v, w, rho, rhow, flux)
    else
      call advect_scalar2D(ncrms, f, u, w, rho, rhow, flux)
    endif

    do icrm = 1 , ncrms
      do k=1,nzm
        fadv(k,icrm)=0.
        do j=1,ny
          do i=1,nx
            fadv(k,icrm)=fadv(k,icrm)+f(i,j,k,icrm)-f0(i,j,k,icrm)
          end do
        end do
      end do
    enddo

  end subroutine advect_scalar

end module advect_scalar_mod
