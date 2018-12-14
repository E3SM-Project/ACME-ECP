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
    real(crm_rknd) tmp
    integer i,j,k,icrm

    !$acc enter data create(f0) async(1)

    if(docolumn) then
      !$acc parallel loop collapse(2) copy(flux) async(1)
      do icrm = 1 , ncrms
        do k = 1 , nz
          flux(k,icrm) = 0.
        enddo
      enddo
      return
    end if

    !$acc parallel loop collapse(4) copyin(f) copy(f0) async(1)
    do icrm = 1 , ncrms
      do k = 1 , nzm
        do j = dimy1_s,dimy2_s
          do i = dimx1_s,dimx2_s
            f0(i,j,k,icrm) = f(i,j,k,icrm)
          enddo
        enddo
      enddo
    enddo

    if(RUN3D) then
      call advect_scalar3D(ncrms, f, u, v, w, rho, rhow, flux)
    else
      call advect_scalar2D(ncrms, f, u, w, rho, rhow, flux)
    endif

    !$acc parallel loop collapse(2) copy(fadv) async(1)
    do icrm = 1 , ncrms
      do k=1,nzm
        fadv(k,icrm)=0.
      enddo
    enddo
    !$acc parallel loop collapse(4) copyin(f,f0) copy(fadv) async(1)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            tmp = f(i,j,k,icrm)-f0(i,j,k,icrm)
            !$acc atomic update
            fadv(k,icrm)=fadv(k,icrm)+tmp
          end do
        end do
      end do
    enddo

    !$acc exit data delete(f0) async(1)

  end subroutine advect_scalar

end module advect_scalar_mod
