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
    use openacc_pool
    implicit none
    integer, intent(in) :: ncrms

    real(crm_rknd) f(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real(crm_rknd) flux(ncrms,nz), fadv(ncrms,nz)
    real(crm_rknd) f2leadv(ncrms,nz),f2legrad(ncrms,nz),fwleadv(ncrms,nz)
    logical doit

    real(crm_rknd), pointer :: df(:,:,:,:)
    integer i,j,k,icrm

    call pool_push(df,(/1,dimx1_s,dimy1_s,1/),(/ncrms,dimx2_s,dimy2_s,nzm/))
    !allocate( df(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) )

    if(docolumn) then
      flux = 0.
      return
    end if

    !call t_startf ('advect_scalars')

    !$acc parallel loop gang vector collapse(4)
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
      call advect_scalar3D(f, u, v, w, rho, rhow, flux, ncrms)
    else
      call advect_scalar2D(f, u, w, rho, rhow, flux, ncrms)
    endif

    !$acc parallel loop gang vector collapse(2)
    do k=1,nzm
      do icrm = 1 , ncrms
        fadv(icrm,k)=0.
        !$acc loop seq
        do j=1,ny
          !$acc loop seq
          do i=1,nx
            fadv(icrm,k)=fadv(icrm,k)+f(icrm,i,j,k)-df(icrm,i,j,k)
          end do
        end do
      end do
    end do

    !call t_stopf ('advect_scalars')
    !deallocate( df )
    call pool_pop()

  end subroutine advect_scalar

end module advect_scalar_mod
