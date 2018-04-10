module diffuse_scalar2D_mod
  implicit none

contains
  subroutine diffuse_scalar2D (dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_z,field,fluxb,fluxt,tkh,rho,rhow,flux,ncrms)

    use grid
    use params
    use openacc_pool
    implicit none
    integer, intent(in) :: ncrms

    ! input
    integer :: dimx1_d,dimx2_d,dimy1_d,dimy2_d
    real(crm_rknd) grdf_x(ncrms,nzm)! grid factor for eddy diffusion in x
    real(crm_rknd) grdf_z(ncrms,nzm)! grid factor for eddy diffusion in z
    real(crm_rknd) field (ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)  ! scalar
    real(crm_rknd) tkh   (ncrms,0:nxp1, 1-YES3D:nyp1, nzm)  ! eddy conductivity
    real(crm_rknd) fluxb (ncrms,nx,ny)    ! bottom flux
    real(crm_rknd) fluxt (ncrms,nx,ny)    ! top flux
    real(crm_rknd) rho   (ncrms,nzm)
    real(crm_rknd) rhow  (ncrms,nz)
    real(crm_rknd) flux  (ncrms,nz)

    ! local
    real(crm_rknd), pointer :: flx (:,:,:,:)
    real(crm_rknd), pointer :: dfdt(:,:,:,:)
    real(crm_rknd) rdx2
    real(crm_rknd) tkx,tkz, tmp1, tmp2
    integer i,j,k,ib,ic,kc,kb, icrm

    if(.not.dosgs.and..not.docolumn) return

    ! allocate(flx (ncrms,0:nx,1,0:nzm))
    ! allocate(dfdt(ncrms,nx,ny,nzm))
    call pool_push(flx,(/1,0,1,0/),(/ncrms,nx,1,nzm/))
    call pool_push(dfdt,(/ncrms,nx,ny,nzm/))

    rdx2=1./(dx*dx)

    j=1

    !$acc parallel loop gang vector collapse(3)
    do k = 1 , nzm
      do i = 1 , nx
        do icrm = 1 , ncrms
          dfdt(icrm,i,j,k) = 0.
        enddo
      enddo
    enddo

    if(dowallx) then

      if(mod(rank,nsubdomains_x).eq.0) then
        !$acc parallel loop gang vector collapse(2)
        do k=1,nzm
          do icrm = 1 , ncrms
            field(icrm,0,j,k) = field(icrm,1,j,k)
          end do
        end do
      end if
      if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
        !$acc parallel loop gang vector collapse(2)
        do k=1,nzm
          do icrm = 1 , ncrms
            field(icrm,nx+1,j,k) = field(icrm,nx,j,k)
          end do
        end do
      end if

    end if

    if(.not.docolumn) then
      !$acc parallel loop gang vector collapse(3)
      do k=1,nzm
        do i=1,nx
          do icrm = 1 , ncrms
            ic=i+1
            ib=i-1
            tmp1 = -0.5*rdx2  *grdf_x(icrm,k)*(tkh(icrm,i ,j,k)+tkh(icrm,ic,j,k))*(field(icrm,ic,j,k)-field(icrm,i ,j,k))
            tmp2 = -0.5*rdx2  *grdf_x(icrm,k)*(tkh(icrm,ib,j,k)+tkh(icrm,i ,j,k))*(field(icrm,i ,j,k)-field(icrm,ib,j,k))
            dfdt(icrm,i,j,k)=dfdt(icrm,i,j,k)-(tmp1-tmp2)
          end do
        end do
      end do
    end if

    !$acc parallel loop gang vector collapse(2)
    do k = 1 , nz
      do icrm = 1 , ncrms
        flux(icrm,k) = 0.
      enddo
    enddo

    !$acc parallel loop gang vector collapse(3)
    do k=1,nzm
      do i=1,nx
        do icrm = 1 , ncrms
          kc=k+1
          if (k <= nzm-1) then
            flx(icrm,i,j,k)=-0.5*grdf_z(icrm,k)/(dz(icrm)*dz(icrm))*(tkh(icrm,i,j,k)+tkh(icrm,i,j,kc))*(field(icrm,i,j,kc)-field(icrm,i,j,k))*rhow(icrm,kc)/adzw(icrm,kc)
            !$acc atomic update
            flux(icrm,kc) = flux(icrm,kc) + flx(icrm,i,j,k)
          elseif (k == nzm) then
            flx(icrm,i,j,0)=fluxb(icrm,i,j)/dz(icrm)*rhow(icrm,1)
            flx(icrm,i,j,nzm)=fluxt(icrm,i,j)/dz(icrm)/adzw(icrm,nz)*rhow(icrm,nz)
            !$acc atomic update
            flux(icrm,1) = flux(icrm,1) + flx(icrm,i,j,0)
          endif
        end do
      end do
    end do

    !$acc parallel loop gang vector collapse(3)
    do k=1,nzm
      do i=1,nx
        do icrm = 1 , ncrms
          kb=k-1
          dfdt(icrm,i,j,k)=dtn*(dfdt(icrm,i,j,k)-(flx(icrm,i,j,k)-flx(icrm,i,j,kb))/(adz(icrm,k)*rho(icrm,k)))
          field(icrm,i,j,k)=field(icrm,i,j,k) + dfdt(icrm,i,j,k)
        end do
      end do
    end do

    ! deallocate(flx )
    ! deallocate(dfdt)
    call pool_pop_multiple(2)

  end subroutine diffuse_scalar2D

end module diffuse_scalar2D_mod
