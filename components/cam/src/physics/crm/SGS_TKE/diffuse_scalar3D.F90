module diffuse_scalar3D_mod
  implicit none

contains

  subroutine diffuse_scalar3D (ncrms,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_y,grdf_z,field,fluxb,fluxt,tkh,rho,rhow,flux)

    use grid
    use params
    use task_util_mod, only: task_rank_to_index
    implicit none
    integer, intent(in) :: ncrms
    ! input
    integer :: dimx1_d,dimx2_d,dimy1_d,dimy2_d
    real(crm_rknd) grdf_x(ncrms,nzm)! grid factor for eddy diffusion in x
    real(crm_rknd) grdf_y(ncrms,nzm)! grid factor for eddy diffusion in y
    real(crm_rknd) grdf_z(ncrms,nzm)! grid factor for eddy diffusion in z
    real(crm_rknd) field(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! scalar
    real(crm_rknd) tkh(ncrms,0:nxp1,1-YES3D:nyp1,nzm) ! eddy conductivity
    real(crm_rknd) fluxb(nx,ny,ncrms)   ! bottom flux
    real(crm_rknd) fluxt(nx,ny,ncrms)   ! top flux
    real(crm_rknd) rho(nzm,ncrms)
    real(crm_rknd) rhow(nz,ncrms)
    real(crm_rknd) flux(nz,ncrms)
    ! local
    real(crm_rknd) flx_x(0:nx,0:ny,0:nzm,ncrms), flx_y(0:nx,0:ny,0:nzm,ncrms), flx_z(0:nx,0:ny,0:nzm,ncrms)
    real(crm_rknd) dfdt(nx,ny,nz,ncrms)
    real(crm_rknd) rdx2,rdy2,rdz2,rdz,rdx5,rdy5,rdz5,tmp
    real(crm_rknd) dxy,dyx,tkx,tky,tkz,rhoi
    integer i,j,k,ib,ic,jb,jc,kc,kb,icrm

    if(.not.dosgs) return

    rdx2=1./(dx*dx)
    rdy2=1./(dy*dy)
    dxy=dx/dy
    dyx=dy/dx

    !$acc enter data create(flx_x,flx_y,flx_z,dfdt) async(asyncid)

    !$acc parallel loop collapse(4) copy(dfdt) async(asyncid)
    do icrm = 1 , ncrms
      do k = 1 , nzm
        do j = 1 , ny
          do i = 1 , nx
            dfdt(i,j,k,icrm)=0.
          enddo
        enddo
      enddo
    enddo

    !-----------------------------------------
    if(dowallx) then
      if(mod(rank,nsubdomains_x).eq.0) then
        !$acc parallel loop collapse(3) copy(field) async(asyncid)
        do icrm = 1 , ncrms
          do k=1,nzm
            do j=1,ny
              field(icrm,0,j,k) = field(icrm,1,j,k)
            enddo
          enddo
        enddo
      endif
      if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
        !$acc parallel loop collapse(3) copy(field) async(asyncid)
        do icrm = 1 , ncrms
          do k=1,nzm
            do j=1,ny
              field(icrm,nx+1,j,k) = field(icrm,nx,j,k)
            enddo
          enddo
        enddo
      endif
    endif

    if(dowally) then
      if(rank.lt.nsubdomains_x) then
        !$acc parallel loop collapse(3) copy(field) async(asyncid)
        do icrm = 1 , ncrms
          do k=1,nzm
            do i=1,nx
              field(icrm,i,1-YES3D,k) = field(icrm,i,1,k)
            enddo
          enddo
        enddo
      endif
      if(rank.gt.nsubdomains-nsubdomains_x-1) then
        !$acc parallel loop collapse(3) copy(field) async(asyncid)
        do icrm = 1 , ncrms
          do k=1,nzm
            do i=1,ny
              field(icrm,i,ny+YES3D,k) = field(icrm,i,ny,k)
            enddo
          enddo
        enddo
      endif
    endif

    if(dowally) then
      !$acc parallel loop collapse(3) copy(field) async(asyncid)
      do icrm = 1 , ncrms
        do k=1,nzm
          do i=1,nx
            field(icrm,i,1-YES3D,k) = field(icrm,i,1,k)
          enddo
        enddo
      enddo
      !$acc parallel loop collapse(3) copy(field) async(asyncid)
      do icrm = 1 , ncrms
        do k=1,nzm
          do i=1,nx
            field(icrm,i,ny+YES3D,k) = field(icrm,i,ny,k)
          enddo
        enddo
      enddo
    endif

    !  Horizontal diffusion:
    !$acc parallel loop collapse(4) copyin(field,tkh,grdf_x,grdf_y) copy(flx_x,flx_y) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=0,ny
          do i=0,nx
            if (j >= 1) then
              ic=i+1
              rdx5=0.5*rdx2  * grdf_x(icrm,k)
              tkx=rdx5*(tkh(icrm,i,j,k)+tkh(icrm,ic,j,k))
              flx_x(i,j,k,icrm)=-tkx*(field(icrm,ic,j,k)-field(icrm,i,j,k))
            endif
            if (i >= 1) then
              jc=j+1
              rdy5=0.5*rdy2  * grdf_y(icrm,k)
              tky=rdy5*(tkh(icrm,i,j,k)+tkh(icrm,i,jc,k))
              flx_y(i,j,k,icrm)=-tky*(field(icrm,i,jc,k)-field(icrm,i,j,k))
            endif
          enddo
        enddo
      enddo
    enddo
    !$acc parallel loop collapse(4) copyin(flx_x,flx_y) copy(dfdt) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            ib=i-1
            dfdt(i,j,k,icrm)=dfdt(i,j,k,icrm)-(flx_x(i,j,k,icrm)-flx_x(ib,j,k,icrm))
            jb=j-1
            dfdt(i,j,k,icrm)=dfdt(i,j,k,icrm)-(flx_y(i,j,k,icrm)-flx_y(i,jb,k,icrm))
          enddo
        enddo
      enddo ! k
    enddo

    !  Vertical diffusion:
    !$acc parallel loop collapse(2) copy(flux) async(asyncid)
    do icrm = 1 , ncrms
      do k = 1 , nzm
        flux(k,icrm) = 0.
      enddo
    enddo

    !$acc parallel loop collapse(4) copyin(rhow,adzw,dz,grdf_z,tkh,field,fluxb,fluxt) copy(flx_z,flux) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            if (k <= nzm-1) then
              kc=k+1
              rhoi = rhow(kc,icrm)/adzw(icrm,kc)
              rdz2=1./(dz(icrm)*dz(icrm))
              rdz5=0.5*rdz2 * grdf_z(icrm,k)
              tkz=rdz5*(tkh(icrm,i,j,k)+tkh(icrm,i,j,kc))
              flx_z(i,j,k,icrm)=-tkz*(field(icrm,i,j,kc)-field(icrm,i,j,k))*rhoi
              !$acc atomic update
              flux(kc,icrm) = flux(kc,icrm) + flx_z(i,j,k,icrm)
            elseif (k == nzm) then
              tmp=1./adzw(icrm,nz)
              rdz=1./dz(icrm)
              flx_z(i,j,0,icrm)=fluxb(i,j,icrm)*rdz*rhow(1,icrm)
              flx_z(i,j,nzm,icrm)=fluxt(i,j,icrm)*rdz*tmp*rhow(nz,icrm)
              !$acc atomic update
              flux(1,icrm) = flux(1,icrm) + flx_z(i,j,0,icrm)
            endif
          enddo
        enddo
      enddo
    enddo

    !$acc parallel loop collapse(4) copyin(adz,rho,flx_z) copy(dfdt,field) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            kb=k-1
            rhoi = 1./(adz(icrm,k)*rho(k,icrm))
            dfdt(i,j,k,icrm)=dtn*(dfdt(i,j,k,icrm)-(flx_z(i,j,k,icrm)-flx_z(i,j,kb,icrm))*rhoi)
            field(icrm,i,j,k)=field(icrm,i,j,k)+dfdt(i,j,k,icrm)
          enddo
        enddo
      enddo
    enddo

    !$acc exit data delete(flx_x,flx_y,flx_z,dfdt) async(asyncid)

  end subroutine diffuse_scalar3D

end module diffuse_scalar3D_mod
