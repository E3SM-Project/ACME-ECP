module bound_exchange_mod
  implicit none

contains

  subroutine bound_exchange(ncrms,f,dimx1,dimx2,dimy1,dimy2,dimz,i_1, i_2, j_1, j_2, id)
    ! periodic boundary exchange
    use grid
    use params, only: crm_rknd
    implicit none
    integer dimx1, dimx2, dimy1, dimy2, dimz, ncrms
    integer i_1, i_2, j_1, j_2
    real(crm_rknd) f(dimx1:dimx2, dimy1:dimy2, dimz, ncrms)
    integer id   ! id of the sent field (dummy variable)
    real(crm_rknd) buffer((nx+ny)*3*nz*ncrms)  ! buffer for sending data
    integer i, j, k, n, icrm
    integer i1, i2, j1, j2

    !$acc enter data create(buffer) async(1)

    i1 = i_1 - 1
    i2 = i_2 - 1
    j1 = j_1 - 1
    j2 = j_2 - 1

    !----------------------------------------------------------------------
    !  Send buffers to neighbors
    !----------------------------------------------------------------------

    if(RUN3D) then
      ! "North" -> "South":
      !$acc parallel loop collapse(4) async(1)
      do icrm = 1 , ncrms
        do k=1,dimz
          do j=ny-j1,ny
            do i=1,nx
              n = (icrm-1)*dimz*(j1+1)*nx + (k-1)*(j1+1)*nx + (j-(ny-j1))*nx + (i-1) + 1
              buffer(n) = f(i,j,k,icrm)
            end do
          end do
        end do
      end do
      !$acc parallel loop collapse(4) async(1)
      do icrm = 1 , ncrms
        do k=1,dimz
          do j=-j1,0
            do i=1,nx
              n = (icrm-1)*dimz*(j1+1)*nx + (k-1)*(j1+1)*nx + (j-(-j1))*nx + (i-1) + 1
              f(i,j,k,icrm) = buffer(n)
            end do
          end do
        end do
      end do

      ! "North-East" -> "South-West":
      !$acc parallel loop collapse(4) async(1)
      do icrm = 1 , ncrms
        do k=1,dimz
          do j=ny-j1,ny
            do i=nx-i1,nx
              n = (icrm-1)*dimz*(j1+1)*(i1+1) + (k-1)*(j1+1)*(i1+1) + (j-(ny-j1))*(i1+1) + (i-(nx-i1)) + 1
              buffer(n) = f(i,j,k,icrm)
            end do
          end do
        end do
      end do
      !$acc parallel loop collapse(4) async(1)
      do icrm = 1 , ncrms
        do k=1,dimz
          do j=-j1,0
            do i=-i1,0
              n = (icrm-1)*dimz*(j1+1)*(i1+1) + (k-1)*(j1+1)*(i1+1) + (j-(-j1))*(i1+1) + (i-(-i1)) + 1
              f(i,j,k,icrm) = buffer(n)
            end do
          end do
        end do
      end do

      ! "South-East" -> "North-West":
      !$acc parallel loop collapse(4) async(1)
      do icrm = 1 , ncrms
        do k=1,dimz
          do j=1,1+j2
            do i=nx-i1,nx
              n = (icrm-1)*dimz*(j2+1)*(i1+1) + (k-1)*(j2+1)*(i1+1) + (j-1)*(i1+1) + (i-(nx-i1)) + 1
              buffer(n) = f(i,j,k,icrm)
            end do
          end do
        end do
      end do
      !$acc parallel loop collapse(4) async(1)
      do icrm = 1 , ncrms
        do k=1,dimz
          do j=nyp1,nyp1+j2
            do i=-i1,0
              n = (icrm-1)*dimz*(j2+1)*(i1+1) + (k-1)*(j2+1)*(i1+1) + (j-nyp1)*(i1+1) + (i-(-i1)) + 1
              f(i,j,k,icrm) = buffer(n)
            end do
          end do
        end do
      end do

      ! "South" -> "North":
      !$acc parallel loop collapse(4) async(1)
      do icrm = 1 , ncrms
        do k=1,dimz
          do j=1,1+j2
            do i=1,nx
              n = (icrm-1)*dimz*(j2+1)*nx + (k-1)*(j2+1)*nx + (j-1)*nx + (i-1) + 1
              buffer(n) = f(i,j,k,icrm)
            end do
          end do
        end do
      end do
      !$acc parallel loop collapse(4) async(1)
      do icrm = 1 , ncrms
        do k=1,dimz
          do j=nyp1,nyp1+j2
            do i=1,nx
              n = (icrm-1)*dimz*(j2+1)*nx + (k-1)*(j2+1)*nx + (j-nyp1)*nx + (i-1) + 1
              f(i,j,k,icrm) = buffer(n)
            end do
          end do
        end do
      end do

      ! "South-West" -> "North-East":
      !$acc parallel loop collapse(4) async(1)
      do icrm = 1 , ncrms
        do k=1,dimz
          do j=1,1+j2
            do i=1,1+i2
              n = (icrm-1)*dimz*(j2+1)*(i2+1) + (k-1)*(j2+1)*(i2+1) + (j-1)*(i2+1) + (i-1) + 1
              buffer(n) = f(i,j,k,icrm)
            end do
          end do
        end do
      end do
      !$acc parallel loop collapse(4) async(1)
      do icrm = 1 , ncrms
        do k=1,dimz
          do j=nyp1,nyp1+j2
            do i=nxp1,nxp1+i2
              n = (icrm-1)*dimz*(j2+1)*(i2+1) + (k-1)*(j2+1)*(i2+1) + (j-nyp1)*(i2+1) + (i-nxp1) + 1
              f(i,j,k,icrm) = buffer(n)
            end do
          end do
        end do
      end do


      ! To "North-West" -> "South-East":
      !$acc parallel loop collapse(4) async(1)
      do icrm = 1 , ncrms
        do k=1,dimz
          do j=ny-j1,ny
            do i=1,1+i2
              n = (icrm-1)*dimz*(j1+1)*(i2+1) + (k-1)*(j1+1)*(i2+1) + (j-(ny-j1))*(i2+1) + (i-1) + 1
              buffer(n) = f(i,j,k,icrm)
            end do
          end do
        end do
      end do
      !$acc parallel loop collapse(4) async(1)
      do icrm = 1 , ncrms
        do k=1,dimz
          do j=-j1,0
            do i=nxp1,nxp1+i2
              n = (icrm-1)*dimz*(j1+1)*(i2+1) + (k-1)*(j1+1)*(i2+1) + (j-(-j1))*(i2+1) + (i-nxp1) + 1
              f(i,j,k,icrm) = buffer(n)
            end do
          end do
        end do
      end do


    endif

    !  "East" -> "West":
    !$acc parallel loop collapse(4) async(1)
    do icrm = 1 , ncrms
      do k=1,dimz
        do j=1,ny
          do i=nx-i1,nx
            n = (icrm-1)*dimz*ny*(i1+1) + (k-1)*ny*(i1+1) + (j-1)*(i1+1) + (i-(nx-i1)) + 1
            buffer(n) = f(i,j,k,icrm)
          end do
        end do
      end do
    end do
    !$acc parallel loop collapse(4) async(1)
    do icrm = 1 , ncrms
      do k=1,dimz
        do j=1,ny
          do i=-i1,0
            n = (icrm-1)*dimz*ny*(i1+1) + (k-1)*ny*(i1+1) + (j-1)*(i1+1) + (i-(-i1)) + 1
            f(i,j,k,icrm) = buffer(n)
          end do
        end do
      end do
    end do

    ! "West" -> "East":
    !$acc parallel loop collapse(4) async(1)
    do icrm = 1 , ncrms
      do k=1,dimz
        do j=1,ny
          do i=1,1+i2
            n = (icrm-1)*dimz*ny*(i2+1) + (k-1)*ny*(i2+1) + (j-1)*(i2+1) + (i-1) + 1
            buffer(n) = f(i,j,k,icrm)
          end do
        end do
      end do
    end do
    !$acc parallel loop collapse(4) async(1)
    do icrm = 1 , ncrms
      do k=1,dimz
        do j=1,ny
          do i=nxp1,nxp1+i2
            n = (icrm-1)*dimz*ny*(i2+1) + (k-1)*ny*(i2+1) + (j-1)*(i2+1) + (i-nxp1) + 1
            f(i,j,k,icrm) = buffer(n)
          end do
        end do
      end do
    end do

    !$acc exit data delete(buffer) async(1)

  end subroutine bound_exchange


end module bound_exchange_mod
