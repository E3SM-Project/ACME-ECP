module bound_exchange_mod
  implicit none

contains

  subroutine bound_exchange(f,dimx1,dimx2,dimy1,dimy2,dimz,i_1, i_2, j_1, j_2, id, ncrms)
    ! periodic boundary exchange
    use grid
    use params, only: crm_rknd
    use openacc_pool
    implicit none
    integer, intent(in) :: ncrms

    integer dimx1, dimx2, dimy1, dimy2, dimz
    integer i_1, i_2, j_1, j_2
    real(crm_rknd) f(ncrms,dimx1:dimx2, dimy1:dimy2, dimz)
    integer id   ! id of the sent field (dummy variable)
    real(crm_rknd), pointer :: buffer(:)	! buffer for sending data

    integer i, j, k, n, icrm
    integer i1, i2, j1, j2

    call pool_push(buffer,(/ncrms*(nx+ny)*3*nz/))
    !allocate(buffer(ncrms*(nx+ny)*3*nz))

    i1 = i_1 - 1
    i2 = i_2 - 1
    j1 = j_1 - 1
    j2 = j_2 - 1

    !----------------------------------------------------------------------
    !  Send buffers to neighbors
    !----------------------------------------------------------------------


    if(RUN3D) then

      ! "North" -> "South":

      n=0
      do k=1,dimz
        do j=ny-j1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              n = n+1
              buffer(n) = f(icrm,i,j,k)
            end do
          end do
        end do
      end do
      n=0
      do k=1,dimz
        do j=-j1,0
          do i=1,nx
            do icrm = 1 , ncrms
              n = n+1
              f(icrm,i,j,k) = buffer(n)
            end do
          end do
        end do
      end do

      ! "North-East" -> "South-West":

      n=0
      do k=1,dimz
        do j=ny-j1,ny
          do i=nx-i1,nx
            do icrm = 1 , ncrms
              n = n+1
              buffer(n) = f(icrm,i,j,k)
            end do
          end do
        end do
      end do
      n=0
      do k=1,dimz
        do j=-j1,0
          do i=-i1,0
            do icrm = 1 , ncrms
              n = n+1
              f(icrm,i,j,k) = buffer(n)
            end do
          end do
        end do
      end do

      ! "South-East" -> "North-West":

      n=0
      do k=1,dimz
        do j=1,1+j2
          do i=nx-i1,nx
            do icrm = 1 , ncrms
              n = n+1
              buffer(n) = f(icrm,i,j,k)
            end do
          end do
        end do
      end do
      n=0
      do k=1,dimz
        do j=nyp1,nyp1+j2
          do i=-i1,0
            do icrm = 1 , ncrms
              n = n+1
              f(icrm,i,j,k) = buffer(n)
            end do
          end do
        end do
      end do

      ! "South" -> "North":

      n=0
      do k=1,dimz
        do j=1,1+j2
          do i=1,nx
            do icrm = 1 , ncrms
              n = n+1
              buffer(n) = f(icrm,i,j,k)
            end do
          end do
        end do
      end do
      n=0
      do k=1,dimz
        do j=nyp1,nyp1+j2
          do i=1,nx
            do icrm = 1 , ncrms
              n = n+1
              f(icrm,i,j,k) = buffer(n)
            end do
          end do
        end do
      end do

      ! "South-West" -> "North-East":

      n=0
      do k=1,dimz
        do j=1,1+j2
          do i=1,1+i2
            do icrm = 1 , ncrms
              n = n+1
              buffer(n) = f(icrm,i,j,k)
            end do
          end do
        end do
      end do
      n=0
      do k=1,dimz
        do j=nyp1,nyp1+j2
          do i=nxp1,nxp1+i2
            do icrm = 1 , ncrms
              n = n+1
              f(icrm,i,j,k) = buffer(n)
            end do
          end do
        end do
      end do


      ! To "North-West" -> "South-East":

      n=0
      do k=1,dimz
        do j=ny-j1,ny
          do i=1,1+i2
            do icrm = 1 , ncrms
              n = n+1
              buffer(n) = f(icrm,i,j,k)
            end do
          end do
        end do
      end do
      n=0
      do k=1,dimz
        do j=-j1,0
          do i=nxp1,nxp1+i2
            do icrm = 1 , ncrms
              n = n+1
              f(icrm,i,j,k) = buffer(n)
            end do
          end do
        end do
      end do


    endif

    !  "East" -> "West":

    !$acc parallel loop gang vector collapse(4)
    do k=1,dimz
      do j=1,ny
        do i=nx-i1,nx
          do icrm = 1 , ncrms
            n = (k-1)*ny*(i1+1)*ncrms + (j-1)*(i1+1)*ncrms + (i-(nx-i1))*ncrms + icrm
            buffer(n) = f(icrm,i,j,k)
          end do
        end do
      end do
    end do
    !$acc parallel loop gang vector collapse(4)
    do k=1,dimz
      do j=1,ny
        do i=-i1,0
          do icrm = 1 , ncrms
            n = (k-1)*ny*(i1+1)*ncrms + (j-1)*(i1+1)*ncrms + (i-(-i1))*ncrms + icrm
            f(icrm,i,j,k) = buffer(n)
          end do
        end do
      end do
    end do

    ! "West" -> "East":

    !$acc parallel loop gang vector collapse(4)
    do k=1,dimz
      do j=1,ny
        do i=1,1+i2
          do icrm = 1 , ncrms
            n = (k-1)*ny*(i2+1)*ncrms + (j-1)*(i2+1)*ncrms + (i-1)*ncrms + icrm
            buffer(n) = f(icrm,i,j,k)
          end do
        end do
      end do
    end do
    !$acc parallel loop gang vector collapse(4)
    do k=1,dimz
      do j=1,ny
        do i=nxp1,nxp1+i2
          do icrm = 1 , ncrms
            n = (k-1)*ny*(i2+1)*ncrms + (j-1)*(i2+1)*ncrms + (i-nxp1)*ncrms + icrm
            f(icrm,i,j,k) = buffer(n)
          end do
        end do
      end do
    end do


    ! deallocate(buffer)
    call pool_pop()


  end subroutine bound_exchange


end module bound_exchange_mod
