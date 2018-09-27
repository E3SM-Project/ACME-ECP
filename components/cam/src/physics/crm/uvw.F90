module uvw_mod
  implicit none

contains

  subroutine uvw(ncrms)
    ! update the velocity field
    use vars
    use params
    implicit none
    integer, intent(in) :: ncrms
    integer :: i, j, k, icrm

    do icrm = 1 , ncrms
      do k = 1 , nzm
        do j = 1 , ny
          do i = 1 , nx
            u(i,j,k,icrm) = dudt(i,j,k,nc(icrm),icrm)
            v(i,j,k,icrm) = dvdt(i,j,k,nc(icrm),icrm)
            w(i,j,k,icrm) = dwdt(i,j,k,nc(icrm),icrm)
          enddo
        enddo
      enddo
    enddo

  end subroutine uvw

end module uvw_mod
