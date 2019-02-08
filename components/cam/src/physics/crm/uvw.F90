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

    !$acc parallel loop collapse(4) copyin(dudt,dvdt,dwdt) copy(u,v,w) async(asyncid)
    do icrm = 1 , ncrms
      do k = 1 , nzm
        do j = 1 , ny
          do i = 1 , nx
            u(icrm,i,j,k) = dudt(i,j,k,nc,icrm)
            v(icrm,i,j,k) = dvdt(i,j,k,nc,icrm)
            w(icrm,i,j,k) = dwdt(i,j,k,nc,icrm)
          enddo
        enddo
      enddo
    enddo

  end subroutine uvw

end module uvw_mod
