module zero_mod
  implicit none

contains

  subroutine zero(ncrms)

    use vars

    implicit none
    integer, intent(in) :: ncrms
    integer :: i,j,k,icrm

    !$acc parallel loop gang vector collapse(4) default(present) async(1)
    do k = 1 , nz
      do j = 1 , nyp1
        do i = 1 , nxp1
          do icrm = 1 , ncrms
            if (j <= ny .and. k <= nzm) dudt(icrm,i,j,k,na) = 0.
            if (i <= nx .and. k <= nzm) dvdt(icrm,i,j,k,na) = 0.
            if (i <= nx .and. j <= ny ) dwdt(icrm,i,j,k,na) = 0.
            if (i <= nx .and. j <= ny ) misc(icrm,i,j,k) = 0.
          enddo
        enddo
      enddo
    enddo

  end

end module zero_mod
