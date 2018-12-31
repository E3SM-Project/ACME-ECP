
module bound_duvdt_mod
  implicit none

contains

  subroutine bound_duvdt(ncrms)
    ! Periodic boundary exchange
    use vars
    implicit none
    integer, intent(in) :: ncrms
    integer i,j,k,icrm

    !$acc parallel loop collapse(3) copy(dudt) async(1)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          dudt(nxp1,j,k,na,icrm) = dudt(1,j,k,na,icrm)
        end do
      end do
    end do

    if(RUN3D) then

      !$acc parallel loop collapse(3) copy(dvdt) async(1)
      do icrm = 1 , ncrms
        do k=1,nzm
          do i=1,nx
            dvdt(i,nyp1,k,na,icrm) = dvdt(i,1,k,na,icrm)
          end do
        end do
      end do

    endif

  end subroutine bound_duvdt

end module bound_duvdt_mod
