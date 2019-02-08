
module bound_duvdt_mod
  use params, only: asyncid
  implicit none

contains

  subroutine bound_duvdt(ncrms)
    ! Periodic boundary exchange
    use vars
    implicit none
    integer, intent(in) :: ncrms
    integer i,j,k,icrm

    !$acc parallel loop collapse(3) copy(dudt) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          dudt(icrm,nxp1,j,k,na) = dudt(icrm,1,j,k,na)
        end do
      end do
    end do

    if(RUN3D) then

      !$acc parallel loop collapse(3) copy(dvdt) async(asyncid)
      do icrm = 1 , ncrms
        do k=1,nzm
          do i=1,nx
            dvdt(icrm,i,nyp1,k,na) = dvdt(icrm,i,1,k,na)
          end do
        end do
      end do

    endif

  end subroutine bound_duvdt

end module bound_duvdt_mod
