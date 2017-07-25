
module cam_abortutils
  implicit none
contains

  subroutine endrun(s)
    implicit none
    character(len=*), intent(in) :: s
    write(*,*) s
    stop
  end subroutine endrun
end module cam_abortutils
