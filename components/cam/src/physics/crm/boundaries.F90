module boundaries_mod
  use periodic_mod
  use task_util_mod
  implicit none

contains

  subroutine boundaries(flag,ncrms)
    use grid, only: dompi
    implicit none
    integer, intent(in) :: ncrms
    integer flag

    !call t_startf ('boundaries')

    if(dompi) then
      call task_boundaries(flag)
    else
      call periodic(flag,ncrms)
    end if

    !call t_stopf ('boundaries')

  end subroutine boundaries
end module boundaries_mod
