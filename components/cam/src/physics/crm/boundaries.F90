module boundaries_mod
  use periodic_mod
  use task_util_mod
  implicit none

contains

  subroutine boundaries(ncrms,icrm,flag)
    use grid, only: dompi
    implicit none
    integer, intent(in) :: ncrms,icrm
    integer flag

    !call t_startf ('boundaries')

    if(dompi) then
      call task_boundaries(flag)
    else
      call periodic(ncrms,icrm,flag)
    end if

    !call t_stopf ('boundaries')

  end subroutine boundaries
end module boundaries_mod
