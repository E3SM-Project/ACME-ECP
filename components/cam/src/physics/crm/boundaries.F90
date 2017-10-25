module boundaries_mod
	use periodic_mod
	use task_util_mod
	implicit none

contains

  subroutine boundaries(flag)
    use grid, only: dompi
    implicit none
    integer flag

    !call t_startf ('boundaries')

    if(dompi) then
      call task_boundaries(flag)
    else
      call periodic(flag)
    end if

    !call t_stopf ('boundaries')

  end subroutine boundaries
end module boundaries_mod
