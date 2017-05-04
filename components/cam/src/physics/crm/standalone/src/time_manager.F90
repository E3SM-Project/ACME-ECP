
module time_manager
  implicit none

contains

  integer function get_nstep()
    !character(len=*), parameter :: sub = 'get_nstep'
    !integer :: rc
    !integer(ESMF_KIND_I8) :: step_no
    !call ESMF_ClockGet(tm_clock, advanceCount=step_no, rc=rc)
    !call chkrc(rc, sub//': error return from ESMF_ClockGet')
    !get_nstep = step_no

    !TODO: I don't think this should need changing.
    get_nstep = 1
  end function get_nstep


end module time_manager
