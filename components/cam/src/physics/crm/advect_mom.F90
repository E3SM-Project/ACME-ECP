module advect_mom_mod
  use advect2_mom_xy_mod
  use advect2_mom_z_mod
  implicit none

contains

  subroutine advect_mom(ncrms,icrm)
    use vars
    use params, only: docolumn
    implicit none
    integer, intent(in) :: ncrms,icrm
    integer i,j,k

    if(docolumn) return

    !call t_startf ('advect_mom')

    call advect2_mom_xy(ncrms,icrm)
    call advect2_mom_z(ncrms,icrm)

    !call t_stopf ('advect_mom')

  end subroutine advect_mom

end module advect_mom_mod
