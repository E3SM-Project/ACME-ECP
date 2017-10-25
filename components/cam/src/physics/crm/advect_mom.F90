module advect_mom_mod
	use advect2_mom_xy_mod
	use advect2_mom_z_mod
	implicit none

contains

  subroutine advect_mom
    use vars
    use params, only: docolumn
    implicit none
    integer i,j,k

    if(docolumn) return

    !call t_startf ('advect_mom')

    call advect2_mom_xy()
    call advect2_mom_z()

    !call t_stopf ('advect_mom')

  end subroutine advect_mom

end module advect_mom_mod
