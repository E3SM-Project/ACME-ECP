module uvw_mod
  implicit none

contains

  subroutine uvw(ncrms,icrm)

    ! update the velocity field

    use vars
    use params
    implicit none
    integer, intent(in) :: ncrms,icrm

    u(icrm,1:nx,1:ny,1:nzm) = dudt(icrm,1:nx,1:ny,1:nzm,nc)
    v(icrm,1:nx,1:ny,1:nzm) = dvdt(icrm,1:nx,1:ny,1:nzm,nc)
    w(icrm,1:nx,1:ny,1:nzm) = dwdt(icrm,1:nx,1:ny,1:nzm,nc)

  end subroutine uvw

end  module uvw_mod
