module uvw_mod
  implicit none

contains

  subroutine uvw(ncrms,icrm)

    ! update the velocity field

    use vars
    use params
    implicit none
    integer, intent(in) :: ncrms,icrm

    u(1:nx,1:ny,1:nzm,icrm) = dudt(1:nx,1:ny,1:nzm,nc,icrm)
    v(1:nx,1:ny,1:nzm,icrm) = dvdt(1:nx,1:ny,1:nzm,nc,icrm)
    w(1:nx,1:ny,1:nzm,icrm) = dwdt(1:nx,1:ny,1:nzm,nc,icrm)

  end subroutine uvw

end module uvw_mod
