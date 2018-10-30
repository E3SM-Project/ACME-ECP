module unittest_mod
  implicit none

  integer, parameter :: nx = 64
  integer, parameter :: ny = 16
  integer, parameter :: nzm = 57
  integer, parameter :: rc = selected_real_kind(12) ! 4 byte real
  integer, parameter :: crm_rknd = selected_real_kind(12) ! 4 byte real
  integer, parameter :: shr_kind_r8 = selected_real_kind(12) ! 4 byte real
  integer, parameter :: index_water_vapor = 1

  integer, parameter :: iulog = 6

  real(crm_rknd), dimension(nx, ny, nzm) :: t ! liquid/ice water static energy
  real(crm_rknd), dimension(nx, ny, nzm) :: qv, qcl, qci
  real(crm_rknd), dimension(nx, ny, nzm, 2) :: micro_field
  real(crm_rknd), dimension(nx, ny, nzm) :: u, v ! liquid/ice water static energy
  real(crm_rknd), dimension(nzm) :: t0, q0, u0, v0

contains
  subroutine endrun(message)
    implicit none
    character(len=*), intent(in) :: message

    write(iulog,*) message
  end subroutine endrun

  ! subroutine initialize()
  !   ! initialize variables
  !   t(:, :, :) = 0.
  !   qv(:, :, :) = 0.
  !   qcl(:, :, :) = 0.
  !   qci(:, :, :) = 0.
  !   micro_field(:, :, :, :) = 0.
  !   u(:, :, :) = 0.
  !   v(:, :, :) = 0.
  !
  !   t0(:) = 0.
  !   q0(:) = 0.
  !   u0(:) = 0.
  !   v0(:) = 0.
  ! end subroutine initialize
end module
