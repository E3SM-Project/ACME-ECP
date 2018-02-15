
module openacc_utils
  use params, only: crm_rknd
  implicit none

contains


  subroutine memzero_real8(a,n)
    real(8), intent(out) :: a(n)
    integer, intent(in ) :: n
    integer :: i
    !$acc parallel loop gang vector
    do i = 1 , n
      a(i) = 0
    enddo
  end subroutine memzero_real8

  subroutine memzero_crm_rknd(a,n)
    real(crm_rknd), intent(out) :: a(n)
    integer       , intent(in ) :: n
    integer :: i
    !$acc parallel loop gang vector
    do i = 1 , n
      a(i) = 0
    enddo
  end subroutine memzero_crm_rknd

  subroutine memzero_integer(a,n)
    integer, intent(out) :: a(n)
    integer, intent(in ) :: n
    integer :: i
    !$acc parallel loop gang vector
    do i = 1 , n
      a(i) = 0
    enddo
  end subroutine memzero_integer

  subroutine memzero_logical(a,n)
    logical, intent(out) :: a(n)
    integer, intent(in ) :: n
    integer :: i
    !$acc parallel loop gang vector
    do i = 1 , n
      a(i) = .false.
    enddo
  end subroutine memzero_logical

  subroutine memcopy_crm_rknd(dest,src,n)
    real(crm_rknd), intent(out) :: dest(n)
    real(crm_rknd), intent(in ) :: src (n)
    integer       , intent(in ) :: n
    integer :: i
    !$acc parallel loop gang vector
    do i = 1 , n
      dest(i) = src(i)
    enddo
  end subroutine memcopy_crm_rknd

  subroutine memcopy_integer(dest,src,n)
    integer, intent(out) :: dest(n)
    integer, intent(in ) :: src (n)
    integer, intent(in ) :: n
    integer :: i
    !$acc parallel loop gang vector
    do i = 1 , n
      dest(i) = src(i)
    enddo
  end subroutine memcopy_integer

  subroutine memcopy_logical(dest,src,n)
    logical, intent(out) :: dest(n)
    logical, intent(in ) :: src (n)
    integer, intent(in ) :: n
    integer :: i
    !$acc parallel loop gang vector
    do i = 1 , n
      dest(i) = src(i)
    enddo
  end subroutine memcopy_logical

end module openacc_utils
