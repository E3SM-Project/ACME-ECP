!>
! Module: mo_gas_optics

! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2015,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!
! Description: Reorder array indecies
!

module mo_util_reorder

  use mo_rte_kind, only: wp

  implicit none

  ! ----- interface for 2D arrays -----
  interface reorder12x21
    module procedure reorder_int_12x21, reorder_wp_12x21
  end interface

  ! ----- interface for 3D arrays -----
  interface reorder123x132
    module procedure reorder_int_123x132, reorder_wp_123x132
  end interface

  interface reorder123x213
    module procedure reorder_int_123x213, reorder_wp_123x213
  end interface

  interface reorder123x231
    module procedure reorder_int_123x231, reorder_wp_123x231
  end interface

  interface reorder123x312
    module procedure reorder_int_123x312, reorder_wp_123x312
  end interface

  interface reorder123x321
    module procedure reorder_int_123x321, reorder_wp_123x321
  end interface

  ! ----- interface for 4D arrays -----
  interface reorder1234x1243
    module procedure reorder_int_1234x1243, reorder_wp_1234x1243
  end interface

  interface reorder1234x1324
    module procedure reorder_int_1234x1324, reorder_wp_1234x1324
  end interface

  interface reorder1234x1342
    module procedure reorder_int_1234x1342, reorder_wp_1234x1342
  end interface

  interface reorder1234x1423
    module procedure reorder_int_1234x1423, reorder_wp_1234x1423
  end interface

  interface reorder1234x1432
    module procedure reorder_int_1234x1432, reorder_wp_1234x1432
  end interface

  interface reorder1234x2134
    module procedure reorder_int_1234x2134, reorder_wp_1234x2134
  end interface

  interface reorder1234x2143
    module procedure reorder_int_1234x2143, reorder_wp_1234x2143
  end interface

  interface reorder1234x2314
    module procedure reorder_int_1234x2314, reorder_wp_1234x2314
  end interface

  interface reorder1234x2341
    module procedure reorder_int_1234x2341, reorder_wp_1234x2341
  end interface

  interface reorder1234x2413
    module procedure reorder_int_1234x2413, reorder_wp_1234x2413
  end interface

  interface reorder1234x2431
    module procedure reorder_int_1234x2431, reorder_wp_1234x2431
  end interface

  interface reorder1234x3124
    module procedure reorder_int_1234x3124, reorder_wp_1234x3124
  end interface

  interface reorder1234x3142
    module procedure reorder_int_1234x3142, reorder_wp_1234x3142
  end interface

  interface reorder1234x3214
    module procedure reorder_int_1234x3214, reorder_wp_1234x3214
  end interface

  interface reorder1234x3241
    module procedure reorder_int_1234x3241, reorder_wp_1234x3241
  end interface

  interface reorder1234x3412
    module procedure reorder_int_1234x3412, reorder_wp_1234x3412
  end interface

  interface reorder1234x3421
    module procedure reorder_int_1234x3421, reorder_wp_1234x3421
  end interface

  interface reorder1234x4123
    module procedure reorder_int_1234x4123, reorder_wp_1234x4123
  end interface

  interface reorder1234x4132
    module procedure reorder_int_1234x4132, reorder_wp_1234x4132
  end interface

  interface reorder1234x4213
    module procedure reorder_int_1234x4213, reorder_wp_1234x4213
  end interface

  interface reorder1234x4231
    module procedure reorder_int_1234x4231, reorder_wp_1234x4231
  end interface

  interface reorder1234x4312
    module procedure reorder_int_1234x4312, reorder_wp_1234x4312
  end interface

  interface reorder1234x4321
    module procedure reorder_int_1234x4321, reorder_wp_1234x4321
  end interface

contains

  ! ----- reorder for 2D arrays -----
  ! reorder the indecies of 4D array
  function reorder_int_12x21(array)
    integer, dimension(:,:), intent(in) :: array
    integer, dimension(size(array,dim=2),size(array,dim=1)) :: reorder_int_12x21
    reorder_int_12x21(:,:) = transpose(array(:,:))
  end function

  ! reorder the indecies of 4D array
  function reorder_wp_12x21(array)
    real(wp), dimension(:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=2),size(array,dim=1)) :: reorder_wp_12x21
    reorder_wp_12x21(:,:) = transpose(array(:,:))
  end function

  ! ----- reorder for 3D arrays -----
  ! reorder the indecies of 4D array
  function reorder_int_123x132(array)
    integer, dimension(:,:,:), intent(in) :: array
    integer, dimension(size(array,dim=1),size(array,dim=3),size(array,dim=2)) :: reorder_int_123x132
    integer :: i2
    do i2 = 1, size(array,dim=2)
      reorder_int_123x132(:,:,i2) = array(:,i2,:)
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_wp_123x132(array)
    real(wp), dimension(:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=1),size(array,dim=3),size(array,dim=2)) :: reorder_wp_123x132
    integer :: i2
    do i2 = 1, size(array,dim=2)
      reorder_wp_123x132(:,:,i2) = array(:,i2,:)
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_int_123x213(array)
    integer, dimension(:,:,:), intent(in) :: array
    integer, dimension(size(array,dim=2),size(array,dim=1),size(array,dim=3)) :: reorder_int_123x213
    integer :: i3
    do i3 = 1, size(array,dim=3)
      reorder_int_123x213(:,:,i3) = transpose(array(:,:,i3))
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_wp_123x213(array)
    real(wp), dimension(:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=2),size(array,dim=1),size(array,dim=3)) :: reorder_wp_123x213
    integer :: i3
    do i3 = 1, size(array,dim=3)
      reorder_wp_123x213(:,:,i3) = transpose(array(:,:,i3))
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_int_123x231(array)
    integer, dimension(:,:,:), intent(in) :: array
    integer, dimension(size(array,dim=2),size(array,dim=3),size(array,dim=1)) :: reorder_int_123x231
    integer :: i1
    do i1 = 1, size(array,dim=1)
      reorder_int_123x231(:,:,i1) = array(i1,:,:)
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_wp_123x231(array)
    real(wp), dimension(:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=2),size(array,dim=3),size(array,dim=1)) :: reorder_wp_123x231
    integer :: i1
    do i1 = 1, size(array,dim=1)
      reorder_wp_123x231(:,:,i1) = array(i1,:,:)
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_int_123x312(array)
    integer, dimension(:,:,:), intent(in) :: array
    integer, dimension(size(array,dim=3),size(array,dim=1),size(array,dim=2)) :: reorder_int_123x312
    integer :: i2
    do i2 = 1, size(array,dim=2)
      reorder_int_123x312(:,:,i2) = transpose(array(:,i2,:))
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_wp_123x312(array)
    real(wp), dimension(:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=3),size(array,dim=1),size(array,dim=2)) :: reorder_wp_123x312
    integer :: i2
    do i2 = 1, size(array,dim=2)
      reorder_wp_123x312(:,:,i2) = transpose(array(:,i2,:))
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_int_123x321(array)
    integer, dimension(:,:,:), intent(in) :: array
    integer, dimension(size(array,dim=3),size(array,dim=2),size(array,dim=1)) :: reorder_int_123x321
    integer :: i1
    do i1 = 1, size(array,dim=1)
      reorder_int_123x321(:,:,i1) = transpose(array(i1,:,:))
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_wp_123x321(array)
    real(wp), dimension(:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=3),size(array,dim=2),size(array,dim=1)) :: reorder_wp_123x321
    integer :: i1
    do i1 = 1, size(array,dim=1)
      reorder_wp_123x321(:,:,i1) = transpose(array(i1,:,:))
    end do
  end function

  ! ----- reorder for 4D arrays -----
  ! reorder the indecies of 4D array
  function reorder_int_1234x1243(array)
    integer, dimension(:,:,:,:), intent(in) :: array
    integer, dimension(size(array,dim=1),size(array,dim=2),size(array,dim=4),size(array,dim=3)) :: reorder_int_1234x1243
    integer :: i4,i3
    do i3 = 1, size(array,dim=3)
      do i4 = 1, size(array,dim=4)
        reorder_int_1234x1243(:,:,i4,i3) = array(:,:,i3,i4)
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_wp_1234x1243(array)
    real(wp), dimension(:,:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=1),size(array,dim=2),size(array,dim=4),size(array,dim=3)) :: reorder_wp_1234x1243
    integer :: i4,i3
    do i3 = 1, size(array,dim=3)
      do i4 = 1, size(array,dim=4)
        reorder_wp_1234x1243(:,:,i4,i3) = array(:,:,i3,i4)
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_int_1234x1324(array)
    integer, dimension(:,:,:,:), intent(in) :: array
    integer, dimension(size(array,dim=1),size(array,dim=3),size(array,dim=2),size(array,dim=4)) :: reorder_int_1234x1324
    integer :: i2,i4
    do i4 = 1, size(array,dim=4)
      do i2 = 1, size(array,dim=2)
        reorder_int_1234x1324(:,:,i2,i4) = array(:,i2,:,i4)
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_wp_1234x1324(array)
    real(wp), dimension(:,:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=1),size(array,dim=3),size(array,dim=2),size(array,dim=4)) :: reorder_wp_1234x1324
    integer :: i2,i4
    do i4 = 1, size(array,dim=4)
      do i2 = 1, size(array,dim=2)
        reorder_wp_1234x1324(:,:,i2,i4) = array(:,i2,:,i4)
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_int_1234x1342(array)
    integer, dimension(:,:,:,:), intent(in) :: array
    integer, dimension(size(array,dim=1),size(array,dim=3),size(array,dim=4),size(array,dim=2)) :: reorder_int_1234x1342
    integer :: i4,i2
    do i2 = 1, size(array,dim=2)
      do i4 = 1, size(array,dim=4)
        reorder_int_1234x1342(:,:,i4,i2) = array(:,i2,:,i4)
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_wp_1234x1342(array)
    real(wp), dimension(:,:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=1),size(array,dim=3),size(array,dim=4),size(array,dim=2)) :: reorder_wp_1234x1342
    integer :: i4,i2
    do i2 = 1, size(array,dim=2)
      do i4 = 1, size(array,dim=4)
        reorder_wp_1234x1342(:,:,i4,i2) = array(:,i2,:,i4)
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_int_1234x1423(array)
    integer, dimension(:,:,:,:), intent(in) :: array
    integer, dimension(size(array,dim=1),size(array,dim=4),size(array,dim=2),size(array,dim=3)) :: reorder_int_1234x1423
    integer :: i2,i3
    do i3 = 1, size(array,dim=3)
      do i2 = 1, size(array,dim=2)
        reorder_int_1234x1423(:,:,i2,i3) = array(:,i2,i3,:)
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_wp_1234x1423(array)
    real(wp), dimension(:,:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=1),size(array,dim=4),size(array,dim=2),size(array,dim=3)) :: reorder_wp_1234x1423
    integer :: i2,i3
    do i3 = 1, size(array,dim=3)
      do i2 = 1, size(array,dim=2)
        reorder_wp_1234x1423(:,:,i2,i3) = array(:,i2,i3,:)
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_int_1234x1432(array)
    integer, dimension(:,:,:,:), intent(in) :: array
    integer, dimension(size(array,dim=1),size(array,dim=4),size(array,dim=3),size(array,dim=2)) :: reorder_int_1234x1432
    integer :: i3,i2
    do i2 = 1, size(array,dim=2)
      do i3 = 1, size(array,dim=3)
        reorder_int_1234x1432(:,:,i3,i2) = array(:,i2,i3,:)
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_wp_1234x1432(array)
    real(wp), dimension(:,:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=1),size(array,dim=4),size(array,dim=3),size(array,dim=2)) :: reorder_wp_1234x1432
    integer :: i3,i2
    do i2 = 1, size(array,dim=2)
      do i3 = 1, size(array,dim=3)
        reorder_wp_1234x1432(:,:,i3,i2) = array(:,i2,i3,:)
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_int_1234x2134(array)
    integer, dimension(:,:,:,:), intent(in) :: array
    integer, dimension(size(array,dim=2),size(array,dim=1),size(array,dim=3),size(array,dim=4)) :: reorder_int_1234x2134
    integer :: i3,i4
    do i4 = 1, size(array,dim=4)
      do i3 = 1, size(array,dim=3)
        reorder_int_1234x2134(:,:,i3,i4) = transpose(array(:,:,i3,i4))
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_wp_1234x2134(array)
    real(wp), dimension(:,:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=2),size(array,dim=1),size(array,dim=3),size(array,dim=4)) :: reorder_wp_1234x2134
    integer :: i3,i4
    do i4 = 1, size(array,dim=4)
      do i3 = 1, size(array,dim=3)
        reorder_wp_1234x2134(:,:,i3,i4) = transpose(array(:,:,i3,i4))
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_int_1234x2143(array)
    integer, dimension(:,:,:,:), intent(in) :: array
    integer, dimension(size(array,dim=2),size(array,dim=1),size(array,dim=4),size(array,dim=3)) :: reorder_int_1234x2143
    integer :: i4,i3
    do i3 = 1, size(array,dim=3)
      do i4 = 1, size(array,dim=4)
        reorder_int_1234x2143(:,:,i4,i3) = transpose(array(:,:,i3,i4))
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_wp_1234x2143(array)
    real(wp), dimension(:,:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=2),size(array,dim=1),size(array,dim=4),size(array,dim=3)) :: reorder_wp_1234x2143
    integer :: i4,i3
    do i3 = 1, size(array,dim=3)
      do i4 = 1, size(array,dim=4)
        reorder_wp_1234x2143(:,:,i4,i3) = transpose(array(:,:,i3,i4))
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_int_1234x2314(array)
    integer, dimension(:,:,:,:), intent(in) :: array
    integer, dimension(size(array,dim=2),size(array,dim=3),size(array,dim=1),size(array,dim=4)) :: reorder_int_1234x2314
    integer :: i1,i4
    do i4 = 1, size(array,dim=4)
      do i1 = 1, size(array,dim=1)
        reorder_int_1234x2314(:,:,i1,i4) = array(i1,:,:,i4)
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_wp_1234x2314(array)
    real(wp), dimension(:,:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=2),size(array,dim=3),size(array,dim=1),size(array,dim=4)) :: reorder_wp_1234x2314
    integer :: i1,i4
    do i4 = 1, size(array,dim=4)
      do i1 = 1, size(array,dim=1)
        reorder_wp_1234x2314(:,:,i1,i4) = array(i1,:,:,i4)
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_int_1234x2341(array)
    integer, dimension(:,:,:,:), intent(in) :: array
    integer, dimension(size(array,dim=2),size(array,dim=3),size(array,dim=4),size(array,dim=1)) :: reorder_int_1234x2341
    integer :: i4,i1
    do i1 = 1, size(array,dim=1)
      do i4 = 1, size(array,dim=4)
        reorder_int_1234x2341(:,:,i4,i1) = array(i1,:,:,i4)
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_wp_1234x2341(array)
    real(wp), dimension(:,:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=2),size(array,dim=3),size(array,dim=4),size(array,dim=1)) :: reorder_wp_1234x2341
    integer :: i4,i1
    do i1 = 1, size(array,dim=1)
      do i4 = 1, size(array,dim=4)
        reorder_wp_1234x2341(:,:,i4,i1) = array(i1,:,:,i4)
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_int_1234x2413(array)
    integer, dimension(:,:,:,:), intent(in) :: array
    integer, dimension(size(array,dim=2),size(array,dim=4),size(array,dim=1),size(array,dim=3)) :: reorder_int_1234x2413
    integer :: i1,i3
    do i3 = 1, size(array,dim=3)
      do i1 = 1, size(array,dim=1)
        reorder_int_1234x2413(:,:,i1,i3) = array(i1,:,i3,:)
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_wp_1234x2413(array)
    real(wp), dimension(:,:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=2),size(array,dim=4),size(array,dim=1),size(array,dim=3)) :: reorder_wp_1234x2413
    integer :: i1,i3
    do i3 = 1, size(array,dim=3)
      do i1 = 1, size(array,dim=1)
        reorder_wp_1234x2413(:,:,i1,i3) = array(i1,:,i3,:)
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_int_1234x2431(array)
    integer, dimension(:,:,:,:), intent(in) :: array
    integer, dimension(size(array,dim=2),size(array,dim=4),size(array,dim=3),size(array,dim=1)) :: reorder_int_1234x2431
    integer :: i3,i1
    do i1 = 1, size(array,dim=1)
      do i3 = 1, size(array,dim=3)
        reorder_int_1234x2431(:,:,i3,i1) = array(i1,:,i3,:)
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_wp_1234x2431(array)
    real(wp), dimension(:,:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=2),size(array,dim=4),size(array,dim=3),size(array,dim=1)) :: reorder_wp_1234x2431
    integer :: i3,i1
    do i1 = 1, size(array,dim=1)
      do i3 = 1, size(array,dim=3)
        reorder_wp_1234x2431(:,:,i3,i1) = array(i1,:,i3,:)
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_int_1234x3124(array)
    integer, dimension(:,:,:,:), intent(in) :: array
    integer, dimension(size(array,dim=3),size(array,dim=1),size(array,dim=2),size(array,dim=4)) :: reorder_int_1234x3124
    integer :: i2,i4
    do i4 = 1, size(array,dim=4)
      do i2 = 1, size(array,dim=2)
        reorder_int_1234x3124(:,:,i2,i4) = transpose(array(:,i2,:,i4))
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_wp_1234x3124(array)
    real(wp), dimension(:,:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=3),size(array,dim=1),size(array,dim=2),size(array,dim=4)) :: reorder_wp_1234x3124
    integer :: i2,i4
    do i4 = 1, size(array,dim=4)
      do i2 = 1, size(array,dim=2)
        reorder_wp_1234x3124(:,:,i2,i4) = transpose(array(:,i2,:,i4))
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_int_1234x3142(array)
    integer, dimension(:,:,:,:), intent(in) :: array
    integer, dimension(size(array,dim=3),size(array,dim=1),size(array,dim=4),size(array,dim=2)) :: reorder_int_1234x3142
    integer :: i4,i2
    do i2 = 1, size(array,dim=2)
      do i4 = 1, size(array,dim=4)
        reorder_int_1234x3142(:,:,i4,i2) = transpose(array(:,i2,:,i4))
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_wp_1234x3142(array)
    real(wp), dimension(:,:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=3),size(array,dim=1),size(array,dim=4),size(array,dim=2)) :: reorder_wp_1234x3142
    integer :: i4,i2
    do i2 = 1, size(array,dim=2)
      do i4 = 1, size(array,dim=4)
        reorder_wp_1234x3142(:,:,i4,i2) = transpose(array(:,i2,:,i4))
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_int_1234x3214(array)
    integer, dimension(:,:,:,:), intent(in) :: array
    integer, dimension(size(array,dim=3),size(array,dim=2),size(array,dim=1),size(array,dim=4)) :: reorder_int_1234x3214
    integer :: i1,i4
    do i4 = 1, size(array,dim=4)
      do i1 = 1, size(array,dim=1)
        reorder_int_1234x3214(:,:,i1,i4) = transpose(array(i1,:,:,i4))
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_wp_1234x3214(array)
    real(wp), dimension(:,:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=3),size(array,dim=2),size(array,dim=1),size(array,dim=4)) :: reorder_wp_1234x3214
    integer :: i1,i4
    do i4 = 1, size(array,dim=4)
      do i1 = 1, size(array,dim=1)
        reorder_wp_1234x3214(:,:,i1,i4) = transpose(array(i1,:,:,i4))
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_int_1234x3241(array)
    integer, dimension(:,:,:,:), intent(in) :: array
    integer, dimension(size(array,dim=3),size(array,dim=2),size(array,dim=4),size(array,dim=1)) :: reorder_int_1234x3241
    integer :: i4,i1
    do i1 = 1, size(array,dim=1)
      do i4 = 1, size(array,dim=4)
        reorder_int_1234x3241(:,:,i4,i1) = transpose(array(i1,:,:,i4))
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_wp_1234x3241(array)
    real(wp), dimension(:,:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=3),size(array,dim=2),size(array,dim=4),size(array,dim=1)) :: reorder_wp_1234x3241
    integer :: i4,i1
    do i1 = 1, size(array,dim=1)
      do i4 = 1, size(array,dim=4)
        reorder_wp_1234x3241(:,:,i4,i1) = transpose(array(i1,:,:,i4))
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_int_1234x3412(array)
    integer, dimension(:,:,:,:), intent(in) :: array
    integer, dimension(size(array,dim=3),size(array,dim=4),size(array,dim=1),size(array,dim=2)) :: reorder_int_1234x3412
    integer :: i1,i2
    do i2 = 1, size(array,dim=2)
      do i1 = 1, size(array,dim=1)
        reorder_int_1234x3412(:,:,i1,i2) = array(i1,i2,:,:)
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_wp_1234x3412(array)
    real(wp), dimension(:,:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=3),size(array,dim=4),size(array,dim=1),size(array,dim=2)) :: reorder_wp_1234x3412
    integer :: i1,i2
    do i2 = 1, size(array,dim=2)
      do i1 = 1, size(array,dim=1)
        reorder_wp_1234x3412(:,:,i1,i2) = array(i1,i2,:,:)
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_int_1234x3421(array)
    integer, dimension(:,:,:,:), intent(in) :: array
    integer, dimension(size(array,dim=3),size(array,dim=4),size(array,dim=2),size(array,dim=1)) :: reorder_int_1234x3421
    integer :: i2,i1
    do i1 = 1, size(array,dim=1)
      do i2 = 1, size(array,dim=2)
        reorder_int_1234x3421(:,:,i2,i1) = array(i1,i2,:,:)
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_wp_1234x3421(array)
    real(wp), dimension(:,:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=3),size(array,dim=4),size(array,dim=2),size(array,dim=1)) :: reorder_wp_1234x3421
    integer :: i2,i1
    do i1 = 1, size(array,dim=1)
      do i2 = 1, size(array,dim=2)
        reorder_wp_1234x3421(:,:,i2,i1) = array(i1,i2,:,:)
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_int_1234x4123(array)
    integer, dimension(:,:,:,:), intent(in) :: array
    integer, dimension(size(array,dim=4),size(array,dim=1),size(array,dim=2),size(array,dim=3)) :: reorder_int_1234x4123
    integer :: i2,i3
    do i3 = 1, size(array,dim=3)
      do i2 = 1, size(array,dim=2)
        reorder_int_1234x4123(:,:,i2,i3) = transpose(array(:,i2,i3,:))
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_wp_1234x4123(array)
    real(wp), dimension(:,:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=4),size(array,dim=1),size(array,dim=2),size(array,dim=3)) :: reorder_wp_1234x4123
    integer :: i2,i3
    do i3 = 1, size(array,dim=3)
      do i2 = 1, size(array,dim=2)
        reorder_wp_1234x4123(:,:,i2,i3) = transpose(array(:,i2,i3,:))
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_int_1234x4132(array)
    integer, dimension(:,:,:,:), intent(in) :: array
    integer, dimension(size(array,dim=4),size(array,dim=1),size(array,dim=3),size(array,dim=2)) :: reorder_int_1234x4132
    integer :: i3,i2
    do i2 = 1, size(array,dim=2)
      do i3 = 1, size(array,dim=3)
        reorder_int_1234x4132(:,:,i3,i2) = transpose(array(:,i2,i3,:))
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_wp_1234x4132(array)
    real(wp), dimension(:,:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=4),size(array,dim=1),size(array,dim=3),size(array,dim=2)) :: reorder_wp_1234x4132
    integer :: i3,i2
    do i2 = 1, size(array,dim=2)
      do i3 = 1, size(array,dim=3)
        reorder_wp_1234x4132(:,:,i3,i2) = transpose(array(:,i2,i3,:))
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_int_1234x4213(array)
    integer, dimension(:,:,:,:), intent(in) :: array
    integer, dimension(size(array,dim=4),size(array,dim=2),size(array,dim=1),size(array,dim=3)) :: reorder_int_1234x4213
    integer :: i1,i3
    do i3 = 1, size(array,dim=3)
      do i1 = 1, size(array,dim=1)
        reorder_int_1234x4213(:,:,i1,i3) = transpose(array(i1,:,i3,:))
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_wp_1234x4213(array)
    real(wp), dimension(:,:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=4),size(array,dim=2),size(array,dim=1),size(array,dim=3)) :: reorder_wp_1234x4213
    integer :: i1,i3
    do i3 = 1, size(array,dim=3)
      do i1 = 1, size(array,dim=1)
        reorder_wp_1234x4213(:,:,i1,i3) = transpose(array(i1,:,i3,:))
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_int_1234x4231(array)
    integer, dimension(:,:,:,:), intent(in) :: array
    integer, dimension(size(array,dim=4),size(array,dim=2),size(array,dim=3),size(array,dim=1)) :: reorder_int_1234x4231
    integer :: i3,i1
    do i1 = 1, size(array,dim=1)
      do i3 = 1, size(array,dim=3)
        reorder_int_1234x4231(:,:,i3,i1) = transpose(array(i1,:,i3,:))
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_wp_1234x4231(array)
    real(wp), dimension(:,:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=4),size(array,dim=2),size(array,dim=3),size(array,dim=1)) :: reorder_wp_1234x4231
    integer :: i3,i1
    do i1 = 1, size(array,dim=1)
      do i3 = 1, size(array,dim=3)
        reorder_wp_1234x4231(:,:,i3,i1) = transpose(array(i1,:,i3,:))
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_int_1234x4312(array)
    integer, dimension(:,:,:,:), intent(in) :: array
    integer, dimension(size(array,dim=4),size(array,dim=3),size(array,dim=1),size(array,dim=2)) :: reorder_int_1234x4312
    integer :: i1,i2
    do i2 = 1, size(array,dim=2)
      do i1 = 1, size(array,dim=1)
        reorder_int_1234x4312(:,:,i1,i2) = transpose(array(i1,i2,:,:))
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_wp_1234x4312(array)
    real(wp), dimension(:,:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=4),size(array,dim=3),size(array,dim=1),size(array,dim=2)) :: reorder_wp_1234x4312
    integer :: i1,i2
    do i2 = 1, size(array,dim=2)
      do i1 = 1, size(array,dim=1)
        reorder_wp_1234x4312(:,:,i1,i2) = transpose(array(i1,i2,:,:))
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_int_1234x4321(array)
    integer, dimension(:,:,:,:), intent(in) :: array
    integer, dimension(size(array,dim=4),size(array,dim=3),size(array,dim=2),size(array,dim=1)) :: reorder_int_1234x4321
    integer :: i2,i1
    do i1 = 1, size(array,dim=1)
      do i2 = 1, size(array,dim=2)
        reorder_int_1234x4321(:,:,i2,i1) = transpose(array(i1,i2,:,:))
      end do
    end do
  end function

  ! reorder the indecies of 4D array
  function reorder_wp_1234x4321(array)
    real(wp), dimension(:,:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=4),size(array,dim=3),size(array,dim=2),size(array,dim=1)) :: reorder_wp_1234x4321
    integer :: i2,i1
    do i1 = 1, size(array,dim=1)
      do i2 = 1, size(array,dim=2)
        reorder_wp_1234x4321(:,:,i2,i1) = transpose(array(i1,i2,:,:))
      end do
    end do
  end function

end module

