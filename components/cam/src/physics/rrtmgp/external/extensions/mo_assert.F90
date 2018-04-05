module mo_assert

   use mo_rte_kind,   only: wp
   use, intrinsic :: ieee_arithmetic, only: &
      isnan => ieee_is_nan, &
      ieee_negative_inf, &
      ieee_positive_inf, &
      ieee_class, &
      operator(==)

   implicit none

   interface assert_valid
      module procedure assert_valid_1d, assert_valid_2d, assert_valid_3d
   end interface

contains

elemental logical function isinf(x)
   real(wp), intent(in) :: x

   isinf = (ieee_positive_inf == ieee_class(x) .or. &
            ieee_negative_inf == ieee_class(x))

end function isinf


function assert_valid_1d(x, message) result(error_message)

   real(wp), intent(in) :: x(:)
   character(len=*), intent(in), optional :: message
   character(len=128) :: error_message

   error_message = ''

   if (any(isnan(x))) then
      error_message = 'assert_valid failed (NaN)'
   else if (any(isinf(x))) then
      error_message = 'assert_valid failed (inf)'
   end if

   if (present(message)) then
      if (len(trim(error_message)) > 0) then
         error_message = trim(error_message) // ': ' // message
      end if
   end if

end function


function assert_valid_2d(x, message) result(error_message)

   real(wp), intent(in) :: x(:,:)
   character(len=*), intent(in), optional :: message
   character(len=128) :: error_message

   error_message = ''

   if (any(isnan(x))) then
      error_message = 'assert_valid failed (NaN)'
   else if (any(isinf(x))) then
      error_message = 'assert_valid failed (inf)'
   end if

   if (present(message)) then
      if (len(trim(error_message)) > 0) then
         error_message = trim(error_message) // ': ' // message
      end if
   end if

end function assert_valid_2d


function assert_valid_3d(x, message) result(error_message)

   real(wp), intent(in) :: x(:,:,:)
   character(len=*), intent(in), optional :: message
   character(len=128) :: error_message

   error_message = ''

   if (any(isnan(x))) then
      error_message = 'assert_valid failed (NaN)'
   else if (any(isinf(x))) then
      error_message = 'assert_valid failed (inf)'
   end if

   if (present(message)) then
      if (len(trim(error_message)) > 0) then
         error_message = trim(error_message) // ': ' // message
      end if
   end if

end function  assert_valid_3d


end module mo_assert
