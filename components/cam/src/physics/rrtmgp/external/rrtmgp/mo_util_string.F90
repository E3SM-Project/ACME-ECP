!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! @brief Module to provide string related utilities
!!
!! @author Andre Wehe (2014)
!!
!

module mo_util_string
  implicit none
  private
  public :: lower_case, string_in_array, string_loc_in_array

  ! List of character for case conversion
  character(len=26), parameter :: LOWER_CASE_CHARS = 'abcdefghijklmnopqrstuvwxyz'
  character(len=26), parameter :: UPPER_CASE_CHARS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

contains
  ! --------------------------------------------------------------------------------------
  pure function lower_case( input_string ) result( output_string )
    character(len=*), intent(in) :: input_string
    character(len=len(input_string)) :: output_string
    integer :: i, n

    ! Copy input string
    output_string = input_string

    ! Convert case character by character
    do i = 1, len(output_string)
      n = index(UPPER_CASE_CHARS, output_string(i:i))
      if ( n /= 0 ) output_string(i:i) = LOWER_CASE_CHARS(n:n)
    end do
  end function
  ! --------------------------------------------------------------------------------------
  !
  ! Is string somewhere in array? 
  !
  pure function string_in_array(string, array) 
    character(len=*),               intent(in) :: string
    character(len=*), dimension(:), intent(in) :: array
    logical                                    :: string_in_array
    
    integer :: i
    character(len=len_trim(string)) :: lc_string  
    
    string_in_array = .false. 
    lc_string = lower_case(trim(string))
    do i = 1, size(array) 
      if(lc_string == lower_case(trim(array(i)))) then 
        string_in_array = .true. 
        exit 
      end if 
    end do 
  end function string_in_array
  ! --------------------------------------------------------------------------------------
  !
  ! Is string somewhere in array? 
  !
  pure function string_loc_in_array(string, array) 
    character(len=*),               intent(in) :: string
    character(len=*), dimension(:), intent(in) :: array
    integer                                    :: string_loc_in_array
    
    integer :: i
    character(len=len_trim(string)) :: lc_string  
    
    string_loc_in_array = -1
    lc_string = lower_case(trim(string))
    do i = 1, size(array) 
      if(lc_string == lower_case(trim(array(i)))) then 
        string_loc_in_array = i
        exit 
      end if 
    end do 
  end function string_loc_in_array
end module
