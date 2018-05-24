!-----------------------------------------------------------------------
!
! title: test for rte_lw
!
! description:
!   Test optical properties operators
!
!-----------------------------------------------------------------------

program test
  use mo_rte_kind, only: wp
  use mo_optical_props
  use mo_test_files_io, only: write_optical_props
  implicit none

  integer, parameter :: ncol = 15, nlay = 5, ngpt = 14, nmom = 10, nbnd = 3
  integer, dimension(2,nbnd), parameter :: band2gpt = reshape((/ 1, 4, 5, 10, 11, 14 /), (/2,nbnd/))
  
  character(len = 128) :: error_msg

  call test_increment_1scalar_by_1scalar
  call test_increment_1scalar_by_2stream
  call test_increment_1scalar_by_nstream
  call test_increment_2stream_by_1scalar
  call test_increment_2stream_by_2stream
  call test_increment_2stream_by_nstream
  call test_increment_nstream_by_1scalar
  call test_increment_nstream_by_2stream
  call test_increment_nstream_by_nstream

  call test_expand_and_increment_1scalar_by_1scalar
  call test_expand_and_increment_1scalar_by_2stream
  call test_expand_and_increment_1scalar_by_nstream
  call test_expand_and_increment_2stream_by_1scalar
  call test_expand_and_increment_2stream_by_2stream
  call test_expand_and_increment_2stream_by_nstream
  call test_expand_and_increment_nstream_by_1scalar
  call test_expand_and_increment_nstream_by_2stream
  call test_expand_and_increment_nstream_by_nstream

  ! all done
  print *, 'test end.'

contains
  ! -----
  subroutine stop_on_err(msg)
    !
    ! Print error message and stop  
    ! 
    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: msg
    if(len_trim(msg) > 0) then 
      write (error_unit,*) trim(msg)
      write (error_unit,*) "test_optical_props stopping"
      stop
    end if 
  end subroutine
  ! -----
  subroutine prep_file(fileName, ncol, nlay)
    ! 
    ! Create a file with col and lay dimensions of the correct size 
    !
    use netcdf 
    character(len=*), intent(in) :: fileName
    integer,          intent(in) :: ncol, nlay
    
    integer :: ncid, dimid 
    
    if (nf90_create(trim(fileName), NF90_CLOBBER, ncid) /= nf90_noerr) & 
      call stop_on_err("Can't create file " // trim(fileName))
    if (nf90_def_dim(ncid, 'col', ncol, dimid) /= nf90_noerr) & 
      call stop_on_err("Can't create dimension col in " // trim(fileName))
    if (nf90_def_dim(ncid, 'lay', nlay, dimid) /= nf90_noerr) & 
      call stop_on_err("Can't create dimension lay in " // trim(fileName))
    dimid = nf90_close(ncid)   
  end subroutine prep_file
  ! -----

  subroutine test_increment_1scalar_by_1scalar
    type(ty_optical_props_1scl) :: a
    type(ty_optical_props_1scl) :: b
    character(len=128) :: fileName = 'out_test_increment_1scalar_by_1scalar.nc'
    call stop_on_err(a%init_1scl(ncol, nlay, ngpt))
    call stop_on_err(b%init_1scl(ncol, nlay, ngpt))
    a%tau = 1._wp
    b%tau = 2._wp
    call stop_on_err(a%increment_by(b))
    call prep_file(fileName, ncol, nlay)
    call write_optical_props(fileName,a)
  end subroutine

  subroutine test_increment_1scalar_by_2stream
    type(ty_optical_props_1scl) :: a
    type(ty_optical_props_2str) :: b
    character(len=128) :: fileName = 'out_test_increment_1scalar_by_2stream.nc'
    call stop_on_err(a%init_1scl(ncol, nlay, ngpt))
    call stop_on_err(b%init_2str(ncol, nlay, ngpt))
    a%tau = 1._wp
    b%tau = 2._wp
    b%ssa = 0.7_wp
    b%g = 0.8_wp
    call stop_on_err(a%increment_by(b))
    call prep_file(fileName, ncol, nlay)
    call write_optical_props(fileName,a)
  end subroutine

  subroutine test_increment_1scalar_by_nstream
    type(ty_optical_props_1scl) :: a
    type(ty_optical_props_nstr) :: b
    character(len=128) :: fileName = 'out_test_increment_1scalar_by_nstream.nc'
    call stop_on_err(a%init_1scl(ncol, nlay, ngpt))
    call stop_on_err(b%init_nstr(nmom, ncol, nlay, ngpt))
    a%tau = 1._wp
    b%tau = 2._wp
    b%ssa = 0.7_wp
    b%p = 0.8_wp
    call stop_on_err(a%increment_by(b))
    call prep_file(fileName, ncol, nlay)
    call write_optical_props(fileName,a)
  end subroutine

  subroutine test_increment_2stream_by_1scalar
    type(ty_optical_props_2str) :: a
    type(ty_optical_props_1scl) :: b
    character(len=128) :: fileName = 'out_test_increment_2stream_by_1scalar.nc'
    call stop_on_err(a%init_2str(ncol, nlay, ngpt))
    call stop_on_err(b%init_1scl(ncol, nlay, ngpt))
    a%tau = 2._wp
    a%ssa = 0.7_wp
    a%g = 0.8_wp
    b%tau = 3._wp
    call stop_on_err(a%increment_by(b))
    call prep_file(fileName, ncol, nlay)
    call write_optical_props(fileName,a)
  end subroutine

  subroutine test_increment_2stream_by_2stream
    type(ty_optical_props_2str) :: a
    type(ty_optical_props_2str) :: b
    character(len=128) :: fileName = 'out_test_increment_2stream_by_2stream.nc'
    call stop_on_err(a%init_2str(ncol, nlay, ngpt))
    call stop_on_err(b%init_2str(ncol, nlay, ngpt))
    a%tau = 2._wp
    a%ssa = 0.7_wp
    a%g = 0.8_wp
    b%tau = 3._wp
    b%ssa = 0.6_wp
    b%g = 0.5_wp
    call stop_on_err(a%increment_by(b))
    call prep_file(fileName, ncol, nlay)
    call write_optical_props(fileName,a)
  end subroutine

  subroutine test_increment_2stream_by_nstream
    type(ty_optical_props_2str) :: a
    type(ty_optical_props_nstr) :: b
    character(len=128) :: fileName = 'out_test_increment_2stream_by_nstream.nc'
    call stop_on_err(a%init_2str(ncol, nlay, ngpt))
    call stop_on_err(b%init_nstr(nmom, ncol, nlay, ngpt))
    a%tau = 2._wp
    a%ssa = 0.4_wp
    a%g = 0.7_wp
    b%tau = 3._wp
    b%ssa = 0.7_wp
    b%p = 0.8_wp
    call stop_on_err(a%increment_by(b))
    call prep_file(fileName, ncol, nlay)
    call write_optical_props(fileName,a)
  end subroutine

  subroutine test_increment_nstream_by_1scalar
    type(ty_optical_props_nstr) :: a
    type(ty_optical_props_1scl) :: b
    character(len=128) :: fileName = 'out_test_increment_nstream_by_1scalar.nc'
    call stop_on_err(a%init_nstr(nmom, ncol, nlay, ngpt))
    call stop_on_err(b%init_1scl(ncol, nlay, ngpt))
    a%tau = 2._wp
    a%ssa = 0.7_wp
    a%p = 0.8_wp
    b%tau = 3._wp
    call stop_on_err(a%increment_by(b))
    call prep_file(fileName, ncol, nlay)
    call write_optical_props(fileName,a)
  end subroutine

  subroutine test_increment_nstream_by_2stream
    type(ty_optical_props_nstr) :: a
    type(ty_optical_props_2str) :: b
    character(len=128) :: fileName = 'out_test_increment_nstream_by_2stream.nc'
    call stop_on_err(a%init_nstr(nmom, ncol, nlay, ngpt))
    call stop_on_err(b%init_2str(ncol, nlay, ngpt))
    a%tau = 3._wp
    a%ssa = 0.7_wp
    a%p = 0.8_wp
    b%tau = 2._wp
    b%ssa = 0.4_wp
    b%g = 0.6_wp
    call stop_on_err(a%increment_by(b))
    call prep_file(fileName, ncol, nlay)
    call write_optical_props(fileName,a)
  end subroutine

  subroutine test_increment_nstream_by_nstream
    type(ty_optical_props_nstr) :: a
    type(ty_optical_props_nstr) :: b
    character(len=128) :: fileName = 'out_test_increment_nstream_by_nstream.nc'
    call stop_on_err(a%init_nstr(nmom, ncol, nlay, ngpt))
    call stop_on_err(b%init_nstr(nmom, ncol, nlay, ngpt))
    a%tau = 3_wp
    a%ssa = 0.7_wp
    a%p = 0.7_wp
    a%p(1,:,:,:) = 0.8_wp
    b%tau = 2._wp
    b%ssa = 0.4_wp
    b%p = 0.6_wp
    b%p(1,:,:,:) = 0.5_wp
    call stop_on_err(a%increment_by(b))
    call prep_file(fileName, ncol, nlay)
    call write_optical_props(fileName,a)
  end subroutine

  ! ---- expand and increment ---------------------------------------------------------------------

  subroutine test_expand_and_increment_1scalar_by_1scalar
    type(ty_optical_props_1scl) :: a
    type(ty_optical_props_1scl) :: b
    character(len=128) :: fileName = 'out_test_expand_and_increment_1scalar_by_1scalar.nc'
    call stop_on_err(a%init_1scl(ncol, nlay, ngpt))
    call stop_on_err(b%init_1scl(ncol, nlay, nbnd))
    a%tau = 1._wp
    b%tau = 2._wp
    call stop_on_err(a%increment_by(b, band2gpt))
    call prep_file(fileName, ncol, nlay)
    call write_optical_props(fileName,a)
  end subroutine

  subroutine test_expand_and_increment_1scalar_by_2stream
    type(ty_optical_props_1scl) :: a
    type(ty_optical_props_2str) :: b
    character(len=128) :: fileName = 'out_test_expand_and_increment_1scalar_by_2stream.nc'
    call stop_on_err(a%init_1scl(ncol, nlay, ngpt))
    call stop_on_err(b%init_2str(ncol, nlay, nbnd))
    a%tau = 1._wp
    b%tau = 2._wp
    b%ssa = 0.7_wp
    b%g = 0.8_wp
    call stop_on_err(a%increment_by(b, band2gpt))
    call prep_file(fileName, ncol, nlay)
    call write_optical_props(fileName,a)
  end subroutine

  subroutine test_expand_and_increment_1scalar_by_nstream
    type(ty_optical_props_1scl) :: a
    type(ty_optical_props_nstr) :: b
    character(len=128) :: fileName = 'out_test_expand_and_increment_1scalar_by_nstream.nc'
    call stop_on_err(a%init_1scl(ncol, nlay, ngpt))
    call stop_on_err(b%init_nstr(nmom, ncol, nlay, nbnd))
    a%tau = 1._wp
    b%tau = 2._wp
    b%ssa = 0.7_wp
    b%p = 0.8_wp
    call stop_on_err(a%increment_by(b, band2gpt))
    call prep_file(fileName, ncol, nlay)
    call write_optical_props(fileName,a)
  end subroutine

  subroutine test_expand_and_increment_2stream_by_1scalar
    type(ty_optical_props_2str) :: a
    type(ty_optical_props_1scl) :: b
    character(len=128) :: fileName = 'out_test_expand_and_increment_2stream_by_1scalar.nc'
    call stop_on_err(a%init_2str(ncol, nlay, ngpt))
    call stop_on_err(b%init_1scl(ncol, nlay, nbnd))
    a%tau = 2._wp
    a%ssa = 0.7_wp
    a%g = 0.8_wp
    b%tau = 3._wp
    call stop_on_err(a%increment_by(b, band2gpt))
    call prep_file(fileName, ncol, nlay)
    call write_optical_props(fileName,a)
  end subroutine

  subroutine test_expand_and_increment_2stream_by_2stream
    type(ty_optical_props_2str) :: a
    type(ty_optical_props_2str) :: b
    character(len=128) :: fileName = 'out_test_expand_and_increment_2stream_by_2stream.nc'
    call stop_on_err(a%init_2str(ncol, nlay, ngpt))
    call stop_on_err(b%init_2str(ncol, nlay, nbnd))
    a%tau = 2._wp
    a%ssa = 0.7_wp
    a%g = 0.8_wp
    b%tau = 3._wp
    b%ssa = 0.6_wp
    b%g = 0.5_wp
    call stop_on_err(a%increment_by(b, band2gpt))
    call prep_file(fileName, ncol, nlay)
    call write_optical_props(fileName,a)
  end subroutine

  subroutine test_expand_and_increment_2stream_by_nstream
    type(ty_optical_props_2str) :: a
    type(ty_optical_props_nstr) :: b
    character(len=128) :: fileName = 'out_test_expand_and_increment_2stream_by_nstream.nc'
    call stop_on_err(a%init_2str(ncol, nlay, ngpt))
    call stop_on_err(b%init_nstr(nmom, ncol, nlay, nbnd))
    a%tau = 2._wp
    a%ssa = 0.4_wp
    a%g = 0.7_wp
    b%tau = 3._wp
    b%ssa = 0.7_wp
    b%p = 0.8_wp
    call stop_on_err(a%increment_by(b, band2gpt))
    call prep_file(fileName, ncol, nlay)
    call write_optical_props(fileName,a)
  end subroutine

  subroutine test_expand_and_increment_nstream_by_1scalar
    type(ty_optical_props_nstr) :: a
    type(ty_optical_props_1scl) :: b
    character(len=128) :: fileName = 'out_test_expand_and_increment_nstream_by_1scalar.nc'
    call stop_on_err(a%init_nstr(nmom, ncol, nlay, ngpt))
    call stop_on_err(b%init_1scl(ncol, nlay, nbnd))
    a%tau = 2._wp
    a%ssa = 0.7_wp
    a%p = 0.8_wp
    b%tau = 3._wp
    call stop_on_err(a%increment_by(b, band2gpt))
    call prep_file(fileName, ncol, nlay)
    call write_optical_props(fileName,a)
  end subroutine

  subroutine test_expand_and_increment_nstream_by_2stream
    type(ty_optical_props_nstr) :: a
    type(ty_optical_props_2str) :: b
    character(len=128) :: fileName = 'out_test_expand_and_increment_nstream_by_2stream.nc'
    call stop_on_err(a%init_nstr(nmom, ncol, nlay, ngpt))
    call stop_on_err(b%init_2str( ncol, nlay, nbnd))
    a%tau = 3._wp
    a%ssa = 0.7_wp
    a%p = 0.8_wp
    b%tau = 2._wp
    b%ssa = 0.4_wp
    b%g = 0.6_wp
    call stop_on_err(a%increment_by(b, band2gpt))
    call prep_file(fileName, ncol, nlay)
    call write_optical_props(fileName,a)
  end subroutine

  subroutine test_expand_and_increment_nstream_by_nstream
    type(ty_optical_props_nstr) :: a
    type(ty_optical_props_nstr) :: b
    character(len=128) :: fileName = 'out_test_expand_and_increment_nstream_by_nstream.nc'
    call stop_on_err(a%init_nstr(nmom, ncol, nlay, ngpt))
    call stop_on_err(b%init_nstr(nmom, ncol, nlay, nbnd))
    a%tau = 3_wp
    a%ssa = 0.7_wp
    a%p = 0.7_wp
    a%p(1,:,:,:) = 0.8_wp
    b%tau = 2._wp
    b%ssa = 0.4_wp
    b%p = 0.6_wp
    b%p(1,:,:,:) = 0.5_wp
    call stop_on_err(a%increment_by(b, band2gpt))
    call prep_file(fileName, ncol, nlay)
    call write_optical_props(fileName,a)
  end subroutine

end program
