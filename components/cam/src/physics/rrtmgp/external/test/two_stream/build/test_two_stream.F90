subroutine stop_on_err(msg)
  !
  ! Print error message and stop
  !
  use iso_fortran_env, only : error_unit
  character(len=*), intent(in) :: msg
  if(len_trim(msg) > 0) then
    write (error_unit,*) trim(msg)
    write (error_unit,*) "test_sw_two_stream stopping"
    stop
  end if
end subroutine
! ----------------------------------------------------------------------------------
program test_two_stream
  use mo_rte_kind,   only: wp
  use mo_optical_props, only: ty_optical_props_arry, &
                              ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
  use mo_rte_solver_kernels, &
                        only: sw_two_stream, lw_two_stream

  use mo_test_files_io, only: read_optical_prop_values, is_sw, is_lw, read_sw_bc, write_two_stream
  implicit none
  ! ----------------------------------------------------------------------------------
  integer :: ncol, nlay, ngpt
  integer :: b, nBlocks, colS, colE
  integer, parameter :: blockSize = 8

  character(len=128) :: fileName = 'rrtmgp-inputs-outputs.nc'

  class(ty_optical_props_arry), allocatable :: atmos

  real(wp), dimension(:    ), allocatable :: mu0, tsi
  real(wp), dimension(:,:  ), allocatable :: gamma1, gamma2
  real(wp), dimension(:,  :), allocatable :: sfc_alb_dir, sfc_alb_dif
  real(wp), dimension(:,:,:), allocatable :: Rdif, Tdif, Rdir, Tdir, Tnoscat
  real(wp)                                :: tsi_scaling = -999._wp
  integer :: i, j, k, ngpts_per_band
  logical :: do_sw
  ! ----------------------------------------------------------------------------------

  call read_optical_prop_values(fileName, atmos)
  do_sw = is_sw(fileName)
  if(do_sw) then
    call read_sw_bc (fileName, mu0, tsi, tsi_scaling, sfc_alb_dir, sfc_alb_dif)
    mu0 = cos(mu0 * acos(-1._wp)/180.)
    ! Some variables aren't needed for this problem
    deallocate(tsi, sfc_alb_dir, sfc_alb_dif)
  end if

  ncol = atmos%get_ncol()
  nlay = atmos%get_nlay()
  ngpt = atmos%get_ngpt()

  allocate(Rdif(ncol,nlay,ngpt), Tdif(ncol,nlay,ngpt))
  if(do_sw) then
    allocate(Rdir(ncol,nlay,ngpt), Tdir(ncol,nlay,ngpt), Tnoscat(ncol,nlay,ngpt))
  else
    allocate(gamma1(blockSize,nlay), gamma2(blockSize,nlay))
  end if

  select type(atmos)
    type is (ty_optical_props_2str)
      !
      ! Loop over subsets of the problem
      !
      nBlocks = ncol/blockSize ! Integer division
      print *, "Doing ", nBlocks, "blocks of size ", blockSize
      do b = 1, nBlocks
        colS = (b-1) * blockSize + 1
        colE =  b    * blockSize
        if(do_sw) then
          do k = 1, ngpt
            call sw_two_stream(colE-colS+1, nlay, &
                               mu0(colS:colE),           atmos%tau(colS:colE,:,k), &
                               atmos%ssa(colS:colE,:,k), atmos%g  (colS:colE,:,k), &
                               Rdif(colS:colE,:,k), Tdif(colS:colE,:,k),           &
                               Rdir(colS:colE,:,k), Tdir(colS:colE,:,k), Tnoscat(colS:colE,:,k))
          end do
        else
          do k = 1, ngpt
            call lw_two_stream(colE-colS+1, nlay, &
                               atmos%tau(colS:colE,:,k), &
                               atmos%ssa(colS:colE,:,k), atmos%g  (colS:colE,:,k), &
                               gamma1, gamma2,           &
                               Rdif(colS:colE,:,k), Tdif(colS:colE,:,k))
          end do
        end if
      end do

      if(mod(ncol, blockSize) /= 0) then
        colS = ncol/blockSize * blockSize + 1  ! Integer arithmetic
        colE = ncol
        print *, "Doing ", colE-colS+1, "extra columns"
        if(do_sw) then
          do k = 1, ngpt
            call sw_two_stream(colE-colS+1, nlay, &
                               mu0(colS:colE),           atmos%tau(colS:colE,:,k), &
                               atmos%ssa(colS:colE,:,k), atmos%g  (colS:colE,:,k), &
                               Rdif(colS:colE,:,k), Tdif(colS:colE,:,k),           &
                               Rdir(colS:colE,:,k), Tdir(colS:colE,:,k), Tnoscat(colS:colE,:,k))
          end do
        else
          do k = 1, ngpt
            call lw_two_stream(colE-colS+1, nlay, &
                               atmos%tau(colS:colE,:,k), &
                               atmos%ssa(colS:colE,:,k), atmos%g  (colS:colE,:,k), &
                               gamma1, gamma2,           &
                               Rdif(colS:colE,:,k), Tdif(colS:colE,:,k))
          end do
        end if
      end if

      if(do_sw) then
        call write_two_stream(fileName, Rdif, Tdif, Rdir, Tdir, Tnoscat)
      else
        call write_two_stream(fileName, Rdif, Tdif)
      end if
    type is (ty_optical_props_nstr)
      call stop_on_err("Haven't implemented two-stream calculations for multi-stream inputs")
    type is (ty_optical_props_1scl)
      call stop_on_err("No ssa, g, information in file")
  end select
end program test_two_stream
