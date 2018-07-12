subroutine stop_on_err(msg)
  !
  ! Print error message and stop
  !
  use iso_fortran_env, only : error_unit
  character(len=*), intent(in) :: msg
  if(len_trim(msg) > 0) then
    write (error_unit,*) trim(msg)
    write (error_unit,*) "test_lw_solver stopping"
    stop
  end if
end subroutine
! ----------------------------------------------------------------------------------
program test_lw_solver
  use mo_rte_kind,         only: wp
  use mo_optical_props,    only: ty_optical_props_arry, &
                                 ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
  use mo_source_functions, only: ty_source_func_lw
  use mo_rte_lw,           only: rte_lw

  use mo_fluxes_bygpoint,  only: ty_fluxes_bygpoint
  use mo_test_files_io,    only: read_optical_prop_values, read_lw_Planck_sources, read_direction, read_lw_bc, read_lw_rt, &
                                 write_gpt_fluxes
  implicit none
  ! ----------------------------------------------------------------------------------
  integer :: ncol, nlay, ngpt
  integer :: nang
  integer :: b, nBlocks, colS, colE
  integer, parameter :: blockSize = 8

  character(len=128) :: fileName = 'rrtmgp-inputs-outputs.nc'

  class(ty_optical_props_arry), allocatable :: atmos_full,   atmos_block
  type (ty_source_func_lw)                  :: sources_full, sources_block

  real(wp), dimension(:,  :), allocatable :: sfc_emis
  real(wp), dimension(:    ), allocatable :: t_sfc
  real(wp), dimension(:,:,:), allocatable, target :: flux_up, flux_dn

  logical :: top_at_1
  type(ty_fluxes_bygpoint) :: fluxes
  ! ----------------------------------------------------------------------------------
  !
  ! In early implementations this called the LW solver at an intermediate level in the
  !   call tree - after some error checking had been done but before deciding e.g.
  !   whether to use no-scattering or two-stream solvers.
  ! The current implementation is functionally the same as compute_fluxes_from_optics but
  !   writes out g-point fluxes
  !

  call read_optical_prop_values (fileName, atmos_full)
  call read_lw_Planck_sources   (fileName, sources_full)
  call read_direction    (fileName, top_at_1)
  call read_lw_bc        (fileName, t_sfc, sfc_emis)
  call read_lw_rt        (fileName, nang)
  if(nang > 1) print *, "  Doing ", nang, "-angle integration"
  ncol = sources_full%get_ncol()
  nlay = sources_full%get_nlay()
  ngpt = sources_full%get_ngpt()


  allocate(flux_up(ncol,nlay+1,ngpt), flux_dn(ncol,nlay+1,ngpt))
  flux_dn(:,MERGE(1, nlay+1, top_at_1),:) = 0._wp

  select type (atmos_full)
    class is (ty_optical_props_1scl)
      allocate(ty_optical_props_1scl::atmos_block)
    class is (ty_optical_props_2str)
      allocate(ty_optical_props_2str::atmos_block)
    class is (ty_optical_props_nstr)
      allocate(ty_optical_props_nstr::atmos_block)
  end select

  !
  ! Loop over subsets of the problem
  !
  nBlocks = ncol/blockSize ! Integer division
  print *, "Doing ", nBlocks, "blocks of size ", blockSize
    do b = 1, nBlocks
      colS = (b-1) * blockSize + 1
      colE =  b    * blockSize

      call stop_on_err(  atmos_full%get_subset(colS, blockSize, atmos_block))
      call stop_on_err(sources_full%get_subset(colS, blockSize, sources_block))
      fluxes%gpt_flux_up => flux_up(colS:colE,:,:)
      fluxes%gpt_flux_dn => flux_dn(colS:colE,:,:)
      call stop_on_err(rte_lw(atmos_block,           &
                              top_at_1,              &
                              sources_block,         &
                              sfc_emis(:,colS:colE), &
                              fluxes, n_gauss_angles = nang))
    end do

    if(mod(ncol, blockSize) /= 0) then
      colS = ncol/blockSize * blockSize + 1  ! Integer arithmetic
      colE = ncol
      print *, "Doing ", colE-colS+1, "extra columns"

      call stop_on_err(  atmos_full%get_subset(colS, colE-colS+1, atmos_block))
      call stop_on_err(sources_full%get_subset(colS, colE-colS+1, sources_block))
      fluxes%gpt_flux_up => flux_up(colS:colE,:,:)
      fluxes%gpt_flux_dn => flux_dn(colS:colE,:,:)
      call stop_on_err(rte_lw(atmos_block,           &
                              top_at_1,              &
                              sources_block,         &
                              sfc_emis(:,colS:colE), &
                              fluxes, n_gauss_angles = nang))
  endif

  call write_gpt_fluxes(fileName, flux_up, flux_dn)
end program test_lw_solver
