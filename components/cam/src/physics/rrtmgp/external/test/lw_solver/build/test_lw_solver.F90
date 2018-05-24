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
  use mo_rte_kind,   only: wp
  use mo_optical_props, only: ty_optical_props_arry, &
                              ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
  use mo_spectral_disc, only: ty_spectral_disc
  use mo_lw_solver,     only: lw_solver_init, lw_solver

  use mo_test_files_io, only: read_optical_props, read_lw_Planck_sources, read_direction, read_lw_bc, read_lw_rt, &
                              read_spectral_disc, write_gpt_fluxes
  implicit none
  ! ----------------------------------------------------------------------------------
  integer :: ncol, nlay, ngpt
  integer :: nang
  integer :: b, nBlocks, colS, colE
  integer, parameter :: blockSize = 8

  character(len=128) :: fileName = 'rrtmgp-inputs-outputs.nc'

  class(ty_optical_props_arry), allocatable :: atmos_full, atmos_block

  real(wp), dimension(:,:,:), allocatable :: lay_source, lev_source_inc, lev_source_dec
  real(wp), dimension(:,  :), allocatable :: sfc_emis, sfc_source
  real(wp), dimension(:    ), allocatable :: t_sfc
  real(wp), dimension(:,:,:), allocatable :: flux_up, flux_dn
  real(wp), dimension(:,  :), allocatable :: sfc_emis_gpt

  integer :: i, ibnd, igpt
  logical :: top_at_1
  type(ty_spectral_disc) :: spectral_disc
  ! ----------------------------------------------------------------------------------

  call read_optical_props(fileName, atmos_full)
  call read_lw_Planck_sources   (fileName, lay_source, lev_source_inc, lev_source_dec, sfc_source)
  call read_direction    (fileName, top_at_1)
  call read_lw_bc        (fileName, t_sfc, sfc_emis)
  call read_lw_rt        (fileName, nang)
  call stop_on_err(lw_solver_init(n_angles=nang))
  if(nang > 1) print *, "  Doing ", nang, "-angle integration"
  ncol = size(lay_source, 1)
  nlay = size(lay_source, 2)
  ngpt = size(lay_source, 3)

   ! Read in gpt2band
  call read_spectral_disc(fileName, spectral_disc)

   ! Set the surface emissivity for each g-point, depending on the band
  allocate(sfc_emis_gpt(ncol,ngpt))
  do igpt = 1, ngpt
    ibnd = spectral_disc%convert_gpt2band(igpt)  ! Get band number for this g-point
    sfc_emis_gpt(1:ncol,igpt) = sfc_emis(ibnd,1:ncol)
  end do

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

      call stop_on_err(atmos_full%get_subset(colS, blockSize, atmos_block))
      call stop_on_err(lw_solver(blockSize, nlay, ngpt, top_at_1, &
                                 atmos_block,                     &
                                 lay_source(colS:colE,:,:),       &
                                 lev_source_inc(colS:colE,:,:),   &
                                 lev_source_dec(colS:colE,:,:),   &
                                 sfc_emis_gpt(colS:colE,  :), &
                                 sfc_source(colS:colE,  :),   &
                                 flux_up(colS:colE,:,:), flux_dn(colS:colE,:,:)) )
    end do

    if(mod(ncol, blockSize) /= 0) then
      colS = ncol/blockSize * blockSize + 1  ! Integer arithmetic
      colE = ncol
      print *, "Doing ", colE-colS+1, "extra columns"

      call stop_on_err(atmos_full%get_subset(colS, colE-colS+1, atmos_block))
      call stop_on_err(lw_solver(colE-colS+1, nlay, ngpt, top_at_1, &
                                 atmos_block,                     &
                                 lay_source(colS:colE,:,:),       &
                                 lev_source_inc(colS:colE,:,:),   &
                                 lev_source_dec(colS:colE,:,:),   &
                                 sfc_emis_gpt(colS:colE,  :), &
                                 sfc_source(colS:colE,  :), &
                                 flux_up(colS:colE,:,:), flux_dn(colS:colE,:,:)) )
  endif

  call write_gpt_fluxes(fileName, flux_up, flux_dn)
end program test_lw_solver
