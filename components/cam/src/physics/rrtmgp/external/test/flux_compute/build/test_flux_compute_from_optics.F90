
subroutine stop_on_err(error_msg)
  use iso_fortran_env, only : error_unit
  character(len=*), intent(in) :: error_msg

  if(error_msg /= "") then
    write (error_unit,*) trim(error_msg)
    write (error_unit,*) "test_flux_compute stopping"
    stop
  end if

end subroutine stop_on_err
!-----------------------------
program flux_compute
  use mo_rte_kind,        only: wp
  use mo_gas_concentrations,       &
                        only: ty_gas_concs
  use mo_optical_props, only: ty_optical_props, ty_optical_props_arry, &
                              ty_optical_props_1scl, ty_optical_props_2str
  use mo_source_functions, &
                        only: ty_source_func_lw
  use mo_fluxes,        only: ty_fluxes
  use mo_fluxes_byband, only: ty_fluxes_byband
  ! ---- RRTMPG driver
  use mo_rte_lw,     only: rte_lw
  use mo_rte_sw,     only: rte_sw
  use mo_heating_rates, only: compute_heating_rate

  ! ---- I/O for test format files.
  use mo_test_files_io,     only: read_atmos, is_lw, is_sw, &
                                  read_lw_bc, read_sw_bc,   &
                                  read_lw_Planck_sources, read_sw_solar_sources, &
                                  read_direction, read_optical_prop_values,      &
                                  read_lw_rt,  &
                                  write_fluxes, write_dir_fluxes, write_heating_rates

  implicit none
  real, parameter :: pi = acos(-1._wp)
  ! ----------------------------------------------------------------------------------
  ! You only need to pressures to compute heating rates, but still ...
  !
  real(wp), dimension(:,:),   allocatable :: p_lay, t_lay, p_lev, t_lev, col_dry
  type(ty_gas_concs)                      :: gas_concs
  ! LW variables
  real(wp), dimension(:,:),   allocatable :: emis_sfc ! First dimension is band
  real(wp), dimension(:),     allocatable :: t_sfc    ! Not needed but included in read() call
  ! Shortwave only
  real(wp), dimension(:),     allocatable :: sza, tsi, mu0
  real(wp), dimension(:,:),   allocatable :: sfc_alb_dir, sfc_alb_dif ! First dimension is band
  real(wp)                                :: tsi_scaling = -999._wp

  ! Source functions
  !   Shortwave
  real(wp), dimension(:,:),    allocatable :: toa_flux


  real(wp), dimension(:,:  ), target, &
                               allocatable ::     flux_up,      flux_dn, &
                                                  flux_net,     flux_dir
  real(wp), dimension(:,:,:), target, &
                               allocatable :: bnd_flux_up,  bnd_flux_dn, &
                                              bnd_flux_net, bnd_flux_dir
  real(wp), dimension(:,:),    allocatable :: heating_rate
  real(wp), dimension(:,:,:),  allocatable :: bnd_heating_rate


  !
  ! Derived types for interacting with RRTMGP
  !
  class(ty_optical_props_arry), allocatable :: atmos, atmos_subset
  type(ty_source_func_lw)                   :: sources, sources_subset
  type(ty_fluxes_byband)                    :: fluxes

  !
  ! Inputs to RRTMGP - arrays, but only extent 1
  !
  logical                                   :: top_at_1

  integer :: ncol, nlay, nbnd, ngpt
  integer :: b, nBlocks, colS, colE, nSubcols, nang
  integer, parameter :: blockSize = 4

  character(len=256) :: input_file
  ! ----------------------------------------------------------------------------------
  !
  ! k-distribution file and input-output files must be paired: LW or SW
  !
  input_file = "rrtmgp-inputs-outputs.nc"
  call read_optical_prop_values(input_file, atmos)
  call read_direction    (input_file, top_at_1)
  call read_atmos(input_file,                 &
                  p_lay, t_lay, p_lev, t_lev, &
                  gas_concs, col_dry)
  deallocate(t_lay, p_lay, t_lev, col_dry)
  !
  ! Problem sizes; allocation of output arrays for full problem
  !
  ncol = atmos%get_ncol()
  nlay = atmos%get_nlay()
  nbnd = atmos%get_nband()
  ngpt = atmos%get_ngpt()

  allocate(atmos_subset, source = atmos)

  !
  ! Outputs for the full problem
  !
  allocate(    flux_up (ncol,nlay+1     ),     flux_dn(ncol,nlay+1     ), &
               flux_net(ncol,nlay+1     ))
  allocate(bnd_flux_up (ncol,nlay+1,nbnd), bnd_flux_dn(ncol,nlay+1,nbnd), &
           bnd_flux_net(ncol,nlay+1,nbnd))
  allocate(heating_rate(ncol,nlay), bnd_heating_rate(ncol,nlay,nbnd))
  if(is_sw(input_file)) &
    allocate(flux_dir(ncol,nlay+1), bnd_flux_dir(ncol,nlay+1,nbnd))

  !
  ! Additional boundary conditions
  !
  if(is_lw(input_file)) then
    call read_lw_bc     (input_file, t_sfc, emis_sfc)
    call read_lw_Planck_sources(input_file, sources)
    ! Number of quadrature angles
    call read_lw_rt(input_file, nang)
  else
    call read_sw_bc     (input_file, sza, tsi, tsi_scaling, sfc_alb_dir, sfc_alb_dif)
    call read_sw_solar_sources(input_file, toa_flux)
    allocate(mu0(size(sza)))
    mu0 = cos(sza * pi/180.)
  end if

  !
  ! Loop over subsets of the problem
  !
  nBlocks = ncol/blockSize ! Integer division
  print *, "Doing ", nBlocks, "blocks of size ", blockSize

  do b = 1, nBlocks
    colS = (b-1) * blockSize + 1
    colE =  b    * blockSize
    nSubcols = colE-colS+1
    call stop_on_err(atmos%get_subset  (colS, nSubcols, atmos_subset))
    if(is_lw(input_file)) call stop_on_err(sources%get_subset(colS, nSubcols, sources_subset))
    fluxes%flux_up      => flux_up(colS:colE,:)
    fluxes%flux_dn      => flux_dn(colS:colE,:)
    fluxes%flux_net     => flux_net(colS:colE,:)
    fluxes%bnd_flux_up  => bnd_flux_up(colS:colE,:,:)
    fluxes%bnd_flux_dn  => bnd_flux_dn(colS:colE,:,:)
    fluxes%bnd_flux_net => bnd_flux_net(colS:colE,:,:)
    if(is_sw(input_file)) then
      fluxes%flux_dn_dir     => flux_dir(colS:colE,:)
      fluxes%bnd_flux_dn_dir => bnd_flux_dir(colS:colE,:,:)
    end if
    if(is_sw(input_file)) then
      if(tsi_scaling > 0.0_wp) toa_flux(:,:) =  toa_flux(:,:) * tsi_scaling
      call stop_on_err(rte_sw(atmos_subset,             &
                                 top_at_1,                 &
                                 mu0(colS:colE),           &
                                 toa_flux(colS:colE, :),   &
                                 sfc_alb_dir(:,colS:colE), &
                                 sfc_alb_dif(:,colS:colE), &
                                 fluxes))
    else
      call stop_on_err(rte_lw(atmos_subset,             &
                                 top_at_1,              &
                                 sources_subset,        &
                                 emis_sfc(:,colS:colE), &
                                 fluxes, n_gauss_angles = nang))
    end if
  end do
  if(mod(ncol, blockSize) /= 0) then
    colS = ncol/blockSize * blockSize + 1  ! Integer arithmetic
    colE = ncol
    print *, "Doing ", colE-colS+1, "extra columns"
    nSubcols = colE-colS+1
    call stop_on_err(atmos%get_subset  (colS, nSubcols, atmos_subset))
    if(is_lw(input_file)) call stop_on_err(sources%get_subset(colS, nSubcols, sources_subset))
    fluxes%flux_up      => flux_up(colS:colE,:)
    fluxes%flux_dn      => flux_dn(colS:colE,:)
    fluxes%flux_net     => flux_net(colS:colE,:)
    fluxes%bnd_flux_up  => bnd_flux_up(colS:colE,:,:)
    fluxes%bnd_flux_dn  => bnd_flux_dn(colS:colE,:,:)
    fluxes%bnd_flux_net => bnd_flux_net(colS:colE,:,:)
    if(is_sw(input_file)) then
      fluxes%flux_dn_dir     => flux_dir(colS:colE,:)
      fluxes%bnd_flux_dn_dir => bnd_flux_dir(colS:colE,:,:)
    end if
    if(is_sw(input_file)) then
      if(tsi_scaling > 0.0_wp) toa_flux(:,:) =  toa_flux(:,:) * tsi_scaling
      call stop_on_err(rte_sw(atmos_subset,             &
                                 top_at_1,                 &
                                 mu0(colS:colE),           &
                                 toa_flux(colS:colE, :),   &
                                 sfc_alb_dir(:,colS:colE), &
                                 sfc_alb_dif(:,colS:colE), &
                                 fluxes))
    else
      call stop_on_err(rte_lw(atmos_subset,             &
                                 top_at_1,              &
                                 sources_subset,        &
                                 emis_sfc(:,colS:colE), &
                                 fluxes, n_gauss_angles = nang))
    end if
  end if

  !
  ! Heating rates
  !
  call stop_on_err(compute_heating_rate(flux_up, flux_dn, p_lev, heating_rate))
  do b = 1, nbnd
      call stop_on_err(compute_heating_rate(bnd_flux_up(:,:,b), bnd_flux_dn(:,:,b), p_lev, bnd_heating_rate(:,:,b)))
  end do

  !
  ! ... and write everything out
  !
  call write_fluxes(input_file, flux_up, flux_dn, flux_net, bnd_flux_up, bnd_flux_dn, bnd_flux_net)
  call write_heating_rates(input_file, heating_rate, bnd_heating_rate)
  if(is_sw(input_file)) &
    call write_dir_fluxes(input_file, flux_dir, bnd_flux_dir)

end program flux_compute
