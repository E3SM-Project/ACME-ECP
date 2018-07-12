
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
  use mo_rte_kind,           only: wp
  use mo_gas_optics,         only: ty_gas_optics_specification
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_optical_props,      only: ty_optical_props, &
                                   ty_optical_props_arry, ty_optical_props_1scl, ty_optical_props_2str
  use mo_source_functions,   only: ty_source_func_lw
  use mo_fluxes_byband,      only: ty_fluxes_byband
  ! ---- RRTMPG driver
  use mo_rte_lw,             only: rte_lw
  use mo_rte_sw,             only: rte_sw
  use mo_heating_rates, only: compute_heating_rate

  ! ---- I/O for test format files.
  use mo_test_files_io,     only: read_atmos, is_lw, is_sw, &
                                  read_lw_bc, read_sw_bc, read_lw_rt,  &
                                  write_fluxes, write_dir_fluxes, write_heating_rates, &
                                  write_spectral_disc
  use mo_load_coefficients, only: load_and_init

  implicit none
  real(wp), parameter :: pi = acos(-1._wp)
  ! ----------------------------------------------------------------------------------

  ! Arrays: dimensions (col, lay)
  real(wp), dimension(:,:),   allocatable :: p_lay, t_lay, p_lev
  real(wp), dimension(:,:), target, &
                              allocatable :: col_dry
  ! Longwave only
  real(wp), dimension(:,:),   target, &
                              allocatable :: t_lev
  real(wp), dimension(:),     allocatable :: t_sfc
  real(wp), dimension(:,:),   allocatable :: emis_sfc ! First dimension is band

  ! Shortwave only
  real(wp), dimension(:),     allocatable :: sza, tsi, mu0
  real(wp), dimension(:,:),   allocatable :: sfc_alb_dir, sfc_alb_dif ! First dimension is band
  real(wp)                                :: tsi_scaling = -999._wp

  ! Source functions
  !   Longwave
  type(ty_source_func_lw)               :: lw_sources
  !   Shortwave
  real(wp), dimension(:,:), allocatable :: toa_flux

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
  type(ty_gas_optics_specification) :: k_dist
  type(ty_gas_concs)                :: gas_concs, gas_concs_subset
  class(ty_optical_props_arry), target, &
                        allocatable :: optical_props
  type(ty_fluxes_byband),       target  &
                                    :: fluxes

  !
  ! Inputs to RRTMGP
  !
  logical :: top_at_1

  integer :: ncol, nlay, nbnd, ngpt
  integer :: b, nBlocks, colS, colE, nSubcols, nang
  integer, parameter :: blockSize = 4

  character(len=256) :: k_dist_file, input_file
  ! ----------------------------------------------------------------------------------
  !
  ! k-distribution file and input-output files must be paired: LW or SW
  !
  k_dist_file = 'coefficients.nc'
  input_file = "rrtmgp-flux-inputs-outputs.nc"

  !
  ! Read temperature, pressure, gas concentrations, then variables specific
  !  to LW or SW problems. Arrays are allocated as they are read
  !
  call read_atmos(input_file,                         &
                   p_lay, t_lay, p_lev, t_lev,        &
                   gas_concs, col_dry)

  if(is_lw(input_file)) then
    call read_lw_bc(input_file, t_sfc, emis_sfc)
    ! Number of quadrature angles
    call read_lw_rt(input_file, nang)
  else
    call read_sw_bc(input_file, sza, tsi, tsi_scaling, sfc_alb_dir, sfc_alb_dif)
    allocate(mu0(size(sza)))
    mu0 = cos(sza * pi/180.)
  end if

  call load_and_init(k_dist, k_dist_file, gas_concs)
  if(k_dist%is_internal_source_present() .neqv. is_lw(input_file)) &
    stop "flux_compute: gas optics and input-output file disagree about SW/LW"

  !
  ! Problem sizes; allocation of output arrays for full problem
  !
  ncol = size(p_lay, 1)
  nlay = size(p_lay, 2)
  nbnd = k_dist%get_nband()
  ngpt = k_dist%get_ngpt()
  top_at_1 = p_lay(1, 1) < p_lay(1, nlay)
  allocate(    flux_up (ncol,nlay+1     ),     flux_dn(ncol,nlay+1     ), &
               flux_net(ncol,nlay+1     ))
  allocate(bnd_flux_up (ncol,nlay+1,nbnd), bnd_flux_dn(ncol,nlay+1,nbnd), &
           bnd_flux_net(ncol,nlay+1,nbnd))
  allocate(heating_rate(ncol,nlay), bnd_heating_rate(ncol,nlay,nbnd))
  if(is_sw(input_file)) &
    allocate(flux_dir(ncol,nlay+1), bnd_flux_dir(ncol,nlay+1,nbnd))

  if(is_sw(input_file)) then
    allocate(ty_optical_props_2str::optical_props)
  else
    allocate(ty_optical_props_1scl::optical_props)
  end if
  call stop_on_err(optical_props%init(k_dist))

  !
  ! Loop over subsets of the problem
  !
  nBlocks = ncol/blockSize ! Integer division
  print *, "Doing ", nBlocks, "blocks of size ", blockSize
  !
  ! Allocation of intermediate arrays
  !
  if(is_sw(input_file)) then
    select type(optical_props)
      class is (ty_optical_props_2str)
        call stop_on_err(optical_props%alloc_2str(blockSize, nlay))
      class default
        call stop_on_err("flux_compute: Trying to do SW calculation but haven't allocated _2str")
    end select

    allocate(toa_flux(blockSize, ngpt))
  else
    select type(optical_props)
      class is (ty_optical_props_1scl)
        call stop_on_err(optical_props%alloc_1scl(blockSize, nlay))
      class default
        call stop_on_err("flux_compute: Trying to do LW calculation but haven't allocated _1scl")
    end select
    call stop_on_err(lw_sources%alloc(k_dist, blockSize, nlay))
  end if

  do b = 1, nBlocks
    print *, "block ", b
    colS = (b-1) * blockSize + 1
    colE =  b    * blockSize
    nSubcols = colE-colS+1
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
    call stop_on_err(gas_concs%get_subset(colS, nSubcols, gas_concs_subset))
    if(is_sw(input_file)) then
      !
      ! Gas optics, including source functions
      !
      if(allocated(col_dry)) then
        call stop_on_err(k_dist%gas_optics(p_lay(colS:colE,:), &
                                           p_lev(colS:colE,:), &
                                           t_lay(colS:colE,:), &
                                           gas_concs_subset,   &
                                           optical_props,      &
                                           toa_flux,           &
                                           col_dry = col_dry(colS:colE,:)))
      else
        call stop_on_err(k_dist%gas_optics(p_lay(colS:colE,:), &
                                           p_lev(colS:colE,:), &
                                           t_lay(colS:colE,:), &
                                           gas_concs_subset,   &
                                           optical_props,      &
                                           toa_flux))
      end if
      if(tsi_scaling > 0.0_wp) toa_flux(:,:) =  toa_flux(:,:) * tsi_scaling
      !
      ! Radiative transfer
      !
      call stop_on_err(rte_sw(optical_props,               &
                                 top_at_1,                 &
                                 mu0(colS:colE),           &
                                 toa_flux,                 &
                                 sfc_alb_dir(:,colS:colE), &
                                 sfc_alb_dif(:,colS:colE), &
                                 fluxes))
    else
      !
      ! Gas optics, including source functions
      !
      if(allocated(col_dry)) then
        call stop_on_err(k_dist%gas_optics(p_lay(colS:colE,:), &
                                           p_lev(colS:colE,:), &
                                           t_lay(colS:colE,:), &
                                           t_sfc(colS:colE  ), &
                                           gas_concs_subset,   &
                                           optical_props,      &
                                           lw_sources,         &
                                           tlev    = t_lev  (colS:colE,:), &
                                           col_dry = col_dry(colS:colE,:)))
      else
        call stop_on_err(k_dist%gas_optics(p_lay(colS:colE,:), &
                                           p_lev(colS:colE,:), &
                                           t_lay(colS:colE,:), &
                                           t_sfc(colS:colE  ), &
                                           gas_concs_subset,   &
                                           optical_props,      &
                                           lw_sources,         &
                                           tlev    = t_lev  (colS:colE,:)))
      end if
      !
      ! Radiative transfer
      !
      call stop_on_err(rte_lw(optical_props,            &
                                 top_at_1,              &
                                 lw_sources,            &
                                 emis_sfc(:,colS:colE), &
                                 fluxes, n_gauss_angles = nang))
    end if
  end do
  if(mod(ncol, blockSize) /= 0) then
    colS = ncol/blockSize * blockSize + 1  ! Integer arithmetic
    colE = ncol
    print *, "Doing ", colE-colS+1, "extra columns"
    nSubcols = colE-colS+1
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
    call stop_on_err(gas_concs%get_subset(colS, nSubcols, gas_concs_subset))

    if(is_sw(input_file)) then
      select type(optical_props)
        class is (ty_optical_props_2str)
          call stop_on_err(optical_props%alloc_2str(nSubcols, nlay))
        class default
          call stop_on_err("flux_compute: Trying to do SW calculation but haven't allocated _2str")
      end select
      if(allocated(toa_flux)) deallocate(toa_flux)
      allocate(toa_flux(nSubcols, ngpt))
    else
      select type(optical_props)
        class is (ty_optical_props_1scl)
          call stop_on_err(optical_props%alloc_1scl(nSubcols, nlay))
        class default
          call stop_on_err("flux_compute: Trying to do LW calculation but haven't allocated _1scl")
      end select
      call stop_on_err(lw_sources%alloc(nSubcols, nlay))
    end if


    if(is_sw(input_file)) then
      !
      ! Gas optics, including source functions
      !
      if(allocated(col_dry)) then
        call stop_on_err(k_dist%gas_optics(p_lay(colS:colE,:), &
                                           p_lev(colS:colE,:), &
                                           t_lay(colS:colE,:), &
                                           gas_concs_subset,   &
                                           optical_props,      &
                                           toa_flux,           &
                                           col_dry = col_dry(colS:colE,:)))
      else
        call stop_on_err(k_dist%gas_optics(p_lay(colS:colE,:), &
                                           p_lev(colS:colE,:), &
                                           t_lay(colS:colE,:), &
                                           gas_concs_subset,   &
                                           optical_props,      &
                                           toa_flux))
      end if
      if(tsi_scaling > 0.0_wp) toa_flux(:,:) =  toa_flux(:,:) * tsi_scaling
      !
      ! Radiative transfer
      !
      call stop_on_err(rte_sw(optical_props,            &
                                 top_at_1,                 &
                                 mu0(colS:colE),           &
                                 toa_flux,                 &
                                 sfc_alb_dir(:,colS:colE), &
                                 sfc_alb_dif(:,colS:colE), &
                                 fluxes))
    else
      !
      ! Gas optics, including source functions
      !
      if(allocated(col_dry)) then
        call stop_on_err(k_dist%gas_optics(p_lay(colS:colE,:), &
                                           p_lev(colS:colE,:), &
                                           t_lay(colS:colE,:), &
                                           t_sfc(colS:colE  ), &
                                           gas_concs_subset,   &
                                           optical_props,      &
                                           lw_sources,         &
                                           tlev    = t_lev  (colS:colE,:), &
                                           col_dry = col_dry(colS:colE,:)))
      else
        call stop_on_err(k_dist%gas_optics(p_lay(colS:colE,:), &
                                           p_lev(colS:colE,:), &
                                           t_lay(colS:colE,:), &
                                           t_sfc(colS:colE  ), &
                                           gas_concs_subset,   &
                                           optical_props,      &
                                           lw_sources,         &
                                           tlev    = t_lev  (colS:colE,:)))
      end if
      !
      ! Radiative transfer
      !
      call stop_on_err(rte_lw(optical_props,         &
                                 top_at_1,              &
                                 lw_sources,         &
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
  call write_spectral_disc(input_file, optical_props)
  call write_fluxes(input_file, flux_up, flux_dn, flux_net, bnd_flux_up, bnd_flux_dn, bnd_flux_net)
  call write_heating_rates(input_file, heating_rate, bnd_heating_rate)
  if(k_dist%is_external_source_present()) &
    call write_dir_fluxes(input_file, flux_dir, bnd_flux_dir)

end program flux_compute
