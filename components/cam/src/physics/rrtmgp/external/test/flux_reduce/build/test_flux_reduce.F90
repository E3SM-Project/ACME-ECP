subroutine stop_on_err(msg)
  !
  ! Print error message and stop
  !
  use iso_fortran_env, only : error_unit
  character(len=*), intent(in) :: msg
  if(len_trim(msg) > 0) then
    write (error_unit,*) trim(msg)
    write (error_unit,*) "test_flux_reduce stopping"
    stop
  end if
end subroutine
! ----------------------------------------------------------------------------------
program test_flux_reduce
  use mo_rte_kind,   only: wp
  use mo_spectral_disc, only: ty_spectral_disc
  use mo_fluxes_byband, only: ty_fluxes_byband

  use mo_test_files_io, only: is_sw, read_gpt_fluxes, read_direction, read_spectral_disc, &
                              write_fluxes, write_dir_fluxes
  implicit none
  ! ----------------------------------------------------------------------------------
  character(len=128) :: fileName    = 'rrtmgp-inputs-outputs.nc'

  type(ty_spectral_disc) :: spectral_disc
  type(ty_fluxes_byband) :: fluxes

  real(wp), dimension(:,  :), allocatable :: toa_src

  real(wp), dimension(:,:,:), allocatable         :: gpt_flux_up, gpt_flux_dn, gpt_flux_dir
  real(wp), dimension(:,:  ), allocatable, target ::     flux_up,     flux_dn,     flux_dir,     flux_net
  real(wp), dimension(:,:,:), allocatable, target :: bnd_flux_up, bnd_flux_dn, bnd_flux_dir, bnd_flux_net

  integer :: ncol, nlev, nbnd, ngpt
  logical :: doing_sw, top_at_1
  ! ----------------------------------------------------------------------------------

  call read_direction(fileName, top_at_1)
  call read_spectral_disc(fileName, spectral_disc)
  doing_sw = is_sw(fileName)
  if(doing_sw) then
    call read_gpt_fluxes(fileName, gpt_flux_up, gpt_flux_dn, gpt_flux_dir)
  else
    call read_gpt_fluxes(fileName, gpt_flux_up, gpt_flux_dn)
  end if

  ncol = size(gpt_flux_up, 1)
  nlev = size(gpt_flux_up, 2)
  nbnd = spectral_disc%get_nband()
  ngpt = spectral_disc%get_ngpt()

  allocate(flux_up    (ncol,nlev     ), flux_dn    (ncol,nlev     ),     flux_net(ncol,nlev))
  allocate(bnd_flux_up(ncol,nlev,nbnd), bnd_flux_dn(ncol,nlev,nbnd), bnd_flux_net(ncol,nlev,nbnd))
  fluxes%ty_fluxes_broadband%flux_up =>     flux_up
  fluxes%ty_fluxes_broadband%flux_dn =>     flux_dn
  fluxes%ty_fluxes_broadband%flux_net=>     flux_net
  fluxes%bnd_flux_up                 => bnd_flux_up
  fluxes%bnd_flux_dn                 => bnd_flux_dn
  fluxes%bnd_flux_net                => bnd_flux_net
  if(doing_sw) then
    allocate(flux_dir(ncol,nlev), bnd_flux_dir(ncol,nlev,nbnd))
    fluxes%ty_fluxes_broadband%flux_dn_dir =>     flux_dir
    fluxes%bnd_flux_dn_dir                 => bnd_flux_dir
  end if

  if(doing_sw) then
    call stop_on_err(fluxes%reduce(gpt_flux_up, gpt_flux_dn, spectral_disc, top_at_1, gpt_flux_dir))
  else
    call stop_on_err(fluxes%reduce(gpt_flux_up, gpt_flux_dn, spectral_disc, top_at_1))
  end if

  call write_fluxes(fileName, flux_up, flux_dn, flux_net, bnd_flux_up, bnd_flux_dn, bnd_flux_net)
  if(doing_sw) call write_dir_fluxes(fileName, flux_dir, bnd_flux_dir)

end program test_flux_reduce
