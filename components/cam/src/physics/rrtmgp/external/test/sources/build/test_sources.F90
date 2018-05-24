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
  use mo_rte_kind,      only: wp, wl
  use mo_optical_props, only: ty_optical_props_arry, &
                              ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
  use mo_spectral_disc, only: ty_spectral_disc
  use mo_rte_solver_kernels, &
                        only: sw_two_stream, sw_source_2str, &
                              lw_two_stream, lw_source_2str, lw_combine_sources, lw_source_noscat, &
                              apply_BC

  use mo_test_files_io, only: read_optical_props, read_spectral_disc, &
                              is_sw, is_lw, read_direction, &
                              read_sw_bc, read_sw_solar_sources,  &
                              read_lw_bc, read_lw_Planck_sources, &
                              write_sources
  implicit none
  ! ----------------------------------------------------------------------------------
  integer :: ncol, nlay, ngpt

  character(len=128) :: fileName = 'rrtmgp-inputs-outputs.nc'

  class(ty_optical_props_arry), allocatable :: atmos

  ! SW-specific
  real(wp), dimension(:    ), allocatable :: mu0, tsi
  real(wp), dimension(:,  :), allocatable :: sfc_alb_dir, sfc_alb_dif, toa_src
  real(wp), dimension(:,  :), allocatable :: Rdir, Tdir, Tnoscat
  real(wp), dimension(:,:,:), allocatable :: flux_dn_dir
  real(wp)                                :: tsi_scaling
  ! LW-specific
  real(wp), dimension(:,:  ), allocatable :: gamma1, gamma2
  real(wp), dimension(:,:  ), allocatable :: sfc_emis, sfc_src, emis_sfc_bnd
  real(wp), dimension(:    ), allocatable :: t_sfc
  real(wp), dimension(:,:,:), allocatable, &
                              target      :: lay_src, lev_src_inc, lev_src_dec
  real(wp), dimension(:,:,:), pointer     :: lev_src_up, lev_src_dn
  ! Generic
  real(wp), dimension(:,  :), allocatable :: Rdif, Tdif
  real(wp), dimension(:,:,:), allocatable :: source_dn, source_up
  real(wp), dimension(:,  :), allocatable :: source_sfc
  real(wp), dimension(:,  :), allocatable :: lev_source

  type(ty_spectral_disc) :: spectral_disc
  integer :: i, j, k
  logical :: do_sw, do_lw, top_is_1
  ! ----------------------------------------------------------------------------------

  call read_optical_props(fileName, atmos)
  ncol = atmos%get_ncol()
  nlay = atmos%get_nlay()
  ngpt = atmos%get_ngpt()
  allocate(source_up(ncol, nlay, ngpt), source_dn(ncol, nlay, ngpt), source_sfc(ncol, ngpt))
  allocate(Rdif(ncol,nlay), Tdif(ncol,nlay))

  call read_spectral_disc(fileName, spectral_disc)
  call read_direction    (fileName, top_is_1)
  do_sw = is_sw(fileName); do_lw = .not. do_sw
  if(do_sw) then
    allocate(flux_dn_dir(ncol, nlay+1, ngpt))
    allocate(Rdir(ncol,nlay), Tdir(ncol,nlay), Tnoscat(ncol,nlay))
    call read_sw_bc (fileName, mu0, tsi, tsi_scaling, sfc_alb_dir, sfc_alb_dif)
    call read_sw_solar_sources(fileName, toa_src)
    mu0 = cos(mu0 * acos(-1._wp)/180.)
    call apply_BC(ncol, nlay, ngpt, logical(top_is_1, wl), toa_src, mu0, flux_dn_dir)
  else
    allocate(gamma1(ncol,nlay), gamma2(ncol,nlay), sfc_emis(ncol, ngpt), lev_source(ncol,nlay+1))
    call read_lw_Planck_sources(fileName, lay_src, lev_src_inc, lev_src_dec, sfc_src)
    call read_lw_bc            (fileName, t_sfc, emis_sfc_bnd)
    do k = 1, ngpt
      sfc_emis(1:ncol, k) = emis_sfc_bnd(spectral_disc%convert_gpt2band(k), :)
    end do
  end if

  select type (atmos)
    type is (ty_optical_props_1scl)
      if(do_lw) then
        if(top_is_1) then
          lev_src_up => lev_src_dec
          lev_src_dn => lev_src_inc
        else
          lev_src_up => lev_src_inc
          lev_src_dn => lev_src_dec
        end if
        do k = 1, ngpt
          call lw_source_noscat(ncol, nlay, &
                                lay_src(:,:,k), lev_src_up(:,:,k), lev_src_dn(:,:,k),     &
                                ! Diffusivity angle would be applied in solver
                                atmos%tau(:,:,k)*1.66_wp, exp(-atmos%tau(:,:,k)*1.66_wp), &
                                source_dn(:,:,k), source_up(:,:,k))
          ! Copied from kernels
          source_sfc(:, k) = sfc_emis(:,k) * sfc_src(:,k)
        end do
      else
        call stop_on_err("Haven't implemented direct beam source for SW")
      end if

    type is (ty_optical_props_2str)
      if(do_lw) then
        do k = 1, ngpt
          call lw_two_stream(ncol, nlay,                                 &
                             atmos%tau (:,:,k), atmos%ssa(:,:,k), atmos%g(:,:,k), &
                             gamma1, gamma2, Rdif, Tdif)
         !
         ! RRTMGP provides source functions at each level using the spectral mapping
         !   of each adjacent layer. Combine these for two-stream calculations
         !
         call lw_combine_sources(ncol, nlay, logical(top_is_1, wl), &
                                 lev_src_inc(:,:,k), lev_src_dec(:,:,k), &
                                 lev_source)
         !
         ! Source function for diffuse radiation
         !
          call lw_source_2str(ncol, nlay, logical(top_is_1, wl), &
                              sfc_emis(:,k), sfc_src(:,k), &
                              lay_src(:,:,k), lev_source, &
                              gamma1, gamma2, Rdif, Tdif, atmos%tau(:,:,k), &
                              source_dn(:,:,k), source_up(:,:,k), source_sfc(:,k))
        end do
      else
        do k = 1, ngpt
          call sw_two_stream(ncol, nlay, mu0,                       &
                             atmos%tau(:,:,k), atmos%ssa(:,:,k), atmos%g(:,:,k), &
                             Rdif, Tdif, Rdir, Tdir, Tnoscat)
          call sw_source_2str(ncol, nlay, logical(top_is_1, wl), Rdir, Tdir, Tnoscat, &
                              sfc_alb_dir(spectral_disc%convert_gpt2band(k), :), &
                              source_up(:,:,k), source_dn(:,:,k), source_sfc(:,k), flux_dn_dir(:,:,k))
        end do
      end if
    type is (ty_optical_props_nstr)
      call stop_on_err("Haven't implemented source calculations for multi-stream inputs")
  end select

  call write_sources(fileName, source_up, source_dn, source_sfc)

end program test_two_stream
