subroutine stop_on_err(msg)
  !
  ! Print error message and stop
  !
  use iso_fortran_env, only : error_unit
  character(len=*), intent(in) :: msg
  if(len_trim(msg) > 0) then
    write (error_unit,*) trim(msg)
    write (error_unit,*) "test_two_stream stopping"
    stop
  end if
end subroutine
! ----------------------------------------------------------------------------------
program test_two_stream
  use mo_rte_kind,      only: wp, wl
  use mo_spectral_disc, only: ty_spectral_disc
  use mo_optical_props, only: ty_optical_props_arry, &
                              ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
  use radiation_two_stream,    &
                        only: calc_reflectance_transmittance_sw, calc_two_stream_gammas_sw, &
                              calc_reflectance_transmittance_lw, calc_two_stream_gammas_lw, &
                              calc_no_scattering_transmittance_lw

  use mo_test_files_io, only: read_optical_props, read_spectral_disc, &
                              is_sw, is_lw, read_direction, &
                              read_sw_bc, read_sw_solar_sources,  &
                              read_lw_bc, read_lw_Planck_sources, &
                              write_two_stream, write_sources
  implicit none
  ! ----------------------------------------------------------------------------------
  integer :: ncol, nlay, ngpt

  character(len=128) :: fileName = 'rrtmgp-inputs-outputs.nc'

  class(ty_optical_props_arry), allocatable :: atmos

  ! SW-specific
  real(wp), dimension(:    ), allocatable :: mu0, tsi
  real(wp), dimension(:,  :), allocatable :: sfc_alb, sfc_alb_dir, sfc_alb_dif, toa_src
  real(wp), dimension(:,:,:), allocatable :: Rdir, Tdir, Tnoscat
  real(wp), dimension(:,:,:), allocatable :: flux_dn_dir
  real(wp)                                :: tsi_scaling = -999._wp
  ! LW-specific
  real(wp), dimension(:,:  ), allocatable :: sfc_emis, sfc_src, emis_sfc_bnd
  real(wp), dimension(:    ), allocatable :: t_sfc
  real(wp), dimension(:,:,:), allocatable, &
                              target      :: lay_src, lev_src_inc, lev_src_dec, planck_hl
  ! Generic
  real(wp), dimension(:,:,:), allocatable :: Rdif, Tdif
  real(wp), dimension(:,:,:), allocatable :: source_dn, source_up
  real(wp), dimension(:,  :), allocatable :: source_sfc

  type(ty_spectral_disc) :: spectral_disc

  integer                                 :: jcol, jlev, jgpt
  real(wp), dimension(:), allocatable     :: gamma1, gamma2, gamma3
  real(wp), parameter :: pi = acos(-1._wp)
  integer :: ngpts_per_band
  logical :: do_sw, top_at_1
  ! ----------------------------------------------------------------------------------
  call read_direction(fileName, top_at_1)
  if(.not. top_at_1) call stop_on_err("atmosphere has to be ordered top to bottom (top_at_1 = .false.) for ECRAD")

  call read_spectral_disc(fileName, spectral_disc)
  call read_optical_props(fileName, atmos)
  do_sw = is_sw(fileName)
  ncol = atmos%get_ncol()
  nlay = atmos%get_nlay()
  ngpt = atmos%get_ngpt()
  allocate(Rdif(ncol,nlay,ngpt), Tdif(ncol,nlay,ngpt), &
           source_up(ncol,nlay,ngpt), source_dn(ncol,nlay,ngpt), source_sfc(ncol, ngpt), &
           gamma1(ncol), gamma2(ncol), gamma3(ncol))

  if(do_sw) then
    allocate(Rdir(ncol,nlay,ngpt), Tdir(ncol,nlay,ngpt), Tnoscat(ncol,nlay,ngpt), flux_dn_dir(ncol, nlay+1, ngpt), &
             sfc_alb(ncol,ngpt))
    call read_sw_bc (fileName, mu0, tsi, tsi_scaling, sfc_alb_dir, sfc_alb_dif)
    call read_sw_solar_sources(fileName, toa_src)
    mu0 = cos(mu0 * acos(-1._wp)/180.)
    do jgpt = 1, ngpt
      sfc_alb(1:ncol, jgpt) = sfc_alb_dir(spectral_disc%convert_gpt2band(jgpt), :)
      flux_dn_dir(1:ncol,1,jgpt) = toa_src(1:ncol,jgpt) * mu0(1:ncol)
    end do
  else
    allocate(sfc_emis(ncol, ngpt), planck_hl(ncol, nlay+1, ngpt))
    call read_lw_Planck_sources(fileName, lay_src, lev_src_inc, lev_src_dec, sfc_src)
    call read_lw_bc            (fileName, t_sfc, emis_sfc_bnd)
    do jgpt = 1, ngpt
      sfc_emis(1:ncol, jgpt) = emis_sfc_bnd(spectral_disc%convert_gpt2band(jgpt), :)
    end do
  end if

  select type(atmos)
    type is (ty_optical_props_2str)
      if(do_sw) then
        do jgpt = 1, ngpt
          do jlev = 1, nlay
            call calc_two_stream_gammas_sw(ncol, 0, &
                 &  mu0, atmos%ssa(:,jlev,jgpt), atmos%g(:,jlev,jgpt), &
                 &  gamma1, gamma2, gamma3)
            call calc_reflectance_transmittance_sw(ncol, &
                 &  mu0, atmos%tau(:,jlev,jgpt), atmos%ssa(:,jlev,jgpt), &
                 &  gamma1, gamma2, gamma3, &
                 &  Rdif(:,jlev,jgpt), Tdif(:,jlev,jgpt), &
                 &  Rdir(:,jlev,jgpt), Tdir(:,jlev,jgpt), &
                 &  Tnoscat(:,jlev,jgpt) )
           !
           ! Solar direct beam and surface source
           !
           source_up(:,jlev,jgpt) =        Rdir(:,jlev,jgpt) * flux_dn_dir(:,jlev,jgpt)
           source_dn(:,jlev,jgpt) =        Tdir(:,jlev,jgpt) * flux_dn_dir(:,jlev,jgpt)
           flux_dn_dir(:,jlev+1,jgpt) = Tnoscat(:,jlev,jgpt) * flux_dn_dir(:,jlev,jgpt)
          end do
        end do
        source_sfc(1:ncol,1:ngpt) = sfc_alb(1:ncol,1:ngpt) * flux_dn_dir(1:ncol,nlay+1,1:jgpt)
        call write_two_stream(fileName, Rdif, Tdif, source_up, source_dn, Tnoscat)
      else
        ! --------------
        ! See lw_combine_sources in mo_rte_solver_kernels
        do jgpt = 1, ngpt
          do jcol = 1, ncol
            jlev = 1
            planck_hl(jcol, jlev, jgpt) =        lev_src_dec(jcol, jlev, jgpt) * pi
            do jlev = 2, nlay
              planck_hl(jcol, jlev, jgpt) = sqrt(lev_src_dec(jcol, jlev, jgpt) * &
                                                 lev_src_inc(jcol, jlev-1, jgpt)) * pi
            end do
            jlev = nlay+1
            planck_hl(jcol, jlev, jgpt) =        lev_src_inc(jcol, jlev-1, jgpt) * pi
          end do
        end do
        do jgpt = 1, ngpt
          do jlev = 1, nlay
            call calc_two_stream_gammas_lw(ncol, 0, &
                 &  atmos%ssa(:,jlev,jgpt), atmos%g(:,jlev,jgpt), &
                 &  gamma1, gamma2)
            call calc_reflectance_transmittance_lw(ncol, &
                 &  atmos%tau(:,jlev,jgpt),              &
                 &  gamma1, gamma2, planck_hl(:,jlev,jgpt), planck_hl(:,jlev+1,jgpt), & ! Planck top, bottom
                 &  Rdif(:,jlev,jgpt), Tdif(:,jlev,jgpt), &
                 &  source_up(:,jlev,jgpt), source_dn(:,jlev,jgpt))
          end do
        end do
        source_sfc(1:ncol,1:ngpt) = sfc_emis(1:ncol,1:ngpt) * sfc_src(1:ncol,1:ngpt) * pi
        call write_two_stream(fileName, Rdif, Tdif)
      end if

    type is (ty_optical_props_1scl)
      if(do_sw) then
        call stop_on_err("No source calculation for direct beam")
      else
        do jgpt = 1, ngpt
          do jlev = 1, nlay
            call calc_no_scattering_transmittance_lw(ncol,       &
                 &  atmos%tau(:,jlev,jgpt),                      & ! diffusivity angle applied in source routine
                 &  lev_src_dec(:,jlev,jgpt), lev_src_inc(:,jlev,jgpt), Tdif(:,jlev,jgpt), & ! Planck top, bottom
                 &  source_up(:,jlev,jgpt), source_dn(:,jlev,jgpt))
          end do
        end do
        source_sfc(1:ncol,1:ngpt) = sfc_emis(1:ncol,1:ngpt) * sfc_src(1:ncol,1:ngpt)
      end if

    type is (ty_optical_props_nstr)
      call stop_on_err("Haven't implemented two-stream calculations for multi-stream inputs")
  end select

  call write_sources(fileName, source_up, source_dn, source_sfc)
end program test_two_stream
