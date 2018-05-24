!-----------------------------------------------------------------------
!
! title: SHDOMPP driver program
!
! description: Reads RRTMGP optical properties from a netcdf file,
!   sets up and calls SOLVE_SHDOMPP, and write the fluxes for each
!   g-point out to the netcdf file.  Due to RRTMGP longwave source
!   function method (layer sources and upwelling and downwelling level
!   sources), in the longwave SHDOMPP uses 2*nlay layers, while using
!   nlay layers for the shortwave.
!
!-----------------------------------------------------------------------

program shdompp_solve
  use mo_rte_kind,   only: wp
  use mo_optical_props, only: ty_optical_props,      &
                              ty_optical_props_arry, &
                              ty_optical_props_1scl, &
                              ty_optical_props_2str, &
                              ty_optical_props_nstr
  use mo_spectral_disc, only: ty_spectral_disc
  use mo_fluxes_byband, only: ty_fluxes_byband
  use mo_test_files_io, only: is_sw, read_direction, &
                              read_optical_props, read_spectral_disc, &
                              read_lw_Planck_sources, read_lw_bc, read_lw_rt, &
                              read_sw_bc, read_sw_solar_sources, &
                              write_gpt_fluxes, write_fluxes, write_dir_fluxes
  use shdompp_rrtm, only: SOLVE_SHDOMPP

  implicit none
  integer :: ncol, nlay, ngpt, nband, nmom, nang
  character(len=128) :: fileName = 'shdompp-rrtmgp-inputs-outputs.nc'

  class(ty_optical_props_arry), allocatable :: atmos
  real(wp), dimension(:,:,:), allocatable :: lay_source, &
                                             lev_source_inc, lev_source_dec
  real(wp), dimension(:,  :), allocatable :: sfc_emis, sfc_source
  real(wp), dimension(:    ), allocatable :: t_sfc
  real(wp), dimension(:    ), allocatable :: mu0, tsi
  real(wp), dimension(:,  :), allocatable :: sfc_alb_dir, sfc_alb_dif
  real(wp), dimension(:,:),   allocatable :: toa_source  ! TOA incident radiation (ncol, ngpt)
  real(wp), dimension(:,:,:), allocatable :: gpt_flux_up, gpt_flux_dn,               gpt_flux_dn_dir
  real(wp), dimension(:,:  ), target, &
                              allocatable ::     flux_up,     flux_dn,     flux_net,     flux_dir
  real(wp), dimension(:,:,:), target, &
                              allocatable :: bnd_flux_up, bnd_flux_dn, bnd_flux_net, bnd_flux_dir
  type(ty_spectral_disc)                :: spectral_disc
  type(ty_fluxes_byband)                :: fluxes
  real(wp)                              :: tsi_scaling = -999._wp

  integer :: icol, igpt, ibnd, i, j, l
  logical :: top_at_1
  integer :: Narg
  character(len=10) :: arg

   ! SHDOMPP arrays and variables:
  real, allocatable :: TAUP(:), ALBEDOP(:), LEGENP(:,:)
  integer, allocatable :: NLEGP(:)
  real, allocatable :: PLANCKSRC(:)
  real, allocatable :: FLUXUP(:), FLUXDN(:), FLUXDIR(:)
  integer :: NLAYSH, NMU, NPHI, ORDINATESET, MAXLEG, MAXIG, MAXITER, ITER
  logical :: DELTAM, ACCELFLAG
  real    :: MAXDELTAU, SPLITACC, SOLACC, SOLCRIT
  real    :: SOLARFLUX, SOLARMU, SFCPLANCK, SFCPARMS(3)
  character(len=1) :: SRCTYPE, SFCTYPE

  ! --------------------------------------------------------------------------

   ! Get the NMU and SPLITACC inputs via command line arguments, as these are
   ! the most important in terms of flux accuracy and running time
  Narg = command_argument_count()
  if (Narg < 1) then
    print *, 'Usage: shdompp_solve NMU [SPLITACC]'
    print *, 'NMU is number of zenith angle streams over both hemispheres (>=4)'
    print *, 'SPLITACC is SHDOMPP layer splitting accuracy (smaller is more accurate, try 0.001 to 0.01)'
    call stop_on_err('Usage: shdompp_solve NMU [SPLITACC]')
  end if
  call get_command_argument(1, arg)
  read (arg,'(I2)') NMU
  if (Narg > 1) then
    call get_command_argument(2, arg)
    read (arg,'(F6.4)') SPLITACC
    SPLITACC = MAX(SPLITACC,0.0001)
  end if


  if (is_sw(fileName)) then
    SRCTYPE = 'S'
  else
    SRCTYPE = 'T'
  end if

  call read_optical_props (fileName, atmos)

   ! Get ncol, nlay, ngpt, and allocate output g-point flux arrays
  ncol = atmos%get_ncol()
  nlay = atmos%get_nlay()
  ngpt = atmos%get_ngpt()
  allocate (gpt_flux_up(ncol,nlay+1,ngpt), gpt_flux_dn(ncol,nlay+1,ngpt), &
            gpt_flux_dn_dir(ncol,nlay+1,ngpt))

   ! Read in gpt2band
  call read_spectral_disc(fileName, spectral_disc)
  nband = spectral_disc%get_nband()

  call read_direction (fileName, top_at_1)

  if (SRCTYPE == 'S') then
    call read_sw_bc (fileName, mu0, tsi, tsi_scaling, sfc_alb_dir, sfc_alb_dif)
    mu0 = cos(mu0 *acos(-1._wp)/180.)
    call read_sw_solar_sources(fileName, toa_source)
  else
    call read_lw_Planck_sources (fileName, lay_source, lev_source_inc, lev_source_dec, &
                                  sfc_source)
    call read_lw_bc (fileName, t_sfc, sfc_emis)
    call read_lw_rt (fileName, nang)
  end if

   ! Print out some useful information about what we are going to do
  print *, 'ncol, nlay, ngpt, nband: ',ncol, nlay, ngpt, nband
  select type (atmos)
    class is (ty_optical_props_1scl)
      print *, 'Scalar optical properties input (tau only)'
      nmom = 0
    class is (ty_optical_props_2str)
      print *, 'Two-stream optical properties input'
      if (SRCTYPE == 'S') then
        print *, 'Rayleigh phase function assumed for layers with zero asymmetry parameter,'
        print *, '  otherwise Henyey-Greenstein phase function used.'
      else
        print *, 'Henyey-Greenstein phase function used'
      end if
      nmom = 0
      DELTAM = .false.
    class is (ty_optical_props_nstr)
      print *, 'N-stream optical properties input'
      nmom = size(atmos%p, 1)
      if (nmom < 4) then
        call stop_on_err('Too few phase function moments: nmom<4')
      else
        print '(A,I2)', 'nmom=',nmom
      end if
      DELTAM = .true.
  end select
  if (SRCTYPE == 'T') &
    print '(A,I3,A)', 'Using ',2*nlay,' SHDOMPP computational layers for longwave'


  ! Start of SHDOMPP related code
    ! max layer optical depth for the "base grid"
  MAXDELTAU = MIN(0.1,MAX(1.0,100.*SPLITACC))
  ORDINATESET = 3   ! 3 for Gaussian quadrature from mu=0 to 1,
                    ! 4 for radiance-to-flux quadrature that RRTM uses for LW
                    ! Warning: only use ORDINATESET=4 for no-scattering longwave

   ! MAXLEG is order of phase function Legendre series:
   !   for two-stream input use NMU terms generated using H-G from asymmetry
   !   parameter, for n-stream input use the nmom terms input in atmos%p
  if (nmom == 0) then
    MAXLEG = NMU
  else
    MAXLEG = nmom
  end if
  if (SRCTYPE == 'S') then
    NPHI = 2*NMU
    NLAYSH = nlay
  else
    NPHI = 1
    NLAYSH = 2*nlay  ! two SHDOMPP computational layers for each RRTMGP layer
  end if
  allocate (TAUP(NLAYSH), ALBEDOP(NLAYSH), NLEGP(NLAYSH), LEGENP(0:MAXLEG,NLAYSH))
  allocate (PLANCKSRC(NLAYSH+1))
  allocate (FLUXUP(NLAYSH+1), FLUXDN(NLAYSH+1), FLUXDIR(NLAYSH+1))
  SFCTYPE = 'L'
  MAXIG = 10000
  SOLACC = 1.0E-6
  MAXITER = 100
  ACCELFLAG = .TRUE.

   ! Loop over columns
  do icol = 1, ncol
    ! Loop over the k-distribution g points for all bands
   do igpt = 1, ngpt
    ibnd = spectral_disc%convert_gpt2band(igpt)
     ! The g-point weights are in the source terms.
    if (SRCTYPE == 'S') then
      SOLARFLUX = mu0(icol)*toa_source(icol,igpt)
      SOLARMU = -mu0(icol)
      SFCPARMS(1) = sfc_alb_dir(ibnd,icol)
      if (sfc_alb_dif(ibnd,icol) /= sfc_alb_dir(ibnd,icol)) &
        print *, 'Warning: ignoring that diffuse surface albedo is different from direct albedo',ibnd,icol
    else
       ! Use average of downwelling and upwelling LW sources for edge levels
       ! and lay_source for middle levels
      if (top_at_1) then
        PLANCKSRC(3:NLAYSH-1:2) = sqrt(lev_source_dec(icol,2:nlay  ,igpt) &
                                      *lev_source_inc(icol,1:nlay-1,igpt))
        PLANCKSRC(1)        = lev_source_dec(icol,1,   igpt)
        PLANCKSRC(NLAYSH+1) = lev_source_inc(icol,nlay,igpt)
        PLANCKSRC(2:NLAYSH:2) = lay_source(icol,1:nlay,igpt)
      else
        PLANCKSRC(3:NLAYSH-1:2) = sqrt(lev_source_dec(icol,nlay  :2:-1,igpt) &
                                      *lev_source_inc(icol,nlay-1:1:-1,igpt))
        PLANCKSRC(1)        = lev_source_inc(icol,nlay,igpt)
        PLANCKSRC(NLAYSH+1) = lev_source_dec(icol,1,   igpt)
        PLANCKSRC(2:NLAYSH:2) = lay_source(icol,nlay:1:-1,igpt)
      end if
      SFCPARMS(1) = 1.0 - sfc_emis(ibnd,icol)
      SFCPLANCK = sfc_source(icol,igpt)
    end if

     ! Use top_at_1 to decide on reversing levels in arrays:
     !   SHDOMPP goes from top down, while RRTMGP can handle either way
     ! Transfer the optical depth profile
    if (top_at_1) then
      if (SRCTYPE == 'S') then
        TAUP(:) = atmos%tau(icol,1:nlay,igpt)
      else
        TAUP(1:NLAYSH-1:2) = 0.5*atmos%tau(icol,1:nlay,igpt)
        TAUP(2:NLAYSH  :2) = 0.5*atmos%tau(icol,1:nlay,igpt)
      end if
    else
      if (SRCTYPE == 'S') then
        TAUP(:) = atmos%tau(icol,nlay:1:-1,igpt)
      else
        TAUP(1:NLAYSH-1:2) = 0.5*atmos%tau(icol,nlay:1:-1,igpt)
        TAUP(2:NLAYSH:2) = 0.5*atmos%tau(icol,nlay:1:-1,igpt)
      end if
    end if
     ! Transfer the single scattering albedo and phase function moments profiles
    select type (atmos)
      class is (ty_optical_props_1scl)
        ALBEDOP(:) = 0.0
        LEGENP(0,:) = 1.0
        LEGENP(1:MAXLEG,:) = 0.0
      class is (ty_optical_props_2str)
        LEGENP(0,:) = 1.0
        do i = 1, nlay
          if (top_at_1) then
            j = i
          else
            j = nlay + 1 - i
          end if
          if (SRCTYPE == 'S') then
            ALBEDOP(i) = atmos%ssa(icol,j,igpt)
             ! If shortwave and g=0 then assume Rayleigh scattering
            if (atmos%g(icol,j,igpt) == 0.0_wp) then
              LEGENP(1,i) = 0.0
              LEGENP(2,i) = 0.5
              LEGENP(3:MAXLEG,i) = 0.0
            else
               ! If only have g, assume Henyey-Greenstein phase function
              do l = 1, MAXLEG
                LEGENP(l,i) = (2*l+1)*atmos%g(icol,j,igpt)**l
              end do
            end if
          else
            ALBEDOP(2*i-1:2*i) = atmos%ssa(icol,j,igpt)
            do l = 1, MAXLEG
              LEGENP(l,2*i-1:2*i) = (2*l+1)*atmos%g(icol,j,igpt)**l
            end do
          end if
        end do
      class is (ty_optical_props_nstr)
        LEGENP(0,:) = 1.0
        do i = 1, nlay
          if (top_at_1) then
            j = i
          else
            j = nlay + 1 - i
          end if
          if (SRCTYPE == 'S') then
            ALBEDOP(i) = atmos%ssa(icol,j,igpt)
             ! Convert from g_l form of moments to Legendre coefficients
            do l = 1, MAXLEG
              LEGENP(l,i) = (2*l+1)*atmos%p(l,icol,j,igpt)
            end do
          else
            ALBEDOP(2*i-1:2*i) = atmos%ssa(icol,j,igpt)
            do l = 1, MAXLEG
              LEGENP(l,2*i-1:2*i) = (2*l+1)*atmos%p(l,icol,j,igpt)
            end do
          end if
        end do
    end select
    NLEGP(:) = MAXLEG

     ! Call SHDOMPP routine for one column
    call SOLVE_SHDOMPP (NLAYSH, MAXLEG, PLANCKSRC, TAUP, ALBEDOP, &
                        NLEGP, LEGENP, DELTAM, SRCTYPE, SOLARFLUX, SOLARMU, &
                        SFCTYPE, SFCPARMS, SFCPLANCK, &
                        NMU, NPHI, ORDINATESET, &
                        MAXIG, MAXDELTAU, SPLITACC, SOLACC, &
                        MAXITER, ACCELFLAG, &
                        SOLCRIT, ITER, FLUXUP, FLUXDN, FLUXDIR)
    if (top_at_1) then
      if (SRCTYPE == 'S') then
        gpt_flux_up(icol,:,igpt) = FLUXUP(1:NLAYSH+1)
        gpt_flux_dn(icol,:,igpt) = FLUXDN(1:NLAYSH+1)
        gpt_flux_dn_dir(icol,:,igpt) = FLUXDIR(1:NLAYSH+1)
      else
        gpt_flux_up(icol,:,igpt) = FLUXUP(1:NLAYSH+1:2)
        gpt_flux_dn(icol,:,igpt) = FLUXDN(1:NLAYSH+1:2)
      end if
    else
      if (SRCTYPE == 'S') then
        gpt_flux_up(icol,:,igpt) = FLUXUP(NLAYSH+1:1:-1)
        gpt_flux_dn(icol,:,igpt) = FLUXDN(NLAYSH+1:1:-1)
        gpt_flux_dn_dir(icol,:,igpt) = FLUXDIR(NLAYSH+1:1:-1)
      else
        gpt_flux_up(icol,:,igpt) = FLUXUP(NLAYSH+1:1:-2)
        gpt_flux_dn(icol,:,igpt) = FLUXDN(NLAYSH+1:1:-2)
      end if
    end if
   end do
  end do

   ! Write out the g-point fluxes
  if (SRCTYPE == 'S') then
    call write_gpt_fluxes (fileName, gpt_flux_up, gpt_flux_dn, gpt_flux_dn_dir)
  else
    call write_gpt_fluxes (fileName, gpt_flux_up, gpt_flux_dn)
  end if

  ! Write out the derived fluxes
  allocate (    flux_up(ncol,nlay+1     ),     flux_dn(ncol,nlay+1     ), &
                flux_net(ncol,nlay+1     ))
  allocate (bnd_flux_up(ncol,nlay+1,nband), bnd_flux_dn(ncol,nlay+1,nband), &
            bnd_flux_net(ncol,nlay+1,nband))
  fluxes%flux_up      => flux_up
  fluxes%flux_dn      => flux_dn
  fluxes%flux_net     => flux_net
  fluxes%bnd_flux_up  => bnd_flux_up
  fluxes%bnd_flux_dn  => bnd_flux_dn
  fluxes%bnd_flux_net => bnd_flux_net
  if (SRCTYPE == 'S') then
    allocate(bnd_flux_dir(ncol,nlay+1,nband), flux_dir(ncol,nlay+1))
    fluxes%flux_dn_dir     => flux_dir
    fluxes%bnd_flux_dn_dir => bnd_flux_dir
    call stop_on_err(fluxes%reduce(gpt_flux_up, gpt_flux_dn, spectral_disc, top_at_1, &
                                          gpt_flux_dn_dir))
  else
    call stop_on_err(fluxes%reduce(gpt_flux_up, gpt_flux_dn, spectral_disc, top_at_1))
  end if

  call write_fluxes(fileName, flux_up, flux_dn, flux_net, bnd_flux_up, bnd_flux_dn, bnd_flux_net)
  if(SRCTYPE == 'S') &
    call write_dir_fluxes(fileName, flux_dir, bnd_flux_dir)
contains

! =-----------------------------------
subroutine stop_on_err(error_msg)
  use iso_fortran_env, only : error_unit
  character(len=*), intent(in) :: error_msg

  if(error_msg /= "") then
    write (error_unit,*) trim(error_msg)
    write (error_unit,*) "shdompp_solve stopping"
    stop
  end if
end subroutine stop_on_err

END program
