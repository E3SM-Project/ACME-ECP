! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2015,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!
! Description: Read in gas optical properties, calculate Mie cloud optical
!   properties for a variety of cloud properties (specified in parameter 
!   statements for ncloud columns), make the cloudy column optical properties
!   for the first ncloudcol (input parameter) columns, add the cloudy and
!   gas optical properties, and output the optical properties in either 
!   two-stream or n-stream format (nmom=2 on command line for two-stream,
!   nmom>=4 for n-stream).  The two-stream cloud optical properties are 
!   delta scaled, but the n-stream ones are not.


program add_cloud_optics
  use mo_rte_kind,   only: wp
  use mo_test_files_io, only: read_atmos, is_lw, & 
                              read_spectral_disc, read_optical_props, & 
                              write_optical_props
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_optical_props, only: ty_optical_props, &
                              ty_optical_props_arry, & 
                              ty_optical_props_2str, &
                              ty_optical_props_nstr
  use mo_spectral_disc, only: ty_spectral_disc
  use MIE_GAMMA_WATER,  only: LIQCLOUD_MIE_CALCS

  implicit none

  real(wp), dimension(:,:), allocatable :: p_lay, t_lay, p_lev, t_lev
  real(wp), dimension(:),   allocatable :: t_sfc
  real(wp), dimension(:,:), allocatable :: emis_sfc
  real(wp), dimension(:,:), allocatable :: col_dry
  type(ty_gas_concs)                    :: gas_concs

  class (ty_optical_props_arry), allocatable :: gas_props, all_props
  class (ty_optical_props_nstr), allocatable :: gas_nstr_props 
  type(ty_spectral_disc)                   :: spectral_disc

  ! dimensions
  integer :: ncol, nlay, ngpt, nbnd

   ! The number of phase function moments needed is equal to the number of
   ! streams in both hemispheres (NMU) for a flux RT calculation with delta-M.
   ! Use nmom=2 for two-stream optical properties output (with g), 
   ! and nmom>=4 for n-stream optical properties output (with p).
   ! nmom is specified as the first command line argument.
  integer :: nmom
  logical :: top_at_1, longwave

   ! Specify the cloud layers for three uniform cloud layers in each column
   ! (0 for no cloud): RRTMGP layer indices increasing away from the surface,
   ! cloud liquid water path (g/m^2), and effective radius (micron) for 
   ! each column.  Note: Fewer r_eff will result in faster Mie table computation.
  integer, parameter :: ncloud=12
  integer, parameter :: cloud_lay(2,3,ncloud) = reshape( &
     [1,1,0,0,0,0,    4,4,0,0,0,0,   24,24,0,0,0,0,  4,8,0,0,0,0, &
      14,16,0,0,0,0,  23,24,0,0,0,0,  4,4,5,5,6,6,   4,6,12,14,0,0, &
      3,9,10,16,17,22, 1,3,12,15,0,0, 4,7,13,15,23,24, 9,11,13,14,0,0 ], &
     [2,3,ncloud] )
  real(wp), parameter :: cloud_lwp(3,ncloud) = reshape( &
      [50.,0.,0.,      100.,0.,0.,     1.0,0.,0.,      600.,0.,0., &      
       100.,0.,0.,      2.0,0.,0.,     100.,200.,300., 300.,100.,0., &
       900.,500.,100.,  100.,300.,0.,  300.,100.,10.,  200.,100.,0. ], &
     [3,ncloud] )
  real(wp), parameter :: cloud_reff(3,ncloud) = reshape( &
      [8.0,0.,0.,       8.0,0.,0.,     25.0,0.,0.,      10.0,0.,0., &
       8.0,0.,0.,       25.0,0.,0.,     8.0,10.0,12.0,  10.0,8.0,0., &
       10.0,12.0,12.0,  8.0,10.0,0.,   12.0,8.0,25.0,   10.0,8.0,0. ], &
     [3,ncloud] )

   ! LIQCLOUD_MIE_CALCS parameters:
  real, parameter :: MAXRADIUS=30.0  ! max cloud droplet radius
  real, parameter :: ALPHA=7.0  ! for effective variance of 0.1 for gamma distributions of liquid cloud droplets
  integer :: NREFF, MAXLEG
  integer, allocatable :: NLEG(:,:)
  real    :: WAVELEN1, WAVELEN2, DELTAWAVE
  real, allocatable :: re(:), REFF(:), EXTINCT(:,:), SSALB(:,:), LEGCOEF(:,:,:)
  character(len=1) :: AVGFLAG = 'A'
   ! Wavelength steps for averaging Mie calculation (microns)
  real, parameter :: deltawave_sw(14) = &
      (/ 0.2,0.05,0.05,0.05,0.01,0.05,0.05,0.05,0.05,0.05,0.05,0.04,0.04,0.03 /)
  real, parameter :: deltawave_lw(16) = &
      (/ 5.0,2.0,1.0,0.2,0.2,0.2,0.2,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1 /)

  integer :: ncloudcol, icolcld, icol, ire, ilay, ilay1, ilay2
  integer :: ibnd, igpt, igpt1, igpt2, l, j, gpt_lims(2) 
  real(wp), allocatable :: tau_cloud(:,:,:), ssa_cloud(:,:,:), p_cloud(:,:,:,:)
  real(wp), allocatable :: fscale(:,:,:)
  real(wp), allocatable :: band_lims_wavenumber(:,:) 
  integer :: Narg
  character(len=8) :: arg1

  character(len=64), parameter :: fileName = 'rrtmgp-inputs-outputs.nc'
  ! ==========================================================

  Narg = command_argument_count()
  if (Narg < 1) then
    print *, 'Usage: add_cloud_optics Nmom [Ncloudcol]'
    print *, 'Nmom is 2 for two-stream gas+cloud optical properties output.'
    print *, 'Nmom is >=4 for N-stream gas+cloud optical properties output.'
    print '(A,I2,A)', ' Ncloudcol is number of cloudy columns made from Ncloud=',ncloud,' clouds available (0 to ncol)'
    stop
  end if
  call get_command_argument(1, arg1)
  read (arg1,'(I2)') nmom
  nmom = max(nmom,2)
  ncloudcol = ncloud
  if (Narg > 1) then
    call get_command_argument(2, arg1)
    read (arg1,'(I2)') ncloudcol
  endif
  print '(A,I2,A,I2)', ' Arguments: Nmom=',nmom, '   Ncloudcol=',ncloudcol

   ! load atmosphere profiles
  call read_atmos (fileName, p_lay, t_lay, p_lev, t_lev, gas_concs, col_dry)

  ncol = size(p_lay, 1)
  nlay = size(p_lay, 2)
  p_lay = p_lay/100._wp
  p_lev = p_lev/100._wp
  top_at_1 = p_lay(1, 1) < p_lay(1, nlay)
  longwave = is_lw(fileName)

   ! read band wavenumber ranges and conversion to g-points
  call read_spectral_disc(fileName, spectral_disc)
  nbnd = spectral_disc%get_nband() 
  allocate(band_lims_wavenumber(2, nbnd)) 
  band_lims_wavenumber = spectral_disc%get_band_lims_wavenumber() 
  
   ! read the gas optical properties
  call read_optical_props (fileName, gas_props)
  ngpt = gas_props%get_ngpt()
  if(ngpt /= spectral_disc%get_ngpt()) call stop_on_err("Number of g-points inconsistent") 
  print *, 'ncol, nlay, ngpt:',ncol, nlay, ngpt

  if (nmom <= 2) then
      ! Make the cloud/total optical properties in the two-stream structure 
    allocate(ty_optical_props_2str::all_props)
  else 
      ! Make the cloud/total optical properties in the nstream structure 
    allocate(ty_optical_props_nstr::all_props)
  endif

  select type (all_props)
    type is (ty_optical_props_2str) ! two-stream calculation
        call stop_on_err(all_props%init_2str(ncol, nlay, ngpt))
        print *, "Initializing output two-stream optical properties" 
    type is (ty_optical_props_nstr) 
        call stop_on_err(all_props%init_nstr(nmom, ncol, nlay, ngpt))
        print *, "Initializing output n-stream optical properties" 
  end select 

  allocate (tau_cloud(ncol,nlay,ngpt), ssa_cloud(ncol,nlay,ngpt))
  allocate (p_cloud(nmom,ncol,nlay,ngpt))
  tau_cloud(:,:,:) = 0.0_wp ; ssa_cloud(:,:,:) = 0.0_wp ; p_cloud(:,:,:,:) = 0.0_wp 

   ! Make a list of unique cloud effective radius in cloud_reff
  allocate (re(3*ncloud))
  NREFF = 1
  re(NREFF) = cloud_reff(1,1)
  do icol = 1, MIN(ncloud,ncloudcol)
    do j = 1, 3
      if (cloud_reff(j,icol) > 0.0) then
        if (ALL(re(1:NREFF) /= cloud_reff(j,icol))) then
          NREFF = NREFF + 1
          re(NREFF) = cloud_reff(j,icol)
        end if
      end if
    end do 
  end do 

  MAXLEG = nmom
  allocate (REFF(NREFF), EXTINCT(NREFF,nbnd), SSALB(NREFF,nbnd))
  allocate (NLEG(NREFF,nbnd), LEGCOEF(0:MAXLEG,NREFF,nbnd))
  REFF(1:NREFF) = re(1:NREFF)

  if (ncloudcol > 0) then
   ! Loop over the bands calling LIQCLOUD_MIE_CALCS for each band
    do ibnd = 1, nbnd
      WAVELEN1 = 10000./band_lims_wavenumber(2,ibnd)
      WAVELEN2 = 10000./band_lims_wavenumber(1,ibnd)
      if (longwave) then 
        DELTAWAVE = deltawave_lw(ibnd)
      else
        DELTAWAVE = deltawave_sw(ibnd)
      endif
      call LIQCLOUD_MIE_CALCS (WAVELEN1, WAVELEN2, AVGFLAG, DELTAWAVE, &
                               ALPHA, NREFF, REFF, MAXRADIUS, MAXLEG, &
                               EXTINCT(:,ibnd), SSALB(:,ibnd), NLEG(:,ibnd), &
                               LEGCOEF(:,:,ibnd))
    end do
  end if

   ! Use the Mie tables to make the cloud optical property profiles for
   ! each uniform cloud in each column.  If ncloudcol > ncloud, then some 
   ! cloud columns in the database are resused.
  do icol = 1, ncloudcol
   do j = 1, 3
    icolcld = MOD(icol-1,ncloud)+1
    ilay1 = cloud_lay(1,j,icolcld)
    ilay2 = cloud_lay(2,j,icolcld)
    if (ilay1 > 0 .and. ilay2 > 0) then
      if (ilay2 < ilay1) stop 'add_cloud_optics: cloud_lay(2,j,icolcld)<cloud_lay(1,j,icolcld)'
       ! Find the index in the Mie table that matches this cloud r_eff
      ire = 1
      do while (REFF(ire) /= cloud_reff(j,icolcld) .and. ire <= NREFF)
        ire = ire + 1
      end do 
      if (ire > NREFF) stop 'add_cloud_optics: cloud_reff not in REFF'

       ! Calculate the layer optical properties from the Mie table for each band.
       !   The cloud LWP*EXTINCT gives the total cloud optical depth, which is
       !   distributed over the layers according to the delta pressure of each.
       !   Convert from Legendre phase function coefficients to "g_l" moments.
      do ibnd = 1, nbnd
        gpt_lims = spectral_disc%convert_band2gpt(ibnd)
        igpt1=gpt_lims(1); igpt2=gpt_lims(2) 
        if (top_at_1) then
          tau_cloud(icol,ilay1:ilay2,igpt1) = &
               0.001*cloud_lwp(j,icolcld)*EXTINCT(ire,ibnd) &
                *abs( (p_lev(icol,nlay+2-ilay2:ilay1:-1) - p_lev(icol,nlay+1-ilay2:ilay1:-1)) &
                    /(p_lev(icol,nlay+1-ilay2) - p_lev(icol,nlay+2-ilay1)) )
        else
          tau_cloud(icol,ilay1:ilay2,igpt1) = &
               0.001*cloud_lwp(j,icolcld)*EXTINCT(ire,ibnd) &
                *abs( (p_lev(icol,ilay1:ilay2) - p_lev(icol,ilay1+1:ilay2+1)) &
                    /(p_lev(icol,ilay1) - p_lev(icol,ilay2+1)) )
        end if
        do igpt = igpt1+1, igpt2
          tau_cloud(icol,ilay1:ilay2,igpt) = tau_cloud(icol,ilay1:ilay2,igpt1)
        end do
        ssa_cloud(icol,ilay1:ilay2,igpt1:igpt2) = SSALB(ire,ibnd)
        do l = 1, nmom
          p_cloud(l,icol,ilay1:ilay2,igpt1:igpt2) = LEGCOEF(l,ire,ibnd)/(2*l+1)
        end do
      end do
    end if
   end do
  end do
  deallocate (REFF, EXTINCT, SSALB, NLEG, LEGCOEF)

   ! Put the cloud optical properties in the desired type (two-stream
   ! or n-stream) output structure. Also deal with the input/output
   ! file having either top down or bottom up arrays.
  allocate (fscale(ncol,nlay,ngpt))
  if (top_at_1) then
    select type (all_props)
      type is (ty_optical_props_2str) ! two-stream output
        all_props%tau(:,nlay:1:-1,:) = tau_cloud(:,1:nlay,:)
        all_props%ssa(:,nlay:1:-1,:) = ssa_cloud(:,1:nlay,:)
        all_props%g(:,nlay:1:-1,:) = p_cloud(1,:,1:nlay,:)
        fscale(:,nlay:1:-1,:) = p_cloud(2,:,1:nlay,:)
      type is (ty_optical_props_nstr) ! n-stream output
        all_props%tau(:,nlay:1:-1,:) = tau_cloud(:,1:nlay,:)
        all_props%ssa(:,nlay:1:-1,:) = ssa_cloud(:,1:nlay,:)
        all_props%p(:,:,nlay:1:-1,:) = p_cloud(:,:,1:nlay,:)
    end select 
  else
    select type (all_props)
      type is (ty_optical_props_2str) ! two-stream output
        all_props%tau(:,1:nlay,:) = tau_cloud(:,1:nlay,:)
        all_props%ssa(:,1:nlay,:) = ssa_cloud(:,1:nlay,:)
        all_props%g(:,1:nlay,:) = p_cloud(1,:,1:nlay,:)
        fscale(:,1:nlay,:) = p_cloud(2,:,1:nlay,:)
      type is (ty_optical_props_nstr) ! n-stream output
        all_props%tau(:,1:nlay,:) = tau_cloud(:,1:nlay,:)
        all_props%ssa(:,1:nlay,:) = ssa_cloud(:,1:nlay,:)
        all_props%p(:,:,1:nlay,:) = p_cloud(:,:,1:nlay,:)
    end select 
  end if

   ! Delta scale the two-stream output cloud optical properties
  select type (all_props)
    type is (ty_optical_props_2str)
      print *, 'Delta-scaling the two-stream cloud optical properties'
       ! Delta scaling with f = g^2:
      ! call stop_on_err (all_props%delta_scale())
       ! Delta scaling with f = p(2):
      call stop_on_err (all_props%delta_scale(fscale))
    type is (ty_optical_props_nstr)
      print *, 'Not delta-scaling the n-stream cloud optical properties'
  end select 
  deallocate (fscale)

   ! Add the gas optical properties to the cloud optical properties.
   !   There is a special case for n-stream output in the shortwave,
   ! where we want to use the proper Rayleigh phase function, which
   ! is not in the two-stream optical properties with g=0.
  select type (all_props)
    type is (ty_optical_props_2str)
      call stop_on_err (all_props%increment_by(gas_props))
    type is (ty_optical_props_nstr)
      if (longwave) then
        call stop_on_err (all_props%increment_by(gas_props))
      else
        print *, 'Assuming Rayleigh phase function for input gas optics in the shortwave.'
         ! Make n-stream optical properties class for the gases and allocate components
        allocate(ty_optical_props_nstr::gas_nstr_props)
        call stop_on_err(gas_nstr_props%init_nstr(nmom,ncol,nlay,ngpt))
         ! Copy the two-stream class to the n-stream class
        call stop_on_err (gas_props%get_subset(1,ncol,gas_nstr_props))
         ! Change all the gas phase functions from isotropic to Rayleigh
        gas_nstr_props%p(2,:,:,:) = 0.1_wp
         ! Finally add the gas optical properties to the cloud optical properties
        call stop_on_err (all_props%increment_by(gas_nstr_props))
      end if
  end select 

    ! Debugging print out
  if (.false.) then
    do icol = 1, ncol
      print *
      print '(A,I2)', 'Band cloud and total optical properties for column: ',icol
      do ibnd = 1, nbnd
        print '(A,I2,2(1X,F8.3))', 'First g-point for band, wavelength ranges: ', &
              ibnd, 10000./band_lims_wavenumber(2:1:-1,ibnd)
        gpt_lims = spectral_disc%convert_band2gpt(ibnd)
        igpt=gpt_lims(1) 
        print *, '                     Cloud                Total            Gas'
        print *, ' Pres  Temp     tau    ssa   asym    tau    ssa   asym     tau'
        do ilay = 1, nlay
          select type (all_props)
           type is (ty_optical_props_2str) 
             print '(1x,f6.1,1x,f5.1,1x,2(1x,f6.2,2(1x,f6.4)),2x,f6.2)', &
               p_lay(icol,ilay), t_lay(icol,ilay), &
               tau_cloud(icol,ilay,igpt), &
               ssa_cloud(icol,ilay,igpt), p_cloud(1,icol,ilay,igpt), &
               all_props%tau(icol,ilay,igpt), &
               all_props%ssa(icol,ilay,igpt), all_props%g(icol,ilay,igpt), &
               gas_props%tau(icol,ilay,igpt)
           type is (ty_optical_props_nstr) 
             print '(1x,f6.1,1x,f5.1,1x,2(1x,f6.2,2(1x,f6.4)),2x,f6.2)', &
               p_lay(icol,ilay), t_lay(icol,ilay), &
               tau_cloud(icol,ilay,igpt), &
               ssa_cloud(icol,ilay,igpt), p_cloud(1,icol,ilay,igpt), &
               all_props%tau(icol,ilay,igpt), &
               all_props%ssa(icol,ilay,igpt), all_props%p(1,icol,ilay,igpt), &
               gas_props%tau(icol,ilay,igpt)
          end select 
        end do
      end do
    end do
  end if
  deallocate (tau_cloud, ssa_cloud, p_cloud)

   ! Write fields out 
  call write_optical_props (fileName, all_props)  
  
contains
! -----------------------------------------------------------------------------------
  subroutine stop_on_err(error_msg) 
    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: error_msg

    if(error_msg /= "") then
      write (error_unit,*) trim(error_msg) 
      write (error_unit,*) "add_cloud_optics stopping" 
      stop  
    end if 

  end subroutine stop_on_err  
! -----------------------------------------------------------------------------------
end program add_cloud_optics
