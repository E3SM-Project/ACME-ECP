subroutine stop_on_err(msg)
  !
  ! Print error message and stop  
  ! 
  use iso_fortran_env, only : error_unit
  character(len=*), intent(in) :: msg
  if(len_trim(msg) > 0) then 
    write (error_unit,*) trim(msg)
    write (error_unit,*) "change_albedo_emis_sza stopping"
    stop
  end if 
end subroutine
! ----------------------------------------------------------------------------------
program change_albedo_emis_sza
  use mo_rte_kind,   only: wp
  use mo_gas_optics, &
                        only: ty_gas_optics_specification
  use mo_test_files_io, only: is_sw, read_sw_bc, read_lw_bc, &
                        read_sfc_test_file, write_sw_surface_albedo, &
                        write_lw_surface_emissivity, write_solar_zenith_angle
  use mo_load_coefficients, only: load_and_init
  implicit none 
  ! ----------------------------------------------------------------------------------
  character(len=128) :: fileName = 'rrtmgp-inputs-outputs.nc'
  character(len=128) :: sfc_test_file

  real(wp), dimension(:,:), allocatable :: sfc_alb, sfc_emis
  real(wp), dimension(:),   allocatable :: sza, tsi 
  real(wp), dimension(:,:), allocatable :: sfc_alb_dir, sfc_alb_dif
  real(wp), dimension(:),   allocatable :: t_sfc
  real(wp), dimension(:,:), allocatable :: emis_sfc 
  real(wp)                              :: tsi_scaling = -999._wp

  integer  :: ncol, nband, nspectra, icol, ispec
  logical  :: using_sfc_test_file, changeSWalbedo, changeLWemis, changeSZA
  logical  :: doing_sw
  real(wp) :: minSWalbedo, maxSWalbedo, minLWemis, maxLWemis, minSZA, maxSZA
  real(wp) :: degrad, minmu0, maxmu0
  integer  :: seed(24)
  real     :: r
  integer  :: Narg
  character(len=128) :: arg1
  character(len=12)  :: arg2, arg3

  ! ----------------------------------------------------------------------------------

    ! Initialize random number generator
  seed(1:24)=(/ 39657, 563173, 999932, 966913, 889471, 909070, 944607, 435667, &
               161840, 780810, 642433, 682461, 782474, 334050, 646162, 281542, &
               712377, 206565,  40657, 905467, 926367, 409936,  44100, 993431 /)
  call random_seed(put=seed)

  Narg = command_argument_count()  
  if (Narg < 1) then
    print *, 'Usage: change_albedo_emis_sza  rrtmgp_sfc_test_file (or none)', &
              '  [sw_albedo_min sw_albedo_max OR lw_emis_min lw_emis_max] [sza_min sza_max]'
    print *, 'If a test file is input, sets each column surface spectral SW albedo or'
    print *, 'LW emissivity to that from a random spectrum in the file.  Otherwise, the'
    print *, 'spectrally uniform albedo or emissivity are randomly chosen between the'
    print *, 'min and max values specified.  The diffuse and direct SW albedo are set to'
    print *, 'the same values.  The solar zenith angles are also randomly set between the'
    print *, 'min and max SZA (uniform in solid angle).'
    stop
  end if

  doing_sw = is_sw(fileName)

  using_sfc_test_file = .false.
  changeSWalbedo = .false.
  changeLWemis = .false.
  changeSZA = .false.
  call get_command_argument(1, arg1)
  read (arg1,'(A)') sfc_test_file
  if (TRIM(sfc_test_file) /= 'none') then
    using_sfc_test_file = .true.
    if (Narg == 3 .and. doing_sw) then
      call get_command_argument(2, arg2)
      read (arg2,*) minSZA
      call get_command_argument(3, arg3)
      read (arg3,*) maxSZA
      changeSZA = .true.
    endif
  else
    if (Narg == 1) then
      call stop_on_err('Arguments specify nothing to do')
    endif
    if (Narg >= 3) then
      call get_command_argument(2, arg2)
      call get_command_argument(3, arg3)
      if (doing_sw) then
        read (arg2,*) minSWalbedo
        read (arg3,*) maxSWalbedo
        changeSWalbedo = .true.
      else
        read (arg2,*) minLWemis
        read (arg3,*) maxLWemis
        changeLWemis = .true.
      endif
    endif
    if (Narg == 5 .and. doing_sw) then
      call get_command_argument(4, arg2)
      read (arg2,*) minSZA
      call get_command_argument(5, arg3)
      read (arg3,*) maxSZA
      changeSZA = .true.
    endif
  endif

  if (using_sfc_test_file) then
    if (doing_sw) then
      print *, 'Setting spectral surface SW albedo from input spectra file.'
      changeSWalbedo = .true.
    else
      print *, 'Setting spectral surface LW emissivity from input spectra file.'
      changeLWemis = .true.
    endif
  else
    if (changeSWalbedo) print *, &
      'Setting spectrally uniform surface SW albedo to random values between min and max input'
    if (changeLWemis) print *, &
      'Setting spectrally uniform surface LW emissivity to random values between min and max input'
  endif
  if (changeSZA) print *, &
     'Setting solar zenith angles to random values (uniform in mu) between min and max input'


  if (doing_sw) then
     ! Shortwave section
    call read_sw_bc(fileName, sza, tsi, tsi_scaling, sfc_alb_dir, sfc_alb_dif)
    nband = size(sfc_alb_dif,1)
    ncol = size(sfc_alb_dif,2)

    if (using_sfc_test_file) then
      call read_sfc_test_file (sfc_test_file, sfc_alb, sfc_emis)
      nspectra = size(sfc_alb,2)
      do icol = 1, ncol
        call random_number(r)
        ispec = int(nspectra*r)+1
        sfc_alb_dif(:,icol) = sfc_alb(:,ispec)
      enddo
    else if (changeSWalbedo) then
      do icol = 1, ncol
        call random_number(r)
        sfc_alb_dif(:,icol) = (1.0-r)*minSWalbedo + r*maxSWalbedo
      enddo
    endif

    if (changeSZA) then
       ! Make random solar zenith angles
      degrad = acos(-1.0)/180.
      minmu0 = cos(degrad*maxSZA)
      maxmu0 = cos(degrad*minSZA)
      do icol = 1, ncol
        call random_number(r)
         ! Uniform distribution in solar zenith angle
        ! sza(icol) = (1.0-r)*minSZA + r*maxSZA
         ! Uniform distribution in solid angle
        sza(icol) = acos((1.0-r)*minmu0 + r*maxmu0)/degrad
      enddo
    endif

  else
     ! Longwave section
    call read_lw_bc (fileName, t_sfc, emis_sfc)
    nband = size(emis_sfc,1)
    ncol = size(emis_sfc,2)

    if (using_sfc_test_file) then
      call read_sfc_test_file (sfc_test_file, sfc_alb, sfc_emis)
      nspectra = size(sfc_emis,2)
      do icol = 1, ncol
        call random_number(r)
        ispec = int(nspectra*r)+1
        emis_sfc(:,icol) = sfc_emis(:,ispec)
      enddo
    else if (changeLWemis) then
      do icol = 1, ncol
        call random_number(r)
        emis_sfc(:,icol) = (1.0-r)*minLWemis + r*maxLWemis
      enddo
    endif
  endif

  if (changeSWalbedo) then
    sfc_alb_dir(:,:) = sfc_alb_dif(:,:)
    call write_sw_surface_albedo (fileName, sfc_alb_dir, sfc_alb_dif)
  endif
  if (changeSZA) then
    call write_solar_zenith_angle (fileName, sza)
  endif
  if (changeLWemis) then
    call write_lw_surface_emissivity (fileName, emis_sfc)
  endif

end program change_albedo_emis_sza
