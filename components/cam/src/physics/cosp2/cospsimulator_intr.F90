module cospsimulator_intr
  ! ######################################################################################
  ! Purpose: CAM interface to
  !         Name:         CFMIP Observational Simulator Package Version 2 (COSP2)
  !         What:         Simulate ISCCP/CloudSat/CALIPSO/MISR/MODIS cloud products from
  !                       GCM inputs
  !         Version:      v2.0 (August 2017)
  !         Authors:      Dustin Swales (dustin.swales@noaa.gov)
  !
  ! Modifications:
  !
  ! ######################################################################################
  use shr_kind_mod,         only: r8 => shr_kind_r8
  use spmd_utils,           only: masterproc
  use ppgrid,               only: pcols, pver, pverp, begchunk, endchunk
  use perf_mod,             only: t_startf, t_stopf
  use cam_abortutils,       only: endrun
  use phys_control,         only: cam_physpkg_is
  use cam_logfile,          only: iulog
  use quickbeam,            only: radar_cfg
  use mod_quickbeam_optics, only: size_distribution
  use mod_cosp,             only: cosp_outputs,cosp_optical_inputs,cosp_column_inputs
  use mod_cosp_config,      only: pres_binCenters, pres_binEdges, tau_binCenters,      &
                                  tau_binEdges,        &
                                  cloudsat_binCenters, &
                                  cloudsat_binEdges,   &
                                  calipso_binCenters,  &
                                  calipso_binEdges,    &
                                  misr_histHgtCenters, &
                                  misr_histHgtEdges,   &
                                  nsza_cosp         => PARASOL_NREFL,       &
                                  PARASOL_SZA,         &
                                  nprs_cosp         => npres,               &
                                  ntau_cosp         => ntau,                &
                                  ntau_cosp_modis   => ntau,                &
                                  ndbze_cosp        => CLOUDSAT_DBZE_BINS,           &
                                  nsr_cosp          => SR_BINS,             &
                                  nhtmisr_cosp      => numMISRHgtBins,      &
                                  nhydro            => N_HYDRO,             &
                                  Nlvgrid, &
                                  R_UNDEF,PARASOL_NREFL,LIDAR_NCAT,LIDAR_NTYPE,SR_BINS, &
                                  N_HYDRO,RTTOV_MAX_CHANNELS,numMISRHgtBins,            &
                                  CLOUDSAT_DBZE_BINS,LIDAR_NTEMP,calipso_histBsct,      &
                                  CFODD_NDBZE,      CFODD_NICOD,                        &
                                  CFODD_BNDRE,      CFODD_NCLASS,                       &
                                  CFODD_DBZE_MIN,   CFODD_DBZE_MAX,                     &
                                  CFODD_ICOD_MIN,   CFODD_ICOD_MAX,                     &
                                  CFODD_DBZE_WIDTH, CFODD_ICOD_WIDTH,                   &
                                  WR_NREGIME,                                           &
                                  numMODISTauBins,numMODISPresBins,         &
                                  numMODISReffIceBins,numMODISReffLiqBins,  &
                                  numISCCPTauBins,numISCCPPresBins,         &
                                  numMISRTauBins,reffICE_binEdges,          &
                                  reffICE_binCenters,reffLIQ_binEdges,      &
                                  reffLIQ_binCenters
   use cosp_kinds, only: wp

  implicit none
  private
  save

  ! Public functions/subroutines
  public :: &
       cospsimulator_intr_readnl,  &
       cospsimulator_intr_register,&
       cospsimulator_intr_init,    &
       cospsimulator_intr_run

  ! ######################################################################################
  ! Public declarations
  ! ######################################################################################
  ! Whether to do COSP calcs and I/O, default is false. If docosp is specified in
  ! the atm_in namelist, this value is overwritten and cosp is run
  logical, public :: docosp = .false.

  ! Frequency at which cosp is called, every cosp_nradsteps radiation timestep
  integer, public :: cosp_nradsteps = 3! CAM namelist variable default, not in COSP namelist

  ! ######################################################################################
  ! Local declarations
  ! ######################################################################################
  real(r8),allocatable, target :: htlim_cosp(:,:)          ! height limits for COSP outputs (Nlvgrid+1)
  real(r8),allocatable, target :: htmid_cosp(:)            ! height midpoints of COSP radar/lidar output (Nlvgrid)
  integer, allocatable, target :: scol_cosp(:)             ! sub-column number (nSubcol)

  ! ######################################################################################
  ! Default namelists
  ! The CAM and COSP namelists defaults are set below.  Some of the COSP namelist
  ! variables are part of the CAM namelist - they all begin with "cosp_" to keep their
  ! names specific to COSP. I set their CAM namelist defaults here, not in namelist_defaults_cam.xml
  !  Variables identified as namelist variables are defined in
  !  ../models/atm/cam/bld/namelist_files/namelist_definition.xml
  ! ######################################################################################
  ! CAM
  logical :: cosp_sample_atrain = .false.    ! CAM namelist variable default, not in COSP namelist
  character(len=256) :: cosp_atrainorbitdata ! CAM namelist variable, no default, need to specify!
  logical :: cosp_amwg             = .false. ! CAM namelist variable default, not in COSP namelist
  logical :: cosp_lite             = .false. ! CAM namelist variable default, not in COSP namelist
  logical :: cosp_passive          = .false. ! CAM namelist variable default, not in COSP namelist
  logical :: cosp_active           = .false. ! CAM namelist variable default, not in COSP namelist
  logical :: cosp_isccp            = .false. ! CAM namelist variable default, not in COSP namelist
  logical :: cosp_cfmip_3hr        = .false. ! CAM namelist variable default, not in COSP namelist
  logical :: cosp_cfmip_da         = .false. ! CAM namelist variable default, not in COSP namelist
  logical :: cosp_cfmip_off        = .false. ! CAM namelist variable default, not in COSP namelist
  logical :: cosp_cfmip_mon        = .false. ! CAM namelist variable default, not in COSP namelist
  logical :: cosp_lradar_sim       = .false. ! CAM namelist variable default
  logical :: cosp_llidar_sim       = .false. ! CAM namelist variable default
  logical :: cosp_lisccp_sim       = .false. ! CAM namelist variable default
  logical :: cosp_lmisr_sim        = .false. ! CAM namelist variable default
  logical :: cosp_lmodis_sim       = .false. ! CAM namelist variable default
  logical :: cosp_histfile_aux     = .false. ! CAM namelist variable default
  logical :: cosp_lfrac_out        = .false. ! CAM namelist variable default
  integer :: cosp_ncolumns         = 50      ! CAM namelist variable default
  integer :: cosp_histfile_num     = 1       ! CAM namelist variable default, not in COSP namelist
  integer :: cosp_histfile_aux_num = -1      ! CAM namelist variable default, not in COSP namelist

  ! COSP
  logical :: cosp_lparasol_sim     = .false.      ! +cosp2
  logical :: cosp_lrttov_sim       = .false.      ! not running rttov, always set to .false.
  logical :: cosp_lgrlidar_sim     = .false.      ! Ground lidar
  logical :: cosp_latlid_sim       = .false.

  ! ######################################################################################
  ! COSP parameters
  ! ######################################################################################
  ! Note: Unless otherwise specified, these are parameters that cannot be set by the CAM namelist.
  integer :: nSubcol = 50                        ! Number of subcolumns in SCOPS (50), can be changed from default by CAM namelist
  integer :: Nlr = 40                            ! Number of levels in statistical outputs
                                                 ! (only used if USE_VGRID=.true.)  (40)
  logical :: use_vgrid = .true.                  ! Use fixed vertical grid for outputs?
                                                 ! (if .true. then define # of levels with nlr)  (.true.)
  logical :: csat_vgrid = .true.                 ! CloudSat vertical grid?
                                                 ! (if .true. then the CloudSat standard grid is used.
                                                 ! If set, overides use_vgrid.) (.true.)
  ! namelist variables for COSP input related to radar simulator
  real(wp) :: radar_freq = 94.0_r8               ! CloudSat radar frequency (GHz) (94.0)
  integer :: surface_radar = 0                   ! surface=1, spaceborne=0 (0)
  integer :: use_mie_tables = 0                  ! use a precomputed lookup table? yes=1,no=0 (0)
  integer :: use_gas_abs = 1                     ! include gaseous absorption? yes=1,no=0 (1)
  integer :: do_ray = 0                          ! calculate/output Rayleigh refl=1, not=0 (0)
  real(wp) :: k2 = -1                            ! |K|^2, -1=use frequency dependent default (-1)
  ! namelist variables for COSP input related to lidar simulator
  integer :: lidar_ice_type = 0                  ! Ice particle shape in lidar calculations
                                                 ! (0=ice-spheres ; 1=ice-non-spherical) (0)
  integer, parameter :: overlap = 3              ! overlap type: 1=max, 2=rand, 3=max/rand (3)

  ! namelist variables for COSP input related to ISCCP simulator
  integer :: isccp_topheight = 1                 ! 1 = adjust top height using both a computed infrared
                                                 ! brightness temperature and the visible
                                                 ! optical depth to adjust cloud top pressure.
                                                 ! Note that this calculation is most appropriate to compare
                                                 ! to ISCCP data during sunlit hours.
                                                 ! 2 = do not adjust top height, that is cloud top pressure
                                                 ! is the actual cloud top pressure in the model
                                                 ! 3 = adjust top height using only the computed infrared
                                                 ! brightness temperature. Note that this calculation is most
                                                 ! appropriate to compare to ISCCP IR only algortihm (i.e.
                                                 ! you can compare to nighttime ISCCP data with this option) (1)
  integer :: isccp_topheight_direction = 2       ! direction for finding atmosphere pressure level with
                                                 ! interpolated temperature equal to the radiance
                                                 ! determined cloud-top temperature
                                                 ! 1 = find the *lowest* altitude (highest pressure) level
                                                 ! with interpolated temperature
                                                 ! equal to the radiance determined cloud-top temperature
                                                 ! 2 = find the *highest* altitude (lowest pressure) level
                                                 ! with interpolated temperature
                                                 ! equal to the radiance determined cloud-top temperature
                                                 ! ONLY APPLICABLE IF top_height EQUALS 1 or 3
                                                 ! 1 = default setting in COSP v1.1, matches all versions of
                                                 ! ISCCP simulator with versions numbers 3.5.1 and lower
                                                 ! 2 = default setting in COSP v1.3. default since V4.0 of ISCCP simulator

  ! ######################################################################################
  ! Other variables
  ! ######################################################################################
  logical,allocatable :: first_run_cosp(:)      !.true. if run_cosp has been populated (allocatable->begchunk:endchunk)
  logical,allocatable :: run_cosp(:,:)          !.true. if cosp should be run by column and
                                                !       chunk (allocatable->1:pcols,begchunk:endchunk)
  ! pbuf indices
  integer :: cld_idx, concld_idx, lsreffrain_idx, lsreffsnow_idx, cvreffliq_idx
  integer :: cvreffice_idx, dpcldliq_idx, dpcldice_idx
  integer :: shcldliq_idx, shcldice_idx, shcldliq1_idx, shcldice1_idx, dpflxprc_idx
  integer :: dpflxsnw_idx, shflxprc_idx, shflxsnw_idx, lsflxprc_idx, lsflxsnw_idx
  integer :: rei_idx, rel_idx
  integer :: crm_qc_idx, crm_qi_idx, crm_qpl_idx, crm_qpi_idx, &
             crm_qr_idx, crm_qs_idx, crm_qg_idx

  ! ######################################################################################
  ! Declarations specific to COSP2
  ! ######################################################################################
  type(radar_cfg)              :: rcfg_cloudsat    ! Radar configuration (Cloudsat)
  type(radar_cfg), allocatable :: rcfg_cs(:)       ! chunked version of rcfg_cloudsat
  type(size_distribution)              :: sd       ! Size distribution used by radar simulator
  type(size_distribution), allocatable :: sd_cs(:) ! chunked version of sd
  character(len=64) :: cloudsat_micro_scheme = 'MMF_v3.5_single_moment'

  integer,parameter :: &
       I_LSCLIQ = 1, & ! Large-scale (stratiform) liquid
       I_LSCICE = 2, & ! Large-scale (stratiform) ice
       I_LSRAIN = 3, & ! Large-scale (stratiform) rain
       I_LSSNOW = 4, & ! Large-scale (stratiform) snow
       I_CVCLIQ = 5, & ! Convective liquid
       I_CVCICE = 6, & ! Convective ice
       I_CVRAIN = 7, & ! Convective rain
       I_CVSNOW = 8, & ! Convective snow
       I_LSGRPL = 9    ! Large-scale (stratiform) groupel

  ! Stratiform and convective clouds in frac_out (scops output).
  integer, parameter :: &
       I_LSC = 1, & ! Large-scale clouds
       I_CVC = 2    ! Convective clouds

  ! Microphysical settings for the precipitation flux to mixing ratio conversion
  real(wp),parameter,dimension(nhydro) :: &
                 !  LSL     LSI         LSR         LSS       CVL     CVI        CVR          CVS          LSG
       N_ax    = (/-1._r8, -1._r8,     8.e6_r8,     3.e6_r8, -1._r8, -1._r8,     8.e6_r8,     3.e6_r8,     4.e6_r8/),&
       N_bx    = (/-1._r8, -1._r8,      0.0_r8,      0.0_r8, -1._r8, -1._r8,      0.0_r8,      0.0_r8,      0.0_r8/),&
       alpha_x = (/-1._r8, -1._r8,      0.0_r8,      0.0_r8, -1._r8, -1._r8,      0.0_r8,      0.0_r8,      0.0_r8/),&
       c_x     = (/-1._r8, -1._r8,    842.0_r8,     4.84_r8, -1._r8, -1._r8,    842.0_r8,     4.84_r8,     94.5_r8/),&
       d_x     = (/-1._r8, -1._r8,      0.8_r8,     0.25_r8, -1._r8, -1._r8,      0.8_r8,     0.25_r8,      0.5_r8/),&
       g_x     = (/-1._r8, -1._r8,      0.5_r8,      0.5_r8, -1._r8, -1._r8,      0.5_r8,      0.5_r8,      0.5_r8/),&
       a_x     = (/-1._r8, -1._r8,    524.0_r8,    52.36_r8, -1._r8, -1._r8,    524.0_r8,    52.36_r8,   209.44_r8/),&
       b_x     = (/-1._r8, -1._r8,      3.0_r8,      3.0_r8, -1._r8, -1._r8,      3.0_r8,      3.0_r8,      3.0_r8/),&
       gamma_1 = (/-1._r8, -1._r8, 17.83725_r8, 8.284701_r8, -1._r8, -1._r8, 17.83725_r8, 8.284701_r8, 11.63230_r8/),&
       gamma_2 = (/-1._r8, -1._r8,      6.0_r8,      6.0_r8, -1._r8, -1._r8,      6.0_r8,      6.0_r8,      6.0_r8/),&
       gamma_3 = (/-1._r8, -1._r8,      2.0_r8,      2.0_r8, -1._r8, -1._r8,      2.0_r8,      2.0_r8,      2.0_r8/),&
       gamma_4 = (/-1._r8, -1._r8,      6.0_r8,      6.0_r8, -1._r8, -1._r8,      6.0_r8,      6.0_r8,      6.0_r8/)

CONTAINS

  ! ######################################################################################
  ! SUBROUTINE setcosp2values
  ! ######################################################################################
  subroutine setcosp2values()
    use mod_cosp,             only: cosp_init
    use mod_cosp_config,      only: vgrid_zl, vgrid_zu, vgrid_z
    use mod_quickbeam_optics, only: hydro_class_init, quickbeam_optics_init
    use crmdims,              only: crm_nx_rad, crm_ny_rad

    ! Local
    logical :: ldouble=.false.
    logical :: lsingle=.true. ! Default is to use single moment
    integer :: i, k

    ! Initialize the distributional parameters for hydrometeors in radar simulator. In COSPv1.4, this was declared in
    ! cosp_defs.f.
    ! TODO: should this stuff go in cospsimulator_intr_init?
    if (cloudsat_micro_scheme == 'MMF_v3.5_two_moment') then
       ldouble = .true.
       lsingle = .false.
    endif
    call hydro_class_init(lsingle,ldouble,sd)
    call quickbeam_optics_init()

    ! DS2017: The setting up of the vertical grid for regridding the CALIPSO and Cloudsat products is
    !         now done in cosp_init, but these fields are stored in cosp_config.F90.
    !         Additionally all static fields used by the individual simulators are set up by calls
    !         to _init functions in cosp_init.
    call COSP_INIT(cosp_lisccp_sim,cosp_lmodis_sim,cosp_lmisr_sim, &
         cosp_lradar_sim,cosp_llidar_sim,cosp_lgrlidar_sim,cosp_latlid_sim, &
         cosp_lparasol_sim,cosp_lrttov_sim, &
         radar_freq,k2,use_gas_abs,do_ray,isccp_topheight,isccp_topheight_direction,    &
         surface_radar,rcfg_cloudsat,use_vgrid,csat_vgrid,Nlr,pver,            &
         cloudsat_micro_scheme)

    ! DJS2017: In COSP2, most of the bin boundaries, centers, and edges are declared in src/cosp_config.F90.
    !          Above I just assign them accordingly in the USE statement. Other bin bounds needed by CAM
    !          are calculated here.
    ! Allocate
    allocate(htlim_cosp(2,Nlvgrid),htmid_cosp(Nlvgrid),scol_cosp(nSubcol))

    ! DJS2017: Just pull from cosp_config
    htmid_cosp      = vgrid_z
    htlim_cosp(1,:) = vgrid_zu
    htlim_cosp(2,:) = vgrid_zl

    scol_cosp(:) = (/(k,k=1,nSubcol)/)

  end subroutine setcosp2values

  ! ######################################################################################
  ! SUBROUTINE cospsimulator_intr_readnl
  !
  ! PURPOSE: to read namelist variables and run setcospvalues subroutine.note: cldfrc_readnl
  ! is a good template in cloud_fraction.F90. Make sure that this routine is reading in a
  ! namelist. models/atm/cam/bld/build-namelist is the perl script to check.
  ! ######################################################################################
  subroutine cospsimulator_intr_readnl(nlfile)
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
#ifdef SPMD
    use mpishorthand,    only: mpicom, mpilog, mpiint, mpichar
#endif
    use crmdims, only: crm_nx_rad, crm_ny_rad

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input  (nlfile=atm_in)

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'cospsimulator_intr_readnl'

    ! this list should include any variable that you might want to include in the namelist
    ! philosophy is to not include COSP output flags but just important COSP settings and cfmip controls.
    namelist /cospsimulator_nl/ docosp, cosp_active, cosp_amwg, cosp_atrainorbitdata, cosp_cfmip_3hr, cosp_cfmip_da, &
         cosp_cfmip_mon, cosp_cfmip_off, cosp_histfile_num, cosp_histfile_aux, cosp_histfile_aux_num, cosp_isccp, cosp_lfrac_out, &
         cosp_lite, cosp_lradar_sim, cosp_llidar_sim, cosp_lisccp_sim,  cosp_lmisr_sim, cosp_lmodis_sim, cosp_ncolumns, &
         cosp_nradsteps, cosp_passive, cosp_sample_atrain

    ! read in the namelist
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )  !! presumably opens the namelist file "nlfile"
       !! position the file to write to the cospsimulator portion of the cam_in namelist
       call find_group_name(unitn, 'cospsimulator_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, cospsimulator_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast(docosp,               1,  mpilog, 0, mpicom)
    !call mpibcast(cosp_atrainorbitdata, len(cosp_atrainorbitdata), mpichar, 0, mpicom)
    !call mpibcast(cosp_sample_atrain,   1,  mpilog, 0, mpicom)
    call mpibcast(cosp_amwg,            1,  mpilog, 0, mpicom)
    call mpibcast(cosp_lite,            1,  mpilog, 0, mpicom)
    call mpibcast(cosp_passive,         1,  mpilog, 0, mpicom)
    call mpibcast(cosp_active,          1,  mpilog, 0, mpicom)
    call mpibcast(cosp_isccp,           1,  mpilog, 0, mpicom)
    call mpibcast(cosp_cfmip_3hr,       1,  mpilog, 0, mpicom)
    call mpibcast(cosp_cfmip_da,        1,  mpilog, 0, mpicom)
    call mpibcast(cosp_cfmip_mon,       1,  mpilog, 0, mpicom)
    call mpibcast(cosp_cfmip_off,       1,  mpilog, 0, mpicom)
    call mpibcast(cosp_lfrac_out,       1,  mpilog, 0, mpicom)
    call mpibcast(cosp_lradar_sim,      1,  mpilog, 0, mpicom)
    call mpibcast(cosp_llidar_sim,      1,  mpilog, 0, mpicom)
    call mpibcast(cosp_lisccp_sim,      1,  mpilog, 0, mpicom)
    call mpibcast(cosp_lmisr_sim,       1,  mpilog, 0, mpicom)
    call mpibcast(cosp_lmodis_sim,      1,  mpilog, 0, mpicom)
    call mpibcast(cosp_ncolumns,        1,  mpiint, 0, mpicom)
    call mpibcast(cosp_histfile_num,    1,  mpiint, 0, mpicom)
    call mpibcast(cosp_histfile_aux_num,1,  mpiint, 0, mpicom)
    call mpibcast(cosp_histfile_aux,    1,  mpilog, 0, mpicom)
    call mpibcast(cosp_nradsteps,       1,  mpiint, 0, mpicom)
#endif

    ! reset COSP namelist variables based on input from cam namelist variables
    !DJS2017: The parasol simulator is now separate from the lidar simulator. To maintain consistency, just
    !         mirror whatever the lidar simulator is doing
    if (cosp_cfmip_3hr) then
       cosp_lradar_sim   = .true.
       cosp_llidar_sim   = .true.
       cosp_lparasol_sim = .true.
       cosp_lisccp_sim   = .true.
    end if
    if (cosp_cfmip_da) then
       cosp_llidar_sim = .true.
       cosp_lparasol_sim = .true.
       cosp_lisccp_sim = .true.
    end if
    if (cosp_cfmip_off) then
       cosp_lradar_sim = .true.
       cosp_llidar_sim = .true.
       cosp_lparasol_sim = .true.
       cosp_lisccp_sim = .true.
    end if
    if (cosp_cfmip_mon) then
       cosp_llidar_sim = .true.
       cosp_lparasol_sim = .true.
       cosp_lisccp_sim = .true.
    end if

    if (cosp_llidar_sim) then
       cosp_lparasol_sim = .true.
    end if

    if (cosp_histfile_aux .and. cosp_histfile_aux_num == -1) then
       cosp_histfile_aux_num = cosp_histfile_num
    end if

    if (cosp_lite) then
       cosp_llidar_sim = .true.
       cosp_lparasol_sim = .true.
       cosp_lisccp_sim = .true.
       cosp_lmisr_sim = .true.
       cosp_lmodis_sim = .true.
       cosp_ncolumns = 10
    end if

    if (cosp_passive) then
       cosp_lisccp_sim = .true.
       cosp_lmisr_sim = .true.
       cosp_lmodis_sim = .true.
       cosp_ncolumns = 10
    end if

    if (cosp_active) then
       cosp_lradar_sim = .true.
       cosp_llidar_sim = .true.
       cosp_lparasol_sim = .true.
       cosp_ncolumns = 10
    end if

    if (cosp_isccp) then
       cosp_lisccp_sim = .true.
       cosp_ncolumns = 10
    end if

    ! If no simulators are turned on at all and docosp is, set cosp_amwg = .true.
    if((docosp) .and. (.not.cosp_lradar_sim) .and. (.not.cosp_llidar_sim) .and. (.not.cosp_lisccp_sim) .and. &
         (.not.cosp_lmisr_sim) .and. (.not.cosp_lmodis_sim)) then
       cosp_amwg = .true.
    end if
    if (cosp_amwg) then
       cosp_lradar_sim = .true.
       cosp_llidar_sim = .true.
       cosp_lparasol_sim = .true.
       cosp_lisccp_sim = .true.
       cosp_lmisr_sim = .true.
       cosp_lmodis_sim = .true.
       cosp_ncolumns = 10
    end if

    ! Reset COSP namelist variables based on input from cam namelist variables
    if (cosp_ncolumns .ne. nSubcol) then
       nSubcol = cosp_ncolumns
    end if

    if (masterproc) then
       if (docosp) then
          write(iulog,*)'COSP configuration:'
          write(iulog,*)'  Number of COSP subcolumns                = ', cosp_ncolumns
          write(iulog,*)'  Frequency at which cosp is called        = ', cosp_nradsteps
          write(iulog,*)'  Enable radar simulator                   = ', cosp_lradar_sim
          write(iulog,*)'  Enable lidar simulator                   = ', cosp_llidar_sim
          write(iulog,*)'  Enable ISCCP simulator                   = ', cosp_lisccp_sim
          write(iulog,*)'  Enable MISR simulator                    = ', cosp_lmisr_sim
          write(iulog,*)'  Enable MODIS simulator                   = ', cosp_lmodis_sim
          write(iulog,*)'  RADAR_SIM microphysics scheme            = ', trim(cloudsat_micro_scheme)
          write(iulog,*)'  Write COSP output to history file        = ', cosp_histfile_num
          write(iulog,*)'  Write COSP input fields                  = ', cosp_histfile_aux
          write(iulog,*)'  Write COSP input fields to history file  = ', cosp_histfile_aux_num
          write(iulog,*)'  Write COSP subcolumn fields              = ', cosp_lfrac_out
       else
          write(iulog,*)'COSP not enabled'
       end if
    end if
  end subroutine cospsimulator_intr_readnl

  ! ######################################################################################
  ! SUBROUTINE cospsimulator_intr_register
  ! ######################################################################################
  subroutine cospsimulator_intr_register()

    use cam_history_support, only: add_hist_coord
    use phys_control, only: phys_getopts
    use crmdims, only: crm_nx_rad, crm_ny_rad
    logical :: use_SPCAM

    ! Reset number subcolumns if using SP/MMF
    call phys_getopts(use_SPCAM_out=use_SPCAM)
    if (use_SPCAM) then
       nSubcol = crm_nx_rad * crm_ny_rad
    end if

    ! Set vertical coordinate, subcolumn, and calculation frequency cosp options based on namelist inputs
    call setcosp2values()

    ! register non-standard variable dimensions
    if (cosp_lisccp_sim .or. cosp_lmodis_sim) then
       call add_hist_coord('cosp_prs', nprs_cosp, 'COSP Mean ISCCP pressure',  &
            'hPa', pres_binCenters, bounds_name='cosp_prs_bnds', bounds=pres_binEdges)
    end if

    if (cosp_lisccp_sim .or. cosp_lmisr_sim) then
       call add_hist_coord('cosp_tau', ntau_cosp,                 &
            'COSP Mean ISCCP optical depth', '1', tau_binCenters, &
            bounds_name='cosp_tau_bnds', bounds=tau_binEdges      )
    end if

    if (cosp_lisccp_sim .or. cosp_llidar_sim .or. cosp_lradar_sim .or. cosp_lmisr_sim) then
       call add_hist_coord('cosp_scol', nSubcol, 'COSP subcolumn',  values=scol_cosp)
    end if

    if (cosp_llidar_sim .or. cosp_lradar_sim) then
       call add_hist_coord('cosp_ht', Nlvgrid,                             &
            'COSP Mean Height for lidar and radar simulator outputs', 'm', &
            htmid_cosp, bounds_name='cosp_ht_bnds', bounds=htlim_cosp,     &
            vertical_coord=.true.)
    end if

    if (cosp_llidar_sim) then
       call add_hist_coord('cosp_sr', nsr_cosp,                                &
            'COSP Mean Scattering Ratio for lidar simulator CFAD output', '1', &
            calipso_binCenters, bounds_name='cosp_sr_bnds', bounds=calipso_binEdges)
    end if

    if (cosp_llidar_sim) then
       call add_hist_coord('cosp_sza', nsza_cosp, 'COSP Parasol SZA', 'degrees', parasol_sza)
    end if

    if (cosp_lradar_sim) then
       call add_hist_coord('cosp_dbze', ndbze_cosp,                  &
            'COSP Mean dBZe for radar simulator CFAD output', 'dBZ', &
            cloudsat_binCenters, bounds_name='cosp_dbze_bnds',       &
            bounds=cloudsat_binEdges)
    end if

    if (cosp_lmisr_sim) then
       call add_hist_coord('cosp_htmisr', nhtmisr_cosp, 'COSP MISR height', &
            'km', misr_histHgtCenters,                                      &
            bounds_name='cosp_htmisr_bnds', bounds=misr_histHgtEdges)
    end if

    if (cosp_lmodis_sim) then
       call add_hist_coord('cosp_tau_modis', ntau_cosp_modis,               &
            'COSP Mean MODIS optical depth', '1', tau_binCenters,           &
            bounds_name='cosp_tau_modis_bnds', bounds=tau_binEdges)
       call add_hist_coord('cosp_reffice',numMODISReffIceBins,                       &
            'COSP Mean MODIS effective radius (ice)', 'microns', reffICE_binCenters, &
            bounds_name='cosp_reffice_bnds',bounds=reffICE_binEdges)
       call add_hist_coord('cosp_reffliq',numMODISReffLiqBins,                          &
            'COSP Mean MODIS effective radius (liquid)', 'microns', reffLIQ_binCenters, &
            bounds_name='cosp_reffliq_bnds',bounds=reffLIQ_binEdges)
    end if

  end subroutine cospsimulator_intr_register

  ! ######################################################################################
  ! SUBROUTINE cospsimulator_intr_init
  ! ######################################################################################
  subroutine cospsimulator_intr_init()

    use cam_history,         only: addfld, add_default, horiz_only
#ifdef SPMD
    use mpishorthand,        only : mpir8, mpiint, mpicom
#endif
    use netcdf,              only : nf90_open, nf90_inq_varid, nf90_get_var, nf90_close, nf90_nowrite
    use error_messages,      only : handle_ncerr, alloc_err
    use physics_buffer,  only: pbuf_get_index
    use phys_control, only : phys_getopts
    use mod_cosp_config,  only : R_UNDEF

    integer :: ncid,latid,lonid,did,hrid,minid,secid, istat
    integer :: i
    logical :: use_SPCAM
    character(len=16) :: SPCAM_microp_scheme

    ! ISCCP OUTPUTS
    if (cosp_lisccp_sim) then
       call addfld('FISCCP1_COSP',(/'cosp_tau','cosp_prs'/),'A','percent', &
            'Grid-box fraction covered by each ISCCP D level cloud type',&
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDTOT_ISCCP', horiz_only,'A','percent', &
            'Total Cloud Fraction Calculated by the ISCCP Simulator ',flag_xyfill=.true., fill_value=R_UNDEF)
       ! Per CFMIP request - weight by ISCCP Total Cloud Fraction (divide by CLDTOT_ISSCP in history file to get weighted average)
       call addfld('MEANCLDALB_ISCCP',horiz_only,'A','1','Mean cloud albedo*CLDTOT_ISCCP',flag_xyfill=.true., fill_value=R_UNDEF)
       ! Per CFMIP request - weight by ISCCP Total Cloud Fraction (divide by CLDTOT_ISSCP in history file to get weighted average)
       call addfld('MEANPTOP_ISCCP',horiz_only,'A','Pa','Mean cloud top pressure*CLDTOT_ISCCP',flag_xyfill=.true., &
            fill_value=R_UNDEF)
       ! For averaging, weight by ISCCP Total Cloud Fraction (divide by CLDTOT_ISSCP in history file to get weighted average)
       call addfld ('MEANTAU_ISCCP',horiz_only,'A','1','Mean optical thickness*CLDTOT_ISCCP',flag_xyfill=.true., &
            fill_value=R_UNDEF)
       call addfld ('MEANTB_ISCCP',horiz_only,'A','K','Mean Infrared Tb from ISCCP simulator',flag_xyfill=.true., &
            fill_value=R_UNDEF)
       call addfld ('MEANTBCLR_ISCCP',horiz_only,'A','K','Mean Clear-sky Infrared Tb from ISCCP simulator',   &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('TAU_ISCCP',(/'cosp_scol'/),'I','1','Optical Depth in each Subcolumn',flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CLDPTOP_ISCCP',(/'cosp_scol'/),'I','Pa','Cloud Top Pressure in each Subcolumn',  &
            flag_xyfill=.true., fill_value=R_UNDEF)

       ! add_default calls for CFMIP experiments or else all fields are added to history file
       !     except those with sub-column dimension
       ! add cfmip-requested variables to two separate cam history files
       if (cosp_cfmip_mon.or.cosp_cfmip_da) then
          if (cosp_cfmip_da) then
             call add_default ('FISCCP1_COSP',2,' ')
             call add_default ('CLDTOT_ISCCP',2,' ')
             call add_default ('MEANCLDALB_ISCCP',2,' ')
             call add_default ('MEANPTOP_ISCCP',2,' ')
          end if
          if (cosp_cfmip_mon) then
             call add_default ('FISCCP1_COSP',1,' ')
             call add_default ('CLDTOT_ISCCP',1,' ')
             call add_default ('MEANCLDALB_ISCCP',1,' ')
             call add_default ('MEANPTOP_ISCCP',1,' ')
          end if
       else
          ! add all isccp outputs to the history file specified by the CAM namelist variable cosp_histfile_num
          call add_default ('FISCCP1_COSP',cosp_histfile_num,' ')
          call add_default ('CLDTOT_ISCCP',cosp_histfile_num,' ')
          call add_default ('MEANCLDALB_ISCCP',cosp_histfile_num,' ')
          call add_default ('MEANPTOP_ISCCP',cosp_histfile_num,' ')
          call add_default ('MEANTAU_ISCCP',cosp_histfile_num,' ')
          call add_default ('MEANTB_ISCCP',cosp_histfile_num,' ')
          call add_default ('MEANTBCLR_ISCCP',cosp_histfile_num,' ')
       end if
    end if

    ! LIDAR SIMULATOR OUTPUTS
    if (cosp_llidar_sim) then
       call addfld('CLDLOW_CAL',horiz_only,'A','percent','Lidar Low-level Cloud Fraction',flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDMED_CAL',horiz_only,'A','percent','Lidar Mid-level Cloud Fraction',flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDHGH_CAL',horiz_only,'A','percent','Lidar High-level Cloud Fraction',flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDTOT_CAL',horiz_only,'A','percent','Lidar Total Cloud Fraction',flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLD_CAL',(/'cosp_ht'/),'A','percent','Lidar Cloud Fraction (532 nm)', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CFAD_SR532_CAL',(/'cosp_sr','cosp_ht'/),'A','fraction',                                    &
            'Lidar Scattering Ratio CFAD (532 nm)',                                                    &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('MOL532_CAL',(/'lev'/),'A','m-1sr-1','Lidar Molecular Backscatter (532 nm) ',              &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('ATB532_CAL',(/'cosp_scol','lev      '/),'I','no_unit_log10(x)',                           &
            'Lidar Attenuated Total Backscatter (532 nm) in each Subcolumn',                        &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLD_CAL_LIQ', (/'cosp_ht'/), 'A','percent', 'Lidar Liquid Cloud Fraction',                 &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLD_CAL_ICE', (/'cosp_ht'/), 'A','percent', 'Lidar Ice Cloud Fraction',                    &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLD_CAL_UN', (/'cosp_ht'/),'A','percent', 'Lidar Undefined-Phase Cloud Fraction',          &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLD_CAL_TMP', (/'cosp_ht'/), 'A','percent', 'NOT SURE WHAT THIS IS Cloud Fraction',        &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLD_CAL_TMPLIQ', (/'cosp_ht'/), 'A','percent', 'NOT SURE WHAT THIS IS Cloud Fraction',     &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLD_CAL_TMPICE', (/'cosp_ht'/), 'A','percent', 'NOT SURE WHAT THIS IS Cloud Fraction',     &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLD_CAL_TMPUN', (/'cosp_ht'/), 'A','percent', 'NOT SURE WHAT THIS IS Cloud Fraction',      &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDTOT_CAL_ICE', horiz_only,'A','percent','Lidar Total Ice Cloud Fraction',    &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDTOT_CAL_LIQ', horiz_only,'A','percent','Lidar Total Liquid Cloud Fraction', &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDTOT_CAL_UN',horiz_only,'A','percent','Lidar Total Undefined-Phase Cloud Fraction',      &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDHGH_CAL_ICE',horiz_only,'A','percent','Lidar High-level Ice Cloud Fraction',            &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDHGH_CAL_LIQ',horiz_only,'A','percent','Lidar High-level Liquid Cloud Fraction',         &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDHGH_CAL_UN',horiz_only,'A','percent','Lidar High-level Undefined-Phase Cloud Fraction', &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDMED_CAL_ICE',horiz_only,'A','percent','Lidar Mid-level Ice Cloud Fraction',             &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDMED_CAL_LIQ',horiz_only,'A','percent','Lidar Mid-level Liquid Cloud Fraction',          &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDMED_CAL_UN',horiz_only,'A','percent','Lidar Mid-level Undefined-Phase Cloud Fraction',  &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDLOW_CAL_ICE',horiz_only,'A','percent','Lidar Low-level Ice Cloud Fraction',             &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDLOW_CAL_LIQ',horiz_only,'A','percent','Lidar Low-level Liquid Cloud Fraction',          &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDLOW_CAL_UN',horiz_only,'A','percent','Lidar Low-level Undefined-Phase Cloud Fraction',  &
            flag_xyfill=.true., fill_value=R_UNDEF)

       ! add_default calls for CFMIP experiments or else all fields are added to history file
       !     except those with sub-column dimension/experimental variables
       if (cosp_cfmip_mon .or. cosp_cfmip_off .or. cosp_cfmip_da .or. cosp_cfmip_3hr) then
          if (cosp_cfmip_da) then
             call add_default ('CLDLOW_CAL',2,' ')
             call add_default ('CLDMED_CAL',2,' ')
             call add_default ('CLDHGH_CAL',2,' ')
             call add_default ('CLDTOT_CAL',2,' ')
             call add_default ('CLD_CAL',2,' ')
          end if
          if (cosp_cfmip_mon.or.cosp_cfmip_off) then
             call add_default ('CLDLOW_CAL',1,' ')
             call add_default ('CLDMED_CAL',1,' ')
             call add_default ('CLDHGH_CAL',1,' ')
             call add_default ('CLDTOT_CAL',1,' ')
             call add_default ('CLD_CAL',1,' ')
          end if
          if (cosp_cfmip_3hr) then
             call add_default ('CFAD_SR532_CAL',3,' ')
             call add_default ('CLDLOW_CAL',3,' ')
             call add_default ('CLDMED_CAL',3,' ')
             call add_default ('CLDHGH_CAL',3,' ')
             call add_default ('CLDTOT_CAL',3,' ')
             call add_default ('CLD_CAL',3,' ')
          end if
          if (cosp_cfmip_off) then
             call add_default ('CFAD_SR532_CAL',1,' ')
          end if
       else
          ! add all lidar outputs to the history file specified by the CAM namelist variable cosp_histfile_num
          call add_default ('CLDLOW_CAL',cosp_histfile_num,' ')
          call add_default ('CLDMED_CAL',cosp_histfile_num,' ')
          call add_default ('CLDHGH_CAL',cosp_histfile_num,' ')
          call add_default ('CLDTOT_CAL',cosp_histfile_num,' ')
          call add_default ('CLD_CAL',cosp_histfile_num,' ')
          call add_default ('CFAD_SR532_CAL',cosp_histfile_num,' ')
          call add_default ('CLD_CAL_LIQ',cosp_histfile_num,' ')
          call add_default ('CLD_CAL_ICE',cosp_histfile_num,' ')
          call add_default ('CLD_CAL_UN',cosp_histfile_num,' ')
          call add_default ('CLDTOT_CAL_ICE',cosp_histfile_num,' ')
          call add_default ('CLDTOT_CAL_LIQ',cosp_histfile_num,' ')
          call add_default ('CLDTOT_CAL_UN',cosp_histfile_num,' ')
          call add_default ('CLDHGH_CAL_ICE',cosp_histfile_num,' ')
          call add_default ('CLDHGH_CAL_LIQ',cosp_histfile_num,' ')
          call add_default ('CLDHGH_CAL_UN',cosp_histfile_num,' ')
          call add_default ('CLDMED_CAL_ICE',cosp_histfile_num,' ')
          call add_default ('CLDMED_CAL_LIQ',cosp_histfile_num,' ')
          call add_default ('CLDMED_CAL_UN',cosp_histfile_num,' ')
          call add_default ('CLDLOW_CAL_ICE',cosp_histfile_num,' ')
          call add_default ('CLDLOW_CAL_LIQ',cosp_histfile_num,' ')
          call add_default ('CLDLOW_CAL_UN',cosp_histfile_num,' ')

          if ((.not.cosp_amwg) .and. (.not.cosp_lite) .and. (.not.cosp_passive) .and. (.not.cosp_active) &
             .and. (.not.cosp_isccp)) then
             call add_default ('MOL532_CAL',cosp_histfile_num,' ')
          end if
       end if
    end if

    if (cosp_lparasol_sim) then
       call addfld ('RFL_PARASOL',(/'cosp_sza'/),'A','fraction','PARASOL-like mono-directional reflectance ',  &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('RFL_PARASOL_PIX',(/'cosp_scol','cosp_sza '/),'I','fraction','PARASOL-like mono-directional reflectance ',  &
            flag_xyfill=.true., fill_value=R_UNDEF)

       ! add_default calls for CFMIP experiments or else all fields are added to history file
       !     except those with sub-column dimension/experimental variables
       if (cosp_cfmip_mon .or. cosp_cfmip_off .or. cosp_cfmip_da .or. cosp_cfmip_3hr) then
          if (cosp_cfmip_da) then
             call add_default ('RFL_PARASOL',2,' ')
          end if
          if (cosp_cfmip_mon.or.cosp_cfmip_off) then
             call add_default ('RFL_PARASOL',1,' ')
          end if
          if (cosp_cfmip_3hr) then
             call add_default ('RFL_PARASOL',3,' ')
          end if
       else
          call add_default ('RFL_PARASOL',cosp_histfile_num,' ')
       end if
    end if
    ! RADAR SIMULATOR OUTPUTS
    if (cosp_lradar_sim) then

       allocate(sd_cs(begchunk:endchunk), rcfg_cs(begchunk:endchunk))
       do i = begchunk, endchunk
          sd_cs(i)   = sd
          rcfg_cs(i) = rcfg_cloudsat
       end do

       call addfld('CFAD_DBZE94_CS',(/'cosp_dbze','cosp_ht  '/),'A','fraction',&
            'Radar Reflectivity Factor CFAD (94 GHz)',&
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CLD_CAL_NOTCS',(/'cosp_ht'/),'A','percent','Cloud occurrence seen by CALIPSO but not CloudSat ',   &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CLDTOT_CALCS',horiz_only,'A','percent',' Lidar and Radar Total Cloud Fraction ',flag_xyfill=.true., &
            fill_value=R_UNDEF)
       call addfld ('CLDTOT_CS',horiz_only,'A','percent',' Radar total cloud amount ',flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CLDTOT_CS2',horiz_only,'A','percent', &
            ' Radar total cloud amount without the data for the first kilometer above surface ', &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('DBZE_CS',(/'cosp_scol','lev      '/),'I','dBZe',' Radar dBZe (94 GHz) in each Subcolumn',&
            flag_xyfill=.true., fill_value=R_UNDEF)

       ! add_default calls for CFMIP experiments or else all fields are added to history file except those with sub-column dimension
       if (cosp_cfmip_off.or.cosp_cfmip_3hr) then
          if (cosp_cfmip_3hr) then
             call add_default ('CFAD_DBZE94_CS',3,' ')
             call add_default ('CLD_CAL_NOTCS',3,' ')
          end if
          if (cosp_cfmip_off) then
             call add_default ('CFAD_DBZE94_CS',1,' ')
             call add_default ('CLD_CAL_NOTCS',1,' ')
          end if
       else
          ! add all radar outputs to the history file specified by the CAM namelist variable cosp_histfile_num
          call add_default ('CFAD_DBZE94_CS',cosp_histfile_num,' ')
          call add_default ('CLD_CAL_NOTCS',cosp_histfile_num,' ')
          call add_default ('CLDTOT_CALCS',cosp_histfile_num,' ')
          call add_default ('CLDTOT_CS',cosp_histfile_num,' ')
          call add_default ('CLDTOT_CS2',cosp_histfile_num,' ')
       end if
    end if

    ! MISR SIMULATOR OUTPUTS
    if (cosp_lmisr_sim) then
       call addfld ('CLD_MISR',(/'cosp_tau   ','cosp_htmisr'/),'A','percent','Cloud Fraction from MISR Simulator',  &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call add_default ('CLD_MISR',cosp_histfile_num,' ')
    end if

    ! MODIS OUTPUT
    if (cosp_lmodis_sim) then
       call addfld ('CLTMODIS',horiz_only,'A','%','MODIS Total Cloud Fraction',flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CLWMODIS',horiz_only,'A','%','MODIS Liquid Cloud Fraction',flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CLIMODIS',horiz_only,'A','%','MODIS Ice Cloud Fraction',flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CLHMODIS',horiz_only,'A','%','MODIS High Level Cloud Fraction',flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CLMMODIS',horiz_only,'A','%','MODIS Mid Level Cloud Fraction',flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CLLMODIS',horiz_only,'A','%','MODIS Low Level Cloud Fraction',flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('TAUTMODIS',horiz_only,'A','1','MODIS Total Cloud Optical Thickness*CLTMODIS',                  &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('TAUWMODIS',horiz_only,'A','1','MODIS Liquid Cloud Optical Thickness*CLWMODIS',                 &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('TAUIMODIS',horiz_only,'A','1','MODIS Ice Cloud Optical Thickness*CLIMODIS',                    &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('TAUTLOGMODIS',horiz_only,'A','1','MODIS Total Cloud Optical Thickness (Log10 Mean)*CLTMODIS',  &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('TAUWLOGMODIS',horiz_only,'A','1','MODIS Liquid Cloud Optical Thickness (Log10 Mean)*CLWMODIS', &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('TAUILOGMODIS',horiz_only,'A','1','MODIS Ice Cloud Optical Thickness (Log10 Mean)*CLIMODIS',    &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('REFFCLWMODIS',horiz_only,'A','m','MODIS Liquid Cloud Particle Size*CLWMODIS',                  &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('REFFCLIMODIS',horiz_only,'A','m','MODIS Ice Cloud Particle Size*CLIMODIS',                     &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('PCTMODIS',horiz_only,'A','Pa','MODIS Cloud Top Pressure*CLTMODIS',flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('LWPMODIS',horiz_only,'A','kg m-2','MODIS Cloud Liquid Water Path*CLWMODIS',                    &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('IWPMODIS',horiz_only,'A','kg m-2','MODIS Cloud Ice Water Path*CLIMODIS',flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CLMODIS',(/'cosp_tau_modis','cosp_prs      '/),'A','%','MODIS Cloud Area Fraction',            &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CLRIMODIS',(/'cosp_tau_modis','cosp_reffice  '/),'A','%','MODIS Cloud Area Fraction',            &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CLRLMODIS',(/'cosp_tau_modis','cosp_reffliq  '/),'A','%','MODIS Cloud Area Fraction',            &
            flag_xyfill=.true., fill_value=R_UNDEF)

       ! add MODIS output to history file specified by the CAM namelist variable cosp_histfile_num
       call add_default ('CLTMODIS',cosp_histfile_num,' ')
       call add_default ('CLWMODIS',cosp_histfile_num,' ')
       call add_default ('CLIMODIS',cosp_histfile_num,' ')
       call add_default ('CLHMODIS',cosp_histfile_num,' ')
       call add_default ('CLMMODIS',cosp_histfile_num,' ')
       call add_default ('CLLMODIS',cosp_histfile_num,' ')
       call add_default ('TAUTMODIS',cosp_histfile_num,' ')
       call add_default ('TAUWMODIS',cosp_histfile_num,' ')
       call add_default ('TAUIMODIS',cosp_histfile_num,' ')
       call add_default ('TAUTLOGMODIS',cosp_histfile_num,' ')
       call add_default ('TAUWLOGMODIS',cosp_histfile_num,' ')
       call add_default ('TAUILOGMODIS',cosp_histfile_num,' ')
       call add_default ('REFFCLWMODIS',cosp_histfile_num,' ')
       call add_default ('REFFCLIMODIS',cosp_histfile_num,' ')
       call add_default ('PCTMODIS',cosp_histfile_num,' ')
       call add_default ('LWPMODIS',cosp_histfile_num,' ')
       call add_default ('IWPMODIS',cosp_histfile_num,' ')
       call add_default ('CLMODIS',cosp_histfile_num,' ')
       call add_default ('CLRIMODIS',cosp_histfile_num,' ')
       call add_default ('CLRLMODIS',cosp_histfile_num,' ')
    end if

    ! SUB-COLUMN OUTPUT
    if (cosp_lfrac_out) then
       call addfld(                                             &
          'SCOPS_OUT', (/'cosp_scol','lev      '/), 'I',        &
          '0 = no cloud, 1 = strat, 2 = conv',                  &
          'SCOPS Subcolumn output',                             &
          flag_xyfill=.true., fill_value=R_UNDEF                &
       )
       call addfld(                                             &
          'COSP_PREC_FRAC', (/'cosp_scol','lev      '/), 'I',   &
          '0 = no prec, 1 = strat, 2 = conv, 3 = both',         &
          'PREC_SCOPS subcolumn output',                        &
          flag_xyfill=.true., fill_value=R_UNDEF                &
       )
       call add_default('SCOPS_OUT', cosp_histfile_num, ' ')
       call add_default('COSP_PREC_FRAC', cosp_histfile_num, ' ')
       if (cosp_lisccp_sim) then
          call add_default('TAU_ISCCP', cosp_histfile_num, ' ')
          call add_default('CLDPTOP_ISCCP', cosp_histfile_num, ' ')
       end if
       if (cosp_llidar_sim) then
          call add_default('ATB532_CAL', cosp_histfile_num, ' ')
       end if
       if (cosp_lradar_sim) then
          call add_default('DBZE_CS', cosp_histfile_num, ' ')
       end if
       if (cosp_lparasol_sim) then
          call add_default('RFL_PARASOL_PIX', cosp_histfile_num, ' ')
       end if
    end if

    ! ADDFLD, ADD_DEFAULT, OUTFLD CALLS FOR COSP OUTPUTS IF RUNNING COSP OFF-LINE
    ! Note: A suggestion was to add all of the CAM variables needed to add to make it possible to run COSP off-line
    ! These fields are available and can be called from the namelist though.  Here, when the cosp_histfile_aux is true
    ! all of the inputs are saved on the cam history file.  This is good de-bugging functionality we should maintain.
    if (cosp_histfile_aux) then
       call addfld ('PS_COSP',         horiz_only,            'I','Pa',     'PS_COSP',                            &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('TS_COSP',         horiz_only,            'I','K',      'TS_COSP',                            &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('P_COSP',          (/            'lev'/), 'I','Pa',     'P_COSP',                             &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('PH_COSP',         (/            'lev'/), 'I','Pa',     'PH_COSP',                            &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('ZLEV_COSP',       (/            'lev'/), 'I','m',      'ZLEV_COSP',                          &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('ZLEV_HALF_COSP',  (/            'lev'/), 'I','m',      'ZLEV_HALF_COSP',                     &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('T_COSP',          (/            'lev'/), 'I','K',      'T_COSP',                             &
            flag_xyfill=.true., fill_value=R_UNDEF)
       ! TODO: What actually gets output for this variable is SPECIFIC HUMIDITY;
       ! change the name here!
       call addfld ('RH_COSP',         (/            'lev'/), 'I','percent','RH_COSP',                            &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('TAU_067',         (/'cosp_scol','lev      '/), 'I','1',      'Subcolumn 0.67micron optical depth', &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('EMISS_11',        (/'cosp_scol','lev      '/), 'I','1',      'Subcolumn 11micron emissivity',      &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('MODIS_fracliq',   (/'cosp_scol','lev      '/), 'I','1',      'Fraction of tau from liquid water',  &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('MODIS_asym',      (/'cosp_scol','lev      '/), 'I','1',      'Assymetry parameter (MODIS)',        &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('MODIS_ssa',       (/'cosp_scol','lev      '/), 'I','1',      'Single-scattering albedo (MODIS)',   &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CAL_betatot',     (/'cosp_scol','lev      '/), 'I','1',      'Backscatter coefficient (CALIPSO)',  &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CAL_betatot_ice', (/'cosp_scol','lev      '/), 'I','1',      'Backscatter coefficient (CALIPSO)',  &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CAL_betatot_liq', (/'cosp_scol','lev      '/), 'I','1',      'Backscatter coefficient (CALIPSO)',  &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CAL_tautot',      (/'cosp_scol','lev      '/), 'I','1',      'Vertically integrated ptical-depth (CALIPSO)', &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CAL_tautot_ice',  (/'cosp_scol','lev      '/), 'I','1',      'Vertically integrated ptical-depth (CALIPSO)', &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CAL_tautot_liq',  (/'cosp_scol','lev      '/), 'I','1',      'Vertically integrated ptical-depth (CALIPSO)', &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CS_z_vol',        (/'cosp_scol','lev      '/), 'I','1',      'Effective reflectivity factor (CLOUDSAT)',     &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CS_kr_vol',       (/'cosp_scol','lev      '/), 'I','1',      'Attenuation coefficient (hydro) (CLOUDSAT)',   &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CS_g_vol',        (/'cosp_scol','lev      '/), 'I','1',      'Attenuation coefficient (gases) (CLOUDSAT)',   &
            flag_xyfill=.true., fill_value=R_UNDEF)

       call add_default ('PS_COSP',         cosp_histfile_aux_num,' ')
       call add_default ('TS_COSP',         cosp_histfile_aux_num,' ')
       call add_default ('P_COSP',          cosp_histfile_aux_num,' ')
       call add_default ('PH_COSP',         cosp_histfile_aux_num,' ')
       call add_default ('ZLEV_COSP',       cosp_histfile_aux_num,' ')
       call add_default ('ZLEV_HALF_COSP',  cosp_histfile_aux_num,' ')
       call add_default ('T_COSP',          cosp_histfile_aux_num,' ')
       call add_default ('RH_COSP',         cosp_histfile_aux_num,' ')
       call add_default ('TAU_067',         cosp_histfile_aux_num,' ')
       call add_default ('EMISS_11',        cosp_histfile_aux_num,' ')
       call add_default ('MODIS_fracliq',   cosp_histfile_aux_num,' ')
       call add_default ('MODIS_asym',      cosp_histfile_aux_num,' ')
       call add_default ('MODIS_ssa',       cosp_histfile_aux_num,' ')
       call add_default ('CAL_betatot',     cosp_histfile_aux_num,' ')
       call add_default ('CAL_betatot_ice', cosp_histfile_aux_num,' ')
       call add_default ('CAL_betatot_liq', cosp_histfile_aux_num,' ')
       call add_default ('CAL_tautot',      cosp_histfile_aux_num,' ')
       call add_default ('CAL_tautot_ice',  cosp_histfile_aux_num,' ')
       call add_default ('CAL_tautot_liq',  cosp_histfile_aux_num,' ')
       call add_default ('CS_z_vol',        cosp_histfile_aux_num,' ')
       call add_default ('CS_kr_vol',       cosp_histfile_aux_num,' ')
       call add_default ('CS_g_vol',        cosp_histfile_aux_num,' ')
    end if

    ! Get pbuf indices just once at init time to avoid unnecessary character
    ! look-ups during runtime
    rei_idx        = pbuf_get_index('REI')
    rel_idx        = pbuf_get_index('REL')
    cld_idx        = pbuf_get_index('CLD')
    concld_idx     = pbuf_get_index('CONCLD')
    lsreffrain_idx = pbuf_get_index('LS_REFFRAIN')
    lsreffsnow_idx = pbuf_get_index('LS_REFFSNOW')
    cvreffliq_idx  = pbuf_get_index('CV_REFFLIQ')
    cvreffice_idx  = pbuf_get_index('CV_REFFICE')
    dpcldliq_idx   = pbuf_get_index('DP_CLDLIQ')
    dpcldice_idx   = pbuf_get_index('DP_CLDICE')
    shcldliq_idx   = pbuf_get_index('SH_CLDLIQ')
    shcldice_idx   = pbuf_get_index('SH_CLDICE')
    shcldliq1_idx  = pbuf_get_index('SH_CLDLIQ1')
    shcldice1_idx  = pbuf_get_index('SH_CLDICE1')
    dpflxprc_idx   = pbuf_get_index('DP_FLXPRC')
    dpflxsnw_idx   = pbuf_get_index('DP_FLXSNW')
    shflxprc_idx   = pbuf_get_index('SH_FLXPRC')
    shflxsnw_idx   = pbuf_get_index('SH_FLXSNW')
    lsflxprc_idx   = pbuf_get_index('LS_FLXPRC')
    lsflxsnw_idx   = pbuf_get_index('LS_FLXSNW')

    ! pbuf indices for CRM fields when using MMF
    call phys_getopts(use_SPCAM_out=use_SPCAM)
    call phys_getopts(SPCAM_microp_scheme_out=SPCAM_microp_scheme)
    if (use_SPCAM) then
       crm_qc_idx = pbuf_get_index('CRM_QC_RAD')
       crm_qi_idx = pbuf_get_index('CRM_QI_RAD')
       crm_qpl_idx = pbuf_get_index('CRM_QPL_RAD')
       crm_qpi_idx = pbuf_get_index('CRM_QPI_RAD')
       if (trim(SPCAM_microp_scheme) == 'm2005') then
          crm_qs_idx = pbuf_get_index('CRM_QS_RAD')
          crm_qr_idx = pbuf_get_index('CRM_QR_RAD')
          crm_qg_idx = pbuf_get_index('CRM_QG_RAD')
       end if
    end if

    allocate(first_run_cosp(begchunk:endchunk))
    first_run_cosp(begchunk:endchunk)=.true.
    allocate(run_cosp(1:pcols,begchunk:endchunk))
    run_cosp(1:pcols,begchunk:endchunk)=.false.

  end subroutine cospsimulator_intr_init

  ! ######################################################################################
  ! SUBROUTINE cospsimulator_intr_run
  ! ######################################################################################
  subroutine cospsimulator_intr_run(state, pbuf, cam_in, emis, coszrs, &
                                    cld_swtau, snow_tau, snow_emis,    &
                                    crm_cld_tau, crm_cld_emis          )
    use physics_types,        only: physics_state
    use physics_buffer,       only: physics_buffer_desc, pbuf_get_field
    use phys_control,         only: phys_getopts
    use camsrfexch,           only: cam_in_t
    use rad_constituents,     only: rad_cnst_get_gas
    use wv_saturation,        only: qsat_water
    use physconst,            only: pi, gravit
    use cam_history,          only: outfld,hist_fld_col_active
    use cam_history_support,  only: max_fieldname_len
    use cmparray_mod,         only: CmpDayNite, ExpDayNite
    use mod_cosp_config,      only: R_UNDEF,parasol_nrefl,Nlvgrid
    use mod_cosp,             only: cosp_simulator
    use mod_quickbeam_optics, only: size_distribution

    ! ######################################################################################
    ! Inputs
    ! ######################################################################################
    type(physics_state), intent(in),target  :: state
    type(physics_buffer_desc),      pointer :: pbuf(:)
    type(cam_in_t),      intent(in)         :: cam_in
    real(r8), intent(in) :: emis(pcols,pver)       ! cloud longwave emissivity
    real(r8), intent(in) :: coszrs(pcols)          ! cosine solar zenith angle (to tell if day or night)
    real(r8), intent(in) :: cld_swtau(pcols,pver)  ! RRTM cld_swtau_in, read in using this variable
    real(r8), intent(inout) :: snow_tau(pcols,pver)   ! RRTM grid-box mean SW snow optical depth, used for CAM5 simulations
    real(r8), intent(inout) :: snow_emis(pcols,pver)  ! RRTM grid-box mean LW snow optical depth, used for CAM5 simulations
    real(r8), intent(in), optional, dimension(:,:,:,:) :: crm_cld_tau, crm_cld_emis

    ! ######################################################################################
    ! Local variables
    ! ######################################################################################
    integer :: lchnk                             ! chunk identifier
    integer :: ncol                              ! number of active atmospheric columns
    integer :: i,k

    ! COSP-related local vars
    type(cosp_outputs)        :: cospOUT                  ! COSP simulator outputs
    type(cosp_optical_inputs) :: cospIN                   ! COSP optical (or derived?) fields needed by simulators
    type(cosp_column_inputs)  :: cospstateIN              ! COSP model fields needed by simulators

    logical :: use_precipitation_fluxes                   ! True if precipitation fluxes are input to the algorithm
    real(r8), parameter :: emsfc_lw = 0.99_r8             ! longwave emissivity of surface at 10.5 microns

    real(r8) :: cam_reff(pcols,pver,nhydro)              ! effective radius for cosp input
    real(r8) :: cam_hydro(pcols,pver,nhydro)             ! hydrometeor mixing ratios and precip fluxes (kg/kg)

    ! ######################################################################################
    ! Simulator output info
    ! ######################################################################################
    integer, parameter :: nf_radar=6                     ! number of radar outputs
    integer, parameter :: nf_lidar=28                    ! number of lidar outputs !!+cosp1.4
    integer, parameter :: nf_isccp=9                     ! number of isccp outputs
    integer, parameter :: nf_misr=1                      ! number of misr outputs
    integer, parameter :: nf_modis=20                    ! number of modis outputs

    ! Cloudsat outputs
    character(len=max_fieldname_len),dimension(nf_radar),parameter ::        &
         fname_radar = (/'CFAD_DBZE94_CS','CLD_CAL_NOTCS ','DBZE_CS       ', &
                         'CLDTOT_CALCS  ','CLDTOT_CS     ','CLDTOT_CS2    '/)
    ! CALIPSO outputs
    character(len=max_fieldname_len),dimension(nf_lidar),parameter :: &
         fname_lidar=(/'CLDLOW_CAL     ','CLDMED_CAL     ','CLDHGH_CAL     ','CLDTOT_CAL     ','CLD_CAL        ',&
                       'RFL_PARASOL    ','CFAD_SR532_CAL ','ATB532_CAL     ','MOL532_CAL     ','CLD_CAL_LIQ    ',&
                       'CLD_CAL_ICE    ','CLD_CAL_UN     ','CLD_CAL_TMP    ','CLD_CAL_TMPLIQ ','CLD_CAL_TMPICE ',&
                       'CLD_CAL_TMPUN  ','CLDTOT_CAL_ICE ','CLDTOT_CAL_LIQ ','CLDTOT_CAL_UN  ','CLDHGH_CAL_ICE ',&
                       'CLDHGH_CAL_LIQ ','CLDHGH_CAL_UN  ','CLDMED_CAL_ICE ','CLDMED_CAL_LIQ ','CLDMED_CAL_UN  ',&
                       'CLDLOW_CAL_ICE ','CLDLOW_CAL_LIQ ','CLDLOW_CAL_UN  '/)  !+COSP1.4 (all same # characters!)
    ! ISCCP outputs
    character(len=max_fieldname_len),dimension(nf_isccp),parameter :: &
         fname_isccp=(/'FISCCP1_COSP    ','CLDTOT_ISCCP    ','MEANCLDALB_ISCCP',&
                       'MEANPTOP_ISCCP  ','TAU_ISCCP       ','CLDPTOP_ISCCP   ','MEANTAU_ISCCP   ',&
                       'MEANTB_ISCCP    ','MEANTBCLR_ISCCP '/)
    ! MISR outputs
    character(len=max_fieldname_len),dimension(nf_misr),parameter :: &
         fname_misr=(/'CLD_MISR '/)
    ! MODIS outputs
    character(len=max_fieldname_len),dimension(nf_modis) :: &
         fname_modis=(/'CLTMODIS    ','CLWMODIS    ','CLIMODIS    ','CLHMODIS    ','CLMMODIS    ',&
                       'CLLMODIS    ','TAUTMODIS   ','TAUWMODIS   ','TAUIMODIS   ','TAUTLOGMODIS',&
                       'TAUWLOGMODIS','TAUILOGMODIS','REFFCLWMODIS','REFFCLIMODIS',&
                       'PCTMODIS    ','LWPMODIS    ','IWPMODIS    ','CLMODIS     ','CLRIMODIS   ',&
                       'CLRLMODIS   '/)

    logical :: run_radar(nf_radar,pcols)                 ! logical telling you if you should run radar simulator
    logical :: run_lidar(nf_lidar,pcols)                 ! logical telling you if you should run lidar simulator
    logical :: run_isccp(nf_isccp,pcols)                 ! logical telling you if you should run isccp simulator
    logical :: run_misr(nf_misr,pcols)                   ! logical telling you if you should run misr simulator
    logical :: run_modis(nf_modis,pcols)                 ! logical telling you if you should run modis simulator

    ! CAM pointers to get variables from radiation interface (get from rad_cnst_get_gas)
    real(r8), pointer, dimension(:,:) :: q               ! specific humidity (kg/kg)
    real(r8), pointer, dimension(:,:) :: o3              ! Mass mixing ratio 03

    ! CAM pointers to get variables from the physics buffer
    real(r8), pointer, dimension(:,:) :: cld             ! cloud fraction, tca - total_cloud_amount (0-1)
    real(r8), pointer, dimension(:,:) :: concld          ! concld fraction, cca - convective_cloud_amount (0-1)

    ! COSPv2 stuff
    character(len=256),dimension(100) :: cosp_status
    integer :: nerror

    ! subcolumn and optics stuff
    real(wp),dimension(:,:,:),  allocatable :: frac_prec
    real(wp),dimension(:,:,:,:),allocatable :: mr_hydro, Reff, Np

    ! For MMF/SPCAM
    logical :: use_SPCAM
    character(len=16) :: SPCAM_microp_scheme

    call t_startf("cosp_init")
    ! ######################################################################################
    ! Initialization
    ! ######################################################################################
    ! Find the chunk and ncol from the state vector
    lchnk = state%lchnk    ! state variable contains a number of columns, one chunk
    ncol  = state%ncol     ! number of columns in the chunk

    ! ######################################################################################
    ! DECIDE WHICH COLUMNS YOU ARE GOING TO RUN COSP ON....
    ! ######################################################################################

    ! run_cosp is set for each column in each chunk in the first timestep of the run
    ! hist_fld_col_active in cam_history.F90 is used to decide if you need to run cosp.
    if (first_run_cosp(lchnk)) then
       run_cosp(1:ncol,lchnk)=.false.
       run_radar(1:nf_radar,1:ncol)=.false.
       run_lidar(1:nf_lidar,1:ncol)=.false.
       run_isccp(1:nf_isccp,1:ncol)=.false.
       run_misr(1:nf_misr,1:ncol)=.false.
       run_modis(1:nf_modis,1:ncol)=.false.

       if (cosp_lradar_sim) then
          do i=1,nf_radar
             run_radar(i,1:pcols)=hist_fld_col_active(fname_radar(i),lchnk,pcols)
          end do
       end if
       if (cosp_llidar_sim) then
          do i=1,nf_lidar
             run_lidar(i,1:pcols)=hist_fld_col_active(fname_lidar(i),lchnk,pcols)
          end do
       end if
       if (cosp_lisccp_sim) then
          do i=1,nf_isccp
             run_isccp(i,1:pcols)=hist_fld_col_active(fname_isccp(i),lchnk,pcols)
          end do
       end if
       if (cosp_lmisr_sim) then
          do i=1,nf_misr
             run_misr(i,1:pcols)=hist_fld_col_active(fname_misr(i),lchnk,pcols)
          end do
       end if
       if (cosp_lmodis_sim) then
          do i=1,nf_modis
             run_modis(i,1:pcols)=hist_fld_col_active(fname_modis(i),lchnk,pcols)
          end do
       end if

       do i=1,ncol
          if ((any(run_radar(:,i))) .or. (any(run_lidar(:,i))) .or. (any(run_isccp(:,i))) &
               .or. (any(run_misr(:,i))) .or. (any(run_modis(:,i)))) then
             run_cosp(i,lchnk)=.true.
          end if
       end do

       first_run_cosp(lchnk)=.false.
    endif

    ! ######################################################################################
    ! GET CAM GEOPHYSICAL VARIABLES NEEDED FOR COSP INPUT
    ! ######################################################################################

    call phys_getopts(use_SPCAM_out=use_SPCAM)
    call phys_getopts(SPCAM_microp_scheme_out=SPCAM_microp_scheme)

    call t_stopf('cosp_init')

    ! initialize local variables
    cam_reff(:,:,:) = 0
    cam_hydro(:,:,:) = 0

    ! ######################################################################################
    ! Construct COSP output derived type.
    ! ######################################################################################
    call t_startf("cosp_construct_cospOUT")
    call construct_cosp_outputs(ncol,nSubcol,pver,Nlvgrid,0,cospOUT)
    call t_stopf("cosp_construct_cospOUT")

    ! ######################################################################################
    ! Construct and populate COSP input types
    ! ######################################################################################
    ! Model state
    call t_startf("cosp_construct_cospstateIN")
    call construct_cospstateIN(ncol,pver,0,cospstateIN)

    ! Convert from coordinate variables radians to degrees_north and degrees_east
    ! Note that -90 < lat < 90; 0 < lon < 360
    cospstateIN%lat = state%lat(1:ncol) * 180._r8 / pi
    cospstateIN%lon = state%lon(1:ncol) * 180._r8 / pi
    cospstateIN%at  = state%t(1:ncol,1:pver)

    ! Radiative constituent interface variables:
    ! specific humidity (q), 03 mass mixing ratio
    call rad_cnst_get_gas(0,'H2O', state, pbuf,  q)
    call rad_cnst_get_gas(0,'O3',  state, pbuf,  o3)
    cospstateIN%qv  = q(1:ncol,1:pver)
    cospstateIN%o3  = o3(1:ncol,1:pver)

    ! Skin temperature comes from cam_in (from land model?)
    cospstateIN%skt = cam_in%ts(1:ncol)

    ! Set sunlit flag
    do i=1,ncol
       if ((coszrs(i) > 0.0_r8)) then
          cospstateIN%sunlit(i) = 1
       else
          cospstateIN%sunlit(i) = 0
       endif
    enddo

    ! Set land mask
    do i = 1,ncol
       if (cam_in%landfrac(i) > 0.01_r8) then
          cospstateIN%land(i) = 1
       else
          cospstateIN%land(i) = 0
       end if
    end do
    cospstateIN%pfull                  = state%pmid(1:ncol,1:pver)
    cospstateIN%phalf(1:ncol,1)        = 0._r8
    cospstateIN%phalf(1:ncol,2:pver+1) = state%pint(1:ncol,2:pver+1)
    ! Add surface height (surface geopotential/gravity) to convert CAM heights based on
    ! geopotential above surface into height above sea level
    ! TODO: This looks like a bug to me; hgt_matrix_half USED to be size nlev,
    ! but now is nlev+1, so should be set appropropriately (here index k is
    ! given index k+1 from state%zi, which is incorrect).
    do k = 1,pver
       cospstateIN%hgt_matrix(1:ncol,k)      = state%zm(1:ncol,k  ) + state%phis(1:ncol) / gravit
       cospstateIN%hgt_matrix_half(1:ncol,k) = state%zi(1:ncol,k+1) + state%phis(1:ncol) / gravit
    end do
    cospstateIN%hgt_matrix_half(1:ncol,pver+1) = 0._r8
    cospstateIN%surfelev(1:ncol) = state%zi(1:ncol,pver+1)
    call t_stopf("cosp_construct_cospstateIN")

    ! Optical inputs
    call t_startf("cosp_construct_cospIN")
    call construct_cospIN(ncol,nSubcol,pver,cospIN)
    cospIN%emsfc_lw = emsfc_lw
    if (cosp_lradar_sim) cospIN%rcfg_cloudsat = rcfg_cs(lchnk)
    call t_stopf("cosp_construct_cospIN")

    ! *NOTE* Fields passed into subsample_and_optics are ordered from TOA-2-SFC.
    call t_startf("cosp_subsample_and_optics")

    ! mixing ratios, effective radii and precipitation fluxes for cloud and precipitation
    allocate(mr_hydro(ncol, nSubcol, pver, nhydro), &
             reff(ncol, nSubcol, pver, nhydro), &
             np(ncol, nSubcol, pver, nhydro), &
             frac_prec(ncol, nSubcol, pver))

    ! If using SPCAM, then do not use large-scale precipitation fluxes (we will
    ! populate subgrid fields directly with precipitation mixing ratios)
    call phys_getopts(use_SPCAM_out=use_SPCAM)
    if (.not. use_SPCAM) then

       ! Use precipitation fluxes instead of mixing ratios (and convert to mixing
       ! ratios later)
       use_precipitation_fluxes = .true.

       ! Get gridbox mean hydrometeor information from GCM fields
       ! Get subcolumn hydrometeor information by subsampling, and calculate
       ! subcolumn optical properties
       call pbuf_get_field(pbuf, cld_idx,    cld   )
       call pbuf_get_field(pbuf, concld_idx, concld)
       call get_cam_hydros(state, pbuf, cld, cam_hydro, cam_reff)
       call cosp_subsample(ncol, pver, nSubcol, nhydro, overlap, &
            use_precipitation_fluxes, cld(1:ncol, 1:pver), concld(1:ncol, 1:pver), &
            cam_hydro(1:ncol,1:pver,1:N_HYDRO), cam_reff(1:ncol,1:pver,1:N_HYDRO), &
            cld_swtau(1:ncol,1:pver), cld_swtau(1:ncol,1:pver), &
            emis(1:ncol,1:pver), emis(1:ncol,1:pver), &
            snow_tau(1:ncol,1:pver), snow_emis(1:ncol,1:pver), &
            state%ps(1:ncol), cospstateIN, cospIN, &
            mr_hydro, Reff, Np, frac_prec)  ! new outputs, need to be added ...
       call calc_cosp_optics( &
          ncol, pver, nSubcol, nhydro, &
          lidar_ice_type, sd_cs(lchnk), &
          mr_hydro, Reff, Np, frac_prec, &
          cospstateIN, cospIN, &
          snow_tau(1:ncol,1:pver), snow_emis(1:ncol,1:pver) &
       )

    else

       ! Do not use precip fluxes for CRM because we have and will use the
       ! precipitation mixing ratios
       use_precipitation_fluxes = .false.

       ! Get subcolumn hydrometeor information from CRM fields, and calculate
       ! subcolumn optical properties
       call t_startf('cosp_get_crm_hydros')
       if (present(crm_cld_tau).and.present(crm_cld_emis)) then
          call get_crm_hydros( &
             state, pbuf, crm_cld_tau, crm_cld_emis, &
             mr_hydro, reff, np, frac_prec, cospIN &
          )
       else
          call endrun('cospsimulator_intr_run: crm_cld_tau or crm_cld_emis not present')
       end if
       call t_stopf('cosp_get_crm_hydros')
       call t_startf('cosp_calc_cosp_optics')
       call calc_cosp_optics( &
          ncol, pver, nSubcol, nhydro, &
          lidar_ice_type, sd_cs(lchnk), &
          mr_hydro, Reff, Np, frac_prec, &
          cospstateIN, cospIN &
       )
       call t_stopf('cosp_calc_cosp_optics')

    end if
    call t_stopf("cosp_subsample_and_optics")

    ! ######################################################################################
    ! Write COSP inputs to output file for offline use.
    ! ######################################################################################
    call t_startf("cosp_histfile_aux")
    if (cosp_histfile_aux) then
       call outfld('SCOPS_OUT',      cospIN%frac_out(1:ncol,:,:),  ncol,lchnk)
       call outfld('COSP_PREC_FRAC', frac_prec(1:ncol,:,:),        ncol,lchnk)
       call outfld('PS_COSP',        state%ps(1:ncol),             ncol,lchnk)
       call outfld('TS_COSP',        cospstateIN%skt,              ncol,lchnk)
       call outfld('P_COSP',         cospstateIN%pfull,            ncol,lchnk)
       call outfld('PH_COSP',        cospstateIN%phalf,            ncol,lchnk)
       call outfld('ZLEV_COSP',      cospstateIN%hgt_matrix,       ncol,lchnk)
       call outfld('ZLEV_HALF_COSP', cospstateIN%hgt_matrix_half,  ncol,lchnk)
       call outfld('T_COSP',         cospstateIN%at,               ncol,lchnk)
       call outfld('RH_COSP',        cospstateIN%qv,               ncol,lchnk)
       call outfld('TAU_067',        cospIN%tau_067,               ncol,lchnk)
       call outfld('EMISS_11',       cospIN%emiss_11,              ncol,lchnk)
       call outfld('MODIS_asym',     cospIN%asym,                  ncol,lchnk)
       call outfld('MODIS_ssa',      cospIN%ss_alb,                ncol,lchnk)
       call outfld('MODIS_fracliq',  cospIN%fracLiq,               ncol,lchnk)
    end if
    call t_stopf("cosp_histfile_aux")

    ! done with these now ...
    deallocate(mr_hydro)
    deallocate(Reff)
    deallocate(Np)
    deallocate(frac_prec)

    ! ######################################################################################
    ! Call COSP
    ! ######################################################################################
    call t_startf("cosp_simulator")
    cosp_status = COSP_SIMULATOR(cospIN, cospstateIN, cospOUT, start_idx=1, stop_idx=ncol, debug=.false.)

    ! Check status flags
    nerror = 0
    do i = 1, ubound(cosp_status, 1)
       if (len_trim(cosp_status(i)) > 0) then
          write(iulog,*) "cosp_simulator: ERROR: "//trim(cosp_status(i))
          nerror = nerror + 1
       end if
    end do
    if (nerror > 0) then
       call endrun('cospsimulator_intr_run: error return from cosp_simulator')
    end if
    call t_stopf("cosp_simulator")

    ! ######################################################################################
    ! Set dark-scenes to fill value. Only done for passive simulators
    ! TODO: we should NOT have to do this, this is supposed to be done
    ! internally in COSP based on the sunlit flag that is passed, but we need to
    ! handle some cases explicitly here because of existing bugs in COSP.
    ! ######################################################################################
    call t_startf("cosp_mask_passive")
    ! ISCCP simulator
    ! TODO: these should not be set undefined for night columns (brightness
    ! temperature is valid at night), but we set these to fillvalues for now
    ! to keep results BFB
    if (cosp_lisccp_sim) then
       where(cospstateIN%sunlit(1:ncol) .eq. 0)
          cospOUT%isccp_meantb(1:ncol)       = R_UNDEF
          cospOUT%isccp_meantbclr(1:ncol)    = R_UNDEF
       end where
    endif

    ! MISR simulator
    ! TODO: the MISR simulator is supposed to handle this for us, but there
    ! is a bug in COSP that leaves these set for night columns, so we need to
    ! explicitly handle sunlit vs not sunlit here for now
    if (cosp_lmisr_sim) then
       do i = 1,ncol
         if (cospstateIN%sunlit(i) == 0) then
            cospOUT%misr_fq(i,:,:) = R_UNDEF
         end if
       end do
    end if

    ! MODIS simulator
    ! TODO: again, this is to workaround a bug in COSP, where these two
    ! variables are not masked properly for night columns
    if (cosp_lmodis_sim) then
       do i = 1,ncol
          if (cospstateIN%sunlit(i) == 0) then
             cospOUT%modis_Optical_Thickness_vs_ReffLiq(i,:,:) = R_UNDEF
             cospOUT%modis_Optical_Thickness_vs_ReffIce(i,:,:) = R_UNDEF
          end if
       end do
    end if
    call t_stopf("cosp_mask_passive")

    call t_startf("cosp_history_output")
    call cosp_history_output(ncol, lchnk, cospIN, cospOUT)
    call t_stopf("cosp_history_output")

    ! ######################################################################################
    ! Clean up
    ! ######################################################################################
    call t_startf("destroy_cospIN")
    call destroy_cospIN(cospIN)
    call t_stopf("destroy_cospIN")
    call t_startf("destroy_cospstateIN")
    call destroy_cospstateIN(cospstateIN)
    call t_stopf("destroy_cospstateIN")
    call t_startf("destroy_cospOUT")
    call destroy_cosp_outputs(cospOUT)
    call t_stopf("destroy_cospOUT")

  end subroutine cospsimulator_intr_run


  subroutine get_cam_hydros(state, pbuf, cld, cam_hydro, cam_reff)

    use physics_types,        only: physics_state
    use physics_buffer,       only: physics_buffer_desc, pbuf_get_field
    use interpolate_data,     only: lininterp_init,lininterp,lininterp_finish,interp_type
    use constituents,         only: cnst_get_ind

    type(physics_state), intent(in) :: state
    type(physics_buffer_desc), pointer :: pbuf(:)
    real(r8), intent(in) :: cld(:,:)
    real(r8), intent(out) :: cam_hydro(:,:,:)
    real(r8), intent(out) :: cam_reff(:,:,:)

    real(r8), pointer, dimension(:,:) :: rel             ! liquid effective drop radius (microns)
    real(r8), pointer, dimension(:,:) :: rei             ! ice effective drop size (microns)
    real(r8), pointer, dimension(:,:) :: ls_reffrain     ! rain effective drop radius (microns)
    real(r8), pointer, dimension(:,:) :: ls_reffsnow     ! snow effective drop size (microns)
    real(r8), pointer, dimension(:,:) :: cv_reffliq      ! convective cld liq effective drop radius (microns)
    real(r8), pointer, dimension(:,:) :: cv_reffice      ! convective cld ice effective drop size (microns)

    ! precip flux pointers (use for cam4 or cam5)
    ! Added pointers;  pbuf in zm_conv_intr.F90, calc in zm_conv.F90
    real(r8), pointer, dimension(:,:) :: dp_flxprc       ! deep interface gbm flux_convective_cloud_rain+snow (kg m^-2 s^-1)
    real(r8), pointer, dimension(:,:) :: dp_flxsnw       ! deep interface gbm flux_convective_cloud_snow (kg m^-2 s^-1)
    ! More pointers;  pbuf in convect_shallow.F90, calc in hk_conv.F90/convect_shallow.F90 (CAM4), uwshcu.F90 (CAM5)
    real(r8), pointer, dimension(:,:) :: sh_flxprc       ! shallow interface gbm flux_convective_cloud_rain+snow (kg m^-2 s^-1)
    real(r8), pointer, dimension(:,:) :: sh_flxsnw       ! shallow interface gbm flux_convective_cloud_snow (kg m^-2 s^-1)
    ! More pointers;  pbuf in stratiform.F90, getting from pbuf here
    ! a) added as output to pcond subroutine in cldwat.F90 and to nmicro_pcond subroutine in cldwat2m_micro.F90
    real(r8), pointer, dimension(:,:) :: ls_flxprc       ! stratiform interface gbm flux_cloud_rain+snow (kg m^-2 s^-1)
    real(r8), pointer, dimension(:,:) :: ls_flxsnw       ! stratiform interface gbm flux_cloud_snow (kg m^-2 s^-1)

    !! cloud mixing ratio pointers (note: large-scale in state)
    ! More pointers;  pbuf in convect_shallow.F90 (cam4) or stratiform.F90 (cam5)
    ! calc in hk_conv.F90 (CAM4 should be 0!), uwshcu.F90 but then affected by micro so values from stratiform.F90 (CAM5)
    real(r8), pointer, dimension(:,:) :: sh_cldliq       ! shallow gbm cloud liquid water (kg/kg)
    real(r8), pointer, dimension(:,:) :: sh_cldice       ! shallow gbm cloud ice water (kg/kg)
    ! More pointers;  pbuf in zm_conv_intr.F90, calc in zm_conv.F90, 0 for CAM4 and CAM5 (same convection scheme)
    real(r8), pointer, dimension(:,:) :: dp_cldliq       ! deep gbm cloud liquid water (kg/kg)
    real(r8), pointer, dimension(:,:) :: dp_cldice       ! deep gmb cloud ice water (kg/kg)

    integer, parameter :: extrap_method = 1  ! sets extrapolation method to boundary value (1)

    ! Local vars related to calculations to go from CAM input to COSP input
    ! cosp convective value includes both deep and shallow convection
    real(r8) :: rain_cv(pcols,pverp)                     ! interface flux_convective_cloud_rain (kg m^-2 s^-1)
    real(r8) :: snow_cv(pcols,pverp)                     ! interface flux_convective_cloud_snow (kg m^-2 s^-1)

    ! Microphysics variables
    integer :: ixcldliq                                   ! cloud liquid amount index for state%q
    integer :: ixcldice                                   ! cloud ice amount index

    type(interp_type)  :: interp_wgts
    integer :: ncol, k, i


    ncol = state%ncol

    ! Convective cloud mixing ratios (use for cam4 and cam5)
    call pbuf_get_field(pbuf, dpcldliq_idx, dp_cldliq)
    call pbuf_get_field(pbuf, dpcldice_idx, dp_cldice)

    ! Added to pbuf in stratiform.F90
    call pbuf_get_field(pbuf, shcldliq1_idx, sh_cldliq)
    call pbuf_get_field(pbuf, shcldice1_idx, sh_cldice)

    ! Precipitation fluxes
    call pbuf_get_field(pbuf, dpflxprc_idx, dp_flxprc)
    call pbuf_get_field(pbuf, dpflxsnw_idx, dp_flxsnw)
    call pbuf_get_field(pbuf, shflxprc_idx, sh_flxprc)
    call pbuf_get_field(pbuf, shflxsnw_idx, sh_flxsnw)
    call pbuf_get_field(pbuf, lsflxprc_idx, ls_flxprc)
    call pbuf_get_field(pbuf, lsflxsnw_idx, ls_flxsnw)

    ! Effective radii
    call pbuf_get_field(pbuf, rel_idx, rel)
    call pbuf_get_field(pbuf, rei_idx, rei)
    call pbuf_get_field(pbuf, lsreffrain_idx, ls_reffrain)
    call pbuf_get_field(pbuf, lsreffsnow_idx, ls_reffsnow)
    call pbuf_get_field(pbuf, cvreffliq_idx,  cv_reffliq )
    call pbuf_get_field(pbuf, cvreffice_idx,  cv_reffice )

    ! add together deep and shallow convection precipitation fluxes, recall *_flxprc variables are rain+snow
    rain_cv(1:ncol,1:pverp) = (sh_flxprc(1:ncol,1:pverp) - sh_flxsnw(1:ncol,1:pverp)) &
                            + (dp_flxprc(1:ncol,1:pverp) - dp_flxsnw(1:ncol,1:pverp))
    snow_cv(1:ncol,1:pverp) =  sh_flxsnw(1:ncol,1:pverp) + dp_flxsnw(1:ncol,1:pverp)

    ! interpolate interface precip fluxes to mid points
    do i=1,ncol
       ! find weights (pressure weighting?)
       call lininterp_init(state%zi(i,1:pverp),pverp,state%zm(i,1:pver),pver,extrap_method,interp_wgts)
       ! interpolate  lininterp1d(arrin, nin, arrout, nout, interp_wgts)
       ! note: lininterp is an interface, contains lininterp1d -- code figures out to use lininterp1d.
       call lininterp(  rain_cv(i,1:pverp), pverp, cam_hydro(i,1:pver,I_CVRAIN), pver, interp_wgts)
       call lininterp(  snow_cv(i,1:pverp), pverp, cam_hydro(i,1:pver,I_CVSNOW), pver, interp_wgts)
       call lininterp(ls_flxprc(i,1:pverp), pverp, cam_hydro(i,1:pver,I_LSRAIN), pver, interp_wgts)
       call lininterp(ls_flxsnw(i,1:pverp), pverp, cam_hydro(i,1:pver,I_LSSNOW), pver, interp_wgts)
       call lininterp_finish(interp_wgts)
       ! ls_flxprc is for rain+snow, find rain_ls_interp by subtracting off snow_ls_interp
       cam_hydro(i,1:pver,I_LSRAIN) = cam_hydro(i,1:pver,I_LSRAIN) - cam_hydro(i,1:pver,I_LSSNOW)
    end do

    ! CAM5 cloud mixing ratio calculations
    ! Note: Although CAM5 has non-zero convective cloud mixing ratios that affect the model state,
    ! Convective cloud water is NOT part of radiation calculations.
    call cnst_get_ind('CLDLIQ',ixcldliq)
    call cnst_get_ind('CLDICE',ixcldice)
    do k=1,pver
       do i=1,ncol
          if (cld(i,k) .gt. 0._r8) then
             ! note: convective mixing ratio is the sum of shallow and deep convective clouds in CAM5
             cam_hydro(i,k,I_CVCLIQ) = sh_cldliq(i,k) + dp_cldliq(i,k)
             cam_hydro(i,k,I_CVCICE) = sh_cldice(i,k) + dp_cldice(i,k)
             cam_hydro(i,k,I_LSCLIQ) = state%q(i,k,ixcldliq)  ! only includes stratiform (kg/kg)
             cam_hydro(i,k,I_LSCICE) = state%q(i,k,ixcldice)  ! only includes stratiform (kg/kg)
          else
             cam_hydro(i,k,I_CVCLIQ) = 0._r8
             cam_hydro(i,k,I_CVCICE) = 0._r8
             cam_hydro(i,k,I_LSCLIQ) = 0._r8
             cam_hydro(i,k,I_LSCICE) = 0._r8
          end if
       end do
    end do

    ! More sizes added to physics buffer in stratiform.F90
    cam_reff(1:ncol,1:pver,I_LSCLIQ) = rel(1:ncol,1:pver)*1.e-6_r8          ! same as effc and effliq in stratiform.F90
    cam_reff(1:ncol,1:pver,I_LSCICE) = rei(1:ncol,1:pver)*1.e-6_r8          ! same as effi and effice in stratiform.F90
    cam_reff(1:ncol,1:pver,I_LSRAIN) = ls_reffrain(1:ncol,1:pver)*1.e-6_r8  ! calculated in cldwat2m_micro.F90, passed to stratiform.F90
    cam_reff(1:ncol,1:pver,I_LSSNOW) = ls_reffsnow(1:ncol,1:pver)*1.e-6_r8  ! calculated in cldwat2m_micro.F90, passed to stratiform.F90
    cam_reff(1:ncol,1:pver,I_CVCLIQ) = cv_reffliq(1:ncol,1:pver)*1.e-6_r8   ! calculated in stratiform.F90, not actually used in radiation
    cam_reff(1:ncol,1:pver,I_CVCICE) = cv_reffice(1:ncol,1:pver)*1.e-6_r8   ! calculated in stratiform.F90, not actually used in radiation
    cam_reff(1:ncol,1:pver,I_CVRAIN) = ls_reffrain(1:ncol,1:pver)*1.e-6_r8  ! same as stratiform per Andrew
    cam_reff(1:ncol,1:pver,I_CVSNOW) = ls_reffsnow(1:ncol,1:pver)*1.e-6_r8  ! same as stratiform per Andrew
    cam_reff(1:ncol,1:pver,I_LSGRPL) = 0._r8                                ! using radar default reff

  end subroutine get_cam_hydros


  subroutine get_crm_hydros(state, pbuf, dtau, dems, mr_hydro, reff, np, frac_prec, cospIN)
    use physics_types,  only: physics_state
    use physics_buffer, only: physics_buffer_desc, pbuf_get_field
    use crmdims,        only: crm_nx_rad, crm_ny_rad, crm_nz
    type(physics_state), intent(in) :: state
    type(physics_buffer_desc), pointer :: pbuf(:)
    real(r8), intent(in ), dimension(:,:,:,:) :: dtau, dems
    real(r8), intent(inout), dimension(:,:,:,:) :: mr_hydro, reff, np
    real(r8), intent(inout), dimension(:,:,:) :: frac_prec
    type(cosp_optical_inputs), intent(inout) :: cospIN

    real(r8), pointer, dimension(:,:) :: rel, rei
    real(r8), pointer, dimension(:,:,:,:) :: crm_qc, crm_qi, crm_qr, crm_qs

    integer :: ncol, icol, isubcol, ilev, ix, iy, iz

    ! Initialize
    mr_hydro = 0
    reff = 0
    np = 0
    frac_prec = 0
    cospIN%tau_067 = 0
    cospIN%emiss_11 = 0

    ! Get fields from pbuf; note this assumes all precipitating ice is snow.
    ! While the sam1mom physics internally makes some assumptions about graupel
    ! vs snow, we have not yet separated that out as a diagnostic output.
    call pbuf_get_field(pbuf, rel_idx, rel)
    call pbuf_get_field(pbuf, rei_idx, rei)
    call pbuf_get_field(pbuf, crm_qc_idx, crm_qc)
    call pbuf_get_field(pbuf, crm_qi_idx, crm_qi)
    call pbuf_get_field(pbuf, crm_qpl_idx, crm_qr)
    call pbuf_get_field(pbuf, crm_qpi_idx, crm_qs)

    ! Set values column by column
    ncol = state%ncol
    do icol = 1,ncol
       do iz = 1,crm_nz
          do iy = 1,crm_ny_rad
             do ix = 1,crm_nx_rad

                ! Find indices
                ilev = pver - iz + 1
                isubcol = (iy - 1) * crm_nx_rad + ix

                ! Set hydrometeor mixing ratios
                mr_hydro(icol,isubcol,ilev,I_LSCLIQ) = crm_qc(icol,ix,iy,iz)
                mr_hydro(icol,isubcol,ilev,I_LSCICE) = crm_qi(icol,ix,iy,iz)
                mr_hydro(icol,isubcol,ilev,I_LSRAIN) = crm_qr(icol,ix,iy,iz)
                mr_hydro(icol,isubcol,ilev,I_LSSNOW) = crm_qs(icol,ix,iy,iz)

                ! Set effective radii for each hydrometeor type
                ! TODO: for now, this is using the effective radii from the
                ! GCM scheme, which I believe is constant effective radii,
                ! different values over land vs ocean. This should do something
                ! more intelligent, but this is also what radiation uses when
                ! using sam1mom microphysics.
                ! TODO: confirm that when reff is zero, defaults are used
                reff(icol,isubcol,ilev,I_LSCLIQ) = rel(icol,ilev) * 1.e-6_r8
                reff(icol,isubcol,ilev,I_LSCICE) = rei(icol,ilev) * 1.e-6_r8
                reff(icol,isubcol,ilev,I_LSRAIN) = 0
                reff(icol,isubcol,ilev,I_LSSNOW) = 0
                reff(icol,isubcol,ilev,I_LSGRPL) = 0

                ! Set precipitation fraction
                if (crm_qr(icol,ix,iy,iz) + crm_qs(icol,ix,iy,iz) > 0._r8) then
                   frac_prec(icol,isubcol,ilev) = 1
                else
                   frac_prec(icol,isubcol,ilev) = 0
                end if

                ! Set optical properties
                cospIN%tau_067(icol,isubcol,ilev) = dtau(icol,ix,iy,iz)
                cospIN%emiss_11(icol,isubcol,ilev) = dems(icol,ix,iy,iz)

             end do
          end do
       end do
    end do

  end subroutine get_crm_hydros


  subroutine cosp_history_output(ncol, lchnk, cospIN, cospOUT)
    use cam_history, only: outfld
    integer, intent(in) :: ncol, lchnk
    type(cosp_optical_inputs), intent(in) :: cospIN
    type(cosp_outputs), intent(in) :: cospOUT

    ! Local variables
    real(r8) :: cfad_dbze94(pcols,ndbze_cosp,Nlvgrid)
    real(r8) :: cfad_lidarsr532(pcols,nsr_cosp,Nlvgrid)
    real(r8) :: meancldalb_isccp(pcols)
    real(r8) :: meanptop_isccp(pcols)
    real(r8) :: cld_cal(pcols,Nlvgrid)
    real(r8) :: cld_cal_liq(pcols,Nlvgrid)
    real(r8) :: cld_cal_ice(pcols,Nlvgrid)
    real(r8) :: cld_cal_un(pcols,Nlvgrid)
    real(r8) :: cld_cal_tmp(pcols,Nlvgrid)
    real(r8) :: cld_cal_tmpliq(pcols,Nlvgrid)
    real(r8) :: cld_cal_tmpice(pcols,Nlvgrid)
    real(r8) :: cld_cal_tmpun(pcols,Nlvgrid)
    real(r8) :: meantau_isccp(pcols)
    real(r8) :: meantb_isccp(pcols)
    real(r8) :: meantbclr_isccp(pcols)
    real(r8) :: cldtot_cs(pcols)
    real(r8) :: cldtot_cs2(pcols)
    real(r8) :: refl_parasol(pcols,nsza_cosp)
    real(r8) :: tautmodis(pcols)
    real(r8) :: tauwmodis(pcols)
    real(r8) :: tauimodis(pcols)
    real(r8) :: tautlogmodis(pcols)
    real(r8) :: tauwlogmodis(pcols)
    real(r8) :: tauilogmodis(pcols)
    real(r8) :: reffclwmodis(pcols)
    real(r8) :: reffclimodis(pcols)
    real(r8) :: pctmodis(pcols)
    real(r8) :: lwpmodis(pcols)
    real(r8) :: iwpmodis(pcols)

    ! ISCCP OUTPUTS
    if (cosp_lisccp_sim) then
       call outfld('FISCCP1_COSP',    cospOUT%isccp_fq(1:ncol,:,:)      , ncol, lchnk)
       call outfld('CLDTOT_ISCCP',    cospOUT%isccp_totalcldarea(1:ncol), ncol, lchnk)
       call outfld('MEANTB_ISCCP',    cospOUT%isccp_meantb(1:ncol),       ncol, lchnk)
       call outfld('MEANTBCLR_ISCCP', cospOUT%isccp_meantbclr(1:ncol),    ncol, lchnk)
       ! The following are weighted by total cloud fraction
       ! TODO: this is probably a bug. meancldalb, meanptop, and meantau should
       ! be set to missing values when totalcldarea = 0. This is done within the
       ! ISCCP simulator, but by not checking against meancldalb, meanptop, and
       ! meantau missing values here, we are mistakenly resetting these missing
       ! values to zero by multiplying by isccp_totalcldarea.
       where (cospOUT%isccp_totalcldarea(1:ncol) .eq. R_UNDEF)
          meancldalb_isccp(:ncol) = R_UNDEF
          meanptop_isccp(:ncol)   = R_UNDEF
          meantau_isccp(:ncol)    = R_UNDEF
       elsewhere
          meancldalb_isccp(:ncol) = cospOUT%isccp_meanalbedocld(:ncol)*cospOUT%isccp_totalcldarea(1:ncol)
          meanptop_isccp(:ncol)   = cospOUT%isccp_meanptop(:ncol)*cospOUT%isccp_totalcldarea(1:ncol)
          meantau_isccp(:ncol)    = cospOUT%isccp_meantaucld(:ncol)*cospOUT%isccp_totalcldarea(1:ncol)
       end where
       call outfld('MEANCLDALB_ISCCP', meancldalb_isccp(:ncol), ncol, lchnk)
       call outfld('MEANPTOP_ISCCP'  , meanptop_isccp(:ncol)  , ncol, lchnk)
       call outfld('MEANTAU_ISCCP'   , meantau_isccp(:ncol)   , ncol, lchnk)
    end if

    ! LIDAR SIMULATOR OUTPUTS
    if (cosp_llidar_sim) then
       call outfld('CLDLOW_CAL',    cospOUT%calipso_cldlayer(:ncol,1), ncol, lchnk)
       call outfld('CLDMED_CAL',    cospOUT%calipso_cldlayer(:ncol,2), ncol, lchnk)
       call outfld('CLDHGH_CAL',    cospOUT%calipso_cldlayer(:ncol,3), ncol, lchnk)
       call outfld('CLDTOT_CAL',    cospOUT%calipso_cldlayer(:ncol,4), ncol, lchnk)
       call outfld('CLDTOT_CAL_ICE',cospOUT%calipso_cldlayerphase(:ncol,4,1), ncol, lchnk)
       call outfld('CLDTOT_CAL_LIQ',cospOUT%calipso_cldlayerphase(:ncol,4,2), ncol, lchnk)
       call outfld('CLDTOT_CAL_UN', cospOUT%calipso_cldlayerphase(:ncol,4,3), ncol, lchnk)
       call outfld('CLDHGH_CAL_ICE',cospOUT%calipso_cldlayerphase(:ncol,3,1), ncol, lchnk)
       call outfld('CLDHGH_CAL_LIQ',cospOUT%calipso_cldlayerphase(:ncol,3,2), ncol, lchnk)
       call outfld('CLDHGH_CAL_UN', cospOUT%calipso_cldlayerphase(:ncol,3,3), ncol, lchnk)
       call outfld('CLDMED_CAL_ICE',cospOUT%calipso_cldlayerphase(:ncol,2,1), ncol, lchnk)
       call outfld('CLDMED_CAL_LIQ',cospOUT%calipso_cldlayerphase(:ncol,2,2), ncol, lchnk)
       call outfld('CLDMED_CAL_UN', cospOUT%calipso_cldlayerphase(:ncol,2,3), ncol, lchnk)
       call outfld('CLDLOW_CAL_ICE',cospOUT%calipso_cldlayerphase(:ncol,1,1), ncol, lchnk)
       call outfld('CLDLOW_CAL_LIQ',cospOUT%calipso_cldlayerphase(:ncol,1,2), ncol, lchnk)
       call outfld('CLDLOW_CAL_UN', cospOUT%calipso_cldlayerphase(:ncol,1,3), ncol, lchnk)
       call outfld('MOL532_CAL', cospOUT%calipso_beta_mol(1:ncol,:), ncol, lchnk)

       ! Need to reset missing values to 0 here; CAM history does not like mix
       ! of missing values and valid values in the vertical dimension. The
       ! missing values here are likely values below sea level.

       ! TODO: reversing vertical dimension looks like a bug here
       cfad_lidarsr532(1:ncol,1:nsr_cosp,1:Nlvgrid) = cospOUT%calipso_cfad_sr(:,:,Nlvgrid:1:-1)
       where (cfad_lidarsr532(1:ncol,:,:) .eq. R_UNDEF) cfad_lidarsr532(1:ncol,:,:) = 0._r8
       call outfld('CFAD_SR532_CAL', cfad_lidarsr532(1:ncol,:,:), ncol, lchnk)

       cld_cal(1:ncol,1:Nlvgrid) = cospOUT%calipso_lidarcld(:,Nlvgrid:1:-1)
       where (cld_cal(:ncol,:Nlvgrid) .eq. R_UNDEF) cld_cal(:ncol,:Nlvgrid) = 0.0_r8
       call outfld('CLD_CAL', cld_cal(:ncol,:), ncol, lchnk)

       cld_cal_liq(1:ncol,1:Nlvgrid) = cospOUT%calipso_lidarcldphase(:,Nlvgrid:1:-1,2)
       where (cld_cal_liq(:ncol,:Nlvgrid) .eq. R_UNDEF) cld_cal_liq(:ncol,:Nlvgrid) = 0.0_r8
       call outfld('CLD_CAL_LIQ', cld_cal_liq(:ncol,:), ncol, lchnk)

       cld_cal_ice(1:ncol,1:Nlvgrid) = cospOUT%calipso_lidarcldphase(:,Nlvgrid:1:-1,1)
       where (cld_cal_ice(:ncol,:Nlvgrid) .eq. R_UNDEF) cld_cal_ice(:ncol,:Nlvgrid) = 0.0_r8
       call outfld('CLD_CAL_ICE', cld_cal_ice(:ncol,:), ncol, lchnk)

       cld_cal_un(1:ncol,1:Nlvgrid) = cospOUT%calipso_lidarcldphase(:,Nlvgrid:1:-1,3)
       where (cld_cal_un(:ncol,:Nlvgrid) .eq. R_UNDEF) cld_cal_un(:ncol,:Nlvgrid) = 0.0_r8
       call outfld('CLD_CAL_UN', cld_cal_un(:ncol,:), ncol, lchnk)

       cld_cal_tmp(1:ncol,1:Nlvgrid) = cospOUT%calipso_lidarcldtmp(:,:,1)
       where (cld_cal_tmp(:ncol,:Nlvgrid) .eq. R_UNDEF) cld_cal_tmp(:ncol,:Nlvgrid) = 0.0_r8
       call outfld('CLD_CAL_TMP', cld_cal_tmp(:ncol,:), ncol, lchnk)

       cld_cal_tmpliq(1:ncol,1:Nlvgrid) = cospOUT%calipso_lidarcldtmp(:,:,2)
       where (cld_cal_tmpliq(:ncol,:Nlvgrid) .eq. R_UNDEF) cld_cal_tmpliq(:ncol,:Nlvgrid) = 0.0_r8
       call outfld('CLD_CAL_TMPLIQ', cld_cal_tmpliq(:ncol,:), ncol, lchnk)

       cld_cal_tmpice(1:ncol,1:Nlvgrid) = cospOUT%calipso_lidarcldtmp(:,:,3)
       where (cld_cal_tmpice(:ncol,:Nlvgrid) .eq. R_UNDEF) cld_cal_tmpice(:ncol,:Nlvgrid) = 0.0_r8
       call outfld('CLD_CAL_TMPICE', cld_cal_tmpice(:ncol,:), ncol, lchnk)

       cld_cal_tmpun(1:ncol,1:Nlvgrid) = cospOUT%calipso_lidarcldtmp(:,:,4)
       where (cld_cal_tmpun(:ncol,:Nlvgrid) .eq. R_UNDEF) cld_cal_tmpun(:ncol,:Nlvgrid) = 0.0_r8
       call outfld('CLD_CAL_TMPUN', cld_cal_tmpun(:ncol,:), ncol, lchnk)
    end if

    if (cosp_lparasol_sim) then
       refl_parasol(1:ncol,1:nsza_cosp) = cospOUT%parasolGrid_refl
       where (refl_parasol(:ncol,:nsza_cosp) .eq. R_UNDEF) refl_parasol(:ncol,:nsza_cosp) = 0
       call outfld('RFL_PARASOL', refl_parasol(:ncol,:), ncol, lchnk)
    end if

    ! RADAR simulator outputs
    if (cosp_lradar_sim) then
       ! Need to reset missing values to zero here to avoid mix of missing and
       ! real values for levels below the surface
       cfad_dbze94(:ncol,:,:) = cospOUT%cloudsat_cfad_ze(:ncol,:,:)
       where (cfad_dbze94(:ncol,:,:) .eq. R_UNDEF) cfad_dbze94(:ncol,:,:) = 0.0_r8
       ! Note levels are flipped, but I think this is maybe a mistake?
       call outfld('CFAD_DBZE94_CS', cfad_dbze94(:ncol,:,Nlvgrid:1:-1), ncol, lchnk)
       call outfld('CLDTOT_CALCS',   cospOUT%radar_lidar_tcc(1:ncol), ncol, lchnk)

       ! TODO: these actually do exist in COSP2 and should be added to the
       ! output, but passing on this for the moment to keep tests BFB.
       cldtot_cs(1:ncol)  = 0  !cospOUT%cloudsat_tcc
       cldtot_cs2(1:ncol) = 0  !cospOUT%cloudsat_tcc2
       call outfld('CLDTOT_CS',     cldtot_cs(1:ncol),  ncol, lchnk)
       call outfld('CLDTOT_CS2',    cldtot_cs2(1:ncol), ncol, lchnk)
       call outfld('CLD_CAL_NOTCS', cospOUT%lidar_only_freq_cloud(1:ncol,:), ncol, lchnk)
    end if

    ! MISR simulator outputs
    if (cosp_lmisr_sim) then
       call outfld('CLD_MISR', cospOUT%misr_fq(1:ncol,:,:), ncol, lchnk)
    end if

    ! MODIS simulator outputs
    if (cosp_lmodis_sim) then

       call outfld('CLTMODIS' , cospOUT%modis_Cloud_Fraction_Total_Mean(1:ncol), ncol, lchnk)
       call outfld('CLWMODIS' , cospOUT%modis_Cloud_Fraction_Water_Mean(1:ncol), ncol, lchnk)
       call outfld('CLIMODIS' , cospOUT%modis_Cloud_Fraction_Ice_Mean  (1:ncol), ncol, lchnk)
       call outfld('CLHMODIS' , cospOUT%modis_Cloud_Fraction_High_Mean (1:ncol), ncol, lchnk)
       call outfld('CLMMODIS' , cospOUT%modis_Cloud_Fraction_Mid_Mean  (1:ncol), ncol, lchnk)
       call outfld('CLLMODIS' , cospOUT%modis_Cloud_Fraction_Low_Mean  (1:ncol), ncol, lchnk)
       call outfld('CLMODIS'  , cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(1:ncol,:,:), ncol, lchnk)
       call outfld('CLRIMODIS', cospOUT%modis_Optical_Thickness_vs_ReffICE(1:ncol,:,:), ncol, lchnk)
       call outfld('CLRLMODIS', cospOUT%modis_Optical_Thickness_vs_ReffLIQ(1:ncol,:,:), ncol, lchnk)

       ! Weight retrievals by cloud fraction
       tautmodis(1:ncol) = masked_product(cospOUT%modis_Optical_Thickness_Total_Mean, cospOUT%modis_Cloud_Fraction_Total_mean(1:ncol), R_UNDEF)
       call outfld('TAUTMODIS', tautmodis(1:ncol), ncol, lchnk)

       tauwmodis(1:ncol) = masked_product(cospOUT%modis_Optical_Thickness_Water_Mean, cospOUT%modis_Cloud_Fraction_Water_mean, R_UNDEF)
       call outfld('TAUWMODIS', tauwmodis(1:ncol), ncol, lchnk)

       tauimodis(1:ncol) = masked_product(cospOUT%modis_Optical_Thickness_Ice_Mean, cospOUT%modis_Cloud_Fraction_Ice_Mean, R_UNDEF)
       call outfld('TAUIMODIS', tauimodis(1:ncol), ncol, lchnk)

       tautlogmodis(1:ncol) = masked_product(cospOUT%modis_Optical_Thickness_Total_LogMean, cospOUT%modis_Cloud_Fraction_Total_Mean, R_UNDEF)
       call outfld('TAUTLOGMODIS', tautlogmodis(1:ncol), ncol, lchnk)

       tauwlogmodis(1:ncol) = masked_product(cospOUT%modis_Optical_Thickness_Water_LogMean, cospOUT%modis_Cloud_Fraction_Water_Mean, R_UNDEF)
       call outfld('TAUWLOGMODIS', tauwlogmodis(1:ncol), ncol, lchnk)

       tauilogmodis(1:ncol) = masked_product(cospOUT%modis_Optical_Thickness_Ice_LogMean, cospOUT%modis_Cloud_Fraction_Ice_Mean, R_UNDEF)
       call outfld('TAUILOGMODIS', tauilogmodis(1:ncol), ncol, lchnk)

       reffclwmodis(1:ncol) = masked_product(cospOUT%modis_Cloud_Particle_Size_Water_Mean, cospOUT%modis_Cloud_Fraction_Water_Mean, R_UNDEF)
       call outfld('REFFCLWMODIS', reffclwmodis(1:ncol), ncol, lchnk)

       reffclimodis(1:ncol) = masked_product(cospOUT%modis_Cloud_Particle_Size_Ice_Mean, cospOUT%modis_Cloud_Fraction_Ice_Mean, R_UNDEF)
       call outfld('REFFCLIMODIS', reffclimodis(1:ncol), ncol, lchnk)

       pctmodis(1:ncol) = masked_product(cospOUT%modis_Cloud_Top_Pressure_Total_Mean, cospOUT%modis_Cloud_Fraction_Total_Mean, R_UNDEF)
       call outfld('PCTMODIS', pctmodis(1:ncol), ncol, lchnk)

       lwpmodis(1:ncol) = masked_product(cospOUT%modis_Liquid_Water_Path_Mean, cospOUT%modis_Cloud_Fraction_Water_Mean, R_UNDEF)
       call outfld('LWPMODIS', lwpmodis(1:ncol), ncol, lchnk)

       iwpmodis(1:ncol) = masked_product(cospOUT%modis_Ice_Water_Path_Mean, cospOUT%modis_Cloud_Fraction_Ice_Mean, R_UNDEF)
       call outfld('IWPMODIS', iwpmodis(1:ncol), ncol, lchnk)

    end if

    ! SUB-COLUMN OUTPUT (fail check_accum if 'A')
    if (cosp_lfrac_out) then
       if (cosp_lisccp_sim) then
          call outfld('TAU_ISCCP'    , cospOUT%isccp_boxtau(1:ncol,:) , ncol, lchnk)
          call outfld('CLDPTOP_ISCCP', cospOUT%isccp_boxptop(1:ncol,:), ncol, lchnk)
       end if
       if (cosp_llidar_sim) then
          call outfld('ATB532_CAL', cospOUT%calipso_beta_tot(1:ncol,:,:), ncol, lchnk)
       end if
       if (cosp_lradar_sim) then
          call outfld('DBZE_CS', cospOUT%cloudsat_Ze_tot(1:ncol,:,:), ncol, lchnk)
       end if
       if (cosp_lparasol_sim) then
          call outfld('RFL_PARASOL_PIX', cospOUT%parasolPix_refl(1:ncol,:,:), ncol, lchnk)
       end if
    end if

end subroutine cosp_history_output


! Provide a function to return a product of two variables when those variables
! are masked with fill/missing values; this is used when weighting COSP output
! variables by cloud fraction
function masked_product(var1, var2, missing_value)
   
   real(r8), intent(in) :: var1(:), var2(:)
   real(r8), intent(in) :: missing_value
   real(r8) :: masked_product(size(var1))

   where ((var1 .ne. missing_value) .and. (var2 .ne. missing_value))
      masked_product = var1 * var2
   elsewhere
      masked_product = missing_value
   end where

   return

end function masked_product


  ! ######################################################################################
  ! SUBROUTINE cosp_subsample
  ! ######################################################################################
  subroutine cosp_subsample(Npoints, nLevels, Ncolumns, nHydro, overlap, &
      use_precipitation_fluxes, tca, cca, &
      gb_hydro, gb_reff, &
      dtau_s, dtau_c, dem_s, dem_c, &
      dtau_s_snow, dem_s_snow, &
      sfcP, cospstateIN, cospIN, &
      mr_hydro, Reff, Np, frac_prec)
   ! Dependencies
   use cosp_kinds,           only: wp
   use mod_rng,              only: rng_state, init_rng
   use mod_scops,            only: scops
   use mod_prec_scops,       only: prec_scops
   use mod_cosp_utils,       only: cosp_precip_mxratio
   use cosp_optics,          only: cosp_simulator_optics

   ! Inputs
   logical,intent(in) :: &
        use_precipitation_fluxes
   integer,intent(in) :: &
        Npoints,      & ! Number of gridpoints
        nLevels,      & ! Number of vertical levels
        Ncolumns,     & ! Number of subcolumns
        nHydro,       & ! Number pf hydrometeor types
        overlap         ! Overlap assumption (1/2/3)
   real(wp),intent(in),dimension(Npoints,nLevels) :: &
        tca,          & ! Total cloud amount (0-1)
        cca,          & ! Convective cloud amount (0-1)
        dtau_s,       & ! 0.67-micron optical depth (stratiform)
        dtau_c,       & ! 0.67-micron optical depth (convective)
        dem_s,        & ! 11-micron emissivity (stratiform)
        dem_c           ! 11-micron emissivity (convective)
   real(wp),intent(inout),dimension(Npoints,nLevels) :: &
        dtau_s_snow,  & ! 0.67-micron optical depth (snow)
        dem_s_snow      ! 11-micron emissivity (snow)
   real(wp),intent(in),dimension(Npoints,nLevels,nHydro) :: &
        gb_hydro,     & ! Gridbox-mean mixing ratios and precip fluxes (kg/kg)
        gb_reff         ! Gridbox-mean effective radii
   real(wp),intent(in),dimension(Npoints) :: &
        sfcP            ! Surface pressure

   ! Outputs
   type(cosp_optical_inputs),intent(inout) :: cospIN
   type(cosp_column_inputs),intent(in)  :: cospstateIN
   real(wp), dimension(Npoints, Ncolumns, nLevels, nHydro), intent(out) :: &
      mr_hydro, Reff, Np
   real(wp), dimension(Npoints, Ncolumns, nLevels), intent(out)  :: frac_prec

   ! Local variables
   integer :: i,j,k
   real(wp),dimension(Npoints,nLevels)      :: column_frac_out,column_prec_out,         &
                                               fl_lsrain,fl_lssnow,fl_lsgrpl,fl_ccrain, &
                                               fl_ccsnow
   type(rng_state),allocatable,dimension(:) :: rngs  ! Seeds for random number generator
   integer,dimension(:),allocatable         :: seed
   real(wp),dimension(:,:),allocatable      :: ls_p_rate,cv_p_rate,frac_ls,frac_cv,     &
                                               prec_ls,prec_cv

   call t_startf("scops")
   if (Ncolumns .gt. 1) then
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Generate subcolumns for clouds (SCOPS) and precipitation type (PREC_SCOPS)
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! RNG used for subcolumn generation
      allocate(rngs(Npoints), seed(Npoints))
      seed = int(sfcP)
      if (Npoints .gt. 1) seed=(sfcP-int(sfcP))*1000000
      call init_rng(rngs, seed)

      ! Call scops
      call scops(Npoints,Nlevels,Ncolumns,rngs,tca,cca,overlap,cospIN%frac_out,0)
      deallocate(seed,rngs)

      ! Sum up precipitation rates. If not using precipitation fluxes, mixing ratios are
      ! stored in _rate variables.
      allocate(ls_p_rate(Npoints,nLevels),cv_p_rate(Npoints,Nlevels))
      ls_p_rate(:,1:nLevels) = gb_hydro(:,:,I_LSRAIN) + gb_hydro(:,:,I_LSSNOW) + gb_hydro(:,:,I_LSGRPL)
      cv_p_rate(:,1:nLevels) = gb_hydro(:,:,I_CVRAIN) + gb_hydro(:,:,I_CVSNOW)

      ! Call PREC_SCOPS
      call prec_scops(Npoints,nLevels,Ncolumns,ls_p_rate,cv_p_rate,cospIN%frac_out,frac_prec)
      deallocate(ls_p_rate,cv_p_rate)

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Compute precipitation fraction in each gridbox
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Allocate
      allocate(frac_ls(Npoints,nLevels),prec_ls(Npoints,nLevels),                       &
               frac_cv(Npoints,nLevels),prec_cv(Npoints,nLevels))

      ! Initialize
      frac_ls(1:Npoints,1:nLevels) = 0._wp
      prec_ls(1:Npoints,1:nLevels) = 0._wp
      frac_cv(1:Npoints,1:nLevels) = 0._wp
      prec_cv(1:Npoints,1:nLevels) = 0._wp
      do j=1,Npoints
         do k=1,nLevels
            do i=1,Ncolumns
               if (cospIN%frac_out(j,i,k)  .eq. 1)  frac_ls(j,k) = frac_ls(j,k)+1._wp
               if (cospIN%frac_out(j,i,k)  .eq. 2)  frac_cv(j,k) = frac_cv(j,k)+1._wp
               if (frac_prec(j,i,k) .eq. 1)         prec_ls(j,k) = prec_ls(j,k)+1._wp
               if (frac_prec(j,i,k) .eq. 2)         prec_cv(j,k) = prec_cv(j,k)+1._wp
               if (frac_prec(j,i,k) .eq. 3)         prec_cv(j,k) = prec_cv(j,k)+1._wp
               if (frac_prec(j,i,k) .eq. 3)         prec_ls(j,k) = prec_ls(j,k)+1._wp
            enddo
            frac_ls(j,k)=frac_ls(j,k)/Ncolumns
            frac_cv(j,k)=frac_cv(j,k)/Ncolumns
            prec_ls(j,k)=prec_ls(j,k)/Ncolumns
            prec_cv(j,k)=prec_cv(j,k)/Ncolumns

            ! Adjust grid-box mean snow properties to local properties
            ! Convert longwave optical depth to longwave emissivity
            if (prec_ls(j,k) .ne. 0._r8 .and. dtau_s_snow(j,k) .gt. 0._r8) then
               dtau_s_snow(j,k) = dtau_s_snow(j,k)/prec_ls(j,k)
            end if
            if (prec_ls(j,k) .ne. 0._r8 .and. dem_s_snow(j,k) .gt. 0._r8) then
               dem_s_snow(j,k) = dem_s_snow(j,k)/prec_ls(j,k)
               dem_s_snow(j,k) = 1._r8 - exp ( -1._r8*dem_s_snow(j,k))
            end if !!+JEK
         enddo
      enddo

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Compute mixing ratios, effective radii and precipitation fluxes for clouds
      ! and precipitation
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Initialize
      mr_hydro(:,:,:,:) = 0._wp
      Reff(:,:,:,:)     = 0._wp
      Np(:,:,:,:)       = 0._wp

      do k=1,Ncolumns
         ! Subcolumn clouds
         column_frac_out = cospIN%frac_out(:,k,:)

         ! LS clouds
         where (column_frac_out == I_LSC)
            mr_hydro(:,k,:,I_LSCLIQ) = gb_hydro(:,:,I_LSCLIQ)
            mr_hydro(:,k,:,I_LSCICE) = gb_hydro(:,:,I_LSCICE)
            Reff(:,k,:,I_LSCLIQ)     = gb_reff(:,:,I_LSCLIQ)
            Reff(:,k,:,I_LSCICE)     = gb_reff(:,:,I_LSCICE)
         ! CONV clouds
         elsewhere (column_frac_out == I_CVC)
            mr_hydro(:,k,:,I_CVCLIQ) = gb_hydro(:,:,I_CVCLIQ)
            mr_hydro(:,k,:,I_CVCICE) = gb_hydro(:,:,I_CVCICE)
            Reff(:,k,:,I_CVCLIQ)     = gb_reff(:,:,I_CVCLIQ)
            Reff(:,k,:,I_CVCICE)     = gb_reff(:,:,I_CVCICE)
         end where

         ! Subcolumn precipitation
         column_prec_out = frac_prec(:,k,:)

         ! LS Precipitation
         where ((column_prec_out == 1) .or. (column_prec_out == 3) )
            Reff(:,k,:,I_LSRAIN) = gb_reff(:,:,I_LSRAIN)
            Reff(:,k,:,I_LSSNOW) = gb_reff(:,:,I_LSSNOW)
            Reff(:,k,:,I_LSGRPL) = gb_reff(:,:,I_LSGRPL)
         ! CONV precipitation
         elsewhere ((column_prec_out == 2) .or. (column_prec_out == 3))
            Reff(:,k,:,I_CVRAIN) = gb_reff(:,:,I_CVRAIN)
            Reff(:,k,:,I_CVSNOW) = gb_reff(:,:,I_CVSNOW)
         end where
      enddo

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Convert the mixing ratio and precipitation fluxes from gridbox mean to
      ! the fraction-based values
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      do k=1,nLevels
         do j=1,Npoints
            ! Clouds
            if (frac_ls(j,k) .ne. 0._r8) then
               mr_hydro(j,:,k,I_LSCLIQ) = mr_hydro(j,:,k,I_LSCLIQ)/frac_ls(j,k)
               mr_hydro(j,:,k,I_LSCICE) = mr_hydro(j,:,k,I_LSCICE)/frac_ls(j,k)
            endif
            if (frac_cv(j,k) .ne. 0._r8) then
               mr_hydro(j,:,k,I_CVCLIQ) = mr_hydro(j,:,k,I_CVCLIQ)/frac_cv(j,k)
               mr_hydro(j,:,k,I_CVCICE) = mr_hydro(j,:,k,I_CVCICE)/frac_cv(j,k)
            endif

            ! Precipitation
            if (use_precipitation_fluxes) then
               if (prec_ls(j,k) .ne. 0._r8) then
                  fl_lsrain(j,k) = gb_hydro(j,k,I_LSRAIN)/prec_ls(j,k)
                  fl_lssnow(j,k) = gb_hydro(j,k,I_LSSNOW)/prec_ls(j,k)
                  fl_lsgrpl(j,k) = gb_hydro(j,k,I_LSGRPL)/prec_ls(j,k)
               endif
               if (prec_cv(j,k) .ne. 0._r8) then
                  fl_ccrain(j,k) = gb_hydro(j,k,I_CVRAIN)/prec_cv(j,k)
                  fl_ccsnow(j,k) = gb_hydro(j,k,I_CVSNOW)/prec_cv(j,k)
               endif
            else
               if (prec_ls(j,k) .ne. 0._r8) then
                  mr_hydro(j,:,k,I_LSRAIN) = mr_hydro(j,:,k,I_LSRAIN)/prec_ls(j,k)
                  mr_hydro(j,:,k,I_LSSNOW) = mr_hydro(j,:,k,I_LSSNOW)/prec_ls(j,k)
                  mr_hydro(j,:,k,I_LSGRPL) = mr_hydro(j,:,k,I_LSGRPL)/prec_ls(j,k)
               endif
               if (prec_cv(j,k) .ne. 0._r8) then
                  mr_hydro(j,:,k,I_CVRAIN) = mr_hydro(j,:,k,I_CVRAIN)/prec_cv(j,k)
                  mr_hydro(j,:,k,I_CVSNOW) = mr_hydro(j,:,k,I_CVSNOW)/prec_cv(j,k)
               endif
            endif
         enddo
      enddo

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Convert precipitation fluxes to mixing ratios
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (use_precipitation_fluxes) then
         ! LS rain
         call cosp_precip_mxratio(Npoints, nLevels, Ncolumns, cospstateIN%pfull,        &
              cospstateIN%at, frac_prec, 1._wp, n_ax(I_LSRAIN), n_bx(I_LSRAIN),         &
              alpha_x(I_LSRAIN), c_x(I_LSRAIN),   d_x(I_LSRAIN),   g_x(I_LSRAIN),       &
              a_x(I_LSRAIN),   b_x(I_LSRAIN),   gamma_1(I_LSRAIN), gamma_2(I_LSRAIN),   &
              gamma_3(I_LSRAIN), gamma_4(I_LSRAIN), fl_lsrain,                          &
              mr_hydro(:,:,:,I_LSRAIN), Reff(:,:,:,I_LSRAIN))
         ! LS snow
         call cosp_precip_mxratio(Npoints, nLevels, Ncolumns, cospstateIN%pfull,        &
              cospstateIN%at, frac_prec, 1._wp,  n_ax(I_LSSNOW),  n_bx(I_LSSNOW),       &
              alpha_x(I_LSSNOW), c_x(I_LSSNOW),  d_x(I_LSSNOW),  g_x(I_LSSNOW),         &
              a_x(I_LSSNOW),   b_x(I_LSSNOW),   gamma_1(I_LSSNOW),  gamma_2(I_LSSNOW),  &
              gamma_3(I_LSSNOW), gamma_4(I_LSSNOW), fl_lssnow,                          &
              mr_hydro(:,:,:,I_LSSNOW), Reff(:,:,:,I_LSSNOW))
         ! CV rain
         call cosp_precip_mxratio(Npoints, nLevels, Ncolumns, cospstateIN%pfull,        &
              cospstateIN%at, frac_prec, 2._wp, n_ax(I_CVRAIN),  n_bx(I_CVRAIN),        &
              alpha_x(I_CVRAIN), c_x(I_CVRAIN),   d_x(I_CVRAIN),   g_x(I_CVRAIN),       &
              a_x(I_CVRAIN),   b_x(I_CVRAIN),   gamma_1(I_CVRAIN), gamma_2(I_CVRAIN),   &
              gamma_3(I_CVRAIN), gamma_4(I_CVRAIN), fl_ccrain,                          &
              mr_hydro(:,:,:,I_CVRAIN), Reff(:,:,:,I_CVRAIN))
         ! CV snow
         call cosp_precip_mxratio(Npoints, nLevels, Ncolumns, cospstateIN%pfull,        &
              cospstateIN%at, frac_prec, 2._wp, n_ax(I_CVSNOW),  n_bx(I_CVSNOW),        &
              alpha_x(I_CVSNOW),  c_x(I_CVSNOW),   d_x(I_CVSNOW),   g_x(I_CVSNOW),      &
              a_x(I_CVSNOW),   b_x(I_CVSNOW),   gamma_1(I_CVSNOW), gamma_2(I_CVSNOW),   &
              gamma_3(I_CVSNOW), gamma_4(I_CVSNOW), fl_ccsnow,                          &
              mr_hydro(:,:,:,I_CVSNOW), Reff(:,:,:,I_CVSNOW))
         ! LS groupel.
         call cosp_precip_mxratio(Npoints, nLevels, Ncolumns, cospstateIN%pfull,        &
              cospstateIN%at, frac_prec, 1._wp, n_ax(I_LSGRPL),  n_bx(I_LSGRPL),        &
              alpha_x(I_LSGRPL), c_x(I_LSGRPL),   d_x(I_LSGRPL),   g_x(I_LSGRPL),       &
              a_x(I_LSGRPL),   b_x(I_LSGRPL),   gamma_1(I_LSGRPL),  gamma_2(I_LSGRPL),  &
              gamma_3(I_LSGRPL), gamma_4(I_LSGRPL), fl_lsgrpl,                          &
              mr_hydro(:,:,:,I_LSGRPL), Reff(:,:,:,I_LSGRPL))
      endif

      ! Subsample input optical properties
      if (allocated(cospIN%emiss_11)) then
         call cosp_simulator_optics(Npoints,nSubcol,nLevels,cospIN%frac_out,dem_c,dem_s,  &
              cospIN%emiss_11)
      end if
      if (allocated(cospIN%tau_067)) then
         call cosp_simulator_optics(Npoints,nSubcol,nLevels,cospIN%frac_out,dtau_c,dtau_s,  &
              cospIN%tau_067)
      end if

   else
      cospIN%frac_out(:,1,:) = 1
      cospIN%emiss_11(:,1,:) = max(dem_s, dem_c)
      cospIN%tau_067 (:,1,:) = max(dtau_s, dtau_c)
      mr_hydro(:,1,:,I_LSCLIQ) = gb_hydro(:,:,I_LSCLIQ)
      mr_hydro(:,1,:,I_LSCICE) = gb_hydro(:,:,I_LSCICE)
      mr_hydro(:,1,:,I_CVCLIQ) = gb_hydro(:,:,I_CVCLIQ)
      mr_hydro(:,1,:,I_CVCICE) = gb_hydro(:,:,I_CVCICE)
      Reff(:,1,:,:)            = gb_reff(:,:,:)
   endif
   call t_stopf("scops")

  end subroutine cosp_subsample

  ! ######################################################################################
  ! SUBROUTINE calc_cosp_optics
  ! ######################################################################################
  subroutine calc_cosp_optics(Npoints, nLevels, nSubcol, nHydro, &
      lidar_ice_type, sd, &
      mr_hydro, Reff, Np, frac_prec, &
      cospstateIN, cospIN, &
      dtau_snow, dem_snow)
   ! Dependencies
   use cosp_kinds, only: wp
   use mod_cosp_config,      only: R_UNDEF
   use mod_quickbeam_optics, only: quickbeam_optics, gases
   use cosp_optics,          only: cosp_simulator_optics, lidar_optics, modis_optics, &
                                   modis_optics_partition
   integer, intent(in) :: &
      Npoints,    & ! Number of gridpoints
      nLevels,    & ! Number of vertical levels
      nSubcol,   & ! Number of subcolumns
      nHydro,     & ! Number pf hydrometeor types
      lidar_ice_type  ! Ice type assumption used by lidar optics
   real(wp), intent(in), dimension(Npoints,nLevels), optional :: &
      dtau_snow,  & ! 0.67-micron optical depth (snow)
      dem_snow      ! 11-micron emissivity (snow)
   real(wp), dimension(Npoints, nSubcol, nLevels, nHydro), intent(inout) :: &
      mr_hydro, Reff, Np
   real(wp), dimension(Npoints, nSubcol, nLevels), intent(in)  :: frac_prec

   type(size_distribution),intent(inout) :: &
      sd
   ! Outputs
   type(cosp_column_inputs),intent(inout)  :: cospstateIN
   type(cosp_optical_inputs),intent(inout) :: cospIN

   ! Local variables
   integer :: i, j, k
   real(wp), dimension(:, :), allocatable :: g_vol
   real(wp),dimension(Npoints, nLevels, nHydro) :: ReffTemp
   logical, parameter :: lground = .false.
   integer, parameter :: lidar_freq = 532
   real(wp),dimension(:,:,:),  allocatable  :: &
                                                MODIS_cloudWater,MODIS_cloudIce,        &
                                                MODIS_watersize,MODIS_iceSize,          &
                                                MODIS_snowSize,MODIS_cloudSnow,         &
                                                MODIS_opticalThicknessLiq,              &
                                                MODIS_opticalThicknessSnow,             &
                                                MODIS_opticalThicknessIce

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! CLOUDSAT RADAR OPTICS
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   call t_startf("cloudsat_optics")
   if (cosp_lradar_sim) then
      ! Compute gaseous absorption (assume identical for each subcolun)
      allocate(g_vol(Npoints,nLevels))
      g_vol(:,:)=0._wp
      do i = 1, Npoints
         do j = 1, nLevels
            if (cospIN%rcfg_cloudsat%use_gas_abs == 1 .or. &
               (cospIN%rcfg_cloudsat%use_gas_abs == 2 .and. j == 1)) then
               g_vol(i,j) = gases(cospstateIN%pfull(i,j), cospstateIN%at(i,j),    &
                                  cospstateIN%qv(i,j), cospIN%rcfg_cloudsat%freq)
            endif
            cospIN%g_vol_cloudsat(i,:,j) = g_vol(i,j)
         end do
      end do

      ! Loop over all subcolumns
      do k=1,nSubcol
         call quickbeam_optics(sd, cospIN%rcfg_cloudsat, Npoints, nLevels, R_UNDEF, &
              mr_hydro(:,k,:,1:nHydro)*1000._wp, Reff(:,k,:,1:nHydro)*1.e6_wp,      &
              Np(:,k,:,1:nHydro), cospstateIN%pfull, cospstateIN%at,                &
              cospstateIN%qv, cospIN%z_vol_cloudsat(1:Npoints,k,:),                 &
              cospIN%kr_vol_cloudsat(1:Npoints,k,:))
      enddo
   endif
   call t_stopf("cloudsat_optics")

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! LIDAR Polarized optics
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   call t_startf("calipso_optics")
   if (cosp_llidar_sim) then
      ! lidar_optics wants gridbox effective radii instead of subcolumn, so back
      ! this out there by taking the max value along the subcolumn dimension so
      ! that we do not have to pass a redundant variable into this routine and
      ! also so that we can use with data only definied on subcolumns (i.e.,
      ! SP/MMF)
      ReffTemp = maxval(Reff, dim=2)
      call lidar_optics(Npoints,nSubcol,nLevels,5,lidar_ice_type,lidar_freq,lground,   &
                        mr_hydro(1:Npoints,1:nSubcol,1:nLevels,I_LSCLIQ),              &
                        mr_hydro(1:Npoints,1:nSubcol,1:nLevels,I_LSCICE),              &
                        mr_hydro(1:Npoints,1:nSubcol,1:nLevels,I_CVCLIQ),              &
                        mr_hydro(1:Npoints,1:nSubcol,1:nLevels,I_CVCICE),              &
                        ReffTemp(1:Npoints,1:nLevels,I_LSCLIQ),                        &
                        ReffTemp(1:Npoints,1:nLevels,I_LSCICE),                        &
                        ReffTemp(1:Npoints,1:nLevels,I_CVCLIQ),                        &
                        ReffTemp(1:Npoints,1:nLevels,I_CVCICE),                        &
                        cospstateIN%pfull(1:Npoints,1:nLevels),                        &
                        cospstateIN%phalf(1:Npoints,1:nLevels+1),                      &
                        cospstateIN%at(1:Npoints,1:nLevels),                           &
                        cospIN%beta_mol_calipso(1:Npoints,1:nLevels),                  &
                        cospIN%betatot_calipso(1:Npoints,1:nSubcol,1:nLevels),         &
                        cospIN%tau_mol_calipso(1:Npoints,1:nLevels),                   &
                        cospIN%tautot_calipso(1:Npoints,1:nSubcol,1:nLevels),          &
                        cospIN%tautot_S_liq(1:Npoints,1:nSubcol),                      &
                        cospIN%tautot_S_ice(1:Npoints,1:nSubcol),                      &
                        cospIN%betatot_ice_calipso(1:Npoints,1:nSubcol,1:nLevels),     &
                        cospIN%betatot_liq_calipso(1:Npoints,1:nSubcol,1:nLevels),     &
                        cospIN%tautot_ice_calipso(1:Npoints,1:nSubcol,1:nLevels),      &
                        cospIN%tautot_liq_calipso(1:Npoints,1:nSubcol,1:nLevels))
   endif
   call t_stopf("calipso_optics")


   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! Compute optical fields for passive simulators (i.e. only sunlit points)
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! 11 micron emissivity (needed by the ISCCP simulator)
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! Add in contributions from snow
   call t_startf("11micron_emissivity")
   if (allocated(cospIN%emiss_11)) then
      if (present(dem_snow)) then
         do j = 1,nSubcol
            where(frac_prec(:,j,:) .eq. 1 .or. frac_prec(:,j,:) .eq. 3)
               cospIN%emiss_11(:,j,:) = 1._wp - (1- cospIN%emiss_11(:,j,:))*(1-dem_snow(:,:))
            endwhere
         end do
      end if
   endif
   call t_stopf("11micron_emissivity")

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! 0.67 micron optical depth (needed by ISCCP, MISR and MODIS simulators)
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! Add in contributions from snow
   call t_startf("067tau")
   if (allocated(cospIN%tau_067)) then
      if (present(dtau_snow)) then
         do j = 1,nSubcol
            where((frac_prec(:,j,:) .eq. 1 .or. frac_prec(:,j,:) .eq. 3) .and. &
                 Reff(:,j,:,I_LSSNOW) .gt. 0._r8 .and. dtau_snow(:,:) .gt. 0._r8)
               cospIN%tau_067(:,j,:)  = cospIN%tau_067(:,j,:)+dtau_snow(:,:)
            endwhere
         end do
      end if
   endif
   call t_stopf("067tau")

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! MODIS optics
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   call t_startf("modis_optics")
   if (cosp_lmodis_sim) then
      allocate(MODIS_cloudWater(Npoints,nSubcol,nLevels),                              &
               MODIS_cloudIce(Npoints,nSubcol,nLevels),                                &
               MODIS_cloudSnow(Npoints,nSubcol,nLevels),                               &
               MODIS_waterSize(Npoints,nSubcol,nLevels),                               &
               MODIS_iceSize(Npoints,nSubcol,nLevels),                                 &
               MODIS_snowSize(Npoints,nSubcol,nLevels),                                &
               MODIS_opticalThicknessLiq(Npoints,nSubcol,nLevels),                     &
               MODIS_opticalThicknessIce(Npoints,nSubcol,nLevels),                     &
               MODIS_opticalThicknessSnow(Npoints,nSubcol,nLevels))

      ! Cloud water
      call cosp_simulator_optics(Npoints,nSubcol,nLevels,cospIN%frac_out,              &
           mr_hydro(:,:,:,I_CVCLIQ),mr_hydro(:,:,:,I_LSCLIQ),MODIS_cloudWater)
      ! Cloud ice
      call cosp_simulator_optics(Npoints,nSubcol,nLevels,cospIN%frac_out,              &
           mr_hydro(:,:,:,I_CVCICE),mr_hydro(:,:,:,I_LSCICE),MODIS_cloudIce)
      ! Cloud water droplet size
      call cosp_simulator_optics(Npoints,nSubcol,nLevels,cospIN%frac_out,              &
           Reff(:,:,:,I_CVCLIQ),Reff(:,:,:,I_LSCLIQ),MODIS_waterSize)
      ! Cloud ice crystal size
      call cosp_simulator_optics(Npoints,nSubcol,nLevels,cospIN%frac_out,              &
           Reff(:,:,:,I_CVCICE),Reff(:,:,:,I_LSCICE),MODIS_iceSize)

!     ! Cloud snow and size; note that these do not appear to be used!
!     ! TODO: do we need to figure out how to incorporate these?
!     if (present(dtau_snow)) then
!        do j=1,nSubcol
!           where((frac_prec(:,j,:) .eq. 1 .or. frac_prec(:,j,:) .eq. 3) .and. &
!                Reff(:,j,:,I_LSSNOW) .gt. 0._r8 .and. dtau_snow(:,:) .gt. 0._r8)
!              MODIS_cloudSnow(:,j,:) = mr_hydro(:,j,:,I_LSSNOW)
!              MODIS_snowSize(:,j,:)  = Reff(:,j,:,I_LSSNOW)
!           elsewhere
!              MODIS_snowSize(:,j,:)  = 0._wp
!              MODIS_cloudSnow(:,j,:) = 0._wp
!           endwhere
!        enddo
!     else
!        MODIS_snowSize(:,:,:)  = Reff(:,:,:,I_LSSNOW)
!        MODIS_cloudSnow(:,:,:) = mr_hydro(:,:,:,I_LSSNOW)
!     end if
      call modis_optics_partition(Npoints, nLevels, nSubcol, MODIS_cloudWater,     &
           MODIS_cloudIce, MODIS_waterSize, MODIS_iceSize,         &
           cospIN%tau_067, MODIS_opticalThicknessLiq,               &
           MODIS_opticalThicknessIce)
      call modis_optics(Npoints, nLevels, nSubcol, MODIS_opticalThicknessLiq,      &
           MODIS_waterSize*1.0e6_wp, MODIS_opticalThicknessIce,                     &
           MODIS_iceSize*1.0e6_wp, cospIN%fracLiq, cospIN%asym, cospIN%ss_alb)
   endif ! MODIS simulator optics
   call t_stopf("modis_optics")
  end subroutine calc_cosp_optics


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE construct_cospIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine construct_cospIN(Npoints,Ncolumns,Nlevels,y)
    ! Inputs
    integer,intent(in) :: &
         Npoints,  & ! Number of horizontal gridpoints
         Ncolumns, & ! Number of subcolumns
         Nlevels     ! Number of vertical levels
    ! Outputs
    type(cosp_optical_inputs),intent(out) :: y

    ! Dimensions
    y%Npoints  = Npoints
    y%Ncolumns = Ncolumns
    y%Nlevels  = Nlevels
    y%Npart    = 4
    y%Nrefl    = PARASOL_NREFL
    allocate(y%frac_out(Npoints,       Ncolumns,Nlevels))

    if (cosp_lmodis_sim .or. cosp_lmisr_sim .or. cosp_lisccp_sim) then
       allocate(y%tau_067(Npoints,        Ncolumns,Nlevels),&
                y%emiss_11(Npoints,       Ncolumns,Nlevels))
    endif
    if (cosp_llidar_sim) then
       allocate(y%betatot_calipso(Npoints,        Ncolumns,Nlevels),&
                y%betatot_ice_calipso(Npoints,    Ncolumns,Nlevels),&
                y%betatot_liq_calipso(Npoints,    Ncolumns,Nlevels),&
                y%tautot_calipso(Npoints,         Ncolumns,Nlevels),&
                y%tautot_ice_calipso(Npoints,     Ncolumns,Nlevels),&
                y%tautot_liq_calipso(Npoints,     Ncolumns,Nlevels),&
                y%beta_mol_calipso(Npoints,                Nlevels),&
                y%tau_mol_calipso(Npoints,                 Nlevels),&
                y%tautot_S_ice(Npoints,   Ncolumns        ),&
                y%tautot_S_liq(Npoints,   Ncolumns        ))
    endif

    if (cosp_lgrlidar_sim) then
       allocate(y%beta_mol_grLidar532(Npoints,          Nlevels),&
                y%betatot_grLidar532(Npoints,  Ncolumns,Nlevels),&
                y%tau_mol_grLidar532(Npoints,           Nlevels),&
                y%tautot_grLidar532(Npoints,   Ncolumns,Nlevels))
    endif

    if (cosp_latlid_sim) then
       allocate(y%beta_mol_atlid(Npoints,             Nlevels),&
                y%betatot_atlid(Npoints,     Ncolumns,Nlevels),&
                y%tau_mol_atlid(Npoints,              Nlevels),&
                y%tautot_atlid(Npoints,      Ncolumns,Nlevels))
    endif

    if (cosp_lradar_sim) then
       allocate(y%z_vol_cloudsat(Npoints,  Ncolumns,Nlevels),&
                y%kr_vol_cloudsat(Npoints, Ncolumns,Nlevels),&
                y%g_vol_cloudsat(Npoints,  Ncolumns,Nlevels),&
                y%fracPrecipIce(Npoints,   Ncolumns))
    endif
    if (cosp_lmodis_sim) then
       allocate(y%fracLiq(Npoints,        Ncolumns,Nlevels),&
                y%asym(Npoints,           Ncolumns,Nlevels),&
                y%ss_alb(Npoints,         Ncolumns,Nlevels))
    endif

  end subroutine construct_cospIN

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE construct_cospstateIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine construct_cospstateIN(Npoints,Nlevels,nchan,y)
    ! Inputs
    integer,intent(in) :: &
         Npoints, & ! Number of horizontal gridpoints
         Nlevels, & ! Number of vertical levels
         nchan      ! Number of channels
    ! Outputs
    type(cosp_column_inputs),intent(out) :: y

    allocate(y%sunlit(Npoints),y%skt(Npoints),y%land(Npoints),y%at(Npoints,Nlevels),     &
             y%pfull(Npoints,Nlevels),y%phalf(Npoints,Nlevels+1),y%qv(Npoints,Nlevels),  &
             y%o3(Npoints,Nlevels),y%hgt_matrix(Npoints,Nlevels),y%u_sfc(Npoints),       &
             y%v_sfc(Npoints),y%lat(Npoints),y%lon(Npoints),y%emis_sfc(nchan),           &
             y%cloudIce(Npoints,Nlevels),y%cloudLiq(Npoints,Nlevels),y%surfelev(Npoints),&
             y%fl_snow(Npoints,Nlevels),y%fl_rain(Npoints,Nlevels),y%seaice(Npoints),    &
             y%tca(Npoints,Nlevels),y%hgt_matrix_half(Npoints,Nlevels+1))

  end subroutine construct_cospstateIN
  ! ######################################################################################
  ! SUBROUTINE construct_cosp_outputs
  !
  ! This subroutine allocates output fields based on input logical flag switches.
  ! ######################################################################################
  subroutine construct_cosp_outputs(Npoints,Ncolumns,Nlevels,Nlvstat,Nchan,x)
    ! Inputs
    integer,intent(in) :: &
         Npoints,         & ! Number of sampled points
         Ncolumns,        & ! Number of subgrid columns
         Nlevels,         & ! Number of model levels
         Nlvstat,         & ! Number of levels in L3 stats computation
         Nchan              ! Number of RTTOV channels

    ! Outputs
    type(cosp_outputs),intent(out) :: &
         x           ! COSP output structure

     ! ISCCP simulator outputs
    if (cosp_lisccp_sim) then
       allocate(x%isccp_boxtau(Npoints,Ncolumns))
       allocate(x%isccp_boxptop(Npoints,Ncolumns))
       allocate(x%isccp_fq(Npoints,numISCCPTauBins,numISCCPPresBins))
       allocate(x%isccp_totalcldarea(Npoints))
       allocate(x%isccp_meanptop(Npoints))
       allocate(x%isccp_meantaucld(Npoints))
       allocate(x%isccp_meantb(Npoints))
       allocate(x%isccp_meantbclr(Npoints))
       allocate(x%isccp_meanalbedocld(Npoints))
    endif

    ! MISR simulator
    if (cosp_lmisr_sim) then
       allocate(x%misr_fq(Npoints,numMISRTauBins,numMISRHgtBins))
       ! *NOTE* These 3 fields are not output, but were part of the v1.4.0 cosp_misr, so
       !        they are still computed. Should probably have a logical to control these
       !        outputs.
       allocate(x%misr_dist_model_layertops(Npoints,numMISRHgtBins))
       allocate(x%misr_meanztop(Npoints))
       allocate(x%misr_cldarea(Npoints))
    endif

    ! MODIS simulator
    if (cosp_lmodis_sim) then
       allocate(x%modis_Cloud_Fraction_Total_Mean(Npoints))
       allocate(x%modis_Cloud_Fraction_Water_Mean(Npoints))
       allocate(x%modis_Cloud_Fraction_Ice_Mean(Npoints))
       allocate(x%modis_Cloud_Fraction_High_Mean(Npoints))
       allocate(x%modis_Cloud_Fraction_Mid_Mean(Npoints))
       allocate(x%modis_Cloud_Fraction_Low_Mean(Npoints))
       allocate(x%modis_Optical_Thickness_Total_Mean(Npoints))
       allocate(x%modis_Optical_Thickness_Water_Mean(Npoints))
       allocate(x%modis_Optical_Thickness_Ice_Mean(Npoints))
       allocate(x%modis_Optical_Thickness_Total_LogMean(Npoints))
       allocate(x%modis_Optical_Thickness_Water_LogMean(Npoints))
       allocate(x%modis_Optical_Thickness_Ice_LogMean(Npoints))
       allocate(x%modis_Cloud_Particle_Size_Water_Mean(Npoints))
       allocate(x%modis_Cloud_Particle_Size_Ice_Mean(Npoints))
       allocate(x%modis_Cloud_Top_Pressure_Total_Mean(Npoints))
       allocate(x%modis_Liquid_Water_Path_Mean(Npoints))
       allocate(x%modis_Ice_Water_Path_Mean(Npoints))
       allocate(x%modis_Optical_Thickness_vs_Cloud_Top_Pressure(Npoints,numModisTauBins,numMODISPresBins))
       allocate(x%modis_Optical_thickness_vs_ReffLIQ(Npoints,numMODISTauBins,numMODISReffLiqBins))
       allocate(x%modis_Optical_Thickness_vs_ReffICE(Npoints,numMODISTauBins,numMODISReffIceBins))
    endif

    ! LIDAR simulator
    if (cosp_llidar_sim) then
       allocate(x%calipso_beta_mol(Npoints,Nlevels))
       allocate(x%calipso_beta_tot(Npoints,Ncolumns,Nlevels))
       allocate(x%calipso_srbval(SR_BINS+1))
       allocate(x%calipso_cfad_sr(Npoints,SR_BINS,Nlvstat))
       allocate(x%calipso_betaperp_tot(Npoints,Ncolumns,Nlevels))
       allocate(x%calipso_lidarcld(Npoints,Nlvstat))
       allocate(x%calipso_cldlayer(Npoints,LIDAR_NCAT))
       allocate(x%calipso_lidarcldphase(Npoints,Nlvstat,6))
       allocate(x%calipso_lidarcldtmp(Npoints,LIDAR_NTEMP,5))
       allocate(x%calipso_cldlayerphase(Npoints,LIDAR_NCAT,6))
       allocate(x%calipso_cldtype(Npoints,LIDAR_NTYPE))
       allocate(x%calipso_cldtypetemp(Npoints,LIDAR_NTYPE))
       allocate(x%calipso_cldtypemeanz(Npoints,2))
       allocate(x%calipso_cldtypemeanzse(Npoints,3))
       allocate(x%calipso_cldthinemis(Npoints))
       allocate(x%calipso_lidarcldtype(Npoints,Nlvstat,LIDAR_NTYPE+1))
       ! These 2 outputs are part of the calipso output type, but are not controlled by an
       ! logical switch in the output namelist, so if all other fields are on, then allocate
       allocate(x%calipso_tau_tot(Npoints,Ncolumns,Nlevels))
       allocate(x%calipso_temp_tot(Npoints,Nlevels))
    end if

    ! GROUND LIDAR @ 532NM simulator
    if (cosp_lgrlidar_sim) then
       allocate(x%grLidar532_beta_mol(Npoints,Nlevels))
       allocate(x%grLidar532_beta_tot(Npoints,Ncolumns,Nlevels))
       allocate(x%grLidar532_srbval(SR_BINS+1))
       allocate(x%grLidar532_cfad_sr(Npoints,SR_BINS,Nlvstat))
       allocate(x%grLidar532_lidarcld(Npoints,Nlvstat))
       allocate(x%grLidar532_cldlayer(Npoints,LIDAR_NCAT))
    end if

    ! ATLID simulator
    if (cosp_latlid_sim) then
       allocate(x%atlid_beta_mol(Npoints,Nlevels))
       allocate(x%atlid_beta_tot(Npoints,Ncolumns,Nlevels))
       allocate(x%atlid_srbval(SR_BINS+1))
       allocate(x%atlid_cfad_sr(Npoints,SR_BINS,Nlvstat))
       allocate(x%atlid_cldlayer(Npoints,LIDAR_NCAT))
    endif

    ! PARASOL
    if (cosp_lparasol_sim) then
       allocate(x%parasolPix_refl(Npoints,Ncolumns,PARASOL_NREFL))
       allocate(x%parasolGrid_refl(Npoints,PARASOL_NREFL))
    endif

    ! Cloudsat simulator
    if (cosp_lradar_sim) then
       allocate(x%cloudsat_Ze_tot(Npoints,Ncolumns,Nlevels))
       allocate(x%cloudsat_cfad_ze(Npoints,cloudsat_DBZE_BINS,Nlvstat))
       allocate(x%cloudsat_precip_cover(Npoints,cloudsat_DBZE_BINS))
       allocate(x%cloudsat_pia(Npoints))
    endif

    ! Combined CALIPSO/CLOUDSAT fields
    if (cosp_lradar_sim .and. cosp_llidar_sim) then
       allocate(x%lidar_only_freq_cloud(Npoints,Nlvstat))
       allocate(x%radar_lidar_tcc(Npoints))
       allocate(x%cloudsat_tcc(Npoints))
       allocate(x%cloudsat_tcc2(Npoints))
    endif

    ! RTTOV
    if (cosp_lrttov_sim) then
      allocate(x%rttov_tbs(Npoints,Nchan))
    end if

    ! Joint MODIS/CloudSat Statistics
    if (cosp_lmodis_sim .and. cosp_lradar_sim) then
       allocate(x%wr_occfreq_ntotal(Npoints,WR_NREGIME))
       allocate(x%cfodd_ntotal(Npoints,CFODD_NDBZE,CFODD_NICOD,CFODD_NCLASS))
    end if

  end subroutine construct_cosp_outputs

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cospIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine destroy_cospIN(y)
    type(cosp_optical_inputs),intent(inout) :: y
    if (allocated(y%tau_067))             deallocate(y%tau_067)
    if (allocated(y%emiss_11))            deallocate(y%emiss_11)
    if (allocated(y%frac_out))            deallocate(y%frac_out)
    if (allocated(y%beta_mol_calipso))    deallocate(y%beta_mol_calipso)
    if (allocated(y%tau_mol_calipso))     deallocate(y%tau_mol_calipso)
    if (allocated(y%betatot_calipso))     deallocate(y%betatot_calipso)
    if (allocated(y%betatot_ice_calipso)) deallocate(y%betatot_ice_calipso)
    if (allocated(y%betatot_liq_calipso)) deallocate(y%betatot_liq_calipso)
    if (allocated(y%tautot_calipso))      deallocate(y%tautot_calipso)
    if (allocated(y%tautot_ice_calipso))  deallocate(y%tautot_ice_calipso)
    if (allocated(y%tautot_liq_calipso))  deallocate(y%tautot_liq_calipso)
    if (allocated(y%tautot_S_liq))        deallocate(y%tautot_S_liq)
    if (allocated(y%tautot_S_ice))        deallocate(y%tautot_S_ice)
    if (allocated(y%z_vol_cloudsat))      deallocate(y%z_vol_cloudsat)
    if (allocated(y%kr_vol_cloudsat))     deallocate(y%kr_vol_cloudsat)
    if (allocated(y%g_vol_cloudsat))      deallocate(y%g_vol_cloudsat)
    if (allocated(y%asym))                deallocate(y%asym)
    if (allocated(y%ss_alb))              deallocate(y%ss_alb)
    if (allocated(y%fracLiq))             deallocate(y%fracLiq)
    if (allocated(y%beta_mol_grLidar532)) deallocate(y%beta_mol_grLidar532)
    if (allocated(y%betatot_grLidar532))  deallocate(y%betatot_grLidar532)
    if (allocated(y%tau_mol_grLidar532))  deallocate(y%tau_mol_grLidar532)
    if (allocated(y%tautot_grLidar532))   deallocate(y%tautot_grLidar532)
    if (allocated(y%beta_mol_atlid))      deallocate(y%beta_mol_atlid)
    if (allocated(y%betatot_atlid))       deallocate(y%betatot_atlid)
    if (allocated(y%tau_mol_atlid))       deallocate(y%tau_mol_atlid)
    if (allocated(y%tautot_atlid))        deallocate(y%tautot_atlid)
    if (allocated(y%fracPrecipIce))      deallocate(y%fracPrecipIce)
  end subroutine destroy_cospIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cospstateIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine destroy_cospstateIN(y)
    type(cosp_column_inputs),intent(inout) :: y

    if (allocated(y%sunlit))          deallocate(y%sunlit)
    if (allocated(y%skt))             deallocate(y%skt)
    if (allocated(y%land))            deallocate(y%land)
    if (allocated(y%at))              deallocate(y%at)
    if (allocated(y%pfull))           deallocate(y%pfull)
    if (allocated(y%phalf))           deallocate(y%phalf)
    if (allocated(y%qv))              deallocate(y%qv)
    if (allocated(y%o3))              deallocate(y%o3)
    if (allocated(y%hgt_matrix))      deallocate(y%hgt_matrix)
    if (allocated(y%u_sfc))           deallocate(y%u_sfc)
    if (allocated(y%v_sfc))           deallocate(y%v_sfc)
    if (allocated(y%lat))             deallocate(y%lat)
    if (allocated(y%lon))             deallocate(y%lon)
    if (allocated(y%emis_sfc))        deallocate(y%emis_sfc)
    if (allocated(y%cloudIce))        deallocate(y%cloudIce)
    if (allocated(y%cloudLiq))        deallocate(y%cloudLiq)
    if (allocated(y%seaice))          deallocate(y%seaice)
    if (allocated(y%fl_rain))         deallocate(y%fl_rain)
    if (allocated(y%fl_snow))         deallocate(y%fl_snow)
    if (allocated(y%tca))             deallocate(y%tca)
    if (allocated(y%hgt_matrix_half)) deallocate(y%hgt_matrix_half)

  end subroutine destroy_cospstateIN

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cosp_outputs
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine destroy_cosp_outputs(y)
     type(cosp_outputs),intent(inout) :: y

     ! Deallocate and nullify
     if (associated(y%calipso_beta_mol)) then
        deallocate(y%calipso_beta_mol)
        nullify(y%calipso_beta_mol)
     endif
     if (associated(y%calipso_temp_tot)) then
        deallocate(y%calipso_temp_tot)
        nullify(y%calipso_temp_tot)
     endif
     if (associated(y%calipso_betaperp_tot)) then
        deallocate(y%calipso_betaperp_tot)
        nullify(y%calipso_betaperp_tot)
     endif
     if (associated(y%calipso_beta_tot)) then
        deallocate(y%calipso_beta_tot)
        nullify(y%calipso_beta_tot)
     endif
     if (associated(y%calipso_tau_tot)) then
        deallocate(y%calipso_tau_tot)
        nullify(y%calipso_tau_tot)
     endif
     if (associated(y%calipso_lidarcldphase)) then
        deallocate(y%calipso_lidarcldphase)
        nullify(y%calipso_lidarcldphase)
     endif
     if (associated(y%calipso_lidarcldtype)) then
        deallocate(y%calipso_lidarcldtype)
        nullify(y%calipso_lidarcldtype)
     endif
     if (associated(y%calipso_cldlayerphase)) then
        deallocate(y%calipso_cldlayerphase)
        nullify(y%calipso_cldlayerphase)
     endif
     if (associated(y%calipso_lidarcldtmp)) then
        deallocate(y%calipso_lidarcldtmp)
        nullify(y%calipso_lidarcldtmp)
     endif
     if (associated(y%calipso_cldlayer)) then
        deallocate(y%calipso_cldlayer)
        nullify(y%calipso_cldlayer)
     endif
     if (associated(y%calipso_cldtype)) then
        deallocate(y%calipso_cldtype)
        nullify(y%calipso_cldtype)
     endif
     if (associated(y%calipso_cldtypetemp)) then
        deallocate(y%calipso_cldtypetemp)
        nullify(y%calipso_cldtypetemp)
     endif
     if (associated(y%calipso_cldtypemeanz)) then
        deallocate(y%calipso_cldtypemeanz)
        nullify(y%calipso_cldtypemeanz)
     endif
     if (associated(y%calipso_cldtypemeanzse)) then
        deallocate(y%calipso_cldtypemeanzse)
        nullify(y%calipso_cldtypemeanzse)
     endif
     if (associated(y%calipso_cldthinemis)) then
        deallocate(y%calipso_cldthinemis)
        nullify(y%calipso_cldthinemis)
     endif
     if (associated(y%calipso_lidarcld)) then
        deallocate(y%calipso_lidarcld)
        nullify(y%calipso_lidarcld)
     endif
     if (associated(y%calipso_srbval)) then
        deallocate(y%calipso_srbval)
        nullify(y%calipso_srbval)
     endif
     if (associated(y%calipso_cfad_sr)) then
        deallocate(y%calipso_cfad_sr)
        nullify(y%calipso_cfad_sr)
     endif
     if (associated(y%grLidar532_beta_mol)) then
        deallocate(y%grLidar532_beta_mol)
        nullify(y%grLidar532_beta_mol)
     endif
     if (associated(y%grLidar532_beta_tot)) then
        deallocate(y%grLidar532_beta_tot)
        nullify(y%grLidar532_beta_tot)
     endif
     if (associated(y%grLidar532_cldlayer)) then
        deallocate(y%grLidar532_cldlayer)
        nullify(y%grLidar532_cldlayer)
     endif
     if (associated(y%grLidar532_lidarcld)) then
        deallocate(y%grLidar532_lidarcld)
        nullify(y%grLidar532_lidarcld)
     endif
     if (associated(y%grLidar532_cfad_sr)) then
        deallocate(y%grLidar532_cfad_sr)
        nullify(y%grLidar532_cfad_sr)
     endif
     if (associated(y%grLidar532_srbval)) then
        deallocate(y%grLidar532_srbval)
        nullify(y%grLidar532_srbval)
     endif
     if (associated(y%atlid_beta_mol)) then
        deallocate(y%atlid_beta_mol)
        nullify(y%atlid_beta_mol)
     endif
     if (associated(y%atlid_beta_tot)) then
        deallocate(y%atlid_beta_tot)
        nullify(y%atlid_beta_tot)
     endif
     if (associated(y%atlid_cldlayer)) then
        deallocate(y%atlid_cldlayer)
        nullify(y%atlid_cldlayer)
     endif
     if (associated(y%atlid_lidarcld)) then
        deallocate(y%atlid_lidarcld)
        nullify(y%atlid_lidarcld)
     endif
     if (associated(y%atlid_cfad_sr)) then
        deallocate(y%atlid_cfad_sr)
        nullify(y%atlid_cfad_sr)
     endif
     if (associated(y%atlid_srbval)) then
        deallocate(y%atlid_srbval)
        nullify(y%atlid_srbval)
     endif
     if (associated(y%parasolPix_refl)) then
        deallocate(y%parasolPix_refl)
        nullify(y%parasolPix_refl)
     endif
     if (associated(y%parasolGrid_refl)) then
        deallocate(y%parasolGrid_refl)
        nullify(y%parasolGrid_refl)
     endif
     if (associated(y%cloudsat_Ze_tot)) then
        deallocate(y%cloudsat_Ze_tot)
        nullify(y%cloudsat_Ze_tot)
     endif
     if (associated(y%cloudsat_cfad_ze)) then
        deallocate(y%cloudsat_cfad_ze)
        nullify(y%cloudsat_cfad_ze)
     endif
     if (associated(y%cloudsat_precip_cover)) then
        deallocate(y%cloudsat_precip_cover)
        nullify(y%cloudsat_precip_cover)
     endif
     if (associated(y%cloudsat_pia)) then
        deallocate(y%cloudsat_pia)
        nullify(y%cloudsat_pia)
     endif
     if (associated(y%cloudsat_tcc)) then
        deallocate(y%cloudsat_tcc)
        nullify(y%cloudsat_tcc)
     endif
     if (associated(y%cloudsat_tcc2)) then
        deallocate(y%cloudsat_tcc2)
        nullify(y%cloudsat_tcc2)
     endif
     if (associated(y%radar_lidar_tcc)) then
        deallocate(y%radar_lidar_tcc)
        nullify(y%radar_lidar_tcc)
     endif
     if (associated(y%cloudsat_tcc)) then
        deallocate(y%cloudsat_tcc)
        nullify(y%cloudsat_tcc)
     endif
     if (associated(y%cloudsat_tcc2)) then
        deallocate(y%cloudsat_tcc2)
        nullify(y%cloudsat_tcc2)
     endif
     if (associated(y%lidar_only_freq_cloud)) then
        deallocate(y%lidar_only_freq_cloud)
        nullify(y%lidar_only_freq_cloud)
     endif
     if (associated(y%isccp_totalcldarea)) then
        deallocate(y%isccp_totalcldarea)
        nullify(y%isccp_totalcldarea)
     endif
     if (associated(y%isccp_meantb)) then
        deallocate(y%isccp_meantb)
        nullify(y%isccp_meantb)
     endif
     if (associated(y%isccp_meantbclr)) then
        deallocate(y%isccp_meantbclr)
        nullify(y%isccp_meantbclr)
     endif
     if (associated(y%isccp_meanptop)) then
        deallocate(y%isccp_meanptop)
        nullify(y%isccp_meanptop)
     endif
     if (associated(y%isccp_meantaucld)) then
        deallocate(y%isccp_meantaucld)
        nullify(y%isccp_meantaucld)
     endif
     if (associated(y%isccp_meanalbedocld)) then
        deallocate(y%isccp_meanalbedocld)
        nullify(y%isccp_meanalbedocld)
     endif
     if (associated(y%isccp_boxtau)) then
        deallocate(y%isccp_boxtau)
        nullify(y%isccp_boxtau)
     endif
     if (associated(y%isccp_boxptop)) then
        deallocate(y%isccp_boxptop)
        nullify(y%isccp_boxptop)
     endif
     if (associated(y%isccp_fq)) then
        deallocate(y%isccp_fq)
        nullify(y%isccp_fq)
     endif
     if (associated(y%misr_fq)) then
        deallocate(y%misr_fq)
        nullify(y%misr_fq)
     endif
     if (associated(y%misr_dist_model_layertops)) then
        deallocate(y%misr_dist_model_layertops)
        nullify(y%misr_dist_model_layertops)
     endif
     if (associated(y%misr_meanztop)) then
        deallocate(y%misr_meanztop)
        nullify(y%misr_meanztop)
     endif
     if (associated(y%misr_cldarea)) then
        deallocate(y%misr_cldarea)
        nullify(y%misr_cldarea)
     endif
     if (associated(y%rttov_tbs)) then
        deallocate(y%rttov_tbs)
        nullify(y%rttov_tbs)
     endif
     if (associated(y%modis_Cloud_Fraction_Total_Mean)) then
        deallocate(y%modis_Cloud_Fraction_Total_Mean)
        nullify(y%modis_Cloud_Fraction_Total_Mean)
     endif
     if (associated(y%modis_Cloud_Fraction_Ice_Mean)) then
        deallocate(y%modis_Cloud_Fraction_Ice_Mean)
        nullify(y%modis_Cloud_Fraction_Ice_Mean)
     endif
     if (associated(y%modis_Cloud_Fraction_Water_Mean)) then
        deallocate(y%modis_Cloud_Fraction_Water_Mean)
        nullify(y%modis_Cloud_Fraction_Water_Mean)
     endif
     if (associated(y%modis_Cloud_Fraction_High_Mean)) then
        deallocate(y%modis_Cloud_Fraction_High_Mean)
        nullify(y%modis_Cloud_Fraction_High_Mean)
     endif
     if (associated(y%modis_Cloud_Fraction_Mid_Mean)) then
        deallocate(y%modis_Cloud_Fraction_Mid_Mean)
        nullify(y%modis_Cloud_Fraction_Mid_Mean)
     endif
     if (associated(y%modis_Cloud_Fraction_Low_Mean)) then
        deallocate(y%modis_Cloud_Fraction_Low_Mean)
        nullify(y%modis_Cloud_Fraction_Low_Mean)
     endif
     if (associated(y%modis_Optical_Thickness_Total_Mean)) then
        deallocate(y%modis_Optical_Thickness_Total_Mean)
        nullify(y%modis_Optical_Thickness_Total_Mean)
     endif
     if (associated(y%modis_Optical_Thickness_Water_Mean)) then
        deallocate(y%modis_Optical_Thickness_Water_Mean)
        nullify(y%modis_Optical_Thickness_Water_Mean)
     endif
     if (associated(y%modis_Optical_Thickness_Ice_Mean)) then
        deallocate(y%modis_Optical_Thickness_Ice_Mean)
        nullify(y%modis_Optical_Thickness_Ice_Mean)
     endif
     if (associated(y%modis_Optical_Thickness_Total_LogMean)) then
        deallocate(y%modis_Optical_Thickness_Total_LogMean)
        nullify(y%modis_Optical_Thickness_Total_LogMean)
     endif
     if (associated(y%modis_Optical_Thickness_Water_LogMean)) then
        deallocate(y%modis_Optical_Thickness_Water_LogMean)
        nullify(y%modis_Optical_Thickness_Water_LogMean)
     endif
     if (associated(y%modis_Optical_Thickness_Ice_LogMean)) then
        deallocate(y%modis_Optical_Thickness_Ice_LogMean)
        nullify(y%modis_Optical_Thickness_Ice_LogMean)
     endif
     if (associated(y%modis_Cloud_Particle_Size_Water_Mean)) then
        deallocate(y%modis_Cloud_Particle_Size_Water_Mean)
        nullify(y%modis_Cloud_Particle_Size_Water_Mean)
     endif
     if (associated(y%modis_Cloud_Particle_Size_Ice_Mean)) then
        deallocate(y%modis_Cloud_Particle_Size_Ice_Mean)
        nullify(y%modis_Cloud_Particle_Size_Ice_Mean)
     endif
     if (associated(y%modis_Cloud_Top_Pressure_Total_Mean)) then
        deallocate(y%modis_Cloud_Top_Pressure_Total_Mean)
        nullify(y%modis_Cloud_Top_Pressure_Total_Mean)
     endif
     if (associated(y%modis_Liquid_Water_Path_Mean)) then
        deallocate(y%modis_Liquid_Water_Path_Mean)
        nullify(y%modis_Liquid_Water_Path_Mean)
     endif
     if (associated(y%modis_Ice_Water_Path_Mean)) then
        deallocate(y%modis_Ice_Water_Path_Mean)
        nullify(y%modis_Ice_Water_Path_Mean)
     endif
     if (associated(y%modis_Optical_Thickness_vs_Cloud_Top_Pressure)) then
        deallocate(y%modis_Optical_Thickness_vs_Cloud_Top_Pressure)
        nullify(y%modis_Optical_Thickness_vs_Cloud_Top_Pressure)
     endif
     if (associated(y%modis_Optical_thickness_vs_ReffLIQ)) then
        deallocate(y%modis_Optical_thickness_vs_ReffLIQ)
        nullify(y%modis_Optical_thickness_vs_ReffLIQ)
     endif
     if (associated(y%modis_Optical_thickness_vs_ReffICE)) then
        deallocate(y%modis_Optical_thickness_vs_ReffICE)
        nullify(y%modis_Optical_thickness_vs_ReffICE)
     endif
     if (associated(y%cfodd_ntotal)) then
        deallocate(y%cfodd_ntotal)
        nullify(y%cfodd_ntotal)
     endif
     if (associated(y%wr_occfreq_ntotal)) then
        deallocate(y%wr_occfreq_ntotal)
        nullify(y%wr_occfreq_ntotal)
     endif

   end subroutine destroy_cosp_outputs

!#######################################################################
end module cospsimulator_intr
