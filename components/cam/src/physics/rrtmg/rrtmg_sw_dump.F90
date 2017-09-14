
#ifdef RRTMG_DUMP
#define _ERR(s,e,l) if (.not. s) then; write(*,*) l,', ',trim(e); stop; endif
#endif

module rrtmg_sw_dump
  implicit none
  !I'm making this a module variable so that it remains persistent across subsequent calls to crm_dump_input() and crm_dump_output(),
  !meaning that if a given time step's input is being use, so should its related output.
  !I'll default it to true so that if CRM_DUMP_RATIO isn't defined, it will output every time step
  logical :: do_dump = .true.
  logical :: first_time = .true.

contains


  subroutine rrtmg_sw_dump_input(igstep,lchnk, Nday, rrtmg_levs, icld,         &
                 pmidmb, pintmb, tlay, tlev, tsfc, &
                 h2ovmr, o3vmr, co2vmr, ch4vmr, o2vmr, n2ovmr, &
                 asdir, asdif, aldir, aldif)
!                 coszrs, eccf, dyofyr, solvar, &
!                 inflgsw, iceflgsw, liqflgsw, &
!                 cld_stosw, tauc_stosw, ssac_stosw, asmc_stosw, fsfc_stosw, &
!                 cicewp_stosw, cliqwp_stosw, rei, rel, &
!                 tau_aer_sw, ssa_aer_sw, asm_aer_sw)
    use mpi
    use dmdf
    use shr_kind_mod, only: r8 => shr_kind_r8
    use ppgrid,          only: pcols
    use parrrsw,         only: nbndsw, ngptsw
    use crmdims
    use microphysics, only: nmicro_fields
    use params,       only: crm_rknd
    
    implicit none
       ! Input arguments
    integer, intent(in) :: igstep
    integer, intent(in) :: lchnk             ! chunk identifier
    integer, intent(in) :: Nday              ! Number of daylight columns
    integer, intent(in) :: rrtmg_levs        ! number of levels rad is applied
    integer, intent(in) :: icld                          ! Flag for cloud overlap method
    real(r8), intent(in) :: pmidmb(pcols,rrtmg_levs)   ! Level pressure (hPa)
    real(r8), intent(in) :: pintmb(pcols,rrtmg_levs+1) ! Model interface pressure (hPa)
    real(r8), intent(in) :: tlay(pcols,rrtmg_levs)     ! mid point temperature
    real(r8), intent(in) :: tlev(pcols,rrtmg_levs+1)   ! interface temperature
    real(r8), intent(in) :: tsfc(pcols)          ! surface temperature
    real(r8), intent(in) :: h2ovmr(pcols,rrtmg_levs)   ! h2o volume mixing ratio
    real(r8), intent(in) :: o3vmr(pcols,rrtmg_levs)    ! o3 volume mixing ratio
    real(r8), intent(in) :: co2vmr(pcols,rrtmg_levs)   ! co2 volume mixing ratio 
    real(r8), intent(in) :: ch4vmr(pcols,rrtmg_levs)   ! ch4 volume mixing ratio 
    real(r8), intent(in) :: o2vmr(pcols,rrtmg_levs)    ! o2  volume mixing ratio 
    real(r8), intent(in) :: n2ovmr(pcols,rrtmg_levs)   ! n2o volume mixing ratio 
    real(r8), intent(in) :: asdir(pcols)     ! 0.2-0.7 micro-meter srfc alb: direct rad
    real(r8), intent(in) :: aldir(pcols)     ! 0.7-5.0 micro-meter srfc alb: direct rad
    real(r8), intent(in) :: asdif(pcols)     ! 0.2-0.7 micro-meter srfc alb: diffuse rad
    real(r8), intent(in) :: aldif(pcols)     ! 0.7-5.0 micro-meter srfc alb: diffuse rad
!    real(r8) :: coszrs(pcols)    ! Cosine solar zenith angle
!    real(r8), intent(in) :: eccf               ! Eccentricity factor (1./earth-sun dist^2)
!    integer :: dyofyr                ! Set to day of year for Earth/Sun distance calculation in
                                    ! rrtmg_sw, or pass in adjustment directly into adjes
!    real(r8) :: solvar              ! solar irradiance variability in each band
!    integer  :: inflgsw               ! flag for cloud parameterization method
!    integer  :: iceflgsw              ! flag for ice cloud parameterization method
!    integer  :: liqflgsw              ! flag for liquid cloud parameterization method
!    integer, parameter :: nsubcsw = ngptsw           ! rrtmg_sw g-point (quadrature point) dimension
!    real(r8) :: cld_stosw(nsubcsw, pcols, rrtmg_levs-1)      ! stochastic cloud fraction
!    real(r8) :: tauc_stosw(nbndsw, pcols, rrtmg_levs-1)         ! cloud optical depth
!    real(r8) :: ssac_stosw(nbndsw, pcols, rrtmg_levs-1)         ! cloud single scat. albedo
!    real(r8) :: asmc_stosw(nbndsw, pcols, rrtmg_levs-1)         ! cloud asymmetry parameter
!    real(r8) :: fsfc_stosw(nbndsw, pcols, rrtmg_levs-1)         ! cloud forward scattering fraction
!    real(r8) :: cicewp_stosw(nsubcsw, pcols, rrtmg_levs-1)   ! stochastic cloud ice water path
!    real(r8) :: cliqwp_stosw(nsubcsw, pcols, rrtmg_levs-1)   ! stochastic cloud liquid wter path
!    real(r8) :: rei(pcols,rrtmg_levs-1)    ! Ice effective drop size (microns)
!    real(r8) :: rel(pcols,rrtmg_levs-1)    ! Liquid effective drop size (microns)
!    real(r8) :: tau_aer_sw(pcols, rrtmg_levs-1, nbndsw)      ! aer optical depth
!    real(r8) :: ssa_aer_sw(pcols, rrtmg_levs-1, nbndsw)      ! aer single scat. albedo
!    real(r8) :: asm_aer_sw(pcols, rrtmg_levs-1, nbndsw)      ! aer asymmetry parameter
    
    integer :: myrank, ierr
    real(crm_rknd) :: myrand    
#ifdef RRTMG_DUMP
        call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)
#ifdef RRTMG_DUMP_RATIO
        !If CRM_DUMP_RATIO is defined, then you only want to output that proportion of the time
        !So sample a random number, and if it's less than the given ratio, do the output; otherwise don't
        !Also note that crm_dump is controlled only here. the output() routine shared the same value that is
        !set by this input routine
        call random_number(myrand)
        do_dump = .false.
        if (myrand < RRTMG_DUMP_RATIO) do_dump = .true.
#endif
        if (do_dump) then
          !The first time we create the file, we have to make sure that we have the needed attributes for ./setup_from_file.sh in the standalone model to work properly.
          if (first_time) then
            first_time = .false.
            call dmdf_write_attr(PLEV     ,myrank,'rrtmg_sw_in','plev'    ); _ERR(success,error_string,__LINE__)
            call dmdf_write_attr(PSUBCOLS ,myrank,'rrtmg_sw_in','psubcols'); _ERR(success,error_string,__LINE__)
            call dmdf_write_attr(PCOLS    ,myrank,'rrtmg_sw_in','pcols'   ); _ERR(success,error_string,__LINE__)
            call dmdf_write_attr(PCNST    ,myrank,'rrtmg_sw_in','pcnst'   ); _ERR(success,error_string,__LINE__)
          endif
          !Don't output the first time because fields are initially uniform, and random noise is added. Keep stochasticity out of the standalone model
          if (igstep > 1) then

          WRITE(*,*),"WRITE FILES IN STARTED", myrank

            call dmdf_write(lchnk            ,myrank,'rrtmg_sw_in','lchnk'                                                             ,.true. ,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(Nday             ,myrank,'rrtmg_sw_in','Nday'                                                              ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(rrtmg_levs       ,myrank,'rrtmg_sw_in','rrtmg_levs'                                                        ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(icld             ,myrank,'rrtmg_sw_in','icld'                                                              ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(pmidmb           ,myrank,'rrtmg_sw_in','pmidmb'          ,(/'pcols','rrtmg_levs'/)                         ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(pintmb           ,myrank,'rrtmg_sw_in','pmitmb'          ,(/'pcols','rrtmg_levs_p1'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(tlay             ,myrank,'rrtmg_sw_in','tlay'            ,(/'pcols','rrtmg_levs'/)                         ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(tlev             ,myrank,'rrtmg_sw_in','tlev'            ,(/'pcols','rrtmg_levs_p1'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(tsfc             ,myrank,'rrtmg_sw_in','tsfc    '        ,(/'pcols'/)                                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(h2ovmr           ,myrank,'rrtmg_sw_in','h2ovmr'          ,(/'pcols','rrtmg_levs'/)                         ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(o3vmr            ,myrank,'rrtmg_sw_in','o3vmr'           ,(/'pcols','rrtmg_levs'/)                         ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(co2vmr           ,myrank,'rrtmg_sw_in','co2vmr'          ,(/'pcols','rrtmg_levs'/)                         ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(ch4vmr           ,myrank,'rrtmg_sw_in','ch4vmr'          ,(/'pcols','rrtmg_levs'/)                         ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(o2vmr            ,myrank,'rrtmg_sw_in','o2vmr'           ,(/'pcols','rrtmg_levs'/)                         ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(n2ovmr           ,myrank,'rrtmg_sw_in','n2ovmr'          ,(/'pcols','rrtmg_levs'/)                         ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(asdir            ,myrank,'rrtmg_sw_in','asdir'          ,(/'pcols'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(aldir            ,myrank,'rrtmg_sw_in','aldir'          ,(/'pcols'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(asdif            ,myrank,'rrtmg_sw_in','asdif'          ,(/'pcols'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(aldif            ,myrank,'rrtmg_sw_in','aldif'          ,(/'pcols'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)        


            WRITE(*,*),"WRITE FILE ENDS",myrank
!            call dmdf_write(coszrs           ,myrank,'rrtmg_sw_in','coszrs'         ,(/'pcols'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
!            call dmdf_write(eccf             ,myrank,'rrtmg_sw_in','eccf'                                                              ,.false.,.false.); _ERR(success,error_string,__LINE__)
!            call dmdf_write(dyofyr           ,myrank,'rrtmg_sw_in','dyofyr'                                                            ,.false.,.false.); _ERR(success,error_string,__LINE__)
!            call dmdf_write(solvar           ,myrank,'rrtmg_sw_in','solvar'                                                            ,.false.,.false.); _ERR(success,error_string,__LINE__)
!            call dmdf_write(inflgsw          ,myrank,'rrtmg_sw_in','inflgsw'                                                           ,.false.,.false.); _ERR(success,error_string,__LINE__)
!            call dmdf_write(iceflgsw         ,myrank,'rrtmg_sw_in','iceflgsw'                                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
!            call dmdf_write(liqflgsw         ,myrank,'rrtmg_sw_in','liqflgsw'                                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
!            call dmdf_write(nsubcsw          ,myrank,'rrtmg_sw_in','nsubcsw'                                                           ,.false.,.false.); _ERR(success,error_string,__LINE__)
!            call dmdf_write(cld_stosw        ,myrank,'rrtmg_sw_in','cld_stosw'      ,(/'nsubcsw', 'pcols', 'rrtmg_levs_m1'/)           ,.false.,.false.); _ERR(success,error_string,__LINE__)
!            call dmdf_write(tauc_stosw          ,myrank,'rrtmg_sw_in','tauc_sw'        ,(/'nsubcsw', 'pcols', 'rrtmg_levs_m1'/)           ,.false.,.false.); _ERR(success,error_string,__LINE__)
!            call dmdf_write(ssac_stosw          ,myrank,'rrtmg_sw_in','ssac_sw'        ,(/'nsubcsw', 'pcols', 'rrtmg_levs_m1'/)           ,.false.,.false.); _ERR(success,error_string,__LINE__)
!            call dmdf_write(asmc_stosw          ,myrank,'rrtmg_sw_in','asmc_sw'        ,(/'nsubcsw', 'pcols', 'rrtmg_levs_m1'/)           ,.false.,.false.); _ERR(success,error_string,__LINE__)
!            call dmdf_write(fsfc_stosw          ,myrank,'rrtmg_sw_in','fsfc_sw'        ,(/'nsubcsw', 'pcols', 'rrtmg_levs_m1'/)           ,.false.,.false.); _ERR(success,error_string,__LINE__)
!            call dmdf_write(cicewp_stosw     ,myrank,'rrtmg_sw_in','cicewp_stosw'   ,(/'nsubcsw', 'pcols', 'rrtmg_levs_m1'/)           ,.false.,.false.); _ERR(success,error_string,__LINE__)
!            call dmdf_write(cliqwp_stosw     ,myrank,'rrtmg_sw_in','cliqwp_stosw'   ,(/'nsubcsw', 'pcols', 'rrtmg_levs_m1'/)           ,.false.,.false.); _ERR(success,error_string,__LINE__)
!            call dmdf_write(rei              ,myrank,'rrtmg_sw_in','rei'            ,(/'pcols', 'rrtmg_levs_m1'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
!            call dmdf_write(rel              ,myrank,'rrtmg_sw_in','rel'            ,(/'pcols', 'rrtmg_levs_m1'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
!            call dmdf_write(tau_aer_sw       ,myrank,'rrtmg_sw_in','tau_aer_sw'     ,(/'pcols', 'rrtmg_levs_m1', 'nbndsw'/)            ,.false.,.false.); _ERR(success,error_string,__LINE__)
!            call dmdf_write(ssa_aer_sw       ,myrank,'rrtmg_sw_in','ssa_aer_sw'     ,(/'pcols', 'rrtmg_levs_m1', 'nbndsw'/)            ,.false.,.false.); _ERR(success,error_string,__LINE__)
!            call dmdf_write(asm_aer_sw       ,myrank,'rrtmg_sw_in','asm_aer_sw'     ,(/'pcols', 'rrtmg_levs_m1', 'nbndsw'/)            ,.false.,.true.); _ERR(success,error_string,__LINE__)


        endif
        endif
#endif

  end subroutine rrtmg_sw_dump_input



  subroutine rrtmg_sw_dump_output(igstep,lchnk, Nday, rrtmg_levs, icld,         &
                 swuflx, swdflx, swhr, swuflxc, swdflxc, swhrc, &
                 dirdnuv, dirdnir, difdnuv, difdnir, ninflx, ninflxc)
    use mpi
    use dmdf
    use shr_kind_mod, only: r8 => shr_kind_r8
    use ppgrid,          only: pcols
    use crmdims
    use microphysics, only: nmicro_fields
#ifdef ECPP
    use ecppvars,  only: NCLASS_CL, ncls_ecpp_in, NCLASS_PR
#endif
    implicit none
       ! Input arguments
    integer, intent(in) :: igstep
    integer, intent(in) :: lchnk             ! chunk identifier
    integer, intent(in) :: Nday              ! Number of daylight columns
    integer, intent(in) :: rrtmg_levs        ! number of levels rad is applied
    integer :: icld                          ! Flag for cloud overlap method
    real(r8) :: swuflx(pcols,rrtmg_levs+1)       ! Total sky shortwave upward flux (W/m2)
    real(r8) :: swdflx(pcols,rrtmg_levs+1)       ! Total sky shortwave downward flux (W/m2)
    real(r8) :: swhr(pcols,rrtmg_levs)           ! Total sky shortwave radiative heating rate (K/d)
    real(r8) :: swuflxc(pcols,rrtmg_levs+1)      ! Clear sky shortwave upward flux (W/m2)
    real(r8) :: swdflxc(pcols,rrtmg_levs+1)      ! Clear sky shortwave downward flux (W/m2)
    real(r8) :: swhrc(pcols,rrtmg_levs)          ! Clear sky shortwave radiative heating rate (K/d)
    real(r8) :: dirdnuv(pcols,rrtmg_levs+1)       ! Direct downward shortwave flux, UV/vis
    real(r8) :: difdnuv(pcols,rrtmg_levs+1)       ! Diffuse downward shortwave flux, UV/vis
    real(r8) :: dirdnir(pcols,rrtmg_levs+1)       ! Direct downward shortwave flux, near-IR
    real(r8) :: difdnir(pcols,rrtmg_levs+1)       ! Diffuse downward shortwave flux, near-IR

   ! Added for net near-IR diagnostic
    real(r8) :: ninflx(pcols,rrtmg_levs+1)        ! Net shortwave flux, near-IR
    real(r8) :: ninflxc(pcols,rrtmg_levs+1)       ! Net clear sky shortwave flux, near-IR

    integer :: myrank, ierr

#ifdef RRTMG_DUMP
    call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)
    if (do_dump) then
      !Don't output the first time because fields are initially uniform, and random noise is added. Keep stochasticity out of the standalone model
      if (igstep > 1) then
            call dmdf_write(lchnk           ,myrank,'rrtmg_sw_out','lchnk'        ,.true. ,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(Nday            ,myrank,'rrtmg_sw_out','Nday'         ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(rrtmg_levs      ,myrank,'rrtmg_sw_out','rrtmg_levs'   ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(icld            ,myrank,'rrtmg_sw_out','icld'         ,.false.,.false.); _ERR(success,error_string,__LINE__)

           call dmdf_write(swuflx           ,myrank,'rrtmg_sw_out',trim('swuflx          '),(/'pcols','rrtmg_levs_p1'/)               ,.false. ,.false.); _ERR(success,error_string,__LINE__)
           call dmdf_write(swdflx           ,myrank,'rrtmg_sw_out',trim('swdflx          '),(/'pcols','rrtmg_levs_p1'/)               ,.false. ,.false.); _ERR(success,error_string,__LINE__)
           call dmdf_write(swhr             ,myrank,'rrtmg_sw_out',trim('swhr            '),(/'pcols','rrtmg_levs'/)                  ,.false. ,.false.); _ERR(success,error_string,__LINE__)
           call dmdf_write(swuflxc          ,myrank,'rrtmg_sw_out',trim('swuflxc         '),(/'pcols','rrtmg_levs_p1'/)               ,.false. ,.false.); _ERR(success,error_string,__LINE__)     
           call dmdf_write(swdflxc          ,myrank,'rrtmg_sw_out',trim('swdflxc         '),(/'pcols','rrtmg_levs_p1'/)               ,.false. ,.false.); _ERR(success,error_string,__LINE__) 
           call dmdf_write(swhrc            ,myrank,'rrtmg_sw_out',trim('swhrc           '),(/'pcols','rrtmg_levs'/)                  ,.false. ,.false.); _ERR(success,error_string,__LINE__)
           call dmdf_write(dirdnuv          ,myrank,'rrtmg_sw_out',trim('dirdnuv         '),(/'pcols','rrtmg_levs_p1'/)               ,.false. ,.false.); _ERR(success,error_string,__LINE__)
           call dmdf_write(difdnuv          ,myrank,'rrtmg_sw_out',trim('difdnuv         '),(/'pcols','rrtmg_levs_p1'/)               ,.false. ,.false.); _ERR(success,error_string,__LINE__)
           call dmdf_write(dirdnir          ,myrank,'rrtmg_sw_out',trim('dirdnir         '),(/'pcols','rrtmg_levs_p1'/)               ,.false. ,.false.); _ERR(success,error_string,__LINE__)
           call dmdf_write(difdnir          ,myrank,'rrtmg_sw_out',trim('difdnir         '),(/'pcols','rrtmg_levs_p1'/)               ,.false. ,.false.); _ERR(success,error_string,__LINE__)
           call dmdf_write(ninflx           ,myrank,'rrtmg_sw_out',trim('ninflx          '),(/'pcols','rrtmg_levs_p1'/)               ,.false. ,.false.); _ERR(success,error_string,__LINE__)
           call dmdf_write(ninflxc          ,myrank,'rrtmg_sw_out',trim('ninflxc         '),(/'pcols','rrtmg_levs_p1'/)               ,.false. ,.true.); _ERR(success,error_string,__LINE__)

     endif
    endif
#endif
  end subroutine rrtmg_sw_dump_output


end module rrtmg_sw_dump
