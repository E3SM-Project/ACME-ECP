
#ifdef RRTMG_DUMP
#define _ERR(s,e,l) if (.not. s) then; write(*,*) l,', ',trim(e); stop; endif
#endif

module rrtmg_lw_dump

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pcols, begchunk, endchunk

!      use parkind, only : jpim, jprb 
  use rrlw_vsn
  use mcica_subcol_gen_lw, only: mcica_subcol_lw
  use rrtmg_lw_cldprmc, only: cldprmc
! Move call to rrtmg_lw_ini and following use association to 
! GCM initialization area
!      use rrtmg_lw_init, only: rrtmg_lw_ini
  use rrtmg_lw_rtrnmc, only: rtrnmc
  use rrtmg_lw_setcoef, only: setcoef
  use rrtmg_lw_taumol, only: taumol

  implicit none
  !I'm making this a module variable so that it remains persistent across subsequent calls to crm_dump_input() and crm_dump_output(),
  !meaning that if a given time step's input is being use, so should its related output.
  !I'll default it to true so that if CRM_DUMP_RATIO isn't defined, it will output every time step
  logical :: do_dump = .true.
  logical :: first_time = .true.

contains


    subroutine rrtmg_lw_dump_input(igstep, &
             lchnk   ,ncol    ,nlay    ,icld    , &
             play    ,plev    ,tlay    ,tlev    ,tsfc    ,h2ovmr  , &
             o3vmr   ,co2vmr  ,ch4vmr  ,o2vmr, n2ovmr  ,cfc11vmr,cfc12vmr, &
             cfc22vmr,ccl4vmr ,emis    ,inflglw ,iceflglw,liqflglw, &
             cldfr   ,taucld  ,cicewp  ,cliqwp  ,reice   ,reliq   , &
             tauaer)

    use mpi
    use dmdf
    use shr_kind_mod, only: r8 => shr_kind_r8
    use parrrtm,         only: nbndlw, ngptlw 
    use crmdims
    use params,       only: crm_rknd
    use microphysics, only: nmicro_fields
    
    implicit none
       ! Input arguments
      integer, intent(in) :: igstep
      integer, intent(in) :: lchnk                      ! chunk identifier
      integer, intent(in) :: ncol                       ! Number of horizontal columns
      integer, intent(in) :: nlay                       ! Number of model layers
      integer, intent(in) :: icld                    ! Cloud overlap method
                                                        !    0: Clear only
                                                        !    1: Random
                                                        !    2: Maximum/random
                                                        !    3: Maximum

      real(kind=r8), intent(in) :: play(ncol,nlay)      ! Layer pressures (hPa, mb)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: plev(ncol,nlay+1)    ! Interface pressures (hPa, mb)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(in) :: tlay(ncol,nlay)      ! Layer temperatures (K)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: tlev(ncol,nlay+1)    ! Interface temperatures (K)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(in) :: tsfc(ncol)           ! Surface temperature (K)
                                                        !    Dimensions: (ncol)
      real(kind=r8), intent(in) :: h2ovmr(ncol,nlay)    ! H2O volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: o3vmr(ncol,nlay)     ! O3 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: co2vmr(ncol,nlay)    ! CO2 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: ch4vmr(ncol,nlay)    ! Methane volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: o2vmr(ncol,nlay)     ! O2 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: n2ovmr(ncol,nlay)    ! Nitrous oxide volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: cfc11vmr(ncol,nlay)  ! CFC11 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: cfc12vmr(ncol,nlay)  ! CFC12 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: cfc22vmr(ncol,nlay)  ! CFC22 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: ccl4vmr(ncol,nlay)   ! CCL4 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: emis(ncol,nbndlw)    ! Surface emissivity
                                                        !    Dimensions: (ncol,nbndlw)

      integer, intent(in) :: inflglw                    ! Flag for cloud optical properties
      integer, intent(in) :: iceflglw                   ! Flag for ice particle specification
      integer, intent(in) :: liqflglw                   ! Flag for liquid droplet specification
 
      real(kind=r8), intent(in) :: cldfr(ncol,nlay)     ! Cloud fraction
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: cicewp(ncol,nlay)    ! Cloud ice water path (g/m2)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: cliqwp(ncol,nlay)    ! Cloud liquid water path (g/m2)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: reice(ncol,nlay)     ! Cloud ice effective radius (microns)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: reliq(ncol,nlay)     ! Cloud water drop effective radius (microns)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: taucld(nbndlw,ncol,nlay)    ! Cloud optical depth
                                                        !    Dimensions: (nbndlw,ncol,nlay)
      real(kind=r8), intent(in) :: tauaer(ncol,nlay,nbndlw)        ! aerosol optical depth
                                                        !   at mid-point of LW spectral bands
                                                        !    Dimensions: (ncol,nlay,nbndlw)
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
        !Don't output the first time because fields are initially uniform, and random noise is added. Keep stochasticity out of the standalone model
          if (igstep > 1) then
            call dmdf_write(lchnk           ,myrank,'rrtmg_lw_in','lchnk'        ,.true. ,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(ncol            ,myrank,'rrtmg_lw_in','ncol'         ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(nlay            ,myrank,'rrtmg_lw_in','nlay'         ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(icld            ,myrank,'rrtmg_lw_in','icld'         ,.false.,.false.); _ERR(success,error_string,__LINE__)

            call dmdf_write(play          ,myrank,'rrtmg_lw_in','play'          ,(/'ncol','nlay'/)                                 ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(plev          ,myrank,'rrtmg_lw_in','plev'          ,(/'ncol','nlay_p1'/)                              ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(tlay          ,myrank,'rrtmg_lw_in','tlay'          ,(/'ncol','nlay'/)                                ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(tlev          ,myrank,'rrtmg_lw_in','tlev'          ,(/'ncol','nlay_p1'/)                             ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(tsfc          ,myrank,'rrtmg_lw_in','tsfc'          ,(/'ncol'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(h2ovmr        ,myrank,'rrtmg_lw_in','h2ovmr'        ,(/'ncol','nlay'/)                                ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(o3vmr         ,myrank,'rrtmg_lw_in','o3vmr'         ,(/'ncol','nlay'/)                                ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(co2vmr        ,myrank,'rrtmg_lw_in','co2vmr'        ,(/'ncol','nlay'/)                                ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(ch4vmr        ,myrank,'rrtmg_lw_in','ch4vmr'        ,(/'ncol','nlay'/)                                ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(o2vmr         ,myrank,'rrtmg_lw_in','o2vmr'         ,(/'ncol','nlay'/)                                ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(n2ovmr        ,myrank,'rrtmg_lw_in','n2ovmr'        ,(/'ncol','nlay'/)                                ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(cfc11vmr      ,myrank,'rrtmg_lw_in','cfc11vmr'      ,(/'ncol','nlay'/)                                ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(cfc12vmr      ,myrank,'rrtmg_lw_in','cfc12vmr'      ,(/'ncol','nlay'/)                                ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(cfc22vmr      ,myrank,'rrtmg_lw_in','cfc22vmr'      ,(/'ncol','nlay'/)                                ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(ccl4vmr       ,myrank,'rrtmg_lw_in','ccl4vmr'       ,(/'ncol','nlay'/)                                ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(emis          ,myrank,'rrtmg_lw_in','emis'          ,(/'ncol','nbndlw'/)                             ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(inflglw       ,myrank,'rrtmg_lw_in','inflglw'                                                           ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(iceflglw      ,myrank,'rrtmg_lw_in','iceflglw'                                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(liqflglw      ,myrank,'rrtmg_lw_in','liqflglw'                                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
            
            call dmdf_write(cldfr         ,myrank,'rrtmg_lw_in','cldfr'         ,(/'ncol','nlay'/)                                 ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(cicewp         ,myrank,'rrtmg_lw_in','cicewp'       ,(/'ncol','nlay'/)                                 ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(cliqwp         ,myrank,'rrtmg_lw_in','cliqwp'       ,(/'ncol','nlay'/)                                 ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(reice          ,myrank,'rrtmg_lw_in','reice'        ,(/'ncol','nlay'/)                                 ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(reliq          ,myrank,'rrtmg_lw_in','reliq'        ,(/'ncol','nlay'/)                                 ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(taucld         ,myrank,'rrtmg_lw_in','taucld'       ,(/'nbndlw','ncol','nlay'/)                        ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(tauaer         ,myrank,'rrtmg_lw_in','tauaer'       ,(/'ncol','nlay','nbndlw'/)                        ,.false.,.true.); _ERR(success,error_string,__LINE__)


        endif
        endif
#endif

  end subroutine rrtmg_lw_dump_input



  subroutine rrtmg_lw_dump_output(igstep, &
                                 lchnk   ,ncol    ,nlay    ,icld,         &
                                 uflx    ,dflx    ,hr      ,uflxc   ,dflxc   ,hrc, &
                                 uflxs   ,dflxs )
    use mpi
    use dmdf
    use shr_kind_mod, only: r8 => shr_kind_r8
    use ppgrid,          only: pcols
    use crmdims
    use microphysics, only: nmicro_fields
    use parrrtm,         only: nbndlw, ngptlw   

 
    implicit none
       ! Input arguments
    integer, intent(in) :: igstep
    integer, intent(in) :: lchnk                      ! chunk identifier
    integer, intent(in) :: ncol                       ! Number of horizontal columns
    integer, intent(in) :: nlay                       ! Number of model layers
    integer, intent(in) :: icld                       ! Cloud overlap method
                                                        !    0: Clear only
                                                        !    1: Random
                                                        !    2: Maximum/random
                                                        !    3: Maximum

    real(kind=r8), intent(in) :: uflx(ncol,nlay+1)           ! Total sky longwave upward flux (W/m2)
                                                        !    Dimensions: (ncol,nlay+1)
    real(kind=r8), intent(in) :: dflx(ncol,nlay+1)           ! Total sky longwave downward flux (W/m2)
                                                        !    Dimensions: (ncol,nlay+1)
    real(kind=r8), intent(in) :: hr(ncol,nlay)          ! Total sky longwave radiative heating rate (K/d)
                                                        !    Dimensions: (ncol,nlay)
    real(kind=r8), intent(in) :: uflxc(ncol,nlay+1)          ! Clear sky longwave upward flux (W/m2)
                                                        !    Dimensions: (ncol,nlay+1)
    real(kind=r8), intent(in) :: dflxc(ncol,nlay+1)          ! Clear sky longwave downward flux (W/m2)
                                                       !    Dimensions: (ncol,nlay+1)
    real(kind=r8), intent(in) :: hrc(ncol,nlay)            ! Clear sky longwave radiative heating rate (K/d)
                                                        !    Dimensions: (ncol,nlay)
    real(kind=r8), intent(in) :: uflxs(nbndlw,ncol,nlay+1)        ! Total sky longwave upward flux spectral (W/m2)
                                                        !    Dimensions: (nbndlw,ncol,nlay+1)
    real(kind=r8), intent(in) :: dflxs(nbndlw,ncol,nlay+1)        ! Total sky longwave downward flux spectral (W/m2)
                                                        !    Dimensions: (nbndlw,ncol,nlay+1)


    integer :: myrank, ierr

#ifdef RRTMG_DUMP
    call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)
    if (do_dump) then
      !Don't output the first time because fields are initially uniform, and random noise is added. Keep stochasticity out of the standalone model
      if (igstep > 1) then

           call dmdf_write(lchnk           ,myrank,'rrtmg_lw_in','lchnk'        ,.true. ,.false.); _ERR(success,error_string,__LINE__)
           call dmdf_write(ncol            ,myrank,'rrtmg_lw_in','ncol'         ,.false.,.false.); _ERR(success,error_string,__LINE__)
           call dmdf_write(nlay            ,myrank,'rrtmg_lw_in','nlay'         ,.false.,.false.); _ERR(success,error_string,__LINE__)
           call dmdf_write(icld            ,myrank,'rrtmg_lw_in','icld'         ,.false.,.false.); _ERR(success,error_string,__LINE__)


           call dmdf_write(uflx            ,myrank,'rrtmg_lw_out',trim('uflx          '),(/'ncol','nlay_p1'/)                 ,.false. ,.false.); _ERR(success,error_string,__LINE__)
           call dmdf_write(uflxc           ,myrank,'rrtmg_lw_out',trim('uflxc         '),(/'ncol','nlay_p1'/)                 ,.false. ,.false.); _ERR(success,error_string,__LINE__)
           call dmdf_write(hr              ,myrank,'rrtmg_lw_out',trim('hr            '),(/'ncol','nlay'/)                    ,.false. ,.false.); _ERR(success,error_string,__LINE__)
           call dmdf_write(uflxc           ,myrank,'rrtmg_lw_out',trim('uflxc         '),(/'ncol','nlay_p1'/)                 ,.false. ,.false.); _ERR(success,error_string,__LINE__)     
           call dmdf_write(dflxc           ,myrank,'rrtmg_lw_out',trim('dflxc         '),(/'ncol','nlay_p1'/)                 ,.false. ,.false.); _ERR(success,error_string,__LINE__) 
           call dmdf_write(hrc             ,myrank,'rrtmg_lw_out',trim('hrc           '),(/'ncol','nlay'/)                    ,.false. ,.true.); _ERR(success,error_string,__LINE__)

           call dmdf_write(uflxs           ,myrank,'rrtmg_lw_out',trim('uflxs         '),(/'nbndlw','ncol','nlay_p1'/)        ,.false. ,.false.); _ERR(success,error_string,__LINE__)
           call dmdf_write(dflxs           ,myrank,'rrtmg_lw_out',trim('dflxs         '),(/'nbndlw','ncol','nlay_p1'/)        ,.false. ,.true.); _ERR(success,error_string,__LINE__)



     endif
    endif
#endif
  end subroutine rrtmg_lw_dump_output


end module rrtmg_lw_dump
