
#ifdef CRM_DUMP
#define _ERR(s,e,l) if (.not. s) then; write(*,*) l,', ',trim(e); stop; endif
#endif

module crm_dump
  implicit none
  !I'm making this a module variable so that it remains persistent across subsequent calls to crm_dump_input() and crm_dump_output(),
  !meaning that if a given time step's input is being use, so should its related output.
  !I'll default it to true so that if CRM_DUMP_RATIO isn't defined, it will output every time step
  logical :: do_dump = .true.
  logical :: first_time = .true.

contains


  subroutine crm_dump_input(igstep,plev,lchnk,icol,latitude0,longitude0,ps,pmid,pdel,phis,zmid,zint,qrad_crm,dt_gl,ocnfrac,tau00,&
                            wndls,bflxls,fluxu00,fluxv00,fluxt00,fluxq00,tl,ql,qccl,qiil,ul,vl, &
#ifdef CLUBB_CRM
                            clubb_buffer    , &
#endif
                            cltot,clhgh,clmed,cllow,u_crm,v_crm,w_crm,t_crm,micro_fields_crm, &
#ifdef m2005
#ifdef MODAL_AERO
                            naermod,vaerosol,hygro , &
#endif
#endif
                            dd_crm,mui_crm,mdi_crm )
    use mpi
    use dmdf
    use shr_kind_mod, only: r8 => shr_kind_r8
    use crmdims
    use microphysics, only: nmicro_fields
    use params,       only: crm_rknd
#ifdef MODAL_AERO
        use modal_aero_data,   only: ntot_amode
#endif
    implicit none
    integer , intent(in) :: igstep
    integer , intent(in) :: lchnk    ! chunk identifier
    integer , intent(in) :: icol     ! column identifier
    real(crm_rknd)    , intent(in) :: latitude0
    real(crm_rknd)    , intent(in) :: longitude0
    integer , intent(in) :: plev     ! number of levels in parent model
    real(r8), intent(in) :: ps ! Global grid surface pressure (Pa)
    real(r8), intent(in) :: pmid(plev) ! Global grid pressure (Pa)
    real(r8), intent(in) :: pdel(plev) ! Layer's pressure thickness (Pa)
    real(r8), intent(in) :: phis ! Global grid surface geopotential (m2/s2)
    real(r8), intent(in) :: zmid(plev) ! Global grid height (m)
    real(r8), intent(in) :: zint(plev+1)! Global grid interface height (m)
    real(r8), intent(in) :: qrad_crm(crm_nx, crm_ny, crm_nz) ! CRM rad. heating
    real(r8), intent(in) :: dt_gl ! global model's time step
    real(r8), intent(in) :: ocnfrac ! area fraction of the ocean
    real(r8), intent(in) :: tau00  ! large-scale surface stress (N/m2)
    real(r8), intent(in) :: wndls  ! large-scale surface wind (m/s)
    real(r8), intent(in) :: bflxls  ! large-scale surface buoyancy flux (K m/s)
    real(r8), intent(in) :: fluxu00  ! surface momenent fluxes [N/m2]
    real(r8), intent(in) :: fluxv00  ! surface momenent fluxes [N/m2]
    real(r8), intent(in) :: fluxt00  ! surface sensible heat fluxes [K Kg/ (m2 s)]
    real(r8), intent(in) :: fluxq00  ! surface latent heat fluxes [ kg/(m2 s)]
    real(r8), intent(in) :: tl(plev) ! Global grid temperature (K)
    real(r8), intent(in) :: ql(plev) ! Global grid water vapor (g/g)
    real(r8), intent(in) :: qccl(plev)! Global grid cloud liquid water (g/g)
    real(r8), intent(in) :: qiil(plev)! Global grid cloud ice (g/g)
    real(r8), intent(in) :: ul(plev) ! Global grid u (m/s)
    real(r8), intent(in) :: vl(plev) ! Global grid v (m/s)
#ifdef CLUBB_CRM
    real(r8), intent(in), target :: clubb_buffer(crm_nx, crm_ny, crm_nz+1,1:nclubbvars)
#endif
    real(r8), intent(in) :: cltot ! shaded cloud fraction
    real(r8), intent(in) :: clhgh ! shaded cloud fraction
    real(r8), intent(in) :: clmed ! shaded cloud fraction
    real(r8), intent(in) :: cllow ! shaded cloud fraction
    real(r8), intent(in) :: u_crm  (crm_nx,crm_ny,crm_nz) ! CRM v-wind component
    real(r8), intent(in) :: v_crm  (crm_nx,crm_ny,crm_nz) ! CRM v-wind component
    real(r8), intent(in) :: w_crm  (crm_nx,crm_ny,crm_nz) ! CRM w-wind component
    real(r8), intent(in) :: t_crm  (crm_nx,crm_ny,crm_nz) ! CRM temperuture
    real(r8), intent(in) :: micro_fields_crm  (crm_nx,crm_ny,crm_nz,nmicro_fields+1) ! CRM total water
#ifdef m2005
#ifdef MODAL_AERO
    real(r8), intent(in) :: naermod(plev, ntot_amode)     ! Aerosol number concentration [/m3]
    real(r8), intent(in) :: vaerosol(plev, ntot_amode)    ! aerosol volume concentration [m3/m3]
    real(r8), intent(in) :: hygro(plev, ntot_amode)       ! hygroscopicity of aerosol mode 
#endif 
#endif
    real(r8), intent(in) :: dd_crm (plev)             ! mass entraiment from downdraft
    real(r8), intent(in) :: mui_crm (plev+1)             ! mass flux up at the interface
    real(r8), intent(in) :: mdi_crm (plev+1)             ! mass flux down at the interface

    integer :: myrank, ierr
    real(crm_rknd) :: myrand

#ifdef CRM_DUMP
        call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)
#ifdef CRM_DUMP_RATIO
        !If CRM_DUMP_RATIO is defined, then you only want to output that proportion of the time
        !So sample a random number, and if it's less than the given ratio, do the output; otherwise don't
        !Also note that crm_dump is controlled only here. the output() routine shared the same value that is
        !set by this input routine
        call random_number(myrand)
        do_dump = .false.
        if (myrand < CRM_DUMP_RATIO) do_dump = .true.
#endif
        if (do_dump) then
          !The first time we create the file, we have to make sure that we have the needed attributes for ./setup_from_file.sh in the standalone model to work properly.
          if (first_time) then
            first_time = .false.
            call dmdf_write_attr(CRM_NX   ,myrank,'crm_in','crm_nx'  ); _ERR(success,error_string,__LINE__)
            call dmdf_write_attr(CRM_NY   ,myrank,'crm_in','crm_ny'  ); _ERR(success,error_string,__LINE__)
            call dmdf_write_attr(CRM_NZ   ,myrank,'crm_in','crm_nz'  ); _ERR(success,error_string,__LINE__)
            call dmdf_write_attr(CRM_DX   ,myrank,'crm_in','crm_dx'  ); _ERR(success,error_string,__LINE__)
            call dmdf_write_attr(CRM_DT   ,myrank,'crm_in','crm_dt'  ); _ERR(success,error_string,__LINE__)
#if   (defined m2005)
            call dmdf_write_attr('m2005'  ,myrank,'crm_in','micro'   ); _ERR(success,error_string,__LINE__)
#elif (defined sam1mom)
            call dmdf_write_attr('sam1mom',myrank,'crm_in','micro'   ); _ERR(success,error_string,__LINE__)
#endif
#if   (defined _UM5)
            call dmdf_write_attr('UM5'    ,myrank,'crm_in','crm_adv' ); _ERR(success,error_string,__LINE__)
#elif (defined _MPDATA)
            call dmdf_write_attr('MPDATA' ,myrank,'crm_in','crm_adv' ); _ERR(success,error_string,__LINE__)
#endif
            call dmdf_write_attr(PLEV     ,myrank,'crm_in','plev'    ); _ERR(success,error_string,__LINE__)
            call dmdf_write_attr(PSUBCOLS ,myrank,'crm_in','psubcols'); _ERR(success,error_string,__LINE__)
            call dmdf_write_attr(PCOLS    ,myrank,'crm_in','pcols'   ); _ERR(success,error_string,__LINE__)
            call dmdf_write_attr(PCNST    ,myrank,'crm_in','pcnst'   ); _ERR(success,error_string,__LINE__)
          endif
          !Don't output the first time because fields are initially uniform, and random noise is added. Keep stochasticity out of the standalone model
          if (igstep > 1) then
            call dmdf_write(lchnk           ,myrank,'crm_in','lchnk'                                                             ,.true. ,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(icol            ,myrank,'crm_in','icol'                                                              ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(latitude0       ,myrank,'crm_in','latitude0'                                                         ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(longitude0      ,myrank,'crm_in','longitude0'                                                        ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(ps              ,myrank,'crm_in','ps'                                                                ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(pmid            ,myrank,'crm_in','pmid'            ,(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(pdel            ,myrank,'crm_in','pdel'            ,(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(phis            ,myrank,'crm_in','phis'                                                              ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(zmid            ,myrank,'crm_in','zmid'            ,(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(zint            ,myrank,'crm_in','zint'            ,(/'plev_p1'/)                                    ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(qrad_crm        ,myrank,'crm_in','qrad_crm'        ,(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(dt_gl           ,myrank,'crm_in','dt_gl'                                                             ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(ocnfrac         ,myrank,'crm_in','ocnfrac'                                                           ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(tau00           ,myrank,'crm_in','tau00'                                                             ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(wndls           ,myrank,'crm_in','wndls'                                                             ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(bflxls          ,myrank,'crm_in','bflxls'                                                            ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(fluxu00         ,myrank,'crm_in','fluxu00'                                                           ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(fluxv00         ,myrank,'crm_in','fluxv00'                                                           ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(fluxt00         ,myrank,'crm_in','fluxt00'                                                           ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(fluxq00         ,myrank,'crm_in','fluxq00'                                                           ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(tl              ,myrank,'crm_in','tl'              ,(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(ql              ,myrank,'crm_in','ql'              ,(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(qccl            ,myrank,'crm_in','qccl'            ,(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(qiil            ,myrank,'crm_in','qiil'            ,(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(ul              ,myrank,'crm_in','ul'              ,(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(vl              ,myrank,'crm_in','vl'              ,(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
#ifdef CLUBB_CRM
            call dmdf_write(clubb_buffer    ,myrank,'crm_in','clubb_buffer'    ,(/'crm_nx','crm_ny','crm_nz_p1','nclubbvars'/)   ,.false.,.false.); _ERR(success,error_string,__LINE__)
#endif
            call dmdf_write(cltot           ,myrank,'crm_in','cltot'                                                             ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(clhgh           ,myrank,'crm_in','clhgh'                                                             ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(clmed           ,myrank,'crm_in','clmed'                                                             ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(cllow           ,myrank,'crm_in','cllow'                                                             ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(u_crm           ,myrank,'crm_in','u_crm'           ,(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(v_crm           ,myrank,'crm_in','v_crm'           ,(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(w_crm           ,myrank,'crm_in','w_crm'           ,(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(t_crm           ,myrank,'crm_in','t_crm'           ,(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(micro_fields_crm,myrank,'crm_in','micro_fields_crm',(/'crm_nx','crm_ny','crm_nz','nmicro_fields_p1'/),.false.,.false.); _ERR(success,error_string,__LINE__)
#ifdef m2005
#ifdef MODAL_AERO
            call dmdf_write(naermod         ,myrank,'crm_in','naermod'         ,(/'plev','ntot_amode'/)                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(vaerosol        ,myrank,'crm_in','vaerosol'        ,(/'plev','ntot_amode'/)                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(hygro           ,myrank,'crm_in','hygro'           ,(/'plev','ntot_amode'/)                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
#endif
#endif
            call dmdf_write(dd_crm          ,myrank,'crm_in','dd_crm'          ,(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(mui_crm         ,myrank,'crm_in','mui_crm'         ,(/'plev_p1'/)                                    ,.false.,.false.); _ERR(success,error_string,__LINE__)
            call dmdf_write(mdi_crm         ,myrank,'crm_in','mdi_crm'         ,(/'plev_p1'/)                                    ,.false.,.true. ); _ERR(success,error_string,__LINE__)
          endif
        endif
#endif
  end subroutine crm_dump_input


  subroutine crm_dump_output(igstep,plev,crm_tk,crm_tkh,cltot,clhgh,clmed,cllow,sltend,u_crm,v_crm,w_crm,t_crm,micro_fields_crm,&
                             qltend,qcltend,qiltend,t_rad,qv_rad,qc_rad,qi_rad,cld_rad,cld3d_crm, &
#ifdef CLUBB_CRM
                             clubb_buffer,crm_cld,clubb_tk,clubb_tkh,relvar,accre_enhan,qclvar , &
#endif
#ifdef CRM3D
                             ultend,vltend , &
#endif
#ifdef m2005
                             nc_rad,ni_rad,qs_rad,ns_rad,wvar_crm,aut_crm,acc_crm,evpc_crm,evpr_crm,mlt_crm,sub_crm,dep_crm,con_crm,aut_crm_a,acc_crm_a,&
                             evpc_crm_a,evpr_crm_a,mlt_crm_a,sub_crm_a,dep_crm_a,con_crm_a,crm_nc,crm_ni,crm_ns,crm_ng,crm_nr, &
#endif
#ifdef ECPP
                             acen,acen_tf,rhcen,qcloudcen,qicecen,qlsinkcen,precrcen,precsolidcen,qlsink_bfcen,qlsink_avgcen,praincen,wwqui_cen,wwqui_cloudy_cen,&
                             abnd,abnd_tf,massflxbnd,wupthresh_bnd,wdownthresh_bnd,wwqui_bnd,wwqui_cloudy_bnd, &
#endif
                             precc,precl,cld,cldtop,gicewp,gliqwp,mc,mcup,mcdn,mcuup,mcudn,crm_qc,crm_qi,crm_qs,crm_qg,crm_qr,mu_crm,md_crm,du_crm,eu_crm,&
                             ed_crm,dd_crm,jt_crm,mx_crm,mui_crm,mdi_crm,flux_qt,fluxsgs_qt,tkez,tkesgsz,tkz,flux_u,flux_v,flux_qp,pflx,qt_ls,qt_trans,   &
                             qp_trans,qp_fall,qp_src,qp_evp,t_ls,prectend,precstend,precsc,precsl,taux_crm,tauy_crm,z0m,timing_factor,qc_crm,qi_crm,qpc_crm,qpi_crm,prec_crm,qtot)
    use mpi
    use dmdf
    use shr_kind_mod, only: r8 => shr_kind_r8
    use crmdims
    use microphysics, only: nmicro_fields
#ifdef ECPP
    use ecppvars,  only: NCLASS_CL, ncls_ecpp_in, NCLASS_PR
#endif
    implicit none
    integer, intent(in) :: igstep
    integer, intent(in) :: plev
#ifdef CLUBB_CRM
    real(r8), intent(in) :: clubb_buffer(crm_nx, crm_ny, crm_nz+1,1:nclubbvars)
    real(r8), intent(in) :: crm_cld(crm_nx, crm_ny, crm_nz+1)
    real(r8), intent(in) :: clubb_tk(crm_nx, crm_ny, crm_nz)
    real(r8), intent(in) :: clubb_tkh(crm_nx, crm_ny, crm_nz)
    real(r8), intent(in) :: relvar(crm_nx, crm_ny, crm_nz) 
    real(r8), intent(in) :: accre_enhan(crm_nx, crm_ny, crm_nz)
    real(r8), intent(in) :: qclvar(crm_nx, crm_ny, crm_nz)
#endif
    real(r8), intent(in) :: crm_tk(crm_nx, crm_ny, crm_nz)
    real(r8), intent(in) :: crm_tkh(crm_nx, crm_ny, crm_nz)
    real(r8), intent(in) :: cltot ! shaded cloud fraction
    real(r8), intent(in) :: clhgh ! shaded cloud fraction
    real(r8), intent(in) :: clmed ! shaded cloud fraction
    real(r8), intent(in) :: cllow ! shaded cloud fraction
#ifdef CRM3D
    real(r8), intent(in) :: ultend(plev) ! tendency of ul
    real(r8), intent(in) :: vltend(plev) ! tendency of vl
#endif
    real(r8), intent(in) :: sltend(plev) ! tendency of static energy
    real(r8), intent(in) :: u_crm  (crm_nx,crm_ny,crm_nz) ! CRM v-wind component
    real(r8), intent(in) :: v_crm  (crm_nx,crm_ny,crm_nz) ! CRM v-wind component
    real(r8), intent(in) :: w_crm  (crm_nx,crm_ny,crm_nz) ! CRM w-wind component
    real(r8), intent(in) :: t_crm  (crm_nx,crm_ny,crm_nz) ! CRM temperuture
    real(r8), intent(in) :: micro_fields_crm  (crm_nx,crm_ny,crm_nz,nmicro_fields+1) ! CRM total water
    real(r8), intent(in) :: qltend(plev) ! tendency of water vapor
    real(r8), intent(in) :: qcltend(plev)! tendency of cloud liquid water
    real(r8), intent(in) :: qiltend(plev)! tendency of cloud ice
    real(r8), intent(in) :: t_rad (crm_nx, crm_ny, crm_nz) ! rad temperuture
    real(r8), intent(in) :: qv_rad(crm_nx, crm_ny, crm_nz) ! rad vapor
    real(r8), intent(in) :: qc_rad(crm_nx, crm_ny, crm_nz) ! rad cloud water
    real(r8), intent(in) :: qi_rad(crm_nx, crm_ny, crm_nz) ! rad cloud ice
    real(r8), intent(in) :: cld_rad(crm_nx, crm_ny, crm_nz) ! rad cloud fraction 
    real(r8), intent(in) :: cld3d_crm(crm_nx, crm_ny, crm_nz) ! instant 3D cloud fraction
#ifdef m2005
    real(r8), intent(in) :: nc_rad(crm_nx, crm_ny, crm_nz) ! rad cloud droplet number (#/kg) 
    real(r8), intent(in) :: ni_rad(crm_nx, crm_ny, crm_nz) ! rad cloud ice crystal number (#/kg)
    real(r8), intent(in) :: qs_rad(crm_nx, crm_ny, crm_nz) ! rad cloud snow (kg/kg)
    real(r8), intent(in) :: ns_rad(crm_nx, crm_ny, crm_nz) ! rad cloud snow crystal number (#/kg)
    real(r8), intent(in) :: wvar_crm(crm_nx, crm_ny, crm_nz) ! vertical velocity variance (m/s)
    real(r8), intent(in) :: aut_crm(crm_nx, crm_ny, crm_nz) ! cloud water autoconversion (1/s)
    real(r8), intent(in) :: acc_crm(crm_nx, crm_ny, crm_nz) ! cloud water accretion (1/s)
    real(r8), intent(in) :: evpc_crm(crm_nx, crm_ny, crm_nz) ! cloud water evaporation (1/s)
    real(r8), intent(in) :: evpr_crm(crm_nx, crm_ny, crm_nz) ! rain evaporation (1/s)
    real(r8), intent(in) :: mlt_crm(crm_nx, crm_ny, crm_nz) ! ice, snow, graupel melting (1/s)
    real(r8), intent(in) :: sub_crm(crm_nx, crm_ny, crm_nz) ! ice, snow, graupel sublimation (1/s)
    real(r8), intent(in) :: dep_crm(crm_nx, crm_ny, crm_nz) ! ice, snow, graupel deposition (1/s)
    real(r8), intent(in) :: con_crm(crm_nx, crm_ny, crm_nz) ! cloud water condensation(1/s)
    real(r8), intent(in) :: aut_crm_a(plev) ! cloud water autoconversion (1/s)
    real(r8), intent(in) :: acc_crm_a(plev) ! cloud water accretion (1/s)
    real(r8), intent(in) :: evpc_crm_a(plev) ! cloud water evaporation (1/s)
    real(r8), intent(in) :: evpr_crm_a(plev) ! rain evaporation (1/s)
    real(r8), intent(in) :: mlt_crm_a(plev) ! ice, snow, graupel melting (1/s)
    real(r8), intent(in) :: sub_crm_a(plev) ! ice, snow, graupel sublimation (1/s)
    real(r8), intent(in) :: dep_crm_a(plev) ! ice, snow, graupel deposition (1/s)
    real(r8), intent(in) :: con_crm_a(plev) ! cloud water condensation(1/s)
#endif
    real(r8), intent(in) :: precc ! convective precip rate (m/s)
    real(r8), intent(in) :: precl ! stratiform precip rate (m/s)
    real(r8), intent(in) :: cld(plev)  ! cloud fraction
    real(r8), intent(in) :: cldtop(plev)  ! cloud top pdf
    real(r8), intent(in) :: gicewp(plev)  ! ice water path
    real(r8), intent(in) :: gliqwp(plev)  ! ice water path
    real(r8), intent(in) :: mc(plev)   ! cloud mass flux
    real(r8), intent(in) :: mcup(plev) ! updraft cloud mass flux
    real(r8), intent(in) :: mcdn(plev) ! downdraft cloud mass flux
    real(r8), intent(in) :: mcuup(plev) ! unsat updraft cloud mass flux
    real(r8), intent(in) :: mcudn(plev) ! unsat downdraft cloud mass flux
    real(r8), intent(in) :: crm_qc(plev)  ! mean cloud water
    real(r8), intent(in) :: crm_qi(plev)  ! mean cloud ice
    real(r8), intent(in) :: crm_qs(plev)  ! mean snow
    real(r8), intent(in) :: crm_qg(plev)  ! mean graupel
    real(r8), intent(in) :: crm_qr(plev)  ! mean rain
#ifdef m2005
    real(r8), intent(in) :: crm_nc(plev)  ! mean cloud water  (#/kg)
    real(r8), intent(in) :: crm_ni(plev)  ! mean cloud ice    (#/kg)
    real(r8), intent(in) :: crm_ns(plev)  ! mean snow         (#/kg)
    real(r8), intent(in) :: crm_ng(plev)  ! mean graupel      (#/kg)
    real(r8), intent(in) :: crm_nr(plev)  ! mean rain         (#/kg)
#endif
    real(r8), intent(in) :: mu_crm (plev)             ! mass flux up
    real(r8), intent(in) :: md_crm (plev)             ! mass flux down
    real(r8), intent(in) :: du_crm (plev)             ! mass detrainment from updraft
    real(r8), intent(in) :: eu_crm (plev)             ! mass entrainment from updraft
    real(r8), intent(in) :: ed_crm (plev)             ! mass detrainment from downdraft
    real(r8), intent(in) :: dd_crm (plev)             ! mass entraiment from downdraft
    real(r8), intent(in) :: jt_crm                    ! index of cloud (convection) top 
    real(r8), intent(in) :: mx_crm                    ! index of cloud (convection) bottom
    real(r8), intent(in) :: mui_crm (plev+1)             ! mass flux up at the interface
    real(r8), intent(in) :: mdi_crm (plev+1)             ! mass flux down at the interface
    real(r8), intent(in) :: flux_qt(plev) ! nonprecipitating water flux           [kg/m2/s]
    real(r8), intent(in) :: fluxsgs_qt(plev) ! sgs nonprecipitating water flux    [kg/m2/s]
    real(r8), intent(in) :: tkez(plev) ! tke profile               [kg/m/s2]
    real(r8), intent(in) :: tkesgsz(plev) ! sgs tke profile        [kg/m/s2]
    real(r8), intent(in) :: tkz(plev)  ! tk profile                [m2/s]
    real(r8), intent(in) :: flux_u(plev) ! x-momentum flux          [m2/s2]
    real(r8), intent(in) :: flux_v(plev) ! y-momentum flux          [m2/s2]
    real(r8), intent(in) :: flux_qp(plev) ! precipitating water flux [kg/m2/s or mm/s]
    real(r8), intent(in) :: pflx(plev)    ! precipitation flux      [m/s]
    real(r8), intent(in) :: qt_ls(plev) ! tendency of nonprec water due to large-scale  [kg/kg/s]
    real(r8), intent(in) :: qt_trans(plev)! tendency of nonprec water due to transport  [kg/kg/s]
    real(r8), intent(in) :: qp_trans(plev) ! tendency of prec water due to transport [kg/kg/s]
    real(r8), intent(in) :: qp_fall(plev) ! tendency of prec water due to fall-out   [kg/kg/s]
    real(r8), intent(in) :: qp_src(plev) ! tendency of prec water due to conversion  [kg/kg/s]
    real(r8), intent(in) :: qp_evp(plev) ! tendency of prec water due to evp         [kg/kg/s]
    real(r8), intent(in) :: t_ls(plev) ! tendency of lwse  due to large-scale        [kg/kg/s] ???
    real(r8), intent(in) :: prectend ! column integrated tendency in precipitating water+ice (kg/m2/s)
    real(r8), intent(in) :: precstend ! column integrated tendency in precipitating ice (kg/m2/s)
    real(r8), intent(in) :: precsc ! convective snow rate (m/s)
    real(r8), intent(in) :: precsl ! stratiform snow rate (m/s)
    real(r8), intent(in) :: taux_crm  ! zonal CRM surface stress perturbation (N/m2)
    real(r8), intent(in) :: tauy_crm  ! merid CRM surface stress perturbation (N/m2)
    real(r8), intent(in) :: z0m ! surface stress (N/m2)
    real(r8), intent(in) :: timing_factor ! crm cpu efficiency
    real(r8), intent(in) :: qc_crm (crm_nx, crm_ny, crm_nz)! CRM cloud water
    real(r8), intent(in) :: qi_crm (crm_nx, crm_ny, crm_nz)! CRM cloud ice
    real(r8), intent(in) :: qpc_crm(crm_nx, crm_ny, crm_nz)! CRM precip water
    real(r8), intent(in) :: qpi_crm(crm_nx, crm_ny, crm_nz)! CRM precip ice
    real(r8), intent(in) :: prec_crm(crm_nx, crm_ny)! CRM precipiation rate
#ifdef ECPP
    real(r8), intent(in) :: acen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)   ! cloud fraction for each sub-sub class for full time period
    real(r8), intent(in) :: acen_tf(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) ! cloud fraction for end-portion of time period
    real(r8), intent(in) :: rhcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! relative humidity (0-1)
    real(r8), intent(in) :: qcloudcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! cloud water (kg/kg)
    real(r8), intent(in) :: qicecen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) ! cloud ice (kg/kg)
    real(r8), intent(in) :: qlsinkcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! cloud water loss rate from precipitation (/s??)
    real(r8), intent(in) :: precrcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)   ! liquid (rain) precipitation rate (kg/m2/s)
    real(r8), intent(in) :: precsolidcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)   ! solid (rain) precipitation rate (kg/m2/s)
    real(r8), intent(in) :: qlsink_bfcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! cloud water loss rate from precipitation calculated 
    real(r8), intent(in) :: qlsink_avgcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! cloud water loss rate from precipitation calculated 
    real(r8), intent(in) :: praincen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! cloud water loss rate from precipitation (kg/kg/s)
    real(r8), intent(in) :: wwqui_cen(plev)                                ! vertical velocity variance in quiescent class (m2/s2)
    real(r8), intent(in) :: wwqui_cloudy_cen(plev)                         ! vertical velocity variance in quiescent, and cloudy class (m2/s2)
    real(r8), intent(in) :: abnd(plev+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)   ! cloud fraction for each sub-sub class for full time period
    real(r8), intent(in) :: abnd_tf(plev+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) ! cloud fraction for end-portion of time period
    real(r8), intent(in) :: massflxbnd(plev+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) ! sub-class vertical mass flux (kg/m2/s) at layer bottom boundary.
    real(r8), intent(in) :: wupthresh_bnd(plev+1)             ! vertical velocity threshold for updraft (m/s)
    real(r8), intent(in) :: wdownthresh_bnd(plev+1)           ! vertical velocity threshold for downdraft (m/s)
    real(r8), intent(in) :: wwqui_bnd(plev+1)                                ! vertical velocity variance in quiescent class (m2/s2)
    real(r8), intent(in) :: wwqui_cloudy_bnd(plev+1)                         ! vertical velocity variance in quiescent, and cloudy class (m2/s2)
#endif
    real(r8), intent(in) :: qtot(20)

    integer :: myrank, ierr

#ifdef CRM_DUMP
    call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)
    if (do_dump) then
      !Don't output the first time because fields are initially uniform, and random noise is added. Keep stochasticity out of the standalone model
      if (igstep > 1) then
        call dmdf_write(crm_tk          ,myrank,'crm_out',trim('crm_tk          '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.true. ,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(crm_tkh         ,myrank,'crm_out',trim('crm_tkh         '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(cltot           ,myrank,'crm_out',trim('cltot           ')                                                     ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(clhgh           ,myrank,'crm_out',trim('clhgh           ')                                                     ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(clmed           ,myrank,'crm_out',trim('clmed           ')                                                     ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(cllow           ,myrank,'crm_out',trim('cllow           ')                                                     ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(sltend          ,myrank,'crm_out',trim('sltend          '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(u_crm           ,myrank,'crm_out',trim('u_crm           '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(v_crm           ,myrank,'crm_out',trim('v_crm           '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(w_crm           ,myrank,'crm_out',trim('w_crm           '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(t_crm           ,myrank,'crm_out',trim('t_crm           '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(micro_fields_crm,myrank,'crm_out',trim('micro_fields_crm'),(/'crm_nx','crm_ny','crm_nz','nmicro_fields_p1'/)   ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(qltend          ,myrank,'crm_out',trim('qltend          '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(qcltend         ,myrank,'crm_out',trim('qcltend         '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(qiltend         ,myrank,'crm_out',trim('qiltend         '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(t_rad           ,myrank,'crm_out',trim('t_rad           '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(qv_rad          ,myrank,'crm_out',trim('qv_rad          '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(qc_rad          ,myrank,'crm_out',trim('qc_rad          '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(qi_rad          ,myrank,'crm_out',trim('qi_rad          '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(cld_rad         ,myrank,'crm_out',trim('cld_rad         '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(cld3d_crm       ,myrank,'crm_out',trim('cld3d_crm       '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
#ifdef CLUBB_CRM
        call dmdf_write(clubb_buffer    ,myrank,'crm_out',trim('clubb_buffer    '),(/'crm_nx','crm_ny','crm_nz_p1,nclubbvars'/)        ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(crm_cld         ,myrank,'crm_out',trim('crm_cld         '),(/'crm_nx','crm_ny','crm_nz_p1'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(clubb_tk        ,myrank,'crm_out',trim('clubb_tk        '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(clubb_tkh       ,myrank,'crm_out',trim('clubb_tkh       '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(relvar          ,myrank,'crm_out',trim('relvar          '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(accre_enhan     ,myrank,'crm_out',trim('accre_enhan     '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(qclvar          ,myrank,'crm_out',trim('qclvar          '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
#endif
#ifdef CRM3D
        call dmdf_write(ultend          ,myrank,'crm_out',trim('ultend          '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(vltend          ,myrank,'crm_out',trim('vltend          '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
#endif
#ifdef m2005
        call dmdf_write(nc_rad          ,myrank,'crm_out',trim('nc_rad          '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(ni_rad          ,myrank,'crm_out',trim('ni_rad          '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(qs_rad          ,myrank,'crm_out',trim('qs_rad          '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(ns_rad          ,myrank,'crm_out',trim('ns_rad          '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(wvar_crm        ,myrank,'crm_out',trim('wvar_crm        '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(aut_crm         ,myrank,'crm_out',trim('aut_crm         '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(acc_crm         ,myrank,'crm_out',trim('acc_crm         '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(evpc_crm        ,myrank,'crm_out',trim('evpc_crm        '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(evpr_crm        ,myrank,'crm_out',trim('evpr_crm        '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(mlt_crm         ,myrank,'crm_out',trim('mlt_crm         '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(sub_crm         ,myrank,'crm_out',trim('sub_crm         '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(dep_crm         ,myrank,'crm_out',trim('dep_crm         '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(con_crm         ,myrank,'crm_out',trim('con_crm         '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(aut_crm_a       ,myrank,'crm_out',trim('aut_crm_a       '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(acc_crm_a       ,myrank,'crm_out',trim('acc_crm_a       '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(evpc_crm_a      ,myrank,'crm_out',trim('evpc_crm_a      '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(evpr_crm_a      ,myrank,'crm_out',trim('evpr_crm_a      '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(mlt_crm_a       ,myrank,'crm_out',trim('mlt_crm_a       '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(sub_crm_a       ,myrank,'crm_out',trim('sub_crm_a       '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(dep_crm_a       ,myrank,'crm_out',trim('dep_crm_a       '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(con_crm_a       ,myrank,'crm_out',trim('con_crm_a       '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
#endif
#ifdef m2005
        call dmdf_write(crm_nc          ,myrank,'crm_out',trim('crm_nc          '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(crm_ni          ,myrank,'crm_out',trim('crm_ni          '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(crm_ns          ,myrank,'crm_out',trim('crm_ns          '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(crm_ng          ,myrank,'crm_out',trim('crm_ng          '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(crm_nr          ,myrank,'crm_out',trim('crm_nr          '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
#endif
#ifdef ECPP
        call dmdf_write(acen            ,myrank,'crm_out',trim('acen            '),(/'plev','NCLASS_CL','ncls_ecpp_in','NCLASS_PR'/)   ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(acen_tf         ,myrank,'crm_out',trim('acen_tf         '),(/'plev','NCLASS_CL','ncls_ecpp_in','NCLASS_PR'/)   ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(rhcen           ,myrank,'crm_out',trim('rhcen           '),(/'plev','NCLASS_CL','ncls_ecpp_in','NCLASS_PR'/)   ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(qcloudcen       ,myrank,'crm_out',trim('qcloudcen       '),(/'plev','NCLASS_CL','ncls_ecpp_in','NCLASS_PR'/)   ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(qicecen         ,myrank,'crm_out',trim('qicecen         '),(/'plev','NCLASS_CL','ncls_ecpp_in','NCLASS_PR'/)   ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(qlsinkcen       ,myrank,'crm_out',trim('qlsinkcen       '),(/'plev','NCLASS_CL','ncls_ecpp_in','NCLASS_PR'/)   ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(precrcen        ,myrank,'crm_out',trim('precrcen        '),(/'plev','NCLASS_CL','ncls_ecpp_in','NCLASS_PR'/)   ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(precsolidcen    ,myrank,'crm_out',trim('precsolidcen    '),(/'plev','NCLASS_CL','ncls_ecpp_in','NCLASS_PR'/)   ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(qlsink_bfcen    ,myrank,'crm_out',trim('qlsink_bfcen    '),(/'plev','NCLASS_CL','ncls_ecpp_in','NCLASS_PR'/)   ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(qlsink_avgcen   ,myrank,'crm_out',trim('qlsink_avgcen   '),(/'plev','NCLASS_CL','ncls_ecpp_in','NCLASS_PR'/)   ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(praincen        ,myrank,'crm_out',trim('praincen        '),(/'plev','NCLASS_CL','ncls_ecpp_in','NCLASS_PR'/)   ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(wwqui_cen       ,myrank,'crm_out',trim('wwqui_cen       '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(wwqui_cloudy_cen,myrank,'crm_out',trim('wwqui_cloudy_cen'),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(abnd            ,myrank,'crm_out',trim('abnd            '),(/'plev_p1','NCLASS_CL','ncls_ecpp_in','NCLASS_PR'/),.false.,.false.); _ERR(success,error_string,__LINE__)  
        call dmdf_write(abnd_tf         ,myrank,'crm_out',trim('abnd_tf         '),(/'plev_p1','NCLASS_CL','ncls_ecpp_in','NCLASS_PR'/),.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(massflxbnd      ,myrank,'crm_out',trim('massflxbnd      '),(/'plev_p1','NCLASS_CL','ncls_ecpp_in','NCLASS_PR'/),.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(wupthresh_bnd   ,myrank,'crm_out',trim('wupthresh_bnd   '),(/'plev_p1'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(wdownthresh_bnd ,myrank,'crm_out',trim('wdownthresh_bnd '),(/'plev_p1'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(wwqui_bnd       ,myrank,'crm_out',trim('wwqui_bnd       '),(/'plev_p1'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(wwqui_cloudy_bnd,myrank,'crm_out',trim('wwqui_cloudy_bnd'),(/'plev_p1'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
#endif
        call dmdf_write(precc           ,myrank,'crm_out',trim('precc           ')                                                     ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(precl           ,myrank,'crm_out',trim('precl           ')                                                     ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(cld             ,myrank,'crm_out',trim('cld             '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(cldtop          ,myrank,'crm_out',trim('cldtop          '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(gicewp          ,myrank,'crm_out',trim('gicewp          '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(gliqwp          ,myrank,'crm_out',trim('gliqwp          '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(mc              ,myrank,'crm_out',trim('mc              '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(mcup            ,myrank,'crm_out',trim('mcup            '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(mcdn            ,myrank,'crm_out',trim('mcdn            '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(mcuup           ,myrank,'crm_out',trim('mcuup           '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(mcudn           ,myrank,'crm_out',trim('mcudn           '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(crm_qc          ,myrank,'crm_out',trim('crm_qc          '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(crm_qi          ,myrank,'crm_out',trim('crm_qi          '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(crm_qs          ,myrank,'crm_out',trim('crm_qs          '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(crm_qg          ,myrank,'crm_out',trim('crm_qg          '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(crm_qr          ,myrank,'crm_out',trim('crm_qr          '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(mu_crm          ,myrank,'crm_out',trim('mu_crm          '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(md_crm          ,myrank,'crm_out',trim('md_crm          '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(du_crm          ,myrank,'crm_out',trim('du_crm          '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(eu_crm          ,myrank,'crm_out',trim('eu_crm          '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(ed_crm          ,myrank,'crm_out',trim('ed_crm          '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(dd_crm          ,myrank,'crm_out',trim('dd_crm          '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(jt_crm          ,myrank,'crm_out',trim('jt_crm          ')                                                     ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(mx_crm          ,myrank,'crm_out',trim('mx_crm          ')                                                     ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(mui_crm         ,myrank,'crm_out',trim('mui_crm         '),(/'plev_p1'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(mdi_crm         ,myrank,'crm_out',trim('mdi_crm         '),(/'plev_p1'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(flux_qt         ,myrank,'crm_out',trim('flux_qt         '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(fluxsgs_qt      ,myrank,'crm_out',trim('fluxsgs_qt      '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(tkez            ,myrank,'crm_out',trim('tkez            '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(tkesgsz         ,myrank,'crm_out',trim('tkesgsz         '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(tkz             ,myrank,'crm_out',trim('tkz             '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(flux_u          ,myrank,'crm_out',trim('flux_u          '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(flux_v          ,myrank,'crm_out',trim('flux_v          '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(flux_qp         ,myrank,'crm_out',trim('flux_qp         '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(pflx            ,myrank,'crm_out',trim('pflx            '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(qt_ls           ,myrank,'crm_out',trim('qt_ls           '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(qt_trans        ,myrank,'crm_out',trim('qt_trans        '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(qp_trans        ,myrank,'crm_out',trim('qp_trans        '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(qp_fall         ,myrank,'crm_out',trim('qp_fall         '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(qp_src          ,myrank,'crm_out',trim('qp_src          '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(qp_evp          ,myrank,'crm_out',trim('qp_evp          '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(t_ls            ,myrank,'crm_out',trim('t_ls            '),(/'plev'/)                                          ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(prectend        ,myrank,'crm_out',trim('prectend        ')                                                     ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(precstend       ,myrank,'crm_out',trim('precstend       ')                                                     ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(precsc          ,myrank,'crm_out',trim('precsc          ')                                                     ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(precsl          ,myrank,'crm_out',trim('precsl          ')                                                     ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(taux_crm        ,myrank,'crm_out',trim('taux_crm        ')                                                     ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(tauy_crm        ,myrank,'crm_out',trim('tauy_crm        ')                                                     ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(z0m             ,myrank,'crm_out',trim('z0m             ')                                                     ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(timing_factor   ,myrank,'crm_out',trim('timing_factor   ')                                                     ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(qc_crm          ,myrank,'crm_out',trim('qc_crm          '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(qi_crm          ,myrank,'crm_out',trim('qi_crm          '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(qpc_crm         ,myrank,'crm_out',trim('qpc_crm         '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(qpi_crm         ,myrank,'crm_out',trim('qpi_crm         '),(/'crm_nx','crm_ny','crm_nz'/)                      ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(prec_crm        ,myrank,'crm_out',trim('prec_crm        '),(/'crm_nx','crm_ny'/)                               ,.false.,.false.); _ERR(success,error_string,__LINE__)
        call dmdf_write(qtot            ,myrank,'crm_out',trim('qtot            '),(/'d20'/)                                           ,.false.,.true. ); _ERR(success,error_string,__LINE__)
      endif
    endif
#endif
  end subroutine crm_dump_output


end module crm_dump
