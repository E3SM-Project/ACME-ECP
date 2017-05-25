
#ifndef CRM_STANDALONE
#define CRM_STANDALONE
#endif

!MRN: Convenient error checking macro for dumping data to file
#define _ERR(s,e,l) if (.not. s) then; write(*,*) l,', ',trim(e); stop; endif

program crm_standalone
  use mpi
  use dmdf
  use crm_dump    , only: crm_dump_output
  use perf_mod    , only: t_startf, t_stopf, t_prf, t_initf, t_finalizef
  use crm_module  , only: crm
  use crmdims     , only: nclubbvars, crm_nx, crm_ny, crm_nz
  use shr_kind_mod, only: r8 => shr_kind_r8
  use microphysics, only: nmicro_fields
  use setparm_mod , only: setparm
#ifdef ECPP
  use ecppvars,  only: NCLASS_CL, ncls_ecpp_in, NCLASS_PR
#endif
  implicit none

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Standalone Variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  character(len=512) :: fname_in, fprefix_out
  integer :: col1, col2, ncols, stat, i, ind, ierr, nsamp_tot
  integer :: rank, nranks, globid
  integer, parameter :: plev = PLEV
  real(r8) :: nper, t1, t2, tmin, tmax, tred(2)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! CRM Variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer , allocatable :: lchnk(:) 
  integer , allocatable :: icol(:) 
#ifdef CRM_STANDALONE
  real    , allocatable :: latitude0(:) 
  real    , allocatable :: longitude0(:) 
#endif
  real(r8), allocatable :: ps(:) 
  real(r8), allocatable :: pmid(:,:) 
  real(r8), allocatable :: pdel(:,:) 
  real(r8), allocatable :: phis(:) 
  real(r8), allocatable :: zmid(:,:) 
  real(r8), allocatable :: zint(:,:)
  real(r8), allocatable :: qrad_crm(:, :, :,:) 
  real(r8), allocatable :: dt_gl(:) 
  real(r8), allocatable :: ocnfrac(:) 
  real(r8), allocatable :: tau00(:) 
  real(r8), allocatable :: wndls(:) 
  real(r8), allocatable :: bflxls(:) 
  real(r8), allocatable :: fluxu00(:) 
  real(r8), allocatable :: fluxv00(:) 
  real(r8), allocatable :: fluxt00(:) 
  real(r8), allocatable :: fluxq00(:) 
  real(r8), allocatable :: tl(:,:) 
  real(r8), allocatable :: ql(:,:) 
  real(r8), allocatable :: qccl(:,:)
  real(r8), allocatable :: qiil(:,:)
  real(r8), allocatable :: ul(:,:) 
  real(r8), allocatable :: vl(:,:) 
#ifdef CLUBB_CRM
  real(r8), allocatable, target :: clubb_buffer(:, :, :,:,:)
  real(r8), allocatable :: crm_cld(:, :, :,:)
  real(r8), allocatable :: clubb_tk(:, :, :,:)
  real(r8), allocatable :: clubb_tkh(:, :, :,:)
  real(r8), allocatable :: relvar(:, :, :,:) 
  real(r8), allocatable :: accre_enhan(:, :, :,:)
  real(r8), allocatable :: qclvar(:, :, :,:)
#endif
  real(r8), allocatable :: crm_tk(:, :, :,:)
  real(r8), allocatable :: crm_tkh(:, :, :,:)
  real(r8), allocatable :: cltot(:) 
  real(r8), allocatable :: clhgh(:) 
  real(r8), allocatable :: clmed(:) 
  real(r8), allocatable :: cllow(:) 
#ifdef CRM3D
  real(r8), allocatable :: ultend(:,:) 
  real(r8), allocatable :: vltend(:,:) 
#endif
  real(r8), allocatable :: sltend(:,:) 
  real(r8), allocatable :: u_crm(:,:,:,:) 
  real(r8), allocatable :: v_crm(:,:,:,:) 
  real(r8), allocatable :: w_crm(:,:,:,:) 
  real(r8), allocatable :: t_crm(:,:,:,:) 
  real(r8), allocatable :: micro_fields_crm(:,:,:,:,:)
  real(r8), allocatable :: qltend(:,:) 
  real(r8), allocatable :: qcltend(:,:)
  real(r8), allocatable :: qiltend(:,:)
  real(r8), allocatable :: t_rad(:, :, :,:) 
  real(r8), allocatable :: qv_rad(:, :, :,:) 
  real(r8), allocatable :: qc_rad(:, :, :,:) 
  real(r8), allocatable :: qi_rad(:, :, :,:) 
  real(r8), allocatable :: cld_rad(:, :, :,:) 
  real(r8), allocatable :: cld3d_crm(:, :, :,:) 
#ifdef m2005
  real(r8), allocatable :: nc_rad(:, :, :,:) 
  real(r8), allocatable :: ni_rad(:, :, :,:) 
  real(r8), allocatable :: qs_rad(:, :, :,:) 
  real(r8), allocatable :: ns_rad(:, :, :,:) 
  real(r8), allocatable :: wvar_crm(:, :, :,:) 
  real(r8), allocatable :: aut_crm(:, :, :,:) 
  real(r8), allocatable :: acc_crm(:, :, :,:) 
  real(r8), allocatable :: evpc_crm(:, :, :,:) 
  real(r8), allocatable :: evpr_crm(:, :, :,:) 
  real(r8), allocatable :: mlt_crm(:, :, :,:) 
  real(r8), allocatable :: sub_crm(:, :, :,:) 
  real(r8), allocatable :: dep_crm(:, :, :,:) 
  real(r8), allocatable :: con_crm(:, :, :,:) 
  real(r8), allocatable :: aut_crm_a(:,:) 
  real(r8), allocatable :: acc_crm_a(:,:) 
  real(r8), allocatable :: evpc_crm_a(:,:) 
  real(r8), allocatable :: evpr_crm_a(:,:) 
  real(r8), allocatable :: mlt_crm_a(:,:) 
  real(r8), allocatable :: sub_crm_a(:,:) 
  real(r8), allocatable :: dep_crm_a(:,:) 
  real(r8), allocatable :: con_crm_a(:,:) 
#endif
  real(r8), allocatable :: precc(:) 
  real(r8), allocatable :: precl(:) 
  real(r8), allocatable :: cld(:,:)  
  real(r8), allocatable :: cldtop(:,:)  
  real(r8), allocatable :: gicewp(:,:)  
  real(r8), allocatable :: gliqwp(:,:)  
  real(r8), allocatable :: mc(:,:)   
  real(r8), allocatable :: mcup(:,:) 
  real(r8), allocatable :: mcdn(:,:) 
  real(r8), allocatable :: mcuup(:,:) 
  real(r8), allocatable :: mcudn(:,:) 
  real(r8), allocatable :: crm_qc(:,:)  
  real(r8), allocatable :: crm_qi(:,:)  
  real(r8), allocatable :: crm_qs(:,:)  
  real(r8), allocatable :: crm_qg(:,:)  
  real(r8), allocatable :: crm_qr(:,:)  
#ifdef m2005
  real(r8), allocatable :: crm_nc(:,:)  
  real(r8), allocatable :: crm_ni(:,:)  
  real(r8), allocatable :: crm_ns(:,:)  
  real(r8), allocatable :: crm_ng(:,:)  
  real(r8), allocatable :: crm_nr(:,:)  
#ifdef MODAL_AERO
  real(r8), allocatable :: naermod(:, :,:)     
  real(r8), allocatable :: vaerosol(:, :,:)    
  real(r8), allocatable :: hygro(:, :,:)       
#endif 
#endif
  real(r8), allocatable :: mu_crm(:,:)             
  real(r8), allocatable :: md_crm(:,:)             
  real(r8), allocatable :: du_crm(:,:)             
  real(r8), allocatable :: eu_crm(:,:)             
  real(r8), allocatable :: ed_crm(:,:)             
  real(r8), allocatable :: dd_crm(:,:)             
  real(r8), allocatable :: jt_crm(:)                   
  real(r8), allocatable :: mx_crm(:)                   
  real(r8), allocatable :: mui_crm(:,:)             
  real(r8), allocatable :: mdi_crm(:,:)             
  real(r8), allocatable :: flux_qt(:,:) 
  real(r8), allocatable :: fluxsgs_qt(:,:) 
  real(r8), allocatable :: tkez(:,:) 
  real(r8), allocatable :: tkesgsz(:,:) 
  real(r8), allocatable :: tkz(:,:)  
  real(r8), allocatable :: flux_u(:,:) 
  real(r8), allocatable :: flux_v(:,:) 
  real(r8), allocatable :: flux_qp(:,:) 
  real(r8), allocatable :: pflx(:,:)    
  real(r8), allocatable :: qt_ls(:,:) 
  real(r8), allocatable :: qt_trans(:,:)
  real(r8), allocatable :: qp_trans(:,:) 
  real(r8), allocatable :: qp_fall(:,:) 
  real(r8), allocatable :: qp_src(:,:) 
  real(r8), allocatable :: qp_evp(:,:) 
  real(r8), allocatable :: t_ls(:,:) 
  real(r8), allocatable :: prectend(:) 
  real(r8), allocatable :: precstend(:) 
  real(r8), allocatable :: precsc(:) 
  real(r8), allocatable :: precsl(:) 
  real(r8), allocatable :: taux_crm(:) 
  real(r8), allocatable :: tauy_crm(:) 
  real(r8), allocatable :: z0m(:) 
  real(r8), allocatable :: timing_factor(:) 
  real(r8), allocatable :: qc_crm(:, :, :,:)
  real(r8), allocatable :: qi_crm(:, :, :,:)
  real(r8), allocatable :: qpc_crm(:, :, :,:)
  real(r8), allocatable :: qpi_crm(:, :, :,:)
  real(r8), allocatable :: prec_crm(:, :,:)
#ifdef ECPP
  real(r8), allocatable :: acen(:,:,:,:,:)   
  real(r8), allocatable :: acen_tf(:,:,:,:,:) 
  real(r8), allocatable :: rhcen(:,:,:,:,:)  
  real(r8), allocatable :: qcloudcen(:,:,:,:,:)  
  real(r8), allocatable :: qicecen(:,:,:,:,:) 
  real(r8), allocatable :: qlsinkcen(:,:,:,:,:)  
  real(r8), allocatable :: precrcen(:,:,:,:,:)   
  real(r8), allocatable :: precsolidcen(:,:,:,:,:)   
  real(r8), allocatable :: qlsink_bfcen(:,:,:,:,:)  
  real(r8), allocatable :: qlsink_avgcen(:,:,:,:,:)  
  real(r8), allocatable :: praincen(:,:,:,:,:)  
  real(r8), allocatable :: wwqui_cen(:,:)                                
  real(r8), allocatable :: wwqui_cloudy_cen(:,:)                         
  real(r8), allocatable :: abnd(:,:,:,:,:)   
  real(r8), allocatable :: abnd_tf(:,:,:,:,:) 
  real(r8), allocatable :: massflxbnd(:,:,:,:,:) 
  real(r8), allocatable :: wupthresh_bnd(:,:)             
  real(r8), allocatable :: wdownthresh_bnd(:,:)           
  real(r8), allocatable :: wwqui_bnd(:,:)                                
  real(r8), allocatable :: wwqui_cloudy_bnd(:,:)                         
#endif
  real(r8), allocatable :: qtot(:,:)

  !Get MPI information
  call MPI_Init( ierr )
  call MPI_Comm_size( MPI_COMM_WORLD , nranks , ierr )
  call MPI_Comm_rank( MPI_COMM_WORLD , rank   , ierr )

  call t_initf('gptl.nl',LogPrint=(rank==0),Mpicom=MPI_COMM_WORLD,MasterTask=(rank==0))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Read in the command line arguments
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call getarg(1,fname_in )
  call getarg(2,fprefix_out )
  !Get the total number of netcdf records
  call dmdf_num_records(fname_in,nsamp_tot)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Compute total number of columns per MPI task
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  globid = rank
  nper = real(nsamp_tot)/nranks
  col1 = nint( nper* globid    )+1
  col2 = nint( nper*(globid+1) )
  ncols = col2 - col1 + 1

  if (rank == 0) write(unit=*,fmt='(2(A,I7),A)') 'Processing: ',ncols,' of ',nsamp_tot,' samples.'

  !I have to call setparm to get the correct value for nmicro_fields
  call setparm()

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Allocate CRM variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate( lchnk(ncols) ) 
  allocate( icol(ncols) ) 
#ifdef CRM_STANDALONE
  allocate( latitude0(ncols) ) 
  allocate( longitude0(ncols) ) 
#endif
  allocate( ps(ncols) ) 
  allocate( pmid(plev,ncols) ) 
  allocate( pdel(plev,ncols) ) 
  allocate( phis(ncols) ) 
  allocate( zmid(plev,ncols) ) 
  allocate( zint(plev+1,ncols) )
  allocate( qrad_crm(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( dt_gl(ncols) ) 
  allocate( ocnfrac(ncols) ) 
  allocate( tau00(ncols) ) 
  allocate( wndls(ncols) ) 
  allocate( bflxls(ncols) ) 
  allocate( fluxu00(ncols) ) 
  allocate( fluxv00(ncols) ) 
  allocate( fluxt00(ncols) ) 
  allocate( fluxq00(ncols) ) 
  allocate( tl(plev,ncols) ) 
  allocate( ql(plev,ncols) ) 
  allocate( qccl(plev,ncols) )
  allocate( qiil(plev,ncols) )
  allocate( ul(plev,ncols) ) 
  allocate( vl(plev,ncols) ) 
#ifdef CLUBB_CRM
  allocate( clubb_buffer(crm_nx, crm_ny, crm_nz+1,1:nclubbvars,ncols) )
  allocate( crm_cld(crm_nx, crm_ny, crm_nz+1,ncols) )
  allocate( clubb_tk(crm_nx, crm_ny, crm_nz,ncols) )
  allocate( clubb_tkh(crm_nx, crm_ny, crm_nz,ncols) )
  allocate( relvar(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( accre_enhan(crm_nx, crm_ny, crm_nz,ncols) )
  allocate( qclvar(crm_nx, crm_ny, crm_nz,ncols) )
#endif
  allocate( crm_tk(crm_nx, crm_ny, crm_nz,ncols) )
  allocate( crm_tkh(crm_nx, crm_ny, crm_nz,ncols) )
  allocate( cltot(ncols) ) 
  allocate( clhgh(ncols) ) 
  allocate( clmed(ncols) ) 
  allocate( cllow(ncols) ) 
#ifdef CRM3D
  allocate( ultend(plev,ncols) ) 
  allocate( vltend(plev,ncols) ) 
#endif
  allocate( sltend(plev,ncols) ) 
  allocate( u_crm(crm_nx,crm_ny,crm_nz,ncols) ) 
  allocate( v_crm(crm_nx,crm_ny,crm_nz,ncols) ) 
  allocate( w_crm(crm_nx,crm_ny,crm_nz,ncols) ) 
  allocate( t_crm(crm_nx,crm_ny,crm_nz,ncols) ) 
  allocate( micro_fields_crm(crm_nx,crm_ny,crm_nz,nmicro_fields+1,ncols) )
  allocate( qltend(plev,ncols) ) 
  allocate( qcltend(plev,ncols) )
  allocate( qiltend(plev,ncols) )
  allocate( t_rad(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( qv_rad(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( qc_rad(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( qi_rad(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( cld_rad(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( cld3d_crm(crm_nx, crm_ny, crm_nz,ncols) ) 
#ifdef m2005
  allocate( nc_rad(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( ni_rad(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( qs_rad(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( ns_rad(crm_nx, crm_ny, crm_nz,ncols) )
  allocate( wvar_crm(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( aut_crm(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( acc_crm(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( evpc_crm(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( evpr_crm(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( mlt_crm(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( sub_crm(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( dep_crm(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( con_crm(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( aut_crm_a(plev,ncols) ) 
  allocate( acc_crm_a(plev,ncols) ) 
  allocate( evpc_crm_a(plev,ncols) ) 
  allocate( evpr_crm_a(plev,ncols) ) 
  allocate( mlt_crm_a(plev,ncols) ) 
  allocate( sub_crm_a(plev,ncols) ) 
  allocate( dep_crm_a(plev,ncols) ) 
  allocate( con_crm_a(plev,ncols) ) 
#endif
  allocate( precc(ncols) ) 
  allocate( precl(ncols) ) 
  allocate( cld(plev,ncols) )  
  allocate( cldtop(plev,ncols) )  
  allocate( gicewp(plev,ncols) )  
  allocate( gliqwp(plev,ncols) )  
  allocate( mc(plev,ncols) )   
  allocate( mcup(plev,ncols) ) 
  allocate( mcdn(plev,ncols) ) 
  allocate( mcuup(plev,ncols) ) 
  allocate( mcudn(plev,ncols) ) 
  allocate( crm_qc(plev,ncols) )  
  allocate( crm_qi(plev,ncols) )  
  allocate( crm_qs(plev,ncols) )  
  allocate( crm_qg(plev,ncols) )  
  allocate( crm_qr(plev,ncols) )  
#ifdef m2005
  allocate( crm_nc(plev,ncols) )  
  allocate( crm_ni(plev,ncols) )  
  allocate( crm_ns(plev,ncols) )  
  allocate( crm_ng(plev,ncols) )  
  allocate( crm_nr(plev,ncols) )  
#ifdef MODAL_AERO
  allocate( naermod(plev, ntot_amode,ncols) )     
  allocate( vaerosol(plev, ntot_amode,ncols) )    
  allocate( hygro(plev, ntot_amode,ncols) )       
#endif 
#endif
  allocate( mu_crm(plev,ncols) )             
  allocate( md_crm(plev,ncols) )             
  allocate( du_crm(plev,ncols) )             
  allocate( eu_crm(plev,ncols) )             
  allocate( ed_crm(plev,ncols) )             
  allocate( dd_crm(plev,ncols) )             
  allocate( jt_crm(ncols) )                   
  allocate( mx_crm(ncols) )                   
  allocate( mui_crm(plev+1,ncols) )             
  allocate( mdi_crm(plev+1,ncols) )             
  allocate( flux_qt(plev,ncols) ) 
  allocate( fluxsgs_qt(plev,ncols) ) 
  allocate( tkez(plev,ncols) ) 
  allocate( tkesgsz(plev,ncols) ) 
  allocate( tkz(plev,ncols) )  
  allocate( flux_u(plev,ncols) ) 
  allocate( flux_v(plev,ncols) ) 
  allocate( flux_qp(plev,ncols) ) 
  allocate( pflx(plev,ncols) )    
  allocate( qt_ls(plev,ncols) ) 
  allocate( qt_trans(plev,ncols) )
  allocate( qp_trans(plev,ncols) ) 
  allocate( qp_fall(plev,ncols) ) 
  allocate( qp_src(plev,ncols) ) 
  allocate( qp_evp(plev,ncols) ) 
  allocate( t_ls(plev,ncols) ) 
  allocate( prectend(ncols) ) 
  allocate( precstend(ncols) ) 
  allocate( precsc(ncols) ) 
  allocate( precsl(ncols) ) 
  allocate( taux_crm(ncols) ) 
  allocate( tauy_crm(ncols) ) 
  allocate( z0m(ncols) ) 
  allocate( timing_factor(ncols) ) 
  allocate( qc_crm(crm_nx, crm_ny, crm_nz,ncols) )
  allocate( qi_crm(crm_nx, crm_ny, crm_nz,ncols) )
  allocate( qpc_crm(crm_nx, crm_ny, crm_nz,ncols) )
  allocate( qpi_crm(crm_nx, crm_ny, crm_nz,ncols) )
  allocate( prec_crm(crm_nx, crm_ny,ncols) )
#ifdef ECPP
  allocate( acen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,ncols) )   
  allocate( acen_tf(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,ncols) ) 
  allocate( rhcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,ncols) )  
  allocate( qcloudcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,ncols) )  
  allocate( qicecen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,ncols) ) 
  allocate( qlsinkcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,ncols) )  
  allocate( precrcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,ncols) )   
  allocate( precsolidcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,ncols) )   
  allocate( qlsink_bfcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,ncols) )  
  allocate( qlsink_avgcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,ncols) )  
  allocate( praincen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,ncols) )  
  allocate( wwqui_cen(plev,ncols) )                                
  allocate( wwqui_cloudy_cen(plev,ncols) )                         
  allocate( abnd(plev+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,ncols) )   
  allocate( abnd_tf(plev+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,ncols) ) 
  allocate( massflxbnd(plev+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,ncols) ) 
  allocate( wupthresh_bnd(plev+1,ncols) )             
  allocate( wdownthresh_bnd(plev+1,ncols) )           
  allocate( wwqui_bnd(plev+1,ncols) )                                
  allocate( wwqui_cloudy_bnd(plev+1,ncols) )                         
#endif
  allocate( qtot(20,ncols) )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Read in Input Data from the dump file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (rank == 0) write(*,*) '*** READING INPUT DATA ***'
  call dmdf_read(lchnk           ,trim(fname_in),'lchnk'           ,col1,col2,.true. ,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(icol            ,trim(fname_in),'icol'            ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(latitude0       ,trim(fname_in),'latitude0'       ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(longitude0      ,trim(fname_in),'longitude0'      ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(ps              ,trim(fname_in),'ps'              ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(pmid            ,trim(fname_in),'pmid'            ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(pdel            ,trim(fname_in),'pdel'            ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(phis            ,trim(fname_in),'phis'            ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(zmid            ,trim(fname_in),'zmid'            ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(zint            ,trim(fname_in),'zint'            ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(qrad_crm        ,trim(fname_in),'qrad_crm'        ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(dt_gl           ,trim(fname_in),'dt_gl'           ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(ocnfrac         ,trim(fname_in),'ocnfrac'         ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(tau00           ,trim(fname_in),'tau00'           ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(wndls           ,trim(fname_in),'wndls'           ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(bflxls          ,trim(fname_in),'bflxls'          ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(fluxu00         ,trim(fname_in),'fluxu00'         ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(fluxv00         ,trim(fname_in),'fluxv00'         ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(fluxt00         ,trim(fname_in),'fluxt00'         ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(fluxq00         ,trim(fname_in),'fluxq00'         ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(tl              ,trim(fname_in),'tl'              ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(ql              ,trim(fname_in),'ql'              ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(qccl            ,trim(fname_in),'qccl'            ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(qiil            ,trim(fname_in),'qiil'            ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(ul              ,trim(fname_in),'ul'              ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(vl              ,trim(fname_in),'vl'              ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
#ifdef CLUBB_CRM
  call dmdf_read(clubb_buffer    ,trim(fname_in),'clubb_buffer'    ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
#endif
  call dmdf_read(cltot           ,trim(fname_in),'cltot'           ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(clhgh           ,trim(fname_in),'clhgh'           ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(clmed           ,trim(fname_in),'clmed'           ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(cllow           ,trim(fname_in),'cllow'           ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(u_crm           ,trim(fname_in),'u_crm'           ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(v_crm           ,trim(fname_in),'v_crm'           ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(w_crm           ,trim(fname_in),'w_crm'           ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(t_crm           ,trim(fname_in),'t_crm'           ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(micro_fields_crm,trim(fname_in),'micro_fields_crm',col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
#ifdef m2005
#ifdef MODAL_AERO
  call dmdf_read(naermod         ,trim(fname_in),'naermod'         ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(vaerosol        ,trim(fname_in),'vaerosol'        ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(hygro           ,trim(fname_in),'hygro'           ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
#endif
#endif
  call dmdf_read(dd_crm          ,trim(fname_in),'dd_crm'          ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(mui_crm         ,trim(fname_in),'mui_crm'         ,col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(mdi_crm         ,trim(fname_in),'mdi_crm'         ,col1,col2,.false.,.true. ); _ERR(success,error_string,__LINE__)

  call MPI_Barrier( MPI_COMM_WORLD , ierr )

  if (rank == 0) write(*,*) '*** FINISHED READING INPUT DATA ***'

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Call the crm routine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !TODO: OpenMP is currently broken in SAM. I suspect it's because they're using module data as common blocks, and
  !TODO: those variables need to be declared threadprivate. It won't be fun to find all of them :-/.
  !!!$omp parallel do
  tmin = 1e10
  tmax = -1.
  do ind = 1 , ncols
    write(unit=*,fmt='(A,I7,A,I6,A)') '*** CALLING CRM ',ind,' from task ',rank,' ***'
    call t_startf('crm_main')
    t1 = MPI_Wtime()
    call crm         ( lchnk(ind), icol(ind), &
#ifdef CRM_STANDALONE
                       latitude0(ind), longitude0(ind), &
#endif
                       tl(:,ind), ql(:,ind), qccl(:,ind), qiil(:,ind), ul(:,ind), vl(:,ind), &
                       ps(ind), pmid(:,ind), pdel(:,ind), phis(ind), &
                       zmid(:,ind), zint(:,ind), dt_gl(ind), plev, &
#ifdef CRM3D
                       ultend(:,ind), vltend(:,ind),          &
#endif
                       qltend(:,ind), qcltend(:,ind), qiltend(:,ind), sltend(:,ind), &
                       u_crm(:,:,:,ind), v_crm(:,:,:,ind), w_crm(:,:,:,ind), t_crm(:,:,:,ind), micro_fields_crm(:,:,:,:,ind), &
                       qrad_crm(:,:,:,ind), &
                       qc_crm(:,:,:,ind), qi_crm(:,:,:,ind), qpc_crm(:,:,:,ind), qpi_crm(:,:,:,ind), prec_crm(:,:,ind), &
                       t_rad(:,:,:,ind), qv_rad(:,:,:,ind), qc_rad(:,:,:,ind), qi_rad(:,:,:,ind), cld_rad(:,:,:,ind), cld3d_crm(:,:,:,ind), &
#ifdef m2005
                       nc_rad(:,:,:,ind), ni_rad(:,:,:,ind), qs_rad(:,:,:,ind), ns_rad(:,:,:,ind), wvar_crm(:,:,:,ind),  &
                       aut_crm(:,:,:,ind), acc_crm(:,:,:,ind), evpc_crm(:,:,:,ind), evpr_crm(:,:,:,ind), mlt_crm(:,:,:,ind), &
                       sub_crm(:,:,:,ind), dep_crm(:,:,:,ind), con_crm(:,:,:,ind), &
                       aut_crm_a(:,ind), acc_crm_a(:,ind), evpc_crm_a(:,ind), evpr_crm_a(:,ind), mlt_crm_a(:,ind), &
                       sub_crm_a(:,ind), dep_crm_a(:,ind), con_crm_a(:,ind), &
#endif
                       precc(ind), precl(ind), precsc(ind), precsl(ind), &
                       cltot(ind), clhgh(ind), clmed(ind), cllow(ind), cld(:,ind), cldtop(:,ind), &
                       gicewp(:,ind), gliqwp(:,ind), &
                       mc(:,ind), mcup(:,ind), mcdn(:,ind), mcuup(:,ind), mcudn(:,ind), &
                       crm_qc(:,ind), crm_qi(:,ind), crm_qs(:,ind), crm_qg(:,ind), crm_qr(:,ind), &
#ifdef m2005
                       crm_nc(:,ind), crm_ni(:,ind), crm_ns(:,ind), crm_ng(:,ind), crm_nr(:,ind), &
#ifdef MODAL_AERO
                       naermod(:,:,ind), vaerosol(:,:,ind), hygro(:,:,ind),     &
#endif 
#endif
#ifdef CLUBB_CRM
                       clubb_buffer(:,:,:,:,ind),                 &
                       crm_cld(:,:,:,ind),                      &
                       clubb_tk(:,:,:,ind), clubb_tkh(:,:,:,ind),          &
                       relvar(:,:,:,ind), accre_enhan(:,:,:,ind), qclvar(:,:,:,ind),  &
#endif
                       crm_tk(:,:,:,ind), crm_tkh(:,:,:,ind),              &
                       mu_crm(:,ind), md_crm(:,ind), du_crm(:,ind), eu_crm(:,ind), ed_crm(:,ind), jt_crm(ind), mx_crm(ind),    &
#ifdef ECPP
                       abnd(:,:,:,:,ind), abnd_tf(:,:,:,:,ind), massflxbnd(:,:,:,:,ind), acen(:,:,:,:,ind), acen_tf(:,:,:,:,ind),           &
                       rhcen(:,:,:,:,ind), qcloudcen(:,:,:,:,ind), qicecen(:,:,:,:,ind), qlsinkcen(:,:,:,:,ind), precrcen(:,:,:,:,ind), precsolidcen(:,:,:,:,ind),  & 
                       qlsink_bfcen(:,:,:,:,ind), qlsink_avgcen(:,:,:,:,ind), praincen(:,:,:,:,ind),     &
                       wupthresh_bnd(:,ind), wdownthresh_bnd(:,ind),   &
                       wwqui_cen(:,ind), wwqui_bnd(:,ind), wwqui_cloudy_cen(:,ind), wwqui_cloudy_bnd(:,ind),   &
#endif
                       tkez(:,ind), tkesgsz(:,ind), tkz(:,ind), flux_u(:,ind), flux_v(:,ind), flux_qt(:,ind), fluxsgs_qt(:,ind),flux_qp(:,ind), &
                       pflx(:,ind), qt_ls(:,ind), qt_trans(:,ind), qp_trans(:,ind), qp_fall(:,ind), &
                       qp_evp(:,ind), qp_src(:,ind), t_ls(:,ind), prectend(ind), precstend(ind), &
                       ocnfrac(ind), wndls(ind), tau00(ind), bflxls(ind), &
                       fluxu00(ind), fluxv00(ind), fluxt00(ind), fluxq00(ind),    &
                       taux_crm(ind), tauy_crm(ind), z0m(ind), timing_factor(ind), qtot(:,ind) )
    t2 = MPI_Wtime()
    tmin = min(tmin,t2-t1)
    tmax = max(tmax,t2-t1)
    call t_stopf('crm_main')
  enddo


  do i = 1 , ncols
    call dmdf_write(crm_tk          (:,:,:,  i),rank,fprefix_out,trim('crm_tk          '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.true. ,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(crm_tkh         (:,:,:,  i),rank,fprefix_out,trim('crm_tkh         '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(cltot           (        i),rank,fprefix_out,trim('cltot           ')                                                  ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(clhgh           (        i),rank,fprefix_out,trim('clhgh           ')                                                  ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(clmed           (        i),rank,fprefix_out,trim('clmed           ')                                                  ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(cllow           (        i),rank,fprefix_out,trim('cllow           ')                                                  ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(sltend          (:,      i),rank,fprefix_out,trim('sltend          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(u_crm           (:,:,:,  i),rank,fprefix_out,trim('u_crm           '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(v_crm           (:,:,:,  i),rank,fprefix_out,trim('v_crm           '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(w_crm           (:,:,:,  i),rank,fprefix_out,trim('w_crm           '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(t_crm           (:,:,:,  i),rank,fprefix_out,trim('t_crm           '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(micro_fields_crm(:,:,:,:,i),rank,fprefix_out,trim('micro_fields_crm'),(/'crm_nx          ','crm_ny          ','crm_nz          ','nmicro_fields_p1'/),.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(qltend          (:,      i),rank,fprefix_out,trim('qltend          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(qcltend         (:,      i),rank,fprefix_out,trim('qcltend         '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(qiltend         (:,      i),rank,fprefix_out,trim('qiltend         '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(t_rad           (:,:,:,  i),rank,fprefix_out,trim('t_rad           '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(qv_rad          (:,:,:,  i),rank,fprefix_out,trim('qv_rad          '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(qc_rad          (:,:,:,  i),rank,fprefix_out,trim('qc_rad          '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(qi_rad          (:,:,:,  i),rank,fprefix_out,trim('qi_rad          '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(cld_rad         (:,:,:,  i),rank,fprefix_out,trim('cld_rad         '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(cld3d_crm       (:,:,:,  i),rank,fprefix_out,trim('cld3d_crm       '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
#ifdef CRM3D
    call dmdf_write(ultend          (:,      i),rank,fprefix_out,trim('ultend          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(vltend          (:,      i),rank,fprefix_out,trim('vltend          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
#endif
#ifdef m2005
    call dmdf_write(nc_rad          (:,:,:,  i),rank,fprefix_out,trim('nc_rad          '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(ni_rad          (:,:,:,  i),rank,fprefix_out,trim('ni_rad          '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(qs_rad          (:,:,:,  i),rank,fprefix_out,trim('qs_rad          '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(ns_rad          (:,:,:,  i),rank,fprefix_out,trim('ns_rad          '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(wvar_crm        (:,:,:,  i),rank,fprefix_out,trim('wvar_crm        '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(aut_crm         (:,:,:,  i),rank,fprefix_out,trim('aut_crm         '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(acc_crm         (:,:,:,  i),rank,fprefix_out,trim('acc_crm         '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(evpc_crm        (:,:,:,  i),rank,fprefix_out,trim('evpc_crm        '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(evpr_crm        (:,:,:,  i),rank,fprefix_out,trim('evpr_crm        '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(mlt_crm         (:,:,:,  i),rank,fprefix_out,trim('mlt_crm         '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(sub_crm         (:,:,:,  i),rank,fprefix_out,trim('sub_crm         '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(dep_crm         (:,:,:,  i),rank,fprefix_out,trim('dep_crm         '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(con_crm         (:,:,:,  i),rank,fprefix_out,trim('con_crm         '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(aut_crm_a       (:,      i),rank,fprefix_out,trim('aut_crm_a       '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(acc_crm_a       (:,      i),rank,fprefix_out,trim('acc_crm_a       '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(evpc_crm_a      (:,      i),rank,fprefix_out,trim('evpc_crm_a      '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(evpr_crm_a      (:,      i),rank,fprefix_out,trim('evpr_crm_a      '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(mlt_crm_a       (:,      i),rank,fprefix_out,trim('mlt_crm_a       '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(sub_crm_a       (:,      i),rank,fprefix_out,trim('sub_crm_a       '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(dep_crm_a       (:,      i),rank,fprefix_out,trim('dep_crm_a       '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(con_crm_a       (:,      i),rank,fprefix_out,trim('con_crm_a       '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(crm_nc          (:,      i),rank,fprefix_out,trim('crm_nc          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(crm_ni          (:,      i),rank,fprefix_out,trim('crm_ni          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(crm_ns          (:,      i),rank,fprefix_out,trim('crm_ns          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(crm_ng          (:,      i),rank,fprefix_out,trim('crm_ng          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(crm_nr          (:,      i),rank,fprefix_out,trim('crm_nr          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
#endif
    call dmdf_write(precc           (        i),rank,fprefix_out,trim('precc           ')                                                  ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(precl           (        i),rank,fprefix_out,trim('precl           ')                                                  ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(cld             (:,      i),rank,fprefix_out,trim('cld             '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(cldtop          (:,      i),rank,fprefix_out,trim('cldtop          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(gicewp          (:,      i),rank,fprefix_out,trim('gicewp          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(gliqwp          (:,      i),rank,fprefix_out,trim('gliqwp          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(mc              (:,      i),rank,fprefix_out,trim('mc              '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(mcup            (:,      i),rank,fprefix_out,trim('mcup            '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(mcdn            (:,      i),rank,fprefix_out,trim('mcdn            '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(mcuup           (:,      i),rank,fprefix_out,trim('mcuup           '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(mcudn           (:,      i),rank,fprefix_out,trim('mcudn           '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(crm_qc          (:,      i),rank,fprefix_out,trim('crm_qc          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(crm_qi          (:,      i),rank,fprefix_out,trim('crm_qi          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(crm_qs          (:,      i),rank,fprefix_out,trim('crm_qs          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(crm_qg          (:,      i),rank,fprefix_out,trim('crm_qg          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(crm_qr          (:,      i),rank,fprefix_out,trim('crm_qr          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(mu_crm          (:,      i),rank,fprefix_out,trim('mu_crm          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(md_crm          (:,      i),rank,fprefix_out,trim('md_crm          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(du_crm          (:,      i),rank,fprefix_out,trim('du_crm          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(eu_crm          (:,      i),rank,fprefix_out,trim('eu_crm          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(ed_crm          (:,      i),rank,fprefix_out,trim('ed_crm          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(dd_crm          (:,      i),rank,fprefix_out,trim('dd_crm          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(jt_crm          (        i),rank,fprefix_out,trim('jt_crm          ')                                                  ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(mx_crm          (        i),rank,fprefix_out,trim('mx_crm          ')                                                  ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(mui_crm         (:,      i),rank,fprefix_out,trim('mui_crm         '),(/'plev_p1'/)                                    ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(mdi_crm         (:,      i),rank,fprefix_out,trim('mdi_crm         '),(/'plev_p1'/)                                    ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(flux_qt         (:,      i),rank,fprefix_out,trim('flux_qt         '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(fluxsgs_qt      (:,      i),rank,fprefix_out,trim('fluxsgs_qt      '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(tkez            (:,      i),rank,fprefix_out,trim('tkez            '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(tkesgsz         (:,      i),rank,fprefix_out,trim('tkesgsz         '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(tkz             (:,      i),rank,fprefix_out,trim('tkz             '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(flux_u          (:,      i),rank,fprefix_out,trim('flux_u          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(flux_v          (:,      i),rank,fprefix_out,trim('flux_v          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(flux_qp         (:,      i),rank,fprefix_out,trim('flux_qp         '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(pflx            (:,      i),rank,fprefix_out,trim('pflx            '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(qt_ls           (:,      i),rank,fprefix_out,trim('qt_ls           '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(qt_trans        (:,      i),rank,fprefix_out,trim('qt_trans        '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(qp_trans        (:,      i),rank,fprefix_out,trim('qp_trans        '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(qp_fall         (:,      i),rank,fprefix_out,trim('qp_fall         '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(qp_src          (:,      i),rank,fprefix_out,trim('qp_src          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(qp_evp          (:,      i),rank,fprefix_out,trim('qp_evp          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(t_ls            (:,      i),rank,fprefix_out,trim('t_ls            '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(prectend        (        i),rank,fprefix_out,trim('prectend        ')                                                  ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(precstend       (        i),rank,fprefix_out,trim('precstend       ')                                                  ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(precsc          (        i),rank,fprefix_out,trim('precsc          ')                                                  ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(precsl          (        i),rank,fprefix_out,trim('precsl          ')                                                  ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(taux_crm        (        i),rank,fprefix_out,trim('taux_crm        ')                                                  ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(tauy_crm        (        i),rank,fprefix_out,trim('tauy_crm        ')                                                  ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(z0m             (        i),rank,fprefix_out,trim('z0m             ')                                                  ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(timing_factor   (        i),rank,fprefix_out,trim('timing_factor   ')                                                  ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(qc_crm          (:,:,:,  i),rank,fprefix_out,trim('qc_crm          '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(qi_crm          (:,:,:,  i),rank,fprefix_out,trim('qi_crm          '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(qpc_crm         (:,:,:,  i),rank,fprefix_out,trim('qpc_crm         '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(qpi_crm         (:,:,:,  i),rank,fprefix_out,trim('qpi_crm         '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(prec_crm        (:,:,    i),rank,fprefix_out,trim('prec_crm        '),(/'crm_nx','crm_ny'/)                            ,.false.,.false.); _ERR(success,error_string,__LINE__)
    call dmdf_write(qtot            (:,      i),rank,fprefix_out,trim('qtot            '),(/'d20'/)                                        ,.false.,.true. ); _ERR(success,error_string,__LINE__)
  enddo

  call mpi_reduce((/tmin/),tred(1),1,MPI_REAL8,MPI_MIN,0,MPI_COMM_WORLD,ierr)
  call mpi_reduce((/tmax/),tred(2),1,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  
  if (rank == 0) write(*,*) 'Time Min , Max , Ratio: ',tred(1),' , ',tred(2),' , ',tred(2)/tred(1)

  call t_prf('timing', MPI_COMM_WORLD)
  call t_finalizef()
  call MPI_Finalize( ierr )

end program crm_standalone

