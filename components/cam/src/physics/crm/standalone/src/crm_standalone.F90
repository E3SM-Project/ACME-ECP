
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
  character(len=512) :: fname_in, fname_out
  integer :: col1, col2, ncols, stat, i, ind, ierr, nsamp_tot
  integer :: rank, nranks, globid
  integer, parameter :: plev = PLEV
  real(r8) :: nper, t1, t2, tmin, tmax

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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! COPIES OF CRM VARIABLES FOR COMPARISONS
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer , allocatable :: lchnk_f(:) 
  integer , allocatable :: icol_f(:) 
#ifdef CRM_STANDALONE
  real    , allocatable :: latitude0_f(:) 
  real    , allocatable :: longitude0_f(:) 
#endif
  real(r8), allocatable :: ps_f(:) 
  real(r8), allocatable :: pmid_f(:,:) 
  real(r8), allocatable :: pdel_f(:,:) 
  real(r8), allocatable :: phis_f(:) 
  real(r8), allocatable :: zmid_f(:,:) 
  real(r8), allocatable :: zint_f(:,:)
  real(r8), allocatable :: qrad_crm_f(:, :, :,:) 
  real(r8), allocatable :: dt_gl_f(:) 
  real(r8), allocatable :: ocnfrac_f(:) 
  real(r8), allocatable :: tau00_f(:) 
  real(r8), allocatable :: wndls_f(:) 
  real(r8), allocatable :: bflxls_f(:) 
  real(r8), allocatable :: fluxu00_f(:) 
  real(r8), allocatable :: fluxv00_f(:) 
  real(r8), allocatable :: fluxt00_f(:) 
  real(r8), allocatable :: fluxq00_f(:) 
  real(r8), allocatable :: tl_f(:,:) 
  real(r8), allocatable :: ql_f(:,:) 
  real(r8), allocatable :: qccl_f(:,:)
  real(r8), allocatable :: qiil_f(:,:)
  real(r8), allocatable :: ul_f(:,:) 
  real(r8), allocatable :: vl_f(:,:) 
#ifdef CLUBB_CRM
  real(r8), allocatable, target :: clubb_buffer_f(:, :, :,:,:)
  real(r8), allocatable :: crm_cld_f(:, :, :,:)
  real(r8), allocatable :: clubb_tk_f(:, :, :,:)
  real(r8), allocatable :: clubb_tkh_f(:, :, :,:)
  real(r8), allocatable :: relvar_f(:, :, :,:) 
  real(r8), allocatable :: accre_enhan_f(:, :, :,:)
  real(r8), allocatable :: qclvar_f(:, :, :,:)
#endif
  real(r8), allocatable :: crm_tk_f(:, :, :,:)
  real(r8), allocatable :: crm_tkh_f(:, :, :,:)
  real(r8), allocatable :: cltot_f(:) 
  real(r8), allocatable :: clhgh_f(:) 
  real(r8), allocatable :: clmed_f(:) 
  real(r8), allocatable :: cllow_f(:) 
#ifdef CRM3D
  real(r8), allocatable :: ultend_f(:,:) 
  real(r8), allocatable :: vltend_f(:,:) 
#endif
  real(r8), allocatable :: sltend_f(:,:) 
  real(r8), allocatable :: u_crm_f(:,:,:,:) 
  real(r8), allocatable :: v_crm_f(:,:,:,:) 
  real(r8), allocatable :: w_crm_f(:,:,:,:) 
  real(r8), allocatable :: t_crm_f(:,:,:,:) 
  real(r8), allocatable :: micro_fields_crm_f(:,:,:,:,:)
  real(r8), allocatable :: qltend_f(:,:) 
  real(r8), allocatable :: qcltend_f(:,:)
  real(r8), allocatable :: qiltend_f(:,:)
  real(r8), allocatable :: t_rad_f(:, :, :,:) 
  real(r8), allocatable :: qv_rad_f(:, :, :,:) 
  real(r8), allocatable :: qc_rad_f(:, :, :,:) 
  real(r8), allocatable :: qi_rad_f(:, :, :,:) 
  real(r8), allocatable :: cld_rad_f(:, :, :,:) 
  real(r8), allocatable :: cld3d_crm_f(:, :, :,:) 
#ifdef m2005
  real(r8), allocatable :: nc_rad_f(:, :, :,:) 
  real(r8), allocatable :: ni_rad_f(:, :, :,:) 
  real(r8), allocatable :: qs_rad_f(:, :, :,:) 
  real(r8), allocatable :: ns_rad_f(:, :, :,:) 
  real(r8), allocatable :: wvar_crm_f(:, :, :,:) 
  real(r8), allocatable :: aut_crm_f(:, :, :,:) 
  real(r8), allocatable :: acc_crm_f(:, :, :,:) 
  real(r8), allocatable :: evpc_crm_f(:, :, :,:) 
  real(r8), allocatable :: evpr_crm_f(:, :, :,:) 
  real(r8), allocatable :: mlt_crm_f(:, :, :,:) 
  real(r8), allocatable :: sub_crm_f(:, :, :,:) 
  real(r8), allocatable :: dep_crm_f(:, :, :,:) 
  real(r8), allocatable :: con_crm_f(:, :, :,:) 
  real(r8), allocatable :: aut_crm_a_f(:,:) 
  real(r8), allocatable :: acc_crm_a_f(:,:) 
  real(r8), allocatable :: evpc_crm_a_f(:,:) 
  real(r8), allocatable :: evpr_crm_a_f(:,:) 
  real(r8), allocatable :: mlt_crm_a_f(:,:) 
  real(r8), allocatable :: sub_crm_a_f(:,:) 
  real(r8), allocatable :: dep_crm_a_f(:,:) 
  real(r8), allocatable :: con_crm_a_f(:,:) 
#endif
  real(r8), allocatable :: precc_f(:) 
  real(r8), allocatable :: precl_f(:) 
  real(r8), allocatable :: cld_f(:,:)  
  real(r8), allocatable :: cldtop_f(:,:)  
  real(r8), allocatable :: gicewp_f(:,:)  
  real(r8), allocatable :: gliqwp_f(:,:)  
  real(r8), allocatable :: mc_f(:,:)   
  real(r8), allocatable :: mcup_f(:,:) 
  real(r8), allocatable :: mcdn_f(:,:) 
  real(r8), allocatable :: mcuup_f(:,:) 
  real(r8), allocatable :: mcudn_f(:,:) 
  real(r8), allocatable :: crm_qc_f(:,:)  
  real(r8), allocatable :: crm_qi_f(:,:)  
  real(r8), allocatable :: crm_qs_f(:,:)  
  real(r8), allocatable :: crm_qg_f(:,:)  
  real(r8), allocatable :: crm_qr_f(:,:)  
#ifdef m2005
  real(r8), allocatable :: crm_nc_f(:,:)  
  real(r8), allocatable :: crm_ni_f(:,:)  
  real(r8), allocatable :: crm_ns_f(:,:)  
  real(r8), allocatable :: crm_ng_f(:,:)  
  real(r8), allocatable :: crm_nr_f(:,:)  
#ifdef MODAL_AERO
  real(r8), allocatable :: naermod_f(:, :,:)     
  real(r8), allocatable :: vaerosol_f(:, :,:)    
  real(r8), allocatable :: hygro_f(:, :,:)       
#endif 
#endif
  real(r8), allocatable :: mu_crm_f(:,:)             
  real(r8), allocatable :: md_crm_f(:,:)             
  real(r8), allocatable :: du_crm_f(:,:)             
  real(r8), allocatable :: eu_crm_f(:,:)             
  real(r8), allocatable :: ed_crm_f(:,:)             
  real(r8), allocatable :: dd_crm_f(:,:)             
  real(r8), allocatable :: jt_crm_f(:)                   
  real(r8), allocatable :: mx_crm_f(:)                   
  real(r8), allocatable :: mui_crm_f(:,:)             
  real(r8), allocatable :: mdi_crm_f(:,:)             
  real(r8), allocatable :: flux_qt_f(:,:) 
  real(r8), allocatable :: fluxsgs_qt_f(:,:) 
  real(r8), allocatable :: tkez_f(:,:) 
  real(r8), allocatable :: tkesgsz_f(:,:) 
  real(r8), allocatable :: tkz_f(:,:)  
  real(r8), allocatable :: flux_u_f(:,:) 
  real(r8), allocatable :: flux_v_f(:,:) 
  real(r8), allocatable :: flux_qp_f(:,:) 
  real(r8), allocatable :: pflx_f(:,:)    
  real(r8), allocatable :: qt_ls_f(:,:) 
  real(r8), allocatable :: qt_trans_f(:,:)
  real(r8), allocatable :: qp_trans_f(:,:) 
  real(r8), allocatable :: qp_fall_f(:,:) 
  real(r8), allocatable :: qp_src_f(:,:) 
  real(r8), allocatable :: qp_evp_f(:,:) 
  real(r8), allocatable :: t_ls_f(:,:) 
  real(r8), allocatable :: prectend_f(:) 
  real(r8), allocatable :: precstend_f(:) 
  real(r8), allocatable :: precsc_f(:) 
  real(r8), allocatable :: precsl_f(:) 
  real(r8), allocatable :: taux_crm_f(:) 
  real(r8), allocatable :: tauy_crm_f(:) 
  real(r8), allocatable :: z0m_f(:) 
  real(r8), allocatable :: timing_factor_f(:) 
  real(r8), allocatable :: qc_crm_f(:, :, :,:)
  real(r8), allocatable :: qi_crm_f(:, :, :,:)
  real(r8), allocatable :: qpc_crm_f(:, :, :,:)
  real(r8), allocatable :: qpi_crm_f(:, :, :,:)
  real(r8), allocatable :: prec_crm_f(:, :,:)
#ifdef ECPP
  real(r8), allocatable :: acen_f(:,:,:,:,:)   
  real(r8), allocatable :: acen_tf_f(:,:,:,:,:) 
  real(r8), allocatable :: rhcen_f(:,:,:,:,:)  
  real(r8), allocatable :: qcloudcen_f(:,:,:,:,:)  
  real(r8), allocatable :: qicecen_f(:,:,:,:,:) 
  real(r8), allocatable :: qlsinkcen_f(:,:,:,:,:)  
  real(r8), allocatable :: precrcen_f(:,:,:,:,:)   
  real(r8), allocatable :: precsolidcen_f(:,:,:,:,:)   
  real(r8), allocatable :: qlsink_bfcen_f(:,:,:,:,:)  
  real(r8), allocatable :: qlsink_avgcen_f(:,:,:,:,:)  
  real(r8), allocatable :: praincen_f(:,:,:,:,:)  
  real(r8), allocatable :: wwqui_cen_f(:,:)                                
  real(r8), allocatable :: wwqui_cloudy_cen_f(:,:)                         
  real(r8), allocatable :: abnd_f(:,:,:,:,:)   
  real(r8), allocatable :: abnd_tf_f(:,:,:,:,:) 
  real(r8), allocatable :: massflxbnd_f(:,:,:,:,:) 
  real(r8), allocatable :: wupthresh_bnd_f(:,:)             
  real(r8), allocatable :: wdownthresh_bnd_f(:,:)           
  real(r8), allocatable :: wwqui_bnd_f(:,:)                                
  real(r8), allocatable :: wwqui_cloudy_bnd_f(:,:)                         
#endif
  real(r8), allocatable :: qtot_f(:,:)

  !Get MPI information
  call MPI_Init( ierr )
  call MPI_Comm_size( MPI_COMM_WORLD , nranks , ierr )
  call MPI_Comm_rank( MPI_COMM_WORLD , rank   , ierr )

  call t_initf('gptl.nl',LogPrint=(rank==0),Mpicom=MPI_COMM_WORLD,MasterTask=(rank==0))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Read in the command line arguments
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call getarg(1,fname_in )
  call getarg(2,fname_out)
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
  write(*,*) globid,nper,col1,col2,ncols

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
  !! Allocate copies of CRM variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate( lchnk_f(ncols) ) 
  allocate( icol_f(ncols) ) 
#ifdef CRM_STANDALONE
  allocate( latitude0_f(ncols) ) 
  allocate( longitude0_f(ncols) ) 
#endif
  allocate( ps_f(ncols) ) 
  allocate( pmid_f(plev,ncols) ) 
  allocate( pdel_f(plev,ncols) ) 
  allocate( phis_f(ncols) ) 
  allocate( zmid_f(plev,ncols) ) 
  allocate( zint_f(plev+1,ncols) )
  allocate( qrad_crm_f(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( dt_gl_f(ncols) ) 
  allocate( ocnfrac_f(ncols) ) 
  allocate( tau00_f(ncols) ) 
  allocate( wndls_f(ncols) ) 
  allocate( bflxls_f(ncols) ) 
  allocate( fluxu00_f(ncols) ) 
  allocate( fluxv00_f(ncols) ) 
  allocate( fluxt00_f(ncols) ) 
  allocate( fluxq00_f(ncols) ) 
  allocate( tl_f(plev,ncols) ) 
  allocate( ql_f(plev,ncols) ) 
  allocate( qccl_f(plev,ncols) )
  allocate( qiil_f(plev,ncols) )
  allocate( ul_f(plev,ncols) ) 
  allocate( vl_f(plev,ncols) ) 
#ifdef CLUBB_CRM
  allocate( clubb_buffer_f(crm_nx, crm_ny, crm_nz+1,1:nclubbvars,ncols) )
  allocate( crm_cld_f(crm_nx, crm_ny, crm_nz+1,ncols) )
  allocate( clubb_tk_f(crm_nx, crm_ny, crm_nz,ncols) )
  allocate( clubb_tkh_f(crm_nx, crm_ny, crm_nz,ncols) )
  allocate( relvar_f(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( accre_enhan_f(crm_nx, crm_ny, crm_nz,ncols) )
  allocate( qclvar_f(crm_nx, crm_ny, crm_nz,ncols) )
#endif
  allocate( crm_tk_f(crm_nx, crm_ny, crm_nz,ncols) )
  allocate( crm_tkh_f(crm_nx, crm_ny, crm_nz,ncols) )
  allocate( cltot_f(ncols) ) 
  allocate( clhgh_f(ncols) ) 
  allocate( clmed_f(ncols) ) 
  allocate( cllow_f(ncols) ) 
#ifdef CRM3D
  allocate( ultend_f(plev,ncols) ) 
  allocate( vltend_f(plev,ncols) ) 
#endif
  allocate( sltend_f(plev,ncols) ) 
  allocate( u_crm_f(crm_nx,crm_ny,crm_nz,ncols) ) 
  allocate( v_crm_f(crm_nx,crm_ny,crm_nz,ncols) ) 
  allocate( w_crm_f(crm_nx,crm_ny,crm_nz,ncols) ) 
  allocate( t_crm_f(crm_nx,crm_ny,crm_nz,ncols) ) 
  allocate( micro_fields_crm_f(crm_nx,crm_ny,crm_nz,nmicro_fields+1,ncols) )
  allocate( qltend_f(plev,ncols) ) 
  allocate( qcltend_f(plev,ncols) )
  allocate( qiltend_f(plev,ncols) )
  allocate( t_rad_f(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( qv_rad_f(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( qc_rad_f(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( qi_rad_f(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( cld_rad_f(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( cld3d_crm_f(crm_nx, crm_ny, crm_nz,ncols) ) 
#ifdef m2005
  allocate( nc_rad_f(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( ni_rad_f(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( qs_rad_f(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( ns_rad_f(crm_nx, crm_ny, crm_nz,ncols) )
  allocate( wvar_crm_f(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( aut_crm_f(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( acc_crm_f(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( evpc_crm_f(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( evpr_crm_f(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( mlt_crm_f(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( sub_crm_f(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( dep_crm_f(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( con_crm_f(crm_nx, crm_ny, crm_nz,ncols) ) 
  allocate( aut_crm_a_f(plev,ncols) ) 
  allocate( acc_crm_a_f(plev,ncols) ) 
  allocate( evpc_crm_a_f(plev,ncols) ) 
  allocate( evpr_crm_a_f(plev,ncols) ) 
  allocate( mlt_crm_a_f(plev,ncols) ) 
  allocate( sub_crm_a_f(plev,ncols) ) 
  allocate( dep_crm_a_f(plev,ncols) ) 
  allocate( con_crm_a_f(plev,ncols) ) 
#endif
  allocate( precc_f(ncols) ) 
  allocate( precl_f(ncols) ) 
  allocate( cld_f(plev,ncols) )  
  allocate( cldtop_f(plev,ncols) )  
  allocate( gicewp_f(plev,ncols) )  
  allocate( gliqwp_f(plev,ncols) )  
  allocate( mc_f(plev,ncols) )   
  allocate( mcup_f(plev,ncols) ) 
  allocate( mcdn_f(plev,ncols) ) 
  allocate( mcuup_f(plev,ncols) ) 
  allocate( mcudn_f(plev,ncols) ) 
  allocate( crm_qc_f(plev,ncols) )  
  allocate( crm_qi_f(plev,ncols) )  
  allocate( crm_qs_f(plev,ncols) )  
  allocate( crm_qg_f(plev,ncols) )  
  allocate( crm_qr_f(plev,ncols) )  
#ifdef m2005
  allocate( crm_nc_f(plev,ncols) )  
  allocate( crm_ni_f(plev,ncols) )  
  allocate( crm_ns_f(plev,ncols) )  
  allocate( crm_ng_f(plev,ncols) )  
  allocate( crm_nr_f(plev,ncols) )  
#ifdef MODAL_AERO
  allocate( naermod_f(plev, ntot_amode,ncols) )     
  allocate( vaerosol_f(plev, ntot_amode,ncols) )    
  allocate( hygro_f(plev, ntot_amode,ncols) )       
#endif 
#endif
  allocate( mu_crm_f(plev,ncols) )             
  allocate( md_crm_f(plev,ncols) )             
  allocate( du_crm_f(plev,ncols) )             
  allocate( eu_crm_f(plev,ncols) )             
  allocate( ed_crm_f(plev,ncols) )             
  allocate( dd_crm_f(plev,ncols) )             
  allocate( jt_crm_f(ncols) )                   
  allocate( mx_crm_f(ncols) )                   
  allocate( mui_crm_f(plev+1,ncols) )             
  allocate( mdi_crm_f(plev+1,ncols) )             
  allocate( flux_qt_f(plev,ncols) ) 
  allocate( fluxsgs_qt_f(plev,ncols) ) 
  allocate( tkez_f(plev,ncols) ) 
  allocate( tkesgsz_f(plev,ncols) ) 
  allocate( tkz_f(plev,ncols) )  
  allocate( flux_u_f(plev,ncols) ) 
  allocate( flux_v_f(plev,ncols) ) 
  allocate( flux_qp_f(plev,ncols) ) 
  allocate( pflx_f(plev,ncols) )    
  allocate( qt_ls_f(plev,ncols) ) 
  allocate( qt_trans_f(plev,ncols) )
  allocate( qp_trans_f(plev,ncols) ) 
  allocate( qp_fall_f(plev,ncols) ) 
  allocate( qp_src_f(plev,ncols) ) 
  allocate( qp_evp_f(plev,ncols) ) 
  allocate( t_ls_f(plev,ncols) ) 
  allocate( prectend_f(ncols) ) 
  allocate( precstend_f(ncols) ) 
  allocate( precsc_f(ncols) ) 
  allocate( precsl_f(ncols) ) 
  allocate( taux_crm_f(ncols) ) 
  allocate( tauy_crm_f(ncols) ) 
  allocate( z0m_f(ncols) ) 
  allocate( timing_factor_f(ncols) ) 
  allocate( qc_crm_f(crm_nx, crm_ny, crm_nz,ncols) )
  allocate( qi_crm_f(crm_nx, crm_ny, crm_nz,ncols) )
  allocate( qpc_crm_f(crm_nx, crm_ny, crm_nz,ncols) )
  allocate( qpi_crm_f(crm_nx, crm_ny, crm_nz,ncols) )
  allocate( prec_crm_f(crm_nx, crm_ny,ncols) )
#ifdef ECPP
  allocate( acen_f(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,ncols) )   
  allocate( acen_tf_f(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,ncols) ) 
  allocate( rhcen_f(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,ncols) )  
  allocate( qcloudcen_f(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,ncols) )  
  allocate( qicecen_f(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,ncols) ) 
  allocate( qlsinkcen_f(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,ncols) )  
  allocate( precrcen_f(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,ncols) )   
  allocate( precsolidcen_f(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,ncols) )   
  allocate( qlsink_bfcen_f(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,ncols) )  
  allocate( qlsink_avgcen_f(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,ncols) )  
  allocate( praincen_f(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,ncols) )  
  allocate( wwqui_cen_f(plev,ncols) )                                
  allocate( wwqui_cloudy_cen_f(plev,ncols) )                         
  allocate( abnd_f(plev+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,ncols) )   
  allocate( abnd_tf_f(plev+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,ncols) ) 
  allocate( massflxbnd_f(plev+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,ncols) ) 
  allocate( wupthresh_bnd_f(plev+1,ncols) )             
  allocate( wdownthresh_bnd_f(plev+1,ncols) )           
  allocate( wwqui_bnd_f(plev+1,ncols) )                                
  allocate( wwqui_cloudy_bnd_f(plev+1,ncols) )                         
#endif
  allocate( qtot_f(20,ncols) )

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
    dmdf_write(u_crm(:,:,:,i),rank,fprefix,vname,dnames,first,last)   !For array values
  enddo

  call MPI_Barrier( MPI_COMM_WORLD , ierr )
  call t_prf('timing', MPI_COMM_WORLD)
  call MPI_Barrier( MPI_COMM_WORLD , ierr )
  call t_finalizef()
  call MPI_Finalize( ierr )

end program crm_standalone

