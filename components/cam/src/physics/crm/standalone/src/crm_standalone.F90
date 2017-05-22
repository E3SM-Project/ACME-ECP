
#ifndef CRM_STANDALONE
#define CRM_STANDALONE
#endif

!MRN: Convenient error checking macro for dumping data to file
#define _ERR(s,e,l) if (.not. s) then; write(*,*) l,', ',trim(e); stop; endif

program crm_standalone
  use mpi
  use omp_lib
  use dmdf
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
  integer :: col1, col2, ncols, stat, i, ind, ierr
  integer, parameter :: plev = PLEV
  integer :: nsamp_tot
  integer, parameter :: nsamp = NUM_SAMPLES
  integer :: rank, nranks, ithr, nthr, globid
  real(r8) :: nper

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! CRM Variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer  :: lchnk   (nsamp) 
  integer  :: icol    (nsamp) 
#ifdef CRM_STANDALONE
  real     :: latitude0(nsamp) 
  real     :: longitude0(nsamp) 
#endif
  real(r8) :: ps(nsamp) 
  real(r8) :: pmid(plev,nsamp) 
  real(r8) :: pdel(plev,nsamp) 
  real(r8) :: phis(nsamp) 
  real(r8) :: zmid(plev,nsamp) 
  real(r8) :: zint(plev+1,nsamp)
  real(r8) :: qrad_crm(crm_nx, crm_ny, crm_nz,nsamp) 
  real(r8) :: dt_gl(nsamp) 
  real(r8) :: ocnfrac(nsamp) 
  real(r8) :: tau00 (nsamp) 
  real(r8) :: wndls (nsamp) 
  real(r8) :: bflxls (nsamp) 
  real(r8) :: fluxu00 (nsamp) 
  real(r8) :: fluxv00 (nsamp) 
  real(r8) :: fluxt00 (nsamp) 
  real(r8) :: fluxq00 (nsamp) 
  real(r8) :: tl(plev,nsamp) 
  real(r8) :: ql(plev,nsamp) 
  real(r8) :: qccl(plev,nsamp)
  real(r8) :: qiil(plev,nsamp)
  real(r8) :: ul(plev,nsamp) 
  real(r8) :: vl(plev,nsamp) 
#ifdef CLUBB_CRM
  real(r8), target :: clubb_buffer(crm_nx, crm_ny, crm_nz+1,1:nclubbvars,nsamp)
  real(r8) :: crm_cld(crm_nx, crm_ny, crm_nz+1,nsamp)
  real(r8) :: clubb_tk(crm_nx, crm_ny, crm_nz,nsamp)
  real(r8) :: clubb_tkh(crm_nx, crm_ny, crm_nz,nsamp)
  real(r8) :: relvar(crm_nx, crm_ny, crm_nz,nsamp) 
  real(r8) :: accre_enhan(crm_nx, crm_ny, crm_nz,nsamp)
  real(r8) :: qclvar(crm_nx, crm_ny, crm_nz,nsamp)
#endif
  real(r8) :: crm_tk(crm_nx, crm_ny, crm_nz,nsamp)
  real(r8) :: crm_tkh(crm_nx, crm_ny, crm_nz,nsamp)
  real(r8) :: cltot(nsamp) 
  real(r8) :: clhgh(nsamp) 
  real(r8) :: clmed(nsamp) 
  real(r8) :: cllow(nsamp) 
#ifdef CRM3D
  real(r8) :: ultend(plev,nsamp) 
  real(r8) :: vltend(plev,nsamp) 
#endif
  real(r8) :: sltend(plev,nsamp) 
  real(r8) :: u_crm  (crm_nx,crm_ny,crm_nz,nsamp) 
  real(r8) :: v_crm  (crm_nx,crm_ny,crm_nz,nsamp) 
  real(r8) :: w_crm  (crm_nx,crm_ny,crm_nz,nsamp) 
  real(r8) :: t_crm  (crm_nx,crm_ny,crm_nz,nsamp) 
  real(r8), allocatable :: micro_fields_crm  (:,:,:,:,:)
  real(r8) :: qltend(plev,nsamp) 
  real(r8) :: qcltend(plev,nsamp)
  real(r8) :: qiltend(plev,nsamp)
  real(r8) :: t_rad (crm_nx, crm_ny, crm_nz,nsamp) 
  real(r8) :: qv_rad(crm_nx, crm_ny, crm_nz,nsamp) 
  real(r8) :: qc_rad(crm_nx, crm_ny, crm_nz,nsamp) 
  real(r8) :: qi_rad(crm_nx, crm_ny, crm_nz,nsamp) 
  real(r8) :: cld_rad(crm_nx, crm_ny, crm_nz,nsamp) 
  real(r8) :: cld3d_crm(crm_nx, crm_ny, crm_nz,nsamp) 
#ifdef m2005
  real(r8) :: nc_rad(crm_nx, crm_ny, crm_nz,nsamp) 
  real(r8) :: ni_rad(crm_nx, crm_ny, crm_nz,nsamp) 
  real(r8) :: qs_rad(crm_nx, crm_ny, crm_nz,nsamp) 
  real(r8) :: ns_rad(crm_nx, crm_ny, crm_nz,nsamp) 
  real(r8) :: wvar_crm(crm_nx, crm_ny, crm_nz,nsamp) 
  real(r8) :: aut_crm(crm_nx, crm_ny, crm_nz,nsamp) 
  real(r8) :: acc_crm(crm_nx, crm_ny, crm_nz,nsamp) 
  real(r8) :: evpc_crm(crm_nx, crm_ny, crm_nz,nsamp) 
  real(r8) :: evpr_crm(crm_nx, crm_ny, crm_nz,nsamp) 
  real(r8) :: mlt_crm(crm_nx, crm_ny, crm_nz,nsamp) 
  real(r8) :: sub_crm(crm_nx, crm_ny, crm_nz,nsamp) 
  real(r8) :: dep_crm(crm_nx, crm_ny, crm_nz,nsamp) 
  real(r8) :: con_crm(crm_nx, crm_ny, crm_nz,nsamp) 
  real(r8) :: aut_crm_a(plev,nsamp) 
  real(r8) :: acc_crm_a(plev,nsamp) 
  real(r8) :: evpc_crm_a(plev,nsamp) 
  real(r8) :: evpr_crm_a(plev,nsamp) 
  real(r8) :: mlt_crm_a(plev,nsamp) 
  real(r8) :: sub_crm_a(plev,nsamp) 
  real(r8) :: dep_crm_a(plev,nsamp) 
  real(r8) :: con_crm_a(plev,nsamp) 
#endif
  real(r8) :: precc(nsamp) 
  real(r8) :: precl(nsamp) 
  real(r8) :: cld(plev,nsamp)  
  real(r8) :: cldtop(plev,nsamp)  
  real(r8) :: gicewp(plev,nsamp)  
  real(r8) :: gliqwp(plev,nsamp)  
  real(r8) :: mc(plev,nsamp)   
  real(r8) :: mcup(plev,nsamp) 
  real(r8) :: mcdn(plev,nsamp) 
  real(r8) :: mcuup(plev,nsamp) 
  real(r8) :: mcudn(plev,nsamp) 
  real(r8) :: crm_qc(plev,nsamp)  
  real(r8) :: crm_qi(plev,nsamp)  
  real(r8) :: crm_qs(plev,nsamp)  
  real(r8) :: crm_qg(plev,nsamp)  
  real(r8) :: crm_qr(plev,nsamp)  
#ifdef m2005
  real(r8) :: crm_nc(plev,nsamp)  
  real(r8) :: crm_ni(plev,nsamp)  
  real(r8) :: crm_ns(plev,nsamp)  
  real(r8) :: crm_ng(plev,nsamp)  
  real(r8) :: crm_nr(plev,nsamp)  
#ifdef MODAL_AERO
  real(r8) :: naermod(plev, ntot_amode,nsamp)     
  real(r8) :: vaerosol(plev, ntot_amode,nsamp)    
  real(r8) :: hygro(plev, ntot_amode,nsamp)       
#endif 
#endif
  real(r8) :: mu_crm (plev,nsamp)             
  real(r8) :: md_crm (plev,nsamp)             
  real(r8) :: du_crm (plev,nsamp)             
  real(r8) :: eu_crm (plev,nsamp)             
  real(r8) :: ed_crm (plev,nsamp)             
  real(r8) :: dd_crm (plev,nsamp)             
  real(r8) :: jt_crm (nsamp)                   
  real(r8) :: mx_crm (nsamp)                   
  real(r8) :: mui_crm (plev+1,nsamp)             
  real(r8) :: mdi_crm (plev+1,nsamp)             
  real(r8) :: flux_qt(plev,nsamp) 
  real(r8) :: fluxsgs_qt(plev,nsamp) 
  real(r8) :: tkez(plev,nsamp) 
  real(r8) :: tkesgsz(plev,nsamp) 
  real(r8) :: tkz(plev,nsamp)  
  real(r8) :: flux_u(plev,nsamp) 
  real(r8) :: flux_v(plev,nsamp) 
  real(r8) :: flux_qp(plev,nsamp) 
  real(r8) :: pflx(plev,nsamp)    
  real(r8) :: qt_ls(plev,nsamp) 
  real(r8) :: qt_trans(plev,nsamp)
  real(r8) :: qp_trans(plev,nsamp) 
  real(r8) :: qp_fall(plev,nsamp) 
  real(r8) :: qp_src(plev,nsamp) 
  real(r8) :: qp_evp(plev,nsamp) 
  real(r8) :: t_ls(plev,nsamp) 
  real(r8) :: prectend(nsamp) 
  real(r8) :: precstend(nsamp) 
  real(r8) :: precsc(nsamp) 
  real(r8) :: precsl(nsamp) 
  real(r8) :: taux_crm (nsamp) 
  real(r8) :: tauy_crm (nsamp) 
  real(r8) :: z0m(nsamp) 
  real(r8) :: timing_factor(nsamp) 
  real(r8) :: qc_crm (crm_nx, crm_ny, crm_nz,nsamp)
  real(r8) :: qi_crm (crm_nx, crm_ny, crm_nz,nsamp)
  real(r8) :: qpc_crm(crm_nx, crm_ny, crm_nz,nsamp)
  real(r8) :: qpi_crm(crm_nx, crm_ny, crm_nz,nsamp)
  real(r8) :: prec_crm(crm_nx, crm_ny,nsamp)
#ifdef ECPP
  real(r8) :: acen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,nsamp)   
  real(r8) :: acen_tf(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,nsamp) 
  real(r8) :: rhcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,nsamp)  
  real(r8) :: qcloudcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,nsamp)  
  real(r8) :: qicecen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,nsamp) 
  real(r8) :: qlsinkcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,nsamp)  
  real(r8) :: precrcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,nsamp)   
  real(r8) :: precsolidcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,nsamp)   
  real(r8) :: qlsink_bfcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,nsamp)  
  real(r8) :: qlsink_avgcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,nsamp)  
  real(r8) :: praincen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,nsamp)  
  real(r8) :: wwqui_cen(plev,nsamp)                                
  real(r8) :: wwqui_cloudy_cen(plev,nsamp)                         
  real(r8) :: abnd(plev+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,nsamp)   
  real(r8) :: abnd_tf(plev+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,nsamp) 
  real(r8) :: massflxbnd(plev+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR,nsamp) 
  real(r8) :: wupthresh_bnd(plev+1,nsamp)             
  real(r8) :: wdownthresh_bnd(plev+1,nsamp)           
  real(r8) :: wwqui_bnd(plev+1,nsamp)                                
  real(r8) :: wwqui_cloudy_bnd(plev+1,nsamp)                         
#endif
  real(r8) :: qtot(20,nsamp)

  real(r8) :: u_crm_f  (crm_nx,crm_ny,crm_nz,nsamp) 
  real(r8) :: v_crm_f  (crm_nx,crm_ny,crm_nz,nsamp) 
  real(r8) :: w_crm_f  (crm_nx,crm_ny,crm_nz,nsamp) 
  real(r8) :: t_crm_f  (crm_nx,crm_ny,crm_nz,nsamp) 

  !Get MPI and OpenMP information
  call MPI_Init( ierr )
  call MPI_Comm_size( MPI_COMM_WORLD , nranks , ierr )
  call MPI_Comm_rank( MPI_COMM_WORLD , rank   , ierr)
  nthr = omp_get_num_threads()
  ithr = omp_get_thread_num()
  if (rank == 0) write(unit=*,fmt='(4(A,I5),A)') 'Taks (Thread): ',rank,' (',ithr,') / ',nranks,' (',nthr,')'


  !I have to call setparm to get the correct value for nmicro_fields
  call setparm()
  allocate( micro_fields_crm  (crm_nx,crm_ny,crm_nz,nmicro_fields+1,nsamp) )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Read in the command line arguments
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call getarg(1,fname_in)
  call getarg(2,fname_out)
  call dmdf_num_records(fname_in,nsamp_tot)
  if (rank == 0) write(unit=*,fmt='(2(A,I7),A)') 'Processing: ',nsamp,' of ',nsamp_tot,' samples.'

  globid = rank*nthr + ithr
  nper = real(nsamp)/(nranks*nthr)
  col1 = nint( nper* globid    )+1
  col2 = nint( nper*(globid+1) )
  ncols = col2 - col1 + 1
  call MPI_Barrier( MPI_COMM_WORLD , ierr )
  write(*,*) col1,col2
  call MPI_Barrier( MPI_COMM_WORLD , ierr )

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
  if (rank == 0) write(*,*) '*** FINISHED READING INPUT DATA ***'

  call MPI_Barrier( MPI_COMM_WORLD , ierr )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Call the crm routine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ind = 1 , ncols
    write(*,*) '*** CALLING CRM ',ind,' ***'
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
  enddo

  call MPI_Barrier( MPI_COMM_WORLD , ierr )

  if (rank == 0) write(*,*) '*** READING OUTPUT DATA ***'
  call dmdf_read(u_crm_f,trim(fname_out),'u_crm',col1,col2,.true. ,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(v_crm_f,trim(fname_out),'v_crm',col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(w_crm_f,trim(fname_out),'w_crm',col1,col2,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(t_crm_f,trim(fname_out),'t_crm',col1,col2,.false.,.true. ); _ERR(success,error_string,__LINE__)
  if (rank == 0) write(*,*) '*** FINISHED READING OUTPUT DATA ***'

  if (rank == 0) write(*,*) '*** CHECKING ANSWERS ***'
  call MPI_Barrier( MPI_COMM_WORLD , ierr )
  do i = 0 , nranks
    if (rank == i) then
      write(*,*) 'Rank: ',i
      write(*,*) sum(abs(u_crm(:,:,:,1:ncols)-u_crm_f(:,:,:,1:ncols))) / sum(abs(u_crm_f(:,:,:,1:ncols)))
      write(*,*) sum(abs(v_crm(:,:,:,1:ncols)-v_crm_f(:,:,:,1:ncols))) / sum(abs(v_crm_f(:,:,:,1:ncols)))
      write(*,*) sum(abs(w_crm(:,:,:,1:ncols)-w_crm_f(:,:,:,1:ncols))) / sum(abs(w_crm_f(:,:,:,1:ncols)))
      write(*,*) sum(abs(t_crm(:,:,:,1:ncols)-t_crm_f(:,:,:,1:ncols))) / sum(abs(t_crm_f(:,:,:,1:ncols)))
      write(*,*)
    endif
    call MPI_Barrier( MPI_COMM_WORLD , ierr )
  enddo

  call MPI_Finalize( ierr )




end program crm_standalone

