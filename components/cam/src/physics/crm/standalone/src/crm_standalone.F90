
#ifndef CRM_STANDALONE
#define CRM_STANDALONE
#endif

!MRN: Convenient error checking macro for dumping data to file
#define _ERR(s,e,l) if (.not. s) then; write(*,*) l,', ',trim(e); stop; endif

program crm_standalone
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
  character(len=32)  :: ind1_str, ind2_str
  integer :: ind1, ind2, stat
  integer, parameter :: plev = PLEV

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! CRM Variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer  :: lchnk    
  integer  :: icol     
#ifdef CRM_STANDALONE
  real     :: latitude0
  real     :: longitude0
#endif
  real(r8) :: ps 
  real(r8) :: pmid(plev) 
  real(r8) :: pdel(plev) 
  real(r8) :: phis 
  real(r8) :: zmid(plev) 
  real(r8) :: zint(plev+1)
  real(r8) :: qrad_crm(crm_nx, crm_ny, crm_nz) 
  real(r8) :: dt_gl 
  real(r8) :: ocnfrac 
  real(r8) :: tau00  
  real(r8) :: wndls  
  real(r8) :: bflxls  
  real(r8) :: fluxu00  
  real(r8) :: fluxv00  
  real(r8) :: fluxt00  
  real(r8) :: fluxq00  
  real(r8) :: tl(plev) 
  real(r8) :: ql(plev) 
  real(r8) :: qccl(plev)
  real(r8) :: qiil(plev)
  real(r8) :: ul(plev) 
  real(r8) :: vl(plev) 
#ifdef CLUBB_CRM
  real(r8), target :: clubb_buffer(crm_nx, crm_ny, crm_nz+1,1:nclubbvars)
  real(r8) :: crm_cld(crm_nx, crm_ny, crm_nz+1)
  real(r8) :: clubb_tk(crm_nx, crm_ny, crm_nz)
  real(r8) :: clubb_tkh(crm_nx, crm_ny, crm_nz)
  real(r8) :: relvar(crm_nx, crm_ny, crm_nz) 
  real(r8) :: accre_enhan(crm_nx, crm_ny, crm_nz)
  real(r8) :: qclvar(crm_nx, crm_ny, crm_nz)
#endif
  real(r8) :: crm_tk(crm_nx, crm_ny, crm_nz)
  real(r8) :: crm_tkh(crm_nx, crm_ny, crm_nz)
  real(r8) :: cltot 
  real(r8) :: clhgh 
  real(r8) :: clmed 
  real(r8) :: cllow 
#ifdef CRM3D
  real(r8) :: ultend(plev) 
  real(r8) :: vltend(plev) 
#endif
  real(r8) :: sltend(plev) 
  real(r8) :: u_crm  (crm_nx,crm_ny,crm_nz) 
  real(r8) :: v_crm  (crm_nx,crm_ny,crm_nz) 
  real(r8) :: w_crm  (crm_nx,crm_ny,crm_nz) 
  real(r8) :: t_crm  (crm_nx,crm_ny,crm_nz) 
  real(r8), allocatable :: micro_fields_crm  (:,:,:,:)
  real(r8) :: qltend(plev) 
  real(r8) :: qcltend(plev)
  real(r8) :: qiltend(plev)
  real(r8) :: t_rad (crm_nx, crm_ny, crm_nz) 
  real(r8) :: qv_rad(crm_nx, crm_ny, crm_nz) 
  real(r8) :: qc_rad(crm_nx, crm_ny, crm_nz) 
  real(r8) :: qi_rad(crm_nx, crm_ny, crm_nz) 
  real(r8) :: cld_rad(crm_nx, crm_ny, crm_nz) 
  real(r8) :: cld3d_crm(crm_nx, crm_ny, crm_nz) 
#ifdef m2005
  real(r8) :: nc_rad(crm_nx, crm_ny, crm_nz) 
  real(r8) :: ni_rad(crm_nx, crm_ny, crm_nz) 
  real(r8) :: qs_rad(crm_nx, crm_ny, crm_nz) 
  real(r8) :: ns_rad(crm_nx, crm_ny, crm_nz) 
  real(r8) :: wvar_crm(crm_nx, crm_ny, crm_nz) 
  real(r8) :: aut_crm(crm_nx, crm_ny, crm_nz) 
  real(r8) :: acc_crm(crm_nx, crm_ny, crm_nz) 
  real(r8) :: evpc_crm(crm_nx, crm_ny, crm_nz) 
  real(r8) :: evpr_crm(crm_nx, crm_ny, crm_nz) 
  real(r8) :: mlt_crm(crm_nx, crm_ny, crm_nz) 
  real(r8) :: sub_crm(crm_nx, crm_ny, crm_nz) 
  real(r8) :: dep_crm(crm_nx, crm_ny, crm_nz) 
  real(r8) :: con_crm(crm_nx, crm_ny, crm_nz) 
  real(r8) :: aut_crm_a(plev) 
  real(r8) :: acc_crm_a(plev) 
  real(r8) :: evpc_crm_a(plev) 
  real(r8) :: evpr_crm_a(plev) 
  real(r8) :: mlt_crm_a(plev) 
  real(r8) :: sub_crm_a(plev) 
  real(r8) :: dep_crm_a(plev) 
  real(r8) :: con_crm_a(plev) 
#endif
  real(r8) :: precc 
  real(r8) :: precl 
  real(r8) :: cld(plev)  
  real(r8) :: cldtop(plev)  
  real(r8) :: gicewp(plev)  
  real(r8) :: gliqwp(plev)  
  real(r8) :: mc(plev)   
  real(r8) :: mcup(plev) 
  real(r8) :: mcdn(plev) 
  real(r8) :: mcuup(plev) 
  real(r8) :: mcudn(plev) 
  real(r8) :: crm_qc(plev)  
  real(r8) :: crm_qi(plev)  
  real(r8) :: crm_qs(plev)  
  real(r8) :: crm_qg(plev)  
  real(r8) :: crm_qr(plev)  
#ifdef m2005
  real(r8) :: crm_nc(plev)  
  real(r8) :: crm_ni(plev)  
  real(r8) :: crm_ns(plev)  
  real(r8) :: crm_ng(plev)  
  real(r8) :: crm_nr(plev)  
#ifdef MODAL_AERO
  real(r8) :: naermod(plev, ntot_amode)     
  real(r8) :: vaerosol(plev, ntot_amode)    
  real(r8) :: hygro(plev, ntot_amode)       
#endif 
#endif
  real(r8) :: mu_crm (plev)             
  real(r8) :: md_crm (plev)             
  real(r8) :: du_crm (plev)             
  real(r8) :: eu_crm (plev)             
  real(r8) :: ed_crm (plev)             
  real(r8) :: dd_crm (plev)             
  real(r8) :: jt_crm                    
  real(r8) :: mx_crm                    
  real(r8) :: mui_crm (plev+1)             
  real(r8) :: mdi_crm (plev+1)             
  real(r8) :: flux_qt(plev) 
  real(r8) :: fluxsgs_qt(plev) 
  real(r8) :: tkez(plev) 
  real(r8) :: tkesgsz(plev) 
  real(r8) :: tkz(plev)  
  real(r8) :: flux_u(plev) 
  real(r8) :: flux_v(plev) 
  real(r8) :: flux_qp(plev) 
  real(r8) :: pflx(plev)    
  real(r8) :: qt_ls(plev) 
  real(r8) :: qt_trans(plev)
  real(r8) :: qp_trans(plev) 
  real(r8) :: qp_fall(plev) 
  real(r8) :: qp_src(plev) 
  real(r8) :: qp_evp(plev) 
  real(r8) :: t_ls(plev) 
  real(r8) :: prectend 
  real(r8) :: precstend 
  real(r8) :: precsc 
  real(r8) :: precsl 
  real(r8) :: taux_crm  
  real(r8) :: tauy_crm  
  real(r8) :: z0m 
  real(r8) :: timing_factor 
  real(r8) :: qc_crm (crm_nx, crm_ny, crm_nz)
  real(r8) :: qi_crm (crm_nx, crm_ny, crm_nz)
  real(r8) :: qpc_crm(crm_nx, crm_ny, crm_nz)
  real(r8) :: qpi_crm(crm_nx, crm_ny, crm_nz)
  real(r8) :: prec_crm(crm_nx, crm_ny)
#ifdef ECPP
  real(r8) :: acen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)   
  real(r8) :: acen_tf(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) 
  real(r8) :: rhcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  
  real(r8) :: qcloudcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  
  real(r8) :: qicecen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) 
  real(r8) :: qlsinkcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  
  real(r8) :: precrcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)   
  real(r8) :: precsolidcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)   
  real(r8) :: qlsink_bfcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  
  real(r8) :: qlsink_avgcen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  
  real(r8) :: praincen(plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  
  real(r8) :: wwqui_cen(plev)                                
  real(r8) :: wwqui_cloudy_cen(plev)                         
  real(r8) :: abnd(plev+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)   
  real(r8) :: abnd_tf(plev+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) 
  real(r8) :: massflxbnd(plev+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) 
  real(r8) :: wupthresh_bnd(plev+1)             
  real(r8) :: wdownthresh_bnd(plev+1)           
  real(r8) :: wwqui_bnd(plev+1)                                
  real(r8) :: wwqui_cloudy_bnd(plev+1)                         
#endif
  real(r8) :: qtot(20)



  real(r8) :: u_crm_f  (crm_nx,crm_ny,crm_nz) 
  real(r8) :: v_crm_f  (crm_nx,crm_ny,crm_nz) 
  real(r8) :: w_crm_f  (crm_nx,crm_ny,crm_nz) 
  real(r8) :: t_crm_f  (crm_nx,crm_ny,crm_nz) 



  !I have to call setparm to get the correct value for nmicro_fields
  call setparm()
  allocate( micro_fields_crm  (crm_nx,crm_ny,crm_nz,nmicro_fields+1) )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Read in the command line arguments
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call getarg(1,fname_in)
  call getarg(2,fname_out)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Read in Input Data from the dump file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call dmdf_read(lchnk           ,trim(fname_in),'lchnk'           ,1,.true. ,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(icol            ,trim(fname_in),'icol'            ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(latitude0       ,trim(fname_in),'latitude0'       ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(longitude0      ,trim(fname_in),'longitude0'      ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(ps              ,trim(fname_in),'ps'              ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(pmid            ,trim(fname_in),'pmid'            ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(pdel            ,trim(fname_in),'pdel'            ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(phis            ,trim(fname_in),'phis'            ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(zmid            ,trim(fname_in),'zmid'            ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(zint            ,trim(fname_in),'zint'            ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(qrad_crm        ,trim(fname_in),'qrad_crm'        ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(dt_gl           ,trim(fname_in),'dt_gl'           ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(ocnfrac         ,trim(fname_in),'ocnfrac'         ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(tau00           ,trim(fname_in),'tau00'           ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(wndls           ,trim(fname_in),'wndls'           ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(bflxls          ,trim(fname_in),'bflxls'          ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(fluxu00         ,trim(fname_in),'fluxu00'         ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(fluxv00         ,trim(fname_in),'fluxv00'         ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(fluxt00         ,trim(fname_in),'fluxt00'         ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(fluxq00         ,trim(fname_in),'fluxq00'         ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(tl              ,trim(fname_in),'tl'              ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(ql              ,trim(fname_in),'ql'              ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(qccl            ,trim(fname_in),'qccl'            ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(qiil            ,trim(fname_in),'qiil'            ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(ul              ,trim(fname_in),'ul'              ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(vl              ,trim(fname_in),'vl'              ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
#ifdef CLUBB_CRM
  call dmdf_read(clubb_buffer    ,trim(fname_in),'clubb_buffer'    ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
#endif
  call dmdf_read(cltot           ,trim(fname_in),'cltot'           ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(clhgh           ,trim(fname_in),'clhgh'           ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(clmed           ,trim(fname_in),'clmed'           ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(cllow           ,trim(fname_in),'cllow'           ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(u_crm           ,trim(fname_in),'u_crm'           ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(v_crm           ,trim(fname_in),'v_crm'           ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(w_crm           ,trim(fname_in),'w_crm'           ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(t_crm           ,trim(fname_in),'t_crm'           ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(micro_fields_crm,trim(fname_in),'micro_fields_crm',1,.false.,.false.); _ERR(success,error_string,__LINE__)
#ifdef m2005
#ifdef MODAL_AERO
  call dmdf_read(naermod         ,trim(fname_in),'naermod'         ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(vaerosol        ,trim(fname_in),'vaerosol'        ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(hygro           ,trim(fname_in),'hygro'           ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
#endif
#endif
  call dmdf_read(dd_crm          ,trim(fname_in),'dd_crm'          ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(mui_crm         ,trim(fname_in),'mui_crm'         ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(mdi_crm         ,trim(fname_in),'mdi_crm'         ,1,.false.,.true. ); _ERR(success,error_string,__LINE__)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Call the crm routine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call crm            (lchnk, icol, &
#ifdef CRM_STANDALONE
                       latitude0, longitude0, &
#endif
                       tl, ql, qccl, qiil, ul, vl, &
                       ps, pmid, pdel, phis, &
                       zmid, zint, dt_gl, plev, &
#ifdef CRM3D
                       ultend, vltend,          &
#endif
                       qltend, qcltend, qiltend, sltend, &
                       u_crm, v_crm, w_crm, t_crm, micro_fields_crm, &
                       qrad_crm, &
                       qc_crm, qi_crm, qpc_crm, qpi_crm, prec_crm, &
                       t_rad, qv_rad, qc_rad, qi_rad, cld_rad, cld3d_crm, &
#ifdef m2005
                       nc_rad, ni_rad, qs_rad, ns_rad, wvar_crm,  &
                       aut_crm, acc_crm, evpc_crm, evpr_crm, mlt_crm, &
                       sub_crm, dep_crm, con_crm, &
                       aut_crm_a, acc_crm_a, evpc_crm_a, evpr_crm_a, mlt_crm_a, &
                       sub_crm_a, dep_crm_a, con_crm_a, &
#endif
                       precc, precl, precsc, precsl, &
                       cltot, clhgh, clmed, cllow, cld, cldtop, &
                       gicewp, gliqwp, &
                       mc, mcup, mcdn, mcuup, mcudn, &
                       crm_qc, crm_qi, crm_qs, crm_qg, crm_qr, &
#ifdef m2005
                       crm_nc, crm_ni, crm_ns, crm_ng, crm_nr, &
#ifdef MODAL_AERO
                       naermod, vaerosol, hygro,     &
#endif 
#endif
#ifdef CLUBB_CRM
                       clubb_buffer,                 &
                       crm_cld,                      &
                       clubb_tk, clubb_tkh,          &
                       relvar, accre_enhan, qclvar,  &
#endif
                       crm_tk, crm_tkh,              &
                       mu_crm, md_crm, du_crm, eu_crm, ed_crm, jt_crm, mx_crm,    &
#ifdef ECPP
                       abnd, abnd_tf, massflxbnd, acen, acen_tf,           &
                       rhcen, qcloudcen, qicecen, qlsinkcen, precrcen, precsolidcen,  & 
                       qlsink_bfcen, qlsink_avgcen, praincen,     &
                       wupthresh_bnd, wdownthresh_bnd,   &
                       wwqui_cen, wwqui_bnd, wwqui_cloudy_cen, wwqui_cloudy_bnd,   &
#endif
                       tkez, tkesgsz, tkz, flux_u, flux_v, flux_qt, fluxsgs_qt,flux_qp, &
                       pflx, qt_ls, qt_trans, qp_trans, qp_fall, &
                       qp_evp, qp_src, t_ls, prectend, precstend, &
                       ocnfrac, wndls, tau00, bflxls, &
                       fluxu00, fluxv00, fluxt00, fluxq00,    &
                       taux_crm, tauy_crm, z0m, timing_factor, qtot)   

  call dmdf_read(u_crm_f           ,trim(fname_out),'u_crm'           ,1,.true. ,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(v_crm_f           ,trim(fname_out),'v_crm'           ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(w_crm_f           ,trim(fname_out),'w_crm'           ,1,.false.,.false.); _ERR(success,error_string,__LINE__)
  call dmdf_read(t_crm_f           ,trim(fname_out),'t_crm'           ,1,.false.,.true. ); _ERR(success,error_string,__LINE__)

write(*,*) sum(abs(u_crm-u_crm_f)) / sum(abs(u_crm_f))
write(*,*) sum(abs(v_crm-v_crm_f)) / sum(abs(v_crm_f))
write(*,*) sum(abs(w_crm-w_crm_f)) / sum(abs(w_crm_f))
write(*,*) sum(abs(t_crm-t_crm_f)) / sum(abs(t_crm_f))




end program crm_standalone

