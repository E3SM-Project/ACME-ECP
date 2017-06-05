
begin

f1 = addfile(fn1,"r")
f2 = addfile(fn2,"r")

vnlist = (/   \
"crm_tk",   \
"crm_tkh",   \
"cltot",   \
"clhgh",   \
"clmed",   \
"cllow",   \
"sltend",   \
"u_crm",   \
"v_crm",   \
"w_crm",   \
"t_crm",   \
"micro_fields_crm",   \
"qltend",   \
"qcltend",   \
"qiltend",   \
"t_rad",   \
"qv_rad",   \
"qc_rad",   \
"qi_rad",   \
"cld_rad",   \
"cld3d_crm",   \
"precc",   \
"precl",   \
"cld",   \
"cldtop",   \
"gicewp",   \
"gliqwp",   \
"mc",   \
"mcup",   \
"mcdn",   \
"mcuup",   \
"mcudn",   \
"crm_qc",   \
"crm_qi",   \
"crm_qs",   \
"crm_qg",   \
"crm_qr",   \
"mu_crm",   \
"md_crm",   \
"du_crm",   \
"eu_crm",   \
"ed_crm",   \
"dd_crm",   \
"jt_crm",   \
"mx_crm",   \
"mui_crm",   \
"mdi_crm",   \
"flux_qt",   \
"fluxsgs_qt",   \
"tkez",   \
"tkesgsz",   \
"tkz",   \
"flux_u",   \
"flux_v",   \
"flux_qp",   \
"pflx",   \
"qt_ls",   \
"qt_trans",   \
"qp_trans",   \
"qp_fall",   \
"qp_src",   \
"qp_evp",   \
"t_ls",   \
"prectend",   \
"precstend",   \
"precsc",   \
"precsl",   \
"taux_crm",   \
"tauy_crm",   \
"z0m",   \
"timing_factor",   \
"qc_crm",   \
"qi_crm",   \
"qpc_crm",   \
"qpi_crm",   \
"prec_crm",   \
"qtot"   /)

;#ifdef CRM3D
;    call dmdf_write(ultend          (:,      i),rank,fprefix_out,trim('ultend          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
;    call dmdf_write(vltend          (:,      i),rank,fprefix_out,trim('vltend          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
;#endif
;#ifdef m2005
;    call dmdf_write(nc_rad          (:,:,:,  i),rank,fprefix_out,trim('nc_rad          '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
;    call dmdf_write(ni_rad          (:,:,:,  i),rank,fprefix_out,trim('ni_rad          '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
;    call dmdf_write(qs_rad          (:,:,:,  i),rank,fprefix_out,trim('qs_rad          '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
;    call dmdf_write(ns_rad          (:,:,:,  i),rank,fprefix_out,trim('ns_rad          '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
;    call dmdf_write(wvar_crm        (:,:,:,  i),rank,fprefix_out,trim('wvar_crm        '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
;    call dmdf_write(aut_crm         (:,:,:,  i),rank,fprefix_out,trim('aut_crm         '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
;    call dmdf_write(acc_crm         (:,:,:,  i),rank,fprefix_out,trim('acc_crm         '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
;    call dmdf_write(evpc_crm        (:,:,:,  i),rank,fprefix_out,trim('evpc_crm        '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
;    call dmdf_write(evpr_crm        (:,:,:,  i),rank,fprefix_out,trim('evpr_crm        '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
;    call dmdf_write(mlt_crm         (:,:,:,  i),rank,fprefix_out,trim('mlt_crm         '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
;    call dmdf_write(sub_crm         (:,:,:,  i),rank,fprefix_out,trim('sub_crm         '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
;    call dmdf_write(dep_crm         (:,:,:,  i),rank,fprefix_out,trim('dep_crm         '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
;    call dmdf_write(con_crm         (:,:,:,  i),rank,fprefix_out,trim('con_crm         '),(/'crm_nx','crm_ny','crm_nz'/)                   ,.false.,.false.); _ERR(success,error_string,__LINE__)
;    call dmdf_write(aut_crm_a       (:,      i),rank,fprefix_out,trim('aut_crm_a       '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
;    call dmdf_write(acc_crm_a       (:,      i),rank,fprefix_out,trim('acc_crm_a       '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
;    call dmdf_write(evpc_crm_a      (:,      i),rank,fprefix_out,trim('evpc_crm_a      '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
;    call dmdf_write(evpr_crm_a      (:,      i),rank,fprefix_out,trim('evpr_crm_a      '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
;    call dmdf_write(mlt_crm_a       (:,      i),rank,fprefix_out,trim('mlt_crm_a       '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
;    call dmdf_write(sub_crm_a       (:,      i),rank,fprefix_out,trim('sub_crm_a       '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
;    call dmdf_write(dep_crm_a       (:,      i),rank,fprefix_out,trim('dep_crm_a       '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
;    call dmdf_write(con_crm_a       (:,      i),rank,fprefix_out,trim('con_crm_a       '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
;    call dmdf_write(crm_nc          (:,      i),rank,fprefix_out,trim('crm_nc          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
;    call dmdf_write(crm_ni          (:,      i),rank,fprefix_out,trim('crm_ni          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
;    call dmdf_write(crm_ns          (:,      i),rank,fprefix_out,trim('crm_ns          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
;    call dmdf_write(crm_ng          (:,      i),rank,fprefix_out,trim('crm_ng          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
;    call dmdf_write(crm_nr          (:,      i),rank,fprefix_out,trim('crm_nr          '),(/'plev'/)                                       ,.false.,.false.); _ERR(success,error_string,__LINE__)
;#endif

nvar = dimsizes(vnlist)

do i = 0 , nvar(0)-1
  v1 = f1->$vnlist(i)$
  v2 = f2->$vnlist(i)$
  
  l1 = sum(abs(v2-v1)) / max( (/ 1e-15 , sum(abs(v1)) /) )
  
  print("Relative L1 Norm for "+vnlist(i)+" is: "+l1)
  delete(v1)
  delete(v2)
end do

end

