;==============================================================================================
;==============================================================================================
undef("strpad")
function strpad (str,n)
local str,n,len,ostr,nstr
begin
    ostr = str
    nstr = dimsizes(str)
    do i = 0,nstr-1
        do while(strlen(ostr(i)).lt.n)
            ostr(i) = ostr(i) + " "
        end do
    end do
    return ostr
end
;=======================================================================
;=======================================================================
undef("printMAM")
procedure printMAM (x)
local x
begin
    if isatt(x,"long_name") then
        print(" "+x@long_name)
    else
        print(" (no name)")
    end if
    print("   min: "+min(x))
    print("   avg: "+avg(x))
    print("   max: "+max(x))
end
;==============================================================================================
;==============================================================================================
begin

    f1 = addfile(fn1,"r")
    f2 = addfile(fn2,"r")

    ;------------------------------------------------------------------------
    ;------------------------------------------------------------------------
    vnlist=(/   \
        "crm_tk","crm_tkh",   \
        "cltot","clhgh","clmed","cllow",   \
        "u_crm","v_crm","w_crm","t_crm",   \
        "qc_crm","qi_crm","qpc_crm","qpi_crm","prec_crm","qtot",   \
        "sltend","qltend","qcltend","qiltend",   \
        "micro_fields_crm",   \
        "t_rad","qv_rad","qc_rad","qi_rad",   \
        "cld_rad","cld3d_crm","cld","cldtop",   \
        "precc","precl","prectend","precstend","precsc","precsl",   \
        "gicewp","gliqwp",   \
        "mc","mcup","mcdn","mcuup","mcudn",   \
        "crm_qc","crm_qi","crm_qs","crm_qg","crm_qr",   \
        "mu_crm","md_crm","du_crm","eu_crm","ed_crm","dd_crm",   \
        "jt_crm","mx_crm",   \
        "flux_qt","fluxsgs_qt","flux_u","flux_v","flux_qp","pflx",   \
        "tkez","tkesgsz","tkz",   \
        "qt_ls","t_ls",   \
        "mui_crm","mdi_crm",   \
        "qt_trans","qp_trans","qp_fall","qp_src","qp_evp",   \
        "taux_crm","tauy_crm",   \
        "z0m","timing_factor"   \
        /)
    nvar = dimsizes(vnlist)
    ;------------------------------------------------------------------------
    ;------------------------------------------------------------------------
    ; aok = False

    do i = 0 , nvar(0)-1
        ;----------------------------------------------
        ;----------------------------------------------
        v1 := f1->$vnlist(i)$
        v2 := f2->$vnlist(i)$

        ;----------------------------------------------
        ;----------------------------------------------
        ; tvar = "mui_crm"
        ; if vnlist(i).eq.tvar then
        ;     ; print(dimsizes(v1))
        ;     ; print(dimsizes(v2))
        ;     ; printMAM(v1)
        ;     ; printMAM(v2)
        ;     ilev = ispan(0,72,1)
        ;     ; print(ilev(::-1)+"      "+dim_avg_n(v1,0)+"     "+dim_avg_n(v2,0))
        ;     tt = 0
        ;     print(ilev(::-1)+"      "+v1(tt,:)+"     "+v2(tt,:))
        ;     exit 
        ; end if
        ;----------------------------------------------
        ; check dimension sizes - adjust zero-th
        ;----------------------------------------------
        dims1 := dimsizes(v1)
        dims2 := dimsizes(v2)
        ndim  := dimsizes(dims1)
        t2 = dims1(0)-1
        if .not.all(dims1.eq.dims2) then 
            if i.eq.0 then 
                print("")
                print("dims:")
                print("  "+dims1+"    "+dims2)   
                print("")
            end if

            t2 = min((/dims1(0),dims2(0)/)) - 1

            if ndim.eq.1 then 
                v1 := v1(:t2)
                v2 := v2(:t2)
            else
                if all(dims1(1:).eq.dims2(1:))
                    if ndim.eq.5 then 
                        v1 := v1(:t2,:,:,:,:)
                        v2 := v2(:t2,:,:,:,:)
                    end if
                    if ndim.eq.4 then 
                        v1 := v1(:t2,:,:,:)
                        v2 := v2(:t2,:,:,:)
                    end if
                    if ndim.eq.3 then 
                        v1 := v1(:t2,:,:)
                        v2 := v2(:t2,:,:)
                    end if
                    if ndim.eq.2 then 
                        v1 := v1(:t2,:)
                        v2 := v2(:t2,:)
                    end if
                else
                    print("ERROR: Dimensions sizes for variable '"+vnlist(i)+"' do not match!")
                    print("  "+dims1+"    "+dims2)
                end if
            end if
        end if
        ;----------------------------------------------
        ; special check for a specific variable
        ;----------------------------------------------
        ; tvar = "micro_fields_crm"
        ; if vnlist(i).eq.tvar then
        ;     v1 := f1->$tvar$
        ;     v2 := f2->$tvar$
        ;     dims1 := dimsizes(v1)
        ;     num_m = dims1(1)
        ;     fmt = "%10.2F"
        ;     do m = 0,num_m-1
        ;         ; print( "m = "+m+"  "+ sprintf(fmt,avg(v1(:,m,:,:,:))) +"  "+ sprintf(fmt,avg(v2(:,m,:,:,:))) )
        ;         ; print( "m = "+m+"  "+ sprintf(fmt,max(v1(:,m,:,:,:))) +"  "+ sprintf(fmt,max(v2(:,m,:,:,:))) )
        ;         print( "m = "+m+"  "+ sprintf(fmt,avg(v1(:t2,m,:,:,:))) +"  "+ sprintf(fmt,avg(v2(:t2,m,:,:,:))) +"  "+ sprintf(fmt,avg(v1(:t2,m,:,:,:)-v2(:t2,m,:,:,:))) )
        ;     end do
        ;     exit 
        ; end if
        ;----------------------------------------------
        ; calculate difference statistics
        ;----------------------------------------------
        if ndim.eq.5 then 
            num_m = dims1(1)
        else
            num_m = 1
        end if

        do m = 0,num_m-1

            if ndim.eq.5 then
                tvar = vnlist(i) + " ( m = "+sprinti("%2.2i",m)+") "
                tv1 := v1(:,m,:,:,:)
                tv2 := v2(:,m,:,:,:)
            else
                tvar = vnlist(i)
                tv1 := v1
                tv2 := v2
            end if

            l1 = sum(abs(tv2-tv1)) / max( (/ 1e-15 , sum(abs(v1)) /) )

            avg_diff = avg( abs(tv2-tv1) )
            max_diff = max( abs(tv2-tv1) )

            fmt1 = "%10.4f"

            l1_str = sprintf(fmt1,l1)
            av_str = sprintf(fmt1,avg_diff)
            xd_str = sprintf(fmt1,max_diff)
            
            theshold = 1e-6
            if l1.lt.theshold       then l1_str = "" end if
            if avg_diff.lt.theshold then av_str = "" end if
            if max_diff.lt.theshold then xd_str = "" end if

            padlen1 = 30
            padlen2 = 20
            
            flag = ""
            if any(isnan_ieee((/tv1,tv2/))) then flag = flag + " !!!! " end if
            if any(isnan_ieee(tv1)) then flag = flag + " 1111 " end if
            if any(isnan_ieee(tv2)) then flag = flag + " 2222 " end if

            print(strpad(tvar,padlen1)\
                  +"  L1 Norm  : "+strpad( l1_str  ,padlen2)  \
                  +"  avg diff : "+strpad( av_str  ,padlen2)  \
                  +"  max diff : "+strpad( xd_str  ,padlen2)  \
                  +flag)
        
        end do
        ;----------------------------------------------
        ;----------------------------------------------
    end do
    ;------------------------------------------------------------------------
    ;------------------------------------------------------------------------
end
;==============================================================================================
;==============================================================================================


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