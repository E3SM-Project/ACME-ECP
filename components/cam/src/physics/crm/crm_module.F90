
module crm_module
!---------------------------------------------------------------
!  Super-parameterization's main driver 
!  Marat Khairoutdinov, 2001-2009
!---------------------------------------------------------------
use setparm_mod, only : setparm

contains

! subroutine crm  (lchnk, icol, &
subroutine crm(lchnk, icol, nvcols, &
!MRN: If this is in standalone mode, lat,lon are passed in directly, not looked up in phys_grid
#ifdef CRM_STANDALONE
                latitude0_in, longitude0_in, &
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
        !---------------------------------------------------------------
    use crm_dump              , only: crm_dump_input, crm_dump_output
    use shr_kind_mod          , only: r8 => shr_kind_r8
    !MRN: phys_grid is a rabbit hole of dependencies I'd rather hijack and avoid.
    !MRN: It's only used to get the longitude and latitude for each call and then
    !MRN: random seed values for perturbations to the initial temprature field on 
    !MRN: the first call to crm(...)

#ifndef CRM_STANDALONE
    use phys_grid             , only: get_rlon_p, get_rlat_p, get_gcol_all_p
#endif
    use ppgrid                , only: pcols
    use vars
    use params
    use microphysics
    use sgs
    use crmtracers

! whannah - Matt Norman added the more specific use statements below - however this was causing problems for the 1-moment configuration
! 
!     use vars                  , only: crm_nx, crm_ny, crm_nz, nx, ny, nz, crm_dx, crm_dt, fcory, fcorzy, latitude, longitude, z, zi, pres, prespot, bet, gamaz, &
!                                       adzw, adz, rho, rhow, u, v, w, tabs, dudt, dvdt, dwdt, p, wsub, cf3d, u0, v0, t0, tabs0, q0, t, qp0, tke0, qv0, qn0,  &
!                                       tke0, ttend, qtend, utend, vtend, ug0, vg0, tg0, qg0, qtotmicro, dt3, precsfc, precssfc, qpfall, precflux, uwle, uwsb, &
!                                       vwle, vwsb, dostatis, nzm, dz, yes3d, qcl, qci, qpl, qpi, qv, fluxbu, fluxbv, fluxbt, fluxbq, fluxtu, fluxtv, fluxtt, &
!                                       fluxtq, fzero, uwle, uwsb, vwle, nstop, dt, nsave3d, nstep, nprint, ncycle, day, day0, time, icycle, dtn, na, dtfactor, &
! #ifdef MODAL_AERO
!                                       naer, vaer, hgaer, &
! #endif
!                                       nc, nb, qsatw_crm
!     use params                , only: latitude0, longitude0, fcor, pi, fcorz, ocean, land, rgas, cp, fac_fus, uhl, vhl, z0, dodamping, dosurface, docoriolis, &
!                                       taux0, tauy0, doclubb, doclubbnoninter, docloud, dosmoke, crm_rknd
!     use microphysics          , only: nmicro_fields, micro_field, cloudliq, aut1a, acc1a, evpc1a, evpr1a, mlt1a, sub1a, dep1a, con1a, mkwsb, mkwle, mkadv, &
!                                       mkdiff, qpsrc, qpevp, ggr, dopredictnc, incl, nc0, fac_cond, fac_sub, iqr, iqs, iqg, inci, ins, wvar, aut1, acc1,    &
!                                       evpc1, evpr1, mlt1, sub1, dep1, con1, inr, ing, iqv, iqci, micro_init, micro_proc
!     use sgs                   , only: tke, tk, tkh, dosgs, sgs_init, sgs_proc, sgs_mom, sgs_scalars

#ifdef MODAL_AERO
    use modal_aero_data       , only: ntot_amode
#endif
#ifdef CLUBB_CRM
    use crmdims               , only: nclubbvars
    use clubb_sgs             , only: advance_clubb_sgs, clubb_sgs_setup, clubb_sgs_cleanup, apply_clubb_sgs_tndcy, apply_clubb_sgs_tndcy_scalars, &
                                      apply_clubb_sgs_tndcy_mom, t2thetal, total_energy
    use clubb_precision       , only: time_precision, core_rknd
    use clubbvars             , only: up2, vp2, wprtp, wpthlp, wp2, wp3, rtp2, thlp2, rtpthlp, upwp, vpwp, cloud_frac, t_tndcy, qc_tndcy, qv_tndcy, &
                                      u_tndcy, v_tndcy, lrestart_clubb, rho_ds_zt, rho_ds_zm, thv_ds_zt, thv_ds_zm, invrs_rho_ds_zt, invrs_rho_ds_zm, &
                                      tracer_tndcy, sclrp2, sclrprtp, sclrpthlp, wpsclrp, relvarg, accre_enhang, qclvarg, edsclr_dim, sclr_dim, rho_ds_zt, &
                                      rho_ds_zm, rtm_spurious_source, thlm_spurious_source
    use fill_holes            , only: vertical_integral
    use numerical_check       , only: calculate_spurious_source
    use grid_class            , only: gr
#endif
#ifdef ECPP
    use ecppvars              , only: qlsink, precr, precsolid, &
                                      area_bnd_final, area_bnd_sum, area_cen_final, area_cen_sum, &
                                      mass_bnd_final, mass_bnd_sum, rh_cen_sum, qcloud_cen_sum, qice_cen_sum, &
                                      qlsink_cen_sum, precr_cen_sum, precsolid_cen_sum, xkhvsum, wup_thresh, wdown_thresh, &
                                      wwqui_cen_sum, wwqui_bnd_sum, wwqui_cloudy_cen_sum, wwqui_cloudy_bnd_sum, &
                                      qlsink_bf_cen_sum, qlsink_avg_cen_sum, prain_cen_sum, qlsink_bf, prain
    use module_ecpp_crm_driver, only: ecpp_crm_stat, ecpp_crm_init, ecpp_crm_cleanup, ntavg1_ss, ntavg2_ss
    use ecppvars              , only: NCLASS_CL, ncls_ecpp_in, NCLASS_PR
#endif
    use cam_abortutils        , only: endrun
    use time_manager          , only: get_nstep

    implicit none
    integer , intent(in   ) :: lchnk                            ! chunk identifier (only for lat/lon and random seed)
    integer , intent(in   ) :: nvcols                           ! Number of "vector" GCM columns to push down into CRM for SIMD vectorization / more threading
    integer , intent(in   ) :: plev                             ! number of levels in parent model
    real(r8), intent(in   ) :: dt_gl                            ! global model's time step
    integer , intent(in   ) :: icol                (nvcols)     ! column identifier (only for lat/lon and random seed)
#ifdef CRM_STANDALONE
    real(crm_rknd)   , intent(in) :: latitude0_in  (nvcols) 
    real(crm_rknd)   , intent(in) :: longitude0_in (nvcols)
#endif
    real(r8), intent(in   ) :: ps                  (nvcols)       ! Global grid surface pressure (Pa)
    real(r8), intent(in   ) :: pmid                (nvcols,plev)  ! Global grid pressure (Pa)
    real(r8), intent(in   ) :: pdel                (nvcols,plev)  ! Layer's pressure thickness (Pa)
    real(r8), intent(in   ) :: phis                (nvcols)       ! Global grid surface geopotential (m2/s2)
    real(r8), intent(in   ) :: zmid                (nvcols,plev)  ! Global grid height (m)
    real(r8), intent(in   ) :: zint                (nvcols,plev+1)! Global grid interface height (m)
    real(r8), intent(in   ) :: qrad_crm            (nvcols,crm_nx, crm_ny, crm_nz) ! CRM rad. heating
    real(r8), intent(in   ) :: ocnfrac             (nvcols)       ! area fraction of the ocean
    real(r8), intent(in   ) :: tau00               (nvcols)       ! large-scale surface stress (N/m2)
    real(r8), intent(in   ) :: wndls               (nvcols)       ! large-scale surface wind (m/s)
    real(r8), intent(in   ) :: bflxls              (nvcols)       ! large-scale surface buoyancy flux (K m/s)
    real(r8), intent(in   ) :: fluxu00             (nvcols)       ! surface momenent fluxes [N/m2]
    real(r8), intent(in   ) :: fluxv00             (nvcols)       ! surface momenent fluxes [N/m2]
    real(r8), intent(in   ) :: fluxt00             (nvcols)       ! surface sensible heat fluxes [K Kg/ (m2 s)]
    real(r8), intent(in   ) :: fluxq00             (nvcols)       ! surface latent heat fluxes [ kg/(m2 s)]
    real(r8), intent(in   ) :: tl                  (nvcols,plev)  ! Global grid temperature (K)
    real(r8), intent(in   ) :: ql                  (nvcols,plev)  ! Global grid water vapor (g/g)
    real(r8), intent(in   ) :: qccl                (nvcols,plev)  ! Global grid cloud liquid water (g/g)
    real(r8), intent(in   ) :: qiil                (nvcols,plev)  ! Global grid cloud ice (g/g)
    real(r8), intent(in   ) :: ul                  (nvcols,plev)  ! Global grid u (m/s)
    real(r8), intent(in   ) :: vl                  (nvcols,plev)  ! Global grid v (m/s)
#ifdef CLUBB_CRM
    real(r8), intent(inout), target :: clubb_buffer(nvcols,crm_nx, crm_ny, crm_nz+1,1:nclubbvars)
    real(r8), intent(  out) :: crm_cld             (nvcols,crm_nx, crm_ny, crm_nz+1)
    real(r8), intent(  out) :: clubb_tk            (nvcols,crm_nx, crm_ny, crm_nz)
    real(r8), intent(  out) :: clubb_tkh           (nvcols,crm_nx, crm_ny, crm_nz)
    real(r8), intent(  out) :: relvar              (nvcols,crm_nx, crm_ny, crm_nz) 
    real(r8), intent(  out) :: accre_enhan         (nvcols,crm_nx, crm_ny, crm_nz)
    real(r8), intent(  out) :: qclvar              (nvcols,crm_nx, crm_ny, crm_nz)
#endif
    real(r8), intent(  out) :: crm_tk              (nvcols,crm_nx, crm_ny, crm_nz)
    real(r8), intent(  out) :: crm_tkh             (nvcols,crm_nx, crm_ny, crm_nz)
    real(r8), intent(inout) :: cltot               (nvcols)                        ! shaded cloud fraction
    real(r8), intent(inout) :: clhgh               (nvcols)                        ! shaded cloud fraction
    real(r8), intent(inout) :: clmed               (nvcols)                        ! shaded cloud fraction
    real(r8), intent(inout) :: cllow               (nvcols)                        ! shaded cloud fraction
#ifdef CRM3D
    real(r8), intent(  out) :: ultend              (nvcols,plev)                   ! tendency of ul
    real(r8), intent(  out) :: vltend              (nvcols,plev)                   ! tendency of vl
#endif
    real(r8), intent(  out) :: sltend              (nvcols,plev)                   ! tendency of static energy
    real(r8), intent(  out) :: qltend              (nvcols,plev)                   ! tendency of water vapor
    real(r8), intent(  out) :: qcltend             (nvcols,plev)                   ! tendency of cloud liquid water
    real(r8), intent(  out) :: qiltend             (nvcols,plev)                   ! tendency of cloud ice
    real(r8), intent(inout) :: u_crm               (nvcols,crm_nx,crm_ny,crm_nz)   ! CRM v-wind component
    real(r8), intent(inout) :: v_crm               (nvcols,crm_nx,crm_ny,crm_nz)   ! CRM v-wind component
    real(r8), intent(inout) :: w_crm               (nvcols,crm_nx,crm_ny,crm_nz)   ! CRM w-wind component
    real(r8), intent(inout) :: t_crm               (nvcols,crm_nx,crm_ny,crm_nz)   ! CRM temperuture
    real(r8), intent(inout) :: micro_fields_crm    (nvcols,crm_nx,crm_ny,crm_nz,nmicro_fields+1) ! CRM total water
    real(r8), intent(  out) :: t_rad               (nvcols,crm_nx, crm_ny, crm_nz) ! rad temperuture
    real(r8), intent(  out) :: qv_rad              (nvcols,crm_nx, crm_ny, crm_nz) ! rad vapor
    real(r8), intent(  out) :: qc_rad              (nvcols,crm_nx, crm_ny, crm_nz) ! rad cloud water
    real(r8), intent(  out) :: qi_rad              (nvcols,crm_nx, crm_ny, crm_nz) ! rad cloud ice
    real(r8), intent(  out) :: cld_rad             (nvcols,crm_nx, crm_ny, crm_nz) ! rad cloud fraction 
    real(r8), intent(  out) :: cld3d_crm           (nvcols,crm_nx, crm_ny, crm_nz) ! instant 3D cloud fraction
#ifdef m2005
    real(r8), intent(  out) :: nc_rad              (nvcols,crm_nx, crm_ny, crm_nz) ! rad cloud droplet number (#/kg) 
    real(r8), intent(  out) :: ni_rad              (nvcols,crm_nx, crm_ny, crm_nz) ! rad cloud ice crystal number (#/kg)
    real(r8), intent(  out) :: qs_rad              (nvcols,crm_nx, crm_ny, crm_nz) ! rad cloud snow (kg/kg)
    real(r8), intent(  out) :: ns_rad              (nvcols,crm_nx, crm_ny, crm_nz) ! rad cloud snow crystal number (#/kg)
    real(r8), intent(  out) :: wvar_crm            (nvcols,crm_nx, crm_ny, crm_nz) ! vertical velocity variance (m/s)
    real(r8), intent(  out) :: aut_crm             (nvcols,crm_nx, crm_ny, crm_nz) ! cloud water autoconversion (1/s)
    real(r8), intent(  out) :: acc_crm             (nvcols,crm_nx, crm_ny, crm_nz) ! cloud water accretion (1/s)
    real(r8), intent(  out) :: evpc_crm            (nvcols,crm_nx, crm_ny, crm_nz) ! cloud water evaporation (1/s)
    real(r8), intent(  out) :: evpr_crm            (nvcols,crm_nx, crm_ny, crm_nz) ! rain evaporation (1/s)
    real(r8), intent(  out) :: mlt_crm             (nvcols,crm_nx, crm_ny, crm_nz) ! ice, snow, graupel melting (1/s)
    real(r8), intent(  out) :: sub_crm             (nvcols,crm_nx, crm_ny, crm_nz) ! ice, snow, graupel sublimation (1/s)
    real(r8), intent(  out) :: dep_crm             (nvcols,crm_nx, crm_ny, crm_nz) ! ice, snow, graupel deposition (1/s)
    real(r8), intent(  out) :: con_crm             (nvcols,crm_nx, crm_ny, crm_nz) ! cloud water condensation(1/s)
    real(r8), intent(  out) :: aut_crm_a           (nvcols,plev)  ! cloud water autoconversion (1/s)
    real(r8), intent(  out) :: acc_crm_a           (nvcols,plev)  ! cloud water accretion (1/s)
    real(r8), intent(  out) :: evpc_crm_a          (nvcols,plev)  ! cloud water evaporation (1/s)
    real(r8), intent(  out) :: evpr_crm_a          (nvcols,plev)  ! rain evaporation (1/s)
    real(r8), intent(  out) :: mlt_crm_a           (nvcols,plev)  ! ice, snow, graupel melting (1/s)
    real(r8), intent(  out) :: sub_crm_a           (nvcols,plev)  ! ice, snow, graupel sublimation (1/s)
    real(r8), intent(  out) :: dep_crm_a           (nvcols,plev)  ! ice, snow, graupel deposition (1/s)
    real(r8), intent(  out) :: con_crm_a           (nvcols,plev)  ! cloud water condensation(1/s)
#endif
    real(r8), intent(  out) :: precc               (nvcols)       ! convective precip rate (m/s)
    real(r8), intent(  out) :: precl               (nvcols)       ! stratiform precip rate (m/s)
    real(r8), intent(  out) :: cld                 (nvcols,plev)  ! cloud fraction
    real(r8), intent(  out) :: cldtop              (nvcols,plev)  ! cloud top pdf
    real(r8), intent(  out) :: gicewp              (nvcols,plev)  ! ice water path
    real(r8), intent(  out) :: gliqwp              (nvcols,plev)  ! ice water path
    real(r8), intent(  out) :: mc                  (nvcols,plev)  ! cloud mass flux
    real(r8), intent(  out) :: mcup                (nvcols,plev)  ! updraft cloud mass flux
    real(r8), intent(  out) :: mcdn                (nvcols,plev)  ! downdraft cloud mass flux
    real(r8), intent(  out) :: mcuup               (nvcols,plev)  ! unsat updraft cloud mass flux
    real(r8), intent(  out) :: mcudn               (nvcols,plev)  ! unsat downdraft cloud mass flux
    real(r8), intent(  out) :: crm_qc              (nvcols,plev)  ! mean cloud water
    real(r8), intent(  out) :: crm_qi              (nvcols,plev)  ! mean cloud ice
    real(r8), intent(  out) :: crm_qs              (nvcols,plev)  ! mean snow
    real(r8), intent(  out) :: crm_qg              (nvcols,plev)  ! mean graupel
    real(r8), intent(  out) :: crm_qr              (nvcols,plev)  ! mean rain
#ifdef m2005
    real(r8), intent(  out) :: crm_nc              (nvcols,plev)  ! mean cloud water  (#/kg)
    real(r8), intent(  out) :: crm_ni              (nvcols,plev)  ! mean cloud ice    (#/kg)
    real(r8), intent(  out) :: crm_ns              (nvcols,plev)  ! mean snow         (#/kg)
    real(r8), intent(  out) :: crm_ng              (nvcols,plev)  ! mean graupel      (#/kg)
    real(r8), intent(  out) :: crm_nr              (nvcols,plev)  ! mean rain         (#/kg)
#ifdef MODAL_AERO
    real(r8), intent(in   )  :: naermod            (nvcols,plev, ntot_amode)    ! Aerosol number concentration [/m3]
    real(r8), intent(in   )  :: vaerosol           (nvcols,plev, ntot_amode)    ! aerosol volume concentration [m3/m3]
    real(r8), intent(in   )  :: hygro              (nvcols,plev, ntot_amode)    ! hygroscopicity of aerosol mode 
#endif 
#endif
    real(r8), intent(  out) :: mu_crm              (nvcols,plev)       ! mass flux up
    real(r8), intent(  out) :: md_crm              (nvcols,plev)       ! mass flux down
    real(r8), intent(  out) :: du_crm              (nvcols,plev)       ! mass detrainment from updraft
    real(r8), intent(  out) :: eu_crm              (nvcols,plev)       ! mass entrainment from updraft
    real(r8), intent(  out) :: ed_crm              (nvcols,plev)       ! mass detrainment from downdraft
    real(r8)                :: dd_crm              (nvcols,plev)       ! mass entraiment from downdraft
    real(r8), intent(  out) :: jt_crm              (nvcols)            ! index of cloud (convection) top 
    real(r8), intent(  out) :: mx_crm              (nvcols)            ! index of cloud (convection) bottom
    real(r8)                :: mui_crm             (nvcols,plev+1)     ! mass flux up at the interface
    real(r8)                :: mdi_crm             (nvcols,plev+1)     ! mass flux down at the interface
    real(r8), intent(  out) :: flux_qt             (nvcols,plev)       ! nonprecipitating water flux           [kg/m2/s]
    real(r8), intent(  out) :: fluxsgs_qt          (nvcols,plev)       ! sgs nonprecipitating water flux    [kg/m2/s]
    real(r8), intent(  out) :: tkez                (nvcols,plev)       ! tke profile               [kg/m/s2]
    real(r8), intent(  out) :: tkesgsz             (nvcols,plev)       ! sgs tke profile        [kg/m/s2]
    real(r8), intent(  out) :: tkz                 (nvcols,plev)       ! tk profile                [m2/s]
    real(r8), intent(  out) :: flux_u              (nvcols,plev)       ! x-momentum flux          [m2/s2]
    real(r8), intent(  out) :: flux_v              (nvcols,plev)       ! y-momentum flux          [m2/s2]
    real(r8), intent(  out) :: flux_qp             (nvcols,plev)       ! precipitating water flux [kg/m2/s or mm/s]
    real(r8), intent(  out) :: pflx                (nvcols,plev)       ! precipitation flux      [m/s]
    real(r8), intent(  out) :: qt_ls               (nvcols,plev)       ! tendency of nonprec water due to large-scale  [kg/kg/s]
    real(r8), intent(  out) :: qt_trans            (nvcols,plev)       ! tendency of nonprec water due to transport  [kg/kg/s]
    real(r8), intent(  out) :: qp_trans            (nvcols,plev)       ! tendency of prec water due to transport [kg/kg/s]
    real(r8), intent(  out) :: qp_fall             (nvcols,plev)       ! tendency of prec water due to fall-out   [kg/kg/s]
    real(r8), intent(  out) :: qp_src              (nvcols,plev)       ! tendency of prec water due to conversion  [kg/kg/s]
    real(r8), intent(  out) :: qp_evp              (nvcols,plev)       ! tendency of prec water due to evp         [kg/kg/s]
    real(r8), intent(  out) :: t_ls                (nvcols,plev)       ! tendency of lwse  due to large-scale        [kg/kg/s] ???
    real(r8), intent(  out) :: prectend            (nvcols)            ! column integrated tendency in precipitating water+ice (kg/m2/s)
    real(r8), intent(  out) :: precstend           (nvcols)            ! column integrated tendency in precipitating ice (kg/m2/s)
    real(r8), intent(  out) :: precsc              (nvcols)            ! convective snow rate (m/s)
    real(r8), intent(  out) :: precsl              (nvcols)            ! stratiform snow rate (m/s)
    real(r8), intent(  out) :: taux_crm            (nvcols)            ! zonal CRM surface stress perturbation (N/m2)
    real(r8), intent(  out) :: tauy_crm            (nvcols)            ! merid CRM surface stress perturbation (N/m2)
    real(r8), intent(  out) :: z0m                 (nvcols)            ! surface stress (N/m2)
    real(r8), intent(  out) :: timing_factor       (nvcols)            ! crm cpu efficiency
    real(r8), intent(  out) :: qc_crm              (nvcols,crm_nx, crm_ny, crm_nz)! CRM cloud water
    real(r8), intent(  out) :: qi_crm              (nvcols,crm_nx, crm_ny, crm_nz)! CRM cloud ice
    real(r8), intent(  out) :: qpc_crm             (nvcols,crm_nx, crm_ny, crm_nz)! CRM precip water
    real(r8), intent(  out) :: qpi_crm             (nvcols,crm_nx, crm_ny, crm_nz)! CRM precip ice
    real(r8), intent(  out) :: prec_crm            (nvcols,crm_nx, crm_ny)        ! CRM precipiation rate at layer center
#ifdef ECPP
    real(r8), intent(  out) :: acen                (nvcols,plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! cloud fraction for each sub-sub class for full time period
    real(r8), intent(  out) :: acen_tf             (nvcols,plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! cloud fraction for end-portion of time period
    real(r8), intent(  out) :: rhcen               (nvcols,plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! relative humidity (0-1)
    real(r8), intent(  out) :: qcloudcen           (nvcols,plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! cloud water (kg/kg)
    real(r8), intent(  out) :: qicecen             (nvcols,plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! cloud ice (kg/kg)
    real(r8), intent(  out) :: qlsinkcen           (nvcols,plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! cloud water loss rate from precipitation (/s??)
    real(r8), intent(  out) :: precrcen            (nvcols,plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! liquid (rain) precipitation rate (kg/m2/s)
    real(r8), intent(  out) :: precsolidcen        (nvcols,plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! solid (rain) precipitation rate (kg/m2/s)
    real(r8), intent(  out) :: qlsink_bfcen        (nvcols,plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! cloud water loss rate from precipitation calculated 
                                                                                                   ! cloud water before precipitatinog (/s)
    real(r8), intent(  out) :: qlsink_avgcen       (nvcols,plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! cloud water loss rate from precipitation calculated 
                                                                                                   ! from praincen and qlcoudcen averaged over 
                                                                                                   ! ntavg1_ss time step (/s??)
    real(r8), intent(  out) :: praincen            (nvcols,plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! cloud water loss rate from precipitation (kg/kg/s)
    real(r8), intent(  out) :: wwqui_cen           (nvcols,plev)                                   ! vertical velocity variance in quiescent class (m2/s2)
    real(r8), intent(  out) :: wwqui_cloudy_cen    (nvcols,plev)                                   ! vertical velocity variance in quiescent, and cloudy class (m2/s2) at layer boundary
    real(r8), intent(  out) :: abnd                (nvcols,plev+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)! cloud fraction for each sub-sub class for full time period
    real(r8), intent(  out) :: abnd_tf             (nvcols,plev+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)! cloud fraction for end-portion of time period
    real(r8), intent(  out) :: massflxbnd          (nvcols,plev+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)! sub-class vertical mass flux (kg/m2/s) at layer bottom boundary.
    real(r8), intent(  out) :: wupthresh_bnd       (nvcols,plev+1)                                 ! vertical velocity threshold for updraft (m/s)
    real(r8), intent(  out) :: wdownthresh_bnd     (nvcols,plev+1)                                 ! vertical velocity threshold for downdraft (m/s)
    real(r8), intent(  out) :: wwqui_bnd           (nvcols,plev+1)                                 ! vertical velocity variance in quiescent class (m2/s2)
    real(r8), intent(  out) :: wwqui_cloudy_bnd    (nvcols,plev+1)                                 ! vertical velocity variance in quiescent, and cloudy class (m2/s2)
#endif

!  Local space:
    real(r8), intent(  out) :: qtot                (nvcols,20)

    !  Local space:
    real(r8),       parameter :: umax = 0.5*crm_dx/crm_dt ! maxumum ampitude of the l.s. wind
    real(r8),       parameter :: wmin = 2.                ! minimum up/downdraft velocity for stat
    real(crm_rknd), parameter :: cwp_threshold = 0.001    ! threshold for cloud condensate for shaded fraction calculation
    real(crm_rknd)  :: dummy(nz), t00(nz)
    real(crm_rknd)  :: fluxbtmp(nx,ny), fluxttmp(nx,ny)    !bloss
    real(crm_rknd)  :: tln(plev), qln(plev), qccln(plev), qiiln(plev), uln(plev), vln(plev)
    real(crm_rknd)  :: cwp(nx,ny), cwph(nx,ny), cwpm(nx,ny), cwpl(nx,ny)
    real(r8)        :: factor_xy, idt_gl
    real(crm_rknd)  :: tmp1, tmp2
    real(crm_rknd)  :: u2z,v2z,w2z
    integer         :: i,j,k,l,ptop,nn,icyc, nstatsteps, vc
    integer         :: kx
    logical         :: flag_top(nx,ny)
    real(crm_rknd)  :: ustar, bflx, wnd, z0_est, qsat, omg
    real(crm_rknd)  :: colprec,colprecs
    real(r8)        :: zs                ! surface elevation
    integer         :: igstep            ! GCM time steps
    integer         :: iseed             ! seed for random perturbation
    integer         :: gcolindex(pcols)  ! array of global latitude indices
    real(crm_rknd)  :: cltemp(nx,ny), cmtemp(nx,ny), chtemp(nx, ny), cttemp(nx, ny)
    real(crm_rknd)  :: ntotal_step
    integer         :: myrank, ierr
#ifdef CLUBB_CRM
    !Array indicies for spurious RTM check
    real(kind=core_rknd) :: rtm_integral_before (nx,ny), rtm_integral_after (nx,ny), rtm_flux_top, rtm_flux_sfc
    real(kind=core_rknd) :: thlm_integral_before(nx,ny), thlm_integral_after(nx,ny), thlm_before(nzm), thlm_after(nzm), thlm_flux_top, thlm_flux_sfc
    real(kind=core_rknd) :: rtm_column(nzm) ! Total water (vapor + liquid)     [kg/kg]
#endif

  !Loop over "vector columns"
  do vc = 1 , nvcols


    !MRN: In standalone mode, we need to pass these things in by parameter, not look them up.
#ifdef CRM_STANDALONE
    latitude0  = latitude0_in (vc)
    longitude0 = longitude0_in(vc)
#else
    latitude0  = get_rlat_p(lchnk, icol(vc)) * 57.296_r8
    longitude0 = get_rlon_p(lchnk, icol(vc)) * 57.296_r8
#endif

    igstep = get_nstep()

    call crm_dump_input( igstep,plev,lchnk,icol(vc),latitude0,longitude0,ps(vc),pmid(vc,:),pdel(vc,:),phis(vc),zmid(vc,:),zint(vc,:),qrad_crm(vc,:,:,:),dt_gl, &
                         ocnfrac(vc),tau00(vc),wndls(vc),bflxls(vc),fluxu00(vc),fluxv00(vc),fluxt00(vc),fluxq00(vc),tl(vc,:),ql(vc,:),qccl(vc,:),qiil(vc,:),   &
                         ul(vc,:),vl(vc,:), &
#ifdef CLUBB_CRM
                         clubb_buffer(vc,:,:,:,:) , &
#endif
                         cltot(vc),clhgh(vc),clmed(vc),cllow(vc),u_crm(vc,:,:,:),v_crm(vc,:,:,:),w_crm(vc,:,:,:),t_crm(vc,:,:,:),micro_fields_crm(vc,:,:,:,:), &
#ifdef m2005
#ifdef MODAL_AERO
                         naermod(vc,:,:),vaerosol(vc,:,:),hygro(vc,:,:) , &
#endif
#endif
                         dd_crm(vc,:),mui_crm(vc,:),mdi_crm(vc,:) )

!-----------------------------------------------

    dostatis  = .false.    ! no statistics are collected. 
    idt_gl    = 1._r8/dt_gl
    ptop      = plev-nzm+1
    factor_xy = 1._r8/dble(nx*ny)
    dummy     = 0.
    t_rad  (vc,:,:,:) = 0.
    qv_rad (vc,:,:,:) = 0.
    qc_rad (vc,:,:,:) = 0.
    qi_rad (vc,:,:,:) = 0.
    cld_rad(vc,:,:,:) = 0.
#ifdef m2005
    nc_rad(vc,:,:,:) = 0.0
    ni_rad(vc,:,:,:) = 0.0
    qs_rad(vc,:,:,:) = 0.0
    ns_rad(vc,:,:,:) = 0.0
#endif
    zs=phis(vc)/ggr
    bflx = bflxls(vc)
    wnd = wndls(vc)

!-----------------------------------------

#ifdef CLUBB_CRM
    if(igstep == 1) then
      lrestart_clubb = .false.
    else
     lrestart_clubb = .true.
    endif
#endif

    call task_init ()
    call setparm()

    if(fcor.eq.-999.) fcor= 4*pi/86400.*sin(latitude0*pi/180.)
    fcorz = sqrt(4.*(2*pi/(3600.*24.))**2-fcor**2)
    fcory(:) = fcor
    fcorzy(:) = fcorz
    do j=1,ny
      do i=1,nx
        latitude (i,j) = latitude0
        longitude(i,j) = longitude0
      end do
    end do

    if(ocnfrac(vc).gt.0.5) then
       OCEAN = .true.
    else
       LAND = .true.
    end if

    ! Create CRM vertical grid and initialize some vertical reference arrays:
    do k = 1, nzm
      z(k) = zmid(vc,plev-k+1) - zint(vc,plev+1)
      zi(k) = zint(vc,plev-k+2)- zint(vc,plev+1)
      pres(k) = pmid(vc,plev-k+1)/100.
      prespot(k)=(1000./pres(k))**(rgas/cp)
      bet(k) = ggr/tl(vc,plev-k+1)
      gamaz(k)=ggr/cp*z(k)
    end do ! k        
   ! zi(nz) =  zint(plev-nz+2)
    zi(nz) = zint(vc,plev-nz+2)-zint(vc,plev+1) !+++mhwang, 2012-02-04

    dz = 0.5*(z(1)+z(2))
    do k=2,nzm
      adzw(k) = (z(k)-z(k-1))/dz
    end do
    adzw(1)  = 1.
    adzw(nz) = adzw(nzm)
    !+++mhwang fix the adz bug. (adz needs to be consistent with zi)
    !2012-02-04 Minghuai Wang (minghuai.wang@pnnl.gov)
    do k=1, nzm
      adz(k)=(zi(k+1)-zi(k))/dz
    end do
    
    do k = 1,nzm
      rho(k) = pdel(vc,plev-k+1)/ggr/(adz(k)*dz)
    end do
    do k=2,nzm
    ! rhow(k) = 0.5*(rho(k)+rho(k-1))
    !+++mhwang fix the rhow bug (rhow needes to be consistent with pmid)
    !2012-02-04 Minghuai Wang (minghuai.wang@pnnl.gov)
      rhow(k) = (pmid(vc,plev-k+2)-pmid(vc,plev-k+1))/ggr/(adzw(k)*dz)
    end do
    rhow(1) = 2.*rhow(2) - rhow(3)
#ifdef CLUBB_CRM /* Fix extrapolation for 30 point grid */
    if (  2.*rhow(nzm) - rhow(nzm-1) > 0. ) then
       rhow(nz)= 2.*rhow(nzm) - rhow(nzm-1)
    else
       rhow(nz)= sqrt( rhow(nzm) )
    endif
#else
    rhow(nz)= 2.*rhow(nzm) - rhow(nzm-1)


#endif /*CLUBB_CRM*/
    colprec=0
    colprecs=0

    !  Initialize:
    ! limit the velocity at the very first step:
    if(u_crm(vc,1,1,1).eq.u_crm(vc,2,1,1).and.u_crm(vc,3,1,2).eq.u_crm(vc,4,1,2)) then
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            u_crm(vc,i,j,k) = min( umax, max(-umax,u_crm(vc,i,j,k)) )
            v_crm(vc,i,j,k) = min( umax, max(-umax,v_crm(vc,i,j,k)) )*YES3D
          enddo
        enddo
      enddo
    endif

    u          (1:nx,1:ny,1:nzm                ) = u_crm           (vc,1:nx,1:ny,1:nzm                )
    v          (1:nx,1:ny,1:nzm                ) = v_crm           (vc,1:nx,1:ny,1:nzm                )*YES3D
    w          (1:nx,1:ny,1:nzm                ) = w_crm           (vc,1:nx,1:ny,1:nzm                )
    tabs       (1:nx,1:ny,1:nzm                ) = t_crm           (vc,1:nx,1:ny,1:nzm                )
    micro_field(1:nx,1:ny,1:nzm,1:nmicro_fields) = micro_fields_crm(vc,1:nx,1:ny,1:nzm,1:nmicro_fields)
#ifdef sam1mom
    qn      (1:nx,1:ny,1:nzm) = micro_fields_crm(vc,1:nx,1:ny,1:nzm,3 )
#endif

#ifdef m2005
    cloudliq(1:nx,1:ny,1:nzm) = micro_fields_crm(vc,1:nx,1:ny,1:nzm,11)
#endif

#ifdef m2005
    do k=1, nzm
#ifdef MODAL_AERO
      ! set aerosol data
      l=plev-k+1
      naer (k, 1:ntot_amode) = naermod (vc,l, 1:ntot_amode)
      vaer (k, 1:ntot_amode) = vaerosol(vc,l, 1:ntot_amode)
      hgaer(k, 1:ntot_amode) = hygro   (vc,l, 1:ntot_amode)
#endif
      do j=1, ny
        do i=1, nx
          if(cloudliq(i,j,k).gt.0) then
            if(dopredictNc) then 
              if( micro_field(i,j,k,incl).eq.0) micro_field(i,j,k,incl) = 1.0e6*Nc0/rho(k)
            endif
          endif
        enddo
      enddo
    enddo
#endif 

    w(:,:,nz)=0.
    wsub (:) = 0.      !used in clubb, +++mhwang
    dudt(1:nx,1:ny,1:nzm,1:3) = 0.
    dvdt(1:nx,1:ny,1:nzm,1:3) = 0.
    dwdt(1:nx,1:ny,1:nz,1:3) = 0.
    tke (1:nx,1:ny,1:nzm) = 0.
    tk  (1:nx,1:ny,1:nzm) = 0.
    tkh (1:nx,1:ny,1:nzm) = 0.
    p   (1:nx,1:ny,1:nzm) = 0.

    CF3D(1:nx,1:ny,1:nzm) = 1.

    call micro_init()

    ! initialize sgs fields
    call sgs_init()
        
    do k=1,nzm
      u0(k)=0.
      v0(k)=0.
      t0(k)=0.
      t00(k)=0.
      tabs0(k)=0.
      q0(k)=0.
      qv0(k)=0.
      !+++mhwang these are not initialized ??
      qn0(k) = 0.0
      qp0(k) = 0.0
      tke0(k) = 0.0
      !---mhwang
      do j=1,ny
        do i=1,nx
          t(i,j,k) = tabs(i,j,k)+gamaz(k) &
                    -fac_cond*qcl(i,j,k)-fac_sub*qci(i,j,k) &
                    -fac_cond*qpl(i,j,k)-fac_sub*qpi(i,j,k)
          colprec=colprec+(qpl(i,j,k)+qpi(i,j,k))*pdel(vc,plev-k+1)
          colprecs=colprecs+qpi(i,j,k)*pdel(vc,plev-k+1)
          u0(k)=u0(k)+u(i,j,k)
          v0(k)=v0(k)+v(i,j,k)
          t0(k)=t0(k)+t(i,j,k)
          t00(k)=t00(k)+t(i,j,k)+fac_cond*qpl(i,j,k)+fac_sub*qpi(i,j,k)
          tabs0(k)=tabs0(k)+tabs(i,j,k)
          q0(k)=q0(k)+qv(i,j,k)+qcl(i,j,k)+qci(i,j,k)
          qv0(k) = qv0(k) + qv(i,j,k)
          qn0(k) = qn0(k) + qcl(i,j,k) + qci(i,j,k)
          qp0(k) = qp0(k) + qpl(i,j,k) + qpi(i,j,k)
          tke0(k)=tke0(k)+tke(i,j,k)
        enddo
      enddo

      u0(k) = u0(k) * factor_xy
      v0(k) = v0(k) * factor_xy
      t0(k) = t0(k) * factor_xy
      t00(k) = t00(k) * factor_xy
      tabs0(k) = tabs0(k) * factor_xy
      q0(k) = q0(k) * factor_xy
      qv0(k) = qv0(k) * factor_xy
      qn0(k) = qn0(k) * factor_xy
      qp0(k) = qp0(k) * factor_xy
      tke0(k) = tke0(k) * factor_xy
#ifdef CLUBB_CRM
      ! Update thetav for CLUBB.  This is needed when we have a higher model top 
      ! than is in the sounding, because we subsequently use tv0 to initialize 
      ! thv_ds_zt/zm, which appear in CLUBB's anelastic buoyancy terms. 
      ! -dschanen UWM 11 Feb 2010
      tv0(k) = tabs0(k)*prespot(k)*(1.+epsv*q0(k))
#endif

      l = plev-k+1
      uln(l) = min( umax, max(-umax,ul(vc,l)) )
      vln(l) = min( umax, max(-umax,vl(vc,l)) )*YES3D
      ttend(k) = (tl(vc,l)+gamaz(k)- fac_cond*(qccl(vc,l)+qiil(vc,l))-fac_fus*qiil(vc,l)-t00(k))*idt_gl
      qtend(k) = (ql(vc,l)+qccl(vc,l)+qiil(vc,l)-q0(k))*idt_gl
      utend(k) = (uln(l)-u0(k))*idt_gl
      vtend(k) = (vln(l)-v0(k))*idt_gl
      ug0(k) = uln(l)
      vg0(k) = vln(l)
      tg0(k) = tl(vc,l)+gamaz(k)-fac_cond*qccl(vc,l)-fac_sub*qiil(vc,l)
      qg0(k) = ql(vc,l)+qccl(vc,l)+qiil(vc,l)

    end do ! k

    uhl = u0(1)
    vhl = v0(1)

! estimate roughness length assuming logarithmic profile of velocity near the surface:

    ustar = sqrt(tau00(vc)/rho(1))
    z0 = z0_est(z(1),bflx,wnd,ustar)
    z0 = max(real(0.00001,crm_rknd),min(real(1.,crm_rknd),z0))

    timing_factor = 0.

    prectend=colprec
    precstend=colprecs


#ifdef CLUBB_CRM
    if(doclubb) then
      fluxbu(:, :) = fluxu00(vc)/rhow(1)
      fluxbv(:, :) = fluxv00(vc)/rhow(1)
      fluxbt(:, :) = fluxt00(vc)/rhow(1)
      fluxbq(:, :) = fluxq00(vc)/rhow(1)
    else
      fluxbu(:, :) = 0.
      fluxbv(:, :) = 0.
      fluxbt(:, :) = 0.
      fluxbq(:, :) = 0.
    endif
#else 
    fluxbu=0.
    fluxbv=0.
    fluxbt=0.
    fluxbq=0.
#endif
    fluxtu=0.
    fluxtv=0.
    fluxtt=0.
    fluxtq=0.
    fzero =0.
    precsfc=0.
    precssfc=0.

!---------------------------------------------------
    cld   (vc,:) = 0.
    cldtop(vc,:) = 0.
    gicewp(vc,:) = 0
    gliqwp(vc,:) = 0
    mc    (vc,:) = 0.
    mcup  (vc,:) = 0.
    mcdn  (vc,:) = 0.
    mcuup (vc,:) = 0.
    mcudn (vc,:) = 0.
    crm_qc(vc,:) = 0.
    crm_qi(vc,:) = 0.
    crm_qs(vc,:) = 0.
    crm_qg(vc,:) = 0.
    crm_qr(vc,:) = 0.
#ifdef m2005
    crm_nc(vc,:) = 0.
    crm_ni(vc,:) = 0.
    crm_ns(vc,:) = 0.
    crm_ng(vc,:) = 0.
    crm_nr(vc,:) = 0.
    ! hm 8/31/11 add new variables
    aut_crm_a (vc,:) = 0.
    acc_crm_a (vc,:) = 0.
    evpc_crm_a(vc,:) = 0.
    evpr_crm_a(vc,:) = 0.
    mlt_crm_a (vc,:) = 0.
    sub_crm_a (vc,:) = 0.
    dep_crm_a (vc,:) = 0.
    con_crm_a (vc,:) = 0.

    ! hm 8/31/11 add new output
    ! these are increments added to calculate gcm-grid and time-step avg
    ! note - these values are also averaged over the icycle loop following
    ! the approach for precsfc
    aut1a = 0.
    acc1a = 0.
    evpc1a = 0.
    evpr1a = 0.
    mlt1a = 0.
    sub1a = 0.
    dep1a = 0.
    con1a = 0.
#endif 

    mu_crm (vc,:) = 0.
    md_crm (vc,:) = 0.
    eu_crm (vc,:) = 0.
    du_crm (vc,:) = 0.
    ed_crm (vc,:) = 0.
    dd_crm (vc,:) = 0.
    jt_crm (vc)   = 0.
    mx_crm (vc)   = 0.

    mui_crm(vc,:)   = 0.
    mdi_crm(vc,:)   = 0.

    flux_qt   (vc,:) = 0.
    flux_u    (vc,:) = 0.
    flux_v    (vc,:) = 0.
    fluxsgs_qt(vc,:) = 0.
    tkez      (vc,:) = 0.
    tkesgsz   (vc,:) = 0.
    tkz       (vc,:) = 0.
    flux_qp   (vc,:) = 0.
    pflx      (vc,:) = 0.
    qt_trans  (vc,:) = 0.
    qp_trans  (vc,:) = 0.
    qp_fall   (vc,:) = 0.
    qp_evp    (vc,:) = 0.
    qp_src    (vc,:) = 0.
    qt_ls     (vc,:) = 0.
    t_ls      (vc,:) = 0.

    uwle = 0.
    uwsb = 0.
    vwle = 0.
    vwsb = 0.
    qpsrc = 0.
    qpevp = 0.
    qpfall = 0.
    precflux = 0.

!--------------------------------------------------
#ifdef sam1mom
    if(doprecip) call precip_init()
#endif

    !MRN: Don't want any stochasticity introduced in the standalone.
    !MRN: Need to make sure the first call to crm(...) is not dumped out
    !MRN: Also want to avoid the rabbit hole of dependencies eminating from get_gcol_all_p in phys_grid!
#ifndef CRM_STANDALONE
    call get_gcol_all_p(lchnk, pcols, gcolindex)
    iseed = gcolindex(icol(vc))
    if(u(1,1,1).eq.u(2,1,1).and.u(3,1,2).eq.u(4,1,2)) &
                call setperturb(iseed)
#endif

#ifndef CLUBB_CRM
    !--------------------------
    ! do a CLUBB sanity check
    if ( doclubb .or. doclubbnoninter ) then
      write(0,*) "Cannot call CLUBB if -DCLUBB is not in FFLAGS"
      call endrun('crm main')
    endif
#endif
#ifdef CLUBB_CRM
    !------------------------------------------------------------------
    ! Do initialization for UWM CLUBB
    !------------------------------------------------------------------
    up2       (1:nx, 1:ny, 1:nz ) = clubb_buffer(vc,1:nx, 1:ny, 1:nz ,  1)
    vp2       (1:nx, 1:ny, 1:nz ) = clubb_buffer(vc,1:nx, 1:ny, 1:nz ,  2)
    wprtp     (1:nx, 1:ny, 1:nz ) = clubb_buffer(vc,1:nx, 1:ny, 1:nz ,  3)
    wpthlp    (1:nx, 1:ny, 1:nz ) = clubb_buffer(vc,1:nx, 1:ny, 1:nz ,  4)
    wp2       (1:nx, 1:ny, 1:nz ) = clubb_buffer(vc,1:nx, 1:ny, 1:nz ,  5)
    wp3       (1:nx, 1:ny, 1:nz ) = clubb_buffer(vc,1:nx, 1:ny, 1:nz ,  6)
    rtp2      (1:nx, 1:ny, 1:nz ) = clubb_buffer(vc,1:nx, 1:ny, 1:nz ,  7)
    thlp2     (1:nx, 1:ny, 1:nz ) = clubb_buffer(vc,1:nx, 1:ny, 1:nz ,  8)
    rtpthlp   (1:nx, 1:ny, 1:nz ) = clubb_buffer(vc,1:nx, 1:ny, 1:nz ,  9)
    upwp      (1:nx, 1:ny, 1:nz ) = clubb_buffer(vc,1:nx, 1:ny, 1:nz , 10)
    vpwp      (1:nx, 1:ny, 1:nz ) = clubb_buffer(vc,1:nx, 1:ny, 1:nz , 11)
    cloud_frac(1:nx, 1:ny, 1:nz ) = clubb_buffer(vc,1:nx, 1:ny, 1:nz , 12)
    t_tndcy   (1:nx, 1:ny, 1:nzm) = clubb_buffer(vc,1:nx, 1:ny, 1:nzm, 13)
    qc_tndcy  (1:nx, 1:ny, 1:nzm) = clubb_buffer(vc,1:nx, 1:ny, 1:nzm, 14)
    qv_tndcy  (1:nx, 1:ny, 1:nzm) = clubb_buffer(vc,1:nx, 1:ny, 1:nzm, 15)
    u_tndcy   (1:nx, 1:ny, 1:nzm) = clubb_buffer(vc,1:nx, 1:ny, 1:nzm, 16)
    v_tndcy   (1:nx, 1:ny, 1:nzm) = clubb_buffer(vc,1:nx, 1:ny, 1:nzm, 17)

    ! since no tracer is carried in the current version of MMF, these 
    ! tracer-related restart varialbes are set to zero. +++mhwang, 2011-08
    tracer_tndcy = 0.0
    sclrp2 = 0.0
    sclrprtp = 0.0
    sclrpthlp = 0.0
    wpsclrp =0.0

    if((doclubb.and.docloud).or.(.not.doclubb .and. .not.docloud)) then
      write(0, *) 'doclubb and docloud can not both be true or be false'
      call endrun('crm_clubb2') 
    endif
    if((doclubb_sfc_fluxes.and.docam_sfc_fluxes)) then
      write(0, *) 'doclubb_sfc_fluxes and dosam_sfc_fluxes can not both be true'
      call endrun('crm_clubb_fluxes')
    endif

    if ( doclubb .or. doclubbnoninter ) then
      call clubb_sgs_setup( real( dt*real( nclubb ), kind=time_precision), &
                            latitude, longitude, z, rho, zi, rhow, tv0, tke )
    endif
#endif

#ifdef ECPP
    !ntavg1_ss = dt_gl/3   ! one third of GCM time step, 10 minutes
    ntavg1_ss = min(600._r8, dt_gl)       ! 10 minutes  or the GCM timestep, whichever smaller
          ! ntavg1_ss = number of seconds to average between computing categories.
    ntavg2_ss = dt_gl   ! GCM time step
          ! ntavg2_ss = number of seconds to average between outputs.
          !    This must be a multiple of ntavgt1_ss.

    !  ecpp_crm_init has to be called after ntavg1_ss and ntavg2_ss are set for
    !  their values are used in ecpp_crm_init. 
    call ecpp_crm_init()

    qlsink = 0.0
    qlsink_bf = 0.0
    prain = 0.0
    precr = 0.0
    precsolid = 0.0
#endif

    !+++mhwangtest
    ! test water conservtion problem
    ntotal_step = 0.0
    qtot(vc,:) = 0.0
    qtotmicro(:) = 0.0
    do k=1, nzm
      l=plev-k+1
      do j=1, ny
        do i=1, nx
#ifdef m2005
          qtot(vc,1) = qtot(vc,1)+((micro_field(i,j,k,iqr)+micro_field(i,j,k,iqs)+micro_field(i,j,k,iqg)) * pdel(vc,l)/ggr)/(nx*ny) 
#endif
#ifdef sam1mom
          qtot(vc,1) = qtot(vc,1)+(qpl(i,j,k)+qpi(i,j,k)) * pdel(vc,l)/ggr/(nx*ny)
#endif
        enddo
      enddo
      qtot(vc,1) = qtot(vc,1) + (ql(vc,l)+qccl(vc,l)+qiil(vc,l)) * pdel(vc,l)/ggr
    enddo
    !---mhwangtest

    nstop = dt_gl/dt
    dt = dt_gl/nstop
    nsave3D = nint(60/dt)
    !if(nint(nsave3D*dt).ne.60)then
    !   print *,'CRM: time step=',dt,' is not divisible by 60 seconds'
    !   print *,'this is needed for output every 60 seconds'
    !   stop
    !endif
    nstep = 0
    nprint = 1
    ncycle = 0
    !nrad = nstop/nrad0
    day=day0

    !------------------------------------------------------------------
    !   Main time loop    
    !------------------------------------------------------------------

    do while (nstep.lt.nstop) 
      nstep = nstep + 1
      time = time + dt
      day = day0 + time/86400.
      timing_factor(vc) = timing_factor(vc)+1
      !------------------------------------------------------------------
      !  Check if the dynamical time step should be decreased 
      !  to handle the cases when the flow being locally linearly unstable
      !------------------------------------------------------------------
      ncycle = 1

      call kurant()

      do icyc=1,ncycle
        icycle = icyc
        dtn = dt/ncycle
        dt3(na) = dtn
        dtfactor = dtn/dt

        !---------------------------------------------
        !  	the Adams-Bashforth scheme in time
        call abcoefs()
 
        !---------------------------------------------
        !  	initialize stuff: 
        call zero()

        !-----------------------------------------------------------
        !       Buoyancy term:
        call buoyancy()

        !+++mhwangtest
        ! test water conservtion problem
        ntotal_step = ntotal_step + 1.
        !---mhwangtest 

        !------------------------------------------------------------
        !       Large-scale and surface forcing:
        call forcing()

        do k=1,nzm
          do j=1,ny
            do i=1,nx
              t(i,j,k) = t(i,j,k) + qrad_crm(vc,i,j,k)*dtn
            enddo
          enddo
        enddo

        !----------------------------------------------------------
        !   	suppress turbulence near the upper boundary (spange):
        if (dodamping) call damping()

        !---------------------------------------------------------
        !   Ice fall-out

#ifdef CLUBB_CRM
        if ( docloud .or. doclubb ) then
          call ice_fall()
        endif
#else
        if(docloud) then
            call ice_fall()
        endif
#endif

        !----------------------------------------------------------
        !     Update scalar boundaries after large-scale processes:
        call boundaries(3)

        !---------------------------------------------------------
        !     Update boundaries for velocities:
        call boundaries(0)

        !-----------------------------------------------
        !     surface fluxes:
        if (dosurface) call crmsurface(bflx)

        !-----------------------------------------------------------
        !  SGS physics:
        if (dosgs) call sgs_proc()
        
#ifdef CLUBB_CRM_OLD   
        !----------------------------------------------------------
        ! Do a timestep with CLUBB if enabled:
        ! -dschanen UWM 16 May 2008

        if ( doclubb .or. doclubbnoninter ) then
          ! In case of ice fall, we recompute qci here for the 
          ! single-moment scheme.  Also, subsidence, diffusion and advection have
          ! been applied to micro_field but not qv/qcl so they must be updated.
          call micro_update()
        endif ! doclubb .or. doclubbnoninter

        if ( doclubb ) then
          ! Calculate the vertical integrals for RTM and THLM so we can later
          ! calculate whether CLUBB is a spurious source or sink of either.
          ! - nielsenb UWM 4 Jun 2010
          do i = 1,nx
            do j = 1,ny
              rtm_column = qv(i,j,1:nzm) + qcl(i,j,1:nzm)
              rtm_integral_before(i,j) = vertical_integral( (nz - 2 + 1), rho_ds_zt(2:nz), & 
                                           rtm_column, gr%invrs_dzt(2:nz) )

              thlm_before = t2thetal( t(i,j,1:nzm), gamaz(1:nzm), &
                                   qcl(i,j,1:nzm), qpl(i,j,1:nzm), &
                                   qci(i,j,1:nzm), qpi(i,j,1:nzm), &
                                   prespot(1:nzm) )

              thlm_integral_before(i,j) = vertical_integral( (nz - 2 + 1), rho_ds_zt(2:nz), &
                                                            thlm_before(1:nzm), gr%invrs_dzt(2:nz) )
            enddo
          enddo
          ! End vertical integral
        endif ! doclubb

        if ( doclubb .or. doclubbnoninter ) then
          ! We call CLUBB here because adjustments to the wind
          ! must occur prior to adams() -dschanen 26 Aug 2008
          ! Here we call clubb only if nstep divides the current timestep,
          ! or we're on the very first timestep
          if ( nstep == 1 .or. mod( nstep, nclubb ) == 0 ) then

            call advance_clubb_sgs( real( dtn*real( nclubb ), kind=time_precision), & ! in
                                    real( 0., kind=time_precision ),         & ! in
                                    real( time, kind=time_precision ),       & ! in
                                    rho, rhow, wsub, u, v, w, qpl, qci, qpi, & ! in
                                    t, qv, qcl ) ! in
          endif ! nstep == 1 .or. mod( nstep, nclubb) == 0
        endif ! doclubb .or. doclubbnoninter

#endif
        !----------------------------------------------------------
        !     Fill boundaries for SGS diagnostic fields:
        call boundaries(4)

        !-----------------------------------------------
        !       advection of momentum:
        call advect_mom()

        !----------------------------------------------------------
        !	SGS effects on momentum:

        if(dosgs) call sgs_mom()
#ifdef CLUBB_CRM_OLD
        if ( doclubb ) then
          !          call apply_clubb_sgs_tndcy_mom &
          !               ( dudt, dvdt ) ! in/out
        endif
#endif /*CLUBB_CRM_OLD*/

        !-----------------------------------------------------------
        !       Coriolis force:
        if (docoriolis) call coriolis()

        !---------------------------------------------------------
        !       compute rhs of the Poisson equation and solve it for pressure. 
        call pressure()

        !---------------------------------------------------------
        !       find velocity field at n+1/2 timestep needed for advection of scalars:
        !  Note that at the end of the call, the velocities are in nondimensional form.
        call adams()

        !----------------------------------------------------------
        !     Update boundaries for all prognostic scalar fields for advection:
        call boundaries(2)

        !---------------------------------------------------------
        !      advection of scalars :
        call advect_all_scalars()

        !-----------------------------------------------------------
        !    Convert velocity back from nondimensional form:
        call uvw()

        !----------------------------------------------------------
        !     Update boundaries for scalars to prepare for SGS effects:
        call boundaries(3)

        !---------------------------------------------------------
        !      SGS effects on scalars :
        if (dosgs) call sgs_scalars()

#ifdef CLUBB_CRM_OLD
        ! Re-compute q/qv/qcl based on values computed in CLUBB
        if ( doclubb ) then
          ! Recalculate q, qv, qcl based on new micro_fields (updated by horizontal diffusion)
          call micro_update()

          ! Then Re-compute q/qv/qcl based on values computed in CLUBB
          call apply_clubb_sgs_tndcy_scalars( real( dtn, kind=time_precision), & ! in
                                              t, qv, qcl) ! in/out

          call micro_adjust( qv, qcl ) ! in

          ! Calculate the vertical integrals for RTM and THLM again so
          ! calculate whether CLUBB is a spurious source or sink of either.
          ! - nielsenb UWM 4 Jun 2010
          do i = 1,nx
            do j = 1,ny
              rtm_flux_top = rho_ds_zm(nz) * wprtp(i,j,nz)
              rtm_flux_sfc = rho_ds_zm(1) * fluxbq(i,j)
              rtm_column = qv(i,j,1:nzm) + qcl(i,j,1:nzm)
              rtm_integral_after(i,j) = vertical_integral( (nz - 2 + 1), rho_ds_zt(2:nz), & 
                                            rtm_column, gr%invrs_dzt(2:nz) )
                                             
              rtm_spurious_source(i,j) = calculate_spurious_source( rtm_integral_after(i,j), &
                                                         rtm_integral_before(i,j), &
                                                         rtm_flux_top, rtm_flux_sfc, &
                                                         0.0_core_rknd, real( dtn, kind=core_rknd) )

              thlm_flux_top = rho_ds_zm(nz) * wpthlp(i,j,nz)
              thlm_flux_sfc = rho_ds_zm(1) * fluxbt(i,j)

              thlm_after = t2thetal( t(i,j,1:nzm), gamaz(1:nzm), &
                                     qcl(i,j,1:nzm), qpl(i,j,1:nzm), &
                                     qci(i,j,1:nzm), qpi(i,j,1:nzm), &
                                     prespot(1:nzm) )

              thlm_integral_after(i,j) = vertical_integral( (nz - 2 + 1), rho_ds_zt(2:nz), &
                                                         thlm_after(1:nzm), gr%invrs_dzt(2:nz))
                                            
              thlm_spurious_source(i,j) = calculate_spurious_source( thlm_integral_after(i,j), &
                                                             thlm_integral_before(i,j), &
                                                             thlm_flux_top, thlm_flux_sfc, &
                                                             0.0_core_rknd, real( dtn, kind=core_rknd ))
            enddo
          enddo
          ! End spurious source calculation
        endif  ! doclubb
#endif /*CLUBB_CRM_OLD*/

        !-----------------------------------------------------------
        !       Cloud condensation/evaporation and precipitation processes:
#ifdef CLUBB_CRM
        if(docloud.or.dosmoke.or.doclubb) call micro_proc()
#else
        if(docloud.or.dosmoke) call micro_proc()
#endif /*CLUBB_CRM*/

        !-----------------------------------------------------------
        !    Compute diagnostics fields:
          call diagnose()

        !----------------------------------------------------------
        ! Rotate the dynamic tendency arrays for Adams-bashforth scheme:
        nn=na
        na=nc
        nc=nb
        nb=nn
      enddo ! icycle	
          
      !----------------------------------------------------------
      !----------------------------------------------------------
#ifdef ECPP
      ! Here ecpp_crm_stat is called every CRM time step (dt), not every subcycle time step (dtn). 
      ! This is what the original MMF model did (t_rad, qv_rad, ...). Do we want to call ecpp_crm_stat
      ! every subcycle time step??? +++mhwang
      call ecpp_crm_stat()   
#endif /*ECPP*/

      cwp  = 0.
      cwph = 0.
      cwpm = 0.
      cwpl = 0.

      flag_top(:,:) = .true.

      cltemp = 0.0; cmtemp = 0.0
      chtemp = 0.0; cttemp = 0.0

      do k=1,nzm
        l = plev-k+1
        do j=1,ny
          do i=1,nx
            ! hm modify 9/7/11 for end of timestep, GCM-grid scale hydrometeor output
            ! instead of time-step-averaged
            ! I also modified this for all q and N variables as well as for sam1mom
            ! for consistency
            !hm           crm_qc(l) = crm_qc(l) + qcl(i,j,k)
            !hm           crm_qi(l) = crm_qi(l) + qci(i,j,k)
            !hm           crm_qr(l) = crm_qr(l) + qpl(i,j,k)
            !hm#ifdef sam1mom
            !hm           omg = max(0.,min(1.,(tabs(i,j,k)-tgrmin)*a_gr))
            !hm           crm_qg(l) = crm_qg(l) + qpi(i,j,k)*omg
            !hm           crm_qs(l) = crm_qs(l) + qpi(i,j,k)*(1.-omg)
            !hm#else
            !           crm_qg(l) = crm_qg(l) + qpi(i,j,k)
            !           crm_qs(l) = crm_qs(l) + 0.     ! temporerary solution
            !hm           crm_qg(l) = crm_qg(l) + micro_field(i,j,k,iqg)
            !hm           crm_qs(l) = crm_qs(l) + micro_field(i,j,k,iqs)   
            
            !hm           crm_nc(l) = crm_nc(l) + micro_field(i,j,k,incl)
            !hm           crm_ni(l) = crm_ni(l) + micro_field(i,j,k,inci)
            !hm           crm_nr(l) = crm_nr(l) + micro_field(i,j,k,inr)
            !hm           crm_ng(l) = crm_ng(l) + micro_field(i,j,k,ing)
            !hm           crm_ns(l) = crm_ns(l) + micro_field(i,j,k,ins)
            
            !hm#endif

            tmp1 = rho(nz-k)*adz(nz-k)*dz*(qcl(i,j,nz-k)+qci(i,j,nz-k))
            cwp(i,j) = cwp(i,j)+tmp1
            cttemp(i,j) = max(CF3D(i,j,nz-k), cttemp(i,j))
            if(cwp(i,j).gt.cwp_threshold.and.flag_top(i,j)) then
                cldtop(vc,k) = cldtop(vc,k) + 1
                flag_top(i,j) = .false.
            endif
            if(pres(nz-k).ge.700.) then
                cwpl(i,j) = cwpl(i,j)+tmp1
                cltemp(i,j) = max(CF3D(i,j,nz-k), cltemp(i,j))
            else if(pres(nz-k).lt.400.) then
                cwph(i,j) = cwph(i,j)+tmp1
                chtemp(i,j) = max(CF3D(i,j,nz-k), chtemp(i,j))
            else
                cwpm(i,j) = cwpm(i,j)+tmp1
                cmtemp(i,j) = max(CF3D(i,j,nz-k), cmtemp(i,j))
            endif

            !     qsat = qsatw_crm(tabs(i,j,k),pres(k))
            !     if(qcl(i,j,k)+qci(i,j,k).gt.min(1.e-5,0.01*qsat)) then
            tmp1 = rho(k)*adz(k)*dz
            if(tmp1*(qcl(i,j,k)+qci(i,j,k)).gt.cwp_threshold) then
                 cld(vc,l) = cld(vc,l) + CF3D(i,j,k)
                 if(w(i,j,k+1)+w(i,j,k).gt.2*wmin) then
                   mcup (vc,l) = mcup (vc,l) + rho(k)*0.5*(w(i,j,k+1)+w(i,j,k)) * CF3D(i,j,k)
                   mcuup(vc,l) = mcuup(vc,l) + rho(k)*0.5*(w(i,j,k+1)+w(i,j,k)) * (1.0 - CF3D(i,j,k))
                 endif
                 if(w(i,j,k+1)+w(i,j,k).lt.-2*wmin) then
                   mcdn (vc,l) = mcdn (vc,l) + rho(k)*0.5*(w(i,j,k+1)+w(i,j,k)) * CF3D(i,j,k)
                   mcudn(vc,l) = mcudn(vc,l) + rho(k)*0.5*(w(i,j,k+1)+w(i,j,k)) * (1. - CF3D(i,j,k))
                 endif
            else 
                 if(w(i,j,k+1)+w(i,j,k).gt.2*wmin) then
                   mcuup(vc,l) = mcuup(vc,l) + rho(k)*0.5*(w(i,j,k+1)+w(i,j,k))
                 endif
                 if(w(i,j,k+1)+w(i,j,k).lt.-2*wmin) then
                   mcudn(vc,l) = mcudn(vc,l) + rho(k)*0.5*(w(i,j,k+1)+w(i,j,k))
                 endif
            endif
            
            t_rad  (vc,i,j,k) = t_rad  (vc,i,j,k)+tabs(i,j,k)
            qv_rad (vc,i,j,k) = qv_rad (vc,i,j,k)+max(real(0.,crm_rknd),qv(i,j,k))
            qc_rad (vc,i,j,k) = qc_rad (vc,i,j,k)+qcl(i,j,k)
            qi_rad (vc,i,j,k) = qi_rad (vc,i,j,k)+qci(i,j,k)
            cld_rad(vc,i,j,k) = cld_rad(vc,i,j,k) +  CF3D(i,j,k)
#ifdef m2005
            nc_rad(vc,i,j,k) = nc_rad(vc,i,j,k)+micro_field(i,j,k,incl)
            ni_rad(vc,i,j,k) = ni_rad(vc,i,j,k)+micro_field(i,j,k,inci)
            qs_rad(vc,i,j,k) = qs_rad(vc,i,j,k)+micro_field(i,j,k,iqs)
            ns_rad(vc,i,j,k) = ns_rad(vc,i,j,k)+micro_field(i,j,k,ins)
#endif 
            gliqwp(vc,l) = gliqwp(vc,l)+qcl(i,j,k)
            gicewp(vc,l) = gicewp(vc,l)+qci(i,j,k)
          enddo
        enddo
      enddo

      ! Diagnose mass fluxes to drive CAM's convective transport of tracers. 
      ! definition of mass fluxes is taken from Xu et al., 2002, QJRMS. 
      do k=1, nzm+1
        l=plev+1-k+1
        do j=1, ny
          do i=1, nx
            if(w(i,j,k).gt.0.) then
              kx=max(1, k-1)
              qsat = qsatw_crm(tabs(i,j,kx),pres(kx))
              if(qcl(i,j,kx)+qci(i,j,kx).gt.min(real(1.e-5,crm_rknd),0.01*qsat)) then
                mui_crm(vc,l) = mui_crm(vc,l)+rhow(k)*w(i,j,k)
              endif
            else if (w(i,j,k).lt.0.) then 
              kx=min(k+1, nzm)
              qsat = qsatw_crm(tabs(i,j,kx),pres(kx))
              if(qcl(i,j,kx)+qci(i,j,kx).gt.min(real(1.e-5,crm_rknd),0.01*qsat)) then
                mdi_crm(vc,l) = mdi_crm(vc,l)+rhow(k)*w(i,j,k)
              else if(qpl(i,j,kx)+qpi(i,j,kx).gt.1.0e-4) then 
                mdi_crm(vc,l) = mdi_crm(vc,l)+rhow(k)*w(i,j,k) 
              endif
            endif
          enddo
        enddo
      enddo

      !do k=1,nzm
      ! radlwup0(k)=radlwup0(k)+radlwup(k)
      ! radlwdn0(k)=radlwdn0(k)+radlwdn(k)
      ! radqrlw0(k)=radqrlw0(k)+radqrlw(k)
      ! radswup0(k)=radswup0(k)+radswup(k)
      ! radswdn0(k)=radswdn0(k)+radswdn(k)
      ! radqrsw0(k)=radqrsw0(k)+radqrsw(k)
      !enddo
        
      do j=1,ny
        do i=1,nx
          !           if(cwp (i,j).gt.cwp_threshold) cltot(vc) = cltot(vc) + 1.
          !           if(cwph(i,j).gt.cwp_threshold) clhgh(vc) = clhgh(vc) + 1.
          !           if(cwpm(i,j).gt.cwp_threshold) clmed(vc) = clmed(vc) + 1.
          !           if(cwpl(i,j).gt.cwp_threshold) cllow(vc) = cllow(vc) + 1.
          !  use maxmimum cloud overlap to calcluate cltot, clhgh, 
          !  cldmed, and cldlow   +++ mhwang
          if(cwp (i,j).gt.cwp_threshold) cltot(vc) = cltot(vc) + cttemp(i,j) 
          if(cwph(i,j).gt.cwp_threshold) clhgh(vc) = clhgh(vc) + chtemp(i,j) 
          if(cwpm(i,j).gt.cwp_threshold) clmed(vc) = clmed(vc) + cmtemp(i,j) 
          if(cwpl(i,j).gt.cwp_threshold) cllow(vc) = cllow(vc) + cltemp(i,j) 
        enddo
      enddo

      !        call stepout()
      !----------------------------------------------------------
    enddo ! main loop
    !----------------------------------------------------------

    tmp1 = 1._r8/ dble(nstop)
    t_rad  (vc,:,:,:) = t_rad  (vc,:,:,:) * tmp1
    qv_rad (vc,:,:,:) = qv_rad (vc,:,:,:) * tmp1
    qc_rad (vc,:,:,:) = qc_rad (vc,:,:,:) * tmp1
    qi_rad (vc,:,:,:) = qi_rad (vc,:,:,:) * tmp1
    cld_rad(vc,:,:,:) = cld_rad(vc,:,:,:) * tmp1
#ifdef m2005
    nc_rad(vc,:,:,:) = nc_rad(vc,:,:,:) * tmp1
    ni_rad(vc,:,:,:) = ni_rad(vc,:,:,:) * tmp1
    qs_rad(vc,:,:,:) = qs_rad(vc,:,:,:) * tmp1
    ns_rad(vc,:,:,:) = ns_rad(vc,:,:,:) * tmp1
#endif

    ! no CRM tendencies above its top
    tln(1:ptop-1) = tl(vc,1:ptop-1)
    qln(1:ptop-1) = ql(vc,1:ptop-1)
    qccln(1:ptop-1)= qccl(vc,1:ptop-1)
    qiiln(1:ptop-1)= qiil(vc,1:ptop-1)
    uln(1:ptop-1) = ul(vc,1:ptop-1)
    vln(1:ptop-1) = vl(vc,1:ptop-1)

    !  Compute tendencies due to CRM:
    tln(ptop:plev) = 0.
    qln(ptop:plev) = 0.
    qccln(ptop:plev)= 0.
    qiiln(ptop:plev)= 0.
    uln(ptop:plev) = 0.
    vln(ptop:plev) = 0.
    
    colprec=0
    colprecs=0
    do k = 1,nzm
      l = plev-k+1
      do i=1,nx
        do j=1,ny
          colprec=colprec+(qpl(i,j,k)+qpi(i,j,k))*pdel(vc,plev-k+1)
          colprecs=colprecs+qpi(i,j,k)*pdel(vc,plev-k+1)
          tln(l) = tln(l)+tabs(i,j,k)
          qln(l) = qln(l)+qv(i,j,k)
          qccln(l)= qccln(l)+qcl(i,j,k)
          qiiln(l)= qiiln(l)+qci(i,j,k)
          uln(l) = uln(l)+u(i,j,k)
          vln(l) = vln(l)+v(i,j,k)
        enddo ! k
      enddo
    enddo ! i

    tln(ptop:plev) = tln(ptop:plev) * factor_xy
    qln(ptop:plev) = qln(ptop:plev) * factor_xy
    qccln(ptop:plev) = qccln(ptop:plev) * factor_xy
    qiiln(ptop:plev) = qiiln(ptop:plev) * factor_xy
    uln(ptop:plev) = uln(ptop:plev) * factor_xy
    vln(ptop:plev) = vln(ptop:plev) * factor_xy

#ifdef SPMOMTRANS
    ! whannah - SP CMT tendencies
    ultend(vc,:) = (uln - ul(vc,:))*idt_gl
    vltend(vc,:) = (vln - vl(vc,:))*idt_gl
#endif

    sltend (vc,:) = cp * (tln   - tl  (vc,:)) * idt_gl
    qltend (vc,:) =      (qln   - ql  (vc,:)) * idt_gl
    qcltend(vc,:) =      (qccln - qccl(vc,:)) * idt_gl
    qiltend(vc,:) =      (qiiln - qiil(vc,:)) * idt_gl
    prectend (vc)=(colprec -prectend (vc))/ggr*factor_xy * idt_gl
    precstend(vc)=(colprecs-precstend(vc))/ggr*factor_xy * idt_gl

    ! don't use CRM tendencies from two crm top levels, 
    ! radiation tendencies are added back after the CRM call (see crm_physics_tend)
    sltend (vc,ptop:ptop+1) = 0.
    qltend (vc,ptop:ptop+1) = 0.
    qcltend(vc,ptop:ptop+1) = 0.
    qiltend(vc,ptop:ptop+1) = 0.
    !-------------------------------------------------------------
    ! 
    ! Save the last step to the permanent core:
    u_crm  (vc,1:nx,1:ny,1:nzm) = u   (1:nx,1:ny,1:nzm)
    v_crm  (vc,1:nx,1:ny,1:nzm) = v   (1:nx,1:ny,1:nzm)
    w_crm  (vc,1:nx,1:ny,1:nzm) = w   (1:nx,1:ny,1:nzm)
    t_crm  (vc,1:nx,1:ny,1:nzm) = tabs(1:nx,1:ny,1:nzm)
    micro_fields_crm(vc,1:nx,1:ny,1:nzm,1:nmicro_fields) = micro_field(1:nx,1:ny,1:nzm,1:nmicro_fields)

#ifdef sam1mom
    micro_fields_crm(vc,1:nx,1:ny,1:nzm,3) = qn(1:nx,1:ny,1:nzm)
#endif
#ifdef m2005
    micro_fields_crm(vc,1:nx,1:ny,1:nzm,11) = cloudliq(1:nx,1:ny,1:nzm)
#endif
    crm_tk   (vc,1:nx,1:ny,1:nzm) = tk  (1:nx, 1:ny, 1:nzm)
    crm_tkh  (vc,1:nx,1:ny,1:nzm) = tkh (1:nx, 1:ny, 1:nzm)
    cld3d_crm(vc,1:nx,1:ny,1:nzm) = CF3D(1:nx, 1:ny, 1:nzm)
#ifdef CLUBB_CRM
    clubb_buffer(vc,1:nx, 1:ny, 1:nz ,  1) = up2       (1:nx, 1:ny, 1:nz )
    clubb_buffer(vc,1:nx, 1:ny, 1:nz ,  2) = vp2       (1:nx, 1:ny, 1:nz )
    clubb_buffer(vc,1:nx, 1:ny, 1:nz ,  3) = wprtp     (1:nx, 1:ny, 1:nz )
    clubb_buffer(vc,1:nx, 1:ny, 1:nz ,  4) = wpthlp    (1:nx, 1:ny, 1:nz )
    clubb_buffer(vc,1:nx, 1:ny, 1:nz ,  5) = wp2       (1:nx, 1:ny, 1:nz )
    clubb_buffer(vc,1:nx, 1:ny, 1:nz ,  6) = wp3       (1:nx, 1:ny, 1:nz )
    clubb_buffer(vc,1:nx, 1:ny, 1:nz ,  7) = rtp2      (1:nx, 1:ny, 1:nz )
    clubb_buffer(vc,1:nx, 1:ny, 1:nz ,  8) = thlp2     (1:nx, 1:ny, 1:nz )
    clubb_buffer(vc,1:nx, 1:ny, 1:nz ,  9) = rtpthlp   (1:nx, 1:ny, 1:nz )
    clubb_buffer(vc,1:nx, 1:ny, 1:nz , 10) = upwp      (1:nx, 1:ny, 1:nz )
    clubb_buffer(vc,1:nx, 1:ny, 1:nz , 11) = vpwp      (1:nx, 1:ny, 1:nz )
    clubb_buffer(vc,1:nx, 1:ny, 1:nz , 12) = cloud_frac(1:nx, 1:ny, 1:nz )
    clubb_buffer(vc,1:nx, 1:ny, 1:nzm, 13) = t_tndcy   (1:nx, 1:ny, 1:nzm)
    clubb_buffer(vc,1:nx, 1:ny, 1:nzm, 14) = qc_tndcy  (1:nx, 1:ny, 1:nzm)
    clubb_buffer(vc,1:nx, 1:ny, 1:nzm, 15) = qv_tndcy  (1:nx, 1:ny, 1:nzm)
    clubb_buffer(vc,1:nx, 1:ny, 1:nzm, 16) = u_tndcy   (1:nx, 1:ny, 1:nzm)
    clubb_buffer(vc,1:nx, 1:ny, 1:nzm, 17) = v_tndcy   (1:nx, 1:ny, 1:nzm)

    crm_cld    (vc,1:nx, 1:ny, 1:nz ) = cloud_frac  (1:nx, 1:ny, 1:nz )
    clubb_tk   (vc,1:nx, 1:ny, 1:nzm) = tk_clubb    (1:nx, 1:ny, 1:nzm)
    clubb_tkh  (vc,1:nx, 1:ny, 1:nzm) = tkh_clubb   (1:nx, 1:ny, 1:nzm)
    relvar     (vc,1:nx, 1:ny, 1:nzm) = relvarg     (1:nx, 1:ny, 1:nzm)
    accre_enhan(vc,1:nx, 1:ny, 1:nzm) = accre_enhang(1:nx, 1:ny, 1:nzm) 
    qclvar     (vc,1:nx, 1:ny, 1:nzm) = qclvarg     (1:nx, 1:ny, 1:nzm)
#endif

    do k=1,nzm
     do j=1,ny
      do i=1,nx
        qc_crm (vc,i,j,k) = qcl(i,j,k)
        qi_crm (vc,i,j,k) = qci(i,j,k)
        qpc_crm(vc,i,j,k) = qpl(i,j,k)
        qpi_crm(vc,i,j,k) = qpi(i,j,k)
#ifdef m2005
        wvar_crm(vc,i,j,k) = wvar (i,j,k)
        aut_crm (vc,i,j,k) = aut1 (i,j,k)
        acc_crm (vc,i,j,k) = acc1 (i,j,k)
        evpc_crm(vc,i,j,k) = evpc1(i,j,k)
        evpr_crm(vc,i,j,k) = evpr1(i,j,k)
        mlt_crm (vc,i,j,k) = mlt1 (i,j,k)
        sub_crm (vc,i,j,k) = sub1 (i,j,k)
        dep_crm (vc,i,j,k) = dep1 (i,j,k)
        con_crm (vc,i,j,k) = con1 (i,j,k)
#endif
        enddo
      enddo
    enddo
    z0m     (vc) = z0 
    taux_crm(vc) = taux0 / dble(nstop)
    tauy_crm(vc) = tauy0 / dble(nstop)

    !---------------------------------------------------------------
    !  Diagnostics:

    ! hm add 9/7/11, change from GCM-time step avg to end-of-timestep
    do k=1,nzm
      l = plev-k+1
      do j=1,ny
        do i=1,nx
          crm_qc(vc,l) = crm_qc(vc,l) + qcl(i,j,k)
          crm_qi(vc,l) = crm_qi(vc,l) + qci(i,j,k)
          crm_qr(vc,l) = crm_qr(vc,l) + qpl(i,j,k)
#ifdef sam1mom
          omg = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(i,j,k)-tgrmin)*a_gr))
          crm_qg(vc,l) = crm_qg(vc,l) + qpi(i,j,k)*omg
          crm_qs(vc,l) = crm_qs(vc,l) + qpi(i,j,k)*(1.-omg)
#else
          !crm_qg(vc,l) = crm_qg(vc,l) + qpi(i,j,k)
          !crm_qs(vc,l) = crm_qs(vc,l) + 0.     ! temporerary solution
          crm_qg(vc,l) = crm_qg(vc,l) + micro_field(i,j,k,iqg)
          crm_qs(vc,l) = crm_qs(vc,l) + micro_field(i,j,k,iqs)

          crm_nc(vc,l) = crm_nc(vc,l) + micro_field(i,j,k,incl)
          crm_ni(vc,l) = crm_ni(vc,l) + micro_field(i,j,k,inci)
          crm_nr(vc,l) = crm_nr(vc,l) + micro_field(i,j,k,inr)
          crm_ng(vc,l) = crm_ng(vc,l) + micro_field(i,j,k,ing)
          crm_ns(vc,l) = crm_ns(vc,l) + micro_field(i,j,k,ins)
#endif
        enddo
      enddo
    enddo

    cld   (vc,:) = min(1._r8,cld   (vc,:)/real(nstop,crm_rknd)*factor_xy)
    cldtop(vc,:) = min(1._r8,cldtop(vc,:)/real(nstop,crm_rknd)*factor_xy)
    gicewp(vc,:) = gicewp(vc,:)*pdel(vc,:)*1000./ggr/real(nstop,crm_rknd)*factor_xy
    gliqwp(vc,:) = gliqwp(vc,:)*pdel(vc,:)*1000./ggr/real(nstop,crm_rknd)*factor_xy
    mcup  (vc,:) = mcup (vc,:) / real(nstop,crm_rknd) * factor_xy
    mcdn  (vc,:) = mcdn (vc,:) / real(nstop,crm_rknd) * factor_xy
    mcuup (vc,:) = mcuup(vc,:) / real(nstop,crm_rknd) * factor_xy
    mcudn (vc,:) = mcudn(vc,:) / real(nstop,crm_rknd) * factor_xy
    mc    (vc,:) = mcup(vc,:) + mcdn(vc,:) + mcuup(vc,:) + mcudn(vc,:)
    ! hm 9/7/11 modify for end-of-timestep instead of timestep-avg output
    !hm        crm_qc = crm_qc / float(nstop) * factor_xy
    !hm        crm_qi = crm_qi / float(nstop) * factor_xy
    !hm        crm_qs = crm_qs / float(nstop) * factor_xy
    !hm        crm_qg = crm_qg / float(nstop) * factor_xy
    !hm        crm_qr = crm_qr / float(nstop) * factor_xy
    !hm#ifdef m2005
    !hm        crm_nc = crm_nc / float(nstop) * factor_xy
    !hm        crm_ni = crm_ni / float(nstop) * factor_xy
    !hm        crm_ns = crm_ns / float(nstop) * factor_xy
    !hm        crm_ng = crm_ng / float(nstop) * factor_xy
    !hm        crm_nr = crm_nr / float(nstop) * factor_xy

    crm_qc(vc,:) = crm_qc(vc,:) * factor_xy
    crm_qi(vc,:) = crm_qi(vc,:) * factor_xy
    crm_qs(vc,:) = crm_qs(vc,:) * factor_xy
    crm_qg(vc,:) = crm_qg(vc,:) * factor_xy
    crm_qr(vc,:) = crm_qr(vc,:) * factor_xy
#ifdef m2005
    crm_nc(vc,:) = crm_nc(vc,:) * factor_xy
    crm_ni(vc,:) = crm_ni(vc,:) * factor_xy
    crm_ns(vc,:) = crm_ns(vc,:) * factor_xy
    crm_ng(vc,:) = crm_ng(vc,:) * factor_xy
    crm_nr(vc,:) = crm_nr(vc,:) * factor_xy

    ! hm 8/31/11 new output, gcm-grid- and time-step avg
    ! add loop over i,j do get horizontal avg, and flip vertical array
    do k=1,nzm
      l = plev-k+1
      do j=1,ny
        do i=1,nx
          aut_crm_a (vc,l) = aut_crm_a (vc,l) + aut1a (i,j,k)
          acc_crm_a (vc,l) = acc_crm_a (vc,l) + acc1a (i,j,k)
          evpc_crm_a(vc,l) = evpc_crm_a(vc,l) + evpc1a(i,j,k)
          evpr_crm_a(vc,l) = evpr_crm_a(vc,l) + evpr1a(i,j,k)
          mlt_crm_a (vc,l) = mlt_crm_a (vc,l) + mlt1a (i,j,k)
          sub_crm_a (vc,l) = sub_crm_a (vc,l) + sub1a (i,j,k)
          dep_crm_a (vc,l) = dep_crm_a (vc,l) + dep1a (i,j,k)
          con_crm_a (vc,l) = con_crm_a (vc,l) + con1a (i,j,k)
        enddo
      enddo
    enddo

    ! note, rates are divded by dt to get mean rate over step
    aut_crm_a (vc,:) = aut_crm_a (vc,:) / dble(nstop) * factor_xy / dt
    acc_crm_a (vc,:) = acc_crm_a (vc,:) / dble(nstop) * factor_xy / dt
    evpc_crm_a(vc,:) = evpc_crm_a(vc,:) / dble(nstop) * factor_xy / dt
    evpr_crm_a(vc,:) = evpr_crm_a(vc,:) / dble(nstop) * factor_xy / dt
    mlt_crm_a (vc,:) = mlt_crm_a (vc,:) / dble(nstop) * factor_xy / dt
    sub_crm_a (vc,:) = sub_crm_a (vc,:) / dble(nstop) * factor_xy / dt
    dep_crm_a (vc,:) = dep_crm_a (vc,:) / dble(nstop) * factor_xy / dt
    con_crm_a (vc,:) = con_crm_a (vc,:) / dble(nstop) * factor_xy / dt

#endif
    precc (vc) = 0.
    precl (vc) = 0.
    precsc(vc) = 0.
    precsl(vc) = 0.
    do j=1,ny
      do i=1,nx
#ifdef sam1mom
        precsfc(i,j) = precsfc(i,j)*dz/dt/dble(nstop)
        precssfc(i,j) = precssfc(i,j)*dz/dt/dble(nstop)
#endif
#ifdef m2005 
        ! precsfc and precssfc from the subroutine of micro_proc in M2005 have a unit mm/s/dz
        !          precsfc(i,j) = precsfc(i,j)*dz/dble(nstop)     !mm/s/dz --> mm/s
        !          precssfc(i,j) = precssfc(i,j)*dz/dble(nstop)   !mm/s/dz --> mm/s
        ! precsfc and precssfc from the subroutine of micro_proc in M2005 have a unit mm/dz
        precsfc(i,j) = precsfc(i,j)*dz/dt/dble(nstop)     !mm/s/dz --> mm/s
        precssfc(i,j) = precssfc(i,j)*dz/dt/dble(nstop)   !mm/s/dz --> mm/s
#endif
        if(precsfc(i,j).gt.10./86400.) then
           precc (vc) = precc (vc) + precsfc(i,j)
           precsc(vc) = precsc(vc) + precssfc(i,j)
        else
           precl (vc) = precl (vc) + precsfc(i,j)
           precsl(vc) = precsl(vc) + precssfc(i,j)
        endif
      enddo
    enddo
    prec_crm(vc,:,:) = precsfc/1000.           !mm/s --> m/s
    precc   (vc)     = precc (vc)*factor_xy/1000.     
    precl   (vc)     = precl (vc)*factor_xy/1000.
    precsc  (vc)     = precsc(vc)*factor_xy/1000.
    precsl  (vc)     = precsl(vc)*factor_xy/1000.

    !+++mhwangtest
    ! test water conservtion problem
    do k=1, nzm
      l=plev-k+1
      do j=1, ny
        do i=1, nx
#ifdef m2005
          qtot(vc,9) = qtot(vc,9)+((micro_field(i,j,k,iqr)+micro_field(i,j,k,iqs)+micro_field(i,j,k,iqg)) * pdel(vc,l)/ggr)/(nx*ny)
          qtot(vc,9) = qtot(vc,9)+((micro_field(i,j,k,iqv)+micro_field(i,j,k,iqci)) * pdel(vc,l)/ggr)/(nx*ny)
#endif
#ifdef sam1mom
          qtot(vc,9) = qtot(vc,9)+((micro_field(i,j,k,1)+micro_field(i,j,k,2)) * pdel(vc,l)/ggr)/(nx*ny)
#endif
        enddo
      enddo
    enddo
    qtot(vc,9) = qtot(vc,9) + (precc(vc)+precl(vc))*1000 * dt_gl

    cltot(vc) = cltot(vc) *factor_xy/nstop
    clhgh(vc) = clhgh(vc) *factor_xy/nstop
    clmed(vc) = clmed(vc) *factor_xy/nstop
    cllow(vc) = cllow(vc) *factor_xy/nstop

    jt_crm(vc) = plev * 1.0
    mx_crm(vc) = 1.0
    do k=1, plev 
      mu_crm(vc,k)=0.5*(mui_crm(vc,k)+mui_crm(vc,k+1))
      md_crm(vc,k)=0.5*(mdi_crm(vc,k)+mdi_crm(vc,k+1))
      mu_crm(vc,k)=mu_crm(vc,k)*ggr/100.          !kg/m2/s --> mb/s
      md_crm(vc,k)=md_crm(vc,k)*ggr/100.          !kg/m2/s --> mb/s
      eu_crm(vc,k) = 0.
      if(mui_crm(vc,k)-mui_crm(vc,k+1).gt.0) then
        eu_crm(vc,k)=(mui_crm(vc,k)-mui_crm(vc,k+1))*ggr/pdel(vc,k)    !/s
      else
        du_crm(vc,k)=-1.0*(mui_crm(vc,k)-mui_crm(vc,k+1))*ggr/pdel(vc,k)   !/s
      endif
      if(mdi_crm(vc,k+1)-mdi_crm(vc,k).lt.0) then
        ed_crm(vc,k)=(mdi_crm(vc,k)-mdi_crm(vc,k+1))*ggr/pdel(vc,k) ! /s
      else
        dd_crm(vc,k)=-1.*(mdi_crm(vc,k)-mdi_crm(vc,k+1))*ggr/pdel(vc,k)   !/s 
      endif
      if(abs(mu_crm(vc,k)).gt.1.0e-15.or.abs(md_crm(vc,k)).gt.1.0e-15) then
        jt_crm(vc) = min(k*1.0_r8, jt_crm(vc))
        mx_crm(vc) = max(k*1.0_r8, mx_crm(vc))
      endif
    enddo
        
    !-------------------------------------------------------------
    !       Fluxes and other stat:
    !-------------------------------------------------------------
    do k=1,nzm
      u2z = 0.
      v2z = 0.
      w2z = 0.
      do j=1,ny
        do i=1,nx
          u2z = u2z+(u(i,j,k)-u0(k))**2
          v2z = v2z+(v(i,j,k)-v0(k))**2
          w2z = w2z+0.5*(w(i,j,k+1)**2+w(i,j,k)**2)
        enddo
      enddo
      !+++mhwang
      ! mkwsb, mkle, mkadv, mkdiff (also flux_u, flux_v) seem not calculted correclty in the spcam3.5 codes. 
      ! Only values at the last time step are calculated, but is averaged over the entire GCM 
      ! time step. 
      !---mhwang

      tmp1 = dz/rhow(k)
      tmp2 = tmp1/dtn                        ! dtn is calculated inside of the icyc loop. 
                                             ! It seems wrong to use it here ???? +++mhwang
      mkwsb (k,:) = mkwsb (k,:) * tmp1*rhow(k) * factor_xy/nstop     !kg/m3/s --> kg/m2/s
      mkwle (k,:) = mkwle (k,:) * tmp2*rhow(k) * factor_xy/nstop     !kg/m3   --> kg/m2/s  
      mkadv (k,:) = mkadv (k,:) * factor_xy*idt_gl     ! kg/kg  --> kg/kg/s
      mkdiff(k,:) = mkdiff(k,:) * factor_xy*idt_gl   ! kg/kg  --> kg/kg/s

      ! qpsrc, qpevp, qpfall in M2005 are calculated in micro_flux. 
      qpsrc   (k) = qpsrc   (k) * factor_xy*idt_gl
      qpevp   (k) = qpevp   (k) * factor_xy*idt_gl
      qpfall  (k) = qpfall  (k) * factor_xy*idt_gl   ! kg/kg in M2005 ---> kg/kg/s 
      precflux(k) = precflux(k) * factor_xy*dz/dt/nstop  !kg/m2/dz in M2005 -->kg/m2/s or mm/s (idt_gl=1/dt/nstop)

      l = plev-k+1
      flux_u    (vc,l) = (uwle(k) + uwsb(k))*tmp1*factor_xy/nstop
      flux_v    (vc,l) = (vwle(k) + vwsb(k))*tmp1*factor_xy/nstop
#ifdef sam1mom
      flux_qt   (vc,l) = mkwle(k,1) + mkwsb(k,1)
      fluxsgs_qt(vc,l) = mkwsb(k,1)
      flux_qp   (vc,l) = mkwle(k,2) + mkwsb(k,2)
      qt_trans  (vc,l) = mkadv(k,1) + mkdiff(k,1)
      qp_trans  (vc,l) = mkadv(k,2) + mkdiff(k,2)
#endif
#ifdef m2005
      flux_qt   (vc,l) = mkwle(k,1   ) + mkwsb(k,1   ) +  &
                         mkwle(k,iqci) + mkwsb(k,iqci)
      fluxsgs_qt(vc,l) = mkwsb(k,1   ) + mkwsb(k,iqci)
      flux_qp   (vc,l) = mkwle(k,iqr) + mkwsb(k,iqr) +  &
                         mkwle(k,iqs) + mkwsb(k,iqs) + mkwle(k,iqg) + mkwsb(k,iqg)
      qt_trans  (vc,l) = mkadv (k,1) + mkadv (k,iqci) + &
                         mkdiff(k,1) + mkdiff(k,iqci) 
      qp_trans  (vc,l) = mkadv (k,iqr) + mkadv (k,iqs) + mkadv (k,iqg) + &
                         mkdiff(k,iqr) + mkdiff(k,iqs) + mkdiff(k,iqg) 
#endif
      tkesgsz   (vc,l)= rho(k)*sum(tke(1:nx,1:ny,k))*factor_xy
      tkez      (vc,l)= rho(k)*0.5*(u2z+v2z*YES3D+w2z)*factor_xy + tkesgsz(vc,l)
      tkz       (vc,l) = sum(tk(1:nx, 1:ny, k)) * factor_xy
      pflx      (vc,l) = precflux(k)/1000.       !mm/s  -->m/s

      qp_fall   (vc,l) = qpfall(k)
      qp_evp    (vc,l) = qpevp(k)
      qp_src    (vc,l) = qpsrc(k)

      qt_ls     (vc,l) = qtend(k)
      t_ls      (vc,l) = ttend(k)
    enddo

#ifdef ECPP
    abnd         (vc,:,:,:,:)=0.0
    abnd_tf      (vc,:,:,:,:)=0.0
    massflxbnd   (vc,:,:,:,:)=0.0
    acen         (vc,:,:,:,:)=0.0
    acen_tf      (vc,:,:,:,:)=0.0
    rhcen        (vc,:,:,:,:)=0.0
    qcloudcen    (vc,:,:,:,:)=0.0
    qicecen      (vc,:,:,:,:)=0.0
    qlsinkcen    (vc,:,:,:,:)=0.0
    precrcen     (vc,:,:,:,:)=0.0
    precsolidcen (vc,:,:,:,:)=0.0
    qlsink_bfcen (vc,:,:,:,:)=0.0
    qlsink_avgcen(vc,:,:,:,:)=0.0
    praincen     (vc,:,:,:,:)=0.0

    wupthresh_bnd   (vc,:)=0.0
    wdownthresh_bnd (vc,:)=0.0
    wwqui_cen       (vc,:)=0.0
    wwqui_bnd       (vc,:)=0.0
    wwqui_cloudy_cen(vc,:)=0.0
    wwqui_cloudy_bnd(vc,:)=0.0

    ! default is clear, non-precipitating, and quiescent class
    abnd   (vc,:,1,1,1)=1.0 
    abnd_tf(vc,:,1,1,1)=1.0
    acen   (vc,:,1,1,1)=1.0
    acen_tf(vc,:,1,1,1)=1.0
    do k=1, nzm
      l=plev-k+1
      acen            (vc,l,:,:,:) = area_cen_sum        (k,:,1:ncls_ecpp_in,:)
      acen_tf         (vc,l,:,:,:) = area_cen_final      (k,:,1:ncls_ecpp_in,:)
      rhcen           (vc,l,:,:,:) = rh_cen_sum          (k,:,1:ncls_ecpp_in,:)
      qcloudcen       (vc,l,:,:,:) = qcloud_cen_sum      (k,:,1:ncls_ecpp_in,:)
      qicecen         (vc,l,:,:,:) = qice_cen_sum        (k,:,1:ncls_ecpp_in,:)
      qlsinkcen       (vc,l,:,:,:) = qlsink_cen_sum      (k,:,1:ncls_ecpp_in,:)
      precrcen        (vc,l,:,:,:) = precr_cen_sum       (k,:,1:ncls_ecpp_in,:)
      precsolidcen    (vc,l,:,:,:) = precsolid_cen_sum   (k,:,1:ncls_ecpp_in,:)
      wwqui_cen       (vc,l)       = wwqui_cen_sum       (k)
      wwqui_cloudy_cen(vc,l)       = wwqui_cloudy_cen_sum(k)
      qlsink_bfcen    (vc,l,:,:,:) = qlsink_bf_cen_sum   (k,:,1:ncls_ecpp_in,:)
      qlsink_avgcen   (vc,l,:,:,:) = qlsink_avg_cen_sum  (k,:,1:ncls_ecpp_in,:)
      praincen        (vc,l,:,:,:) = prain_cen_sum       (k,:,1:ncls_ecpp_in,:)
    enddo
    do k=1, nzm+1
      l=plev+1-k+1
      abnd            (vc,l,:,:,:) = area_bnd_sum        (k,:,1:ncls_ecpp_in,:)
      abnd_tf         (vc,l,:,:,:) = area_bnd_final      (k,:,1:ncls_ecpp_in,:)
      massflxbnd      (vc,l,:,:,:) = mass_bnd_sum        (k,:,1:ncls_ecpp_in,:)
      wupthresh_bnd   (vc,l)       = wup_thresh          (k)
      wdownthresh_bnd (vc,l)       = wdown_thresh        (k)
      wwqui_bnd       (vc,l)       = wwqui_bnd_sum       (k)
      wwqui_cloudy_bnd(vc,l)       = wwqui_cloudy_bnd_sum(k)
    enddo
#endif
        
    timing_factor(vc) = timing_factor(vc) / nstop

#ifdef CLUBB_CRM
    ! Deallocate CLUBB variables, etc.
    ! -UWM
    if ( doclubb .or. doclubbnoninter ) call clubb_sgs_cleanup( )
#endif
#ifdef ECPP
    ! Deallocate ECPP variables
    call ecpp_crm_cleanup ()
#endif

    call crm_dump_output( igstep,plev,crm_tk(vc,:,:,:),crm_tkh(vc,:,:,:),cltot(vc),clhgh(vc),clmed(vc),cllow(vc),sltend(vc,:),u_crm(vc,:,:,:),v_crm(vc,:,:,:),&
                          w_crm(vc,:,:,:),t_crm(vc,:,:,:),micro_fields_crm(vc,:,:,:,:),qltend(vc,:),qcltend(vc,:),qiltend(vc,:),t_rad(vc,:,:,:),qv_rad(vc,:,:,:),&
                          qc_rad(vc,:,:,:),qi_rad(vc,:,:,:),cld_rad(vc,:,:,:),cld3d_crm(vc,:,:,:), &
#ifdef CLUBB_CRM
                          clubb_buffer(vc,:,:,:,:),crm_cld(vc,:,:,:),clubb_tk(vc,:,:,:),clubb_tkh(vc,:,:,:),relvar(vc,:,:,:),accre_enhan(vc,:,:,:),qclvar(vc,:,:,:) , &
#endif
#ifdef CRM3D
                          ultend(vc,:),vltend(vc,:) , &
#endif
#ifdef m2005
                          nc_rad(vc,:,:,:),ni_rad(vc,:,:,:),qs_rad(vc,:,:,:),ns_rad(vc,:,:,:),wvar_crm(vc,:,:,:),aut_crm(vc,:,:,:),acc_crm(vc,:,:,:),evpc_crm(vc,:,:,:), &
                          evpr_crm(vc,:,:,:),mlt_crm(vc,:,:,:),sub_crm(vc,:,:,:),dep_crm(vc,:,:,:),con_crm(vc,:,:,:),aut_crm_a(vc,:),acc_crm_a(vc,:),evpc_crm_a(vc,:), &
                          evpr_crm_a(vc,:),mlt_crm_a(vc,:),sub_crm_a(vc,:),dep_crm_a(vc,:),con_crm_a(vc,:),crm_nc(vc,:),crm_ni(vc,:),crm_ns(vc,:),crm_ng(vc,:),crm_nr(vc,:), &
#endif
#ifdef ECPP
                          acen(vc,:,:,:,:),acen_tf(vc,:,:,:,:),rhcen(vc,:,:,:,:),qcloudcen(vc,:,:,:,:),qicecen(vc,:,:,:,:),qlsinkcen(vc,:,:,:,:),precrcen(vc,:,:,:,:),&
                          precsolidcen(vc,:,:,:,:),qlsink_bfcen(vc,:,:,:,:),qlsink_avgcen(vc,:,:,:,:),praincen(vc,:,:,:,:),wwqui_cen(vc,:),wwqui_cloudy_cen(vc,:), &
                          abnd(vc,:,:,:,:),abnd_tf(vc,:,:,:,:),massflxbnd(vc,:,:,:,:),wupthresh_bnd(vc,:),wdownthresh_bnd(vc,:),wwqui_bnd(vc,:),wwqui_cloudy_bnd(vc,:), &
#endif
                          precc(vc),precl(vc),cld(vc,:),cldtop(vc,:),gicewp(vc,:),gliqwp(vc,:),mc(vc,:),mcup(vc,:),mcdn(vc,:),mcuup(vc,:),mcudn(vc,:),crm_qc(vc,:), &
                          crm_qi(vc,:),crm_qs(vc,:),crm_qg(vc,:),crm_qr(vc,:),mu_crm(vc,:),md_crm(vc,:),du_crm(vc,:),eu_crm(vc,:),ed_crm(vc,:),dd_crm(vc,:),jt_crm(vc), &
                          mx_crm(vc),mui_crm(vc,:),mdi_crm(vc,:),flux_qt(vc,:),fluxsgs_qt(vc,:),tkez(vc,:),tkesgsz(vc,:),tkz(vc,:),flux_u(vc,:),flux_v(vc,:),flux_qp(vc,:), &
                          pflx(vc,:),qt_ls(vc,:),qt_trans(vc,:),qp_trans(vc,:),qp_fall(vc,:),qp_src(vc,:),qp_evp(vc,:),t_ls(vc,:),prectend(vc),precstend(vc),precsc(vc), &
                          precsl(vc),taux_crm(vc),tauy_crm(vc),z0m(vc),timing_factor(vc),qc_crm(vc,:,:,:),qi_crm(vc,:,:,:),qpc_crm(vc,:,:,:),qpi_crm(vc,:,:,:), &
                          prec_crm(vc,:,:),qtot(vc,:) )
  enddo

  end subroutine crm


end module crm_module
