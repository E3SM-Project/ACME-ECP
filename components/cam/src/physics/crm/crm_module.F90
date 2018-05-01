
module crm_module
  use openacc_utils
  use task_init_mod, only: task_init
  use abcoefs_mod, only: abcoefs
  use kurant_mod, only: kurant
  use setperturb_mod, only: setperturb
  use boundaries_mod, only: boundaries
  use forcing_mod, only: forcing
  use advect_mom_mod, only: advect_mom
  use adams_mod, only: adams
  use advect_all_scalars_mod, only: advect_all_scalars
  use sat_mod
  use crmsurface_mod
  use perf_mod
#ifdef sam1mom
  use precip_init_mod
  use micro_params
#endif
  use zero_mod
  use buoyancy_mod
  use pressure_mod
  use uvw_mod
  use diagnose_mod
  use damping_mod
  use ice_fall_mod
  use coriolis_mod
  !---------------------------------------------------------------
  !  Super-parameterization's main driver
  !  Marat Khairoutdinov, 2001-2009
  !---------------------------------------------------------------
  use setparm_mod, only : setparm

contains

! subroutine crm  (lchnk, icol, &
subroutine crm(lchnk, icol, ncrms, is_first_step , &
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
    use phys_grid             , only: get_rlon_p, get_rlat_p, get_gcol_p  !, get_gcol_all_p
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
  use crmdims               , only: nclubbvars, crm_nx_rad, crm_ny_rad
#ifdef CLUBB_CRM
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
    integer , intent(in   ) :: ncrms                            ! Number of crm instances
    logical , intent(in   ) :: is_first_step                    ! flag to indicate first CRM integration - used to call setperturb()
    integer , intent(in   ) :: plev                             ! number of levels in parent model
    real(r8), intent(in   ) :: dt_gl                            ! global model's time step
    integer , intent(in   ) :: icol              (ncrms)     ! column identifier (only for lat/lon and random seed)
#ifdef CRM_STANDALONE
  real(crm_rknd)   , intent(in) :: latitude0_in  (ncrms)
  real(crm_rknd)   , intent(in) :: longitude0_in (ncrms)
#endif
  real(r8), intent(in   ) :: ps                  (ncrms)       ! Global grid surface pressure (Pa)
  real(r8), intent(in   ) :: pmid                (ncrms,plev)  ! Global grid pressure (Pa)
  real(r8), intent(in   ) :: pdel                (ncrms,plev)  ! Layer's pressure thickness (Pa)
  real(r8), intent(in   ) :: phis                (ncrms)       ! Global grid surface geopotential (m2/s2)
  real(r8), intent(in   ) :: zmid                (ncrms,plev)  ! Global grid height (m)
  real(r8), intent(in   ) :: zint                (ncrms,plev+1)! Global grid interface height (m)
  real(r8), intent(in   ) :: qrad_crm            (ncrms,crm_nx_rad, crm_ny_rad, crm_nz) ! CRM rad. heating
  real(r8), intent(in   ) :: ocnfrac             (ncrms)       ! area fraction of the ocean
  real(r8), intent(in   ) :: tau00               (ncrms)       ! large-scale surface stress (N/m2)
  real(r8), intent(in   ) :: wndls               (ncrms)       ! large-scale surface wind (m/s)
  real(r8), intent(in   ) :: bflxls              (ncrms)       ! large-scale surface buoyancy flux (K m/s)
  real(r8), intent(in   ) :: fluxu00             (ncrms)       ! surface momenent fluxes [N/m2]
  real(r8), intent(in   ) :: fluxv00             (ncrms)       ! surface momenent fluxes [N/m2]
  real(r8), intent(in   ) :: fluxt00             (ncrms)       ! surface sensible heat fluxes [K Kg/ (m2 s)]
  real(r8), intent(in   ) :: fluxq00             (ncrms)       ! surface latent heat fluxes [ kg/(m2 s)]
  real(r8), intent(in   ) :: tl                  (ncrms,plev)  ! Global grid temperature (K)
  real(r8), intent(in   ) :: ql                  (ncrms,plev)  ! Global grid water vapor (g/g)
  real(r8), intent(in   ) :: qccl                (ncrms,plev)  ! Global grid cloud liquid water (g/g)
  real(r8), intent(in   ) :: qiil                (ncrms,plev)  ! Global grid cloud ice (g/g)
  real(r8), intent(in   ) :: ul                  (ncrms,plev)  ! Global grid u (icrm,m/s)
  real(r8), intent(in   ) :: vl                  (ncrms,plev)  ! Global grid v (icrm,m/s)
#ifdef CLUBB_CRM
  real(r8), intent(inout), target :: clubb_buffer(ncrms,crm_nx, crm_ny, crm_nz+1,1:nclubbvars)
  real(r8), intent(  out) :: crm_cld             (ncrms,crm_nx, crm_ny, crm_nz+1)
  real(r8), intent(  out) :: clubb_tk            (ncrms,crm_nx, crm_ny, crm_nz)
  real(r8), intent(  out) :: clubb_tkh           (ncrms,crm_nx, crm_ny, crm_nz)
  real(r8), intent(  out) :: relvar              (ncrms,crm_nx, crm_ny, crm_nz)
  real(r8), intent(  out) :: accre_enhan         (ncrms,crm_nx, crm_ny, crm_nz)
  real(r8), intent(  out) :: qclvar              (ncrms,crm_nx, crm_ny, crm_nz)
#endif
  real(r8), intent(  out) :: crm_tk              (ncrms,crm_nx, crm_ny, crm_nz)
  real(r8), intent(  out) :: crm_tkh             (ncrms,crm_nx, crm_ny, crm_nz)
  real(r8), intent(inout) :: cltot               (ncrms)                        ! shaded cloud fraction
  real(r8), intent(inout) :: clhgh               (ncrms)                        ! shaded cloud fraction
  real(r8), intent(inout) :: clmed               (ncrms)                        ! shaded cloud fraction
  real(r8), intent(inout) :: cllow               (ncrms)                        ! shaded cloud fraction
#ifdef CRM3D
  real(r8), intent(  out) :: ultend              (ncrms,plev)                   ! tendency of ul
  real(r8), intent(  out) :: vltend              (ncrms,plev)                   ! tendency of vl
#endif
  real(r8), intent(  out) :: sltend              (ncrms,plev)                   ! tendency of static energy
  real(r8), intent(  out) :: qltend              (ncrms,plev)                   ! tendency of water vapor
  real(r8), intent(  out) :: qcltend             (ncrms,plev)                   ! tendency of cloud liquid water
  real(r8), intent(  out) :: qiltend             (ncrms,plev)                   ! tendency of cloud ice
  real(r8), intent(inout) :: u_crm               (ncrms,crm_nx,crm_ny,crm_nz)   ! CRM u-wind component
  real(r8), intent(inout) :: v_crm               (ncrms,crm_nx,crm_ny,crm_nz)   ! CRM v-wind component
  real(r8), intent(inout) :: w_crm               (ncrms,crm_nx,crm_ny,crm_nz)   ! CRM w-wind component
  real(r8), intent(inout) :: t_crm               (ncrms,crm_nx,crm_ny,crm_nz)   ! CRM temperuture
  real(r8), intent(inout) :: micro_fields_crm    (ncrms,crm_nx,crm_ny,crm_nz,nmicro_fields+1) ! CRM total water
  real(r8), intent(  out) :: cld3d_crm           (ncrms,crm_nx, crm_ny, crm_nz) ! instant 3D cloud fraction
  ! real(r8), intent(  out) :: t_rad               (ncrms,crm_nx, crm_ny, crm_nz) ! rad temperuture
  ! real(r8), intent(  out) :: qv_rad              (ncrms,crm_nx, crm_ny, crm_nz) ! rad vapor
  ! real(r8), intent(  out) :: qc_rad              (ncrms,crm_nx, crm_ny, crm_nz) ! rad cloud water
  ! real(r8), intent(  out) :: qi_rad              (ncrms,crm_nx, crm_ny, crm_nz) ! rad cloud ice
  ! real(r8), intent(  out) :: cld_rad             (ncrms,crm_nx, crm_ny, crm_nz) ! rad cloud fraction
  real(r8), intent(  out) :: t_rad               (ncrms,crm_nx_rad, crm_ny_rad, crm_nz) ! rad temperuture
  real(r8), intent(  out) :: qv_rad              (ncrms,crm_nx_rad, crm_ny_rad, crm_nz) ! rad vapor
  real(r8), intent(  out) :: qc_rad              (ncrms,crm_nx_rad, crm_ny_rad, crm_nz) ! rad cloud water
  real(r8), intent(  out) :: qi_rad              (ncrms,crm_nx_rad, crm_ny_rad, crm_nz) ! rad cloud ice
  real(r8), intent(  out) :: cld_rad             (ncrms,crm_nx_rad, crm_ny_rad, crm_nz) ! rad cloud fraction
#ifdef m2005
  ! real(r8), intent(  out) :: nc_rad              (ncrms,crm_nx, crm_ny, crm_nz) ! rad cloud droplet number (#/kg)
  ! real(r8), intent(  out) :: ni_rad              (ncrms,crm_nx, crm_ny, crm_nz) ! rad cloud ice crystal number (#/kg)
  ! real(r8), intent(  out) :: qs_rad              (ncrms,crm_nx, crm_ny, crm_nz) ! rad cloud snow (kg/kg)
  ! real(r8), intent(  out) :: ns_rad              (ncrms,crm_nx, crm_ny, crm_nz) ! rad cloud snow crystal number (#/kg)
  real(r8), intent(  out) :: nc_rad              (ncrms,crm_nx_rad, crm_ny_rad, crm_nz) ! rad cloud droplet number (#/kg)
  real(r8), intent(  out) :: ni_rad              (ncrms,crm_nx_rad, crm_ny_rad, crm_nz) ! rad cloud ice crystal number (#/kg)
  real(r8), intent(  out) :: qs_rad              (ncrms,crm_nx_rad, crm_ny_rad, crm_nz) ! rad cloud snow (kg/kg)
  real(r8), intent(  out) :: ns_rad              (ncrms,crm_nx_rad, crm_ny_rad, crm_nz) ! rad cloud snow crystal number (#/kg)
  real(r8), intent(  out) :: wvar_crm            (ncrms,crm_nx, crm_ny, crm_nz) ! vertical velocity variance (m/s)
  real(r8), intent(  out) :: aut_crm             (ncrms,crm_nx, crm_ny, crm_nz) ! cloud water autoconversion (1/s)
  real(r8), intent(  out) :: acc_crm             (ncrms,crm_nx, crm_ny, crm_nz) ! cloud water accretion (1/s)
  real(r8), intent(  out) :: evpc_crm            (ncrms,crm_nx, crm_ny, crm_nz) ! cloud water evaporation (1/s)
  real(r8), intent(  out) :: evpr_crm            (ncrms,crm_nx, crm_ny, crm_nz) ! rain evaporation (1/s)
  real(r8), intent(  out) :: mlt_crm             (ncrms,crm_nx, crm_ny, crm_nz) ! ice, snow, graupel melting (1/s)
  real(r8), intent(  out) :: sub_crm             (ncrms,crm_nx, crm_ny, crm_nz) ! ice, snow, graupel sublimation (1/s)
  real(r8), intent(  out) :: dep_crm             (ncrms,crm_nx, crm_ny, crm_nz) ! ice, snow, graupel deposition (1/s)
  real(r8), intent(  out) :: con_crm             (ncrms,crm_nx, crm_ny, crm_nz) ! cloud water condensation(1/s)
  real(r8), intent(  out) :: aut_crm_a           (ncrms,plev)  ! cloud water autoconversion (1/s)
  real(r8), intent(  out) :: acc_crm_a           (ncrms,plev)  ! cloud water accretion (1/s)
  real(r8), intent(  out) :: evpc_crm_a          (ncrms,plev)  ! cloud water evaporation (1/s)
  real(r8), intent(  out) :: evpr_crm_a          (ncrms,plev)  ! rain evaporation (1/s)
  real(r8), intent(  out) :: mlt_crm_a           (ncrms,plev)  ! ice, snow, graupel melting (1/s)
  real(r8), intent(  out) :: sub_crm_a           (ncrms,plev)  ! ice, snow, graupel sublimation (1/s)
  real(r8), intent(  out) :: dep_crm_a           (ncrms,plev)  ! ice, snow, graupel deposition (1/s)
  real(r8), intent(  out) :: con_crm_a           (ncrms,plev)  ! cloud water condensation(1/s)
#endif
  real(r8), intent(  out) :: precc               (ncrms)       ! convective precip rate (m/s)
  real(r8), intent(  out) :: precl               (ncrms)       ! stratiform precip rate (m/s)
  real(r8), intent(  out) :: cld                 (ncrms,plev)  ! cloud fraction
  real(r8), intent(  out) :: cldtop              (ncrms,plev)  ! cloud top pdf
  real(r8), intent(  out) :: gicewp              (ncrms,plev)  ! ice water path
  real(r8), intent(  out) :: gliqwp              (ncrms,plev)  ! ice water path
  real(r8), intent(  out) :: mc                  (ncrms,plev)  ! cloud mass flux
  real(r8), intent(  out) :: mcup                (ncrms,plev)  ! updraft cloud mass flux
  real(r8), intent(  out) :: mcdn                (ncrms,plev)  ! downdraft cloud mass flux
  real(r8), intent(  out) :: mcuup               (ncrms,plev)  ! unsat updraft cloud mass flux
  real(r8), intent(  out) :: mcudn               (ncrms,plev)  ! unsat downdraft cloud mass flux
  real(r8), intent(  out) :: crm_qc              (ncrms,plev)  ! mean cloud water
  real(r8), intent(  out) :: crm_qi              (ncrms,plev)  ! mean cloud ice
  real(r8), intent(  out) :: crm_qs              (ncrms,plev)  ! mean snow
  real(r8), intent(  out) :: crm_qg              (ncrms,plev)  ! mean graupel
    real(r8), intent(  out) :: crm_qr              (ncrms,plev)  ! mean rain
#ifdef m2005
  real(r8), intent(  out) :: crm_nc              (ncrms,plev)  ! mean cloud water  (#/kg)
  real(r8), intent(  out) :: crm_ni              (ncrms,plev)  ! mean cloud ice    (#/kg)
  real(r8), intent(  out) :: crm_ns              (ncrms,plev)  ! mean snow         (#/kg)
  real(r8), intent(  out) :: crm_ng              (ncrms,plev)  ! mean graupel      (#/kg)
  real(r8), intent(  out) :: crm_nr              (ncrms,plev)  ! mean rain         (#/kg)
#ifdef MODAL_AERO
  real(r8), intent(in   )  :: naermod            (ncrms,plev, ntot_amode)    ! Aerosol number concentration [/m3]
  real(r8), intent(in   )  :: vaerosol           (ncrms,plev, ntot_amode)    ! aerosol volume concentration [m3/m3]
  real(r8), intent(in   )  :: hygro              (ncrms,plev, ntot_amode)    ! hygroscopicity of aerosol mode
#endif
#endif
  real(r8), intent(  out) :: mu_crm              (ncrms,plev)       ! mass flux up
  real(r8), intent(  out) :: md_crm              (ncrms,plev)       ! mass flux down
  real(r8), intent(  out) :: du_crm              (ncrms,plev)       ! mass detrainment from updraft
  real(r8), intent(  out) :: eu_crm              (ncrms,plev)       ! mass entrainment from updraft
  real(r8), intent(  out) :: ed_crm              (ncrms,plev)       ! mass detrainment from downdraft
  real(r8)                :: dd_crm              (ncrms,plev)       ! mass entraiment from downdraft
  real(r8), intent(  out) :: jt_crm              (ncrms)            ! index of cloud (convection) top
  real(r8), intent(  out) :: mx_crm              (ncrms)            ! index of cloud (convection) bottom
  real(r8)                :: mui_crm             (ncrms,plev+1)     ! mass flux up at the interface
  real(r8)                :: mdi_crm             (ncrms,plev+1)     ! mass flux down at the interface
  real(r8), intent(  out) :: flux_qt             (ncrms,plev)       ! nonprecipitating water flux           [kg/m2/s]
  real(r8), intent(  out) :: fluxsgs_qt          (ncrms,plev)       ! sgs nonprecipitating water flux    [kg/m2/s]
  real(r8), intent(  out) :: tkez                (ncrms,plev)       ! tke profile               [kg/m/s2]
  real(r8), intent(  out) :: tkesgsz             (ncrms,plev)       ! sgs tke profile        [kg/m/s2]
  real(r8), intent(  out) :: tkz                 (ncrms,plev)       ! tk profile                [m2/s]
  real(r8), intent(  out) :: flux_u              (ncrms,plev)       ! x-momentum flux          [m2/s2]
  real(r8), intent(  out) :: flux_v              (ncrms,plev)       ! y-momentum flux          [m2/s2]
  real(r8), intent(  out) :: flux_qp             (ncrms,plev)       ! precipitating water flux [kg/m2/s or mm/s]
  real(r8), intent(  out) :: pflx                (ncrms,plev)       ! precipitation flux      [m/s]
  real(r8), intent(  out) :: qt_ls               (ncrms,plev)       ! tendency of nonprec water due to large-scale  [kg/kg/s]
  real(r8), intent(  out) :: qt_trans            (ncrms,plev)       ! tendency of nonprec water due to transport  [kg/kg/s]
  real(r8), intent(  out) :: qp_trans            (ncrms,plev)       ! tendency of prec water due to transport [kg/kg/s]
  real(r8), intent(  out) :: qp_fall             (ncrms,plev)       ! tendency of prec water due to fall-out   [kg/kg/s]
  real(r8), intent(  out) :: qp_src              (ncrms,plev)       ! tendency of prec water due to conversion  [kg/kg/s]
  real(r8), intent(  out) :: qp_evp              (ncrms,plev)       ! tendency of prec water due to evp         [kg/kg/s]
  real(r8), intent(  out) :: t_ls                (ncrms,plev)       ! tendency of lwse  due to large-scale        [kg/kg/s] ???
  real(r8), intent(  out) :: prectend            (ncrms)            ! column integrated tendency in precipitating water+ice (kg/m2/s)
  real(r8), intent(  out) :: precstend           (ncrms)            ! column integrated tendency in precipitating ice (kg/m2/s)
  real(r8), intent(  out) :: precsc              (ncrms)            ! convective snow rate (m/s)
  real(r8), intent(  out) :: precsl              (ncrms)            ! stratiform snow rate (m/s)
  real(r8), intent(  out) :: taux_crm            (ncrms)            ! zonal CRM surface stress perturbation (N/m2)
  real(r8), intent(  out) :: tauy_crm            (ncrms)            ! merid CRM surface stress perturbation (N/m2)
  real(r8), intent(  out) :: z0m                 (ncrms)            ! surface stress (N/m2)
  real(r8), intent(  out) :: timing_factor       (ncrms)            ! crm cpu efficiency
  real(r8), intent(  out) :: qc_crm              (ncrms,crm_nx, crm_ny, crm_nz)! CRM cloud water
  real(r8), intent(  out) :: qi_crm              (ncrms,crm_nx, crm_ny, crm_nz)! CRM cloud ice
  real(r8), intent(  out) :: qpc_crm             (ncrms,crm_nx, crm_ny, crm_nz)! CRM precip water
  real(r8), intent(  out) :: qpi_crm             (ncrms,crm_nx, crm_ny, crm_nz)! CRM precip ice
  real(r8), intent(  out) :: prec_crm            (ncrms,crm_nx, crm_ny)        ! CRM precipiation rate at layer center
#ifdef ECPP
  real(r8), intent(  out) :: acen                (ncrms,plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! cloud fraction for each sub-sub class for full time period
  real(r8), intent(  out) :: acen_tf             (ncrms,plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! cloud fraction for end-portion of time period
  real(r8), intent(  out) :: rhcen               (ncrms,plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! relative humidity (0-1)
  real(r8), intent(  out) :: qcloudcen           (ncrms,plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! cloud water (kg/kg)
  real(r8), intent(  out) :: qicecen             (ncrms,plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! cloud ice (kg/kg)
  real(r8), intent(  out) :: qlsinkcen           (ncrms,plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! cloud water loss rate from precipitation (/s??)
  real(r8), intent(  out) :: precrcen            (ncrms,plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! liquid (rain) precipitation rate (kg/m2/s)
  real(r8), intent(  out) :: precsolidcen        (ncrms,plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! solid (rain) precipitation rate (kg/m2/s)
  real(r8), intent(  out) :: qlsink_bfcen        (ncrms,plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! cloud water loss rate from precipitation calculated
                                                                                                 ! cloud water before precipitatinog (/s)
  real(r8), intent(  out) :: qlsink_avgcen       (ncrms,plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! cloud water loss rate from precipitation calculated
                                                                                                 ! from praincen and qlcoudcen averaged over
                                                                                                 ! ntavg1_ss time step (/s??)
  real(r8), intent(  out) :: praincen            (ncrms,plev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)  ! cloud water loss rate from precipitation (kg/kg/s)
  real(r8), intent(  out) :: wwqui_cen           (ncrms,plev)                                   ! vertical velocity variance in quiescent class (m2/s2)
  real(r8), intent(  out) :: wwqui_cloudy_cen    (ncrms,plev)                                   ! vertical velocity variance in quiescent, and cloudy class (m2/s2) at layer boundary
  real(r8), intent(  out) :: abnd                (ncrms,plev+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)! cloud fraction for each sub-sub class for full time period
  real(r8), intent(  out) :: abnd_tf             (ncrms,plev+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)! cloud fraction for end-portion of time period
  real(r8), intent(  out) :: massflxbnd          (ncrms,plev+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR)! sub-class vertical mass flux (kg/m2/s) at layer bottom boundary.
  real(r8), intent(  out) :: wupthresh_bnd       (ncrms,plev+1)                                 ! vertical velocity threshold for updraft (m/s)
  real(r8), intent(  out) :: wdownthresh_bnd     (ncrms,plev+1)                                 ! vertical velocity threshold for downdraft (m/s)
  real(r8), intent(  out) :: wwqui_bnd           (ncrms,plev+1)                                 ! vertical velocity variance in quiescent class (m2/s2)
  real(r8), intent(  out) :: wwqui_cloudy_bnd    (ncrms,plev+1)                                 ! vertical velocity variance in quiescent, and cloudy class (m2/s2)
#endif
  real(r8), intent(  out) :: qtot                (ncrms,20)

  !Parameters
  real(r8),       parameter :: umax = 0.5*crm_dx/crm_dt ! maxumum ampitude of the l.s. wind
  real(r8),       parameter :: wmin = 2.                ! minimum up/downdraft velocity for stat
  real(crm_rknd), parameter :: cwp_threshold = 0.001    ! threshold for cloud condensate for shaded fraction calculation

  !Arrays
  real(crm_rknd), allocatable :: t00      (:)
  real(crm_rknd), allocatable :: tln      (:,:)
  real(crm_rknd), allocatable :: qln      (:,:)
  real(crm_rknd), allocatable :: qccln    (:,:)
  real(crm_rknd), allocatable :: qiiln    (:,:)
  real(crm_rknd), allocatable :: uln      (:,:)
  real(crm_rknd), allocatable :: vln      (:,:)
  real(crm_rknd), allocatable :: cwp      (:,:,:)
  real(crm_rknd), allocatable :: cwph     (:,:,:)
  real(crm_rknd), allocatable :: cwpm     (:,:,:)
  real(crm_rknd), allocatable :: cwpl     (:,:,:)
  logical       , allocatable :: flag_top (:,:,:)
  real(crm_rknd), allocatable :: bflx     (:)
  real(crm_rknd), allocatable :: wnd      (:)
  real(crm_rknd), allocatable :: colprec  (:)
  real(crm_rknd), allocatable :: colprecs (:)
  integer       , allocatable :: gcolindex(:)  ! array of global latitude indices
  real(crm_rknd), allocatable :: cltemp   (:,:,:)
  real(crm_rknd), allocatable :: cmtemp   (:,:,:)
  real(crm_rknd), allocatable :: chtemp   (:,:,:)
  real(crm_rknd), allocatable :: cttemp   (:,:,:)
#ifdef CLUBB_CRM
  !Array indicies for spurious RTM check
  real(kind=core_rknd), allocatable :: rtm_integral_before (:,:)
  real(kind=core_rknd), allocatable :: rtm_integral_after  (:,:)
  real(kind=core_rknd), allocatable :: thlm_integral_before(:,:)
  real(kind=core_rknd), allocatable :: thlm_integral_after (:,:)
  real(kind=core_rknd), allocatable :: thlm_before         (:)
  real(kind=core_rknd), allocatable :: thlm_after          (:)
  real(kind=core_rknd), allocatable :: rtm_column          (:) ! Total water (vapor + liquid)     [kg/kg]
#endif

  !Scalars
#ifdef CLUBB_CRM
  real(kind=core_rknd) :: rtm_flux_top, rtm_flux_sfc
  real(kind=core_rknd) :: thlm_flux_top, thlm_flux_sfc
#endif
  real(r8)        :: factor_xy, idt_gl, tmp
  real(crm_rknd)  :: tmp1, tmp2
  real(crm_rknd)  :: u2z,v2z,w2z
  integer         :: i,j,k,l,ptop,nn,icyc, nstatsteps, icrm
  integer         :: kx
  real(crm_rknd)  :: qsat, omg
  ! real(r8)        :: zs                ! surface elevation
  integer         :: igstep            ! GCM time steps
  integer         :: iseed             ! seed for random perturbation
  ! real(crm_rknd)  :: ntotal_step
  integer         :: myrank, ierr
  real(crm_rknd)  :: fcorz      ! Vertical Coriolis parameter
  real(crm_rknd)  :: fcor     ! Coriolis parameter

  ! whannah - variables for new radiation group method
  real(crm_rknd) :: crm_nx_rad_fac
  real(crm_rknd) :: crm_ny_rad_fac
  integer        :: i_rad
  integer        :: j_rad

  call t_startf('initial allocates and zeros')

  !Allocate local arrays
  allocate( t00      (ncrms)      )
  allocate( tln      (ncrms,plev)       )
  allocate( qln      (ncrms,plev)       )
  allocate( qccln    (ncrms,plev)       )
  allocate( qiiln    (ncrms,plev)       )
  allocate( uln      (ncrms,plev) )
  allocate( vln      (ncrms,plev) )
  allocate( cwp      (ncrms,nx,ny)      )
  allocate( cwph     (ncrms,nx,ny)      )
  allocate( cwpm     (ncrms,nx,ny)      )
  allocate( cwpl     (ncrms,nx,ny)      )
  allocate( flag_top (ncrms,nx,ny)      )
  allocate( bflx     (ncrms)      )
  allocate( wnd      (ncrms)      )
  allocate( colprec  (ncrms)      )
  allocate( colprecs (ncrms)      )
  allocate( gcolindex(pcols)      )
  allocate( cltemp   (ncrms,nx,ny)      )
  allocate( cmtemp   (ncrms,nx,ny)      )
  allocate( chtemp   (ncrms,nx,ny)      )
  allocate( cttemp   (ncrms,nx,ny)      )
#ifdef CLUBB_CRM
  allocate( rtm_integral_before (nx,ny) )
  allocate( rtm_integral_after  (nx,ny) )
  allocate( thlm_integral_before(nx,ny) )
  allocate( thlm_integral_after (nx,ny) )
  allocate( thlm_before         (nzm)   )
  allocate( thlm_after          (nzm)   )
  allocate( rtm_column          (nzm)   )
#endif

  !Initialize local arrays to zero
  t00        = 0
  tln        = 0
  qln        = 0
  qccln      = 0
  qiiln      = 0
  uln        = 0
  vln        = 0
  cwp        = 0
  cwph       = 0
  cwpm       = 0
  cwpl       = 0
  flag_top   = 0
  bflx       = 0
  wnd        = 0
  colprec    = 0
  colprecs   = 0
  gcolindex  = 0
  cltemp     = 0
  cmtemp     = 0
  chtemp     = 0
  cttemp     = 0
#ifdef CLUBB_CRM
  rtm_integral_before  = 0
  rtm_integral_after   = 0
  thlm_integral_before = 0
  thlm_integral_after  = 0
  thlm_before          = 0
  thlm_after           = 0
  rtm_column           = 0
#endif

  call allocate_grid(ncrms)
  call allocate_params(ncrms)
  call allocate_vars(ncrms)
  call allocate_microphysics(ncrms)
  call allocate_tracers(ncrms)
  call allocate_sgs(ncrms)
#ifdef sam1mom
  call allocate_micro_params(ncrms)
#endif
  call t_stopf('initial allocates and zeros')

  call t_startf('before time step loop')
  !MRN: In standalone mode, we need to pass these things in by parameter, not look them up.
#ifdef CRM_STANDALONE
  latitude0 (:) = latitude0_in (:)
  longitude0(:) = longitude0_in(:)
#else
  do icrm = 1 , ncrms
    latitude0 (icrm) = get_rlat_p(lchnk, icol(icrm)) * 57.296_r8
    longitude0(icrm) = get_rlon_p(lchnk, icol(icrm)) * 57.296_r8
  enddo
#endif

  igstep = get_nstep()

#ifdef CRM_DUMP
  do icrm = 1 , ncrms
    call crm_dump_input( igstep,plev,lchnk,icol(icrm),latitude0(icrm),longitude0(icrm),ps(icrm),pmid(icrm,:),pdel(icrm,:),phis(icrm),zmid(icrm,:),zint(icrm,:),qrad_crm(icrm,:,:,:),dt_gl, &
                         ocnfrac(icrm),tau00(icrm),wndls(icrm),bflxls(icrm),fluxu00(icrm),fluxv00(icrm),fluxt00(icrm),fluxq00(icrm),tl(icrm,:),ql(icrm,:),qccl(icrm,:),qiil(icrm,:),   &
                         ul(icrm,:),vl(icrm,:), &
#ifdef CLUBB_CRM
                         clubb_buffer(icrm,:,:,:,:) , &
#endif
                         cltot(icrm),clhgh(icrm),clmed(icrm),cllow(icrm),u_crm(icrm,:,:,:),v_crm(icrm,:,:,:),w_crm(icrm,:,:,:),t_crm(icrm,:,:,:),micro_fields_crm(icrm,:,:,:,:), &
#ifdef m2005
#ifdef MODAL_AERO
                         naermod(icrm,:,:),vaerosol(icrm,:,:),hygro(icrm,:,:) , &
#endif
#endif
                         dd_crm(icrm,:),mui_crm(icrm,:),mdi_crm(icrm,:) )
  enddo
#endif

!-----------------------------------------------

  dostatis  = .false.    ! no statistics are collected.
  idt_gl    = 1._r8/dt_gl
  ptop      = plev-nzm+1
  factor_xy = 1._r8/dble(nx*ny)
  t_rad    = 0
  qv_rad   = 0
  qc_rad   = 0
  qi_rad   = 0
  cld_rad  = 0
#ifdef m2005
  nc_rad  = 0
  ni_rad  = 0
  qs_rad  = 0
  ns_rad  = 0
#endif
  ! zs=phis(icrm)/ggr
  bflx(:) = bflxls(:)
  wnd(:) = wndls(:)

!-----------------------------------------

#ifdef CLUBB_CRM
  if(igstep == 1) then
    lrestart_clubb = .false.
  else
   lrestart_clubb = .true.
  endif
#endif

  call task_init()
  call setparm()

  do icrm = 1 , ncrms
    fcor= 4*pi/86400.*sin(latitude0(icrm)*pi/180.)
    fcorz = sqrt(4.*(2*pi/(3600.*24.))**2-fcor**2)
    fcory(icrm,:) = fcor
    fcorzy(icrm,:) = fcorz
  enddo
  do j=1,ny
    do i=1,nx
      do icrm = 1 , ncrms
        latitude (icrm,i,j) = latitude0 (icrm)
        longitude(icrm,i,j) = longitude0(icrm)
      end do
    end do
  end do

  ! if(ocnfrac(icrm).gt.0.5) then
  !    OCEAN(icrm) = .true.
  ! else
  !    LAND(icrm) = .true.
  ! end if

  ! Create CRM vertical grid and initialize some vertical reference arrays:
  do k = 1, nzm
    do icrm = 1 , ncrms
      z(icrm,k) = zmid(icrm,plev-k+1) - zint(icrm,plev+1)
      zi(icrm,k) = zint(icrm,plev-k+2)- zint(icrm,plev+1)
      pres(icrm,k) = pmid(icrm,plev-k+1)/100.
      prespot(icrm,k)=(1000./pres(icrm,k))**(rgas/cp)
      bet(icrm,k) = ggr/tl(icrm,plev-k+1)
      gamaz(icrm,k)=ggr/cp*z(icrm,k)
    end do ! k
  end do ! k
  do icrm = 1 , ncrms
    zi(icrm,nz) = zint(icrm,plev-nz+2)-zint(icrm,plev+1) !+++mhwang, 2012-02-04
    dz(icrm) = 0.5*(z(icrm,1)+z(icrm,2))
    do k=2,nzm
      adzw(icrm,k) = (z(icrm,k)-z(icrm,k-1))/dz(icrm)
    end do
    adzw(icrm,1)  = 1.
    adzw(icrm,nz) = adzw(icrm,nzm)
  enddo

  !+++mhwang fix the adz bug. (adz needs to be consistent with zi)
  !2012-02-04 Minghuai Wang (minghuai.wang@pnnl.gov)
  do k=1, nzm
    do icrm = 1 , ncrms
      adz(icrm,k)=(zi(icrm,k+1)-zi(icrm,k))/dz(icrm)
      rho(icrm,k) = pdel(icrm,plev-k+1)/ggr/(adz(icrm,k)*dz(icrm))
      if (k >= 2) then
        rhow(icrm,k) = (pmid(icrm,plev-k+2)-pmid(icrm,plev-k+1))/ggr/(adzw(icrm,k)*dz(icrm))
      endif
    end do
  end do

  do icrm = 1 , ncrms
    rhow(icrm,1) = 2.*rhow(icrm,2) - rhow(icrm,3)
#ifdef CLUBB_CRM /* Fix extrapolation for 30 point grid */
    if (  2.*rhow(icrm,nzm) - rhow(icrm,nzm-1) > 0. ) then
       rhow(icrm,nz)= 2.*rhow(icrm,nzm) - rhow(icrm,nzm-1)
    else
       rhow(icrm,nz)= sqrt( rhow(icrm,nzm) )
    endif
#else
    rhow(icrm,nz)= 2.*rhow(icrm,nzm) - rhow(icrm,nzm-1)
#endif /*CLUBB_CRM*/
  enddo

  !  Initialize:
  ! limit the velocity at the very first step:
  do k=1,nzm
    do j=1,ny
      do i=1,nx
        do icrm = 1 , ncrms
          if(u_crm(icrm,1,1,1).eq.u_crm(icrm,2,1,1).and.u_crm(icrm,3,1,2).eq.u_crm(icrm,4,1,2)) then
            u_crm(icrm,i,j,k) = min( umax, max(-umax,u_crm(icrm,i,j,k)) )
            v_crm(icrm,i,j,k) = min( umax, max(-umax,v_crm(icrm,i,j,k)) )*YES3D
          endif
        enddo
      enddo
    enddo
  enddo

  do k = 1 , nzm
    do j = 1 , ny
      do i = 1 , nx
        do icrm = 1 , ncrms
          u          (icrm,i,j,k                ) = u_crm           (icrm,i,j,k                )
          v          (icrm,i,j,k                ) = v_crm           (icrm,i,j,k                )*YES3D
          w          (icrm,i,j,k                ) = w_crm           (icrm,i,j,k                )
          tabs       (icrm,i,j,k                ) = t_crm           (icrm,i,j,k                )
          micro_field(icrm,i,j,k,1:nmicro_fields) = micro_fields_crm(icrm,i,j,k,1:nmicro_fields)
#ifdef sam1mom
          qn         (icrm,i,j,k)                 = micro_fields_crm(icrm,i,j,k,3 )
#endif
#ifdef m2005
          cloudliq   (icrm,i,j,k)                 = micro_fields_crm(icrm,i,j,k,11)
#endif
        enddo
      enddo
    enddo
  enddo


#ifdef m2005
  do j = 1 , ntot_amode
    do k=1, nzm
      do icrm = 1 , ncrms
#ifdef MODAL_AERO
        ! set aerosol data
        l=plev-k+1
        naer (icrm,k,j) = naermod (icrm,l,j)
        vaer (icrm,k,j) = vaerosol(icrm,l,j)
        hgaer(icrm,k,j) = hygro   (icrm,l,j)
#endif
      enddo
    enddo
  enddo
  do k=1, nzm
    do j=1, ny
      do i=1, nx
        do icrm = 1 , ncrms
          if(cloudliq(icrm,i,j,k).gt.0) then
            if(dopredictNc) then
              if( micro_field(icrm,i,j,k,incl).eq.0) micro_field(icrm,i,j,k,incl) = 1.0e6*Nc0/rho(icrm,k)
            endif
          endif
        enddo
      enddo
    enddo
  enddo
#endif

  do k = 1 , nzm
    do j = 1 , ny
      do i = 1 , nx
        do icrm = 1 , ncrms
          CF3D(icrm,i,j,k) = 1.
        enddo
      enddo
    enddo
  enddo

  call micro_init(ncrms)

  ! initialize sgs fields
  call sgs_init(ncrms)

  do icrm = 1 , ncrms
    do k=1,nzm
      u0   (icrm,k) = 0.
      v0   (icrm,k) = 0.
      t0   (icrm,k) = 0.
      t00  (icrm  ) = 0.
      tabs0(icrm,k) = 0.
      q0   (icrm,k) = 0.
      qv0  (icrm,k) = 0.
      qn0  (icrm,k) = 0.
      qp0  (icrm,k) = 0.
      tke0 (icrm,k) = 0.
      do j=1,ny
        do i=1,nx
          t(icrm,i,j,k) = tabs(icrm,i,j,k)+gamaz(icrm,k) &
                    -fac_cond*qcl(icrm,i,j,k)-fac_sub*qci(icrm,i,j,k) &
                    -fac_cond*qpl(icrm,i,j,k)-fac_sub*qpi(icrm,i,j,k)
          colprec (icrm  )=colprec (icrm  )+(qpl(icrm,i,j,k)+qpi(icrm,i,j,k))*pdel(icrm,plev-k+1)
          colprecs(icrm  )=colprecs(icrm  )+qpi(icrm,i,j,k)*pdel(icrm,plev-k+1)
          u0      (icrm,k)=u0      (icrm,k)+u(icrm,i,j,k)
          v0      (icrm,k)=v0      (icrm,k)+v(icrm,i,j,k)
          t0      (icrm,k)=t0      (icrm,k)+t(icrm,i,j,k)
          t00     (icrm  )=t00     (icrm  )+t(icrm,i,j,k)+fac_cond*qpl(icrm,i,j,k)+fac_sub*qpi(icrm,i,j,k)
          tabs0   (icrm,k)=tabs0   (icrm,k)+tabs(icrm,i,j,k)
          q0      (icrm,k)=q0      (icrm,k)+qv(icrm,i,j,k)+qcl(icrm,i,j,k)+qci(icrm,i,j,k)
          qv0     (icrm,k)=qv0     (icrm,k)+qv(icrm,i,j,k)
          qn0     (icrm,k)=qn0     (icrm,k)+qcl(icrm,i,j,k) + qci(icrm,i,j,k)
          qp0     (icrm,k)=qp0     (icrm,k)+qpl(icrm,i,j,k) + qpi(icrm,i,j,k)
          tke0    (icrm,k)=tke0    (icrm,k)+tke(icrm,i,j,k)
        enddo
      enddo
      u0(icrm,k) = u0(icrm,k) * factor_xy
      v0(icrm,k) = v0(icrm,k) * factor_xy
      t0(icrm,k) = t0(icrm,k) * factor_xy
      t00(icrm) = t00(icrm) * factor_xy
      tabs0(icrm,k) = tabs0(icrm,k) * factor_xy
      q0(icrm,k) = q0(icrm,k) * factor_xy
      qv0(icrm,k) = qv0(icrm,k) * factor_xy
      qn0(icrm,k) = qn0(icrm,k) * factor_xy
      qp0(icrm,k) = qp0(icrm,k) * factor_xy
      tke0(icrm,k) = tke0(icrm,k) * factor_xy
#ifdef CLUBB_CRM
      ! Update thetav for CLUBB.  This is needed when we have a higher model top
      ! than is in the sounding, because we subsequently use tv0 to initialize
      ! thv_ds_zt/zm, which appear in CLUBB's anelastic buoyancy terms.
      ! -dschanen UWM 11 Feb 2010
      tv0(icrm,k) = tabs0(icrm,k)*prespot(icrm,k)*(1.+epsv*q0(icrm,k))
#endif
      l = plev-k+1
      uln(icrm,l) = min( umax, max(-umax,ul(icrm,l)) )
      vln(icrm,l) = min( umax, max(-umax,vl(icrm,l)) )*YES3D
      ttend(icrm,k) = (tl(icrm,l)+gamaz(icrm,k)- fac_cond*(qccl(icrm,l)+qiil(icrm,l))-fac_fus*qiil(icrm,l)-t00(icrm))*idt_gl
      qtend(icrm,k) = (ql(icrm,l)+qccl(icrm,l)+qiil(icrm,l)-q0(icrm,k))*idt_gl
      utend(icrm,k) = (uln(icrm,l)-u0(icrm,k))*idt_gl
      vtend(icrm,k) = (vln(icrm,l)-v0(icrm,k))*idt_gl
      ug0(icrm,k) = uln(icrm,l)
      vg0(icrm,k) = vln(icrm,l)
      tg0(icrm,k) = tl(icrm,l)+gamaz(icrm,k)-fac_cond*qccl(icrm,l)-fac_sub*qiil(icrm,l)
      qg0(icrm,k) = ql(icrm,l)+qccl(icrm,l)+qiil(icrm,l)
    end do ! k
    uhl(icrm) = u0(icrm,1)
    vhl(icrm) = v0(icrm,1)
    prectend(icrm)=colprec(icrm)
    precstend(icrm)=colprecs(icrm)
  end do ! k

! estimate roughness length assuming logarithmic profile of velocity near the surface:

  !MRN: This is not on the GPU because it gives the wrong answer for some reason
  do icrm = 1 , ncrms
    z0(icrm) = z0_est(z(icrm,1),bflx(icrm),wnd(icrm),sqrt(tau00(icrm)/rho(icrm,1)))
    z0(icrm) = max(real(0.00001,crm_rknd),min(real(1.,crm_rknd),z0(icrm)))
  enddo

  timing_factor = 0.


#ifdef CLUBB_CRM
  if(doclubb) then
    fluxbu(:,:, :) = fluxu00(:)/rhow(:,1)
    fluxbv(:,:, :) = fluxv00(:)/rhow(:,1)
    fluxbt(:,:, :) = fluxt00(:)/rhow(:,1)
    fluxbq(:,:, :) = fluxq00(:)/rhow(:,1)
  else
  endif
#else
#endif

!--------------------------------------------------
#ifdef sam1mom

  !$acc data copy(rho,tabs0,pres,accrsi,accrsc,coefice,evaps1,evaps2,accrgi,accrgc,evapg1,evapg2,accrrc,evapr1,evapr2)
  if(doprecip) call precip_init(ncrms)
  !$acc wait(1)
  !$acc end data
#endif

  !MRN: Don't want any stochasticity introduced in the standalone.
  !MRN: Need to make sure the first call to crm(...) is not dumped out
  !MRN: Also want to avoid the rabbit hole of dependencies eminating from get_gcol_all_p in phys_grid!
#ifndef CRM_STANDALONE
    if (is_first_step) then
        iseed = get_gcol_p(lchnk,icol(icrm))
        call setperturb(iseed,ncrms)
    end if
#endif

  !--------------------------
  ! whannah - sanity check for new method to calculate radiation
  ! over averaged groups of columns instead of each individually
  if ( mod(nx,crm_nx_rad)==0 .or. mod(nx,crm_nx_rad)==0  ) then
    crm_nx_rad_fac = real(crm_nx_rad,crm_rknd)/real(nx,crm_rknd)
    crm_ny_rad_fac = real(crm_ny_rad,crm_rknd)/real(ny,crm_rknd)
  else
    write(0,*) "crm_nx_rad and crm_ny_rad need to be divisible by nx and ny"
    call endrun('crm main')
  end if

#ifndef CLUBB_CRM
  !--------------------------
  ! do a CLUBB sanity check
  if ( doclubb .or. doclubbnoninter ) then
    write(0,*) "Cannot call CLUBB if -DCLUBB is not in FFLAGS"
    call endrun('crm main')
  endif
#endif
#ifdef CLUBB_CRM
  do icrm = 1 , ncrms
    !------------------------------------------------------------------
    ! Do initialization for UWM CLUBB
    !------------------------------------------------------------------
    up2       (1:nx, 1:ny, 1:nz ) = clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  1)
    vp2       (1:nx, 1:ny, 1:nz ) = clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  2)
    wprtp     (1:nx, 1:ny, 1:nz ) = clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  3)
    wpthlp    (1:nx, 1:ny, 1:nz ) = clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  4)
    wp2       (1:nx, 1:ny, 1:nz ) = clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  5)
    wp3       (1:nx, 1:ny, 1:nz ) = clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  6)
    rtp2      (1:nx, 1:ny, 1:nz ) = clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  7)
    thlp2     (1:nx, 1:ny, 1:nz ) = clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  8)
    rtpthlp   (1:nx, 1:ny, 1:nz ) = clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  9)
    upwp      (1:nx, 1:ny, 1:nz ) = clubb_buffer(icrm,1:nx, 1:ny, 1:nz , 10)
    vpwp      (1:nx, 1:ny, 1:nz ) = clubb_buffer(icrm,1:nx, 1:ny, 1:nz , 11)
    cloud_frac(1:nx, 1:ny, 1:nz ) = clubb_buffer(icrm,1:nx, 1:ny, 1:nz , 12)
    t_tndcy   (1:nx, 1:ny, 1:nzm) = clubb_buffer(icrm,1:nx, 1:ny, 1:nzm, 13)
    qc_tndcy  (1:nx, 1:ny, 1:nzm) = clubb_buffer(icrm,1:nx, 1:ny, 1:nzm, 14)
    qv_tndcy  (1:nx, 1:ny, 1:nzm) = clubb_buffer(icrm,1:nx, 1:ny, 1:nzm, 15)
    u_tndcy   (1:nx, 1:ny, 1:nzm) = clubb_buffer(icrm,1:nx, 1:ny, 1:nzm, 16)
    v_tndcy   (1:nx, 1:ny, 1:nzm) = clubb_buffer(icrm,1:nx, 1:ny, 1:nzm, 17)

    ! since no tracer is carried in the current version of MMF, these
    ! tracer-related restart varialbes are set to zero. +++mhwang, 2011-08
    tracer_tndcy = 0.0
    sclrp2       = 0.0
    sclrprtp     = 0.0
    sclrpthlp    = 0.0
    wpsclrp      = 0.0

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
  enddo
#endif

#ifdef ECPP
  do icrm = 1 , ncrms
    !ntavg1_ss = dt_gl/3   ! one third of GCM time step, 10 minutes
    ntavg1_ss = min(600._r8, dt_gl)       ! 10 minutes  or the GCM timestep, whichever smaller
          ! ntavg1_ss = number of seconds to average between computing categories.
    ntavg2_ss = dt_gl   ! GCM time step
          ! ntavg2_ss = number of seconds to average between outputs.
          !    This must be a multiple of ntavgt1_ss.

    ! ecpp_crm_init has to be called after ntavg1_ss and ntavg2_ss
    ! are set for their values are used in ecpp_crm_init.
    call ecpp_crm_init()

    qlsink    = 0.0
    qlsink_bf = 0.0
    prain     = 0.0
    precr     = 0.0
    precsolid = 0.0
  enddo
#endif

  !+++mhwangtest
  ! test water conservtion problem
  ! ntotal_step = 0.0
  qtot(:,:) = 0.0
  qtotmicro(:,:) = 0.0
  do k=1, nzm
    l=plev-k+1
    do j=1, ny
      do i=1, nx
#ifdef m2005
        qtot(:,1) = qtot(:,1)+((micro_field(:,i,j,k,iqr)+micro_field(:,i,j,k,iqs)+micro_field(:,i,j,k,iqg)) * pdel(:,l)/ggr)/(nx*ny)
#endif
#ifdef sam1mom
        qtot(:,1) = qtot(:,1)+(qpl(:,i,j,k)+qpi(:,i,j,k)) * pdel(:,l)/ggr/(nx*ny)
#endif
      enddo
    enddo
    qtot(:,1) = qtot(:,1) + (ql(:,l)+qccl(:,l)+qiil(:,l)) * pdel(:,l)/ggr
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
  nprint = 1
  ncycle = 0

  nstep  = 0
  !nrad = nstop/nrad0
  day=day0
  call t_stopf('before time step loop')

  !------------------------------------------------------------------
  !   Main time loop
  !------------------------------------------------------------------
  call t_startf('main time loop')
  !$acc data copy(t00, tln, qln, qccln, qiiln, uln, vln, cwp, cwph, cwpm, cwpl, flag_top, bflx, wnd, colprec, colprecs, gcolindex, cltemp, cmtemp, chtemp, cttemp, &
  !$acc&          z, pres, zi, presi, adz, adzw, dz, latitude0, longitude0, z0, uhl, &
  !$acc&          vhl, taux0, tauy0, u, v, w, t, p, tabs, qv, qcl, qpl, qci, qpi, tke2, tk2, dudt, dvdt, dwdt, misc, fluxbu, fluxbv, fluxbt, fluxbq, fluxtu, fluxtv, fluxtt, fluxtq, fzero, &
  !$acc&          precsfc, precssfc, t0, q0, qv0, tabs0, tl0, tv0, u0, v0, tg0, qg0, ug0, vg0, p0, tke0, t01, q01, qp0, qn0, prespot, rho, rhow, bet, gamaz, wsub, qtend, ttend, utend, &
  !$acc&          vtend, sstxy, fcory, fcorzy, latitude, longitude, prec_xy, pw_xy, cw_xy, iw_xy, cld_xy, u200_xy, usfc_xy, v200_xy, vsfc_xy, w500_xy, twle, twsb, precflux, uwle, uwsb, &
  !$acc&          vwle, vwsb, tkelediss, t2leadv, t2legrad, t2lediff, t2lediss, twleadv, twlediff, tadv, tdiff, tlat, tlatqi, qifall, qpfall, w_max, u_max, total_water_before, &
  !$acc&          total_water_after, total_water_evap, total_water_prec, total_water_ls, total_water_clubb, total_energy_before, total_energy_after, total_energy_evap, total_energy_prec, &
  !$acc&          total_energy_ls, total_energy_clubb, total_energy_rad, qtotmicro, CF3D, u850_xy, v850_xy, psfc_xy, swvp_xy, cloudtopheight, echotopheight, cloudtoptemp, &
  !$acc&          micro_field, fluxbmk, fluxtmk, mkwle, mkwsb, mkadv, mklsadv, mkdiff, mstor, qn, qpsrc, qpevp, tracer, fluxbtr, fluxttr, trwle, trwsb, tradv, trdiff, trphys, &
  !$acc&          sgs_field, sgs_field_diag, fluxbsgs, fluxtsgs, sgswle, sgswsb, sgsadv, sgslsadv, sgsdiff, grdf_x, grdf_y, grdf_z, tkesbbuoy, tkesbshear, tkesbdiss, tkesbdiff, accrsc, &
  !$acc&          accrsi, accrrc, coefice, accrgc, accrgi, evaps1, evaps2, evapr1, evapr2, evapg1, evapg2, &
  !$acc&          icol, latitude0_in, longitude0_in, ps, pmid, pdel, phis, zmid, zint, qrad_crm, ocnfrac, tau00, wndls, bflxls, fluxu00, fluxv00, fluxt00, fluxq00, tl, ql, qccl, qiil, ul, vl, &
  !$acc&          crm_tk, crm_tkh, cltot, clhgh, clmed, cllow, sltend, qltend, qcltend, qiltend, &
  !$acc&          u_crm, v_crm, w_crm, t_crm, micro_fields_crm, cld3d_crm, t_rad, qv_rad, qc_rad, qi_rad, cld_rad, t_rad, qv_rad, qc_rad, qi_rad, cld_rad, &
  !$acc&          precc, precl, cld, cldtop, gicewp, gliqwp, mc, mcup, mcdn, mcuup, mcudn, crm_qc, crm_qi, crm_qs, crm_qg, crm_qr, &
  !$acc&          mu_crm, md_crm, du_crm, eu_crm, ed_crm, jt_crm, mx_crm, flux_qt, fluxsgs_qt, tkez, tkesgsz, tkz, flux_u, flux_v, flux_qp, pflx, qt_ls, qt_trans, &
  !$acc&          qp_trans, qp_fall, qp_src, qp_evp, t_ls, prectend, precstend, precsc, precsl, taux_crm, tauy_crm, z0m, timing_factor, qc_crm, qi_crm, qpc_crm, qpi_crm, prec_crm, &
  !$acc&          qtot, dt3, mui_crm, mdi_crm)
  do while (nstep.lt.nstop)
    nstep = nstep + 1
    time = time + dt
    day = day0 + time/86400.

    !$acc parallel loop gang vector default(present) async(1)
    do icrm = 1 , ncrms
      timing_factor(icrm) = timing_factor(icrm)+1
    enddo
    !------------------------------------------------------------------
    !  Check if the dynamical time step should be decreased
    !  to handle the cases when the flow being locally linearly unstable
    !------------------------------------------------------------------
    !TODO: Find a better way to handle kurant
    ncycle = 1 ! call kurant(ncrms)
    do icyc=1,ncycle

      icycle = icyc
      dtn = dt/ncycle
      dt3(na) = dtn
      !$acc update device(dt3) async(1)
      dtfactor = dtn/dt

      !---------------------------------------------
      !  	the Adams-Bashforth scheme in time
      call abcoefs()

      !---------------------------------------------
      !  	initialize stuff:
      call zero(ncrms)

      !-----------------------------------------------------------
      !       Buoyancy term:
      call buoyancy(ncrms)

      !------------------------------------------------------------
      !       Large-scale and surface forcing:
      call forcing(ncrms)

      !$acc parallel loop gang vector collapse(4) default(present) async(1)
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              i_rad = ceiling( real(i,crm_rknd) * crm_nx_rad_fac )
              j_rad = ceiling( real(j,crm_rknd) * crm_ny_rad_fac )
              t(icrm,i,j,k) = t(icrm,i,j,k) + qrad_crm(icrm,i_rad,j_rad,k)*dtn
            enddo
          enddo
        enddo
      enddo

      !----------------------------------------------------------
      !   	suppress turbulence near the upper boundary (spange):
      if (dodamping) call damping(ncrms)

      !---------------------------------------------------------
      !   Ice fall-out
#ifdef CLUBB_CRM
      if ( docloud .or. doclubb ) then
        call ice_fall(ncrms)
      endif
#else
      if(docloud) then
        call ice_fall(ncrms)
      endif
#endif

      !----------------------------------------------------------
      !     Update scalar boundaries after large-scale processes:
      call boundaries(3,ncrms)

      !---------------------------------------------------------
      !     Update boundaries for velocities:
      call boundaries(0,ncrms)

      !-----------------------------------------------
      !     surface fluxes:
      if (dosurface) then
        call crmsurface(bflx,ncrms)
      endif

      !-----------------------------------------------------------
      !  SGS physics:
      if (dosgs) then
        call sgs_proc(ncrms)
      endif

      !----------------------------------------------------------
      !     Fill boundaries for SGS diagnostic fields:

      call boundaries(4,ncrms)

      !-----------------------------------------------
      !       advection of momentum:
      call advect_mom(ncrms)

      !----------------------------------------------------------
      !	SGS effects on momentum:

      if(dosgs) then
        call sgs_mom(ncrms)
      endif

      !-----------------------------------------------------------
      !       Coriolis force:
      if (docoriolis) then
        call coriolis(ncrms)
      endif

      !---------------------------------------------------------
      !       compute rhs of the Poisson equation and solve it for pressure.
      call pressure(ncrms)

      !---------------------------------------------------------
      !       find velocity field at n+1/2 timestep needed for advection of scalars:
      !  Note that at the end of the call, the velocities are in nondimensional form.
      call adams(ncrms)

      !----------------------------------------------------------
      !     Update boundaries for all prognostic scalar fields for advection:
      call boundaries(2,ncrms)

      !---------------------------------------------------------
      !      advection of scalars :
      call advect_all_scalars(ncrms)

      !-----------------------------------------------------------
      !    Convert velocity back from nondimensional form:
      call uvw(ncrms)

      !----------------------------------------------------------
      !     Update boundaries for scalars to prepare for SGS effects:
      call boundaries(3,ncrms)

      !---------------------------------------------------------
      !      SGS effects on scalars :
      if (dosgs) then
        call sgs_scalars(ncrms)
      endif

      !-----------------------------------------------------------
      !       Cloud condensation/evaporation and precipitation processes:
#ifdef CLUBB_CRM
      if(docloud.or.dosmoke.or.doclubb) call micro_proc(ncrms)
#else
      if(docloud.or.dosmoke) call micro_proc(ncrms)
#endif /*CLUBB_CRM*/

      !-----------------------------------------------------------
      !    Compute diagnostics fields:
      call diagnose(ncrms)

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
    do icrm = 1 , ncrms
      call ecpp_crm_stat(ncrms,icrm)
    enddo
#endif /*ECPP*/

    !$acc parallel loop gang vector collapse(3) default(present) async(1)
    do j = 1 , ny
      do i = 1 , nx
        do icrm = 1 , ncrms
          cwp     (icrm,i,j) = 0.
          cwph    (icrm,i,j) = 0.
          cwpm    (icrm,i,j) = 0.
          cwpl    (icrm,i,j) = 0.
          flag_top(icrm,i,j) = .true.
          cltemp  (icrm,i,j) = 0.
          cmtemp  (icrm,i,j) = 0.
          chtemp  (icrm,i,j) = 0.
          cttemp  (icrm,i,j) = 0.
        enddo
      enddo
    enddo

    !$acc parallel loop gang vector collapse(3) default(present) async(1)
    do j=1,ny
      do i=1,nx
        do icrm = 1 , ncrms
          do k=1,nzm
            l = plev-k+1
            tmp1 = rho(icrm,nz-k)*adz(icrm,nz-k)*dz(icrm)*(qcl(icrm,i,j,nz-k)+qci(icrm,i,j,nz-k))
            cwp(icrm,i,j) = cwp(icrm,i,j)+tmp1
            cttemp(icrm,i,j) = max(CF3D(icrm,i,j,nz-k), cttemp(icrm,i,j))
            if(cwp(icrm,i,j).gt.cwp_threshold.and.flag_top(icrm,i,j)) then
              !$acc atomic update
              cldtop(icrm,k) = cldtop(icrm,k) + 1
              flag_top(icrm,i,j) = .false.
            endif
            if(pres(icrm,nz-k).ge.700.) then
              cwpl(icrm,i,j) = cwpl(icrm,i,j)+tmp1
              cltemp(icrm,i,j) = max(CF3D(icrm,i,j,nz-k), cltemp(icrm,i,j))
            else if(pres(icrm,nz-k).lt.400.) then
              cwph(icrm,i,j) = cwph(icrm,i,j)+tmp1
              chtemp(icrm,i,j) = max(CF3D(icrm,i,j,nz-k), chtemp(icrm,i,j))
            else
              cwpm(icrm,i,j) = cwpm(icrm,i,j)+tmp1
              cmtemp(icrm,i,j) = max(CF3D(icrm,i,j,nz-k), cmtemp(icrm,i,j))
            endif

            tmp1 = rho(icrm,k)*adz(icrm,k)*dz(icrm)
            if(tmp1*(qcl(icrm,i,j,k)+qci(icrm,i,j,k)).gt.cwp_threshold) then
              !$acc atomic update
              cld(icrm,l) = cld(icrm,l) + CF3D(icrm,i,j,k)
              if(w(icrm,i,j,k+1)+w(icrm,i,j,k).gt.2*wmin) then
                tmp = rho(icrm,k)*0.5*(w(icrm,i,j,k+1)+w(icrm,i,j,k)) * CF3D(icrm,i,j,k)
                !$acc atomic update
                mcup (icrm,l) = mcup (icrm,l) + tmp
                tmp = rho(icrm,k)*0.5*(w(icrm,i,j,k+1)+w(icrm,i,j,k)) * (1.0 - CF3D(icrm,i,j,k))
                !$acc atomic update
                mcuup(icrm,l) = mcuup(icrm,l) + tmp
              endif
              if(w(icrm,i,j,k+1)+w(icrm,i,j,k).lt.-2*wmin) then
                tmp = rho(icrm,k)*0.5*(w(icrm,i,j,k+1)+w(icrm,i,j,k)) * CF3D(icrm,i,j,k)
                !$acc atomic update
                mcdn (icrm,l) = mcdn (icrm,l) + tmp
                tmp = rho(icrm,k)*0.5*(w(icrm,i,j,k+1)+w(icrm,i,j,k)) * (1. - CF3D(icrm,i,j,k))
                !$acc atomic update
                mcudn(icrm,l) = mcudn(icrm,l) + tmp
              endif
            else
              if(w(icrm,i,j,k+1)+w(icrm,i,j,k).gt.2*wmin) then
                tmp = rho(icrm,k)*0.5*(w(icrm,i,j,k+1)+w(icrm,i,j,k))
                !$acc atomic update
                mcuup(icrm,l) = mcuup(icrm,l) + tmp
              endif
              if(w(icrm,i,j,k+1)+w(icrm,i,j,k).lt.-2*wmin) then
                tmp = rho(icrm,k)*0.5*(w(icrm,i,j,k+1)+w(icrm,i,j,k))
                !$acc atomic update
                mcudn(icrm,l) = mcudn(icrm,l) + tmp
              endif
            endif
            !$acc atomic update
            gliqwp(icrm,l) = gliqwp(icrm,l) + qcl(icrm,i,j,k)
            !$acc atomic update
            gicewp(icrm,l) = gicewp(icrm,l) + qci(icrm,i,j,k)
          enddo
        enddo
      enddo
    enddo

    !$acc parallel loop gang vector collapse(4) default(present) async(1)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            !!! whannah - new method allows for fewer radiation calculation over column groups
            i_rad = ceiling( real(i,crm_rknd) * crm_nx_rad_fac )
            j_rad = ceiling( real(j,crm_rknd) * crm_ny_rad_fac )

            !$acc atomic update
            t_rad  (icrm,i_rad,j_rad,k) = t_rad  (icrm,i_rad,j_rad,k) + tabs(icrm,i,j,k)
            tmp = max(real(0.,crm_rknd),qv(icrm,i,j,k))
            !$acc atomic update
            qv_rad (icrm,i_rad,j_rad,k) = qv_rad (icrm,i_rad,j_rad,k) + tmp
            !$acc atomic update
            qc_rad (icrm,i_rad,j_rad,k) = qc_rad (icrm,i_rad,j_rad,k) + qcl(icrm,i,j,k)
            !$acc atomic update
            qi_rad (icrm,i_rad,j_rad,k) = qi_rad (icrm,i_rad,j_rad,k) + qci(icrm,i,j,k)
            !$acc atomic update
            cld_rad(icrm,i_rad,j_rad,k) = cld_rad(icrm,i_rad,j_rad,k) + CF3D(icrm,i,j,k)
#ifdef m2005
            !$acc atomic update
            nc_rad(icrm,i_rad,j_rad,k) = nc_rad(icrm,i_rad,j_rad,k) + micro_field(icrm,i,j,k,incl)
            !$acc atomic update
            ni_rad(icrm,i_rad,j_rad,k) = ni_rad(icrm,i_rad,j_rad,k) + micro_field(icrm,i,j,k,inci)
            !$acc atomic update
            qs_rad(icrm,i_rad,j_rad,k) = qs_rad(icrm,i_rad,j_rad,k) + micro_field(icrm,i,j,k,iqs)
            !$acc atomic update
            ns_rad(icrm,i_rad,j_rad,k) = ns_rad(icrm,i_rad,j_rad,k) + micro_field(icrm,i,j,k,ins)
#endif
          enddo
        enddo
      enddo
    enddo

    ! Diagnose mass fluxes to drive CAM's convective transport of tracers.
    ! definition of mass fluxes is taken from Xu et al., 2002, QJRMS.
    !$acc parallel loop gang vector collapse(4) default(present) async(1)
    do k=1, nzm+1
      do j=1, ny
        do i=1, nx
          do icrm = 1 , ncrms
            l=plev+1-k+1
            if(w(icrm,i,j,k).gt.0.) then
              kx=max(1, k-1)
              qsat = qsatw_crm(tabs(icrm,i,j,kx),pres(icrm,kx))
              if(qcl(icrm,i,j,kx)+qci(icrm,i,j,kx).gt.min(real(1.e-5,crm_rknd),0.01*qsat)) then
                tmp = rhow(icrm,k)*w(icrm,i,j,k)
                !$acc atomic update
                mui_crm(icrm,l) = mui_crm(icrm,l)+tmp
              endif
            else if (w(icrm,i,j,k).lt.0.) then
              kx=min(k+1, nzm)
              qsat = qsatw_crm(tabs(icrm,i,j,kx),pres(icrm,kx))
              if(qcl(icrm,i,j,kx)+qci(icrm,i,j,kx).gt.min(real(1.e-5,crm_rknd),0.01*qsat)) then
                tmp = rhow(icrm,k)*w(icrm,i,j,k)
                !$acc atomic update
                mdi_crm(icrm,l) = mdi_crm(icrm,l)+tmp
              else if(qpl(icrm,i,j,kx)+qpi(icrm,i,j,kx).gt.1.0e-4) then
                tmp = rhow(icrm,k)*w(icrm,i,j,k)
                !$acc atomic update
                mdi_crm(icrm,l) = mdi_crm(icrm,l)+tmp
              endif
            endif
          enddo
        enddo
      enddo
    enddo

    !$acc parallel loop gang vector collapse(3) default(present) async(1)
    do j=1,ny
      do i=1,nx
        do icrm = 1 , ncrms
          if(cwp (icrm,i,j).gt.cwp_threshold) then
            !$acc atomic update
            cltot(icrm) = cltot(icrm) + cttemp(icrm,i,j)
          endif
          if(cwph(icrm,i,j).gt.cwp_threshold) then
            !$acc atomic update
            clhgh(icrm) = clhgh(icrm) + chtemp(icrm,i,j)
          endif
          if(cwpm(icrm,i,j).gt.cwp_threshold) then
            !$acc atomic update
            clmed(icrm) = clmed(icrm) + cmtemp(icrm,i,j)
          endif
          if(cwpl(icrm,i,j).gt.cwp_threshold) then
            !$acc atomic update
            cllow(icrm) = cllow(icrm) + cltemp(icrm,i,j)
          endif
        enddo
      enddo
    enddo
    !----------------------------------------------------------
  enddo ! main loop
  call t_stopf('main time loop')

  !$acc parallel loop gang vector collapse(4) default(present) async(1)
  do k = 1 , crm_nz
    do j = 1 , crm_ny_rad
      do i = 1 , crm_nx_rad
        do icrm = 1 , ncrms
          tmp1 = crm_nx_rad_fac * crm_ny_rad_fac / real(nstop,crm_rknd)

          t_rad  (icrm,i,j,k) = t_rad  (icrm,i,j,k) * tmp1
          qv_rad (icrm,i,j,k) = qv_rad (icrm,i,j,k) * tmp1
          qc_rad (icrm,i,j,k) = qc_rad (icrm,i,j,k) * tmp1
          qi_rad (icrm,i,j,k) = qi_rad (icrm,i,j,k) * tmp1
          cld_rad(icrm,i,j,k) = cld_rad(icrm,i,j,k) * tmp1
#ifdef m2005
          nc_rad (icrm,i,j,k) = nc_rad (icrm,i,j,k) * tmp1
          ni_rad (icrm,i,j,k) = ni_rad (icrm,i,j,k) * tmp1
          qs_rad (icrm,i,j,k) = qs_rad (icrm,i,j,k) * tmp1
          ns_rad (icrm,i,j,k) = ns_rad (icrm,i,j,k) * tmp1
#endif
        enddo
      enddo
    enddo
  enddo

  !$acc parallel loop gang vector collapse(2) default(present) async(1)
  do k = 1 , plev
    do icrm = 1 , ncrms
      if (k <= ptop-1) then
        ! no CRM tendencies above its top
        tln  (icrm,k) =   tl(icrm,k)
        qln  (icrm,k) =   ql(icrm,k)
        qccln(icrm,k) = qccl(icrm,k)
        qiiln(icrm,k) = qiil(icrm,k)
        uln  (icrm,k) =   ul(icrm,k)
        vln  (icrm,k) =   vl(icrm,k)
      else
        !  Compute tendencies due to CRM:
        tln  (icrm,k) = 0.
        qln  (icrm,k) = 0.
        qccln(icrm,k) = 0.
        qiiln(icrm,k) = 0.
        uln  (icrm,k) = 0.
        vln  (icrm,k) = 0.
      endif
      if (k == 1) then
        colprec (icrm)=0
        colprecs(icrm)=0
      endif
    enddo
  enddo

  !$acc parallel loop gang vector collapse(4) default(present) async(1)
  do k = 1,nzm
    do i=1,nx
      do j=1,ny
        do icrm = 1 , ncrms
          l = plev-k+1
          tmp = (qpl(icrm,i,j,k)+qpi(icrm,i,j,k))*pdel(icrm,plev-k+1)
          !$acc atomic update
          colprec (icrm)=colprec (icrm)+tmp
          tmp = qpi(icrm,i,j,k)*pdel(icrm,plev-k+1)
          !$acc atomic update
          colprecs(icrm)=colprecs(icrm)+tmp
          !$acc atomic update
          tln  (icrm,l) = tln  (icrm,l)+tabs(icrm,i,j,k)
          !$acc atomic update
          qln  (icrm,l) = qln  (icrm,l)+qv  (icrm,i,j,k)
          !$acc atomic update
          qccln(icrm,l) = qccln(icrm,l)+qcl (icrm,i,j,k)
          !$acc atomic update
          qiiln(icrm,l) = qiiln(icrm,l)+qci (icrm,i,j,k)
          !$acc atomic update
          uln  (icrm,l) = uln  (icrm,l)+u   (icrm,i,j,k)
          !$acc atomic update
          vln  (icrm,l) = vln  (icrm,l)+v   (icrm,i,j,k)
        enddo ! k
      enddo
    enddo ! i
  enddo ! i

  !$acc parallel loop gang vector collapse(2) default(present) async(1)
  do k = ptop , plev
    do icrm = 1 , ncrms
      tln  (icrm,k) = tln  (icrm,k) * factor_xy
      qln  (icrm,k) = qln  (icrm,k) * factor_xy
      qccln(icrm,k) = qccln(icrm,k) * factor_xy
      qiiln(icrm,k) = qiiln(icrm,k) * factor_xy
      uln  (icrm,k) = uln  (icrm,k) * factor_xy
      vln  (icrm,k) = vln  (icrm,k) * factor_xy
      if (k == ptop) then
        prectend (icrm)=(colprec(icrm) -prectend (icrm))/ggr*factor_xy * idt_gl
        precstend(icrm)=(colprecs(icrm)-precstend(icrm))/ggr*factor_xy * idt_gl
      endif
    enddo
  enddo

  !$acc parallel loop gang vector collapse(2) default(present) async(1)
  do icrm = 1 , ncrms
    do k = 1 , plev
      sltend (icrm,k) = cp * (tln  (icrm,k) - tl  (icrm,k)) * idt_gl
      qltend (icrm,k) =      (qln  (icrm,k) - ql  (icrm,k)) * idt_gl
      qcltend(icrm,k) =      (qccln(icrm,k) - qccl(icrm,k)) * idt_gl
      qiltend(icrm,k) =      (qiiln(icrm,k) - qiil(icrm,k)) * idt_gl
    enddo
  enddo

  !$acc parallel loop gang vector collapse(2) default(present) async(1)
  do icrm = 1 , ncrms
    do k = ptop,ptop+1
      ! don't use CRM tendencies from two crm top levels,
      ! radiation tendencies are added back after the CRM call (see crm_physics_tend)
      sltend (icrm,k) = 0.
      qltend (icrm,k) = 0.
      qcltend(icrm,k) = 0.
      qiltend(icrm,k) = 0.
    enddo
  enddo

  !$acc parallel loop gang vector collapse(4) default(present) present(tkh,tk) async(1)
  do k = 1 , nz
    do j = 1 , ny
      do i = 1 , nx
        do icrm = 1 , ncrms
          if (k <= nzm) then
            !-------------------------------------------------------------
            ! Save the last step to the permanent core:
            u_crm  (icrm,i,j,k) = u   (icrm,i,j,k)
            v_crm  (icrm,i,j,k) = v   (icrm,i,j,k)
            w_crm  (icrm,i,j,k) = w   (icrm,i,j,k)
            t_crm  (icrm,i,j,k) = tabs(icrm,i,j,k)
            micro_fields_crm(icrm,i,j,k,1:nmicro_fields) = micro_field(icrm,i,j,k,1:nmicro_fields)

#ifdef sam1mom
            micro_fields_crm(icrm,i,j,k,3) = qn(icrm,i,j,k)
#endif
#ifdef m2005
            micro_fields_crm(icrm,i,j,k,11) = cloudliq(icrm,i,j,k)
#endif
            crm_tk   (icrm,i,j,k) = tk  (icrm,i,j,k)
            crm_tkh  (icrm,i,j,k) = tkh (icrm,i,j,k)
            cld3d_crm(icrm,i,j,k) = CF3D(icrm,i,j,k)
#ifdef CLUBB_CRM
            clubb_buffer(icrm,i,j,k, 13) = t_tndcy     (i,j,k)
            clubb_buffer(icrm,i,j,k, 14) = qc_tndcy    (i,j,k)
            clubb_buffer(icrm,i,j,k, 15) = qv_tndcy    (i,j,k)
            clubb_buffer(icrm,i,j,k, 16) = u_tndcy     (i,j,k)
            clubb_buffer(icrm,i,j,k, 17) = v_tndcy     (i,j,k)
            clubb_tk    (icrm,i,j,k)     = tk_clubb    (i,j,k)
            clubb_tkh   (icrm,i,j,k)     = tkh_clubb   (i,j,k)
            relvar      (icrm,i,j,k)     = relvarg     (i,j,k)
            accre_enhan (icrm,i,j,k)     = accre_enhang(i,j,k)
            qclvar      (icrm,i,j,k)     = qclvarg     (i,j,k)
#endif
            qc_crm (icrm,i,j,k) = qcl(icrm,i,j,k)
            qi_crm (icrm,i,j,k) = qci(icrm,i,j,k)
            qpc_crm(icrm,i,j,k) = qpl(icrm,i,j,k)
            qpi_crm(icrm,i,j,k) = qpi(icrm,i,j,k)
#ifdef m2005
            wvar_crm(icrm,i,j,k) = wvar (icrm,i,j,k)
            aut_crm (icrm,i,j,k) = aut1 (icrm,i,j,k)
            acc_crm (icrm,i,j,k) = acc1 (icrm,i,j,k)
            evpc_crm(icrm,i,j,k) = evpc1(icrm,i,j,k)
            evpr_crm(icrm,i,j,k) = evpr1(icrm,i,j,k)
            mlt_crm (icrm,i,j,k) = mlt1 (icrm,i,j,k)
            sub_crm (icrm,i,j,k) = sub1 (icrm,i,j,k)
            dep_crm (icrm,i,j,k) = dep1 (icrm,i,j,k)
            con_crm (icrm,i,j,k) = con1 (icrm,i,j,k)
#endif
        endif
#ifdef CLUBB_CRM
          clubb_buffer(icrm,i,j,k,  1) = up2       (i,j,k)
          clubb_buffer(icrm,i,j,k,  2) = vp2       (i,j,k)
          clubb_buffer(icrm,i,j,k,  3) = wprtp     (i,j,k)
          clubb_buffer(icrm,i,j,k,  4) = wpthlp    (i,j,k)
          clubb_buffer(icrm,i,j,k,  5) = wp2       (i,j,k)
          clubb_buffer(icrm,i,j,k,  6) = wp3       (i,j,k)
          clubb_buffer(icrm,i,j,k,  7) = rtp2      (i,j,k)
          clubb_buffer(icrm,i,j,k,  8) = thlp2     (i,j,k)
          clubb_buffer(icrm,i,j,k,  9) = rtpthlp   (i,j,k)
          clubb_buffer(icrm,i,j,k, 10) = upwp      (i,j,k)
          clubb_buffer(icrm,i,j,k, 11) = vpwp      (i,j,k)
          clubb_buffer(icrm,i,j,k, 12) = cloud_frac(i,j,k)
          crm_cld     (icrm,i,j,k)     = cloud_frac(i,j,k)
#endif
        enddo
      enddo
    enddo
  enddo

  !$acc wait(1)
  !$acc end data

  call t_startf('after time step loop')

  do icrm = 1 , ncrms
    z0m     (icrm) = z0(icrm)
    taux_crm(icrm) = taux0(icrm) / dble(nstop)
    tauy_crm(icrm) = tauy0(icrm) / dble(nstop)

    !---------------------------------------------------------------
    !  Diagnostics:

    ! hm add 9/7/11, change from GCM-time step avg to end-of-timestep
    do k=1,nzm
      l = plev-k+1
      do j=1,ny
        do i=1,nx
          crm_qc(icrm,l) = crm_qc(icrm,l) + qcl(icrm,i,j,k)
          crm_qi(icrm,l) = crm_qi(icrm,l) + qci(icrm,i,j,k)
          crm_qr(icrm,l) = crm_qr(icrm,l) + qpl(icrm,i,j,k)
#ifdef sam1mom
          omg = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(icrm,i,j,k)-tgrmin)*a_gr))
          crm_qg(icrm,l) = crm_qg(icrm,l) + qpi(icrm,i,j,k)*omg
          crm_qs(icrm,l) = crm_qs(icrm,l) + qpi(icrm,i,j,k)*(1.-omg)
#else
          !crm_qg(icrm,l) = crm_qg(icrm,l) + qpi(icrm,i,j,k)
          !crm_qs(icrm,l) = crm_qs(icrm,l) + 0.     ! temporerary solution
          crm_qg(icrm,l) = crm_qg(icrm,l) + micro_field(icrm,i,j,k,iqg)
          crm_qs(icrm,l) = crm_qs(icrm,l) + micro_field(icrm,i,j,k,iqs)

          crm_nc(icrm,l) = crm_nc(icrm,l) + micro_field(icrm,i,j,k,incl)
          crm_ni(icrm,l) = crm_ni(icrm,l) + micro_field(icrm,i,j,k,inci)
          crm_nr(icrm,l) = crm_nr(icrm,l) + micro_field(icrm,i,j,k,inr)
          crm_ng(icrm,l) = crm_ng(icrm,l) + micro_field(icrm,i,j,k,ing)
          crm_ns(icrm,l) = crm_ns(icrm,l) + micro_field(icrm,i,j,k,ins)
#endif
        enddo
      enddo
    enddo

    cld   (icrm,:) = min(1._r8,cld   (icrm,:)/real(nstop,crm_rknd)*factor_xy)
    cldtop(icrm,:) = min(1._r8,cldtop(icrm,:)/real(nstop,crm_rknd)*factor_xy)
    gicewp(icrm,:) = gicewp(icrm,:)*pdel(icrm,:)*1000./ggr/real(nstop,crm_rknd)*factor_xy
    gliqwp(icrm,:) = gliqwp(icrm,:)*pdel(icrm,:)*1000./ggr/real(nstop,crm_rknd)*factor_xy
    mcup  (icrm,:) = mcup (icrm,:) / real(nstop,crm_rknd) * factor_xy
    mcdn  (icrm,:) = mcdn (icrm,:) / real(nstop,crm_rknd) * factor_xy
    mcuup (icrm,:) = mcuup(icrm,:) / real(nstop,crm_rknd) * factor_xy
    mcudn (icrm,:) = mcudn(icrm,:) / real(nstop,crm_rknd) * factor_xy
    mc    (icrm,:) = mcup(icrm,:) + mcdn(icrm,:) + mcuup(icrm,:) + mcudn(icrm,:)

    crm_qc(icrm,:) = crm_qc(icrm,:) * factor_xy
    crm_qi(icrm,:) = crm_qi(icrm,:) * factor_xy
    crm_qs(icrm,:) = crm_qs(icrm,:) * factor_xy
    crm_qg(icrm,:) = crm_qg(icrm,:) * factor_xy
    crm_qr(icrm,:) = crm_qr(icrm,:) * factor_xy
#ifdef m2005
    crm_nc(icrm,:) = crm_nc(icrm,:) * factor_xy
    crm_ni(icrm,:) = crm_ni(icrm,:) * factor_xy
    crm_ns(icrm,:) = crm_ns(icrm,:) * factor_xy
    crm_ng(icrm,:) = crm_ng(icrm,:) * factor_xy
    crm_nr(icrm,:) = crm_nr(icrm,:) * factor_xy

    ! add loop over i,j do get horizontal avg, and flip vertical array
    do k=1,nzm
      l = plev-k+1
      do j=1,ny
        do i=1,nx
          aut_crm_a (icrm,l) = aut_crm_a (icrm,l) + aut1a (icrm,i,j,k)
          acc_crm_a (icrm,l) = acc_crm_a (icrm,l) + acc1a (icrm,i,j,k)
          evpc_crm_a(icrm,l) = evpc_crm_a(icrm,l) + evpc1a(icrm,i,j,k)
          evpr_crm_a(icrm,l) = evpr_crm_a(icrm,l) + evpr1a(icrm,i,j,k)
          mlt_crm_a (icrm,l) = mlt_crm_a (icrm,l) + mlt1a (icrm,i,j,k)
          sub_crm_a (icrm,l) = sub_crm_a (icrm,l) + sub1a (icrm,i,j,k)
          dep_crm_a (icrm,l) = dep_crm_a (icrm,l) + dep1a (icrm,i,j,k)
          con_crm_a (icrm,l) = con_crm_a (icrm,l) + con1a (icrm,i,j,k)
        enddo
      enddo
    enddo

    ! note, rates are divded by dt to get mean rate over step
    aut_crm_a (icrm,:) = aut_crm_a (icrm,:) / dble(nstop) * factor_xy / dt
    acc_crm_a (icrm,:) = acc_crm_a (icrm,:) / dble(nstop) * factor_xy / dt
    evpc_crm_a(icrm,:) = evpc_crm_a(icrm,:) / dble(nstop) * factor_xy / dt
    evpr_crm_a(icrm,:) = evpr_crm_a(icrm,:) / dble(nstop) * factor_xy / dt
    mlt_crm_a (icrm,:) = mlt_crm_a (icrm,:) / dble(nstop) * factor_xy / dt
    sub_crm_a (icrm,:) = sub_crm_a (icrm,:) / dble(nstop) * factor_xy / dt
    dep_crm_a (icrm,:) = dep_crm_a (icrm,:) / dble(nstop) * factor_xy / dt
    con_crm_a (icrm,:) = con_crm_a (icrm,:) / dble(nstop) * factor_xy / dt

#endif
    precc (icrm) = 0.
    precl (icrm) = 0.
    precsc(icrm) = 0.
    precsl(icrm) = 0.
    do j=1,ny
      do i=1,nx
#ifdef sam1mom
        precsfc(icrm,i,j) = precsfc(icrm,i,j)*dz(icrm)/dt/dble(nstop)
        precssfc(icrm,i,j) = precssfc(icrm,i,j)*dz(icrm)/dt/dble(nstop)
#endif
#ifdef m2005
        ! precsfc and precssfc from the subroutine of micro_proc in M2005 have a unit mm/dz
        precsfc(icrm,i,j) = precsfc(icrm,i,j)*dz(icrm)/dt/dble(nstop)     !mm/s/dz --> mm/s
        precssfc(icrm,i,j) = precssfc(icrm,i,j)*dz(icrm)/dt/dble(nstop)   !mm/s/dz --> mm/s
#endif
        if(precsfc(icrm,i,j).gt.10./86400.) then
           precc (icrm) = precc (icrm) + precsfc(icrm,i,j)
           precsc(icrm) = precsc(icrm) + precssfc(icrm,i,j)
        else
           precl (icrm) = precl (icrm) + precsfc(icrm,i,j)
           precsl(icrm) = precsl(icrm) + precssfc(icrm,i,j)
        endif
      enddo
    enddo
    prec_crm(icrm,:,:) = precsfc(icrm,:,:)/1000.           !mm/s --> m/s
    precc   (icrm)     = precc (icrm)*factor_xy/1000.
    precl   (icrm)     = precl (icrm)*factor_xy/1000.
    precsc  (icrm)     = precsc(icrm)*factor_xy/1000.
    precsl  (icrm)     = precsl(icrm)*factor_xy/1000.

    !+++mhwangtest
    ! test water conservtion problem
    do k=1, nzm
      l=plev-k+1
      do j=1, ny
        do i=1, nx
#ifdef m2005
          qtot(icrm,9) = qtot(icrm,9)+((micro_field(icrm,i,j,k,iqr)+micro_field(icrm,i,j,k,iqs)+micro_field(icrm,i,j,k,iqg)) * pdel(icrm,l)/ggr)/(nx*ny)
          qtot(icrm,9) = qtot(icrm,9)+((micro_field(icrm,i,j,k,iqv)+micro_field(icrm,i,j,k,iqci)) * pdel(icrm,l)/ggr)/(nx*ny)
#endif
#ifdef sam1mom
          qtot(icrm,9) = qtot(icrm,9)+((micro_field(icrm,i,j,k,1)+micro_field(icrm,i,j,k,2)) * pdel(icrm,l)/ggr)/(nx*ny)
#endif
        enddo
      enddo
    enddo
    qtot(icrm,9) = qtot(icrm,9) + (precc(icrm)+precl(icrm))*1000 * dt_gl

    cltot(icrm) = cltot(icrm) *factor_xy/nstop
    clhgh(icrm) = clhgh(icrm) *factor_xy/nstop
    clmed(icrm) = clmed(icrm) *factor_xy/nstop
    cllow(icrm) = cllow(icrm) *factor_xy/nstop

    jt_crm(icrm) = plev * 1.0
    mx_crm(icrm) = 1.0
    do k=1, plev
      mu_crm(icrm,k)=0.5*(mui_crm(icrm,k)+mui_crm(icrm,k+1))
      md_crm(icrm,k)=0.5*(mdi_crm(icrm,k)+mdi_crm(icrm,k+1))
      mu_crm(icrm,k)=mu_crm(icrm,k)*ggr/100.          !kg/m2/s --> mb/s
      md_crm(icrm,k)=md_crm(icrm,k)*ggr/100.          !kg/m2/s --> mb/s
      eu_crm(icrm,k) = 0.
      if(mui_crm(icrm,k)-mui_crm(icrm,k+1).gt.0) then
        eu_crm(icrm,k)=(mui_crm(icrm,k)-mui_crm(icrm,k+1))*ggr/pdel(icrm,k)    !/s
      else
        du_crm(icrm,k)=-1.0*(mui_crm(icrm,k)-mui_crm(icrm,k+1))*ggr/pdel(icrm,k)   !/s
      endif
      if(mdi_crm(icrm,k+1)-mdi_crm(icrm,k).lt.0) then
        ed_crm(icrm,k)=(mdi_crm(icrm,k)-mdi_crm(icrm,k+1))*ggr/pdel(icrm,k) ! /s
      else
        dd_crm(icrm,k)=-1.*(mdi_crm(icrm,k)-mdi_crm(icrm,k+1))*ggr/pdel(icrm,k)   !/s
      endif
      if(abs(mu_crm(icrm,k)).gt.1.0e-15.or.abs(md_crm(icrm,k)).gt.1.0e-15) then
        jt_crm(icrm) = min(k*1.0_r8, jt_crm(icrm))
        mx_crm(icrm) = max(k*1.0_r8, mx_crm(icrm))
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
          u2z = u2z+(u(icrm,i,j,k)-u0(icrm,k))**2
          v2z = v2z+(v(icrm,i,j,k)-v0(icrm,k))**2
          w2z = w2z+0.5*(w(icrm,i,j,k+1)**2+w(icrm,i,j,k)**2)
        enddo
      enddo
      !+++mhwang
      ! mkwsb, mkle, mkadv, mkdiff (icrm,also flux_u, flux_v) seem not calculted correclty in the spcam3.5 codes.
      ! Only values at the last time step are calculated, but is averaged over the entire GCM
      ! time step.
      !---mhwang

      tmp1 = dz(icrm)/rhow(icrm,k)
      tmp2 = tmp1/dtn                        ! dtn is calculated inside of the icyc loop.
                                             ! It seems wrong to use it here ???? +++mhwang
      mkwsb (icrm,k,:) = mkwsb (icrm,k,:) * tmp1*rhow(icrm,k) * factor_xy/nstop     !kg/m3/s --> kg/m2/s
      mkwle (icrm,k,:) = mkwle (icrm,k,:) * tmp2*rhow(icrm,k) * factor_xy/nstop     !kg/m3   --> kg/m2/s
      mkadv (icrm,k,:) = mkadv (icrm,k,:) * factor_xy*idt_gl     ! kg/kg  --> kg/kg/s
      mkdiff(icrm,k,:) = mkdiff(icrm,k,:) * factor_xy*idt_gl   ! kg/kg  --> kg/kg/s

      ! qpsrc, qpevp, qpfall in M2005 are calculated in micro_flux.
      qpsrc   (icrm,k) = qpsrc   (icrm,k) * factor_xy*idt_gl
      qpevp   (icrm,k) = qpevp   (icrm,k) * factor_xy*idt_gl
      qpfall  (icrm,k) = qpfall  (icrm,k) * factor_xy*idt_gl   ! kg/kg in M2005 ---> kg/kg/s
      precflux(icrm,k) = precflux(icrm,k) * factor_xy*dz(icrm)/dt/nstop  !kg/m2/dz in M2005 -->kg/m2/s or mm/s (idt_gl=1/dt/nstop)

      l = plev-k+1
      flux_u    (icrm,l) = (uwle(icrm,k) + uwsb(icrm,k))*tmp1*factor_xy/nstop
      flux_v    (icrm,l) = (vwle(icrm,k) + vwsb(icrm,k))*tmp1*factor_xy/nstop
#ifdef sam1mom
      flux_qt   (icrm,l) = mkwle(icrm,k,1) + mkwsb(icrm,k,1)
      fluxsgs_qt(icrm,l) = mkwsb(icrm,k,1)
      flux_qp   (icrm,l) = mkwle(icrm,k,2) + mkwsb(icrm,k,2)
      qt_trans  (icrm,l) = mkadv(icrm,k,1) + mkdiff(icrm,k,1)
      qp_trans  (icrm,l) = mkadv(icrm,k,2) + mkdiff(icrm,k,2)
#endif
#ifdef m2005
      flux_qt   (icrm,l) = mkwle(icrm,k,1   ) + mkwsb(icrm,k,1   ) +  &
                         mkwle(icrm,k,iqci) + mkwsb(icrm,k,iqci)
      fluxsgs_qt(icrm,l) = mkwsb(icrm,k,1   ) + mkwsb(icrm,k,iqci)
      flux_qp   (icrm,l) = mkwle(icrm,k,iqr) + mkwsb(icrm,k,iqr) +  &
                         mkwle(icrm,k,iqs) + mkwsb(icrm,k,iqs) + mkwle(icrm,k,iqg) + mkwsb(icrm,k,iqg)
      qt_trans  (icrm,l) = mkadv (icrm,k,1) + mkadv (icrm,k,iqci) + &
                         mkdiff(icrm,k,1) + mkdiff(icrm,k,iqci)
      qp_trans  (icrm,l) = mkadv (icrm,k,iqr) + mkadv (icrm,k,iqs) + mkadv (icrm,k,iqg) + &
                         mkdiff(icrm,k,iqr) + mkdiff(icrm,k,iqs) + mkdiff(icrm,k,iqg)
#endif
      tkesgsz   (icrm,l)= rho(icrm,k)*sum(tke(icrm,1:nx,1:ny,k))*factor_xy
      tkez      (icrm,l)= rho(icrm,k)*0.5*(u2z+v2z*YES3D+w2z)*factor_xy + tkesgsz(icrm,l)
      tkz       (icrm,l) = sum(tk(icrm,1:nx, 1:ny, k)) * factor_xy
      pflx      (icrm,l) = precflux(icrm,k)/1000.       !mm/s  -->m/s

      qp_fall   (icrm,l) = qpfall(icrm,k)
      qp_evp    (icrm,l) = qpevp(icrm,k)
      qp_src    (icrm,l) = qpsrc(icrm,k)

      qt_ls     (icrm,l) = qtend(icrm,k)
      t_ls      (icrm,l) = ttend(icrm,k)
    enddo

#ifdef ECPP
    abnd         (icrm,:,:,:,:)=0.0
    abnd_tf      (icrm,:,:,:,:)=0.0
    massflxbnd   (icrm,:,:,:,:)=0.0
    acen         (icrm,:,:,:,:)=0.0
    acen_tf      (icrm,:,:,:,:)=0.0
    rhcen        (icrm,:,:,:,:)=0.0
    qcloudcen    (icrm,:,:,:,:)=0.0
    qicecen      (icrm,:,:,:,:)=0.0
    qlsinkcen    (icrm,:,:,:,:)=0.0
    precrcen     (icrm,:,:,:,:)=0.0
    precsolidcen (icrm,:,:,:,:)=0.0
    qlsink_bfcen (icrm,:,:,:,:)=0.0
    qlsink_avgcen(icrm,:,:,:,:)=0.0
    praincen     (icrm,:,:,:,:)=0.0

    wupthresh_bnd   (icrm,:)=0.0
    wdownthresh_bnd (icrm,:)=0.0
    wwqui_cen       (icrm,:)=0.0
    wwqui_bnd       (icrm,:)=0.0
    wwqui_cloudy_cen(icrm,:)=0.0
    wwqui_cloudy_bnd(icrm,:)=0.0

    ! default is clear, non-precipitating, and quiescent class
    abnd   (icrm,:,1,1,1)=1.0
    abnd_tf(icrm,:,1,1,1)=1.0
    acen   (icrm,:,1,1,1)=1.0
    acen_tf(icrm,:,1,1,1)=1.0
    do k=1, nzm
      l=plev-k+1
      acen            (icrm,l,:,:,:) = area_cen_sum        (k,:,1:ncls_ecpp_in,:)
      acen_tf         (icrm,l,:,:,:) = area_cen_final      (k,:,1:ncls_ecpp_in,:)
      rhcen           (icrm,l,:,:,:) = rh_cen_sum          (k,:,1:ncls_ecpp_in,:)
      qcloudcen       (icrm,l,:,:,:) = qcloud_cen_sum      (k,:,1:ncls_ecpp_in,:)
      qicecen         (icrm,l,:,:,:) = qice_cen_sum        (k,:,1:ncls_ecpp_in,:)
      qlsinkcen       (icrm,l,:,:,:) = qlsink_cen_sum      (k,:,1:ncls_ecpp_in,:)
      precrcen        (icrm,l,:,:,:) = precr_cen_sum       (k,:,1:ncls_ecpp_in,:)
      precsolidcen    (icrm,l,:,:,:) = precsolid_cen_sum   (k,:,1:ncls_ecpp_in,:)
      wwqui_cen       (icrm,l)       = wwqui_cen_sum       (k)
      wwqui_cloudy_cen(icrm,l)       = wwqui_cloudy_cen_sum(k)
      qlsink_bfcen    (icrm,l,:,:,:) = qlsink_bf_cen_sum   (k,:,1:ncls_ecpp_in,:)
      qlsink_avgcen   (icrm,l,:,:,:) = qlsink_avg_cen_sum  (k,:,1:ncls_ecpp_in,:)
      praincen        (icrm,l,:,:,:) = prain_cen_sum       (k,:,1:ncls_ecpp_in,:)
    enddo
    do k=1, nzm+1
      l=plev+1-k+1
      abnd            (icrm,l,:,:,:) = area_bnd_sum        (k,:,1:ncls_ecpp_in,:)
      abnd_tf         (icrm,l,:,:,:) = area_bnd_final      (k,:,1:ncls_ecpp_in,:)
      massflxbnd      (icrm,l,:,:,:) = mass_bnd_sum        (k,:,1:ncls_ecpp_in,:)
      wupthresh_bnd   (icrm,l)       = wup_thresh          (k)
      wdownthresh_bnd (icrm,l)       = wdown_thresh        (k)
      wwqui_bnd       (icrm,l)       = wwqui_bnd_sum       (k)
      wwqui_cloudy_bnd(icrm,l)       = wwqui_cloudy_bnd_sum(k)
    enddo
#endif

    timing_factor(icrm) = timing_factor(icrm) / nstop

#ifdef CLUBB_CRM
    ! Deallocate CLUBB variables, etc.
    ! -UWM
    if ( doclubb .or. doclubbnoninter ) call clubb_sgs_cleanup( )
#endif
#ifdef ECPP
    ! Deallocate ECPP variables
    call ecpp_crm_cleanup ()
#endif

#ifdef CRM_DUMP
    call crm_dump_output( igstep,plev,crm_tk(icrm,:,:,:),crm_tkh(icrm,:,:,:),cltot(icrm),clhgh(icrm),clmed(icrm),cllow(icrm),sltend(icrm,:),u_crm(icrm,:,:,:),v_crm(icrm,:,:,:),&
                          w_crm(icrm,:,:,:),t_crm(icrm,:,:,:),micro_fields_crm(icrm,:,:,:,:),qltend(icrm,:),qcltend(icrm,:),qiltend(icrm,:),t_rad(icrm,:,:,:),qv_rad(icrm,:,:,:),&
                          qc_rad(icrm,:,:,:),qi_rad(icrm,:,:,:),cld_rad(icrm,:,:,:),cld3d_crm(icrm,:,:,:), &
#ifdef CLUBB_CRM
                          clubb_buffer(icrm,:,:,:,:),crm_cld(icrm,:,:,:),clubb_tk(icrm,:,:,:),clubb_tkh(icrm,:,:,:),relvar(icrm,:,:,:),accre_enhan(icrm,:,:,:),qclvar(icrm,:,:,:) , &
#endif
#ifdef CRM3D
                          ultend(icrm,:),vltend(icrm,:) , &
#endif
#ifdef m2005
                          nc_rad(icrm,:,:,:),ni_rad(icrm,:,:,:),qs_rad(icrm,:,:,:),ns_rad(icrm,:,:,:),wvar_crm(icrm,:,:,:),aut_crm(icrm,:,:,:),acc_crm(icrm,:,:,:),evpc_crm(icrm,:,:,:), &
                          evpr_crm(icrm,:,:,:),mlt_crm(icrm,:,:,:),sub_crm(icrm,:,:,:),dep_crm(icrm,:,:,:),con_crm(icrm,:,:,:),aut_crm_a(icrm,:),acc_crm_a(icrm,:),evpc_crm_a(icrm,:), &
                          evpr_crm_a(icrm,:),mlt_crm_a(icrm,:),sub_crm_a(icrm,:),dep_crm_a(icrm,:),con_crm_a(icrm,:),crm_nc(icrm,:),crm_ni(icrm,:),crm_ns(icrm,:),crm_ng(icrm,:),crm_nr(icrm,:), &
#endif
#ifdef ECPP
                          acen(icrm,:,:,:,:),acen_tf(icrm,:,:,:,:),rhcen(icrm,:,:,:,:),qcloudcen(icrm,:,:,:,:),qicecen(icrm,:,:,:,:),qlsinkcen(icrm,:,:,:,:),precrcen(icrm,:,:,:,:),&
                          precsolidcen(icrm,:,:,:,:),qlsink_bfcen(icrm,:,:,:,:),qlsink_avgcen(icrm,:,:,:,:),praincen(icrm,:,:,:,:),wwqui_cen(icrm,:),wwqui_cloudy_cen(icrm,:), &
                          abnd(icrm,:,:,:,:),abnd_tf(icrm,:,:,:,:),massflxbnd(icrm,:,:,:,:),wupthresh_bnd(icrm,:),wdownthresh_bnd(icrm,:),wwqui_bnd(icrm,:),wwqui_cloudy_bnd(icrm,:), &
#endif
                          precc(icrm),precl(icrm),cld(icrm,:),cldtop(icrm,:),gicewp(icrm,:),gliqwp(icrm,:),mc(icrm,:),mcup(icrm,:),mcdn(icrm,:),mcuup(icrm,:),mcudn(icrm,:),crm_qc(icrm,:), &
                          crm_qi(icrm,:),crm_qs(icrm,:),crm_qg(icrm,:),crm_qr(icrm,:),mu_crm(icrm,:),md_crm(icrm,:),du_crm(icrm,:),eu_crm(icrm,:),ed_crm(icrm,:),dd_crm(icrm,:),jt_crm(icrm), &
                          mx_crm(icrm),mui_crm(icrm,:),mdi_crm(icrm,:),flux_qt(icrm,:),fluxsgs_qt(icrm,:),tkez(icrm,:),tkesgsz(icrm,:),tkz(icrm,:),flux_u(icrm,:),flux_v(icrm,:),flux_qp(icrm,:), &
                          pflx(icrm,:),qt_ls(icrm,:),qt_trans(icrm,:),qp_trans(icrm,:),qp_fall(icrm,:),qp_src(icrm,:),qp_evp(icrm,:),t_ls(icrm,:),prectend(icrm),precstend(icrm),precsc(icrm), &
                          precsl(icrm),taux_crm(icrm),tauy_crm(icrm),z0m(icrm),timing_factor(icrm),qc_crm(icrm,:,:,:),qi_crm(icrm,:,:,:),qpc_crm(icrm,:,:,:),qpi_crm(icrm,:,:,:), &
                          prec_crm(icrm,:,:),qtot(icrm,:) )
#endif
  enddo

  call deallocate_grid()
  call deallocate_params()
  call deallocate_vars()
  call deallocate_microphysics()
  call deallocate_tracers()
  call deallocate_sgs()
#ifdef sam1mom
  call deallocate_micro_params()
#endif

  deallocate( t00       )
  deallocate( tln       )
  deallocate( qln       )
  deallocate( qccln     )
  deallocate( qiiln     )
  deallocate( uln       )
  deallocate( vln       )
  deallocate( cwp       )
  deallocate( cwph      )
  deallocate( cwpm      )
  deallocate( cwpl      )
  deallocate( flag_top  )
  deallocate( bflx      )
  deallocate( wnd       )
  deallocate( colprec   )
  deallocate( colprecs  )
  deallocate( gcolindex )
  deallocate( cltemp    )
  deallocate( cmtemp    )
  deallocate( chtemp    )
  deallocate( cttemp    )
#ifdef CLUBB_CRM
  deallocate( rtm_integral_before  )
  deallocate( rtm_integral_after   )
  deallocate( thlm_integral_before )
  deallocate( thlm_integral_after  )
  deallocate( thlm_before          )
  deallocate( thlm_after           )
  deallocate( rtm_column           )
#endif

  call t_stopf('after time step loop')

  end subroutine crm


end module crm_module
