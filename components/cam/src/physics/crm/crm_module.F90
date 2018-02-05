
module crm_module
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
#ifdef sam1mom
  use precip_init_mod
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
subroutine crm(lchnk, icol, ncrms, &
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
    integer , intent(in   ) :: ncrms                            ! Number of CRM instances
    integer , intent(in   ) :: plev                             ! number of levels in parent model
    real(r8), intent(in   ) :: dt_gl                            ! global model's time step
    integer , intent(in   ) :: icol                (ncrms)     ! column identifier (only for lat/lon and random seed)
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

!  Local space:
    real(r8), intent(  out) :: qtot                (ncrms,20)

    !  Local space:
    real(r8),       parameter :: umax = 0.5*crm_dx/crm_dt ! maxumum ampitude of the l.s. wind
    real(r8),       parameter :: wmin = 2.                ! minimum up/downdraft velocity for stat
    real(crm_rknd), parameter :: cwp_threshold = 0.001    ! threshold for cloud condensate for shaded fraction calculation
    real(crm_rknd)  :: t00(ncrms)
    real(crm_rknd)  :: tln(plev), qln(plev), qccln(plev), qiiln(plev), uln(ncrms,plev), vln(ncrms,plev)
    real(crm_rknd)  :: cwp(nx,ny), cwph(nx,ny), cwpm(nx,ny), cwpl(nx,ny)
    real(r8)        :: factor_xy, idt_gl
    real(crm_rknd)  :: tmp1, tmp2
    real(crm_rknd)  :: u2z,v2z,w2z
    integer         :: i,j,k,l,ptop,nn,icyc, nstatsteps, icrm
    integer         :: kx
    logical         :: flag_top(nx,ny)
    real(crm_rknd)  :: bflx(ncrms), wnd(ncrms), qsat, omg
    real(crm_rknd)  :: colprec(ncrms),colprecs(ncrms)
    ! real(r8)        :: zs                ! surface elevation
    integer         :: igstep            ! GCM time steps
    integer         :: iseed             ! seed for random perturbation
    integer         :: gcolindex(pcols)  ! array of global latitude indices
    real(crm_rknd)  :: cltemp(nx,ny), cmtemp(nx,ny), chtemp(nx, ny), cttemp(nx, ny)
    ! real(crm_rknd)  :: ntotal_step
    integer         :: myrank, ierr
    real(crm_rknd)  :: fcorz      ! Vertical Coriolis parameter
    real(crm_rknd)  :: fcor     ! Coriolis parameter
#ifdef CLUBB_CRM
    !Array indicies for spurious RTM check
    real(kind=core_rknd) :: rtm_integral_before (nx,ny), rtm_integral_after (nx,ny), rtm_flux_top, rtm_flux_sfc
    real(kind=core_rknd) :: thlm_integral_before(nx,ny), thlm_integral_after(nx,ny), thlm_before(nzm), thlm_after(nzm), thlm_flux_top, thlm_flux_sfc
    real(kind=core_rknd) :: rtm_column(nzm) ! Total water (vapor + liquid)     [kg/kg]
#endif

  ! whannah - variables for new radiation group method
  real(crm_rknd) :: crm_nx_rad_fac
  real(crm_rknd) :: crm_ny_rad_fac
  integer        :: i_rad
  integer        :: j_rad

  call allocate_grid(ncrms)
  call allocate_params(ncrms)
  call allocate_vars(ncrms)
  call allocate_microphysics(ncrms)
  call allocate_tracers(ncrms)
  call allocate_sgs(ncrms)

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
  t_rad  (:,:,:,:) = 0.
  qv_rad (:,:,:,:) = 0.
  qc_rad (:,:,:,:) = 0.
  qi_rad (:,:,:,:) = 0.
  cld_rad(:,:,:,:) = 0.
#ifdef m2005
  nc_rad(:,:,:,:) = 0.0
  ni_rad(:,:,:,:) = 0.0
  qs_rad(:,:,:,:) = 0.0
  ns_rad(:,:,:,:) = 0.0
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
      latitude (:,i,j) = latitude0 (:)
      longitude(:,i,j) = longitude0(:)
    end do
  end do

  ! if(ocnfrac(icrm).gt.0.5) then
  !    OCEAN(icrm) = .true.
  ! else
  !    LAND(icrm) = .true.
  ! end if

  ! Create CRM vertical grid and initialize some vertical reference arrays:
  do k = 1, nzm
    z(:,k) = zmid(:,plev-k+1) - zint(:,plev+1)
    zi(:,k) = zint(:,plev-k+2)- zint(:,plev+1)
    pres(:,k) = pmid(:,plev-k+1)/100.
    prespot(:,k)=(1000./pres(:,k))**(rgas/cp)
    bet(:,k) = ggr/tl(:,plev-k+1)
    gamaz(:,k)=ggr/cp*z(:,k)
  end do ! k
  ! zi(icrm,nz) =  zint(plev-nz+2)
  zi(:,nz) = zint(:,plev-nz+2)-zint(:,plev+1) !+++mhwang, 2012-02-04

  dz(:) = 0.5*(z(:,1)+z(:,2))
  do k=2,nzm
    adzw(:,k) = (z(:,k)-z(:,k-1))/dz(:)
  end do
  adzw(:,1)  = 1.
  adzw(:,nz) = adzw(:,nzm)
  !+++mhwang fix the adz bug. (adz needs to be consistent with zi)
  !2012-02-04 Minghuai Wang (minghuai.wang@pnnl.gov)
  do k=1, nzm
    adz(:,k)=(zi(:,k+1)-zi(:,k))/dz(:)
  end do

  do k = 1,nzm
    rho(:,k) = pdel(:,plev-k+1)/ggr/(adz(:,k)*dz(:))
  end do
  do k=2,nzm
    ! rhow(icrm,k) = 0.5*(rho(icrm,k)+rho(icrm,k-1))
    !+++mhwang fix the rhow bug (rhow needes to be consistent with pmid)
    !2012-02-04 Minghuai Wang (minghuai.wang@pnnl.gov)
    rhow(:,k) = (pmid(:,plev-k+2)-pmid(:,plev-k+1))/ggr/(adzw(:,k)*dz(:))
  end do
  rhow(:,1) = 2.*rhow(:,2) - rhow(:,3)
  do icrm = 1 , ncrms
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
  colprec=0
  colprecs=0

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

  u          (:,1:nx,1:ny,1:nzm                ) = u_crm           (:,1:nx,1:ny,1:nzm                )
  v          (:,1:nx,1:ny,1:nzm                ) = v_crm           (:,1:nx,1:ny,1:nzm                )*YES3D
  w          (:,1:nx,1:ny,1:nzm                ) = w_crm           (:,1:nx,1:ny,1:nzm                )
  tabs       (:,1:nx,1:ny,1:nzm                ) = t_crm           (:,1:nx,1:ny,1:nzm                )
  micro_field(:,1:nx,1:ny,1:nzm,1:nmicro_fields) = micro_fields_crm(:,1:nx,1:ny,1:nzm,1:nmicro_fields)
#ifdef sam1mom
  qn      (:,1:nx,1:ny,1:nzm) = micro_fields_crm(:,1:nx,1:ny,1:nzm,3 )
#endif

#ifdef m2005
  cloudliq(1:nx,1:ny,1:nzm) = micro_fields_crm(:,1:nx,1:ny,1:nzm,11)
#endif

#ifdef m2005
  do k=1, nzm
#ifdef MODAL_AERO
    ! set aerosol data
    l=plev-k+1
    naer (:,k, 1:ntot_amode) = naermod (:,l, 1:ntot_amode)
    vaer (:,k, 1:ntot_amode) = vaerosol(:,l, 1:ntot_amode)
    hgaer(:,k, 1:ntot_amode) = hygro   (:,l, 1:ntot_amode)
#endif
    do j=1, ny
      do i=1, nx
        if(cloudliq(i,j,k).gt.0) then
          if(dopredictNc) then
            if( micro_field(:,i,j,k,incl).eq.0) micro_field(:,i,j,k,incl) = 1.0e6*Nc0/rho(:,k)
          endif
        endif
      enddo
    enddo
  enddo
#endif

  w(:,:,:,nz)=0.
  wsub (:,:) = 0.      !used in clubb, +++mhwang
  dudt(:,1:nx,1:ny,1:nzm,1:3) = 0.
  dvdt(:,1:nx,1:ny,1:nzm,1:3) = 0.
  dwdt(:,1:nx,1:ny,1:nz,1:3) = 0.
  tke (:,1:nx,1:ny,1:nzm) = 0.
  tk  (:,1:nx,1:ny,1:nzm) = 0.
  tkh (:,1:nx,1:ny,1:nzm) = 0.
  p   (:,1:nx,1:ny,1:nzm) = 0.

  CF3D(:,1:nx,1:ny,1:nzm) = 1.

  call micro_init(ncrms)

  ! initialize sgs fields
  call sgs_init(ncrms)

  do k=1,nzm
    u0(:,k)=0.
    v0(:,k)=0.
    t0(:,k)=0.
    t00(:)=0.
    tabs0(:,k)=0.
    q0(:,k)=0.
    qv0(:,k)=0.
    !+++mhwang these are not initialized ??
    qn0(:,k) = 0.0
    qp0(:,k) = 0.0
    tke0(:,k) = 0.0
    !---mhwang
    do j=1,ny
      do i=1,nx
        t(:,i,j,k) = tabs(:,i,j,k)+gamaz(:,k) &
                  -fac_cond*qcl(:,i,j,k)-fac_sub*qci(:,i,j,k) &
                  -fac_cond*qpl(:,i,j,k)-fac_sub*qpi(:,i,j,k)
        colprec(:)=colprec(:)+(qpl(:,i,j,k)+qpi(:,i,j,k))*pdel(:,plev-k+1)
        colprecs(:)=colprecs(:)+qpi(:,i,j,k)*pdel(:,plev-k+1)
        u0(:,k)=u0(:,k)+u(:,i,j,k)
        v0(:,k)=v0(:,k)+v(:,i,j,k)
        t0(:,k)=t0(:,k)+t(:,i,j,k)
        t00(:)=t00(:)+t(:,i,j,k)+fac_cond*qpl(:,i,j,k)+fac_sub*qpi(:,i,j,k)
        tabs0(:,k)=tabs0(:,k)+tabs(:,i,j,k)
        q0(:,k)=q0(:,k)+qv(:,i,j,k)+qcl(:,i,j,k)+qci(:,i,j,k)
        qv0(:,k) = qv0(:,k) + qv(:,i,j,k)
        qn0(:,k) = qn0(:,k) + qcl(:,i,j,k) + qci(:,i,j,k)
        qp0(:,k) = qp0(:,k) + qpl(:,i,j,k) + qpi(:,i,j,k)
        tke0(:,k)=tke0(:,k)+tke(:,i,j,k)
      enddo
    enddo

    u0(:,k) = u0(:,k) * factor_xy
    v0(:,k) = v0(:,k) * factor_xy
    t0(:,k) = t0(:,k) * factor_xy
    t00(:) = t00(:) * factor_xy
    tabs0(:,k) = tabs0(:,k) * factor_xy
    q0(:,k) = q0(:,k) * factor_xy
    qv0(:,k) = qv0(:,k) * factor_xy
    qn0(:,k) = qn0(:,k) * factor_xy
    qp0(:,k) = qp0(:,k) * factor_xy
    tke0(:,k) = tke0(:,k) * factor_xy
#ifdef CLUBB_CRM
    ! Update thetav for CLUBB.  This is needed when we have a higher model top
    ! than is in the sounding, because we subsequently use tv0 to initialize
    ! thv_ds_zt/zm, which appear in CLUBB's anelastic buoyancy terms.
    ! -dschanen UWM 11 Feb 2010
    tv0(:,k) = tabs0(:,k)*prespot(:,k)*(1.+epsv*q0(:,k))
#endif

    l = plev-k+1
    uln(:,l) = min( umax, max(-umax,ul(:,l)) )
    vln(:,l) = min( umax, max(-umax,vl(:,l)) )*YES3D
    ttend(:,k) = (tl(:,l)+gamaz(:,k)- fac_cond*(qccl(:,l)+qiil(:,l))-fac_fus*qiil(:,l)-t00(:))*idt_gl
    qtend(:,k) = (ql(:,l)+qccl(:,l)+qiil(:,l)-q0(:,k))*idt_gl
    utend(:,k) = (uln(:,l)-u0(:,k))*idt_gl
    vtend(:,k) = (vln(:,l)-v0(:,k))*idt_gl
    ug0(:,k) = uln(:,l)
    vg0(:,k) = vln(:,l)
    tg0(:,k) = tl(:,l)+gamaz(:,k)-fac_cond*qccl(:,l)-fac_sub*qiil(:,l)
    qg0(:,k) = ql(:,l)+qccl(:,l)+qiil(:,l)

  end do ! k

  uhl(:) = u0(:,1)
  vhl(:) = v0(:,1)

! estimate roughness length assuming logarithmic profile of velocity near the surface:

  do icrm = 1 , ncrms
    z0(icrm) = z0_est(z(icrm,1),bflx(icrm),wnd(icrm),sqrt(tau00(icrm)/rho(icrm,1)))
    z0(icrm) = max(real(0.00001,crm_rknd),min(real(1.,crm_rknd),z0(icrm)))
  enddo

  timing_factor = 0.

  prectend(:)=colprec(:)
  precstend(:)=colprecs(:)


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

!---------------------------------------------------
  cld   (:,:) = 0.
  cldtop(:,:) = 0.
  gicewp(:,:) = 0
  gliqwp(:,:) = 0
  mc    (:,:) = 0.
  mcup  (:,:) = 0.
  mcdn  (:,:) = 0.
  mcuup (:,:) = 0.
  mcudn (:,:) = 0.
  crm_qc(:,:) = 0.
  crm_qi(:,:) = 0.
  crm_qs(:,:) = 0.
  crm_qg(:,:) = 0.
  crm_qr(:,:) = 0.
#ifdef m2005
  crm_nc(:,:) = 0.
  crm_ni(:,:) = 0.
  crm_ns(:,:) = 0.
  crm_ng(:,:) = 0.
  crm_nr(:,:) = 0.
  ! hm 8/31/11 add new variables
  aut_crm_a (:,:) = 0.
  acc_crm_a (:,:) = 0.
  evpc_crm_a(:,:) = 0.
  evpr_crm_a(:,:) = 0.
  mlt_crm_a (:,:) = 0.
  sub_crm_a (:,:) = 0.
  dep_crm_a (:,:) = 0.
  con_crm_a (:,:) = 0.

  ! hm 8/31/11 add new output
  ! these are increments added to calculate gcm-grid and time-step avg
  ! note - these values are also averaged over the icycle loop following
  ! the approach for precsfc
  aut1a  = 0.
  acc1a  = 0.
  evpc1a = 0.
  evpr1a = 0.
  mlt1a  = 0.
  sub1a  = 0.
  dep1a  = 0.
  con1a  = 0.
#endif

  mu_crm (:,:) = 0.
  md_crm (:,:) = 0.
  eu_crm (:,:) = 0.
  du_crm (:,:) = 0.
  ed_crm (:,:) = 0.
  dd_crm (:,:) = 0.
  jt_crm (:)   = 0.
  mx_crm (:)   = 0.

  mui_crm(:,:) = 0.
  mdi_crm(:,:) = 0.

  flux_qt   (:,:) = 0.
  flux_u    (:,:) = 0.
  flux_v    (:,:) = 0.
  fluxsgs_qt(:,:) = 0.
  tkez      (:,:) = 0.
  tkesgsz   (:,:) = 0.
  tkz       (:,:) = 0.
  flux_qp   (:,:) = 0.
  pflx      (:,:) = 0.
  qt_trans  (:,:) = 0.
  qp_trans  (:,:) = 0.
  qp_fall   (:,:) = 0.
  qp_evp    (:,:) = 0.
  qp_src    (:,:) = 0.
  qt_ls     (:,:) = 0.
  t_ls      (:,:) = 0.

  uwle(:,:)     = 0.
  uwsb(:,:)     = 0.
  vwle(:,:)     = 0.
  vwsb(:,:)     = 0.
  qpsrc(:,:)    = 0.
  qpevp(:,:)    = 0.
  qpfall  (:,:) = 0.
  precflux(:,:) = 0.

!--------------------------------------------------
#ifdef sam1mom
  do icrm = 1 , ncrms
    if(doprecip) call precip_init(ncrms,icrm)
  enddo
#endif

    !MRN: Don't want any stochasticity introduced in the standalone.
    !MRN: Need to make sure the first call to crm(...) is not dumped out
    !MRN: Also want to avoid the rabbit hole of dependencies eminating from get_gcol_all_p in phys_grid!
#ifndef CRM_STANDALONE
  do icrm = 1 , ncrms
    call get_gcol_all_p(lchnk, pcols, gcolindex)
    iseed = gcolindex(icol(icrm))
    if(u(icrm,1,1,1).eq.u(icrm,2,1,1).and.u(icrm,3,1,2).eq.u(icrm,4,1,2)) call setperturb(iseed,ncrms,icrm)
  enddo
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

  !------------------------------------------------------------------
  !   Main time loop
  !------------------------------------------------------------------
  do while (nstep.lt.nstop)
    nstep = nstep + 1
    time = time + dt
    day = day0 + time/86400.

    timing_factor(:) = timing_factor(:)+1
    !------------------------------------------------------------------
    !  Check if the dynamical time step should be decreased
    !  to handle the cases when the flow being locally linearly unstable
    !------------------------------------------------------------------
    call kurant(ncrms)

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
      call zero(ncrms)

      !-----------------------------------------------------------
      !       Buoyancy term:
      call buoyancy(ncrms)

      !------------------------------------------------------------
      !       Large-scale and surface forcing:
      call forcing(ncrms)

      do k=1,nzm
        do j=1,ny
          do i=1,nx
            i_rad = ceiling( real(i,crm_rknd) * crm_nx_rad_fac )
            j_rad = ceiling( real(j,crm_rknd) * crm_ny_rad_fac )
            t(:,i,j,k) = t(:,i,j,k) + qrad_crm(:,i_rad,j_rad,k)*dtn
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
      do icrm = 1 , ncrms
        call advect_mom(ncrms,icrm)
      enddo

      !----------------------------------------------------------
      !	SGS effects on momentum:

      if(dosgs) then
        do icrm = 1 , ncrms
          call sgs_mom(ncrms,icrm)
        enddo
      endif

      !-----------------------------------------------------------
      !       Coriolis force:
      if (docoriolis) then
        do icrm = 1 , ncrms
          call coriolis(ncrms,icrm)
        enddo
      endif

      !---------------------------------------------------------
      !       compute rhs of the Poisson equation and solve it for pressure.
      do icrm = 1 , ncrms
        call pressure(ncrms,icrm)
      enddo

      !---------------------------------------------------------
      !       find velocity field at n+1/2 timestep needed for advection of scalars:
      !  Note that at the end of the call, the velocities are in nondimensional form.
      do icrm = 1 , ncrms
        call adams(ncrms,icrm)
      enddo

      !----------------------------------------------------------
      !     Update boundaries for all prognostic scalar fields for advection:
      call boundaries(2,ncrms)

      !---------------------------------------------------------
      !      advection of scalars :
      call advect_all_scalars(ncrms)

      !-----------------------------------------------------------
      !    Convert velocity back from nondimensional form:
      do icrm = 1 , ncrms
        call uvw(ncrms,icrm)
      enddo

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
      if(docloud.or.dosmoke.or.doclubb) then
        do icrm = 1 , ncrms
          call micro_proc(ncrms,icrm)
        enddo
      endif
#else
      if(docloud.or.dosmoke) then
        do icrm = 1 , ncrms
          call micro_proc(ncrms,icrm)
        enddo
      endif
#endif /*CLUBB_CRM*/

      !-----------------------------------------------------------
      !    Compute diagnostics fields:
      do icrm = 1 , ncrms
        call diagnose(ncrms,icrm)
      enddo

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

    do icrm = 1 , ncrms
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
            tmp1 = rho(icrm,nz-k)*adz(icrm,nz-k)*dz(icrm)*(qcl(icrm,i,j,nz-k)+qci(icrm,i,j,nz-k))
            cwp(i,j) = cwp(i,j)+tmp1
            cttemp(i,j) = max(CF3D(icrm,i,j,nz-k), cttemp(i,j))
            if(cwp(i,j).gt.cwp_threshold.and.flag_top(i,j)) then
                cldtop(icrm,k) = cldtop(icrm,k) + 1
                flag_top(i,j) = .false.
            endif
            if(pres(icrm,nz-k).ge.700.) then
                cwpl(i,j) = cwpl(i,j)+tmp1
                cltemp(i,j) = max(CF3D(icrm,i,j,nz-k), cltemp(i,j))
            else if(pres(icrm,nz-k).lt.400.) then
                cwph(i,j) = cwph(i,j)+tmp1
                chtemp(i,j) = max(CF3D(icrm,i,j,nz-k), chtemp(i,j))
            else
                cwpm(i,j) = cwpm(i,j)+tmp1
                cmtemp(i,j) = max(CF3D(icrm,i,j,nz-k), cmtemp(i,j))
            endif

            tmp1 = rho(icrm,k)*adz(icrm,k)*dz(icrm)
            if(tmp1*(qcl(icrm,i,j,k)+qci(icrm,i,j,k)).gt.cwp_threshold) then
                 cld(icrm,l) = cld(icrm,l) + CF3D(icrm,i,j,k)
                 if(w(icrm,i,j,k+1)+w(icrm,i,j,k).gt.2*wmin) then
                   mcup (icrm,l) = mcup (icrm,l) + rho(icrm,k)*0.5*(w(icrm,i,j,k+1)+w(icrm,i,j,k)) * CF3D(icrm,i,j,k)
                   mcuup(icrm,l) = mcuup(icrm,l) + rho(icrm,k)*0.5*(w(icrm,i,j,k+1)+w(icrm,i,j,k)) * (1.0 - CF3D(icrm,i,j,k))
                 endif
                 if(w(icrm,i,j,k+1)+w(icrm,i,j,k).lt.-2*wmin) then
                   mcdn (icrm,l) = mcdn (icrm,l) + rho(icrm,k)*0.5*(w(icrm,i,j,k+1)+w(icrm,i,j,k)) * CF3D(icrm,i,j,k)
                   mcudn(icrm,l) = mcudn(icrm,l) + rho(icrm,k)*0.5*(w(icrm,i,j,k+1)+w(icrm,i,j,k)) * (1. - CF3D(icrm,i,j,k))
                 endif
            else
                 if(w(icrm,i,j,k+1)+w(icrm,i,j,k).gt.2*wmin) then
                   mcuup(icrm,l) = mcuup(icrm,l) + rho(icrm,k)*0.5*(w(icrm,i,j,k+1)+w(icrm,i,j,k))
                 endif
                 if(w(icrm,i,j,k+1)+w(icrm,i,j,k).lt.-2*wmin) then
                   mcudn(icrm,l) = mcudn(icrm,l) + rho(icrm,k)*0.5*(w(icrm,i,j,k+1)+w(icrm,i,j,k))
                 endif
            endif

            !!! whannah - new method allows for fewer radiation calculation over column groups
            i_rad = ceiling( real(i,crm_rknd) * crm_nx_rad_fac )
            j_rad = ceiling( real(j,crm_rknd) * crm_ny_rad_fac )

            t_rad  (icrm,i_rad,j_rad,k) = t_rad  (icrm,i_rad,j_rad,k) + tabs(icrm,i,j,k)
            qv_rad (icrm,i_rad,j_rad,k) = qv_rad (icrm,i_rad,j_rad,k) + max(real(0.,crm_rknd),qv(icrm,i,j,k))
            qc_rad (icrm,i_rad,j_rad,k) = qc_rad (icrm,i_rad,j_rad,k) + qcl(icrm,i,j,k)
            qi_rad (icrm,i_rad,j_rad,k) = qi_rad (icrm,i_rad,j_rad,k) + qci(icrm,i,j,k)
            cld_rad(icrm,i_rad,j_rad,k) = cld_rad(icrm,i_rad,j_rad,k) + CF3D(icrm,i,j,k)
#ifdef m2005
            nc_rad(icrm,i_rad,j_rad,k) = nc_rad(icrm,i_rad,j_rad,k) + micro_field(icrm,i,j,k,incl)
            ni_rad(icrm,i_rad,j_rad,k) = ni_rad(icrm,i_rad,j_rad,k) + micro_field(icrm,i,j,k,inci)
            qs_rad(icrm,i_rad,j_rad,k) = qs_rad(icrm,i_rad,j_rad,k) + micro_field(icrm,i,j,k,iqs)
            ns_rad(icrm,i_rad,j_rad,k) = ns_rad(icrm,i_rad,j_rad,k) + micro_field(icrm,i,j,k,ins)
#endif
            gliqwp(icrm,l) = gliqwp(icrm,l) + qcl(icrm,i,j,k)
            gicewp(icrm,l) = gicewp(icrm,l) + qci(icrm,i,j,k)
          enddo
        enddo
      enddo

      ! Diagnose mass fluxes to drive CAM's convective transport of tracers.
      ! definition of mass fluxes is taken from Xu et al., 2002, QJRMS.
      do k=1, nzm+1
        l=plev+1-k+1
        do j=1, ny
          do i=1, nx
            if(w(icrm,i,j,k).gt.0.) then
              kx=max(1, k-1)
              qsat = qsatw_crm(tabs(icrm,i,j,kx),pres(icrm,kx))
              if(qcl(icrm,i,j,kx)+qci(icrm,i,j,kx).gt.min(real(1.e-5,crm_rknd),0.01*qsat)) then
                mui_crm(icrm,l) = mui_crm(icrm,l)+rhow(icrm,k)*w(icrm,i,j,k)
              endif
            else if (w(icrm,i,j,k).lt.0.) then
              kx=min(k+1, nzm)
              qsat = qsatw_crm(tabs(icrm,i,j,kx),pres(icrm,kx))
              if(qcl(icrm,i,j,kx)+qci(icrm,i,j,kx).gt.min(real(1.e-5,crm_rknd),0.01*qsat)) then
                mdi_crm(icrm,l) = mdi_crm(icrm,l)+rhow(icrm,k)*w(icrm,i,j,k)
              else if(qpl(icrm,i,j,kx)+qpi(icrm,i,j,kx).gt.1.0e-4) then
                mdi_crm(icrm,l) = mdi_crm(icrm,l)+rhow(icrm,k)*w(icrm,i,j,k)
              endif
            endif
          enddo
        enddo
      enddo

      do j=1,ny
        do i=1,nx
          if(cwp (i,j).gt.cwp_threshold) cltot(icrm) = cltot(icrm) + cttemp(i,j)
          if(cwph(i,j).gt.cwp_threshold) clhgh(icrm) = clhgh(icrm) + chtemp(i,j)
          if(cwpm(i,j).gt.cwp_threshold) clmed(icrm) = clmed(icrm) + cmtemp(i,j)
          if(cwpl(i,j).gt.cwp_threshold) cllow(icrm) = cllow(icrm) + cltemp(i,j)
        enddo
      enddo

      !        call stepout()
      !----------------------------------------------------------
    enddo
    !----------------------------------------------------------
  enddo ! main loop

  do icrm = 1 , ncrms
    tmp1 = crm_nx_rad_fac * crm_ny_rad_fac / real(nstop,crm_rknd)

    t_rad  (icrm,:,:,:) = t_rad  (icrm,:,:,:) * tmp1
    qv_rad (icrm,:,:,:) = qv_rad (icrm,:,:,:) * tmp1
    qc_rad (icrm,:,:,:) = qc_rad (icrm,:,:,:) * tmp1
    qi_rad (icrm,:,:,:) = qi_rad (icrm,:,:,:) * tmp1
    cld_rad(icrm,:,:,:) = cld_rad(icrm,:,:,:) * tmp1
#ifdef m2005
    nc_rad(icrm,:,:,:) = nc_rad(icrm,:,:,:) * tmp1
    ni_rad(icrm,:,:,:) = ni_rad(icrm,:,:,:) * tmp1
    qs_rad(icrm,:,:,:) = qs_rad(icrm,:,:,:) * tmp1
    ns_rad(icrm,:,:,:) = ns_rad(icrm,:,:,:) * tmp1
#endif

    ! no CRM tendencies above its top
    tln  (1:ptop-1) =   tl(icrm,1:ptop-1)
    qln  (1:ptop-1) =   ql(icrm,1:ptop-1)
    qccln(1:ptop-1) = qccl(icrm,1:ptop-1)
    qiiln(1:ptop-1) = qiil(icrm,1:ptop-1)
    uln  (icrm,1:ptop-1) =   ul(icrm,1:ptop-1)
    vln  (icrm,1:ptop-1) =   vl(icrm,1:ptop-1)

    !  Compute tendencies due to CRM:
    tln (ptop:plev)  = 0.
    qln (ptop:plev)  = 0.
    qccln(ptop:plev) = 0.
    qiiln(ptop:plev) = 0.
    uln (icrm,ptop:plev)  = 0.
    vln (icrm,ptop:plev)  = 0.

    colprec(icrm)=0
    colprecs(icrm)=0
    do k = 1,nzm
      l = plev-k+1
      do i=1,nx
        do j=1,ny
          colprec(icrm)=colprec(icrm)+(qpl(icrm,i,j,k)+qpi(icrm,i,j,k))*pdel(icrm,plev-k+1)
          colprecs(icrm)=colprecs(icrm)+qpi(icrm,i,j,k)*pdel(icrm,plev-k+1)
          tln(l) = tln(l)+tabs(icrm,i,j,k)
          qln(l) = qln(l)+qv(icrm,i,j,k)
          qccln(l)= qccln(l)+qcl(icrm,i,j,k)
          qiiln(l)= qiiln(l)+qci(icrm,i,j,k)
          uln(icrm,l) = uln(icrm,l)+u(icrm,i,j,k)
          vln(icrm,l) = vln(icrm,l)+v(icrm,i,j,k)
        enddo ! k
      enddo
    enddo ! i

    tln(ptop:plev) = tln(ptop:plev) * factor_xy
    qln(ptop:plev) = qln(ptop:plev) * factor_xy
    qccln(ptop:plev) = qccln(ptop:plev) * factor_xy
    qiiln(ptop:plev) = qiiln(ptop:plev) * factor_xy
    uln(icrm,ptop:plev) = uln(icrm,ptop:plev) * factor_xy
    vln(icrm,ptop:plev) = vln(icrm,ptop:plev) * factor_xy

#ifdef SPMOMTRANS
    ! whannah - SP CMT tendencies
    ultend(icrm,:) = (ulnicrm, - ul(icrm,:))*idt_gl
    vltend(icrm,:) = (vlnicrm, - vl(icrm,:))*idt_gl
#endif

    sltend (icrm,:) = cp * (tln   - tl  (icrm,:)) * idt_gl
    qltend (icrm,:) =      (qln   - ql  (icrm,:)) * idt_gl
    qcltend(icrm,:) =      (qccln - qccl(icrm,:)) * idt_gl
    qiltend(icrm,:) =      (qiiln - qiil(icrm,:)) * idt_gl
    prectend (icrm)=(colprec(icrm) -prectend (icrm))/ggr*factor_xy * idt_gl
    precstend(icrm)=(colprecs(icrm)-precstend(icrm))/ggr*factor_xy * idt_gl

    ! don't use CRM tendencies from two crm top levels,
    ! radiation tendencies are added back after the CRM call (see crm_physics_tend)
    sltend (icrm,ptop:ptop+1) = 0.
    qltend (icrm,ptop:ptop+1) = 0.
    qcltend(icrm,ptop:ptop+1) = 0.
    qiltend(icrm,ptop:ptop+1) = 0.
    !-------------------------------------------------------------
    !
    ! Save the last step to the permanent core:
    u_crm  (icrm,1:nx,1:ny,1:nzm) = u   (icrm,1:nx,1:ny,1:nzm)
    v_crm  (icrm,1:nx,1:ny,1:nzm) = v   (icrm,1:nx,1:ny,1:nzm)
    w_crm  (icrm,1:nx,1:ny,1:nzm) = w   (icrm,1:nx,1:ny,1:nzm)
    t_crm  (icrm,1:nx,1:ny,1:nzm) = tabs(icrm,1:nx,1:ny,1:nzm)
    micro_fields_crm(icrm,1:nx,1:ny,1:nzm,1:nmicro_fields) = micro_field(icrm,1:nx,1:ny,1:nzm,1:nmicro_fields)

#ifdef sam1mom
    micro_fields_crm(icrm,1:nx,1:ny,1:nzm,3) = qn(icrm,1:nx,1:ny,1:nzm)
#endif
#ifdef m2005
    micro_fields_crm(icrm,1:nx,1:ny,1:nzm,11) = cloudliq(1:nx,1:ny,1:nzm)
#endif
    crm_tk   (icrm,1:nx,1:ny,1:nzm) = tk  (icrm,1:nx, 1:ny, 1:nzm)
    crm_tkh  (icrm,1:nx,1:ny,1:nzm) = tkh (icrm,1:nx, 1:ny, 1:nzm)
    cld3d_crm(icrm,1:nx,1:ny,1:nzm) = CF3D(icrm,1:nx, 1:ny, 1:nzm)
#ifdef CLUBB_CRM
    clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  1) = up2       (1:nx, 1:ny, 1:nz )
    clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  2) = vp2       (1:nx, 1:ny, 1:nz )
    clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  3) = wprtp     (1:nx, 1:ny, 1:nz )
    clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  4) = wpthlp    (1:nx, 1:ny, 1:nz )
    clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  5) = wp2       (1:nx, 1:ny, 1:nz )
    clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  6) = wp3       (1:nx, 1:ny, 1:nz )
    clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  7) = rtp2      (1:nx, 1:ny, 1:nz )
    clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  8) = thlp2     (1:nx, 1:ny, 1:nz )
    clubb_buffer(icrm,1:nx, 1:ny, 1:nz ,  9) = rtpthlp   (1:nx, 1:ny, 1:nz )
    clubb_buffer(icrm,1:nx, 1:ny, 1:nz , 10) = upwp      (1:nx, 1:ny, 1:nz )
    clubb_buffer(icrm,1:nx, 1:ny, 1:nz , 11) = vpwp      (1:nx, 1:ny, 1:nz )
    clubb_buffer(icrm,1:nx, 1:ny, 1:nz , 12) = cloud_frac(1:nx, 1:ny, 1:nz )
    clubb_buffer(icrm,1:nx, 1:ny, 1:nzm, 13) = t_tndcy   (1:nx, 1:ny, 1:nzm)
    clubb_buffer(icrm,1:nx, 1:ny, 1:nzm, 14) = qc_tndcy  (1:nx, 1:ny, 1:nzm)
    clubb_buffer(icrm,1:nx, 1:ny, 1:nzm, 15) = qv_tndcy  (1:nx, 1:ny, 1:nzm)
    clubb_buffer(icrm,1:nx, 1:ny, 1:nzm, 16) = u_tndcy   (1:nx, 1:ny, 1:nzm)
    clubb_buffer(icrm,1:nx, 1:ny, 1:nzm, 17) = v_tndcy   (1:nx, 1:ny, 1:nzm)

    crm_cld    (icrm,1:nx, 1:ny, 1:nz ) = cloud_frac  (1:nx, 1:ny, 1:nz )
    clubb_tk   (icrm,1:nx, 1:ny, 1:nzm) = tk_clubb    (1:nx, 1:ny, 1:nzm)
    clubb_tkh  (icrm,1:nx, 1:ny, 1:nzm) = tkh_clubb   (1:nx, 1:ny, 1:nzm)
    relvar     (icrm,1:nx, 1:ny, 1:nzm) = relvarg     (1:nx, 1:ny, 1:nzm)
    accre_enhan(icrm,1:nx, 1:ny, 1:nzm) = accre_enhang(1:nx, 1:ny, 1:nzm)
    qclvar     (icrm,1:nx, 1:ny, 1:nzm) = qclvarg     (1:nx, 1:ny, 1:nzm)
#endif

    do k=1,nzm
     do j=1,ny
      do i=1,nx
        qc_crm (icrm,i,j,k) = qcl(icrm,i,j,k)
        qi_crm (icrm,i,j,k) = qci(icrm,i,j,k)
        qpc_crm(icrm,i,j,k) = qpl(icrm,i,j,k)
        qpi_crm(icrm,i,j,k) = qpi(icrm,i,j,k)
#ifdef m2005
        wvar_crm(icrm,i,j,k) = wvar (i,j,k)
        aut_crm (icrm,i,j,k) = aut1 (i,j,k)
        acc_crm (icrm,i,j,k) = acc1 (i,j,k)
        evpc_crm(icrm,i,j,k) = evpc1(i,j,k)
        evpr_crm(icrm,i,j,k) = evpr1(i,j,k)
        mlt_crm (icrm,i,j,k) = mlt1 (i,j,k)
        sub_crm (icrm,i,j,k) = sub1 (i,j,k)
        dep_crm (icrm,i,j,k) = dep1 (i,j,k)
        con_crm (icrm,i,j,k) = con1 (i,j,k)
#endif
        enddo
      enddo
    enddo
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
          aut_crm_a (icrm,l) = aut_crm_a (icrm,l) + aut1a (i,j,k)
          acc_crm_a (icrm,l) = acc_crm_a (icrm,l) + acc1a (i,j,k)
          evpc_crm_a(icrm,l) = evpc_crm_a(icrm,l) + evpc1a(i,j,k)
          evpr_crm_a(icrm,l) = evpr_crm_a(icrm,l) + evpr1a(i,j,k)
          mlt_crm_a (icrm,l) = mlt_crm_a (icrm,l) + mlt1a (i,j,k)
          sub_crm_a (icrm,l) = sub_crm_a (icrm,l) + sub1a (i,j,k)
          dep_crm_a (icrm,l) = dep_crm_a (icrm,l) + dep1a (i,j,k)
          con_crm_a (icrm,l) = con_crm_a (icrm,l) + con1a (i,j,k)
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

  end subroutine crm


end module crm_module
