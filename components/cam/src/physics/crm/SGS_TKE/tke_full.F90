module tke_full_mod

  use shear_prod2D_mod
  use shear_prod3D_mod
  use sat_mod

  implicit none

contains

subroutine tke_full(ncrms,dimx1_d, dimx2_d, dimy1_d, dimy2_d,   &
                    grdf_x, grdf_y, grdf_z, dosmagor,     &
                    tkesbdiss, tkesbshear, tkesbbuoy,     &
                    tke, tk, tkh)
  !-----------------------------------------------------------------------
  ! Purpose: solve the TKE equation
  !-----------------------------------------------------------------------

  use vars
  use params

  implicit none
  integer, intent(in) :: ncrms
  !-----------------------------------------------------------------------
  !!! Interface Arguments
  integer       , intent(in)                 :: dimx1_d     ! scalar dimension parameter
  integer       , intent(in)                 :: dimx2_d     ! scalar dimension parameter
  integer       , intent(in)                 :: dimy1_d     ! scalar dimension parameter
  integer       , intent(in)                 :: dimy2_d     ! scalar dimension parameter
  real(crm_rknd), intent(in), dimension(ncrms,nzm) :: grdf_x      ! grid length in x direction
  real(crm_rknd), intent(in), dimension(ncrms,nzm) :: grdf_y      ! grid length in y direction
  real(crm_rknd), intent(in), dimension(ncrms,nzm) :: grdf_z      ! grid length in z direction
  logical       , intent(in)                 :: dosmagor    ! flag for diagnostic smagorinsky scheme

  real(crm_rknd), intent(out), dimension(nz,ncrms) :: tkesbdiss   ! TKE dissipation
  real(crm_rknd), intent(out), dimension(nz,ncrms) :: tkesbshear  ! TKE production by shear
  real(crm_rknd), intent(out), dimension(nz,ncrms) :: tkesbbuoy   ! TKE production by buoyancy

  real(crm_rknd), intent(  out) :: tke(ncrms, dimx1_s:dimx2_s , dimy1_s:dimy2_s , nzm )   ! SGS TKE
  real(crm_rknd), intent(inout) :: tk(ncrms, dimx1_d:dimx2_d , dimy1_d:dimy2_d , nzm )   ! SGS eddy viscosity
  real(crm_rknd), intent(inout) :: tkh(ncrms, dimx1_d:dimx2_d , dimy1_d:dimy2_d , nzm )   ! SGS eddy conductivity

  !-----------------------------------------------------------------------
  !!! Local Variables
  real(crm_rknd), dimension(nx,ny,nzm,ncrms) :: def2
  real(crm_rknd), dimension(nx,ny,0:nzm,ncrms)     :: buoy_sgs_vert
  real(crm_rknd), dimension(nx,ny,0:nzm,ncrms)     :: a_prod_bu_vert
  real(crm_rknd) :: grd           !
  real(crm_rknd) :: betdz         !
  real(crm_rknd) :: Ck            !
  real(crm_rknd) :: Ce            !
  real(crm_rknd) :: Ces           !
  real(crm_rknd) :: Ce1           !
  real(crm_rknd) :: Ce2           !
  real(crm_rknd) :: smix          !
  real(crm_rknd) :: Pr            ! Prandtl number
  real(crm_rknd) :: Cee           !
  real(crm_rknd) :: Cs            !
  real(crm_rknd) :: buoy_sgs      !
  real(crm_rknd) :: ratio         !
  real(crm_rknd) :: a_prod_sh     ! shear production of TKE
  real(crm_rknd) :: a_prod_bu     ! buoyant production of TKE
  real(crm_rknd) :: a_diss        ! TKE dissipation
  real(crm_rknd) :: lstarn        !
  real(crm_rknd) :: lstarp        !
  real(crm_rknd) :: bbb           !
  real(crm_rknd) :: omn           !
  real(crm_rknd) :: omp           !
  real(crm_rknd) :: qsatt         !
  real(crm_rknd) :: dqsat         !
  real(crm_rknd) :: cx            ! correction factor for eddy visc CFL criteria
  real(crm_rknd) :: cy            ! correction factor for eddy visc CFL criteria
  real(crm_rknd) :: cz            ! correction factor for eddy visc CFL criteria
  real(crm_rknd) :: tkmax         ! Maximum TKE (CFL limiter)

  integer :: i,j,k,icrm
  integer :: kc      ! = k+1
  integer :: kb      ! = k-1

  real(crm_rknd) :: tabs_interface    ! tabs interpolated to interfaces
  real(crm_rknd) :: qp_interface      ! qp   interpolated to interfaces
  real(crm_rknd) :: qtot_interface    ! qtot interpolated to interfaces
  real(crm_rknd) :: qsat_check        ! used to check for cloud
  real(crm_rknd) :: qctot             ! total condensate

  real(crm_rknd) :: tk_min_value      ! min value for eddy viscosity (TK)
  real(crm_rknd) :: tk_min_depth      ! near-surface depth to apply tk_min (meters)
  real(crm_rknd) :: tmp

  !$acc enter data create(def2,buoy_sgs_vert,a_prod_bu_vert) async(asyncid)

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

  tk_min_value = 0.05
  tk_min_depth = 500.

  Cs  = 0.15
  Ck  = 0.1
  Ce  = Ck**3/Cs**4
  Ces = Ce/0.7*3.0
  Pr  = 1.

  if(RUN3D) then
    call shear_prod3D(ncrms,def2)
  else
    call shear_prod2D(ncrms,def2)
  endif

  !!! initialize surface and top buoyancy flux to zero
  !$acc parallel loop collapse(3) copy(buoy_sgs_vert,a_prod_bu_vert) async(asyncid)
  do icrm = 1 , ncrms
    do j = 1 , ny
      do i = 1 , nx
        a_prod_bu_vert(i,j,0,icrm) = 0
        buoy_sgs_vert(i,j,0,icrm) = 0
        a_prod_bu_vert(i,j,nzm,icrm) = 0
        buoy_sgs_vert(i,j,nzm,icrm) = 0
      enddo
    enddo
  enddo

  !-----------------------------------------------------------------------
  ! compute SGS buoyancy flux at w-levels and SGS quantities in interior of domain
  !-----------------------------------------------------------------------
  !!! compute subgrid buoyancy flux assuming clear conditions
  !!! we will over-write this later if conditions are cloudy
  !$acc parallel loop collapse(4) copyin(tkh,tabs,bet,qv,qpi,qcl,qpl,t,qci,adzw,presi,dz) copy(buoy_sgs_vert,a_prod_bu_vert) async(asyncid)
  do icrm = 1 , ncrms
    do k = 1,nzm-1
      do j = 1,ny
        do i = 1,nx
          if (k.lt.nzm) then
            kc = k+1
            kb = k
          else
            kc = k
            kb = k-1
          end if
          !!! first compute subgrid buoyancy flux at interface above this level.
          !!! average betdz to w-levels
          betdz = 0.5*(bet(icrm,kc)+bet(icrm,kb))/dz(icrm)/adzw(icrm,k+1)

          !!! compute temperature of mixture between two grid levels if all cloud
          !!!   were evaporated and sublimated
          tabs_interface = &
               0.5*( tabs(icrm,i,j,kc) + fac_cond*qcl(icrm,i,j,kc) + fac_sub*qci(icrm,i,j,kc) &
                   + tabs(icrm,i,j,kb) + fac_cond*qcl(icrm,i,j,kb) + fac_sub*qci(icrm,i,j,kb) )

          !!! similarly for water vapor if all cloud evaporated/sublimated
          qtot_interface = &
               0.5*( qv(icrm,i,j,kc) + qcl(icrm,i,j,kc) + qci(icrm,i,j,kc) &
                   + qv(icrm,i,j,kb) + qcl(icrm,i,j,kb) + qci(icrm,i,j,kb) )

          qp_interface = 0.5*( qpl(icrm,i,j,kc) + qpi(icrm,i,j,kc) + qpl(icrm,i,j,kb) + qpi(icrm,i,j,kb) )

          bbb = 1.+epsv*qtot_interface - qp_interface
          buoy_sgs=betdz*( bbb*(t(icrm,i,j,kc)-t(icrm,i,j,kb)) &
               +epsv*tabs_interface* &
               (qv(icrm,i,j,kc)+qcl(icrm,i,j,kc)+qci(icrm,i,j,kc)-qv(icrm,i,j,kb)-qcl(icrm,i,j,kb)-qci(icrm,i,j,kb)) &
               +(bbb*fac_cond-tabs_interface)*(qpl(icrm,i,j,kc)-qpl(icrm,i,j,kb)) &
               +(bbb*fac_sub -tabs_interface)*(qpi(icrm,i,j,kc)-qpi(icrm,i,j,kb)) )

          buoy_sgs_vert(i,j,k,icrm) = buoy_sgs
          a_prod_bu_vert(i,j,k,icrm) = -0.5*(tkh(icrm,i,j,kc)+tkh(icrm,i,j,kb)+0.002)*buoy_sgs

          !-----------------------------------------------------------------------
          ! now go back and check for cloud
          !-----------------------------------------------------------------------

          !!! if there's any cloud in the grid cells above or below, check to see if
          !!! the mixture between the two levels is also cloudy
          qctot = qcl(icrm,i,j,kc)+qci(icrm,i,j,kc)+qcl(icrm,i,j,kb)+qci(icrm,i,j,kb)
          if(qctot .gt. 0.) then

            !!! figure out the fraction of condensate that's liquid
            omn = (qcl(icrm,i,j,kc)+qcl(icrm,i,j,kb))/(qctot+1.e-20)

            !!! compute temperature of mixture between two grid levels
            !!! if all cloud were evaporated and sublimated
            tabs_interface = &
                 0.5*( tabs(icrm,i,j,kc) + fac_cond*qcl(icrm,i,j,kc) + fac_sub*qci(icrm,i,j,kc) &
                     + tabs(icrm,i,j,kb) + fac_cond*qcl(icrm,i,j,kb) + fac_sub*qci(icrm,i,j,kb) )

            !!! similarly for total water (vapor + cloud) mixing ratio
            qtot_interface = &
                 0.5*( qv(icrm,i,j,kc) + qcl(icrm,i,j,kc) + qci(icrm,i,j,kc) &
                     + qv(icrm,i,j,kb) + qcl(icrm,i,j,kb) + qci(icrm,i,j,kb) )

            !!! compute saturation mixing ratio at this temperature
            qsat_check =      omn*qsatw_crm(tabs_interface,presi(k+1,icrm)) &
                        +(1.-omn)*qsati_crm(tabs_interface,presi(k+1,icrm))

            !!! check to see if the total water exceeds this saturation mixing ratio.
            !!! if so, apply the cloudy relations for subgrid buoyancy flux
            if(qtot_interface.gt.qsat_check) then
              !!! apply cloudy relations for buoyancy flux, use the liquid-ice breakdown computed above.
              lstarn = fac_cond+(1.-omn)*fac_fus
              !!! use the average values of T from the two levels to compute qsat, dqsat
              !!! and the multipliers for the subgrid buoyancy fluxes.  Note that the
              !!! interface is halfway between neighboring levels, so that the potential
              !!! energy cancels out.  This is approximate and neglects the effects of
              !!! evaporation/condensation with mixing.  Hopefully good enough.
              tabs_interface = 0.5*( tabs(icrm,i,j,kc) + tabs(icrm,i,j,kb) )

              qp_interface = 0.5*( qpl(icrm,i,j,kc) + qpi(icrm,i,j,kc) + qpl(icrm,i,j,kb) + qpi(icrm,i,j,kb) )

              dqsat =     omn*dtqsatw_crm(tabs_interface,presi(k+1,icrm)) + &
                     (1.-omn)*dtqsati_crm(tabs_interface,presi(k+1,icrm))
              qsatt =       omn*qsatw_crm(tabs_interface,presi(k+1,icrm)) + &
                       (1.-omn)*qsati_crm(tabs_interface,presi(k+1,icrm))

              !!! condensate loading term
              bbb = 1. + epsv*qsatt &
                   + qsatt - qtot_interface - qp_interface &
                   +1.61*tabs_interface*dqsat
              bbb = bbb / (1.+lstarn*dqsat)

              buoy_sgs = betdz*(bbb*(t(icrm,i,j,kc)-t(icrm,i,j,kb)) &
                   +(bbb*lstarn - (1.+lstarn*dqsat)*tabs_interface)* &
                   (qv(icrm,i,j,kc)+qcl(icrm,i,j,kc)+qci(icrm,i,j,kc)-qv(icrm,i,j,kb)-qcl(icrm,i,j,kb)-qci(icrm,i,j,kb)) &
                   + ( bbb*fac_cond-(1.+fac_cond*dqsat)*tabs(icrm,i,j,k) ) * ( qpl(icrm,i,j,kc)-qpl(icrm,i,j,kb) )  &
                   + ( bbb*fac_sub -(1.+fac_sub *dqsat)*tabs(icrm,i,j,k) ) * ( qpi(icrm,i,j,kc)-qpi(icrm,i,j,kb) ) )

              buoy_sgs_vert(i,j,k,icrm) = buoy_sgs
              a_prod_bu_vert(i,j,k,icrm) = -0.5*(tkh(icrm,i,j,kc)+tkh(icrm,i,j,kb)+0.002)*buoy_sgs
            end if ! if saturated at interface
          end if ! if either layer is cloudy
        end do ! i
      end do ! j
    enddo !k
  enddo !icrm

  !$acc parallel loop collapse(2) copy(tkelediss,tkesbshear,tkesbdiss,tkesbbuoy) async(asyncid)
  do icrm = 1 , ncrms
    do k = 1,nzm-1
      tkelediss(k,icrm)  = 0.
      tkesbdiss(k,icrm)  = 0.
      tkesbshear(k,icrm) = 0.
      tkesbbuoy(k,icrm)  = 0.
    enddo
  enddo

  !$acc parallel loop collapse(4) copyin(z,buoy_sgs_vert,def2,a_prod_bu_vert,dz,grdf_x,adz,grdf_y,adzw,grdf_z) copy(tkelediss,tkesbbuoy,tkesbshear,tkh,tk,tke,tkesbdiss) async(asyncid)
  do icrm = 1 , ncrms
    do k = 1,nzm-1
      do j = 1,ny
        do i = 1,nx
          grd = dz(icrm)*adz(icrm,k)
          Ce1 = Ce/0.7*0.19
          Ce2 = Ce/0.7*0.51
          !!! compute correction factors for eddy visc/cond not to acceed 3D stability
          cx = dx**2/dt/grdf_x(icrm,k)
          cy = dy**2/dt/grdf_y(icrm,k)
          cz = (dz(icrm)*min(adzw(icrm,k),adzw(icrm,k+1)))**2/dt/grdf_z(icrm,k)
          !!! maximum value of eddy visc/cond
          tkmax = 0.09/(1./cx+1./cy+1./cz)
          buoy_sgs = 0.5*( buoy_sgs_vert(i,j,k-1,icrm) + buoy_sgs_vert(i,j,k,icrm) )
          if(buoy_sgs.le.0.) then
            smix = grd
          else
            smix = min(grd,max(0.1*grd, sqrt(0.76*tk(icrm,i,j,k)/Ck/sqrt(buoy_sgs+1.e-10))))
          end if
          ratio = smix/grd
          Cee = Ce1+Ce2*ratio
          if(dosmagor) then
            tk(icrm,i,j,k) = sqrt(Ck**3/Cee*max(0._crm_rknd,def2(i,j,k,icrm)-Pr*buoy_sgs))*smix**2
#if defined( SP_TK_LIM )
              !!! put a hard lower limit on near-surface tk
              if ( z(icrm,k).lt.tk_min_depth ) then
                tk(icrm,i,j,k) = max( tk(icrm,i,j,k), tk_min_value )
              end if
#endif
            tke(icrm,i,j,k) = (tk(icrm,i,j,k)/(Ck*smix))**2
            a_prod_sh = (tk(icrm,i,j,k)+0.001)*def2(i,j,k,icrm)
            ! a_prod_bu=-(tk(icrm,i,j,k)+0.001)*Pr*buoy_sgs
            a_prod_bu = 0.5*( a_prod_bu_vert(i,j,k-1,icrm) + a_prod_bu_vert(i,j,k,icrm) )
            a_diss = a_prod_sh+a_prod_bu
          else
            tke(icrm,i,j,k) = max(real(0.,crm_rknd),tke(icrm,i,j,k))
            a_prod_sh = (tk(icrm,i,j,k)+0.001)*def2(i,j,k,icrm)
            a_prod_bu = 0.5*( a_prod_bu_vert(i,j,k-1,icrm) + a_prod_bu_vert(i,j,k,icrm) )
            !!! cap the diss rate (useful for large time steps)
            a_diss = min(tke(icrm,i,j,k)/(4.*dt),Cee/smix*tke(icrm,i,j,k)**1.5)
            tke(icrm,i,j,k) = max(real(0.,crm_rknd),tke(icrm,i,j,k)+dtn*(max(0._crm_rknd,a_prod_sh+a_prod_bu)-a_diss))
            tk(icrm,i,j,k)  = Ck*smix*sqrt(tke(icrm,i,j,k))
          end if
          tk(icrm,i,j,k)  = min(tk(icrm,i,j,k),tkmax)
          tkh(icrm,i,j,k) = Pr*tk(icrm,i,j,k)

          tmp = a_prod_sh/float(nx*ny)
          !$acc atomic update
          tkelediss(k,icrm)  = tkelediss(k,icrm) - tmp
          !$acc atomic update
          tkesbdiss(k,icrm)  = tkesbdiss(k,icrm) + a_diss
          !$acc atomic update
          tkesbshear(k,icrm) = tkesbshear(k,icrm)+ a_prod_sh
          !$acc atomic update
          tkesbbuoy(k,icrm)  = tkesbbuoy(k,icrm) + a_prod_bu
        end do ! i
      end do ! j
    end do ! k
  enddo !icrm

  !$acc exit data delete(def2,buoy_sgs_vert,a_prod_bu_vert) async(asyncid)

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

end subroutine tke_full

end module tke_full_mod
