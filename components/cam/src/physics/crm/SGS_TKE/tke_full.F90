module tke_full_mod

  use shear_prod2D_mod
  use shear_prod3D_mod
  use sat_mod

  implicit none

contains

subroutine tke_full(ncrms,icrm,dimx1_d, dimx2_d, dimy1_d, dimy2_d,   &
                    grdf_x, grdf_y, grdf_z, dosmagor,     &
                    tkesbdiss, tkesbshear, tkesbbuoy,     &
                    tke, tk, tkh)
  !-----------------------------------------------------------------------
  ! Purpose: solve the TKE equation
  !-----------------------------------------------------------------------

  use vars
  use params
  
  implicit none
  integer, intent(in) :: ncrms,icrm
  !-----------------------------------------------------------------------
  !!! Interface Arguments
  integer       , intent(in)                 :: dimx1_d     ! scalar dimension parameter
  integer       , intent(in)                 :: dimx2_d     ! scalar dimension parameter
  integer       , intent(in)                 :: dimy1_d     ! scalar dimension parameter
  integer       , intent(in)                 :: dimy2_d     ! scalar dimension parameter
  real(crm_rknd), intent(in), dimension(nzm) :: grdf_x      ! grid length in x direction
  real(crm_rknd), intent(in), dimension(nzm) :: grdf_y      ! grid length in y direction
  real(crm_rknd), intent(in), dimension(nzm) :: grdf_z      ! grid length in z direction
  logical       , intent(in)                 :: dosmagor    ! flag for diagnostic smagorinsky scheme

  real(crm_rknd), intent(out), dimension(nz) :: tkesbdiss   ! TKE dissipation
  real(crm_rknd), intent(out), dimension(nz) :: tkesbshear  ! TKE production by shear
  real(crm_rknd), intent(out), dimension(nz) :: tkesbbuoy   ! TKE production by buoyancy

  real(crm_rknd), intent(out), dimension(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) :: tke   ! SGS TKE
  real(crm_rknd), intent(out), dimension(dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm) :: tk    ! SGS eddy viscosity
  real(crm_rknd), intent(out), dimension(dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm) :: tkh   ! SGS eddy conductivity
  
  !-----------------------------------------------------------------------
  !!! Local Variables
  real(crm_rknd), dimension(nx,ny,nzm) :: def2
  real(crm_rknd), dimension(nx,ny)     :: buoy_sgs_below
  real(crm_rknd), dimension(nx,ny)     :: buoy_sgs_above
  real(crm_rknd), dimension(nx,ny)     :: a_prod_bu_below
  real(crm_rknd), dimension(nx,ny)     :: a_prod_bu_above
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
  
  integer :: i,j,k
  integer :: kc      ! = k+1
  integer :: kb      ! = k-1

  real(crm_rknd) :: tabs_interface    ! tabs interpolated to interfaces
  real(crm_rknd) :: qp_interface      ! qp   interpolated to interfaces
  real(crm_rknd) :: qtot_interface    ! qtot interpolated to interfaces
  real(crm_rknd) :: qsat_check        ! used to check for cloud
  real(crm_rknd) :: qctot             ! total condensate

  real(crm_rknd) :: tk_min_value      ! min value for eddy viscosity (TK)
  real(crm_rknd) :: tk_min_depth      ! near-surface depth to apply tk_min (meters)

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
    call shear_prod3D(ncrms,icrm,def2)
  else
    call shear_prod2D(ncrms,icrm,def2)
  endif
  
  !!! initialize surface buoyancy flux to zero
  a_prod_bu_below(:,:) = real(0.0,crm_rknd)
  buoy_sgs_below(:,:) = real(0.0,crm_rknd)

  !-----------------------------------------------------------------------
  ! compute SGS buoyancy flux at w-levels and SGS quantities in interior of domain
  !-----------------------------------------------------------------------
  do k = 1,nzm

    if (k.lt.nzm) then
      kc = k+1
      kb = k
    else
      kc = k
      kb = k-1
    end if

    !!! first compute subgrid buoyancy flux at interface above this level.

    !!! average betdz to w-levels
    betdz = 0.5*(bet(kc,icrm)+bet(kb,icrm))/dz/adzw(k+1)

    !!! compute subgrid buoyancy flux assuming clear conditions
    !!! we will over-write this later if conditions are cloudy
    do j = 1,ny
      do i = 1,nx

        !!! compute temperature of mixture between two grid levels if all cloud
        !!!   were evaporated and sublimated
        tabs_interface = &
             0.5*( tabs(i,j,kc,icrm) + fac_cond*qcl(i,j,kc,icrm) + fac_sub*qci(i,j,kc,icrm) &
                 + tabs(i,j,kb,icrm) + fac_cond*qcl(i,j,kb,icrm) + fac_sub*qci(i,j,kb,icrm) )

        !!! similarly for water vapor if all cloud evaporated/sublimated
        qtot_interface = &
             0.5*( qv(i,j,kc,icrm) + qcl(i,j,kc,icrm) + qci(i,j,kc,icrm) &
                 + qv(i,j,kb,icrm) + qcl(i,j,kb,icrm) + qci(i,j,kb,icrm) )

        qp_interface = 0.5*( qpl(i,j,kc,icrm) + qpi(i,j,kc,icrm) + qpl(i,j,kb,icrm) + qpi(i,j,kb,icrm) )

        bbb = 1.+epsv*qtot_interface - qp_interface
        buoy_sgs=betdz*( bbb*(t(i,j,kc)-t(i,j,kb)) &
             +epsv*tabs_interface* &
             (qv(i,j,kc,icrm)+qcl(i,j,kc,icrm)+qci(i,j,kc,icrm)-qv(i,j,kb,icrm)-qcl(i,j,kb,icrm)-qci(i,j,kb,icrm)) &
             +(bbb*fac_cond-tabs_interface)*(qpl(i,j,kc,icrm)-qpl(i,j,kb,icrm)) &
             +(bbb*fac_sub -tabs_interface)*(qpi(i,j,kc,icrm)-qpi(i,j,kb,icrm)) )

        buoy_sgs_above(i,j) = buoy_sgs
        a_prod_bu_above(i,j) = -0.5*(tkh(i,j,kc)+tkh(i,j,kb)+0.002)*buoy_sgs

      end do ! i
    end do ! j

    !-----------------------------------------------------------------------
    ! now go back and check for cloud
    !-----------------------------------------------------------------------
    do j = 1,ny
      do i = 1,nx
        !!! if there's any cloud in the grid cells above or below, check to see if
        !!! the mixture between the two levels is also cloudy
        qctot = qcl(i,j,kc,icrm)+qci(i,j,kc,icrm)+qcl(i,j,kb,icrm)+qci(i,j,kb,icrm)
        if(qctot .gt. 0.) then
        
          !!! figure out the fraction of condensate that's liquid
          omn = (qcl(i,j,kc,icrm)+qcl(i,j,kb,icrm))/(qctot+1.e-20)
          
          !!! compute temperature of mixture between two grid levels 
          !!! if all cloud were evaporated and sublimated
          tabs_interface = &
               0.5*( tabs(i,j,kc,icrm) + fac_cond*qcl(i,j,kc,icrm) + fac_sub*qci(i,j,kc,icrm) &
                   + tabs(i,j,kb,icrm) + fac_cond*qcl(i,j,kb,icrm) + fac_sub*qci(i,j,kb,icrm) )

          !!! similarly for total water (vapor + cloud) mixing ratio
          qtot_interface = &
               0.5*( qv(i,j,kc,icrm) + qcl(i,j,kc,icrm) + qci(i,j,kc,icrm) &
                   + qv(i,j,kb,icrm) + qcl(i,j,kb,icrm) + qci(i,j,kb,icrm) )

          !!! compute saturation mixing ratio at this temperature
          qsat_check =      omn*qsatw_crm(tabs_interface,presi(k+1)) &
                      +(1.-omn)*qsati_crm(tabs_interface,presi(k+1))

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
            tabs_interface = 0.5*( tabs(i,j,kc,icrm) + tabs(i,j,kb,icrm) )

            qp_interface = 0.5*( qpl(i,j,kc,icrm) + qpi(i,j,kc,icrm) + qpl(i,j,kb,icrm) + qpi(i,j,kb,icrm) )

            dqsat =     omn*dtqsatw_crm(tabs_interface,presi(k+1)) + &
                   (1.-omn)*dtqsati_crm(tabs_interface,presi(k+1))
            qsatt =       omn*qsatw_crm(tabs_interface,presi(k+1)) + &
                     (1.-omn)*qsati_crm(tabs_interface,presi(k+1))

            !!! condensate loading term
            bbb = 1. + epsv*qsatt &
                 + qsatt - qtot_interface - qp_interface & 
                 +1.61*tabs_interface*dqsat
            bbb = bbb / (1.+lstarn*dqsat)

            buoy_sgs = betdz*(bbb*(t(i,j,kc)-t(i,j,kb)) &
                 +(bbb*lstarn - (1.+lstarn*dqsat)*tabs_interface)* &
                 (qv(i,j,kc,icrm)+qcl(i,j,kc,icrm)+qci(i,j,kc,icrm)-qv(i,j,kb,icrm)-qcl(i,j,kb,icrm)-qci(i,j,kb,icrm)) & 
                 + ( bbb*fac_cond-(1.+fac_cond*dqsat)*tabs(i,j,k,icrm) ) * ( qpl(i,j,kc,icrm)-qpl(i,j,kb,icrm) )  &
                 + ( bbb*fac_sub -(1.+fac_sub *dqsat)*tabs(i,j,k,icrm) ) * ( qpi(i,j,kc,icrm)-qpi(i,j,kb,icrm) ) )

            buoy_sgs_above(i,j) = buoy_sgs
            a_prod_bu_above(i,j) = -0.5*(tkh(i,j,kc)+tkh(i,j,kb)+0.002)*buoy_sgs

            ! buoy_sgs_above(i,j) = -0.5*(tkh(i,j,kc)+tkh(i,j,kb)*buoy_sgs    !bloss: Should we add the offset to tkh or not??

          end if ! if saturated at interface

        end if ! if either layer is cloudy

      end do ! i
    end do ! j

    if(k.eq.nzm) then
      !!! zero out the subgrid buoyancy flux at the top level
      buoy_sgs_above(:,:)  = 0.
      a_prod_bu_above(:,:) = 0.
    end if

    grd = dz*adz(k)

    Ce1 = Ce/0.7*0.19
    Ce2 = Ce/0.7*0.51


    tkelediss(k)  = 0.
    tkesbdiss(k)  = 0.
    tkesbshear(k) = 0.
    tkesbbuoy(k)  = 0.

    !!! compute correction factors for eddy visc/cond not to acceed 3D stability
    cx = dx**2/dt/grdf_x(k)
    cy = dy**2/dt/grdf_y(k)
    cz = (dz*min(adzw(k),adzw(k+1)))**2/dt/grdf_z(k)

    !!! maximum value of eddy visc/cond
    tkmax = 0.09/(1./cx+1./cy+1./cz)  

    do j = 1,ny
      do i = 1,nx

        buoy_sgs = 0.5*( buoy_sgs_below(i,j) + buoy_sgs_above(i,j) )

        if(buoy_sgs.le.0.) then
          smix = grd
        else
          smix = min(grd,max(0.1*grd, sqrt(0.76*tk(i,j,k)/Ck/sqrt(buoy_sgs+1.e-10))))
        end if


        ratio = smix/grd
        Cee = Ce1+Ce2*ratio

        if(dosmagor) then

          tk(i,j,k) = sqrt(Ck**3/Cee*max(0.,def2(i,j,k)-Pr*buoy_sgs))*smix**2

#if defined( SP_TK_LIM )
            !!! put a hard lower limit on near-surface tk
            if ( z(k).lt.tk_min_depth ) then
              tk(i,j,k) = max( tk(i,j,k), tk_min_value ) 
            end if
#endif

          tke(i,j,k) = (tk(i,j,k)/(Ck*smix))**2
          a_prod_sh = (tk(i,j,k)+0.001)*def2(i,j,k)
          ! a_prod_bu=-(tk(i,j,k)+0.001)*Pr*buoy_sgs
          a_prod_bu = 0.5*( a_prod_bu_below(i,j) + a_prod_bu_above(i,j) )
          a_diss = a_prod_sh+a_prod_bu

        else

          tke(i,j,k) = max(real(0.,crm_rknd),tke(i,j,k))
          a_prod_sh = (tk(i,j,k)+0.001)*def2(i,j,k)
          a_prod_bu = 0.5*( a_prod_bu_below(i,j) + a_prod_bu_above(i,j) )
          !!! cap the diss rate (useful for large time steps)
          a_diss = min(tke(i,j,k)/(4.*dt),Cee/smix*tke(i,j,k)**1.5)               
          tke(i,j,k) = max(real(0.,crm_rknd),tke(i,j,k)+dtn*(max(0.,a_prod_sh+a_prod_bu)-a_diss))
          tk(i,j,k)  = Ck*smix*sqrt(tke(i,j,k))

        end if

        tk(i,j,k)  = min(tk(i,j,k),tkmax)
        tkh(i,j,k) = Pr*tk(i,j,k)

        tkelediss(k)  = tkelediss(k) - a_prod_sh
        tkesbdiss(k)  = tkesbdiss(k) + a_diss
        tkesbshear(k) = tkesbshear(k)+ a_prod_sh
        tkesbbuoy(k)  = tkesbbuoy(k) + a_prod_bu

      end do ! i
    end do ! j

    tkelediss(k) = tkelediss(k)/float(nx*ny)

    buoy_sgs_below(:,:) = buoy_sgs_above(:,:)
    a_prod_bu_below(:,:) = a_prod_bu_above(:,:)

  end do ! k

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

end subroutine tke_full

end module tke_full_mod
