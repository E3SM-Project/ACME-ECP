module tke_full_mod
  use shear_prod2D_mod
  use shear_prod3D_mod
  use sat_mod
  implicit none

contains

  subroutine tke_full(tkesbdiss, tkesbshear, tkesbbuoy, tke, tk, tkh, &
                      dimx1_d, dimx2_d, dimy1_d, dimy2_d, &
                      grdf_x, grdf_y, grdf_z, &
                      dosmagor)
    !-----------------------------------------------------------------------
    ! Purpose: solve the TKE equation
    !-----------------------------------------------------------------------

    use vars
    use params
    
    implicit none

    !-----------------------------------------------------------------------
    !!! Interface Arguments
    real(crm_rknd), intent(out), dimension(nz) :: tkesbdiss
    real(crm_rknd), intent(out), dimension(nz) :: tkesbshear
    real(crm_rknd), intent(out), dimension(nz) :: tkesbbuoy
    real(crm_rknd), intent(out), dimension(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) :: tke   ! SGS TKE
    real(crm_rknd), intent(out), dimension(dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm) :: tk    ! SGS eddy viscosity
    real(crm_rknd), intent(out), dimension(dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm) :: tkh   ! SGS eddy conductivity
    integer       , intent(in)                 :: dimx1_d
    integer       , intent(in)                 :: dimx2_d
    integer       , intent(in)                 :: dimy1_d
    integer       , intent(in)                 :: dimy2_d
    real(crm_rknd), intent(in), dimension(nzm) :: grdf_x
    real(crm_rknd), intent(in), dimension(nzm) :: grdf_y
    real(crm_rknd), intent(in), dimension(nzm) :: grdf_z
    logical       , intent(in) :: dosmagor
    
    !-----------------------------------------------------------------------
    !!! Local Variables
    real(crm_rknd), dimension(nx,ny,nzm) :: def2
    real(crm_rknd), dimension(nx,ny)     :: buoy_sgs_below
    real(crm_rknd), dimension(nx,ny)     :: buoy_sgs_above
    real(crm_rknd), dimension(nx,ny)     :: a_prod_bu_below
    real(crm_rknd), dimension(nx,ny)     :: a_prod_bu_above
    real(crm_rknd) grd,betdz,Ck,Ce,Ces,Ce1,Ce2,smix,Pr,Cee,Cs
    real(crm_rknd) buoy_sgs,ratio,a_prod_sh,a_prod_bu,a_diss
    real(crm_rknd) lstarn, lstarp, bbb, omn, omp
    real(crm_rknd) qsatt,dqsat
    real(crm_rknd) cx,cy,cz,tkmax
    integer i,j,k,kc,kb

    real(crm_rknd) tabs_interface, qp_interface, qtot_interface, qsat_check, qctot

    real(crm_rknd) tk_min_value   ! min value for eddy viscosity (TK)
    real(crm_rknd) tk_min_depth   ! near-surface depth to apply tk_min (meters)

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
      call shear_prod3D(def2)
    else
      call shear_prod2D(def2)
    endif

    !!! whannah - not sure what to do about this first time step block - brought over from SAM6.11
    !!! the SP_TK_LIM option already guarantees non-zero tke at the surface - so I don't think we need Blossey's code

    ! if((nstep.eq.1).AND.(icycle.eq.1)) then
    !   !bloss(2016-05-09): At start of simulation, make sure that subgrid TKE
    !   !   is non-zero at surface if surface buoyancy fluxes are positive.
    !   !   If they are, compute the TKE implied by local equilibrium between
    !   !   turbulence production by the surface fluxes and dissipation.  Take
    !   !   the initial TKE in the lowest level to be the larger of that and the 
    !   !   initial value.  Since the present values of TKE, eddy viscosity and eddy
    !   !   diffusivity are used in computing the new values, it is
    !   !   important for them not to be zero initially if the buoyancy
    !   !   flux is non-zero initially.

    !   k = 1
    !   do j = 1,ny
    !     do i = 1,nx
    !       !bloss: compute suface buoyancy flux
    !       bbb = 1.+epsv*qv(i,j,k)
    !       a_prod_bu_below(i,j) = bbb*bet(k)*fluxbt(i,j) + bet(k)*epsv*(t00+sstxy(i,j))*fluxbq(i,j) 
      
    !       grd=dz*adz(k)
    !       Pr=1. 
    !       Ce1=Ce/0.7*0.19
    !       Ce2=Ce/0.7*0.51
    !       Cee=Ce1+Ce2

    !       ! Choose the subgrid TKE to be the larger of the initial value or
    !       !   that which satisfies local equilibrium, buoyant production = dissipation
    !       !   or a_prod_bu = Cee/grd * tke^(3/2).
    !       !   NOTE: We're ignoring shear production here.
    !       tke(i,j,1) = MAX( tke(i,j,1), &
    !                         ( grd/Cee * MAX(1.e-20, 0.5*a_prod_bu_below(i,j) ) )**(2./3.) )

    !       ! eddy viscosity = Ck*grd * sqrt(tke) --- analogous for Smagorinksy.
    !       tk(i,j,1) = Ck*grd * SQRT( tke(i,j,1) )

    !       ! eddy diffusivity = Pr * eddy viscosity
    !       tkh(i,j,1) = Pr*tk(i,j,1)
    !     end do
    !   end do

    ! end if ! if(nstep.eq.1.AND.icycle.eq.1)

    !!! whannah - also commenting this out because we don't deal with surface fluxes directly for SP
    !!! instead we'll just initialize surface buoyancy flux to zero
    a_prod_bu_below(:,:) = real(0.0,crm_rknd)
    buoy_sgs_below(:,:) = real(0.0,crm_rknd)

    !-----------------------------------------------------------------------
    ! compute subgrid buoyancy flux at w-levels, starting with surface buoyancy flux
    !-----------------------------------------------------------------------
    ! k = 1
    ! do j = 1,ny
    !   do i = 1,nx
    !     !bloss: Use SST along with surface vapor mixing ratio.  This is slightly inconsistent, 
    !     !  but the error is small, and it's cheaper than another saturation mixing ratio computation.
    !     bbb = 1.+epsv*qv(i,j,k)
    !     a_prod_bu_below(i,j) = bbb*bet(k)*fluxbt(i,j) + bet(k)*epsv*(t00+sstxy(i,j))*fluxbq(i,j) 
    !     ! back buoy_sgs out from buoyancy flux, a_prod_bu = - (tkh(i,j,k)+0.001)*buoy_sgs
    !     buoy_sgs_below(i,j) = - a_prod_bu_below(i,j)/(tkh(i,j,k)+0.001)
    !   end do
    ! end do

    !-----------------------------------------------------------------------
    ! now compute SGS quantities in interior of domain.
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
      betdz = 0.5*(bet(kc)+bet(kb))/dz/adzw(k+1)

      !!! compute subgrid buoyancy flux assuming clear conditions
      !!! we will over-write this later if conditions are cloudy
      do j = 1,ny
        do i = 1,nx

          !!! compute temperature of mixture between two grid levels if all cloud
          !!!   were evaporated and sublimated
          tabs_interface = &
               0.5*( tabs(i,j,kc) + fac_cond*qcl(i,j,kc) + fac_sub*qci(i,j,kc) &
                   + tabs(i,j,kb) + fac_cond*qcl(i,j,kb) + fac_sub*qci(i,j,kb) )

          !!! similarly for water vapor if all cloud evaporated/sublimated
          qtot_interface = &
               0.5*( qv(i,j,kc) + qcl(i,j,kc) + qci(i,j,kc) &
                   + qv(i,j,kb) + qcl(i,j,kb) + qci(i,j,kb) )

          qp_interface = 0.5*( qpl(i,j,kc) + qpi(i,j,kc) + qpl(i,j,kb) + qpi(i,j,kb) )

          bbb = 1.+epsv*qtot_interface - qp_interface
          buoy_sgs=betdz*( bbb*(t(i,j,kc)-t(i,j,kb)) &
               +epsv*tabs_interface* &
               (qv(i,j,kc)+qcl(i,j,kc)+qci(i,j,kc)-qv(i,j,kb)-qcl(i,j,kb)-qci(i,j,kb)) &
               +(bbb*fac_cond-tabs_interface)*(qpl(i,j,kc)-qpl(i,j,kb)) &
               +(bbb*fac_sub -tabs_interface)*(qpi(i,j,kc)-qpi(i,j,kb)) )

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
          qctot = qcl(i,j,kc)+qci(i,j,kc)+qcl(i,j,kb)+qci(i,j,kb)
          if(qctot .gt. 0.) then
          
            !!! figure out the fraction of condensate that's liquid
            omn = (qcl(i,j,kc)+qcl(i,j,kb))/(qctot+1.e-20)
            
            !!! compute temperature of mixture between two grid levels if all cloud were evaporated and sublimated
            tabs_interface = &
                 0.5*( tabs(i,j,kc) + fac_cond*qcl(i,j,kc) + fac_sub*qci(i,j,kc) &
                     + tabs(i,j,kb) + fac_cond*qcl(i,j,kb) + fac_sub*qci(i,j,kb) )

            !!! similarly for total water (vapor + cloud) mixing ratio
            qtot_interface = &
                 0.5*( qv(i,j,kc) + qcl(i,j,kc) + qci(i,j,kc) &
                     + qv(i,j,kb) + qcl(i,j,kb) + qci(i,j,kb) )

            !!! compute saturation mixing ratio at this temperature
            qsat_check = omn*qsatw_crm(tabs_interface,presi(k+1))+(1.-omn)*qsati_crm(tabs_interface,presi(k+1))

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
              tabs_interface = 0.5*( tabs(i,j,kc) + tabs(i,j,kb) )

              qp_interface = 0.5*( qpl(i,j,kc) + qpi(i,j,kc) + qpl(i,j,kb) + qpi(i,j,kb) )

              dqsat = omn*dtqsatw_crm(tabs_interface,presi(k+1))+ &
                   (1.-omn)*dtqsati_crm(tabs_interface,presi(k+1))
              qsatt = omn*qsatw_crm(tabs_interface,presi(k+1))+(1.-omn)*qsati_crm(tabs_interface,presi(k+1))

              bbb = 1. + epsv*qsatt &
                   + qsatt - qtot_interface - qp_interface & ! condensate loading term
                   +1.61*tabs_interface*dqsat
              bbb = bbb / (1.+lstarn*dqsat)

              buoy_sgs = betdz*(bbb*(t(i,j,kc)-t(i,j,kb)) &
                   +(bbb*lstarn - (1.+lstarn*dqsat)*tabs_interface)* &
                   (qv(i,j,kc)+qcl(i,j,kc)+qci(i,j,kc)-qv(i,j,kb)-qcl(i,j,kb)-qci(i,j,kb)) & 
                   + (bbb*fac_cond - (1.+fac_cond*dqsat)*tabs(i,j,k))*(qpl(i,j,kc)-qpl(i,j,kb))  &
                   + (bbb*fac_sub  - (1.+fac_sub *dqsat)*tabs(i,j,k))*(qpi(i,j,kc)-qpi(i,j,kb)) )

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


      tkelediss(k) = 0.
      tkesbdiss(k) = 0.
      tkesbshear(k)= 0.
      tkesbbuoy(k) = 0.

      !!! compute correction factors for eddy visc/cond not to acceed 3D stability
      cx = dx**2/dt/grdf_x(k)
      cy = dy**2/dt/grdf_y(k)
      cz = (dz*min(adzw(k),adzw(k+1)))**2/dt/grdf_z(k)
      tkmax = 0.09/(1./cx+1./cy+1./cz)  ! maximum value of eddy visc/cond

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
