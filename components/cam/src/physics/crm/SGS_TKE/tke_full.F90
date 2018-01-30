module tke_full_mod
  use shear_prod2D_mod
  use shear_prod3D_mod
  use sat_mod
  implicit none

contains

  subroutine tke_full(tkesbdiss, tkesbshear, tkesbbuoy, tke, tk, tkh, dimx1_d, dimx2_d, dimy1_d, dimy2_d, dosmagor, ncrms, icrm)

    !	this subroutine solves the TKE equation

    use vars
    use params
    implicit none
    integer, intent(in) :: ncrms,icrm
    logical :: dosmagor
    integer :: dimx1_d, dimx2_d, dimy1_d, dimy2_d
    real(crm_rknd) tkesbbuoy(nz), tkesbshear(nz), tkesbdiss(nz)
    real(crm_rknd) tke(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)   ! SGS TKE
    real(crm_rknd) tk  (dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm)  ! SGS eddy viscosity
    real(crm_rknd) tkh (dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm)  ! SGS eddy conductivity

    real(crm_rknd) def2(nx,ny,nzm)
    real(crm_rknd) grd,betdz,Ck,Ce,Ces,Ce1,Ce2,smix,Pr,Cee,Cs
    real(crm_rknd) buoy_sgs,ratio,a_prod_sh,a_prod_bu,a_diss
    real(crm_rknd) lstarn, lstarp, bbb, omn, omp
    real(crm_rknd) qsatt,dqsat
    integer i,j,k,kc,kb

    real(crm_rknd) tk_min_value   ! whannah - min value for eddy viscosity (TK)
    real(crm_rknd) tk_min_depth   ! whannah - near-surface depth to apply tk_min (meters)

    tk_min_value = 0.05
    tk_min_depth = 500.

    !call t_startf('tke_full')

    !Cs = 0.1944
    Cs = 0.15
    Ck=0.1
    Ce=Ck**3/Cs**4
    Ces=Ce/0.7*3.0

    if(RUN3D) then
      call shear_prod3D(def2, ncrms, icrm)
    else
      call shear_prod2D(def2, ncrms, icrm)
    endif

    do k=1,nzm
      kb=k-1
      kc=k+1

      grd=dz*adz(k)

      betdz=bet(k)/dz/(adzw(kc)+adzw(k))
      Ce1=Ce/0.7*0.19
      Ce2=Ce/0.7*0.51
      if(k.eq.1) then
        kb=1
        kc=2
        betdz=bet(k)/dz/adzw(kc)
        Ce1=Ces/0.7*0.19
        Ce2=Ces/0.7*0.51
      end if
      if(k.eq.nzm) then
        kb=nzm-1
        kc=nzm
        betdz=bet(k)/dz/adzw(k)
        Ce1=Ces/0.7*0.19
        Ce2=Ces/0.7*0.51
      end if
      tkelediss(k) = 0.
      tkesbdiss(k) = 0.
      tkesbshear(k)= 0.
      tkesbbuoy(k) = 0.
      do j=1,ny
        do i=1,nx
          !  SGS buoyancy flux

          !bloss: removed temperature diagnostics for omn.
          !         - use mauss weighted qsat, dqsat and latent heat for cloud
          !         - separate buoyancy contributions for precipitating water and ice.


          if(qcl(i,j,k)+qci(i,j,k) .gt. 0.) then

            omn = qcl(i,j,k)/(qcl(i,j,k)+qci(i,j,k)+1.e-20)
            lstarn = fac_cond+(1.-omn)*fac_fus

            dqsat = omn*dtqsatw_crm(tabs(icrm,i,j,k),pres(k))+ &
            (1.-omn)*dtqsati_crm(tabs(icrm,i,j,k),pres(k))
            qsatt = omn*qsatw_crm(tabs(icrm,i,j,k),pres(k))+(1.-omn)*qsati_crm(tabs(icrm,i,j,k),pres(k))
            bbb = 1. + epsv*qsatt-qcl(i,j,k)-qci(i,j,k) -qpl(i,j,k)-qpi(i,j,k)+1.61*tabs(icrm,i,j,k)*dqsat
            bbb = bbb / (1.+lstarn*dqsat)
            buoy_sgs=betdz*(bbb*(t(icrm,i,j,kc)-t(icrm,i,j,kb)) &
            +(bbb*lstarn - (1.+lstarn*dqsat)*tabs(icrm,i,j,k))* &
            (qv(i,j,kc)+qcl(i,j,kc)+qci(i,j,kc)-qv(i,j,kb)-qcl(i,j,kb)-qci(i,j,kb)) &
            + (bbb*fac_cond - (1.+fac_cond*dqsat)*tabs(icrm,i,j,k))*(qpl(i,j,kc)-qpl(i,j,kb))  &
            + (bbb*fac_sub  - (1.+fac_sub *dqsat)*tabs(icrm,i,j,k))*(qpi(i,j,kc)-qpi(i,j,kb)) )
            !bloss   +(bbb*lstarp - (1.+lstarp*dqsat)*tabs(icrm,i,j,k))* &
            !bloss            (qpl(i,j,kc)+qpi(i,j,kc)-qpl(i,j,kb)-qpi(i,j,kb)) )
          else

            bbb = 1.+epsv*qv(i,j,k)-qpl(i,j,k)-qpi(i,j,k)
            buoy_sgs=betdz*( bbb*(t(icrm,i,j,kc)-t(icrm,i,j,kb)) &
            +epsv*tabs(icrm,i,j,k)* &
            (qv(i,j,kc)+qcl(i,j,kc)+qci(i,j,kc)-qv(i,j,kb)-qcl(i,j,kb)-qci(i,j,kb)) &
            +(bbb*fac_cond-tabs(icrm,i,j,k))*(qpl(i,j,kc)-qpl(i,j,kb)) &
            +(bbb*fac_sub -tabs(icrm,i,j,k))*(qpi(i,j,kc)-qpi(i,j,kb)) )
            !bloss    +(bbb*lstarp-tabs(icrm,i,j,k))* &
            !bloss         (qpl(i,j,kc)+qpi(i,j,kc)-qpl(i,j,kb)-qpi(i,j,kb)) )
          end if

          if(buoy_sgs.le.0.) then
            smix=grd
          else
            smix=min(grd,max(0.1*grd, sqrt(0.76*tk(i,j,k)/Ck/sqrt(buoy_sgs+1.e-10))))
          end if


          ratio=smix/grd
          Pr=1.
          !   Pr=1. +2.*ratio
          Cee=Ce1+Ce2*ratio

          if(dosmagor) then

            tk(i,j,k)=sqrt(Ck**3/Cee*max(real(0.,crm_rknd),def2(i,j,k)-Pr*buoy_sgs))*smix**2

#ifdef SP_TK_LIM
            ! whannah - put a hard limit on near-surface tk
            if ( z(k).lt.tk_min_depth ) then
              tk(i,j,k) = max( tk(i,j,k), tk_min_value ) 
            end if
#endif

            ! tk(i,j,k)=sqrt(Ck**3/Cee*max(real(0.,crm_rknd),def2(i,j,k)-Pr*buoy_sgs))*smix**2
            tke(i,j,k) = (tk(i,j,k)/(Ck*smix))**2
            a_prod_sh=(tk(i,j,k)+0.001)*def2(i,j,k)
            a_prod_bu=-(tk(i,j,k)+0.001)*Pr*buoy_sgs
            a_diss=a_prod_sh+a_prod_bu

          else

            tke(i,j,k)=max(real(0.,crm_rknd),tke(i,j,k))
            a_prod_sh=(tk(i,j,k)+0.001)*def2(i,j,k)
            a_prod_bu=-(tk(i,j,k)+0.001)*Pr*buoy_sgs
            a_diss=min(tke(i,j,k)/(4.*dt),Cee/smix*tke(i,j,k)**1.5) ! cap the diss rate (useful for large time steps
            tke(i,j,k)=max(real(0.,crm_rknd),tke(i,j,k)+dtn*(max(real(0.,crm_rknd),a_prod_sh+a_prod_bu)-a_diss))
            tk(i,j,k)=Ck*smix*sqrt(tke(i,j,k))

          end if

          tkh(i,j,k)=Pr*tk(i,j,k)

          tkelediss(k) = tkelediss(k) - a_prod_sh
          tkesbdiss(k) = tkesbdiss(k) + a_diss
          tkesbshear(k)= tkesbshear(k)+ a_prod_sh
          tkesbbuoy(k) = tkesbbuoy(k) + a_prod_bu

        end do ! i
      end do ! j

      tkelediss(k) = tkelediss(k)/real(nx*ny,crm_rknd)


    end do ! k

    !call t_stopf('tke_full')

  end

end module tke_full_mod
