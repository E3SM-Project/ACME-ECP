module tke_full_mod
  use shear_prod2D_mod
  use shear_prod3D_mod
  use sat_mod
  implicit none

contains

  subroutine tke_full(tkesbdiss, tkesbshear, tkesbbuoy, tke, tk, tkh, dimx1_d, dimx2_d, dimy1_d, dimy2_d, dosmagor, ncrms)

    !	this subroutine solves the TKE equation

    use vars
    use params
    implicit none
    integer, intent(in) :: ncrms
    logical :: dosmagor
    integer :: dimx1_d, dimx2_d, dimy1_d, dimy2_d
    real(crm_rknd) tkesbbuoy(ncrms,nz), tkesbshear(ncrms,nz), tkesbdiss(ncrms,nz)
    real(crm_rknd) tke (ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)   ! SGS TKE
    real(crm_rknd) tk  (ncrms,dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm)  ! SGS eddy viscosity
    real(crm_rknd) tkh (ncrms,dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm)  ! SGS eddy conductivity

    real(crm_rknd) def2(ncrms,nx,ny,nzm)
    real(crm_rknd) grd,betdz,Ck,Ce,Ces,Ce1,Ce2,smix,Pr,Cee,Cs
    real(crm_rknd) buoy_sgs,ratio,a_prod_sh,a_prod_bu,a_diss
    real(crm_rknd) lstarn, lstarp, bbb, omn, omp
    real(crm_rknd) qsatt,dqsat
    integer i,j,k,kc,kb,icrm

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
      call shear_prod3D(def2, ncrms)
    else
      call shear_prod2D(def2, ncrms)
    endif

    do icrm = 1 , ncrms

    do k=1,nzm
      kb=k-1
      kc=k+1

      grd=dz(icrm)*adz(icrm,k)

      betdz=bet(icrm,k)/dz(icrm)/(adzw(icrm,kc)+adzw(icrm,k))
      Ce1=Ce/0.7*0.19
      Ce2=Ce/0.7*0.51
      if(k.eq.1) then
        kb=1
        kc=2
        betdz=bet(icrm,k)/dz(icrm)/adzw(icrm,kc)
        Ce1=Ces/0.7*0.19
        Ce2=Ces/0.7*0.51
      end if
      if(k.eq.nzm) then
        kb=nzm-1
        kc=nzm
        betdz=bet(icrm,k)/dz(icrm)/adzw(icrm,k)
        Ce1=Ces/0.7*0.19
        Ce2=Ces/0.7*0.51
      end if
      tkelediss(icrm,k) = 0.
      tkesbdiss(icrm,k) = 0.
      tkesbshear(icrm,k)= 0.
      tkesbbuoy(icrm,k) = 0.
      do j=1,ny
        do i=1,nx
          !  SGS buoyancy flux

          !bloss: removed temperature diagnostics for omn.
          !         - use mauss weighted qsat, dqsat and latent heat for cloud
          !         - separate buoyancy contributions for precipitating water and ice.


          if(qcl(icrm,i,j,k)+qci(icrm,i,j,k) .gt. 0.) then

            omn = qcl(icrm,i,j,k)/(qcl(icrm,i,j,k)+qci(icrm,i,j,k)+1.e-20)
            lstarn = fac_cond+(1.-omn)*fac_fus

            dqsat = omn*dtqsatw_crm(tabs(icrm,i,j,k),pres(icrm,k))+ &
            (1.-omn)*dtqsati_crm(tabs(icrm,i,j,k),pres(icrm,k))
            qsatt = omn*qsatw_crm(tabs(icrm,i,j,k),pres(icrm,k))+(1.-omn)*qsati_crm(tabs(icrm,i,j,k),pres(icrm,k))
            bbb = 1. + epsv*qsatt-qcl(icrm,i,j,k)-qci(icrm,i,j,k) -qpl(icrm,i,j,k)-qpi(icrm,i,j,k)+1.61*tabs(icrm,i,j,k)*dqsat
            bbb = bbb / (1.+lstarn*dqsat)
            buoy_sgs=betdz*(bbb*(t(icrm,i,j,kc)-t(icrm,i,j,kb)) &
            +(bbb*lstarn - (1.+lstarn*dqsat)*tabs(icrm,i,j,k))* &
            (qv(icrm,i,j,kc)+qcl(icrm,i,j,kc)+qci(icrm,i,j,kc)-qv(icrm,i,j,kb)-qcl(icrm,i,j,kb)-qci(icrm,i,j,kb)) &
            + (bbb*fac_cond - (1.+fac_cond*dqsat)*tabs(icrm,i,j,k))*(qpl(icrm,i,j,kc)-qpl(icrm,i,j,kb))  &
            + (bbb*fac_sub  - (1.+fac_sub *dqsat)*tabs(icrm,i,j,k))*(qpi(icrm,i,j,kc)-qpi(icrm,i,j,kb)) )
            !bloss   +(bbb*lstarp - (1.+lstarp*dqsat)*tabs(icrm,i,j,k))* &
            !bloss            (qpl(icrm,i,j,kc)+qpi(icrm,i,j,kc)-qpl(icrm,i,j,kb)-qpi(icrm,i,j,kb)) )
          else

            bbb = 1.+epsv*qv(icrm,i,j,k)-qpl(icrm,i,j,k)-qpi(icrm,i,j,k)
            buoy_sgs=betdz*( bbb*(t(icrm,i,j,kc)-t(icrm,i,j,kb)) &
            +epsv*tabs(icrm,i,j,k)* &
            (qv(icrm,i,j,kc)+qcl(icrm,i,j,kc)+qci(icrm,i,j,kc)-qv(icrm,i,j,kb)-qcl(icrm,i,j,kb)-qci(icrm,i,j,kb)) &
            +(bbb*fac_cond-tabs(icrm,i,j,k))*(qpl(icrm,i,j,kc)-qpl(icrm,i,j,kb)) &
            +(bbb*fac_sub -tabs(icrm,i,j,k))*(qpi(icrm,i,j,kc)-qpi(icrm,i,j,kb)) )
            !bloss    +(bbb*lstarp-tabs(icrm,i,j,k))* &
            !bloss         (qpl(icrm,i,j,kc)+qpi(icrm,i,j,kc)-qpl(icrm,i,j,kb)-qpi(icrm,i,j,kb)) )
          end if

          if(buoy_sgs.le.0.) then
            smix=grd
          else
            smix=min(grd,max(0.1*grd, sqrt(0.76*tk(icrm,i,j,k)/Ck/sqrt(buoy_sgs+1.e-10))))
          end if


          ratio=smix/grd
          Pr=1.
          !   Pr=1. +2.*ratio
          Cee=Ce1+Ce2*ratio

          if(dosmagor) then

            tk(icrm,i,j,k)=sqrt(Ck**3/Cee*max(real(0.,crm_rknd),def2(icrm,i,j,k)-Pr*buoy_sgs))*smix**2

#ifdef SP_TK_LIM
            ! whannah - put a hard limit on near-surface tk
            if ( z(icrm,k).lt.tk_min_depth ) then
              tk(icrm,i,j,k) = max( tk(icrm,i,j,k), tk_min_value )
            end if
#endif

            ! tk(icrm,i,j,k)=sqrt(Ck**3/Cee*max(real(0.,crm_rknd),def2(icrm,i,j,k)-Pr*buoy_sgs))*smix**2
            tke(icrm,i,j,k) = (tk(icrm,i,j,k)/(Ck*smix))**2
            a_prod_sh=(tk(icrm,i,j,k)+0.001)*def2(icrm,i,j,k)
            a_prod_bu=-(tk(icrm,i,j,k)+0.001)*Pr*buoy_sgs
            a_diss=a_prod_sh+a_prod_bu

          else

            tke(icrm,i,j,k)=max(real(0.,crm_rknd),tke(icrm,i,j,k))
            a_prod_sh=(tk(icrm,i,j,k)+0.001)*def2(icrm,i,j,k)
            a_prod_bu=-(tk(icrm,i,j,k)+0.001)*Pr*buoy_sgs
            a_diss=min(tke(icrm,i,j,k)/(4.*dt),Cee/smix*tke(icrm,i,j,k)**1.5) ! cap the diss rate (useful for large time steps
            tke(icrm,i,j,k)=max(real(0.,crm_rknd),tke(icrm,i,j,k)+dtn*(max(real(0.,crm_rknd),a_prod_sh+a_prod_bu)-a_diss))
            tk(icrm,i,j,k)=Ck*smix*sqrt(tke(icrm,i,j,k))

          end if

          tkh(icrm,i,j,k)=Pr*tk(icrm,i,j,k)

          tkelediss(icrm,k) = tkelediss(icrm,k) - a_prod_sh
          tkesbdiss(icrm,k) = tkesbdiss(icrm,k) + a_diss
          tkesbshear(icrm,k)= tkesbshear(icrm,k)+ a_prod_sh
          tkesbbuoy(icrm,k) = tkesbbuoy(icrm,k) + a_prod_bu

        end do ! i
      end do ! j

      tkelediss(icrm,k) = tkelediss(icrm,k)/real(nx*ny,crm_rknd)


    end do ! k
  enddo

    !call t_stopf('tke_full')

  end subroutine tke_full

end module tke_full_mod
