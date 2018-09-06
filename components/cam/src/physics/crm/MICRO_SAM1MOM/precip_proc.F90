module precip_proc_mod
  implicit none

contains

  subroutine precip_proc(ncrms,icrm,qpsrc,qpevp,qp,q,qn)

    use vars
    use micro_params
    use params
    use sat_mod

    implicit none
    integer, intent(in) :: ncrms,icrm
    real(crm_rknd) q(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)   ! total nonprecipitating water
    real(crm_rknd) qp(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)  ! total precipitating water
    real(crm_rknd) qn(nx,ny,nzm,ncrms)  ! cloud condensate (liquid + ice)
    real(crm_rknd) qpsrc(nz,ncrms)  ! source of precipitation microphysical processes
    real(crm_rknd) qpevp(nz,ncrms)  ! sink of precipitating water due to evaporation

    integer i,j,k
    real(crm_rknd) autor, autos, accrr, accris, accrcs, accrig, accrcg
    real(crm_rknd) dq, omn, omp, omg, qsatt
    real(crm_rknd) pows1, pows2, powg1, powg2, powr1, powr2, tmp
    real(crm_rknd) qii, qcc, qrr, qss, qgg

    powr1 = (3 + b_rain) / 4.
    powr2 = (5 + b_rain) / 8.
    pows1 = (3 + b_snow) / 4.
    pows2 = (5 + b_snow) / 8.
    powg1 = (3 + b_grau) / 4.
    powg2 = (5 + b_grau) / 8.

    !call t_startf ('precip_proc')

    do k=1,nzm
      qpsrc(k,icrm)=0.
      qpevp(k,icrm)=0.
      do j=1,ny
        do i=1,nx

          !-------     Autoconversion/accretion

          if(qn(i,j,k,icrm)+qp(i,j,k).gt.0.) then


            omn = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(i,j,k,icrm)-tbgmin)*a_bg))
            omp = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(i,j,k,icrm)-tprmin)*a_pr))
            omg = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(i,j,k,icrm)-tgrmin)*a_gr))

            if(qn(i,j,k,icrm).gt.0.) then

              qcc = qn(i,j,k,icrm) * omn
              qii = qn(i,j,k,icrm) * (1.-omn)

              if(qcc .gt. qcw0) then
                autor = alphaelq
              else
                autor = 0.
              endif

              if(qii .gt. qci0) then
                autos = betaelq*coefice(k,icrm)
              else
                autos = 0.
              endif

              accrr = 0.
              if(omp.gt.0.001) then
                qrr = qp(i,j,k) * omp
                accrr = accrrc(k,icrm) * qrr ** powr1
              end if
              accrcs = 0.
              accris = 0.
              if(omp.lt.0.999.and.omg.lt.0.999) then
                qss = qp(i,j,k) * (1.-omp)*(1.-omg)
                tmp = qss ** pows1
                accrcs = accrsc(k,icrm) * tmp
                accris = accrsi(k,icrm) * tmp
              end if
              accrcg = 0.
              accrig = 0.
              if(omp.lt.0.999.and.omg.gt.0.001) then
                qgg = qp(i,j,k) * (1.-omp)*omg
                tmp = qgg ** powg1
                accrcg = accrgc(k,icrm) * tmp
                accrig = accrgi(k,icrm) * tmp
              endif
              qcc = (qcc+dtn*autor*qcw0)/(1.+dtn*(accrr+accrcs+accrcg+autor))
              qii = (qii+dtn*autos*qci0)/(1.+dtn*(accris+accrig+autos))
              dq = dtn *(accrr*qcc + autor*(qcc-qcw0)+ &
              (accris+accrig)*qii + (accrcs+accrcg)*qcc + autos*(qii-qci0))
              dq = min(dq,qn(i,j,k,icrm))
              qp(i,j,k) = qp(i,j,k) + dq
              q(i,j,k) = q(i,j,k) - dq
              qn(i,j,k,icrm) = qn(i,j,k,icrm) - dq
              qpsrc(k,icrm) = qpsrc(k,icrm) + dq

            elseif(qp(i,j,k).gt.qp_threshold.and.qn(i,j,k,icrm).eq.0.) then

              qsatt = 0.
              if(omn.gt.0.001) qsatt = qsatt + omn*qsatw_crm(tabs(i,j,k,icrm),pres(k,icrm))
              if(omn.lt.0.999) qsatt = qsatt + (1.-omn)*qsati_crm(tabs(i,j,k,icrm),pres(k,icrm))
              dq = 0.
              if(omp.gt.0.001) then
                qrr = qp(i,j,k) * omp
                dq = dq + evapr1(k,icrm)*sqrt(qrr) + evapr2(k,icrm)*qrr**powr2
              end if
              if(omp.lt.0.999.and.omg.lt.0.999) then
                qss = qp(i,j,k) * (1.-omp)*(1.-omg)
                dq = dq + evaps1(k,icrm)*sqrt(qss) + evaps2(k,icrm)*qss**pows2
              end if
              if(omp.lt.0.999.and.omg.gt.0.001) then
                qgg = qp(i,j,k) * (1.-omp)*omg
                dq = dq + evapg1(k,icrm)*sqrt(qgg) + evapg2(k,icrm)*qgg**powg2
              end if
              dq = dq * dtn * (q(i,j,k) /qsatt-1.)
              dq = max(-0.5*qp(i,j,k),dq)
              qp(i,j,k) = qp(i,j,k) + dq
              q(i,j,k) = q(i,j,k) - dq
              qpevp(k,icrm) = qpevp(k,icrm) + dq

            else

              q(i,j,k) = q(i,j,k) + qp(i,j,k)
              qpevp(k,icrm) = qpevp(k,icrm) - qp(i,j,k)
              qp(i,j,k) = 0.

            endif

          endif

          dq = qp(i,j,k)
          qp(i,j,k)=max(real(0.,crm_rknd),qp(i,j,k))
          q(i,j,k) = q(i,j,k) + (dq-qp(i,j,k))

        end do
      enddo
    enddo



    !call t_stopf ('precip_proc')

  end subroutine precip_proc

end module precip_proc_mod
