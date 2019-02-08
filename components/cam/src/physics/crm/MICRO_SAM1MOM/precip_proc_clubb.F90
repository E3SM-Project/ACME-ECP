#define CLDFRAC
#ifdef CLDFRAC
module precip_proc_clubb_mod
  implicit none

contains

  subroutine precip_proc_clubb()

#ifdef CLUBB_CRM
    use vars
    use microphysics
    use micro_params
    use params
    use vars, only: CF3D

    implicit none
    integer i,j,k
    real(crm_rknd) autor, autos, accrr, accris, accrcs, accrig, accrcg
    real(crm_rknd) dq, omn, omp, omg, qsatt
    real(crm_rknd) pows1, pows2, powg1, powg2, powr1, powr2, tmp
    real(crm_rknd) qii, qcc, qrr, qss, qgg

    real(crm_rknd) cld3d(nx, ny, nzm), cldmax(nx, ny, nzm)
    real(crm_rknd) cld3d_temp(nx, ny, nzm)
    real(crm_rknd) cloud_frac_thresh
    real(crm_rknd) qclr
    real(crm_rknd) dqpsrc, dqpevp

    powr1 = (3 + b_rain) / 4.
    powr2 = (5 + b_rain) / 8.
    pows1 = (3 + b_snow) / 4.
    pows2 = (5 + b_snow) / 8.
    powg1 = (3 + b_grau) / 4.
    powg2 = (5 + b_grau) / 8.

    !call t_startf ('precip_proc_clubb')

    ! Get cloud fraction of non-precipitating condensate
    ! and precipitating condensate
    cloud_frac_thresh = 0.005
    do j=1, ny
      do i=1, nx
        do k=nzm, 1, -1
          cld3d(i, j, k) = CF3D(i,j,k,icrm)
          cld3d_temp(i, j, k) = min(0.999, max(CF3D(i,j,k,icrm), cloud_frac_thresh))
        end do
        cldmax(i,j,nzm)=cld3d_temp(i,j,nzm)

        do k=nzm-1, 1, -1
          ! if precipitating condensate is smaller than threshold, set cldmax
          ! to cloud fraction at current level
          if(qp(icrm,i, j, k+1).ge.qp_threshold) then
            cldmax(i,j,k) = max(cldmax(i,j,k+1), cld3d_temp(i,j,k))
          else
            cldmax(i,j,k) = cld3d_temp(i,j,k)
          end if

          !    if(cld3d(i,j,k).le.cloud_frac_thresh .and. qp(icrm,i,j,k).gt.qp_threshold) then
          !       if(cldmax(i,j,k).lt.0.1) then
          !         cldmax(i,j,k) = 0.50
          !       end if
          !    end if
        end do
        ! test: assume precipitating hydrometer fill the whole grid box
        !  cldmax(i,j,:) = 0.999

      end do
    end do


    do k=1,nzm
      qpsrc(k,icrm)=0.
      qpevp(k,icrm)=0.
      do j=1,ny
        do i=1,nx
          dqpsrc = 0.0
          dqpevp = 0.0

          !-------     Autoconversion/accretion

          if(qn(i,j,k,icrm)+qp(icrm,i,j,k).gt.0.) then


            omn = max(0.,min(1.,(tabs(i,j,k,icrm)-tbgmin)*a_bg))
            omp = max(0.,min(1.,(tabs(i,j,k,icrm)-tprmin)*a_pr))
            omg = max(0.,min(1.,(tabs(i,j,k,icrm)-tgrmin)*a_gr))

            !	 if(qn(i,j,k,icrm).gt.0.) then
            if(cld3d(i,j,k).gt.0.) then  ! the generation of precipitating condensate

              qcc = qn(i,j,k,icrm) * omn /cld3d_temp(i,j,k)
              qii = qn(i,j,k,icrm) * (1.-omn)/cld3d_temp(i,j,k)

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
                qrr = qp(icrm,i,j,k) * omp / cldmax(i,j,k)
                accrr = accrrc(k,icrm) * qrr ** powr1
              end if
              accrcs = 0.
              accris = 0.
              if(omp.lt.0.999.and.omg.lt.0.999) then
                qss = qp(icrm,i,j,k) * (1.-omp)*(1.-omg) / cldmax(i,j,k)
                tmp = qss ** pows1
                accrcs = accrsc(k,icrm) * tmp
                accris = accrsi(k,icrm) * tmp
              end if
              accrcg = 0.
              accrig = 0.
              if(omp.lt.0.999.and.omg.gt.0.001) then
                qgg = qp(icrm,i,j,k) * (1.-omp)*omg / cldmax(i,j,k)
                tmp = qgg ** powg1
                accrcg = accrgc(k,icrm) * tmp
                accrig = accrgi(k,icrm) * tmp
              endif
              qcc = (qcc+dtn*autor*qcw0)/(1.+dtn*(accrr+accrcs+accrcg+autor))
              qii = (qii+dtn*autos*qci0)/(1.+dtn*(accris+accrig+autos))
              dq = dtn *(accrr*qcc + autor*(qcc-qcw0)+ &
              (accris+accrig)*qii + (accrcs+accrcg)*qcc + autos*(qii-qci0))

              dq = dq * cld3d(i,j,k)  ! convert fro the in-cloud value to grid-mean value

              dq = min(dq,qn(i,j,k,icrm))
              !           qp(icrm,i,j,k) = qp(icrm,i,j,k) + dq
              !           q(icrm,i,j,k) = q(icrm,i,j,k) - dq
              !           qn(i,j,k,icrm) = qn(i,j,k,icrm) - dq
              dqpsrc = dq
              qpsrc(k,icrm) = qpsrc(k,icrm) + dq

            end if

            !elseif(qp(icrm,i,j,k).gt.qp_threshold.and.qn(i,j,k,icrm).eq.0.) then
            ! Evaporation is only allowed when cldmax exceeds cld3d_temp
            !         if(qp(icrm,i,j,k).gt.qp_threshold.and.cldmax(i,j,k).gt.cld3d_temp(i,j,k)) then
            if(qp(icrm,i,j,k).gt.qp_threshold.and.qn(i,j,k,icrm).eq.0.) then

              qsatt = 0.
              if(omn.gt.0.001) qsatt = qsatt + omn*qsatw_crm(tabs(i,j,k,icrm),pres(k,icrm))
              if(omn.lt.0.999) qsatt = qsatt + (1.-omn)*qsati_crm(tabs(i,j,k,icrm),pres(k,icrm))
              dq = 0.
              if(omp.gt.0.001) then
                qrr = qp(icrm,i,j,k) * omp /cldmax(i,j,k)
                dq = dq + evapr1(k,icrm)*sqrt(qrr) + evapr2(k,icrm)*qrr**powr2
              end if
              if(omp.lt.0.999.and.omg.lt.0.999) then
                qss = qp(icrm,i,j,k) * (1.-omp)*(1.-omg) / cldmax(i,j,k)
                dq = dq + evaps1(k,icrm)*sqrt(qss) + evaps2(k,icrm)*qss**pows2
              end if
              if(omp.lt.0.999.and.omg.gt.0.001) then
                qgg = qp(icrm,i,j,k) * (1.-omp)*omg /cldmax(i,j,k)
                dq = dq + evapg1(k,icrm)*sqrt(qgg) + evapg2(k,icrm)*qgg**powg2
              end if

              !           dq = dq * dtn * (q(icrm,i,j,k) /qsatt-1.)
              qclr = max(0., (q(icrm,i,j,k)-qn(i,j,k,icrm)-qsatt * cld3d(i,j,k)))/max(0.001, (1-cld3d(i,j,k)))
              qclr = min(qclr, qsatt)
              dq = dq * dtn * (qclr/qsatt-1.)
              dq = dq * (cldmax(i,j,k) - cld3d_temp(i,j,k))  ! convert this to the grid-mean value

              dq = max(-0.5*qp(icrm,i,j,k),dq)
              !           qp(icrm,i,j,k) = qp(icrm,i,j,k) + dq
              !           q(icrm,i,j,k) = q(icrm,i,j,k) - dq
              dqpevp = dq
              qpevp(k,icrm) = qpevp(k,icrm) + dq

            end if

            if(qp(icrm,i,j,k).le.qp_threshold .and. cld3d(i,j,k).le.0) then
              !           q(icrm,i,j,k) = q(icrm,i,j,k) + qp(icrm,i,j,k)
              dqpevp = dqpevp - qp(icrm,i,j,k)
              qpevp(k,icrm) = qpevp(k,icrm) - qp(icrm,i,j,k)
              !           qp(icrm,i,j,k) = 0.
            endif

          endif

          qp(icrm,i,j,k) = qp(icrm,i,j,k) + dqpsrc + dqpevp
          q(icrm,i,j,k) = q(icrm,i,j,k) - dqpsrc - dqpevp
          qn(i,j,k,icrm) = qn(i,j,k,icrm) - dqpsrc

          dq = qp(icrm,i,j,k)
          qp(icrm,i,j,k)=max(0.,qp(icrm,i,j,k))
          q(icrm,i,j,k) = q(icrm,i,j,k) + (dq-qp(icrm,i,j,k))

        end do
      enddo
    enddo

    !call t_stopf ('precip_proc_clubb')

#endif /*CLUBB_CRM*/
  end subroutine precip_proc_clubb

end module precip_proc_clubb_mod
#endif
