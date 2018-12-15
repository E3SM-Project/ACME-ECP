module cloud_mod
  implicit none

contains

  subroutine cloud(ncrms,q,qp,qn)

    !  Condensation of cloud water/cloud ice.

    use vars
    use micro_params
    use params
    use sat_mod

    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) :: q (dimx1_s:dimx2_s, dimy1_s:dimy2_s, 1:nzm, 1:ncrms)
    real(crm_rknd) :: qp(dimx1_s:dimx2_s, dimy1_s:dimy2_s, 1:nzm, 1:ncrms)
    real(crm_rknd) qn(nx,ny,nzm,ncrms)  ! cloud condensate (liquid + ice)

    integer i,j,k, kb, kc,icrm
    real(crm_rknd) dtabs, tabs1, an, bn, ap, bp, om, ag, omp
    real(crm_rknd) fac1,fac2
    real(crm_rknd) fff,dfff,qsatt,dqsat
    real(crm_rknd) lstarn,dlstarn,lstarp,dlstarp
    integer niter

    an = 1./(tbgmax-tbgmin)
    bn = tbgmin * an
    ap = 1./(tprmax-tprmin)
    bp = tprmin * ap
    fac1 = fac_cond+(1+bp)*fac_fus
    fac2 = fac_fus*ap
    ag = 1./(tgrmax-tgrmin)

    !$acc enter data copyin(t,gamaz,q,qp,tabs,qn,pres,qsatt,dtabs,dqsat) async(1)

    !$acc parallel loop collapse(4) default(present) async(1)
    do icrm = 1 , ncrms
      do k = 1, nzm
        do j = 1, ny
          do i = 1, nx

            q(i,j,k,icrm)=max(real(0.,crm_rknd),q(i,j,k,icrm))


            ! Initail guess for temperature assuming no cloud water/ice:


            tabs(i,j,k,icrm) = t(i,j,k,icrm)-gamaz(k,icrm)
            tabs1=(tabs(i,j,k,icrm)+fac1*qp(i,j,k,icrm))/(1.+fac2*qp(i,j,k,icrm))

            ! Warm cloud:

            if(tabs1.ge.tbgmax) then

              tabs1=tabs(i,j,k,icrm)+fac_cond*qp(i,j,k,icrm)
              qsatt = qsatw_crm(tabs1,pres(k,icrm))

              ! Ice cloud:

            elseif(tabs1.le.tbgmin) then

              tabs1=tabs(i,j,k,icrm)+fac_sub*qp(i,j,k,icrm)
              qsatt = qsati_crm(tabs1,pres(k,icrm))

              ! Mixed-phase cloud:

            else

              om = an*tabs1-bn
              qsatt = om*qsatw_crm(tabs1,pres(k,icrm))+(1.-om)*qsati_crm(tabs1,pres(k,icrm))

            endif


            !  Test if condensation is possible:


            if(q(i,j,k,icrm).gt.qsatt) then

              niter=0
              dtabs = 100.
              do while(abs(dtabs).gt.0.01.and.niter.lt.10)
                if(tabs1.ge.tbgmax) then
                  om=1.
                  lstarn=fac_cond
                  dlstarn=0.
                  qsatt=qsatw_crm(tabs1,pres(k,icrm))
                  dqsat=dtqsatw_crm(tabs1,pres(k,icrm))
                else if(tabs1.le.tbgmin) then
                  om=0.
                  lstarn=fac_sub
                  dlstarn=0.
                  qsatt=qsati_crm(tabs1,pres(k,icrm))
                  dqsat=dtqsati_crm(tabs1,pres(k,icrm))
                else
                  om=an*tabs1-bn
                  lstarn=fac_cond+(1.-om)*fac_fus
                  dlstarn=an*fac_fus
                  qsatt=om*qsatw_crm(tabs1,pres(k,icrm))+(1.-om)*qsati_crm(tabs1,pres(k,icrm))
                  dqsat=om*dtqsatw_crm(tabs1,pres(k,icrm))+(1.-om)*dtqsati_crm(tabs1,pres(k,icrm))
                endif
                if(tabs1.ge.tprmax) then
                  omp=1.
                  lstarp=fac_cond
                  dlstarp=0.
                else if(tabs1.le.tprmin) then
                  omp=0.
                  lstarp=fac_sub
                  dlstarp=0.
                else
                  omp=ap*tabs1-bp
                  lstarp=fac_cond+(1.-omp)*fac_fus
                  dlstarp=ap*fac_fus
                endif
                fff = tabs(i,j,k,icrm)-tabs1+lstarn*(q(i,j,k,icrm)-qsatt)+lstarp*qp(i,j,k,icrm)
                dfff=dlstarn*(q(i,j,k,icrm)-qsatt)+dlstarp*qp(i,j,k,icrm)-lstarn*dqsat-1.
                dtabs=-fff/dfff
                niter=niter+1
                tabs1=tabs1+dtabs
              end do

              qsatt = qsatt + dqsat * dtabs
              qn(i,j,k,icrm) = max(real(0.,crm_rknd),q(i,j,k,icrm)-qsatt)

            else

              qn(i,j,k,icrm) = 0.

            endif

            tabs(i,j,k,icrm) = tabs1
            qp(i,j,k,icrm) = max(real(0.,crm_rknd),qp(i,j,k,icrm)) ! just in case

          end do
        end do
      end do
    end do

    !$acc exit data copyout(t,gamaz,q,qp,tabs,qn,pres,qsatt,dtabs,dqsat) async(1)

  end subroutine cloud

end module cloud_mod
