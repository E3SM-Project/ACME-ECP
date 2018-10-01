module diagnose_mod
  use sat_mod
  use task_util_mod
  implicit none

contains

  subroutine diagnose(ncrms)
    ! Diagnose some useful stuff
    use vars
    use params
    use sgs, only: sgs_diagnose
    implicit none
    integer, intent(in) :: ncrms
    integer i,j,k,kb,kc,k200(ncrms),k500(ncrms),k850(ncrms),icrm
    real(8) coef, coef1
    real(crm_rknd) tmp_lwp

    coef = 1./real(nx*ny,crm_rknd)

    do icrm = 1 , ncrms
      k200(icrm) = nzm
      k500(icrm) = nzm
      k850(icrm) = nzm
      !$acc loop seq
      do k=1,nzm
        u0(k,icrm)=0.
        v0(k,icrm)=0.
        t01(k,icrm) = tabs0(k,icrm)
        q01(k,icrm) = q0(k,icrm)
        t0(k,icrm)=0.
        tabs0(k,icrm)=0.
        q0(k,icrm)=0.
        qn0(k,icrm)=0.
        qp0(k,icrm)=0.
        p0(k,icrm)=0.
        kc=min(nzm,k+1)
        kb=max(1,k-1)
        if(pres(kc,icrm).le.200..and.pres(kb,icrm).gt.200.) k200(icrm)=k
        if(pres(kc,icrm).le.500..and.pres(kb,icrm).gt.500.) k500(icrm)=k
        if(pres(kc,icrm).le.850..and.pres(kb,icrm).gt.850.) k850(icrm)=k
      enddo
    enddo
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            coef1 = rho(k,icrm)*dz(icrm)*adz(k,icrm)*dtfactor
            tabs(i,j,k,icrm) = t(i,j,k,icrm)-gamaz(k,icrm)+ fac_cond * (qcl(i,j,k,icrm)+qpl(i,j,k,icrm)) +&
            fac_sub *(qci(i,j,k,icrm) + qpi(i,j,k,icrm))
            u0(k,icrm)=u0(k,icrm)+u(i,j,k,icrm)
            v0(k,icrm)=v0(k,icrm)+v(i,j,k,icrm)
            p0(k,icrm)=p0(k,icrm)+p(i,j,k,icrm)
            t0(k,icrm)=t0(k,icrm)+t(i,j,k,icrm)
            tabs0(k,icrm)=tabs0(k,icrm)+tabs(i,j,k,icrm)
            q0(k,icrm)=q0(k,icrm)+qv(i,j,k,icrm)+qcl(i,j,k,icrm)+qci(i,j,k,icrm)
            qn0(k,icrm) = qn0(k,icrm) + qcl(i,j,k,icrm) + qci(i,j,k,icrm)
            qp0(k,icrm) = qp0(k,icrm) + qpl(i,j,k,icrm) + qpi(i,j,k,icrm)
            pw_xy(i,j,icrm) = pw_xy(i,j,icrm)+qv(i,j,k,icrm)*coef1
            cw_xy(i,j,icrm) = cw_xy(i,j,icrm)+qcl(i,j,k,icrm)*coef1
            iw_xy(i,j,icrm) = iw_xy(i,j,icrm)+qci(i,j,k,icrm)*coef1
          enddo
        enddo
      enddo
    enddo
    do icrm = 1 , ncrms
      do k=1,nzm
        u0(k,icrm)=u0(k,icrm)*coef
        v0(k,icrm)=v0(k,icrm)*coef
        t0(k,icrm)=t0(k,icrm)*coef
        tabs0(k,icrm)=tabs0(k,icrm)*coef
        q0(k,icrm)=q0(k,icrm)*coef
        qn0(k,icrm)=qn0(k,icrm)*coef
        qp0(k,icrm)=qp0(k,icrm)*coef
        p0(k,icrm)=p0(k,icrm)*coef
      enddo ! k
    enddo

    do icrm = 1 , ncrms
      do j=1,ny
        do i=1,nx
          usfc_xy(i,j,icrm) = usfc_xy(i,j,icrm) + u(i,j,1,icrm)*dtfactor
          vsfc_xy(i,j,icrm) = vsfc_xy(i,j,icrm) + v(i,j,1,icrm)*dtfactor
          u200_xy(i,j,icrm) = u200_xy(i,j,icrm) + u(i,j,k200(icrm),icrm)*dtfactor
          v200_xy(i,j,icrm) = v200_xy(i,j,icrm) + v(i,j,k200(icrm),icrm)*dtfactor
          w500_xy(i,j,icrm) = w500_xy(i,j,icrm) + w(i,j,k500(icrm),icrm)*dtfactor
        enddo
      enddo
    enddo

    do icrm = 1 , ncrms
      qv0(:,icrm) = q0(:,icrm) - qn0(:,icrm)
    enddo

    !=====================================================
    ! UW ADDITIONS
    ! FIND VERTICAL INDICES OF 850MB, COMPUTE SWVP
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            coef1 = rho(k,icrm)*dz(icrm)*adz(k,icrm)*dtfactor
            ! Saturated water vapor path with respect to water. Can be used
            ! with water vapor path (= pw) to compute column-average
            ! relative humidity.
            swvp_xy(i,j,icrm) = swvp_xy(i,j,icrm)+qsatw_crm(tabs(i,j,k,icrm),pres(k,icrm))*coef1
          enddo
        enddo
      enddo ! k
    enddo

    ! ACCUMULATE AVERAGES OF TWO-DIMENSIONAL STATISTICS
    do icrm = 1 , ncrms
      do j=1,ny
        do i=1,nx
          psfc_xy(i,j,icrm) = psfc_xy(i,j,icrm) + (100.*pres(1,icrm) + p(i,j,1,icrm))*dtfactor
          ! 850 mbar horizontal winds
          u850_xy(i,j,icrm) = u850_xy(i,j,icrm) + u(i,j,k850(icrm),icrm)*dtfactor
          v850_xy(i,j,icrm) = v850_xy(i,j,icrm) + v(i,j,k850(icrm),icrm)*dtfactor
        enddo
      enddo
    enddo

    ! COMPUTE CLOUD/ECHO HEIGHTS AS WELL AS CLOUD TOP TEMPERATURE
    ! WHERE CLOUD TOP IS DEFINED AS THE HIGHEST MODEL LEVEL WITH A
    ! CONDENSATE PATH OF 0.01 kg/m2 ABOVE.  ECHO TOP IS THE HIGHEST LEVEL
    ! WHERE THE PRECIPITATE MIXING RATIO > 0.001 G/KG.
    ! initially, zero out heights and set cloudtoptemp to SST
    do icrm = 1 , ncrms
      cloudtopheight(:,:,icrm) = 0.
      cloudtoptemp(:,:,icrm) = sstxy(1:nx,1:ny,icrm)
      echotopheight(:,:,icrm) = 0.
    enddo
    do icrm = 1 , ncrms
      do j = 1,ny
        do i = 1,nx
          ! FIND CLOUD TOP HEIGHT
          tmp_lwp = 0.
          do k = nzm,1,-1
            tmp_lwp = tmp_lwp + (qcl(i,j,k,icrm)+qci(i,j,k,icrm))*rho(k,icrm)*dz(icrm)*adz(k,icrm)
            if (tmp_lwp.gt.0.01) then
              cloudtopheight(i,j,icrm) = z(k,icrm)
              cloudtoptemp(i,j,icrm) = tabs(i,j,k,icrm)
              cld_xy(i,j,icrm) = cld_xy(i,j,icrm) + dtfactor
              EXIT
            endif
          enddo
          ! FIND ECHO TOP HEIGHT
          do k = nzm,1,-1
            if (qpl(i,j,k,icrm)+qpi(i,j,k,icrm).gt.1.e-6) then
              echotopheight(i,j,icrm) = z(k,icrm)
              EXIT
            endif
          enddo
        enddo
      enddo
    enddo
    ! END UW ADDITIONS
    !=====================================================

    ! compute some sgs diagnostics:
    !call sgs_diagnose()

  end subroutine diagnose

end module diagnose_mod
