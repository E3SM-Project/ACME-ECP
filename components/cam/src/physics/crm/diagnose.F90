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
    real(crm_rknd) tmp_lwp, tmp

    coef = 1./real(nx*ny,crm_rknd)

    !!$acc enter data create(k200,k500,k850) async(1)

    !!$acc parallel loop collapse(2) copy(u0,v0,t01,q01,t0,tabs0,q0,qn0,qp0,p0) async(1)
    do icrm = 1 , ncrms
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
      enddo
    enddo

    !!$acc parallel loop copy(k200,k500,k850) copyin(pres) async(1)
    do icrm = 1 , ncrms
      k200(icrm) = nzm
      k500(icrm) = nzm
      k850(icrm) = nzm
      do k=1,nzm
        kc=min(nzm,k+1)
        kb=max(1,k-1)
        if(pres(kc,icrm).le.200..and.pres(kb,icrm).gt.200.) k200(icrm)=k
        if(pres(kc,icrm).le.500..and.pres(kb,icrm).gt.500.) k500(icrm)=k
        if(pres(kc,icrm).le.850..and.pres(kb,icrm).gt.850.) k850(icrm)=k
      enddo
    enddo


    !!$acc parallel loop collapse(4) copyin(gamaz,qcl,qpl,p,u,qv,v,t,adz,qpi,qci,rho,dz) copy(p0,tabs,iw_xy,tabs0,cw_xy,qp0,q0,t0,u0,v0,pw_xy,qn0) async(1)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            coef1 = rho(k,icrm)*dz(icrm)*adz(k,icrm)*dtfactor
            tabs(i,j,k,icrm) = t(i,j,k,icrm)-gamaz(k,icrm)+ fac_cond * (qcl(i,j,k,icrm)+qpl(i,j,k,icrm)) + fac_sub *(qci(i,j,k,icrm) + qpi(i,j,k,icrm))
            !$acc atomic update
            u0(k,icrm)=u0(k,icrm)+u(i,j,k,icrm)
            !$acc atomic update
            v0(k,icrm)=v0(k,icrm)+v(i,j,k,icrm)
            !$acc atomic update
            p0(k,icrm)=p0(k,icrm)+p(i,j,k,icrm)
            !$acc atomic update
            t0(k,icrm)=t0(k,icrm)+t(i,j,k,icrm)
            !$acc atomic update
            tabs0(k,icrm)=tabs0(k,icrm)+tabs(i,j,k,icrm)
            tmp = qv(i,j,k,icrm)+qcl(i,j,k,icrm)+qci(i,j,k,icrm)
            !$acc atomic update
            q0(k,icrm)=q0(k,icrm)+tmp
            tmp = qcl(i,j,k,icrm) + qci(i,j,k,icrm)
            !$acc atomic update
            qn0(k,icrm) = qn0(k,icrm) + tmp
            tmp = qpl(i,j,k,icrm) + qpi(i,j,k,icrm)
            !$acc atomic update
            qp0(k,icrm) = qp0(k,icrm) + tmp
            tmp = qv(i,j,k,icrm)*coef1
            !$acc atomic update
            pw_xy(i,j,icrm) = pw_xy(i,j,icrm)+tmp
            tmp = qcl(i,j,k,icrm)*coef1
            !$acc atomic update
            cw_xy(i,j,icrm) = cw_xy(i,j,icrm)+tmp
            tmp = qci(i,j,k,icrm)*coef1
            !$acc atomic update
            iw_xy(i,j,icrm) = iw_xy(i,j,icrm)+tmp
          enddo
        enddo
      enddo
    enddo
    !!$acc parallel loop collapse(2) copy(qn0,q0,p0,t0,qp0,v0,u0,tabs0) async(1)
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

    !!$acc parallel loop collapse(3) copyin(k500,u,v,k200,w) copy(u200_xy,v200_xy,usfc_xy,vsfc_xy,w500_xy) async(1)
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

    !!$acc parallel loop collapse(2) copyin(qn0,q0) copy(qv0) async(1)
    do icrm = 1 , ncrms
      do k = 1 , nzm
        qv0(k,icrm) = q0(k,icrm) - qn0(k,icrm)
      enddo
    enddo

    !=====================================================
    ! UW ADDITIONS
    ! FIND VERTICAL INDICES OF 850MB, COMPUTE SWVP
    !!$acc parallel loop collapse(4) copyin(rho,adz,dz) copy(pres,tabs,swvp_xy) async(1)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            coef1 = rho(k,icrm)*dz(icrm)*adz(k,icrm)*dtfactor
            ! Saturated water vapor path with respect to water. Can be used
            ! with water vapor path (= pw) to compute column-average
            ! relative humidity.
            tmp = qsatw_crm(tabs(i,j,k,icrm),pres(k,icrm))*coef1
            !$acc atomic update
            swvp_xy(i,j,icrm) = swvp_xy(i,j,icrm)+tmp
          enddo
        enddo
      enddo ! k
    enddo

    ! ACCUMULATE AVERAGES OF TWO-DIMENSIONAL STATISTICS
    !!$acc parallel loop collapse(3) copyin(k850,p,pres,u,v) copy(psfc_xy,u850_xy,v850_xy) async(1)
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
    !!$acc parallel loop collapse(3) copyin(sstxy) copy(cloudtopheight,cloudtoptemp,echotopheight) async(1)
    do icrm = 1 , ncrms
      do j = 1 , ny
        do i = 1 , nx
          cloudtopheight(i,j,icrm) = 0.
          cloudtoptemp  (i,j,icrm) = sstxy(i,j,icrm)
          echotopheight (i,j,icrm) = 0.
        enddo
      enddo
    enddo
    !!$acc parallel loop collapse(3) copyin(qcl,qpi,qci,qpl,adz,rho,z,tabs,dz) copy(echotopheight,cloudtoptemp,cloudtopheight,cld_xy) async(1)
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

    !!$acc exit data delete(k200,k500,k850) async(1)

  end subroutine diagnose

end module diagnose_mod
