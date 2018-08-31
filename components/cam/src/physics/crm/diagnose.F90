module diagnose_mod
  use sat_mod
  use task_util_mod
  implicit none

contains

  subroutine diagnose(ncrms,icrm)

    ! Diagnose some useful stuff

    use vars
    use params
    use sgs, only: sgs_diagnose
    implicit none
    integer, intent(in) :: ncrms,icrm
    integer i,j,k,kb,kc,k200,k500,k850
    real(8) coef, coef1, buffer(nzm,9), buffer1(nzm,8)
    real(crm_rknd) omn, omp, tmp_lwp

    coef = 1./real(nx*ny,crm_rknd)


    k200 = nzm

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
      if(pres(kc).le.200..and.pres(kb).gt.200.) k200=k
      coef1 = rho(k,icrm)*dz*adz(k)*dtfactor
      do j=1,ny
        do i=1,nx
          tabs(i,j,k,icrm) = t(i,j,k)-gamaz(k,icrm)+ fac_cond * (qcl(i,j,k,icrm)+qpl(i,j,k,icrm)) +&
          fac_sub *(qci(i,j,k,icrm) + qpi(i,j,k,icrm))
          u0(k,icrm)=u0(k,icrm)+u(i,j,k)
          v0(k,icrm)=v0(k,icrm)+v(i,j,k)
          p0(k,icrm)=p0(k,icrm)+p(i,j,k,icrm)
          t0(k,icrm)=t0(k,icrm)+t(i,j,k)
          tabs0(k,icrm)=tabs0(k,icrm)+tabs(i,j,k,icrm)
          q0(k,icrm)=q0(k,icrm)+qv(i,j,k,icrm)+qcl(i,j,k,icrm)+qci(i,j,k,icrm)
          qn0(k,icrm) = qn0(k,icrm) + qcl(i,j,k,icrm) + qci(i,j,k,icrm)
          qp0(k,icrm) = qp0(k,icrm) + qpl(i,j,k,icrm) + qpi(i,j,k,icrm)

          pw_xy(i,j,icrm) = pw_xy(i,j,icrm)+qv(i,j,k,icrm)*coef1
          cw_xy(i,j,icrm) = cw_xy(i,j,icrm)+qcl(i,j,k,icrm)*coef1
          iw_xy(i,j,icrm) = iw_xy(i,j,icrm)+qci(i,j,k,icrm)*coef1

        end do
      end do
      u0(k,icrm)=u0(k,icrm)*coef
      v0(k,icrm)=v0(k,icrm)*coef
      t0(k,icrm)=t0(k,icrm)*coef
      tabs0(k,icrm)=tabs0(k,icrm)*coef
      q0(k,icrm)=q0(k,icrm)*coef
      qn0(k,icrm)=qn0(k,icrm)*coef
      qp0(k,icrm)=qp0(k,icrm)*coef
      p0(k,icrm)=p0(k,icrm)*coef

    end do ! k

    k500 = nzm
    do k = 1,nzm
      kc=min(nzm,k+1)
      if((pres(kc).le.500.).and.(pres(k).gt.500.)) then
        if ((500.-pres(kc)).lt.(pres(k)-500.))then
          k500=kc
        else
          k500=k
        end if
      end if
    end do


    do j=1,ny
      do i=1,nx
        usfc_xy(i,j,icrm) = usfc_xy(i,j,icrm) + u(i,j,1)*dtfactor
        vsfc_xy(i,j,icrm) = vsfc_xy(i,j,icrm) + v(i,j,1)*dtfactor
        u200_xy(i,j,icrm) = u200_xy(i,j,icrm) + u(i,j,k200)*dtfactor
        v200_xy(i,j,icrm) = v200_xy(i,j,icrm) + v(i,j,k200)*dtfactor
        w500_xy(i,j,icrm) = w500_xy(i,j,icrm) + w(i,j,k500)*dtfactor
      end do
    end do

    if(dompi) then

      coef1 = 1./real(nsubdomains,crm_rknd)
      do k=1,nzm
        buffer(k,1) = u0(k,icrm)
        buffer(k,2) = v0(k,icrm)
        buffer(k,3) = t0(k,icrm)
        buffer(k,4) = q0(k,icrm)
        buffer(k,5) = p0(k,icrm)
        buffer(k,6) = tabs0(k,icrm)
        buffer(k,7) = qn0(k,icrm)
        buffer(k,8) = qp0(k,icrm)
      end do
      call task_sum_real8(buffer,buffer1,nzm*8)
      do k=1,nzm
        u0(k,icrm)=buffer1(k,1)*coef1
        v0(k,icrm)=buffer1(k,2)*coef1
        t0(k,icrm)=buffer1(k,3)*coef1
        q0(k,icrm)=buffer1(k,4)*coef1
        p0(k,icrm)=buffer1(k,5)*coef1
        tabs0(k,icrm)=buffer1(k,6)*coef1
        qn0(k,icrm)=buffer1(k,7)*coef1
        qp0(k,icrm)=buffer1(k,8)*coef1
      end do

    end if ! dompi

    qv0(:,icrm) = q0(:,icrm) - qn0(:,icrm)

    !=====================================================
    ! UW ADDITIONS

    ! FIND VERTICAL INDICES OF 850MB, COMPUTE SWVP
    k850 = 1
    do k = 1,nzm
      if(pres(k).le.850.) then
        k850 = k
        EXIT
      end if
    end do

    do k=1,nzm
      coef1 = rho(k,icrm)*dz*adz(k)*dtfactor
      do j=1,ny
        do i=1,nx

          ! Saturated water vapor path with respect to water. Can be used
          ! with water vapor path (= pw) to compute column-average
          ! relative humidity.
          swvp_xy(i,j,icrm) = swvp_xy(i,j,icrm)+qsatw_crm(tabs(i,j,k,icrm),pres(k))*coef1
        end do
      end do
    end do ! k

    ! ACCUMULATE AVERAGES OF TWO-DIMENSIONAL STATISTICS
    do j=1,ny
      do i=1,nx
        psfc_xy(i,j,icrm) = psfc_xy(i,j,icrm) + (100.*pres(1) + p(i,j,1,icrm))*dtfactor

        ! 850 mbar horizontal winds
        u850_xy(i,j,icrm) = u850_xy(i,j,icrm) + u(i,j,k850)*dtfactor
        v850_xy(i,j,icrm) = v850_xy(i,j,icrm) + v(i,j,k850)*dtfactor

      end do
    end do

    ! COMPUTE CLOUD/ECHO HEIGHTS AS WELL AS CLOUD TOP TEMPERATURE
    ! WHERE CLOUD TOP IS DEFINED AS THE HIGHEST MODEL LEVEL WITH A
    ! CONDENSATE PATH OF 0.01 kg/m2 ABOVE.  ECHO TOP IS THE HIGHEST LEVEL
    ! WHERE THE PRECIPITATE MIXING RATIO > 0.001 G/KG.

    ! initially, zero out heights and set cloudtoptemp to SST
    cloudtopheight(:,:,icrm) = 0.
    cloudtoptemp(:,:,icrm) = sstxy(1:nx,1:ny,icrm)
    echotopheight(:,:,icrm) = 0.
    do j = 1,ny
      do i = 1,nx
        ! FIND CLOUD TOP HEIGHT
        tmp_lwp = 0.
        do k = nzm,1,-1
          tmp_lwp = tmp_lwp + (qcl(i,j,k,icrm)+qci(i,j,k,icrm))*rho(k,icrm)*dz*adz(k)
          if (tmp_lwp.gt.0.01) then
            cloudtopheight(i,j,icrm) = z(k)
            cloudtoptemp(i,j,icrm) = tabs(i,j,k,icrm)
            cld_xy(i,j,icrm) = cld_xy(i,j,icrm) + dtfactor
            EXIT
          end if
        end do
        ! FIND ECHO TOP HEIGHT
        do k = nzm,1,-1
          if (qpl(i,j,k,icrm)+qpi(i,j,k,icrm).gt.1.e-6) then
            echotopheight(i,j,icrm) = z(k)
            EXIT
          end if
        end do
      end do
    end do

    ! END UW ADDITIONS
    !=====================================================

    !-----------------
    ! compute some sgs diagnostics:

    call sgs_diagnose()

    !-----------------

    ! recompute pressure levels, except at restart (saved levels are used).
    !if(dtfactor.ge.0.) call pressz()   ! recompute pressure levels

  end subroutine diagnose

end module diagnose_mod
