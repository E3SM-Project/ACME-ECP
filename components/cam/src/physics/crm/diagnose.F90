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

    integer :: i,j,k,kb,kc,icrm
    real(8) coef, coef1(ncrms), tmp
    real(crm_rknd) omn, omp, tmp_lwp
    integer, allocatable :: k200(:),k500(:),k850(:)

    allocate(k200(ncrms))
    allocate(k500(ncrms))
    allocate(k850(ncrms))

    coef = 1./real(nx*ny,crm_rknd)

    k200(:) = nzm
    !$acc parallel loop gang vector collapse(2)
    do k=1,nzm
      do icrm = 1 , ncrms
        u0(icrm,k)=0.
        v0(icrm,k)=0.
        t01(icrm,k) = tabs0(icrm,k)
        q01(icrm,k) = q0(icrm,k)
        t0(icrm,k)=0.
        tabs0(icrm,k)=0.
        q0(icrm,k)=0.
        qn0(icrm,k)=0.
        qp0(icrm,k)=0.
        p0(icrm,k)=0.
        kc=min(nzm,k+1)
        kb=max(1,k-1)
        if(pres(icrm,kc).le.200..and.pres(icrm,kb).gt.200.) k200(icrm)=k
        coef1(icrm) = rho(icrm,k)*dz(icrm)*adz(icrm,k)*dtfactor
      enddo
    enddo
    !$acc parallel loop gang vector collapse(4)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            tabs(icrm,i,j,k) = t(icrm,i,j,k)-gamaz(icrm,k)+ fac_cond * (qcl(icrm,i,j,k)+qpl(icrm,i,j,k)) +&
            fac_sub *(qci(icrm,i,j,k) + qpi(icrm,i,j,k))
            !$acc atomic update
            u0(icrm,k)=u0(icrm,k)+u(icrm,i,j,k)
            !$acc atomic update
            v0(icrm,k)=v0(icrm,k)+v(icrm,i,j,k)
            !$acc atomic update
            p0(icrm,k)=p0(icrm,k)+p(icrm,i,j,k)
            !$acc atomic update
            t0(icrm,k)=t0(icrm,k)+t(icrm,i,j,k)
            !$acc atomic update
            tabs0(icrm,k)=tabs0(icrm,k)+tabs(icrm,i,j,k)
            tmp = qv(icrm,i,j,k)+qcl(icrm,i,j,k)+qci(icrm,i,j,k)
            !$acc atomic update
            q0(icrm,k)=q0(icrm,k)+tmp
            tmp = qcl(icrm,i,j,k) + qci(icrm,i,j,k)
            !$acc atomic update
            qn0(icrm,k) = qn0(icrm,k) + tmp
            tmp = qpl(icrm,i,j,k) + qpi(icrm,i,j,k)
            !$acc atomic update
            qp0(icrm,k) = qp0(icrm,k) + tmp

            tmp = qv (icrm,i,j,k)*coef1(icrm)
            !$acc atomic update
            pw_xy(icrm,i,j) = pw_xy(icrm,i,j)+tmp
            tmp = qcl(icrm,i,j,k)*coef1(icrm)
            !$acc atomic update
            cw_xy(icrm,i,j) = cw_xy(icrm,i,j)+tmp
            tmp = qci(icrm,i,j,k)*coef1(icrm)
            !$acc atomic update
            iw_xy(icrm,i,j) = iw_xy(icrm,i,j)+tmp
          end do
        end do
      end do
    end do
    !$acc parallel loop gang vector collapse(2)
    do k=1,nzm
      do icrm = 1 , ncrms
        u0(icrm,k)=u0(icrm,k)*coef
        v0(icrm,k)=v0(icrm,k)*coef
        t0(icrm,k)=t0(icrm,k)*coef
        tabs0(icrm,k)=tabs0(icrm,k)*coef
        q0(icrm,k)=q0(icrm,k)*coef
        qn0(icrm,k)=qn0(icrm,k)*coef
        qp0(icrm,k)=qp0(icrm,k)*coef
        p0(icrm,k)=p0(icrm,k)*coef
      end do
    end do ! k

    !$acc parallel loop gang vector
    do icrm = 1 , ncrms
      k500(icrm) = nzm
      !$acc loop seq
      do k = 1,nzm
        kc=min(nzm,k+1)
        if((pres(icrm,kc).le.500.).and.(pres(icrm,k).gt.500.)) then
          if ((500.-pres(icrm,kc)).lt.(pres(icrm,k)-500.))then
            k500(icrm)=kc
          else
            k500(icrm)=k
          end if
        end if
      end do
    end do

    !$acc parallel loop gang vector collapse(3) copyin(w)
    do j=1,ny
      do i=1,nx
        do icrm = 1 , ncrms
          usfc_xy(icrm,i,j) = usfc_xy(icrm,i,j) + u(icrm,i,j,1)*dtfactor
          vsfc_xy(icrm,i,j) = vsfc_xy(icrm,i,j) + v(icrm,i,j,1)*dtfactor
          u200_xy(icrm,i,j) = u200_xy(icrm,i,j) + u(icrm,i,j,k200(icrm))*dtfactor
          v200_xy(icrm,i,j) = v200_xy(icrm,i,j) + v(icrm,i,j,k200(icrm))*dtfactor
          w500_xy(icrm,i,j) = w500_xy(icrm,i,j) + w(icrm,i,j,k500(icrm))*dtfactor
        end do
      end do
    end do

    !$acc parallel loop gang vector collapse(2)
    do k = 1 , nzm
      do icrm = 1 , ncrms
        qv0(icrm,k) = q0 (icrm,k)- qn0(icrm,k)
      enddo
    enddo

    !=====================================================
    ! UW ADDITIONS

    ! FIND VERTICAL INDICES OF 850MB, COMPUTE SWVP
    !$acc parallel loop gang vector
    do icrm = 1 , ncrms
      k850(icrm) = 1
      do k = 1,nzm
        if(pres(icrm,k).le.850.) then
          k850(icrm) = k
          EXIT
        end if
      end do
    end do

    !$acc parallel loop gang vector collapse(4)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            tmp = rho(icrm,k)*dz(icrm)*adz(icrm,k)*dtfactor
            ! Saturated water vapor path with respect to water. Can be used
            ! with water vapor path (= pw) to compute column-average
            ! relative humidity.
            swvp_xy(icrm,i,j) = swvp_xy(icrm,i,j)+qsatw_crm(tabs(icrm,i,j,k),pres(icrm,k))*tmp
          end do
        end do
      end do
    end do ! k

    ! ACCUMULATE AVERAGES OF TWO-DIMENSIONAL STATISTICS
    !$acc parallel loop gang vector collapse(3) copy(psfc_xy,u850_xy,v850_xy) copyin(u,v,pres,p)
    do j=1,ny
      do i=1,nx
        do icrm = 1 , ncrms
          psfc_xy(icrm,i,j) = psfc_xy(icrm,i,j) + (100.*pres(icrm,1) + p(icrm,i,j,1))*dtfactor

          ! 850 mbar horizontal winds
          u850_xy(icrm,i,j) = u850_xy(icrm,i,j) + u(icrm,i,j,k850(icrm))*dtfactor
          v850_xy(icrm,i,j) = v850_xy(icrm,i,j) + v(icrm,i,j,k850(icrm))*dtfactor

        end do
      end do
    end do

    ! COMPUTE CLOUD/ECHO HEIGHTS AS WELL AS CLOUD TOP TEMPERATURE
    ! WHERE CLOUD TOP IS DEFINED AS THE HIGHEST MODEL LEVEL WITH A
    ! CONDENSATE PATH OF 0.01 kg/m2 ABOVE.  ECHO TOP IS THE HIGHEST LEVEL
    ! WHERE THE PRECIPITATE MIXING RATIO > 0.001 G/KG.

    ! initially, zero out heights and set cloudtoptemp to SST
    !$acc parallel loop gang vector collapse(3)
    do j = 1 , ny
      do i = 1 , nx
        do icrm = 1 , ncrms
          cloudtopheight(icrm,i,j) = 0.
          cloudtoptemp(icrm,i,j) = sstxy(icrm,i,j)
          echotopheight(icrm,i,j) = 0.
          ! FIND CLOUD TOP HEIGHT
          tmp_lwp = 0.
          do k = nzm,1,-1
            tmp_lwp = tmp_lwp + (qcl(icrm,i,j,k)+qci(icrm,i,j,k))*rho(icrm,k)*dz(icrm)*adz(icrm,k)
            if (tmp_lwp.gt.0.01) then
              cloudtopheight(icrm,i,j) = z(icrm,k)
              cloudtoptemp(icrm,i,j) = tabs(icrm,i,j,k)
              cld_xy(icrm,i,j) = cld_xy(icrm,i,j) + dtfactor
              EXIT
            end if
          end do
          ! FIND ECHO TOP HEIGHT
          do k = nzm,1,-1
            if (qpl(icrm,i,j,k)+qpi(icrm,i,j,k).gt.1.e-6) then
              echotopheight(icrm,i,j) = z(icrm,k)
              EXIT
            end if
          end do
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

    deallocate(k200)
    deallocate(k500)
    deallocate(k850)
  end subroutine diagnose

end module diagnose_mod
