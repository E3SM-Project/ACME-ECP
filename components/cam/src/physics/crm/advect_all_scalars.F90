module advect_all_scalars_mod
  use advect_scalar_mod
  implicit none

contains

  subroutine advect_all_scalars(ncrms)

    use vars
    use microphysics
    use sgs
    use crmtracers
#ifdef CLUBB_CRM
    use params, only: dotracers, doclubb, doclubbnoninter
#else
    use params, only: dotracers
#endif
    use scalar_momentum_mod
    implicit none
    integer, intent(in) :: ncrms
    integer k,icrm, i, j, kk
    real(crm_rknd) :: micro_field_tmp(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, ncrms)
    real(crm_rknd) :: adv_tmp(nz,ncrms)
    real(crm_rknd) :: wle_tmp(nz,ncrms)
    real(crm_rknd) :: esmt_offset    ! whannah - offset for advecting scalar momentum tracers

    do icrm = 1 , ncrms

    !      advection of scalars :
    call advect_scalar(ncrms,icrm,t(:,:,:,icrm),tadv(:,icrm),twle(:,icrm))

    !    Advection of microphysics prognostics:
    do k = 1,nmicro_fields
      if(   k.eq.index_water_vapor             &! transport water-vapor variable no metter what
#ifdef CLUBB_CRM
      !Added preprocessor directives. - nielsenb UWM 30 July 2008
      .or. ( docloud .or. doclubb .or. doclubbnoninter ) .and.flag_precip(k).ne.1    & ! transport non-precipitation vars
#else
      .or. docloud.and.flag_precip(k).ne.1    & ! transport non-precipitation vars
#endif
      .or. doprecip.and.flag_precip(k).eq.1 ) then
        do kk = 1 , nzm
          do j = dimy1_s,dimy2_s
            do i = dimx1_s,dimx2_s
              micro_field_tmp(i,j,kk,icrm) = micro_field(i,j,kk,k,icrm)
            enddo
          enddo
        enddo
        do kk = 1 , nz
          adv_tmp(kk,icrm) = mkadv(kk,k,icrm)
          wle_tmp(kk,icrm) = mkwle(kk,k,icrm)
        enddo
        call advect_scalar(ncrms,icrm,micro_field_tmp(:,:,:,icrm),adv_tmp(:,icrm),wle_tmp(:,icrm))
        do kk = 1 , nzm
          do j = dimy1_s,dimy2_s
            do i = dimx1_s,dimx2_s
              micro_field(i,j,kk,k,icrm) = micro_field_tmp(i,j,kk,icrm)
            enddo
          enddo
        enddo
        do kk = 1 , nz
          mkadv(kk,k,icrm) = adv_tmp(kk,icrm)
          mkwle(kk,k,icrm) = wle_tmp(kk,icrm)
        enddo
      endif
    end do

    !    Advection of sgs prognostics:
    if(dosgs.and.advect_sgs) then
      do k = 1,nsgs_fields
        do kk = 1 , nz
          adv_tmp(kk,icrm) = sgsadv(kk,k,icrm)
          wle_tmp(kk,icrm) = sgswle(kk,k,icrm)
        enddo
        call advect_scalar(ncrms,icrm,sgs_field(:,:,:,icrm,k),adv_tmp(:,icrm),wle_tmp(:,icrm))
        do kk = 1 , nz
          sgsadv(kk,k,icrm) = adv_tmp(kk,icrm)
          sgswle(kk,k,icrm) = wle_tmp(kk,icrm)
        enddo
      end do
    end if

    !   Precipitation fallout:
    if(doprecip) then
      total_water_prec(icrm) = total_water_prec(icrm) + total_water(ncrms,icrm)
      call micro_precip_fall(ncrms,icrm)
      total_water_prec(icrm) = total_water_prec(icrm) - total_water(ncrms,icrm)
    end if

    ! advection of tracers:
    !There aren't any of these
    !if(dotracers) then
    !  do k = 1,ntracers
    !    call advect_scalar(ncrms,icrm,tracer(:,:,:,k,icrm),tradv(:,k,icrm),trwle(:,k,icrm))
    !  end do
    !end if

#if defined(SP_ESMT)
    ! whannah - the esmt_offset simply ensures that the scalar momentum
    ! tracers are positive definite during the advection calculation
    esmt_offset = abs( minval( (/ minval(u_esmt(:,:,:,icrm)), minval(v_esmt(:,:,:,icrm)) /) ) ) + 50.
    u_esmt(:,:,:,icrm) = u_esmt(:,:,:,icrm) + esmt_offset
    v_esmt(:,:,:,icrm) = v_esmt(:,:,:,icrm) + esmt_offset
    ! advection of scalar momentum tracers
    call advect_scalar(ncrms,icrm,u_esmt(:,:,:,icrm),u_esmt_adv(:,icrm),u_esmt_wle(:,icrm))
    call advect_scalar(ncrms,icrm,v_esmt(:,:,:,icrm),v_esmt_adv(:,icrm),v_esmt_wle(:,icrm))
    u_esmt(:,:,:,icrm) = u_esmt(:,:,:,icrm) - esmt_offset
    v_esmt(:,:,:,icrm) = v_esmt(:,:,:,icrm) - esmt_offset
#endif
    enddo

  end subroutine advect_all_scalars

end module advect_all_scalars_mod
