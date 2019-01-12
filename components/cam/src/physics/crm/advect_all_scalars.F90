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

    real(crm_rknd) :: adv_tmp(nz,ncrms)
    real(crm_rknd) :: wle_tmp(nz,ncrms)
    real(crm_rknd) :: esmt_offset(ncrms)    ! whannah - offset for advecting scalar momentum tracers

    !$acc enter data create(adv_tmp,wle_tmp) async(1)

    !      advection of scalars :
    call advect_scalar(ncrms,t,tadv,twle)

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
        call advect_scalar(ncrms,micro_field(:,:,:,:,k),adv_tmp,wle_tmp)
        !$acc parallel loop collapse(2) copyin(adv_tmp,wle_tmp) copy(mkadv,mkwle) async(1)
        do icrm = 1 , ncrms
          do kk = 1 , nz
            mkadv(kk,k,icrm) = adv_tmp(kk,icrm)
            mkwle(kk,k,icrm) = wle_tmp(kk,icrm)
          enddo
        enddo
      endif
    end do

    !    Advection of sgs prognostics:
    if(dosgs.and.advect_sgs) then
      do k = 1,nsgs_fields
        call advect_scalar(ncrms,sgs_field(:,:,:,:,k),adv_tmp,wle_tmp)
        !$acc parallel loop collapse(2) copyin(adv_tmp,wle_tmp) copy(sgsadv,sgswle) async(1)
        do icrm = 1 , ncrms
          do kk = 1 , nz
            sgsadv(kk,k,icrm) = adv_tmp(kk,icrm)
            sgswle(kk,k,icrm) = wle_tmp(kk,icrm)
          enddo
        enddo
      end do
    end if

    !   Precipitation fallout:
    if(doprecip) then
      !do icrm = 1 , ncrms
      !  total_water_prec(icrm) = total_water_prec(icrm) + total_water(ncrms,icrm)
      !enddo
      call micro_precip_fall(ncrms)
      !do icrm = 1 , ncrms
      !  total_water_prec(icrm) = total_water_prec(icrm) - total_water(ncrms,icrm)
      !enddo
    end if

    !$acc exit data delete(adv_tmp,wle_tmp) async(1)

    ! advection of tracers:
    !There aren't any of these. We need to delete crmtracers.F90 too at some point
    !if(dotracers) then
    !  do k = 1,ntracers
    !    call advect_scalar(ncrms,icrm,tracer(:,:,:,k,icrm),tradv(:,k,icrm),trwle(:,k,icrm))
    !  end do
    !end if

#if defined(SP_ESMT)
    ! whannah - the esmt_offset simply ensures that the scalar momentum
    ! tracers are positive definite during the advection calculation
    do icrm = 1 , ncrms
      esmt_offset(icrm) = abs( minval( (/ minval(u_esmt(:,:,:,icrm)), minval(v_esmt(:,:,:,icrm)) /) ) ) + 50.
      u_esmt(:,:,:,icrm) = u_esmt(:,:,:,icrm) + esmt_offset(icrm)
      v_esmt(:,:,:,icrm) = v_esmt(:,:,:,icrm) + esmt_offset(icrm)
    enddo
    ! advection of scalar momentum tracers
    call advect_scalar(ncrms,u_esmt,u_esmt_adv,u_esmt_wle)
    call advect_scalar(ncrms,v_esmt,v_esmt_adv,v_esmt_wle)
    do icrm = 1 , ncrms
      u_esmt(:,:,:,icrm) = u_esmt(:,:,:,icrm) - esmt_offset(icrm)
      v_esmt(:,:,:,icrm) = v_esmt(:,:,:,icrm) - esmt_offset(icrm)
    enddo
#endif

  end subroutine advect_all_scalars

end module advect_all_scalars_mod
