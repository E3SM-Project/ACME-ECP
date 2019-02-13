module advect_all_scalars_mod
  use advect_scalar_mod
  implicit none

contains

  subroutine advect_all_scalars(ncrms)

    use vars
    use microphysics
    use sgs
    use crmtracers
    use params, only: dotracers
    use scalar_momentum_mod
    implicit none
    integer, intent(in) :: ncrms
    integer k,icrm, i, j, kk

    real(crm_rknd) :: esmt_offset(ncrms)    ! whannah - offset for advecting scalar momentum tracers

    !      advection of scalars :
    call advect_scalar(ncrms,t,tadv,twle)

    !    Advection of microphysics prognostics:
    do k = 1,nmicro_fields
      if(   k.eq.index_water_vapor             &! transport water-vapor variable no metter what
      .or. docloud.and.flag_precip(k).ne.1    & ! transport non-precipitation vars
      .or. doprecip.and.flag_precip(k).eq.1 ) then
        call advect_scalar(ncrms,micro_field(:,:,:,:,k),mkadv(:,:,k),mkwle(:,:,k))
      endif
    end do

    !    Advection of sgs prognostics:
    if(dosgs.and.advect_sgs) then
      do k = 1,nsgs_fields
        call advect_scalar(ncrms,sgs_field(:,:,:,:,k),sgsadv(:,:,k),sgswle(:,:,k))
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
      esmt_offset(icrm) = abs( minval( (/ minval(u_esmt(icrm,:,:,:)), minval(v_esmt(icrm,:,:,:)) /) ) ) + 50.
      u_esmt(icrm,:,:,:) = u_esmt(icrm,:,:,:) + esmt_offset(icrm)
      v_esmt(icrm,:,:,:) = v_esmt(icrm,:,:,:) + esmt_offset(icrm)
    enddo
    ! advection of scalar momentum tracers
    call advect_scalar(ncrms,u_esmt,u_esmt_adv,u_esmt_wle)
    call advect_scalar(ncrms,v_esmt,v_esmt_adv,v_esmt_wle)
    do icrm = 1 , ncrms
      u_esmt(icrm,:,:,:) = u_esmt(icrm,:,:,:) - esmt_offset(icrm)
      v_esmt(icrm,:,:,:) = v_esmt(icrm,:,:,:) - esmt_offset(icrm)
    enddo
#endif

  end subroutine advect_all_scalars

end module advect_all_scalars_mod
