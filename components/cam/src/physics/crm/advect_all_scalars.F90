module advect_all_scalars_mod
  use advect_scalar_mod
  implicit none

contains

  subroutine advect_all_scalars(ncrms,icrm)

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
    integer, intent(in) :: ncrms,icrm
    ! real dummy(nz)
    real(crm_rknd) dummy(nz)
    integer k

    real(crm_rknd) esmt_offset    ! whannah - offset for advecting scalar momentum tracers


    !---------------------------------------------------------
    !      advection of scalars :

    call advect_scalar(ncrms,icrm,t,tadv,twle(:,icrm),t2leadv,t2legrad,twleadv,.true.)

    !
    !    Advection of microphysics prognostics:
    !

    do k = 1,nmicro_fields
      if(   k.eq.index_water_vapor             &! transport water-vapor variable no metter what
#ifdef CLUBB_CRM
      !Added preprocessor directives. - nielsenb UWM 30 July 2008
      .or. ( docloud .or. doclubb .or. doclubbnoninter ) .and.flag_precip(k).ne.1    & ! transport non-precipitation vars
#else
      .or. docloud.and.flag_precip(k).ne.1    & ! transport non-precipitation vars
#endif
      .or. doprecip.and.flag_precip(k).eq.1 ) &
      call advect_scalar(ncrms,icrm,micro_field(:,:,:,k),mkadv(:,k),mkwle(:,k),dummy,dummy,dummy,.false.)
    end do

    !
    !    Advection of sgs prognostics:
    !

    if(dosgs.and.advect_sgs) then
      do k = 1,nsgs_fields
        call advect_scalar(ncrms,icrm,sgs_field(:,:,:,k),sgsadv(:,k),sgswle(:,k),dummy,dummy,dummy,.false.)
      end do
    end if


    !
    !   Precipitation fallout:
    !
    if(doprecip) then

      total_water_prec = total_water_prec + total_water(ncrms,icrm)

      call micro_precip_fall(ncrms,icrm)

      total_water_prec = total_water_prec - total_water(ncrms,icrm)


    end if

    ! advection of tracers:

    if(dotracers) then

      do k = 1,ntracers
        call advect_scalar(ncrms,icrm,tracer(:,:,:,k),tradv(:,k),trwle(:,k),dummy,dummy,dummy,.false.)
      end do

    end if

#if defined(SP_ESMT)
    
    ! whannah - the esmt_offset simply ensures that the scalar momentum  
    ! tracers are positive definite during the advection calculation
    ! esmt_offset = 1000.

    esmt_offset = abs( minval( (/ minval(u_esmt), minval(v_esmt) /) ) ) + 50.

    u_esmt(:,:,:) = u_esmt(:,:,:) + esmt_offset
    v_esmt(:,:,:) = v_esmt(:,:,:) + esmt_offset

    ! advection of scalar momentum tracers
    call advect_scalar(ncrms,icrm,u_esmt,u_esmt_adv,u_esmt_wle,dummy,dummy,dummy,.false.)
    call advect_scalar(ncrms,icrm,v_esmt,v_esmt_adv,v_esmt_wle,dummy,dummy,dummy,.false.)

    u_esmt(:,:,:) = u_esmt(:,:,:) - esmt_offset
    v_esmt(:,:,:) = v_esmt(:,:,:) - esmt_offset

#endif

  end subroutine advect_all_scalars

end module advect_all_scalars_mod
