module periodic_mod
  use bound_exchange_mod
  use task_util_mod
  implicit none

contains

  subroutine periodic(flag,ncrms)

    use vars
    use microphysics
    use sgs
    use params, only: dotracers, dosgs
    use crmtracers
#ifdef CLUBB_CRM
    use params, only: doclubb, doclubbnoninter
#endif
    implicit none
    integer, intent(in) :: ncrms

    integer flag, i, j, icrm

    if(flag.eq.0) then

      call bound_exchange(u(:,:,:,:),dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm,1,1,1,1,1,ncrms)
      call bound_exchange(v(:,:,:,:),dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm,1,1,1,1,2,ncrms)
      ! use w at the top level  - 0s anyway - to exchange the sst boundaries (for
      ! surface fluxes call
      !$acc parallel loop gang vector collapse(3)
      do j = 1 , ny
        do i = 1 , nx
          do icrm = 1 , ncrms
            w(icrm,i,j,nz) = sstxy(icrm,i,j)
          enddo
        enddo
      enddo
      call bound_exchange(w(:,:,:,:),dimx1_w,dimx2_w,dimy1_w,dimy2_w,nz,1,1,1,1,3,ncrms)
      !$acc parallel loop gang vector collapse(3)
      do j = 1-YES3D,ny+YES3D
        do i = 1 , nx+1
          do icrm = 1 , ncrms
            if (i <= nx .and. j <= ny) sstxy(icrm,i,j) = w(icrm,i,j,nz)
            w(icrm,i,j,nz) = 0.
          enddo
        enddo
      enddo

    endif


    if(flag.eq.2) then

      call bound_exchange(u(:,:,:,:),dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm,2,3,2+NADV,2+NADV,1,ncrms)
      call bound_exchange(v(:,:,:,:),dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm,2+NADV,2+NADV,2,3,2,ncrms)
      call bound_exchange(w(:,:,:,:),dimx1_w,dimx2_w,dimy1_w,dimy2_w,nz,2+NADV,2+NADV,2+NADV,2+NADV,3,ncrms)

      call bound_exchange(t(:,:,:,:),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3+NADVS,3+NADVS,3+NADVS,3+NADVS,4,ncrms)
      do i = 1,nsgs_fields
        if(dosgs.and.advect_sgs) &
        call bound_exchange(sgs_field(:,:,:,:,i),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3+NADVS,3+NADVS,3+NADVS,3+NADVS,4+i,ncrms)
      end do
      do i = 1,nmicro_fields
        if(   i.eq.index_water_vapor             &
#ifdef CLUBB_CRM
        ! Vince Larson (UWM) changed so that bound_exchange is called even if
        !     docloud = .false. and doclubb = .true.    11 Nov 2007
        .or. (docloud.or.doclubb.or.doclubbnoninter) .and.flag_precip(i).ne.1    &
#else
        .or. docloud.and.flag_precip(i).ne.1    &
#endif
        .or. doprecip.and.flag_precip(i).eq.1 ) &
        call bound_exchange(micro_field(:,:,:,:,i),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3+NADVS,3+NADVS,3+NADVS,3+NADVS,4+nsgs_fields+nsgs_fields_diag+i,ncrms)
      end do
      if(dotracers) then
        do i=1,ntracers
          call bound_exchange(tracer(:,:,:,:,i),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3+NADVS,3+NADVS,3+NADVS,3+NADVS,4+nsgs_fields+nsgs_fields_diag+nmicro_fields+i,ncrms)
        end do
      end if

    endif

    if(flag.eq.3) then

      call bound_exchange(t(:,:,:,:),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,4,ncrms)
      do i = 1,nsgs_fields
        if(dosgs.and.advect_sgs) &
        call bound_exchange(sgs_field(:,:,:,:,i),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,4+i,ncrms)
      end do
      do i = 1,nmicro_fields
        if(   i.eq.index_water_vapor             &
#ifdef CLUBB_CRM
        ! Vince Larson (UWM) changed so that bound_exchange is called even if
        !     docloud = .false. and doclubb = .true.    11 Nov 2007
        .or. (docloud.or.doclubb.or.doclubbnoninter) .and.flag_precip(i).ne.1    &
#else
        .or. docloud.and.flag_precip(i).ne.1    &
#endif
        .or. doprecip.and.flag_precip(i).eq.1 ) &
        call bound_exchange(micro_field(:,:,:,:,i),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,4+nsgs_fields+nsgs_fields_diag+i,ncrms)
      end do
      if(dotracers) then
        do i=1,ntracers
          call bound_exchange(tracer(:,:,:,:,i),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,4+nsgs_fields+nsgs_fields_diag+nmicro_fields+i,ncrms)
        end do
      end if

    endif

    if(flag.eq.4) then

      do i = 1,nsgs_fields_diag
        if(dosgs.and.do_sgsdiag_bound) &
        call bound_exchange(sgs_field_diag(:,:,:,:,i),dimx1_d,dimx2_d,dimy1_d,dimy2_d,nzm,1+dimx1_d,dimx2_d-nx,YES3D+dimy1_d,1-YES3D+dimy2_d-ny,4+nsgs_fields+i,ncrms)
      end do

    end if




  end subroutine periodic

end module periodic_mod
