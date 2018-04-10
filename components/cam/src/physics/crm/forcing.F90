module forcing_mod
  implicit none

contains

  subroutine forcing(ncrms)

    use vars
    use params
    use microphysics, only: micro_field, index_water_vapor, total_water
    use openacc_pool
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) :: coef
    real(crm_rknd), pointer :: qneg(:,:),qpoz(:,:), factor(:,:)
    integer i,j,k,icrm

    ! allocate(qneg  (ncrms,nzm))
    ! allocate(qpoz  (ncrms,nzm))
    ! allocate(factor(ncrms,nzm))
    call pool_push(qneg,(/ncrms,nzm/))
    call pool_push(qpoz,(/ncrms,nzm/))
    call pool_push(factor,(/ncrms,nzm/))

    coef = 1./3600.

    !$acc parallel loop gang vector collapse(4)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            t(icrm,i,j,k)=t(icrm,i,j,k) + ttend(icrm,k) * dtn
            micro_field(icrm,i,j,k,index_water_vapor)=micro_field(icrm,i,j,k,index_water_vapor) + qtend(icrm,k) * dtn
            dudt(icrm,i,j,k,na)=dudt(icrm,i,j,k,na) + utend(icrm,k)
            dvdt(icrm,i,j,k,na)=dvdt(icrm,i,j,k,na) + vtend(icrm,k)
            if (j == 1 .and. i == 1) then
              qpoz(icrm,k) = 0.
              qneg(icrm,k) = 0.
            endif
          end do
        end do
      end do
    enddo

    !$acc parallel loop gang vector collapse(2)
    do k=1,nzm
      do icrm = 1 , ncrms
        !$acc loop seq
        do j=1,ny
          !$acc loop seq
          do i=1,nx
            qneg(icrm,k) = qneg(icrm,k) + min( 0._crm_rknd , micro_field(icrm,i,j,k,index_water_vapor) )
            qpoz(icrm,k) = qpoz(icrm,k) + max( 0._crm_rknd , micro_field(icrm,i,j,k,index_water_vapor) )
          enddo
        enddo
        factor(icrm,k) = 1. + qneg(icrm,k)/qpoz(icrm,k)
      enddo
    enddo

    !$acc parallel loop gang vector collapse(4)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            if(qpoz(icrm,k)+qneg(icrm,k).gt.0.) then
              micro_field(icrm,i,j,k,index_water_vapor) = max(real(0.,crm_rknd),micro_field(icrm,i,j,k,index_water_vapor)*factor(icrm,k))
            end if
          end do
        end do
      end do
    end do

    ! deallocate(qneg  )
    ! deallocate(qpoz  )
    ! deallocate(factor)
    call pool_pop_multiple(3)

  end subroutine forcing

end module forcing_mod
