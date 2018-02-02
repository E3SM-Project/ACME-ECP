module forcing_mod
  implicit none

contains

  subroutine forcing(ncrms)

    use vars
    use params
    use microphysics, only: micro_field, index_water_vapor, total_water

    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) :: coef,qneg(ncrms),qpoz(ncrms), factor(ncrms)
    integer i,j,k,nneg(ncrms), icrm

    coef = 1./3600.

    do k=1,nzm
      qpoz(:) = 0.
      qneg(:) = 0.
      nneg(:) = 0

      do j=1,ny
        do i=1,nx
          t(:,i,j,k)=t(:,i,j,k) + ttend(:,k) * dtn
          micro_field(:,i,j,k,index_water_vapor)=micro_field(:,i,j,k,index_water_vapor) + qtend(:,k) * dtn
          qneg(:) = qneg(:) + min( 0._crm_rknd , micro_field(:,i,j,k,index_water_vapor) )
          qpoz(:) = qpoz(:) + max( 0._crm_rknd , micro_field(:,i,j,k,index_water_vapor) )
          dudt(:,i,j,k,na)=dudt(:,i,j,k,na) + utend(:,k)
          dvdt(:,i,j,k,na)=dvdt(:,i,j,k,na) + vtend(:,k)
        end do
      end do

      factor(:) = 1. + qneg(:)/qpoz(:)
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            if(qpoz(icrm)+qneg(icrm).gt.0.) then
              micro_field(icrm,i,j,k,index_water_vapor) = max(real(0.,crm_rknd),micro_field(icrm,i,j,k,index_water_vapor)*factor(icrm))
            end if
          end do
        end do
      end do
    end do

  end subroutine forcing

end module forcing_mod
