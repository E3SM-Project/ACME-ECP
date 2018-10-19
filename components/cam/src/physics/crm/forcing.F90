
module forcing_mod
  implicit none

contains

  subroutine forcing(ncrms)
    use vars
    use params
    use microphysics, only: micro_field, index_water_vapor, total_water
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) coef,qneg(nzm,ncrms),qpoz(nzm,ncrms), factor
    integer i,j,k,nneg(nzm,ncrms),icrm

    coef = 1./3600.

    do icrm = 1 , ncrms
      do k=1,nzm
        qpoz(k,icrm) = 0.
        qneg(k,icrm) = 0.
        nneg(k,icrm) = 0
      enddo
    enddo

    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            t(i,j,k,icrm)=t(i,j,k,icrm) + ttend(k,icrm) * dtn
            micro_field(i,j,k,index_water_vapor,icrm)=micro_field(i,j,k,index_water_vapor,icrm) + qtend(k,icrm) * dtn
            if(micro_field(i,j,k,index_water_vapor,icrm).lt.0.) then
              nneg(k,icrm) = nneg(k,icrm) + 1
              qneg(k,icrm) = qneg(k,icrm) + micro_field(i,j,k,index_water_vapor,icrm)
            else
              qpoz(k,icrm) = qpoz(k,icrm) + micro_field(i,j,k,index_water_vapor,icrm)
            end if
            dudt(i,j,k,na,icrm)=dudt(i,j,k,na,icrm) + utend(k,icrm)
            dvdt(i,j,k,na,icrm)=dvdt(i,j,k,na,icrm) + vtend(k,icrm)
          end do
        end do
      end do
    end do

    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            if(nneg(k,icrm).gt.0.and.qpoz(k,icrm)+qneg(k,icrm).gt.0.) then
              factor = 1. + qneg(k,icrm)/qpoz(k,icrm)
              micro_field(i,j,k,index_water_vapor,icrm) = max(real(0.,crm_rknd),micro_field(i,j,k,index_water_vapor,icrm)*factor)
            end if
          end do
        end do
      end do
    end do

  end subroutine forcing

end module forcing_mod
