module forcing_mod
  implicit none

contains

  subroutine forcing(ncrms,icrm)
    use vars
    use params
    use microphysics, only: micro_field, index_water_vapor, total_water
    implicit none
    integer, intent(in) :: ncrms,icrm
    real(crm_rknd) coef,qneg,qpoz, factor
    integer i,j,k,nneg

    coef = 1./3600.

    do k=1,nzm

      qpoz = 0.
      qneg = 0.
      nneg = 0

      do j=1,ny
        do i=1,nx
          t(i,j,k)=t(i,j,k) + ttend(k) * dtn
          micro_field(i,j,k,index_water_vapor)=micro_field(i,j,k,index_water_vapor) + qtend(k) * dtn
          if(micro_field(i,j,k,index_water_vapor).lt.0.) then
            nneg = nneg + 1
            qneg = qneg + micro_field(i,j,k,index_water_vapor)
          else
            qpoz = qpoz + micro_field(i,j,k,index_water_vapor)
          end if
          dudt(i,j,k,na,icrm)=dudt(i,j,k,na,icrm) + utend(k)
          dvdt(i,j,k,na,icrm)=dvdt(i,j,k,na,icrm) + vtend(k)
        end do
      end do

      if(nneg.gt.0.and.qpoz+qneg.gt.0.) then
        factor = 1. + qneg/qpoz
        do j=1,ny
          do i=1,nx
            micro_field(i,j,k,index_water_vapor) = max(real(0.,crm_rknd),micro_field(i,j,k,index_water_vapor)*factor)
          end do
        end do
      end if

    end do

  end subroutine forcing

end module forcing_mod
