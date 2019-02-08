
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

    !$acc enter data create(qneg,qpoz,nneg) async(asyncid)

    coef = 1./3600.

    !$acc parallel loop collapse(2) copyout(qpoz,qneg,nneg) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        qpoz(k,icrm) = 0.
        qneg(k,icrm) = 0.
        nneg(k,icrm) = 0
      enddo
    enddo

    !$acc parallel loop collapse(4) copy(t,micro_field,nneg,qneg,qpoz,dudt,dvdt) copyin(ttend,qtend,utend,vtend) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            t(icrm,i,j,k)=t(icrm,i,j,k) + ttend(k,icrm) * dtn
            micro_field(icrm,i,j,k,index_water_vapor)=micro_field(icrm,i,j,k,index_water_vapor) + qtend(k,icrm) * dtn
            if(micro_field(icrm,i,j,k,index_water_vapor).lt.0.) then
              !$acc atomic update
              nneg(k,icrm) = nneg(k,icrm) + 1
              !$acc atomic update
              qneg(k,icrm) = qneg(k,icrm) + micro_field(icrm,i,j,k,index_water_vapor)
            else
              !$acc atomic update
              qpoz(k,icrm) = qpoz(k,icrm) + micro_field(icrm,i,j,k,index_water_vapor)
            end if
            dudt(icrm,i,j,k,na)=dudt(icrm,i,j,k,na) + utend(k,icrm)
            dvdt(icrm,i,j,k,na)=dvdt(icrm,i,j,k,na) + vtend(k,icrm)
          end do
        end do
      end do
    end do

    !$acc parallel loop collapse(4) private(factor) copy(micro_field) copyin(qneg,qpoz,nneg) async(asyncid)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            if(nneg(k,icrm).gt.0.and.qpoz(k,icrm)+qneg(k,icrm).gt.0.) then
              factor = 1. + qneg(k,icrm)/qpoz(k,icrm)
              micro_field(icrm,i,j,k,index_water_vapor) = max(real(0.,crm_rknd),micro_field(icrm,i,j,k,index_water_vapor)*factor)
            end if
          end do
        end do
      end do
    end do

    !$acc exit data delete(qneg,qpoz,nneg) async(asyncid)

  end subroutine forcing

end module forcing_mod
