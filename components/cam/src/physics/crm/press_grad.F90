module press_grad_mod
  use task_util_mod
  use bound_duvdt_mod
  implicit none

contains

  subroutine press_grad(ncrms)
    !       pressure term of the momentum equations
    use vars
    use params, only: dowallx, dowally
    implicit none
    integer, intent(in) :: ncrms
    real *8 rdx,rdy,rdz
    integer i,j,k,kb,jb,ib, icrm

    rdx=1./dx
    rdy=1./dy

    !$acc parallel loop collapse(4) copyin(p) copy(dudt,dvdt,dwdt) async(1)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            kb=max(1,k-1)
            rdz = 1./(dz(icrm)*adzw(k,icrm))
            jb=j-YES3D
            ib=i-1
            dudt(i,j,k,na,icrm)=dudt(i,j,k,na,icrm)-(p(i,j,k,icrm)-p(ib,j,k,icrm))*rdx
            dvdt(i,j,k,na,icrm)=dvdt(i,j,k,na,icrm)-(p(i,j,k,icrm)-p(i,jb,k,icrm))*rdy
            dwdt(i,j,k,na,icrm)=dwdt(i,j,k,na,icrm)-(p(i,j,k,icrm)-p(i,j,kb,icrm))*rdz
          end do ! i
        end do ! j
      end do ! k
    enddo

    !$acc parallel loop collapse(4) copy(p) copyin(rho) async(1)
    do icrm = 1 , ncrms
      do k=1,nzm
        do j=1-YES3D,ny !bloss: 0,n* fixes computation of dp/d* in stats.
          do i=0,nx
            p(i,j,k,icrm)=p(i,j,k,icrm)*rho(k,icrm)  ! convert p'/rho to p'
          end do
        end do
      end do
    enddo

    if(dowallx.and.mod(rank,nsubdomains_x).eq.0) then

      !$acc parallel loop collapse(3) async(1)
      do icrm = 1 , ncrms
        do k=1,nzm
          do j=1,ny
            dudt(1,j,k,na,icrm) = 0.
          end do
        end do
      enddo

    end if

    if(dowally.and.RUN3D.and.rank.lt.nsubdomains_x) then

      !$acc parallel loop collapse(3) async(1)
      do icrm = 1 , ncrms
        do k=1,nzm
          do i=1,nx
            dvdt(i,1,k,na,icrm) = 0.
          end do
        end do
      enddo

    end if

    call bound_duvdt(ncrms)

  end subroutine press_grad

end module press_grad_mod
