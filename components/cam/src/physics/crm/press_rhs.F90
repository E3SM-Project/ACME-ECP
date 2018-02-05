module press_rhs_mod
  use task_util_mod
  use bound_duvdt_mod
  implicit none

contains

  subroutine press_rhs(ncrms)

    !       right-hand-side of the Poisson equation for pressure

    use vars
    use params, only: dowallx, dowally

    implicit none
    integer, intent(in) :: ncrms


    real *8 dta,rdx,rdy,rdz(ncrms),btat,ctat,rup(ncrms),rdn(ncrms)
    integer i,j,k,ic,jc,kc,icrm

    if(dowallx.and.mod(rank,nsubdomains_x).eq.0) then

      do k=1,nzm
        do j=1,ny
          dudt(:,1,j,k,na) = 0.
        end do
      end do

    end if

    if(dowally.and.RUN3D.and.rank.lt.nsubdomains_x) then

      do k=1,nzm
        do i=1,nx
          dvdt(:,i,1,k,na) = 0.
        end do
      end do

    end if

    call bound_duvdt(ncrms)

    dta=1./dt3(na)/at
    rdx=1./dx
    rdy=1./dy
    btat=bt/at
    ctat=ct/at

    if(RUN3D) then

      do k=1,nzm
        kc=k+1
        rdz(:)=1./(adz(:,k)*dz(:))
        rup(:) = rhow(:,kc)/rho(:,k)*rdz(:)
        rdn(:) = rhow(:,k)/rho(:,k)*rdz(:)
        do j=1,ny
          jc=j+1
          do i=1,nx
            do icrm = 1 , ncrms
              ic=i+1
              p(icrm,i,j,k)=(rdx*(u(icrm,ic,j,k)-u(icrm,i,j,k))+ &
              rdy*(v(icrm,i,jc,k)-v(icrm,i,j,k))+ &
              (w(icrm,i,j,kc)*rup(icrm)-w(icrm,i,j,k)*rdn(icrm)) )*dta + &
              (rdx*(dudt(icrm,ic,j,k,na)-dudt(icrm,i,j,k,na))+ &
              rdy*(dvdt(icrm,i,jc,k,na)-dvdt(icrm,i,j,k,na))+ &
              (dwdt(icrm,i,j,kc,na)*rup(icrm)-dwdt(icrm,i,j,k,na)*rdn(icrm)) ) + &
              btat*(rdx*(dudt(icrm,ic,j,k,nb)-dudt(icrm,i,j,k,nb))+ &
              rdy*(dvdt(icrm,i,jc,k,nb)-dvdt(icrm,i,j,k,nb))+ &
              (dwdt(icrm,i,j,kc,nb)*rup(icrm)-dwdt(icrm,i,j,k,nb)*rdn(icrm)) ) + &
              ctat*(rdx*(dudt(icrm,ic,j,k,nc)-dudt(icrm,i,j,k,nc))+ &
              rdy*(dvdt(icrm,i,jc,k,nc)-dvdt(icrm,i,j,k,nc))+ &
              (dwdt(icrm,i,j,kc,nc)*rup(icrm)-dwdt(icrm,i,j,k,nc)*rdn(icrm)) )
              p(icrm,i,j,k)=p(icrm,i,j,k)*rho(icrm,k)
            end do
          end do
        end do
      end do


    else

      j=1

      do k=1,nzm
        kc=k+1
        rdz(:)=1./(adz(:,k)*dz(:))
        rup(:) = rhow(:,kc)/rho(:,k)*rdz(:)
        rdn(:) = rhow(:,k)/rho(:,k)*rdz(:)
        do i=1,nx
          do icrm = 1 , ncrms
            ic=i+1
            p(icrm,i,j,k)=(rdx*(u(icrm,ic,j,k)-u(icrm,i,j,k))+ &
            (w(icrm,i,j,kc)*rup(icrm)-w(icrm,i,j,k)*rdn(icrm)) )*dta + &
            (rdx*(dudt(icrm,ic,j,k,na)-dudt(icrm,i,j,k,na))+ &
            (dwdt(icrm,i,j,kc,na)*rup(icrm)-dwdt(icrm,i,j,k,na)*rdn(icrm)) ) + &
            btat*(rdx*(dudt(icrm,ic,j,k,nb)-dudt(icrm,i,j,k,nb))+ &
            (dwdt(icrm,i,j,kc,nb)*rup(icrm)-dwdt(icrm,i,j,k,nb)*rdn(icrm)) ) + &
            ctat*(rdx*(dudt(icrm,ic,j,k,nc)-dudt(icrm,i,j,k,nc))+ &
            (dwdt(icrm,i,j,kc,nc)*rup(icrm)-dwdt(icrm,i,j,k,nc)*rdn(icrm)) )
            p(icrm,i,j,k)=p(icrm,i,j,k)*rho(icrm,k)
          end do
        end do
      end do


    endif

    call task_barrier()

  end subroutine press_rhs

end module press_rhs_mod
