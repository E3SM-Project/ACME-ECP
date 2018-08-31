module press_rhs_mod
  use task_util_mod
  use bound_duvdt_mod
  implicit none

contains

  subroutine press_rhs(ncrms,icrm)
    !       right-hand-side of the Poisson equation for pressure
    use vars
    use params, only: dowallx, dowally
    implicit none
    integer, intent(in) :: ncrms, icrm
    real *8 dta,rdx,rdy,rdz,btat,ctat,rup,rdn
    integer i,j,k,ic,jc,kc

    if(dowallx.and.mod(rank,nsubdomains_x).eq.0) then

      do k=1,nzm
        do j=1,ny
          dudt(1,j,k,na,icrm) = 0.
        end do
      end do

    end if

    if(dowally.and.RUN3D.and.rank.lt.nsubdomains_x) then

      do k=1,nzm
        do i=1,nx
          dvdt(i,1,k,na,icrm) = 0.
        end do
      end do

    end if


    if(dompi) then
      call task_bound_duvdt()
    else
      call bound_duvdt(ncrms,icrm)
    endif

    dta=1./dt3(na)/at
    rdx=1./dx
    rdy=1./dy
    btat=bt/at
    ctat=ct/at

    if(RUN3D) then

      do k=1,nzm
        kc=k+1
        rdz=1./(adz(k)*dz(icrm))
        rup = rhow(kc,icrm)/rho(k,icrm)*rdz
        rdn = rhow(k,icrm)/rho(k,icrm)*rdz
        do j=1,ny
          jc=j+1
          do i=1,nx
            ic=i+1
            p(i,j,k,icrm)=(rdx*(u(ic,j,k)-u(i,j,k))+ &
            rdy*(v(i,jc,k)-v(i,j,k))+ &
            (w(i,j,kc)*rup-w(i,j,k)*rdn) )*dta + &
            (rdx*(dudt(ic,j,k,na,icrm)-dudt(i,j,k,na,icrm))+ &
            rdy*(dvdt(i,jc,k,na,icrm)-dvdt(i,j,k,na,icrm))+ &
            (dwdt(i,j,kc,na,icrm)*rup-dwdt(i,j,k,na,icrm)*rdn) ) + &
            btat*(rdx*(dudt(ic,j,k,nb,icrm)-dudt(i,j,k,nb,icrm))+ &
            rdy*(dvdt(i,jc,k,nb,icrm)-dvdt(i,j,k,nb,icrm))+ &
            (dwdt(i,j,kc,nb,icrm)*rup-dwdt(i,j,k,nb,icrm)*rdn) ) + &
            ctat*(rdx*(dudt(ic,j,k,nc,icrm)-dudt(i,j,k,nc,icrm))+ &
            rdy*(dvdt(i,jc,k,nc,icrm)-dvdt(i,j,k,nc,icrm))+ &
            (dwdt(i,j,kc,nc,icrm)*rup-dwdt(i,j,k,nc,icrm)*rdn) )
            p(i,j,k,icrm)=p(i,j,k,icrm)*rho(k,icrm)
          end do
        end do
      end do


    else

      j=1

      do k=1,nzm
        kc=k+1
        rdz=1./(adz(k)*dz(icrm))
        rup = rhow(kc,icrm)/rho(k,icrm)*rdz
        rdn = rhow(k,icrm)/rho(k,icrm)*rdz
        do i=1,nx
          ic=i+1
          p(i,j,k,icrm)=(rdx*(u(ic,j,k)-u(i,j,k))+ &
          (w(i,j,kc)*rup-w(i,j,k)*rdn) )*dta + &
          (rdx*(dudt(ic,j,k,na,icrm)-dudt(i,j,k,na,icrm))+ &
          (dwdt(i,j,kc,na,icrm)*rup-dwdt(i,j,k,na,icrm)*rdn) ) + &
          btat*(rdx*(dudt(ic,j,k,nb,icrm)-dudt(i,j,k,nb,icrm))+ &
          (dwdt(i,j,kc,nb,icrm)*rup-dwdt(i,j,k,nb,icrm)*rdn) ) + &
          ctat*(rdx*(dudt(ic,j,k,nc,icrm)-dudt(i,j,k,nc,icrm))+ &
          (dwdt(i,j,kc,nc,icrm)*rup-dwdt(i,j,k,nc,icrm)*rdn) )
          p(i,j,k,icrm)=p(i,j,k,icrm)*rho(k,icrm)
        end do
      end do


    endif

    call task_barrier()

  end subroutine press_rhs

end module press_rhs_mod
