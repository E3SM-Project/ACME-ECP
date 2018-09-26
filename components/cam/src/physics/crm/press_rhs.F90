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
          dudt(1,j,k,na(icrm),icrm) = 0.
        end do
      end do

    end if

    if(dowally.and.RUN3D.and.rank.lt.nsubdomains_x) then

      do k=1,nzm
        do i=1,nx
          dvdt(i,1,k,na(icrm),icrm) = 0.
        end do
      end do

    end if


    if(dompi) then
      call task_bound_duvdt()
    else
      call bound_duvdt(ncrms,icrm)
    endif

    dta=1./dt3(na(icrm),icrm)/at
    rdx=1./dx
    rdy=1./dy
    btat=bt/at
    ctat=ct/at

    if(RUN3D) then

      do k=1,nzm
        kc=k+1
        rdz=1./(adz(k,icrm)*dz(icrm))
        rup = rhow(kc,icrm)/rho(k,icrm)*rdz
        rdn = rhow(k,icrm)/rho(k,icrm)*rdz
        do j=1,ny
          jc=j+1
          do i=1,nx
            ic=i+1
            p(i,j,k,icrm)=(rdx*(u(ic,j,k,icrm)-u(i,j,k,icrm))+ &
            rdy*(v(i,jc,k,icrm)-v(i,j,k,icrm))+ &
            (w(i,j,kc,icrm)*rup-w(i,j,k,icrm)*rdn) )*dta + &
            (rdx*(dudt(ic,j,k,na(icrm),icrm)-dudt(i,j,k,na(icrm),icrm))+ &
            rdy*(dvdt(i,jc,k,na(icrm),icrm)-dvdt(i,j,k,na(icrm),icrm))+ &
            (dwdt(i,j,kc,na(icrm),icrm)*rup-dwdt(i,j,k,na(icrm),icrm)*rdn) ) + &
            btat*(rdx*(dudt(ic,j,k,nb(icrm),icrm)-dudt(i,j,k,nb(icrm),icrm))+ &
            rdy*(dvdt(i,jc,k,nb(icrm),icrm)-dvdt(i,j,k,nb(icrm),icrm))+ &
            (dwdt(i,j,kc,nb(icrm),icrm)*rup-dwdt(i,j,k,nb(icrm),icrm)*rdn) ) + &
            ctat*(rdx*(dudt(ic,j,k,nc(icrm),icrm)-dudt(i,j,k,nc(icrm),icrm))+ &
            rdy*(dvdt(i,jc,k,nc(icrm),icrm)-dvdt(i,j,k,nc(icrm),icrm))+ &
            (dwdt(i,j,kc,nc(icrm),icrm)*rup-dwdt(i,j,k,nc(icrm),icrm)*rdn) )
            p(i,j,k,icrm)=p(i,j,k,icrm)*rho(k,icrm)
          end do
        end do
      end do


    else

      j=1

      do k=1,nzm
        kc=k+1
        rdz=1./(adz(k,icrm)*dz(icrm))
        rup = rhow(kc,icrm)/rho(k,icrm)*rdz
        rdn = rhow(k,icrm)/rho(k,icrm)*rdz
        do i=1,nx
          ic=i+1
          p(i,j,k,icrm)=(rdx*(u(ic,j,k,icrm)-u(i,j,k,icrm))+ &
          (w(i,j,kc,icrm)*rup-w(i,j,k,icrm)*rdn) )*dta + &
          (rdx*(dudt(ic,j,k,na(icrm),icrm)-dudt(i,j,k,na(icrm),icrm))+ &
          (dwdt(i,j,kc,na(icrm),icrm)*rup-dwdt(i,j,k,na(icrm),icrm)*rdn) ) + &
          btat*(rdx*(dudt(ic,j,k,nb(icrm),icrm)-dudt(i,j,k,nb(icrm),icrm))+ &
          (dwdt(i,j,kc,nb(icrm),icrm)*rup-dwdt(i,j,k,nb(icrm),icrm)*rdn) ) + &
          ctat*(rdx*(dudt(ic,j,k,nc(icrm),icrm)-dudt(i,j,k,nc(icrm),icrm))+ &
          (dwdt(i,j,kc,nc(icrm),icrm)*rup-dwdt(i,j,k,nc(icrm),icrm)*rdn) )
          p(i,j,k,icrm)=p(i,j,k,icrm)*rho(k,icrm)
        end do
      end do


    endif

    call task_barrier()

  end subroutine press_rhs

end module press_rhs_mod
