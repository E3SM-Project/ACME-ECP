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
    real *8 dta,rdx,rdy,rdz,btat,ctat,rup,rdn
    integer i,j,k,ic,jc,kc, icrm

    if(dowallx.and.mod(rank,nsubdomains_x).eq.0) then
      !$acc parallel loop collapse(3) copy(dudt) async(asyncid)
      do icrm = 1 , ncrms
        do k=1,nzm
          do j=1,ny
            dudt(icrm,1,j,k,na) = 0.
          end do
        end do
      end do
    end if

    if(dowally.and.RUN3D.and.rank.lt.nsubdomains_x) then
      !$acc parallel loop collapse(3) copy(dvdt) async(asyncid)
      do icrm = 1 , ncrms
        do k=1,nzm
          do i=1,nx
            dvdt(icrm,i,1,k,na) = 0.
          end do
        end do
      end do
    end if

    call bound_duvdt(ncrms)

    rdx=1./dx
    rdy=1./dy
    btat=bt/at
    ctat=ct/at

    if(RUN3D) then

      !$acc parallel loop collapse(4) copyin(rhow,u,v,w,dt3,dvdt,dudt,adz,rho,dwdt,dz) copy(p) async(asyncid)
      do icrm = 1 , ncrms
        do k=1,nzm
          do j=1,ny
            do i=1,nx
              kc=k+1
              rdz=1./(adz(k,icrm)*dz(icrm))
              rup = rhow(kc,icrm)/rho(k,icrm)*rdz
              rdn = rhow(k,icrm)/rho(k,icrm)*rdz
              jc=j+1
              ic=i+1
              dta=1./dt3(na)/at
              p(i,j,k,icrm)=(rdx*(u(icrm,ic,j,k)-u(icrm,i,j,k))+ &
              rdy*(v(icrm,i,jc,k)-v(icrm,i,j,k))+ &
              (w(icrm,i,j,kc)*rup-w(icrm,i,j,k)*rdn) )*dta + &
              (rdx*(dudt(icrm,ic,j,k,na)-dudt(icrm,i,j,k,na))+ &
              rdy*(dvdt(icrm,i,jc,k,na)-dvdt(icrm,i,j,k,na))+ &
              (dwdt(icrm,i,j,kc,na)*rup-dwdt(icrm,i,j,k,na)*rdn) ) + &
              btat*(rdx*(dudt(icrm,ic,j,k,nb)-dudt(icrm,i,j,k,nb))+ &
              rdy*(dvdt(icrm,i,jc,k,nb)-dvdt(icrm,i,j,k,nb))+ &
              (dwdt(icrm,i,j,kc,nb)*rup-dwdt(icrm,i,j,k,nb)*rdn) ) + &
              ctat*(rdx*(dudt(icrm,ic,j,k,nc)-dudt(icrm,i,j,k,nc))+ &
              rdy*(dvdt(icrm,i,jc,k,nc)-dvdt(icrm,i,j,k,nc))+ &
              (dwdt(icrm,i,j,kc,nc)*rup-dwdt(icrm,i,j,k,nc)*rdn) )
              p(i,j,k,icrm)=p(i,j,k,icrm)*rho(k,icrm)
            end do
          end do
        end do
      end do

    else

      j=1
      !$acc parallel loop collapse(3) copyin(rhow,u,w,dt3,dudt,adz,rho,dwdt,dz) copy(p) async(asyncid)
      do icrm = 1 , ncrms
        do k=1,nzm
          do i=1,nx
            kc=k+1
            rdz=1./(adz(k,icrm)*dz(icrm))
            rup = rhow(kc,icrm)/rho(k,icrm)*rdz
            rdn = rhow(k,icrm)/rho(k,icrm)*rdz
            ic=i+1
            dta=1./dt3(na)/at
            p(i,j,k,icrm)=(rdx*(u(icrm,ic,j,k)-u(icrm,i,j,k))+ &
            (w(icrm,i,j,kc)*rup-w(icrm,i,j,k)*rdn) )*dta + &
            (rdx*(dudt(icrm,ic,j,k,na)-dudt(icrm,i,j,k,na))+ &
            (dwdt(icrm,i,j,kc,na)*rup-dwdt(icrm,i,j,k,na)*rdn) ) + &
            btat*(rdx*(dudt(icrm,ic,j,k,nb)-dudt(icrm,i,j,k,nb))+ &
            (dwdt(icrm,i,j,kc,nb)*rup-dwdt(icrm,i,j,k,nb)*rdn) ) + &
            ctat*(rdx*(dudt(icrm,ic,j,k,nc)-dudt(icrm,i,j,k,nc))+ &
            (dwdt(icrm,i,j,kc,nc)*rup-dwdt(icrm,i,j,k,nc)*rdn) )
            p(i,j,k,icrm)=p(i,j,k,icrm)*rho(k,icrm)
          end do
        end do
      end do

    endif

  end subroutine press_rhs

end module press_rhs_mod
