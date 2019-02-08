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
            dudt(1,j,k,na,icrm) = 0.
          end do
        end do
      end do
    end if

    if(dowally.and.RUN3D.and.rank.lt.nsubdomains_x) then
      !$acc parallel loop collapse(3) copy(dvdt) async(asyncid)
      do icrm = 1 , ncrms
        do k=1,nzm
          do i=1,nx
            dvdt(i,1,k,na,icrm) = 0.
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
            (rdx*(dudt(ic,j,k,na,icrm)-dudt(i,j,k,na,icrm))+ &
            (dwdt(i,j,kc,na,icrm)*rup-dwdt(i,j,k,na,icrm)*rdn) ) + &
            btat*(rdx*(dudt(ic,j,k,nb,icrm)-dudt(i,j,k,nb,icrm))+ &
            (dwdt(i,j,kc,nb,icrm)*rup-dwdt(i,j,k,nb,icrm)*rdn) ) + &
            ctat*(rdx*(dudt(ic,j,k,nc,icrm)-dudt(i,j,k,nc,icrm))+ &
            (dwdt(i,j,kc,nc,icrm)*rup-dwdt(i,j,k,nc,icrm)*rdn) )
            p(i,j,k,icrm)=p(i,j,k,icrm)*rho(k,icrm)
          end do
        end do
      end do

    endif

  end subroutine press_rhs

end module press_rhs_mod
