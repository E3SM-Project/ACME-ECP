
module kurant_mod
   use task_util_mod
   implicit none

   contains

   subroutine kurant(ncrms)
      use vars
      use sgs, only: kurant_sgs
      use params, only: crm_rknd
      implicit none
      integer, intent(in) :: ncrms
      integer i, j, k, ncycle1(1),ncycle2(1),icrm
      real(crm_rknd) wm(nz,ncrms)  ! maximum vertical wind velocity
      real(crm_rknd) uhm(nz,ncrms) ! maximum horizontal wind velocity
      real(crm_rknd) cfl, cfl_sgs, tmp

      ncycle = 1
      !$acc parallel loop collapse(2) copyout(wm,uhm) async(1)
      do icrm = 1 , ncrms
        do k = 1 , nz
          wm (k,icrm) = 0.
          uhm(k,icrm) = 0.
        enddo
      enddo

      !$acc parallel loop collapse(4) private(tmp) copy(wm,w_max,uhm,u_max) copyin(u,v,w) async(1)
      do icrm = 1 , ncrms
        do k = 1,nzm
          do j = 1 , ny
            do i = 1 , nx
              tmp = abs(w(i,j,k,icrm))
              !$acc atomic update
              wm(k,icrm) = max( wm(k,icrm) , tmp )
              !$acc atomic update
              w_max(icrm) = max( w_max(icrm) , tmp )

              tmp = sqrt(u(i,j,k,icrm)**2+YES3D*v(i,j,k,icrm)**2)
              !$acc atomic update
              uhm(k,icrm) = max( uhm(k,icrm) , tmp )
              !$acc atomic update
              u_max(icrm) = max( u_max(icrm) , tmp )
            enddo
          enddo
        enddo
      enddo

      cfl = 0.
      !$acc parallel loop collapse(2) private(tmp) copyin(wm,uhm,dz,adzw) copy(cfl) async(1)
      do icrm = 1 , ncrms
        do k=1,nzm
          tmp = max( uhm(k,icrm)*dt*sqrt((1./dx)**2+YES3D*(1./dy)**2) , max(wm(k,icrm),wm(k+1,icrm))*dt/(dz(icrm)*adzw(k,icrm)) )
          !$acc atomic update
          cfl = max( cfl , tmp )
        end do
      end do

      call kurant_sgs(ncrms,cfl)
      !$acc wait(1)
      ncycle = max(ncycle,max(1,ceiling(cfl/0.7)))

      if(ncycle.gt.4) then
        if(masterproc) print *,'kurant() - the number of cycles exceeded 4.'
        do icrm = 1 , ncrms
          write(0, 5550) cfl, cfl_sgs, latitude(1,1,icrm), longitude(1,1,icrm)
          do k=1, nzm
            write(0, 5551) k, wm(k,icrm), uhm(k,icrm), tabs(1,1,k,icrm)
          end do
        end do
        call task_abort()
      end if

5550 format('kurant() - cfl: ',f12.2,'  cfl_sgs: ',f12.2,'  lat: ',f6.2,'  lon: ',f6.2)
5551 format('k: ',i5,'  wm: ',f10.2,'  uhm: ',f10.2,'  tabs: ',f8.2)

   end subroutine kurant

end module kurant_mod
