
#include "directives.inc"

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

      !_dir _enter_data _dcreate(wm,uhm) _async(1)

      !_dir _par _loop _gang _vector collapse(2) _kout(wm,uhm) _async(1)
      do icrm = 1 , ncrms
        do k = 1 , nz
          wm (k,icrm) = 0.
          uhm(k,icrm) = 0.
        enddo
      enddo

      !_dir _par _loop _gang _vector collapse(4) private(tmp) _kinout(wm,w_max,uhm,u_max) _kin(u,v,w) _async(1)
      do icrm = 1 , ncrms
        do k = 1,nzm
          do j = 1 , ny
            do i = 1 , nx
              tmp = abs(w(i,j,k,icrm))
              !_dir atomic update
              wm(k,icrm) = max( wm(k,icrm) , tmp )
              !_dir atomic update
              w_max(icrm) = max( w_max(icrm) , tmp )

              tmp = sqrt(u(i,j,k,icrm)**2+YES3D*v(i,j,k,icrm)**2)
              !_dir atomic update
              uhm(k,icrm) = max( uhm(k,icrm) , tmp )
              !_dir atomic update
              u_max(icrm) = max( u_max(icrm) , tmp )
            enddo
          enddo
        enddo
      enddo

      cfl = 0.
      !_dir _par _loop _gang _vector collapse(2) private(tmp) _kin(uhm,wm,dz,adzw) _kinout(cfl) _async(1)
      do icrm = 1 , ncrms
        do k=1,nzm
          tmp = max( uhm(k,icrm)*dt*sqrt((1./dx)**2+YES3D*(1./dy)**2) , max(wm(k,icrm),wm(k+1,icrm))*dt/(dz(icrm)*adzw(k,icrm)) )
          !_dir atomic update
          cfl = max( cfl , tmp )
        end do
      end do

      call kurant_sgs(ncrms,cfl)
      !_dir _wait(1)
      ncycle = max(ncycle,max(1,ceiling(cfl/0.7)))

      if(ncycle.gt.4) then
        !$_dir _update_host(wm,uhm,tabs,latitude,longitude)
        if(masterproc) print *,'kurant() - the number of cycles exceeded 4.'
        do icrm = 1 , ncrms
          write(0, 5550) cfl, cfl_sgs, latitude(1,1,icrm), longitude(1,1,icrm)
          do k=1, nzm
            write(0, 5551) k, wm(k,icrm), uhm(k,icrm), tabs(1,1,k,icrm)
          end do
        end do
        call task_abort()
      end if

      !_dir _exit_data _ddelete(wm,uhm) _async(1)

5550 format('kurant() - cfl: ',f12.2,'  cfl_sgs: ',f12.2,'  lat: ',f6.2,'  lon: ',f6.2)
5551 format('k: ',i5,'  wm: ',f10.2,'  uhm: ',f10.2,'  tabs: ',f8.2)

   end subroutine kurant

end module kurant_mod
