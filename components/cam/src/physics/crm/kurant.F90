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
      real(crm_rknd) wm(nz)  ! maximum vertical wind velocity
      real(crm_rknd) uhm(nz) ! maximum horizontal wind velocity
      real(crm_rknd) cfl, cfl_sgs
      ncycle = 1

      do icrm = 1 , ncrms
        wm(nz)=0.
        do k = 1,nzm
          wm(k) = maxval(abs(w(1:nx,1:ny,k,icrm)))
          uhm(k) = sqrt(maxval(u(1:nx,1:ny,k,icrm)**2+YES3D*v(1:nx,1:ny,k,icrm)**2))
        end do
        w_max(icrm)=max( w_max(icrm), real(maxval(w(1:nx,1:ny,1:nz,icrm)),kind(w_max(icrm))) )
        u_max(icrm)=max( u_max(icrm), real(maxval(uhm(1:nzm))       ,kind(u_max(icrm))) )

        cfl = 0.
        do k=1,nzm
          cfl = max(cfl,uhm(k)*dt*sqrt((1./dx)**2+YES3D*(1./dy)**2), &
          max(wm(k),wm(k+1))*dt/(dz(icrm)*adzw(k,icrm)) )
        end do

        call kurant_sgs(ncrms,icrm,cfl_sgs)
        cfl = max(cfl,cfl_sgs)

        ncycle = max(ncycle,max(1,ceiling(cfl/0.7)))

        if(ncycle.gt.4) then
          if(masterproc) print *,'kurant() - the number of cycles exceeded 4.'
          write(0, 5550) cfl, cfl_sgs, latitude(1,1,icrm), longitude(1,1,icrm)
          do k=1, nzm
            write(0, 5551) k, wm(k), uhm(k), tabs(1,1,k,icrm)
          end do

          call task_abort()
        end if
      enddo

5550 format('kurant() - cfl: ',f12.2,'  cfl_sgs: ',f12.2,'  lat: ',f6.2,'  lon: ',f6.2)
5551 format('k: ',i5,'  wm: ',f10.2,'  uhm: ',f10.2,'  tabs: ',f8.2)

   end subroutine kurant

end module kurant_mod
