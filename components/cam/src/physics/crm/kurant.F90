
subroutine kurant

use vars
use sgs, only: kurant_sgs
use params, only: crm_rknd

implicit none

integer i, j, k, ncycle1(1),ncycle2(1)
real(crm_rknd) wm(nz)  ! maximum vertical wind velocity
real(crm_rknd) uhm(nz) ! maximum horizontal wind velocity
real(crm_rknd) cfl, cfl_sgs

ncycle = 1
	
wm(nz)=0.
do k = 1,nzm
 wm(k) = maxval(abs(w(1:nx,1:ny,k)))
 uhm(k) = sqrt(maxval(u(1:nx,1:ny,k)**2+YES3D*v(1:nx,1:ny,k)**2))
end do
w_max=max( w_max, real(maxval(w(1:nx,1:ny,1:nz)),kind(w_max)) )
u_max=max( u_max, real(maxval(uhm(1:nzm))       ,kind(u_max)) )

cfl = 0.
do k=1,nzm
  cfl = max(cfl,uhm(k)*dt*sqrt((1./dx)**2+YES3D*(1./dy)**2), &
                   max(wm(k),wm(k+1))*dt/(dz*adzw(k)) )
end do

call kurant_sgs(cfl_sgs)
cfl = max(cfl,cfl_sgs)
	
ncycle = max(1,ceiling(cfl/0.7))

if(dompi) then
  ncycle1(1)=ncycle
  call task_max_integer(ncycle1,ncycle2,1)
  ncycle=ncycle2(1)
end if
if(ncycle.gt.4) then
   if(masterproc) print *,'the number of cycles exceeded 4.'
!+++ test +++mhwang
   write(0, *) 'cfl', cfl, cfl_sgs, latitude(1, 1), longitude(1,1)
   do k=1, nzm
      write(0, *) 'k=', k, wm(k), uhm(k)
   end do
   do i=1, nx
     write(0, *) 'i=', i,  u(i, 1, 4), v(i, 1, 4), tabs(i,1,4)
   end do
!---mhwang
   call task_abort()
end if

end subroutine kurant	
