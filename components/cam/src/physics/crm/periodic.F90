
subroutine periodic(flag)

use vars
use microphysics
use crmtracers
implicit none

integer flag, i

if(flag.eq.0) then

  call bound_exchange(u,dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm,1,1,1,1,1)
  call bound_exchange(v,dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm,1,1,1,1,2)
   ! use w at the top level  - 0s anyway - to exchange the sst boundaries (for
   ! surface fluxes call
  w(1:nx,1:ny,nz) = sstxy(1:nx,1:ny)
  call bound_exchange(w,dimx1_w,dimx2_w,dimy1_w,dimy2_w,nz,1,1,1,1,3)
  sstxy(0:nx,1-YES3D:ny) = w(0:nx,1-YES3D:ny,nz)
  w(0:nx+1,1-YES3D:ny+YES3D,nz) = 0.

endif

if(flag.eq.1) then

  call bound_exchange(u,dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm,2,3,2,2,1)
  call bound_exchange(v,dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm,2,2,2,3,2)
  call bound_exchange(w,dimx1_w,dimx2_w,dimy1_w,dimy2_w,nz,2,2,2,2,3)

endif

if(flag.eq.2) then

 call bound_exchange(t,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,4)
 if(dosgs.and..not.dosmagor.or.doscalar) &
     call bound_exchange(tke,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,5)
 call bound_exchange(tk,0,nxp1,1-YES3D,nyp1,nzm,1,1,1,1,6)
 call bound_exchange(tkh,0,nxp1,1-YES3D,nyp1,nzm,1,1,1,1,7)
#ifdef CLUBB_CRM     
 ! Vince Larson (UWM) changed so that bound_exchange is called even if
 !     docloud = .false. and doclubb = .true.    11 Nov 2007
#endif 
 do i = 1,nmicro_fields
    if(   i.eq.index_water_vapor             &
#ifdef CLUBB_CRM
     .or. (docloud.or.doclubb.or.doclubbnoninter) .and.flag_precip(i).ne.1    &
#else
     .or. docloud.and.flag_precip(i).ne.1    &
#endif
     .or. doprecip.and.flag_precip(i).eq.1 ) &
     call bound_exchange(micro_field(:,:,:,i),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,7+i)
 end do
 if(dotracers) then
   do i=1,ntracers
     call bound_exchange(tracer(:,:,:,i),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,7+nmicro_fields+i)
   end do
 end if

endif
        
if(flag.eq.3) then
        
 call bound_exchange(t,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,4)
 if(dosgs.and..not.dosmagor.or.doscalar) &
     call bound_exchange(tke,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,5)
#ifdef CLUBB_CRM     
 ! Vince Larson (UWM) changed so that bound_exchange is called even if
 !     docloud = .false. and doclubb = .true.    11 Nov 2007
#endif 
 do i = 1,nmicro_fields
    if(   i.eq.index_water_vapor             &
#ifdef CLUBB_CRM     
     .or. (docloud.or.doclubb.or.doclubbnoninter) .and.flag_precip(i).ne.1    &
#else
     .or. docloud.and.flag_precip(i).ne.1    &
#endif
     .or. doprecip.and.flag_precip(i).eq.1 ) &
     call bound_exchange(micro_field(:,:,:,i),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,7+i)
 end do
 if(dotracers) then
   do i=1,ntracers
     call bound_exchange(tracer(:,:,:,i),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,7+nmicro_fields+i)
   end do
 end if
        
endif
        
        
end subroutine periodic

