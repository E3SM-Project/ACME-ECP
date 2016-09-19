subroutine stepout(nstatsteps)

use vars
!use rad, only: qrad
use crmtracers
use microphysics, only: micro_print
use params
implicit none	
	
integer i,j,k,ic,jc,nstatsteps
real div, divmax, divmin
real rdx, rdy, rdz, coef
integer im,jm,km
real wmax, qnmax(1), qnmax1(1)
real(8) xbuffer(5), xbuffer1(5)



!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Print stuff out:

call t_startf ('print_out')

if(masterproc) print *,'NSTEP = ',nstep,'    NCYCLE=',ncycle

if(mod(nstep,nprint).eq.0) then
	

 divmin=1.e20
 divmax=-1.e20
	 
 rdx = 1./dx
 rdy = 1./dy

 wmax=0.
 do k=1,nzm
  coef = rho(k)*adz(k)*dz
  rdz = 1./coef
  if(ny.ne.1) then
   do j=1,ny-1*YES3D
    jc = j+1*YES3D
    do i=1,nx-1
     ic = i+1
     div = (u(ic,j,k)-u(i,j,k))*rdx + (v(i,jc,k)-v(i,j,k))*rdy + &
		  (w(i,j,k+1)*rhow(k+1)-w(i,j,k)*rhow(k))*rdz
     divmax = max(divmax,div)
     divmin = min(divmin,div)
     if(w(i,j,k).gt.wmax) then
	wmax=w(i,j,k)
	im=i
	jm=j
	km=k
     endif
    end do
   end do
  else
    j = 1
    do i=1,nx-1
    ic = i+1
     div = (u(ic,j,k)-u(i,j,k))*rdx +(w(i,j,k+1)*rhow(k+1)-w(i,j,k)*rhow(k))*rdz
     divmax = max(divmax,div)
     divmin = min(divmin,div)
     if(w(i,j,k).gt.wmax) then
	wmax=w(i,j,k)
	im=i
	jm=j
	km=k
     endif
    end do
  endif
 end do

 if(dompi) then
   xbuffer(1) = total_water_before
   xbuffer(2) = total_water_after
   xbuffer(3) = total_water_evap
   xbuffer(4) = total_water_prec
   xbuffer(5) = total_water_ls
   call task_sum_real8(xbuffer, xbuffer1,5)
   total_water_before = xbuffer1(1)
   total_water_after = xbuffer1(2)
   total_water_evap = xbuffer1(3)
   total_water_prec = xbuffer1(4)
   total_water_ls = xbuffer1(5)
 end if

!print*,rank,minval(u(1:nx,1:ny,:)),maxval(u(1:nx,1:ny,:))
!print*,rank,'min:',minloc(u(1:nx,1:ny,:))
!print*,rank,'max:',maxloc(u(1:nx,1:ny,:))

!if(rank.eq.2) then

!print*,'p:'
!write(6,'(16f7.2)')((p(i,13,k),i=1,16),k=30,1,-1)
!print*,'u:'
!write(6,'(16f7.2)')((u(i,13,k),i=1,16),k=30,1,-1)
!print*,'v:'
!write(6,'(16f7.2)')((v(i,13,k),i=1,16),k=30,1,-1)
!print*,'w:'
!write(6,'(16f7.2)')((w(i,13,k),i=1,16),k=30,1,-1)
!print*,'qcl:'
!write(6,'(16f7.2)')((qcl(i,13,k)*1000.,i=1,16),k=30,1,-1)
!print*,'qpl:'
!write(6,'(16f7.2)')((qpl(i,13,k)*1000.,i=1,16),k=30,1,-1)
!print*,'qrad:'
!write(6,'(16f7.2)')((qrad(i,13,k)*3600.,i=1,16),k=30,1,-1)
!print*,'qv:'
!write(6,'(16f7.2)')((qv(i,13,k)*1000.,i=1,16),k=30,1,-1)
!print*,'tabs:'
!write(6,'(16f7.2)')((tabs(i,13,k),i=1,16),k=30,1,-1)
!
!end if

!--------------------------------------------------------
 if(masterproc) then
	
    print*,'DAY = ',day	
    write(6,*) 'NSTEP=',nstep
    write(6,*) 'div:',divmax,divmin
    write(6,*) 'SST=',tabs_s, '  surface pressure=',pres0

 endif

 call fminmax_print('u:',u,dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm)
 call fminmax_print('v:',v,dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm-5)
 call fminmax_print('w:',w,dimx1_w,dimx2_w,dimy1_w,dimy2_w,nz)
 call fminmax_print('p:',p,0,nx,1-YES3D,ny,nzm)
 call fminmax_print('t:',t,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
 call fminmax_print('tabs:',tabs,1,nx,1,ny,nzm)
 call fminmax_print('qv:',qv,1,nx,1,ny,nzm)
 call fminmax_print('tke:',tke,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
 call fminmax_print('tk:',tk,0,nxp1,1-YES3D,nyp1,nzm)
 call fminmax_print('tkh:',tkh,0,nxp1,1-YES3D,nyp1,nzm)
#ifdef CLUBB_CRM
 if(docloud.or.doclubb) then
#else
 if(docloud) then
#endif /*CLUBB_CRM*/
   call fminmax_print('qcl:',qcl,1,nx,1,ny,nzm)
   call fminmax_print('qci:',qci,1,nx,1,ny,nzm)
   call micro_print()
 end if
 if(doprecip) then
   call fminmax_print('qpl:',qpl,1,nx,1,ny,nzm)
   call fminmax_print('qpi:',qpi,1,nx,1,ny,nzm)
 end if
! if(dolongwave.or.doshortwave) call fminmax_print('qrad(K/day):',qrad*86400.,1,nx,1,ny,nzm)
 if(dotracers) then
   do k=1,ntracers
     call fminmax_print(trim(tracername(k))//':',tracer(:,:,:,k),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
   end do
 end if
 call fminmax_print('shf:',fluxbt*cp*rho(1),1,nx,1,ny,1)
 call fminmax_print('lhf:',fluxbq*lcond*rho(1),1,nx,1,ny,1)
 call fminmax_print('uw:',fluxbu,1,nx,1,ny,1)
 call fminmax_print('vw:',fluxbv,1,nx,1,ny,1)
 call fminmax_print('sst:',sstxy,0,nx,1-YES3D,ny,1)

end if ! (mod(nstep,nprint).eq.0)

call t_stopf ('print_out')

end
