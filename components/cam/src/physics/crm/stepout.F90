module stepout_mod
  use utils, only: fminmax_print
  implicit none

contains

  subroutine stepout(nstatsteps,ncrms,icrm)

    use vars
    !use rad, only: qrad
    use sgs, only: tk, sgs_print
    use crmtracers
    use microphysics, only: micro_print
    use params
    implicit none
    integer, intent(in) :: ncrms,icrm

    integer i,j,k,ic,jc,nstatsteps
    integer n
    real(crm_rknd) div, divmax, divmin
    real(crm_rknd) rdx, rdy, rdz, coef
    integer im,jm,km
    real(crm_rknd) wmax, qnmax(1), qnmax1(1)
    real(8) buffer(6), buffer1(6)
    real(8) qi0(nzm)

#ifdef CLUBB_CRM
    real(8) buffer_e(7), buffer1_e(7)
#endif



    !------------------------------------------------------------------------------
    !------------------------------------------------------------------------------
    ! Print stuff out:

    !call t_startf ('print_out')

    if(masterproc) print *,'NSTEP = ',nstep(icrm),'    NCYCLE=',ncycle(icrm)

    if(mod(nstep(icrm),nprint).eq.0) then


      divmin=1.e20
      divmax=-1.e20

      rdx = 1./dx
      rdy = 1./dy

      wmax=0.
      do k=1,nzm
        coef = rho(icrm,k)*adz(icrm,k)*dz(icrm)
        rdz = 1./coef
        if(ny.ne.1) then
          do j=1,ny-1*YES3D
            jc = j+1*YES3D
            do i=1,nx-1
              ic = i+1
              div = (u(icrm,ic,j,k)-u(icrm,i,j,k))*rdx + (v(icrm,i,jc,k)-v(icrm,i,j,k))*rdy + &
              (w(icrm,i,j,k+1)*rhow(icrm,k+1)-w(icrm,i,j,k)*rhow(icrm,k))*rdz
              divmax = max(divmax,div)
              divmin = min(divmin,div)
              if(w(icrm,i,j,k).gt.wmax) then
                wmax=w(icrm,i,j,k)
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
            div = (u(icrm,ic,j,k)-u(icrm,i,j,k))*rdx +(w(icrm,i,j,k+1)*rhow(icrm,k+1)-w(icrm,i,j,k)*rhow(icrm,k))*rdz
            divmax = max(divmax,div)
            divmin = min(divmin,div)
            if(w(icrm,i,j,k).gt.wmax) then
              wmax=w(icrm,i,j,k)
              im=i
              jm=j
              km=k
            endif
          end do
        endif
      end do

      if(dompi) then
        buffer(1) = total_water_before(icrm)
        buffer(2) = total_water_after (icrm)
        buffer(3) = total_water_evap  (icrm)
        buffer(4) = total_water_prec  (icrm)
        buffer(5) = total_water_ls    (icrm)
#ifdef CLUBB_CRM
        buffer(6) = total_water_clubb (icrm)

        buffer_e(1) = total_energy_before(icrm)
        buffer_e(2) = total_energy_after(icrm)
        buffer_e(3) = total_energy_evap(icrm)
        buffer_e(4) = total_energy_prec(icrm)
        buffer_e(5) = total_energy_ls(icrm)
        buffer_e(6) = total_energy_clubb(icrm)
        buffer_e(7) = total_energy_rad(icrm)
#endif
        call task_sum_real8(buffer, buffer1,6)
        total_water_before(icrm) = buffer1(1)
        total_water_after(icrm) = buffer1(2)
        total_water_evap(icrm) = buffer1(3)
        total_water_prec(icrm) = buffer1(4)
        total_water_ls(icrm) = buffer1(5)
#ifdef CLUBB_CRM
        total_water_clubb(icrm) = buffer1(6)

        call task_sum_real8(buffer_e, buffer1_e,7)
        total_energy_before(icrm) = buffer1_e(1)
        total_energy_after(icrm) = buffer1_e(2)
        total_energy_evap(icrm) = buffer1_e(3)
        total_energy_prec(icrm) = buffer1_e(4)
        total_energy_ls(icrm) = buffer1_e(5)
        total_energy_clubb(icrm) = buffer1_e(6)
        total_energy_rad(icrm) = buffer1_e(7)
#endif
      end if

      !print*,rank,minval(u(icrm,1:nx,1:ny,:)),maxval(u(icrm,1:nx,1:ny,:))
      !print*,rank,'min:',minloc(u(icrm,1:nx,1:ny,:))
      !print*,rank,'max:',maxloc(u(icrm,1:nx,1:ny,:))

      !if(masterproc) then

      !print*,'--->',tk(icrm,27,1,1)
      !print*,'tk->:'
      !write(6,'(16f7.2)')((tk(icrm,i,1,k),i=1,16),k=nzm,1,-1)
      !print*,'p->:'
      !write(6,'(16f7.2)')((p(icrm,i,1,k),i=1,16),k=nzm,1,-1)
      !print*,'u->:'
      !write(6,'(16f7.2)')((u(icrm,i,1,k),i=1,16),k=nzm,1,-1)
      !print*,'v->:'
      !write(6,'(16f7.2)')((v(icrm,i,1,k),i=1,16),k=nzm,1,-1)
      !print*,'w->:'
      !write(6,'(16f7.2)')((w(icrm,i,1,k),i=1,16),k=nzm,1,-1)
      !print*,'qcl:'
      !write(6,'(16f7.2)')((qcl(icrm,i,13,k)*1000.,i=1,16),k=30,1,-1)
      !print*,'qpl:'
      !write(6,'(16f7.2)')((qpl(icrm,i,13,k)*1000.,i=1,16),k=30,1,-1)
      !print*,'qrad:'
      !write(6,'(16f7.2)')((qrad(i,13,k)*3600.,i=1,16),k=30,1,-1)
      !print*,'qv:'
      !write(6,'(16f7.2)')((qv(icrm,i,13,k)*1000.,i=1,16),k=30,1,-1)
      !print*,'tabs:'
      !write(6,'(16f7.2)')((tabs(icrm,i,13,k),i=1,16),k=30,1,-1)
      !
      !end if

      !--------------------------------------------------------
      if(masterproc) then

        print*,'DAY = ',day(icrm)
        write(6,*) 'NSTEP=',nstep(icrm)
        write(6,*) 'div:',divmax,divmin
        !if(.not.dodynamicocean) write(6,*) 'SST=',tabs_s
        ! write(6,*) 'surface pressure=',pres0

      endif

      call fminmax_print('u:',u(icrm,:,:,:),dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm)
      call fminmax_print('v:',v(icrm,:,:,:),dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm-5)
      call fminmax_print('w:',w(icrm,:,:,:),dimx1_w,dimx2_w,dimy1_w,dimy2_w,nz)
      call fminmax_print('p:',p(icrm,:,:,:),0,nx,1-YES3D,ny,nzm)
      call fminmax_print('t:',t(icrm,:,:,:),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
      call fminmax_print('tabs:',tabs(icrm,:,:,:),1,nx,1,ny,nzm)
      call fminmax_print('qv:',qv(icrm,:,:,:),1,nx,1,ny,nzm)
      if(dosgs) call sgs_print()
#ifdef CLUBB_CRM
      if(docloud.or.doclubb) then
#else
        if(docloud) then
#endif /*CLUBB_CRM*/
          call fminmax_print('qcl:',qcl(icrm,:,:,:),1,nx,1,ny,nzm)
          call fminmax_print('qci:',qci(icrm,:,:,:),1,nx,1,ny,nzm)
          call micro_print()
        end if
        if(doprecip) then
          call fminmax_print('qpl:',qpl(icrm,:,:,:),1,nx,1,ny,nzm)
          call fminmax_print('qpi:',qpi(icrm,:,:,:),1,nx,1,ny,nzm)
        end if
        ! if(dolongwave.or.doshortwave) call fminmax_print('qrad(K/day):',qrad*86400.,1,nx,1,ny,nzm)
        if(dotracers) then
          do k=1,ntracers
            call fminmax_print(trim(tracername(k))//':',tracer(icrm,:,:,:,k),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
          end do
        end if
        call fminmax_print('shf:',fluxbt(icrm,:,:)*cp*rhow(icrm,1),1,nx,1,ny,1)
        call fminmax_print('lhf:',fluxbq(icrm,:,:)*lcond*rhow(icrm,1),1,nx,1,ny,1)
        call fminmax_print('uw:',fluxbu(icrm,:,:),1,nx,1,ny,1)
        call fminmax_print('vw:',fluxbv(icrm,:,:),1,nx,1,ny,1)
        call fminmax_print('sst:',sstxy(icrm,:,:),0,nx,1-YES3D,ny,1)

      end if ! (mod(nstep(icrm),nprint).eq.0)

      !call t_stopf ('print_out')

    end

  end module stepout_mod
