module microphysics

! module for original SAM bulk microphysics
! Marat Khairoutdinov, 2006

use grid, only: nx,ny,nzm,nz, dimx1_s,dimx2_s,dimy1_s,dimy2_s ! subdomain grid information 
use micro_params
implicit none

!----------------------------------------------------------------------
!!! required definitions:

integer, parameter :: nmicro_fields = 2   ! total number of prognostic water vars

!!! microphysics prognostic variables are storred in this array:

real micro_field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, nmicro_fields)

integer, parameter :: flag_wmass(nmicro_fields) = (/1,1/)
integer, parameter :: index_water_vapor = 1 ! index for variable that has water vapor
integer, parameter :: index_cloud_ice = 1   ! index for cloud ice (sedimentation)
integer, parameter :: flag_precip(nmicro_fields) = (/0,1/)

! both variables correspond to mass, not number
integer, parameter :: flag_number(nmicro_fields) = (/0,0/)

! SAM1MOM 3D microphysical fields are output by default.
integer, parameter :: flag_micro3Dout(nmicro_fields) = (/0,0/)

real fluxbmk (nx, ny, 1:nmicro_fields) ! surface flux of tracers
real fluxtmk (nx, ny, 1:nmicro_fields) ! top boundary flux of tracers

!!! these arrays are needed for output statistics:

real mkwle(nz,1:nmicro_fields)  ! resolved vertical flux
real mkwsb(nz,1:nmicro_fields)  ! SGS vertical flux
real mkadv(nz,1:nmicro_fields)  ! tendency due to vertical advection
real mklsadv(nz,1:nmicro_fields)  ! tendency due to large-scale vertical advection
real mkdiff(nz,1:nmicro_fields)  ! tendency due to vertical diffusion

!======================================================================
! UW ADDITIONS

!bloss: arrays with names/units for microphysical outputs in statistics.
character*3, dimension(nmicro_fields) :: mkname
character*80, dimension(nmicro_fields) :: mklongname
character*10, dimension(nmicro_fields) :: mkunits
real, dimension(nmicro_fields) :: mkoutputscale

! END UW ADDITIONS
!======================================================================

!------------------------------------------------------------------
! Optional (internal) definitions)

! make aliases for prognostic variables:
! note that the aliases should be local to microphysics

real q(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)   ! total nonprecipitating water
real qp(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)  ! total precipitating water
equivalence (q(dimx1_s,dimy1_s,1),micro_field(dimx1_s,dimy1_s,1,1))
equivalence (qp(dimx1_s,dimy1_s,1),micro_field(dimx1_s,dimy1_s,1,2))

real qn(nx,ny,nzm)  ! cloud condensate (liquid + ice)

real qpsrc(nz)  ! source of precipitation microphysical processes
real qpevp(nz)  ! sink of precipitating water due to evaporation

real vrain, vsnow, vgrau, crain, csnow, cgrau  ! precomputed coefs for precip terminal velocity

CONTAINS

! required microphysics subroutines and function:
!----------------------------------------------------------------------
!!! Read microphysics options from prm file

subroutine micro_setparm()
  ! no user-definable options in SAM1MOM microphysics.
end subroutine micro_setparm

!----------------------------------------------------------------------
!!! Initialize microphysics:


subroutine micro_init()

#ifdef CLUBB_CRM  
  use vars, only: q0, docloud, doprecip, nrestart, dosmoke, &
    doclubb, doclubbnoninter ! dschanen UWM 21 May 2008
  use grid, only: nclubb
#else
  use vars, only: q0, docloud, doprecip, nrestart, dosmoke
#endif
  integer k
#ifdef CLUBB_CRM  
  if ( nclubb /= 1 ) then
    write(0,*) "The namelist parameter nclubb is not equal to 1,",  &
      " but SAM single moment microphysics is enabled."
    write(0,*) "This will create unrealistic results in subsaturated grid boxes. ", &
      "Exiting..."
    call task_abort()
  end if
#endif

  a_bg = 1./(tbgmax-tbgmin)
  a_pr = 1./(tprmax-tprmin)
  a_gr = 1./(tgrmax-tgrmin)

  if(nrestart.eq.0) then

#ifndef CRM
     micro_field = 0.
     do k=1,nzm
      q(:,:,k) = q0(k)
     end do
     qn = 0.
#endif

     fluxbmk = 0.
     fluxtmk = 0.

#ifdef CLUBB_CRM
     if ( docloud .or. doclubb ) then
#else
     if(docloud) then
#endif
#ifndef CRM
       call cloud()
#endif
       call micro_diagnose()
     end if
     if(dosmoke) then
       call micro_diagnose()
     end if

  end if

  mkwle = 0.
  mkwsb = 0.
  mkadv = 0.
  mkdiff = 0.

  qpsrc = 0.
  qpevp = 0.

  mkname(1) = 'QT'
  mklongname(1) = 'TOTAL WATER (VAPOR + CONDENSATE)'
  mkunits(1) = 'g/kg'
  mkoutputscale(1) = 1.e3

  mkname(2) = 'QP'
  mklongname(2) = 'PRECIPITATING WATER'
  mkunits(2) = 'g/kg'
  mkoutputscale(2) = 1.e3

end subroutine micro_init

!----------------------------------------------------------------------
!!! fill-in surface and top boundary fluxes:
!
subroutine micro_flux()

  use vars, only: fluxbq, fluxtq

#ifdef CLUBB_CRM
  ! Added by dschanen UWM
  use grid, only: doclubb, doclubb_sfc_fluxes, docam_sfc_fluxes
  if ( .not. (doclubb_sfc_fluxes .or. docam_sfc_fluxes) ) then
    fluxbmk(:,:,index_water_vapor) = fluxbq(:,:)
    fluxtmk(:,:,index_water_vapor) = fluxtq(:,:)
  else
    ! Add this in later
    fluxbmk(:,:,index_water_vapor) = 0.0
    fluxtmk(:,:,index_water_vapor) = 0.0
  end if
#else
  fluxbmk(:,:,index_water_vapor) = fluxbq(:,:)
  fluxtmk(:,:,index_water_vapor) = fluxtq(:,:)
#endif /*CLUBB_CRM*/

end subroutine micro_flux

!----------------------------------------------------------------------
!!! compute local microphysics processes (bayond advection and SGS diffusion):
!
subroutine micro_proc()

#ifdef CLUBB_CRM
   use vars, only: nstep,dt,icycle,docloud,doprecip,dosmoke &
     , doclubb, doclubbnoninter ! dschanen UWM 21 May 2008
#else
   use vars, only: nstep,dt,icycle,docloud,doprecip,dosmoke
#endif /*CLUBB_CRM*/

   ! Update bulk coefficient
   if(doprecip.and.icycle.eq.1) call precip_init() 

   if(docloud) then
     call cloud()
     if(doprecip) call precip_proc()
     call micro_diagnose()
   end if
   if(dosmoke) then
     call micro_diagnose()
   end if
#ifdef CLUBB_CRM
   if ( doclubb ) then ! -dschanen UWM 21 May 2008
     if(doprecip) call precip_proc()
     call micro_diagnose()
   end if
#endif /*CLUBB_CRM*/

end subroutine micro_proc

!----------------------------------------------------------------------
!!! Diagnose arrays nessesary for dynamical core and statistics:
!
subroutine micro_diagnose()
 
   use vars

   real omn, omp
   integer i,j,k

   do k=1,nzm
    do j=1,ny
     do i=1,nx
       qv(i,j,k) = q(i,j,k) - qn(i,j,k)
       omn = max(0.,min(1.,(tabs(i,j,k)-tbgmin)*a_bg))
       qcl(i,j,k) = qn(i,j,k)*omn
       qci(i,j,k) = qn(i,j,k)*(1.-omn)
       omp = max(0.,min(1.,(tabs(i,j,k)-tprmin)*a_pr))
       qpl(i,j,k) = qp(i,j,k)*omp
       qpi(i,j,k) = qp(i,j,k)*(1.-omp)
     end do
    end do
   end do
       


end subroutine micro_diagnose

#ifdef CLUBB_CRM
!---------------------------------------------------------------------
subroutine micro_update()

! Description:
! This subroutine essentially does what micro_proc does but does not
! call any microphysics subroutines.  We need this so that CLUBB gets a
! properly updated value of ice fed in. 
!
! dschanen UWM 7 Jul 2008
!---------------------------------------------------------------------

   call cloud()
   call micro_diagnose()

end subroutine micro_update

!---------------------------------------------------------------------
subroutine micro_adjust( new_qv, new_qc )
! Description:
! Adjust vapor and liquid water.
! Microphysical variables are stored separately in 
!    SAM's dynamics + CLUBB ( e.g. qv, qcl, qci) and
!    SAM's microphysics. (e.g. q and qn).
! This subroutine stores values of qv, qcl updated by CLUBB
!   in the single-moment microphysical variables q and qn.
!
! dschanen UWM 20 May 2008
!---------------------------------------------------------------------

  use vars, only: qci

  implicit none

  real, dimension(nx,ny,nzm), intent(in) :: &
  new_qv, & ! Water vapor mixing ratio that has been adjusted by CLUBB [kg/kg]
  new_qc    ! Cloud water mixing ratio that has been adjusted by CLUBB [kg/kg]

  q(1:nx,1:ny,1:nzm) = new_qv + new_qc + qci ! Vapor + Liquid + Ice
  qn(1:nx,1:ny,1:nzm) = new_qc + qci ! Liquid + Ice

  return
end subroutine micro_adjust
#endif /*CLUBB_CRM*/

!----------------------------------------------------------------------
!!! function to compute terminal velocity for precipitating variables:
! In this particular case there is only one precipitating variable.

real function term_vel_qp(i,j,k,ind)
  
  use vars
  integer, intent(in) :: i,j,k,ind
  real wmax, omp, omg, qrr, qss, qgg

  term_vel_qp = 0.
  if(qp(i,j,k).gt.qp_threshold) then
    omp = max(0.,min(1.,(tabs(i,j,k)-tprmin)*a_pr))
    if(omp.eq.1.) then
       term_vel_qp = vrain*(rho(k)*qp(i,j,k))**crain
    elseif(omp.eq.0.) then
       omg = max(0.,min(1.,(tabs(i,j,k)-tgrmin)*a_gr))
       qgg=omg*qp(i,j,k)
       qss=qp(i,j,k)-qgg
       term_vel_qp = (omg*vgrau*(rho(k)*qgg)**cgrau &
                                 +(1.-omg)*vsnow*(rho(k)*qss)**csnow)
    else
       omg = max(0.,min(1.,(tabs(i,j,k)-tgrmin)*a_gr))
       qrr=omp*qp(i,j,k)
       qss=qp(i,j,k)-qrr
       qgg=omg*qss
       qss=qss-qgg
       term_vel_qp = (omp*vrain*(rho(k)*qrr)**crain &
                     +(1.-omp)*(omg*vgrau*(rho(k)*qgg)**cgrau &
                          +(1.-omg)*vsnow*(rho(k)*qss)**csnow))
    endif
  end if  
end function term_vel_qp

!----------------------------------------------------------------------
!!! compute sedimentation 
!
subroutine micro_precip_fall()
  
  use vars
  use params, only : pi

  real omega(nx,ny,nzm)
  integer ind
  integer i,j,k

  crain = b_rain / 4.
  csnow = b_snow / 4.
  cgrau = b_grau / 4.
  vrain = a_rain * gamr3 / 6. / (pi * rhor * nzeror) ** crain
  vsnow = a_snow * gams3 / 6. / (pi * rhos * nzeros) ** csnow
  vgrau = a_grau * gamg3 / 6. / (pi * rhog * nzerog) ** cgrau

 do k=1,nzm
  do j=1,ny
   do i=1,nx
       omega(i,j,k) = max(0.,min(1.,(tabs(i,j,k)-tprmin)*a_pr))
   end do
  end do
 end do

 call precip_fall(qp, term_vel_qp, 2, omega, ind)

 do j=1,ny
   do i=1,nx
     if(qp(i,j,1).gt.1.e-6) s_ar=s_ar+dtfactor
   end do
 end do


end subroutine micro_precip_fall

!----------------------------------------------------------------------
! called when stepout() called

subroutine micro_print()
end subroutine micro_print

!-----------------------------------------------------------------------
! Supply function that computes total water in a domain:
!
real(8) function total_water()

  use vars, only : nstep,nprint,adz,dz,rho
  real(8) tmp
  integer i,j,k,m

  total_water = 0.
  do m=1,nmicro_fields
   if(flag_wmass(m).eq.1) then
    do k=1,nzm
      tmp = 0.
      do j=1,ny
        do i=1,nx
          tmp = tmp + micro_field(i,j,k,m)
        end do
      end do
      total_water = total_water + tmp*adz(k)*dz*rho(k)
    end do
   end if
  end do

end function total_water


end module microphysics



