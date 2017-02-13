module params

use grid, only: nzm
#ifdef CLUBB_CRM
! Use the CLUBB values for these constants for consistency
use constants_clubb, only: Cp, ggr => grav, lcond => Lv, lfus => Lf, &
  lsub => Ls, Rv, rgas => Rd 

#else

#ifdef CRM
use shr_const_mod, only: shr_const_rdair, shr_const_cpdair, shr_const_latvap, &
                           shr_const_latice, shr_const_latsub, shr_const_rgas, &
                           shr_const_mwwv, shr_const_stebol, shr_const_tkfrz, &
                           shr_const_mwdair, shr_const_g, shr_const_karman, &
                           shr_const_rhofw
#endif /*CRM*/

#endif

implicit none

!   Constants:

#ifdef CLUBB_CRM
! Define Cp, ggr, etc. in module constants_clubb
#else
#ifndef CRM
real, parameter :: cp = 1004.             ! Specific heat of air, J/kg/K
real, parameter :: ggr = 9.81             ! Gravity acceleration, m/s2
real, parameter :: lcond = 2.5104e+06     ! Latent heat of condensation, J/kg
real, parameter :: lfus = 0.3336e+06      ! Latent heat of fusion, J/kg
real, parameter :: lsub = 2.8440e+06      ! Latent heat of sublimation, J/kg
real, parameter :: rv = 461.              ! Gas constant for water vapor, J/kg/K
real, parameter :: rgas = 287.            ! Gas constant for dry air, J/kg/K
#else
real, parameter :: cp = shr_const_cpdair
real, parameter :: ggr = shr_const_g
real, parameter :: lcond = shr_const_latvap 
real, parameter :: lfus = shr_const_latice 
real, parameter :: lsub = lcond + lfus
real, parameter :: rv = shr_const_rgas/shr_const_mwwv 
real, parameter :: rgas = shr_const_rdair
#endif
#endif
real, parameter :: diffelq = 2.21e-05     ! Diffusivity of water vapor, m2/s
real, parameter :: therco = 2.40e-02      ! Thermal conductivity of air, J/m/s/K
real, parameter :: muelq = 1.717e-05      ! Dynamic viscosity of air

real, parameter :: fac_cond = lcond/cp 
real, parameter :: fac_fus = lfus/cp
real, parameter :: fac_sub = lsub/cp


!  Variables:

            
real  pres0      ! Reference surface pressure, Pa
real  ug	 ! Velocity of the Domain's drift in x direction
real  vg	 ! Velocity of the Domain's drift in y direction
real  fcor       ! Coriolis parameter	
real  fcorz      ! Vertical Coriolis parameter
real  pi         ! 3.1415...

real longitude0  ! latitude of the domain's center 
real latitude0   ! longitude of the domain's center 
real coszrs      ! solar zenith angle


!  Surface stuff:   	

real   tabs_s	! surface temperature,K
real   fluxt0   ! surface sensible flux, Km/s
real   fluxq0   ! surface latent flux, m/s
real   tau0     ! surface stress, m2/s2
real   z0	! roughness length
real   soil_wetness ! wetness coeff for soil (from 0 to 1.)
real   epsv     ! = (1-eps)/eps, where eps= Rv/Ra, or =0. if dosmoke=.true.
integer ocean_type ! type of SST forcing
 
!  Large-scale stuff

real  timelargescale ! time to start large-scale forcing
real  tauls	     ! nudging-to-large-scaler-profile time-scale

real uhl        ! current large-scale velocity in x near sfc
real vhl        ! current large-scale velocity in y near sfc
real   taux0    ! surface stress in x, m2/s2
real   tauy0    ! surface stress in y, m2/s2


end module params
