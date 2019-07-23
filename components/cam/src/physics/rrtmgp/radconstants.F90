module radconstants

! This module contains constants that are specific to the radiative transfer
! code used in the RRTMG model.
!
! TODO: This all needs to be checked for RRTMGP implementation! Some of this
! might change; band mapping in the shortwave has definitely changed, and this
! module should reflect those changes.
!
! TODO: Should this data be handled in a more robust way? Much of this contains
! explicit mappings to indices, which would probably be better handled with get_
! functions. I.e., get_nswbands() could query the kdist objects in case of
! RRTMGP, and the diag indices could look up the actual bands used in the kdist
! objects as well. On that note, this module should probably go away if
! possible in the future, and we should provide more robust access to the
! radiation interface.

use shr_kind_mod,         only: r8 => shr_kind_r8
use cam_abortutils,       only: endrun
use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp

implicit none
private
save

type(ty_gas_optics_rrtmgp), public :: k_dist_sw, k_dist_lw
integer, public :: nlwgpts
integer, public :: nswgpts

! SHORTWAVE DATA

! number of shorwave spectral intervals
integer, parameter, public :: nswbands = 14

! These are indices to the band for diagnostic output
integer, parameter, public :: idx_sw_diag = 11 ! index to sw visible band
integer, parameter, public :: idx_nir_diag = 9 ! index to sw near infrared (778-1240 nm) band
integer, parameter, public :: idx_uv_diag = 12 ! index to sw uv (345-441 nm) band

integer, parameter, public :: rrtmg_sw_cloudsim_band = 10  ! rrtmg band for .67 micron

! Number of evenly spaced intervals in rh
! The globality of this mesh may not be necessary
! Perhaps it could be specific to the aerosol
! But it is difficult to see how refined it must be
! for lookup.  This value was found to be sufficient
! for Sulfate and probably necessary to resolve the
! high variation near rh = 1.  Alternative methods
! were found to be too slow.
! Optimal approach would be for cam to specify size of aerosol
! based on each aerosol's characteristics.  Radiation 
! should know nothing about hygroscopic growth!
integer, parameter, public :: nrh = 1000  

! LONGWAVE DATA

! These are indices to the band for diagnostic output
integer, parameter, public :: idx_lw_diag = 7 ! index to (H20 window) LW band

integer, parameter, public :: rrtmg_lw_cloudsim_band = 6  ! rrtmg band for 10.5 micron

! number of lw bands
integer, parameter, public :: nlwbands = 16

!These can go away when old camrt disappears
! Index of volc. abs., H2O non-window
integer, public, parameter :: idx_LW_H2O_NONWND=1
! Index of volc. abs., H2O window
integer, public, parameter :: idx_LW_H2O_WINDOW=2
! Index of volc. cnt. abs. 0500--0650 cm-1
integer, public, parameter :: idx_LW_0500_0650=3
! Index of volc. cnt. abs. 0650--0800 cm-1
integer, public, parameter :: idx_LW_0650_0800=4
! Index of volc. cnt. abs. 0800--1000 cm-1
integer, public, parameter :: idx_LW_0800_1000=5
! Index of volc. cnt. abs. 1000--1200 cm-1
integer, public, parameter :: idx_LW_1000_1200=6
! Index of volc. cnt. abs. 1200--2000 cm-1
integer, public, parameter :: idx_LW_1200_2000=7

! GASES TREATED BY RADIATION (line spectrae)

! gasses required by radiation
integer, public, parameter :: gasnamelength = 5
integer, public, parameter :: nradgas = 8
character(len=gasnamelength), public, parameter :: gaslist(nradgas) &
   = (/'H2O  ','O3   ', 'O2   ', 'CO2  ', 'N2O  ', 'CH4  ', 'CFC11', 'CFC12'/)

! what is the minimum mass mixing ratio that can be supported by radiation implementation?
real(r8), public, parameter :: minmmr(nradgas) &
   = epsilon(1._r8)

! Length of "optics type" string specified in optics files.
integer, parameter, public :: ot_length = 32

public :: rad_gas_index

public :: get_number_sw_bands, &
          get_sw_spectral_boundaries, &
          get_lw_spectral_boundaries, &
          get_sw_spectral_midpoints, &
          get_lw_spectral_midpoints

contains

!------------------------------------------------------------------------------

subroutine get_number_sw_bands(number_of_bands)

   ! number of solar (shortwave) bands in the rrtmg code
   integer, intent(out) :: number_of_bands

   number_of_bands = nswbands

end subroutine get_number_sw_bands

!------------------------------------------------------------------------------

subroutine get_lw_spectral_boundaries(low_boundaries, high_boundaries, units)

   ! Provide spectral boundaries of each longwave band

   real(r8), intent(out) :: low_boundaries(nlwbands), high_boundaries(nlwbands)
   character(*), intent(in) :: units ! requested units
   real(r8) :: wavenumber_bounds(2,nlwbands)

   wavenumber_bounds = k_dist_lw%get_band_lims_wavenumber()

   select case (units)
   case ('inv_cm','cm^-1','cm-1')
      low_boundaries  = wavenumber_bounds(1,:)
      high_boundaries = wavenumber_bounds(2,:)
   case('m','meter','meters')
      low_boundaries  = 1.e-2_r8/wavenumber_bounds(2,:)
      high_boundaries = 1.e-2_r8/wavenumber_bounds(1,:)
   case('nm','nanometer','nanometers')
      low_boundaries  = 1.e7_r8/wavenumber_bounds(2,:)
      high_boundaries = 1.e7_r8/wavenumber_bounds(1,:)
   case('um','micrometer','micrometers','micron','microns')
      low_boundaries  = 1.e4_r8/wavenumber_bounds(2,:)
      high_boundaries = 1.e4_r8/wavenumber_bounds(1,:)
   case('cm','centimeter','centimeters')
      low_boundaries  = 1._r8/wavenumber_bounds(2,:)
      high_boundaries = 1._r8/wavenumber_bounds(1,:)
   case default
      call endrun('get_lw_spectral_boundaries: spectral units not acceptable'//units)
   end select

end subroutine get_lw_spectral_boundaries

!------------------------------------------------------------------------------

subroutine get_lw_spectral_midpoints(band_midpoints, units)
   real(r8), intent(out) :: band_midpoints(nlwbands)
   character(len=*), intent(in) :: units
   real(r8) :: lower_boundaries(nlwbands)
   real(r8) :: upper_boundaries(nlwbands)
   integer :: iband

   ! Get band boundaries
   call get_lw_spectral_boundaries(lower_boundaries, upper_boundaries, units)

   ! Get band midpoints
   do iband = 1,nlwbands
      band_midpoints(iband) = 0.5 * ( &
         lower_boundaries(iband) + upper_boundaries(iband) &
      )
   end do
end subroutine get_lw_spectral_midpoints

!------------------------------------------------------------------------------

subroutine get_sw_spectral_boundaries(low_boundaries, high_boundaries, units)

   ! Provide spectral boundaries of each shortwave band

   real(r8), intent(out) :: low_boundaries(nswbands), high_boundaries(nswbands)
   character(*), intent(in) :: units ! requested units
   real(r8) :: wavenumber_bounds(2,nswbands)

   wavenumber_bounds = k_dist_sw%get_band_lims_wavenumber()

   select case (units)
   case ('inv_cm','cm^-1','cm-1')
      low_boundaries = wavenumber_bounds(1,:)
      high_boundaries = wavenumber_bounds(2,:)
   case('m','meter','meters')
      low_boundaries = 1.e-2_r8/wavenumber_bounds(2,:)
      high_boundaries = 1.e-2_r8/wavenumber_bounds(1,:)
   case('nm','nanometer','nanometers')
      low_boundaries = 1.e7_r8/wavenumber_bounds(2,:)
      high_boundaries = 1.e7_r8/wavenumber_bounds(1,:)
   case('um','micrometer','micrometers','micron','microns')
      low_boundaries = 1.e4_r8/wavenumber_bounds(2,:)
      high_boundaries = 1.e4_r8/wavenumber_bounds(1,:)
   case('cm','centimeter','centimeters')
      low_boundaries  = 1._r8/wavenumber_bounds(2,:)
      high_boundaries = 1._r8/wavenumber_bounds(1,:)
   case default
      call endrun('radconstants.F90: spectral units not acceptable'//units)
   end select

end subroutine get_sw_spectral_boundaries

!------------------------------------------------------------------------------

subroutine get_sw_spectral_midpoints(band_midpoints, units)
   real(r8), intent(out) :: band_midpoints(nswbands)
   character(len=*), intent(in) :: units
   real(r8) :: lower_boundaries(nswbands)
   real(r8) :: upper_boundaries(nswbands)
   integer :: iband

   ! Get band boundaries
   call get_sw_spectral_boundaries(lower_boundaries, upper_boundaries, units)

   ! Get band midpoints
   do iband = 1,nswbands
      band_midpoints(iband) = 0.5 * ( &
         lower_boundaries(iband) + upper_boundaries(iband) &
      )
   end do
end subroutine get_sw_spectral_midpoints

!------------------------------------------------------------------------------

integer function rad_gas_index(gasname)

   ! return the index in the gaslist array of the specified gasname

   character(len=*),intent(in) :: gasname
   integer :: igas

   rad_gas_index = -1
   do igas = 1, nradgas
      if (trim(gaslist(igas)).eq.trim(gasname)) then
         rad_gas_index = igas
         return
      endif
   enddo
   call endrun ("rad_gas_index: can not find gas with name "//gasname)
end function rad_gas_index

end module radconstants
