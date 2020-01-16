!---------------------------------------------------------------------------
! $Id$
!===============================================================================
module array_index

  ! Description:
  ! Contains indices to variables in larger arrays.
  ! Note that the 'ii' is necessary because 'i' is used in
  ! statistics to track locations in the zt/zm/sfc derived types.

  ! References:
  !   None
  !-------------------------------------------------------------------------

  use clubb_precision, only: &
    core_rknd      ! Precision

  implicit none

  ! Variables
  ! Microphysics mixing ratios
  integer, public :: &
    iirr,   & ! Hydrometeor array index for rain water mixing ratio, rr
    iirs,   & ! Hydrometeor array index for snow mixing ratio, rs
    iiri,   & ! Hydrometeor array index for ice mixing ratio, ri
    iirg      ! Hydrometeor array index for graupel mixing ratio, rg

  ! Microphysics concentrations
  integer, public :: &
    iiNr,   & ! Hydrometeor array index for rain drop concentration, Nr
    iiNs,   & ! Hydrometeor array index for snow concentration, Ns
    iiNi,   & ! Hydrometeor array index for ice concentration, Ni
    iiNg      ! Hydrometeor array index for graupel concentration, Ng

  ! Scalar quantities
  integer, public :: & 
    iisclr_rt, iisclr_thl, iisclr_CO2, & ! [kg/kg]/[K]/[1e6 mol/mol]
    iiedsclr_rt, iiedsclr_thl, iiedsclr_CO2 ! "    "

  ! Logical fields
  logical, dimension(:), allocatable, public :: &
    l_frozen_hm, & ! if true, then the hydrometeor is frozen; otherwise liquid
    l_mix_rat_hm   ! if true, then the quantity is a hydrometeor mixing ratio

  character(len=10), dimension(:), allocatable, public :: & 
    hydromet_list


  real( kind = core_rknd ), dimension(:), allocatable, public :: &
    hydromet_tol    ! Tolerance values for all hydrometeors    [units vary]


  ! Latin hypercube indices / Correlation array indices
  integer, public :: &
    iiPDF_chi = -1, &
    iiPDF_eta = -1, &
    iiPDF_w   = -1

  integer, public :: &
   iiPDF_rr = -1, &
   iiPDF_rs = -1, &
   iiPDF_ri = -1, &
   iiPDF_rg = -1

  integer, public :: &
   iiPDF_Nr  = -1, &
   iiPDF_Ns  = -1, &
   iiPDF_Ni  = -1, &
   iiPDF_Ng  = -1, &
   iiPDF_Ncn = -1

  private ! Default Scope

!===============================================================================

end module array_index
