module crm_rad_module
   ! Module to encapsulate data and methods specific to radiation. This exists
   ! as a separate module with separate derived types because radiation is
   ! handled in a special way when interacting between the GCM and CRM.
   ! The radiation calculations may be on a separate grid than the CRM, and the
   ! quantities may be aggregated and updated between calls in a different way than
   ! the other diagnostic quantities. This module should also contain methods to
   ! update the radiation in the future, should we choose to put the radiation
   ! calculations on the CRM.
   use params,       only: crm_rknd
   use shr_kind_mod, only: r8 => shr_kind_r8

   implicit none

   public crm_rad_type

   type crm_rad_type
      ! Quantities used by the radiation code. Note that these are strange in that they are 
      ! time-averages, but spatially-resolved.
      real(crm_rknd), pointer :: t_rad  (:,:,:,:) ! rad temperature
      real(crm_rknd), pointer :: qv_rad (:,:,:,:) ! rad vapor
      real(crm_rknd), pointer :: qc_rad (:,:,:,:) ! rad cloud water
      real(crm_rknd), pointer :: qi_rad (:,:,:,:) ! rad cloud ice
      real(crm_rknd), pointer :: cld_rad(:,:,:,:) ! rad cloud fraction

      ! Only relevant when using 2-moment microphysics
      real(crm_rknd), pointer :: nc_rad (:,:,:,:) ! rad cloud droplet number (#/kg)
      real(crm_rknd), pointer :: ni_rad (:,:,:,:) ! rad cloud ice crystal number (#/kg)
      real(crm_rknd), pointer :: qs_rad (:,:,:,:) ! rad cloud snow (kg/kg)
      real(crm_rknd), pointer :: ns_rad (:,:,:,:) ! rad cloud snow crystal number (#/kg)
   contains
      procedure, public :: initialize=>crm_state_initialize
      procedure, public :: finalize=>crm_state_finalize
   end type crm_rad_type

contains

   !------------------------------------------------------------------------------------------------
   ! Type-bound procedures for crm_state_type
   subroutine crm_rad_initialize(this)
      class(crm_rad_type), intent(inout) :: this

      ! Nullify pointers
      this%t_rad => null()
      this%qv_rad => null()
      this%qi_rad => null()
      this%cld_rad => null()

      this%nc_rad => null()
      this%ni_rad => null()
      this%qs_rad => null()
      this%ns_rad => null()

   end subroutine crm_rad_initialize
   !------------------------------------------------------------------------------------------------
   subroutine crm_rad_finalize(this)
      class(crm_rad_type), intent(inout) :: this

      ! Nullify pointers
      this%t_rad => null()
      this%qv_rad => null()
      this%qi_rad => null()
      this%cld_rad => null()

      this%nc_rad => null()
      this%ni_rad => null()
      this%qs_rad => null()
      this%ns_rad => null()

   end subroutine crm_rad_finalize
   !------------------------------------------------------------------------------------------------

end module crm_rad_module
