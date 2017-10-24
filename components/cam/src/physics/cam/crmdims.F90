module crmdims
! whannah - non-SP not compiling without this #ifdef - says it can't find params.mod - not sure how to fix
#ifdef CRM 
    use params, only: crm_rknd
#endif
    implicit none
! whannah - non-SP configuration needs to know the value of crm_rknd
#ifndef CRM 
    integer, parameter :: crm_rknd = selected_real_kind(12) ! 8 byte real  
#endif
    integer, parameter :: nclubbvars = 17

    integer, parameter ::  crm_nx=CRM_NX
    integer, parameter ::  crm_ny=CRM_NY
    integer, parameter ::  crm_nz=CRM_NZ

    integer, parameter ::  crm_nx_rad=CRM_NX_RAD
    integer, parameter ::  crm_ny_rad=CRM_NY_RAD

    real(crm_rknd), parameter :: crm_dx=CRM_DX
    real(crm_rknd), parameter :: crm_dy=crm_dx
    real(crm_rknd), parameter :: crm_dt=CRM_DT

end module crmdims