module crmdims
! whannah - non-SP not compiling without this #ifdef - says it can't find params.mod - not sure how to fix
#ifdef CRM 
    use params, only: crm_rknd
#endif
    implicit none
! whannah - non-SP not compiling without this #ifdef - says it can't find params.mod - not sure how to fix
#ifndef CRM 
    integer, parameter :: crm_rknd = 8  
#endif
    integer, parameter :: nclubbvars = 17

    integer, parameter ::  crm_nx=CRM_NX, crm_ny=CRM_NY, crm_nz=CRM_NZ
    real(crm_rknd), parameter :: crm_dx=CRM_DX, crm_dy=crm_dx, crm_dt=CRM_DT

end module crmdims