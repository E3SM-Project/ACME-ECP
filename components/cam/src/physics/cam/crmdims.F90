#ifdef CRM
module crmdims
	use params, only: crm_rknd
    implicit none

       integer, parameter :: nclubbvars = 17

       integer, parameter ::  crm_nx=CRM_NX, crm_ny=CRM_NY, crm_nz=CRM_NZ
       real(crm_rknd), parameter :: crm_dx=CRM_DX, crm_dy=crm_dx, crm_dt=CRM_DT

end module crmdims
#endif