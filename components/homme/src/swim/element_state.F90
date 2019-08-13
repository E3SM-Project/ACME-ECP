#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module element_state

  use kinds,                  only: real_kind, long_kind, int_kind
  use dimensions_mod,         only: np, npsq, nlev, nlevp

  implicit none
  private
  integer, public, parameter :: timelevels = 3


  type, public :: elem_state_t

    ! prognostic variables for shallow-water solver
     real (kind=real_kind) :: p(np,np,nlev,timelevels)
     real (kind=real_kind) :: ps(np,np)                               ! surface geopotential
     real (kind=real_kind) :: gradps(np,np,2)                         ! gradient of surface geopotential
     real (kind=real_kind) :: v(np,np,2,nlev,timelevels)              ! contravarient comp

    ! Variables for incoming dynamics state used to calculate 
    ! dynamics tendencies to physics when fv_nphys>0
    real (kind=real_kind) :: T_in  (np,np,nlev)         ! temperature
    real (kind=real_kind) :: V_in  (np,np,2,nlev)       ! velocity
    real (kind=real_kind) :: Q_in  (np,np,nlev,qsize_d) ! Tracer concentration
    real (kind=real_kind) :: dp_in (np,np,nlev)         ! pressure thickness

  end type elem_state_t

  !___________________________________________________________________
  type, public :: derived_state_t
  end type derived_state_t

  type, public :: elem_accum_t
  end type elem_accum_t


contains
end module 
