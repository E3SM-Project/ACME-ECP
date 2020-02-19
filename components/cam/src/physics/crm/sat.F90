module sat_mod
  use params, only: crm_rknd
  implicit none
  public :: esatw_crm
  public :: qsatw_crm
  public :: dtesatw_crm
  public :: dtqsatw_crm
  public :: esati_crm
  public :: qsati_crm
  public :: dtesati_crm
  public :: dtqsati_crm
contains
  ! Saturation vapor pressure and mixing ratio.
  ! Based on Flatau et.al, (JAM, 1992:1507) - valid for T > -80C
  ! sat. vapor over ice below -80C - used Murphy and Koop (2005)
  ! For water below -80C simply assumed esw/esi = 2.
  ! des/dT below -80C computed as a finite difference of es
  function esatw_crm(t)
#if defined(_OPENACC)
    !$acc routine seq
#elif defined(_OPENMP)
    !!!$omp declare target
#endif
    implicit none
    real(crm_rknd) :: esatw_crm
    real(crm_rknd), intent(in) :: t  ! temperature (K)
    real(crm_rknd) :: a0,a1,a2,a3,a4,a5,a6,a7,a8
    real(crm_rknd) :: dt
    a0 = 6.105851
    a1 = 0.4440316
    a2 = 0.1430341e-1
    a3 = 0.2641412e-3
    a4 = 0.2995057e-5
    a5 = 0.2031998e-7
    a6 = 0.6936113e-10
    a7 = 0.2564861e-13
    a8 = -0.3704404e-15
    dt = t-273.16
    if(dt.gt.-80.) then
      esatw_crm = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
    else
      esatw_crm = 2.*0.01*exp(9.550426 - 5723.265/t + 3.53068*Log(t) - 0.00728332*t)
    end if
  end function esatw_crm


  function qsatw_crm(t,p)
#if defined(_OPENACC)
    !$acc routine seq
#elif defined(_OPENMP)
    !!!$omp declare target
#endif
    implicit none
    real(crm_rknd) :: qsatw_crm
    real(crm_rknd), intent(in) :: t  ! temperature (K)
    real(crm_rknd), intent(in) :: p  ! pressure    (mb)
    real(crm_rknd) :: esat_crm
    esat_crm = esatw_crm(t)
    qsatw_crm = 0.622 * esat_crm/max(esat_crm,p-esat_crm)
  end function qsatw_crm


  function dtesatw_crm(t)
#if defined(_OPENACC)
    !$omp routine seq
#elif defined(_OPENMP)
    !!!$omp declare target
#endif
    implicit none
    real(crm_rknd) :: dtesatw_crm
    real(crm_rknd), intent(in) :: t  ! temperature (K)
    real(crm_rknd) :: a0,a1,a2,a3,a4,a5,a6,a7,a8
    real(crm_rknd) :: dt

    a0 = 0.443956472
    a1 = 0.285976452e-1
    a2 = 0.794747212e-3
    a3 = 0.121167162e-4
    a4 = 0.103167413e-6
    a5 = 0.385208005e-9
    a6 = -0.604119582e-12
    a7 = -0.792933209e-14
    a8 = -0.599634321e-17
    dt = t-273.16
    if(dt.gt.-80.) then
      dtesatw_crm = a0 + dt* (a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
    else
      dtesatw_crm = esatw_crm(t+1)-esatw_crm(t)
    end if
  end function dtesatw_crm


  function dtqsatw_crm(t,p)
#if defined(_OPENACC)
    !$acc routine seq
#elif defined(_OPENMP)
    !!!$omp declare target
#endif
    implicit none
    real(crm_rknd) :: dtqsatw_crm
    real(crm_rknd), intent(in) :: t  ! temperature (K)
    real(crm_rknd), intent(in) :: p  ! pressure    (mb)
    dtqsatw_crm = 0.622*dtesatw_crm(t)/p
  end function dtqsatw_crm


  function esati_crm(t)
#if defined(_OPENACC)
    !$acc routine seq
#elif defined(_OPENMP)
    !!!$omp declare target
#endif
    implicit none
    real(crm_rknd) :: esati_crm
    real(crm_rknd), intent(in) :: t  ! temperature (K)
    real(crm_rknd) :: a0,a1,a2,a3,a4,a5,a6,a7,a8
    real(crm_rknd) :: dt

    a0 = 6.11147274
    a1 = 0.503160820
    a2 = 0.188439774e-1
    a3 = 0.420895665e-3
    a4 = 0.615021634e-5
    a5 = 0.602588177e-7
    a6 = 0.385852041e-9
    a7 = 0.146898966e-11
    a8 = 0.252751365e-14
    
    dt = t-273.16
    if(dt.gt.-80.) then
      esati_crm = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
    else
      esati_crm = 0.01*exp(9.550426 - 5723.265/t + 3.53068*Log(t) - 0.00728332*t)
    end if
  end function esati_crm


  function qsati_crm(t,p)
#if defined(_OPENACC)
    !$acc routine seq
#elif defined(_OPENMP)
    !!!$omp declare target
#endif
    implicit none
    real(crm_rknd) :: qsati_crm
    real(crm_rknd), intent(in) :: t  ! temperature (K)
    real(crm_rknd), intent(in) :: p  ! pressure    (mb)
    real(crm_rknd) :: esat_crm
    esat_crm=esati_crm(t)
    qsati_crm=0.622 * esat_crm/max(esat_crm,p-esat_crm)
  end function qsati_crm


  function dtesati_crm(t)
#if defined(_OPENACC)
    !$acc routine seq
#elif defined(_OPENMP)
    !!!$omp declare target
#endif
    implicit none
    real(crm_rknd) :: dtesati_crm
    real(crm_rknd), intent(in) :: t  ! temperature (K)
    real(crm_rknd) :: a0,a1,a2,a3,a4,a5,a6,a7,a8
    real(crm_rknd) :: dt

    a0 = 0.503223089
    a1 = 0.377174432e-1
    a2 = 0.126710138e-2
    a3 = 0.249065913e-4
    a4 = 0.312668753e-6
    a5 = 0.255653718e-8
    a6 = 0.132073448e-10
    a7 = 0.390204672e-13
    a8 = 0.497275778e-16

    dt = t-273.16
    if(dt.gt.-80.) then
      dtesati_crm = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
    else
      dtesati_crm = esati_crm(t+1.)-esati_crm(t)
    end if
  end function dtesati_crm


  function dtqsati_crm(t,p)
#if defined(_OPENACC)
    !$acc routine seq
#elif defined(_OPENMP)
    !!!$omp declare target
#endif
    implicit none
    real(crm_rknd) :: dtqsati_crm
    real(crm_rknd), intent(in) :: t  ! temperature (K)
    real(crm_rknd), intent(in) :: p  ! pressure    (mb)
    dtqsati_crm=0.622*dtesati_crm(t)/p
  end function dtqsati_crm

end module sat_mod
