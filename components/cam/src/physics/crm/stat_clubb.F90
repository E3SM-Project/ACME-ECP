! $Id: stat_clubb.F90 688 2010-03-02 16:46:36Z dschanen@uwm.edu $
module stat_clubb
#ifdef CLUBB_CRM

  implicit none

  public :: stats_clubb_update

  private :: LIN_INT

  private

  contains 

  subroutine stats_clubb_update( upwp, vpwp, up2, vp2, wprtp, wpthlp, &
    wp2, wp3, rtp2, thlp2, rtpthlp, cloud_frac, rcm, um, vm )

  use grid, only: nx, ny, nzm, nz, dimx1_s, dimx2_s, dimy1_s, dimy2_s, &
    zi, z ! Interface and pressure levels, respectively.    [m]

#ifndef CRM
  use hbuffer, only: hbuf_put, hbuf_avg_put
#endif

  implicit none

  real, dimension(nx, ny, nz), intent(in) :: &
    upwp,        &! u'w'                          [m^2/s^2]
    vpwp,        &! u'w'                          [m^2/s^2]
    up2,         &! u'^2                          [m^2/s^2]
    vp2,         &! v'^2                          [m^2/s^2]
    wprtp,       &! w' r_t'                       [(m kg)/(s kg)]
    wpthlp,      &! w' th_l'                      [(m K)/s]
    wp2,         &! w'^2                          [m^2/s^2]
    rtp2,        &! r_t'^2                        [(kg/kg)^2]
    thlp2,       &! th_l'^2                       [K^2]
    rtpthlp,     &! r_t' th_l'                    [(kg K)/kg]
    cloud_frac,  &! Cloud Fraction                [-]
    rcm           ! Cloud water                   [kg/kg]

  ! w'^3 is requires additional ghost points on the x and y dimension
  real, dimension(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nz), intent(in) :: &
    wp3,&    ! w'^3                       [m^3/s^3]
    um, &    ! x-wind                     [m/s]
    vm       ! y-wind                     [m/s]

  ! Local variables
  real, dimension(nzm) :: &
    upwp_avg,   &
    vpwp_avg,   &
    up2_avg,    &
    vp2_avg,    &
    wprtp_avg,  &
    wpthlp_avg, &
    wp2_avg,    &
    thlp2_avg,  &
    rtp2_avg,   &
    rtpthlp_avg,&
    sigma_sqd_w_avg, &
    Kh_zt_avg,  &
    tau_zm_avg

  real :: factor_xy

  integer :: i, j, k

  !---------------------------------------------------------
  ! CLUBB variables
  ! Notes: The variables located on the vertical velocity levels 
  ! must be interpolated for the stats grid, which is on the pressure levels.
  ! -dschanen 21 Jul 2008
  factor_xy = 1. / real( nx*ny )

  upwp_avg   = 0.0
  vpwp_avg   = 0.0
  vp2_avg    = 0.0
  up2_avg    = 0.0
  wprtp_avg  = 0.0
  wpthlp_avg = 0.0
  wp2_avg    = 0.0

  thlp2_avg   = 0.0
  rtp2_avg    = 0.0
  rtpthlp_avg = 0.0

  do i = 1, nx
    do j = 1, ny
      do k = 1, nzm
        upwp_avg(k) = upwp_avg(k) &
          + lin_int( upwp(i,j,k+1), upwp(i,j,k), zi(k+1), zi(k), z(k) )
        vpwp_avg(k) = vpwp_avg(k) &
          + lin_int( vpwp(i,j,k+1), vpwp(i,j,k), zi(k+1), zi(k), z(k) )
        vp2_avg(k) = vp2_avg(k) &
          + lin_int( vp2(i,j,k+1), vp2(i,j,k), zi(k+1), zi(k), z(k) )
        up2_avg(k) = up2_avg(k) &
          + lin_int( up2(i,j,k+1), up2(i,j,k), zi(k+1), zi(k), z(k) )
        wprtp_avg(k) = wprtp_avg(k) &
          + lin_int( wprtp(i,j,k+1), wprtp(i,j,k), zi(k+1), zi(k), z(k) )
        wpthlp_avg(k) = wpthlp_avg(k) &
          + lin_int( wpthlp(i,j,k+1), wpthlp(i,j,k), zi(k+1), zi(k), z(k) )
        wp2_avg(k) = wp2_avg(k) &
          + lin_int( wp2(i,j,k+1), wp2(i,j,k), zi(k+1), zi(k), z(k) )
        rtp2_avg(k) = rtp2_avg(k) &
          + lin_int( rtp2(i,j,k+1), rtp2(i,j,k), zi(k+1), zi(k), z(k) )
        thlp2_avg(k) = thlp2_avg(k) &
          + lin_int( thlp2(i,j,k+1), thlp2(i,j,k), zi(k+1), zi(k), z(k) )
        rtpthlp_avg(k) = rtpthlp_avg(k) &
          + lin_int( rtpthlp(i,j,k+1), rtpthlp(i,j,k), zi(k+1), zi(k), z(k) )

      end do
    end do
  end do

#ifndef CRM
  ! Velocity grid variables
  call hbuf_put('UPWP', upwp_avg, factor_xy)
  call hbuf_put('VPWP', vpwp_avg, factor_xy)
  call hbuf_put('VP2', vp2_avg, factor_xy)
  call hbuf_put('UP2', up2_avg, factor_xy)
  call hbuf_put('WPRTP', wprtp_avg, factor_xy)
  call hbuf_put('WPTHLP', wpthlp_avg, factor_xy)
  call hbuf_put('WP2', wp2_avg, factor_xy)
  call hbuf_put('RTP2', rtp2_avg, factor_xy)
  call hbuf_put('THLP2', thlp2_avg, factor_xy)
  call hbuf_put('RTPTHLP', rtpthlp_avg, factor_xy)

  ! CLUBB thermodynamic grid varibles (SAM pressure levels + ghost point)
  call hbuf_avg_put('CLD_FRAC', cloud_frac(1:nx,1:ny,2:nz), 1,nx, 1,ny, nzm, 1.)
  call hbuf_avg_put('RCM', rcm(1:nx,1:ny,2:nz), 1,nx, 1,ny, nzm, 1.)
  call hbuf_avg_put('UM', um(1:nx,1:ny,2:nz), 1,nx, 1,ny, nzm, 1.)
  call hbuf_avg_put('VM', vm(1:nx,1:ny,2:nz), 1,nx, 1,ny, nzm, 1.)
  call hbuf_avg_put('WP3', wp3(1:nx,1:ny,2:nz), 1,nx, 1,ny, nzm, 1.)
#endif /*CRM*/

  return
  end subroutine stats_clubb_update

!-------------------------------------------------------------------------------
FUNCTION LIN_INT( var_high, var_low, height_high, height_low, height_int )

! This function computes a linear interpolation of the value of variable.
! Given two known values of a variable at two height values, the value
! of that variable at a height between those two height levels (rather 
! than a height outside of those two height levels) is computed.
!
! Here is a diagram:
!
!  ################################ Height high, know variable value
!
!
!
!  -------------------------------- Height to be interpolated to; linear interpolation
!
!
!
!
!
!  ################################ Height low, know variable value
!
!
! FORMULA:
!
! variable(@ Height interpolation) =
!
! [ (variable(@ Height high) - variable(@ Height low)) / (Height high - Height low) ]
! * (Height interpolation - Height low)  +  variable(@ Height low)

! Author: Brian Griffin, UW-Milwaukee

! References: None

IMPLICIT NONE

! Input Variables
REAL, INTENT(IN):: var_high
REAL, INTENT(IN):: var_low
REAL, INTENT(IN):: height_high
REAL, INTENT(IN):: height_low
REAL, INTENT(IN):: height_int

! Output Variable
REAL:: lin_int

lin_int = ( var_high - var_low ) / ( height_high - height_low ) &
         * ( height_int - height_low ) + var_low


END FUNCTION LIN_INT

#endif /* CLUBB_CRM */
end module stat_clubb
