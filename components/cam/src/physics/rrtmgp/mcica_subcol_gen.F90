module mcica_subcol_gen

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2006-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

!----------------------------------------------------------------------------------------
! 
! Purpose: Create McICA stochastic arrays for cloud optical properties.
! Input cloud optical properties directly: cloud optical depth, single
! scattering albedo and asymmetry parameter.  Output will be stochastic
! arrays of these variables.  (longwave scattering is not yet available)
!
! Original code: From RRTMG based on Raisanen et al., QJRMS, 2004.
! 
! Uses the KISS random number generator.
!
! Overlap assumption: maximum-random.
! 
!----------------------------------------------------------------------------------------

use shr_kind_mod,     only: r8 => shr_kind_r8
use shr_RandNum_mod,  only: ShrKissRandGen
use kissvec_mod, only: kissvec
use assertions

implicit none
private
save

public :: mcica_subcol_mask

!========================================================================================
contains
!========================================================================================

subroutine mcica_subcol_mask( &
   ngpt, ncol, pver, changeseed, &
   pmid, cldfrac, iscloudy)

   ! Arrays use CAM vertical index convention: index increases from top to bottom.
   ! This index ordering is assumed in the maximum-random overlap algorithm which starts
   ! at the top of a column and marches down, with each layer depending on the state
   ! of the layer above it.
   !
   ! For GCM mode, changeseed must be offset between LW and SW by at least the
   ! number of subcolumns

   ! arguments
   integer, intent(in) :: ngpt              ! number of subcolumns (g-point intervals)
   integer, intent(in) :: ncol              ! number of columns
   integer, intent(in) :: pver              ! number of levels
   integer, intent(in) :: changeseed        ! if the subcolumn generator is called multiple times, 
                                            ! permute the seed between each call.
   real(r8), intent(in) :: pmid(:,:)        ! layer pressures (Pa)
   real(r8), intent(in) :: cldfrac(:,:)     ! layer cloud fraction
   logical, intent(out) :: iscloudy(:,:,:)  ! flag that says whether a gridbox is cloudy

   ! Local vars
   integer :: iter, icol, ilev, igpt
   real(r8), parameter :: cldmin = 1.0e-80_r8  ! min cloud fraction
   !real(r8) :: cldf(ncol,pver)       ! cloud fraction clipped to cldmin
   !integer  :: kiss_seed(ncol,4)
   !real(r8) :: rand_num(ncol)   ! random number (kissvec)
   !real(r8) :: cdf(ngpt,ncol,pver)   ! random numbers
   real(r8), allocatable :: cldf(:,:)       ! cloud fraction clipped to cldmin
   integer , allocatable :: kiss_seed(:,:)
   real(r8), allocatable :: rand_num(:)   ! random number (kissvec)
   real(r8), allocatable :: cdf(:,:,:)   ! random numbers
   
   character(len=128) :: sub_name = 'mcica_subcol_gen_mask'
   !------------------------------------------------------------------------------------------ 

   ! Make sure inputs are TOA to surface
   call assert(pmid(1,1) < pmid(1,2), trim(sub_name) // ': inputs not TOA to SFC')
   call assert(size(iscloudy, 3) == ngpt, trim(sub_name) // ': wrong dimension for iscloudy')

   allocate(cldf(ncol,pver), kiss_seed(ncol,4), rand_num(ncol), cdf(ncol,pver,ngpt))

   ! clip cloud fraction
   !!$acc parallel loop collapse(2)
   do ilev = 1,pver
      do icol = 1,ncol
         if (cldfrac(icol,ilev) >= cldmin) then
            cldf(icol,ilev) = cldfrac(icol,ilev)
         else
            cldf(icol,ilev) = 0
         end if
      end do
   end do

   ! Create a seed that depends on the state of the columns.
   ! Use pmid from bottom four layers. 
   !!$acc parallel loop
   do icol = 1,ncol
      kiss_seed(icol,1) = (pmid(icol,pver)   - int(pmid(icol,pver)))    * 1000000000
      kiss_seed(icol,2) = (pmid(icol,pver-1) - int(pmid(icol,pver-1)))  * 1000000000
      kiss_seed(icol,3) = (pmid(icol,pver-2) - int(pmid(icol,pver-2)))  * 1000000000
      kiss_seed(icol,4) = (pmid(icol,pver-3) - int(pmid(icol,pver-3)))  * 1000000000
   end do

   ! Advance random number generator by changeseed values
   do iter = 1,changeseed
      call kissvec( &
         kiss_seed(:,1), kiss_seed(:,2), kiss_seed(:,3), kiss_seed(:,4), &
         rand_num, ncol &
      )
   end do

   ! Generate random numbers in each subcolumn at every level; this is done on
   ! the CPU to keep bit for bit reproducibility
   do igpt = 1,ngpt
      do ilev = 1,pver
         call kissvec( &
            kiss_seed(:,1), kiss_seed(:,2), kiss_seed(:,3), kiss_seed(:,4), &
            cdf(:,ilev,igpt), ncol &
         )
      end do
   end do

   ! Maximum-Random overlap
   ! i) pick a random number for top layer.
   ! ii) walk down the column: 
   !    - if the layer above is cloudy, use the same random number as in the layer above
   !    - if the layer above is clear, use a new random number 
   ! NOTE: loop carried dependence prevents parallelization over pver
   !!$acc parallel loop
   do igpt = 1,ngpt
      do ilev = 2,pver
         do icol = 1,ncol
            if (cdf(icol,ilev-1,igpt) > 1._r8 - cldf(icol,ilev-1)) then
               cdf(icol,ilev,igpt) = cdf(icol,ilev-1,igpt)
            else
               cdf(icol,ilev,igpt) = cdf(icol,ilev,igpt) * (1._r8 - cldf(icol,ilev-1))
            end if
         end do
      end do
   end do
 
   ! Set binary cloudy flag
   !!$acc parallel loop collapse(3)
   do igpt = 1,ngpt
      do ilev = 1,pver
         do icol = 1,ncol
            iscloudy(icol,ilev,igpt) = (cdf(icol,ilev,igpt) >= 1._r8 - cldf(icol,ilev))
         end do
      end do
   end do

   deallocate(cldf, kiss_seed, rand_num, cdf)

end subroutine mcica_subcol_mask

end module mcica_subcol_gen
