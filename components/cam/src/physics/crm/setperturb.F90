module setperturb_mod
   use random_mod
   implicit none

contains

   subroutine setperturb(iseed)

      ! Add random noise near the surface to help turbulence develop

      !  This surboutine has been updated for SPCAM5 (Minghuai.Wang@pnnl.gov, April, 2012).
      !  Now the random generator is seeded based on the global column id, which gets rid
      !  of the dependence of the SPCAM results on pcols.

      ! This module was updated to use a Mersenne Twister algorithm, because compiler deependent
      ! issues were identified with the intrinisic random number routines (e.g. random_number())
      ! Walter Hannah - LLNL - Mar 2018

      use vars
      use grid,   only: pres
      use sgs,    only: setperturb_sgs
      use params, only: crm_rknd
      use RNG_MT

      implicit none

      integer, intent(in) :: iseed

      integer i,j,k
      real(crm_rknd) :: rand_perturb
      real(crm_rknd) :: t02 
      real(crm_rknd) :: factor_k
      real(crm_rknd) :: factor_xy
      integer        :: perturb_num_layers      ! number of layers to add perturbations

      integer, parameter :: perturb_t_magnitude = 1.0       ! perturbation t amplitube (max at bottom of the perturbed region)   [K]

      factor_xy = 1./real((nx*ny),crm_rknd)

      call setperturb_sgs(0)  ! set sgs fields

      !!! set the seed
      call RNG_MT_set_seed(iseed)

      !!! find number of layers under some pressure level
      do k = 1,nzm
         if ( pres(k) > 700. ) then
            perturb_num_layers = k
         else
            exit
         end if
      end do

      !--------------------------------------------------------
      ! Apply random liquid static energy (LSE) perturbations
      !--------------------------------------------------------
      do k = 1,perturb_num_layers

         !!! set factor_k so that perturbation tapers upwards
         factor_k = real( perturb_num_layers+1-k ,crm_rknd) / real( perturb_num_layers+1 ,crm_rknd)

         t02 = 0.0
         do j = 1,ny
            do i = 1,nx

               !!! Generate a uniform random number in interval (0,1)
               call RNG_MT_gen_rand(rand_perturb)

               !!! convert perturbation range from (0,1) to (-1,1)
               rand_perturb = 1.-2.*rand_perturb

               !!! apply perturbation to temperature field
               t(i,j,k) = t(i,j,k) + perturb_t_magnitude * rand_perturb * factor_k
               
               t02 = t02 + t (i,j,k)*factor_xy

            end do ! i
         end do ! j

         !!! enforce energy conservation
         do j = 1,ny
            do i = 1,nx
               t (i,j,k) = t (i,j,k) *  t0(k)/t02
            end do ! i
         end do ! j

      end do ! k
      !--------------------------------------------------------
      !--------------------------------------------------------

   end subroutine setperturb

end module setperturb_mod
