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
      use sgs, only: setperturb_sgs
      use params, only: crm_rknd
      use RNG_MT

      implicit none

      integer, intent(in) :: iseed

      integer i,j,k
      real(crm_rknd) :: rand_perturb
      real(crm_rknd) :: t02(nzm)

      integer, parameter :: num_perturb_layers = 5   ! number of layers to add noise

      call setperturb_sgs(0)  ! set sgs fields

      !!! set the seed (based on the global physics column index)
      call RNG_MT_set_seed(iseed)

      t02 = 0.0
      do k=1,nzm

         do j=1,ny
            do i=1,nx

               !!! Generate a uniform random number in interval (0,1)
               call RNG_MT_gen_rand(rand_perturb)

               !!! convert perturbation range from (0,1) to (-1,1)
               rand_perturb = 1.-2.*rand_perturb

               !!! apply perturbation to temperature field
               if(k.le.num_perturb_layers) then
                  t(i,j,k)=t(i,j,k)+0.02*rand_perturb*(6-k)
               endif
               t02(k) = t02(k) + t(i,j,k)/(nx*ny)

            end do ! i
         end do ! j

         !!! enforce energy conservation
         do j=1, ny
            do i=1, nx
               if(k.le.num_perturb_layers) then
                  t(i,j,k) = t(i,j,k) * t0(k)/t02(k)
               end if
            end do ! i
         end do ! j

      end do ! k

   end subroutine setperturb

end module setperturb_mod
