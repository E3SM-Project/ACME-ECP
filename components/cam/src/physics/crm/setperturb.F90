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
    ! integer, allocatable :: rndm_seed(:)
    ! integer :: rndm_seed_sz
    real(crm_rknd) :: rrr
    real(crm_rknd) :: t02(nzm)

    !!! old method using intrinic RNG
    ! !call ranset_(30*rank)
    ! call random_seed(size=rndm_seed_sz)
    ! allocate(rndm_seed(rndm_seed_sz))

    ! rndm_seed = iseed
    ! call random_seed(put=rndm_seed)

    call setperturb_sgs(0)  ! set sgs fields

    !!! set the seed (based on the global physics column index)
    call RNG_MT_set_seed(iseed)

    t02 = 0.0
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          !!! old method using intrinic RNG
          ! rrr=1.-2.*ranf_()

          !!! Generate a uniform random number in interval (0,1)
          call RNG_MT_gen_rand(rrr)
          rrr = 1.-2.*rrr

          if(k.le.5) then
            t(i,j,k)=t(i,j,k)+0.02*rrr*(6-k)
          endif
          t02(k) = t02(k) + t(i,j,k)/(nx*ny)
        end do
      end do

      ! enforce energy conservation +++mhwang (2012-06)
      do j=1, ny
        do i=1, nx
          if(k.le.5) then
            t(i,j,k) = t(i,j,k) * t0(k)/t02(k)
          end if
        end do
      end do
    end do

    !!! old method using intrinic RNG
    ! deallocate(rndm_seed)

  end

end module setperturb_mod
