module zero_mod
  implicit none

contains

  subroutine zero(ncrms,icrm)

    use vars
    use microphysics, only : total_water

    implicit none
    integer, intent(in) :: ncrms,icrm

    integer k

    dudt(icrm,:,:,:,na(icrm)) = 0.
    dvdt(icrm,:,:,:,na(icrm)) = 0.
    dwdt(icrm,:,:,:,na(icrm)) = 0.
    misc(icrm,:,:,:) = 0.

  end

end module zero_mod
