module zero_mod
  implicit none

contains

  subroutine zero(ncrms,icrm)
    use vars
    use microphysics, only : total_water
    implicit none
    integer, intent(in) :: ncrms,icrm
    integer k
    dudt(:,:,:,na,icrm) = 0.
    dvdt(:,:,:,na,icrm) = 0.
    dwdt(:,:,:,na,icrm) = 0.
    misc(:,:,:) = 0.
  end

end module zero_mod
