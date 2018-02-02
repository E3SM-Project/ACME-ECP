module zero_mod
  implicit none

contains

  subroutine zero(ncrms)

    use vars
    use microphysics, only : total_water

    implicit none
    integer, intent(in) :: ncrms

    dudt(:,:,:,:,na) = 0.
    dvdt(:,:,:,:,na) = 0.
    dwdt(:,:,:,:,na) = 0.
    misc(:,:,:,:) = 0.

  end

end module zero_mod
