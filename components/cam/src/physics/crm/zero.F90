
subroutine zero
	
use vars
use microphysics, only : total_water

implicit none
	
integer k
	
dudt(:,:,:,na) = 0.
dvdt(:,:,:,na) = 0.
dwdt(:,:,:,na) = 0.
misc(:,:,:) = 0.

end
