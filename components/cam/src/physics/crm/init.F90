module init_mod
	implicit none

contains

  ! Initialize some arrays from vars module:
  subroutine init(ncrms,icrm)
    use vars
    use crmtracers
    implicit none
    integer, intent(in) :: ncrms,icrm

    fzero(icrm,:,:) = 0.

    ttend(icrm,:) = 0.
    qtend(icrm,:) = 0.
    wsub(icrm,:) = 0.
    !unudge = 0.
    !vnudge = 0.
    !tnudge = 0.
    !qnudge = 0.
    !qlsvadv = 0.
    !tlsvadv = 0.
    !ulsvadv = 0.
    !vlsvadv = 0.
    !qstor = 0.
    !tstor = 0.
    !ustor = 0.
    !vstor = 0.
    !qtostor = 0.

    !radlwup = 0.
    !radlwdn = 0.
    !radswup = 0.
    !radswdn = 0.
    !radqrlw = 0.
    !radqrsw = 0.

    tlat(icrm,:) = 0.
    tlatqi(icrm,:) = 0.
    tadv(icrm,:) = 0.
    tdiff(icrm,:) = 0.
    qifall(icrm,:) = 0.
    qpfall(icrm,:) = 0.

    trwle (icrm,:,:) = 0.
    trwsb (icrm,:,:) = 0.
    tradv (icrm,:,:) = 0.
    trdiff(icrm,:,:) = 0.
    trphys(icrm,:,:) = 0.

    !gamt0 = 0.
    !gamq0 = 0.

  end subroutine init



end module init_mod
