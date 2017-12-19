! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2015,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!
!
! This module is for packaging output quantities from RRTMGP based on spectral flux profiles
!    This implementation provides band-by-band flux profiles
!
module mo_fluxes_byband
  use mo_rte_kind,   only: wp
  use mo_fluxes,        only: ty_fluxes
  use mo_spectral_disc, only: ty_spectral_disc
  implicit none

  ! Output from radiation calculations
  !   Data components are pointers so results can be written directly into memory
  !   reduce() function accepts spectral flux profiles
  type, extends(ty_fluxes) :: ty_fluxes_byband
    real(wp), dimension(:,:,:), pointer :: bnd_flux_up => NULL(), & ! Band-by-band fluxes
                                           bnd_flux_dn => NULL()    ! (ncol, nlev, nband)
    real(wp), dimension(:,:,:), pointer :: bnd_flux_net => NULL()   ! Net (down - up)
    real(wp), dimension(:,:,:), pointer :: bnd_flux_dn_dir => NULL() ! Direct flux down 
  contains
    procedure :: reduce => reduce_byband
    procedure :: are_desired => are_desired_byband
  end type ty_fluxes_byband
contains
  ! --------------------------------------------------------------------------------------
  function reduce_byband(this, gpt_flux_up, gpt_flux_dn, spectral_disc, top_at_1, gpt_flux_dn_dir) result(error_msg)
    class(ty_fluxes_byband),           intent(inout) :: this
    real(kind=wp), dimension(:,:,:),   intent(in   ) :: gpt_flux_up ! Fluxes by gpoint [W/m2](ncol, nlay+1, ngpt)
    real(kind=wp), dimension(:,:,:),   intent(in   ) :: gpt_flux_dn ! Fluxes by gpoint [W/m2](ncol, nlay+1, ngpt)
    class(ty_spectral_disc),            intent(in   ) :: spectral_disc  !< derived type with spectral information
    logical,                           intent(in   ) :: top_at_1
    real(kind=wp), dimension(:,:,:), optional, & 
                                       intent(in   ) :: gpt_flux_dn_dir! Direct flux down 
    character(len=128)                               :: error_msg
    ! ------
    integer :: ncol, nlev, ngpt, nbnd
    integer :: bnd, gpt
    integer :: gpt_lim(2)
    logical :: up_lims_ok = .false., dn_lims_ok = .false.
    ! ------
    ncol = size(gpt_flux_up, DIM=1)
    nlev = size(gpt_flux_up, DIM=2)
    ngpt = size(gpt_flux_up, DIM=3)
    nbnd = spectral_disc%get_nband()


    ! Compute broadband fluxes
    !   This also checks that input arrays are consistently sized
    ! 
    error_msg = this%ty_fluxes%reduce(gpt_flux_up, gpt_flux_dn, spectral_disc, top_at_1, gpt_flux_dn_dir)
    if(error_msg /= '') return

    ! Check sizes of output arrays 
    if(associated(this%bnd_flux_up)) then
      if(any([size(this%bnd_flux_up, 1) /= ncol,  &
              size(this%bnd_flux_up, 2) /= nlev,  &
              size(this%bnd_flux_up, 3) /= nbnd])) then 
        error_msg = "reduce: bnd_flux_up array incorrectly sized (can't compute net flux either)"
        return
      end if
    end if
    if(associated(this%bnd_flux_dn)) then
      if(any([size(this%bnd_flux_dn, 1) /= ncol,  &
              size(this%bnd_flux_dn, 2) /= nlev,  &
              size(this%bnd_flux_dn, 3) /= nbnd])) then 
        error_msg = "reduce: bnd_flux_dn array incorrectly sized (can't compute net flux either)"
        return
      end if
    end if
    if(associated(this%bnd_flux_dn_dir)) then
      if(any([size(this%bnd_flux_dn_dir, 1) /= ncol,  &
              size(this%bnd_flux_dn_dir, 2) /= nlev,  &
              size(this%bnd_flux_dn_dir, 3) /= nbnd])) then 
        error_msg = "reduce: bnd_flux_dn_dir array incorrectly sized"
        return
      end if
    end if
    if(associated(this%bnd_flux_net)) then
      if(any([size(this%bnd_flux_net, 1) /= ncol,  &
              size(this%bnd_flux_net, 2) /= nlev,  &
              size(this%bnd_flux_net, 3) /= nbnd])) then 
        error_msg = "reduce: bnd_flux_net array incorrectly sized (can't compute net flux either)"
        return
      end if
    end if
    !
    ! Self-consistency -- shouldn't be asking for direct beam flux if it isn't supplied 
    if(associated(this%bnd_flux_dn_dir) .and. .not. present(gpt_flux_dn_dir)) then 
      error_msg = "reduce: requesting bnd_flux_dn_dir but direct flux hasn't been supplied"
      return
    end if 
        
    ! -------
    ! Band-by-band fluxes
    ! Up flux
    if(associated(this%bnd_flux_up)) then
      do bnd = 1, nbnd
        gpt_lim(:) = spectral_disc%convert_band2gpt(bnd)
        this%bnd_flux_up(:,:,bnd) = sum(gpt_flux_up(:,:,gpt_lim(1):gpt_lim(2)), dim=3)
      end do
    end if 

    ! -------
    ! Down flux
    if(associated(this%bnd_flux_dn)) then
      do bnd = 1, nbnd
        gpt_lim(:) = spectral_disc%convert_band2gpt(bnd)
        this%bnd_flux_dn(:,:,bnd) = sum(gpt_flux_dn(:,:,gpt_lim(1):gpt_lim(2)), dim=3)
      end do  
    end if

    if(associated(this%bnd_flux_dn_dir)) then
      do bnd = 1, nbnd
        gpt_lim(:) = spectral_disc%convert_band2gpt(bnd)
        this%bnd_flux_dn_dir(:,:,bnd) = sum(gpt_flux_dn_dir(:,:,gpt_lim(1):gpt_lim(2)), dim=3)
      end do  
    end if

    ! -------
    ! Net flux
    !
    if(associated(this%bnd_flux_net)) then
      ! 
      !  Reuse down and up results if possible
      ! 
      if(associated(this%bnd_flux_dn) .and. associated(this%bnd_flux_up)) then 
        this%bnd_flux_net = this%bnd_flux_dn - this%bnd_flux_up 
      else 
        do bnd = 1, nbnd
          gpt_lim(:) = spectral_disc%convert_band2gpt(bnd)
          this%bnd_flux_net(:,:,bnd) = sum(gpt_flux_dn(:,:,gpt_lim(1):gpt_lim(2)), dim=3) - & 
                                       sum(gpt_flux_up(:,:,gpt_lim(1):gpt_lim(2)), dim=3)
        end do  
      end if 
    end if 
  end function reduce_byband
  ! --------------------------------------------------------------------------------------
  ! Are any fluxes desired from this set of g-point fluxes? We can tell because memory will
  !   be allocated for output
  !
  function are_desired_byband(this)
    class(ty_fluxes_byband), intent(in   ) :: this
    logical                                :: are_desired_byband

    are_desired_byband = any([associated(this%bnd_flux_up),     &
                              associated(this%bnd_flux_dn),     &
                              associated(this%bnd_flux_dn_dir), &
                              associated(this%bnd_flux_net),    &
                              this%ty_fluxes%are_desired()])
  end function are_desired_byband

end module mo_fluxes_byband
