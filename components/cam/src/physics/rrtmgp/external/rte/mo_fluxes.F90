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
! This module is for packaging output quantities from RRTMGP based on spectral flux profiles 
!    This implementation is the most basic and computes only broadband fluxes. 
!    The intent is for users to extend it as required (see mo_flxues_byband) 
!
module mo_fluxes
  use mo_rte_kind,   only: wp
  use mo_spectral_disc, only: ty_spectral_disc
  implicit none

  !----------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------
  ! Output from radiation calculations 
  !   Data components are pointers so results can be written directly into memory 
  !   reduce() function accepts spectral flux profiles 
  type :: ty_fluxes
    real(wp), dimension(:,:), pointer :: flux_up => NULL(), flux_dn => NULL()
    real(wp), dimension(:,:), pointer :: flux_net => NULL()    ! Net (down - up)  
    real(wp), dimension(:,:), pointer :: flux_dn_dir => NULL() ! Direct flux down 
  contains
    procedure, public :: reduce      => reduce_base
    procedure, public :: are_desired => are_desired_base
  end type ty_fluxes
contains
  ! --------------------------------------------------------------------------------------
  ! Broadband fluxes -- simply sum over the spectral dimension
  !
  function reduce_base(this, gpt_flux_up, gpt_flux_dn, spectral_disc, top_at_1, gpt_flux_dn_dir) result(error_msg)
    class(ty_fluxes),                  intent(inout) :: this
    real(kind=wp), dimension(:,:,:),   intent(in   ) :: gpt_flux_up ! Fluxes by gpoint [W/m2](ncol, nlay+1, ngpt)
    real(kind=wp), dimension(:,:,:),   intent(in   ) :: gpt_flux_dn ! Fluxes by gpoint [W/m2](ncol, nlay+1, ngpt)
    class(ty_spectral_disc),           intent(in   ) :: spectral_disc  !< derived type with spectral information
    logical,                           intent(in   ) :: top_at_1
    real(kind=wp), dimension(:,:,:), optional, & 
                                       intent(in   ) :: gpt_flux_dn_dir! Direct flux down 
    character(len=128)                               :: error_msg
    ! ------
    integer :: ncol, nlev, ngpt
    
    ! ------
    ncol = size(gpt_flux_up, DIM=1)
    nlev = size(gpt_flux_up, DIM=2)
    ngpt = size(gpt_flux_up, DIM=3)
    error_msg = ""
 
    ! 
    ! Check array sizes
    !  Input arrays 
    !
    if(any([size(gpt_flux_dn, 1) /= ncol, & 
            size(gpt_flux_dn, 2) /= nlev, & 
            size(gpt_flux_dn, 3) /= ngpt])) then 
      error_msg = "reduce: gpt_flux_dn array incorrectly sized"
      return
    end if
    if(present(gpt_flux_dn_dir)) then 
    if(any([size(gpt_flux_dn_dir, 1) /= ncol, & 
            size(gpt_flux_dn_dir, 2) /= nlev, & 
            size(gpt_flux_dn_dir, 3) /= ngpt])) then 
        error_msg = "reduce: gpt_flux_dn_dir array incorrectly sized"
        return
      end if
    end if
    !
    ! Output arrays 
    !
    if(associated(this%flux_up)) then 
      if(any([size(this%flux_up, 1), size(this%flux_up, 2)] /= [ncol,nlev])) then 
        error_msg = 'reduce: flux_up array incorrectly sized'
        return
      end if 
    end if 
    if(associated(this%flux_dn)) then 
      if(any([size(this%flux_dn, 1), size(this%flux_dn, 2)] /= [ncol,nlev])) then 
        error_msg = 'reduce: flux_dn array incorrectly sized'
        return
      end if 
    end if 
    if(associated(this%flux_net)) then 
      if(any([size(this%flux_net, 1), size(this%flux_net, 2)] /= [ncol,nlev])) then 
        error_msg = 'reduce: flux_net array incorrectly sized'
        return
      end if 
    end if 
    if(associated(this%flux_dn_dir)) then 
      if(any([size(this%flux_dn_dir, 1), size(this%flux_dn_dir, 2)] /= [ncol,nlev])) then 
        error_msg = 'reduce: flux_dn_dir array incorrectly sized'
        return
      end if 
    end if 
    !
    ! Self-consistency -- shouldn't be asking for direct beam flux if it isn't supplied 
    if(associated(this%flux_dn_dir) .and. .not. present(gpt_flux_dn_dir)) then 
      error_msg = "reduce: requesting direct downward flux but this hasn't been supplied"
      return
    end if 
    
    !
    ! Broadband fluxes
    !
    if(associated(this%flux_up)) then 
        this%flux_up = sum(gpt_flux_up, DIM=3)
    end if

    if(associated(this%flux_dn)) then 
        this%flux_dn = sum(gpt_flux_dn, DIM=3)
    end if
    
    if(associated(this%flux_dn_dir)) then 
      this%flux_dn_dir = sum(gpt_flux_dn_dir, DIM=3)
    end if

    if(associated(this%flux_net)) then
      ! 
      !  Reuse down and up results if possible
      ! 
      if(associated(this%flux_dn) .and. associated(this%flux_up)) then 
        this%flux_net = this%flux_dn - this%flux_up 
      else 
        this%flux_net = sum(gpt_flux_dn - gpt_flux_up, dim = 3) 
      end if 
    end if 
  end function reduce_base
  ! --------------------------------------------------------------------------------------
  ! Are any fluxes desired from this set of g-point fluxes? We can tell because memory will 
  !   be allocated for output 
  !
  function are_desired_base(this) 
    class(ty_fluxes), intent(in   ) :: this
    logical                         :: are_desired_base
    
    are_desired_base = any([associated(this%flux_up),     & 
                            associated(this%flux_dn),     &
                            associated(this%flux_dn_dir), & 
                            associated(this%flux_net)])
  end function are_desired_base
end module mo_fluxes

