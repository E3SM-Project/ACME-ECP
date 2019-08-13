!---------------------------------------------------------------------------------------------------
! This module contains routines for mapping between dynamics and physics columns
! when dynamics is on the spectral element grid (i.e. GLL) and the physics is on
! a quasi-equal area finite volume grid that evenly divides the dyn element
! 
! The current mapping method utilizes a simple piece-wise linear approach where
! the physics tendencies are uniformly copied to the underlying GLL nodes. This
! comes with the restriction of only using 1 or 2x2 phsyics cells per element, 
! but it also carries the advantage of not needing any information outside of 
! the element. Higher order mapping method can be implemented, but they will 
! require special consideration of the model's vertical coordinate, which is 
! currently based on moist pressure. In either case, dyn_to_fv_phys() would stay
! the same, with simple weighted averaging over element subcells.
! 
! Author: Walter Hannah (LLNL)
!---------------------------------------------------------------------------------------------------
module fv_physics_coupling_mod
  use element_mod,    only: element_t
  use shr_kind_mod,   only: r8=>shr_kind_r8, i4=>shr_kind_i4
  use constituents,   only: pcnst, cnst_name
  use dimensions_mod, only: np, npsq, nelemd, nlev
  use dyn_grid,       only: fv_nphys, fv_physgrid
  use ppgrid,         only: pcols, pver, pverp
  use cam_abortutils, only: endrun
  
  private

  ! These method encapsulate the coupling method for the fv phys grid,
  ! they are used by d_p_coupling() and p_d_coupling() in dp_coupling.F90 
  public :: fv_phys_to_dyn
  public :: fv_phys_to_dyn_topo
  public :: dyn_to_fv_phys

contains
  !=================================================================================================
  !=================================================================================================
  subroutine fv_phys_to_dyn(elem,T_tmp,uv_tmp,q_tmp)
    ! Purpose: Copy physics state to dynamics grid
    use control_mod,    only: ftype
    use dyn_comp,       only: TimeLevel, hvcoord
    use derivative_mod, only: subcell_integration
    implicit none
    !---------------------------------------------------------------------------
    ! interface arguments
    type(element_t), intent(inout) :: elem(:)         ! dynamics element structure
    real(r8),        intent(inout) :: T_tmp (:,:,:)   ! temp array to hold T
    real(r8),        intent(inout) :: uv_tmp(:,:,:,:) ! temp array to hold u and v
    real(r8),        intent(inout) :: q_tmp (:,:,:,:) ! temp to hold advected constituents
    ! local variables
    integer(i4) :: ie, m, i, j, icol, ilyr            ! loop iterators
    integer(i4) :: ii, jj, gi, gj                     ! GLL loop iterator and indices for pg2
    integer     :: tl_f
    real(r8), dimension(fv_nphys*fv_nphys,pver,pcnst) :: qo_phys       ! reconstructed initial physics state 
    real(r8), dimension(np,np)                        :: tmp_area_gll  ! area for weighting
    real(r8), dimension(fv_nphys*fv_nphys)            :: inv_area_fvm  ! inverse area for weighting
    real(r8), dimension(np,np,pver)                   :: dp_gll
    real(r8), dimension(fv_nphys*fv_nphys,pver)       :: dp_fvm_sum
    !---------------------------------------------------------------------------
    ! Copy tendencies on the physics grid over to the dynamics grid (GLL)
    !---------------------------------------------------------------------------
    tl_f = TimeLevel%n0
    tmp_area_gll(:,:) = 1.0_r8
    do ie = 1,nelemd
      !-------------------------------------------------------------------------
      ! Recalculate state that was previously sent to physics 
      !-------------------------------------------------------------------------
      if (ftype==2.or.ftype==4) then
        inv_area_fvm(:) = 1.0_r8/RESHAPE( subcell_integration(tmp_area_gll,np,fv_nphys,elem(ie)%metdet(:,:)), (/fv_nphys*fv_nphys/) )
        do ilyr = 1,pver
          dp_gll(:,:,ilyr) = elem(ie)%state%dp3d(:,:,ilyr,tl_f)
          dp_fvm_sum(:,ilyr) = RESHAPE( subcell_integration(dp_gll(:,:,ilyr),np,fv_nphys,elem(ie)%metdet(:,:)), (/fv_nphys*fv_nphys/) )
          do m = 1,pcnst
            qo_phys(:,ilyr,m)  = RESHAPE( subcell_integration(                    &
                                  elem(ie)%state%Q(:,:,ilyr,m)*dp_gll(:,:,ilyr),  &
                                  np, fv_nphys, elem(ie)%metdet(:,:) ),           &
                                  (/fv_nphys*fv_nphys/) ) / dp_fvm_sum(:,ilyr)
          end do ! m
        end do ! ilyr
      end if 
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      icol = 0
      do j = 1,fv_nphys
        do i = 1,fv_nphys 
          icol = icol + 1
          do ilyr = 1,pver
            !-------------------------------------------------------------------
            !-------------------------------------------------------------------
            ! the pg1 case is simple, just copy to all GLL nodes in the element
            if (fv_nphys == 1) then
              elem(ie)%derived%FT(:,:,  ilyr)   = T_tmp (icol,  ilyr,ie)
              elem(ie)%derived%FM(:,:,1,ilyr)   = uv_tmp(icol,1,ilyr,ie)
              elem(ie)%derived%FM(:,:,2,ilyr)   = uv_tmp(icol,2,ilyr,ie)
              do m = 1,pcnst
                if ( ftype==2 .or. ftype==4 ) then
                  ! subtract initial phys state and add previous dyn state
                  elem(ie)%derived%FQ(:,:,ilyr,m) = ( q_tmp(icol,ilyr,m,ie)   &
                                                     -qo_phys(icol,ilyr,m) )  &
                                                    *dp_fvm_sum(icol,ilyr)        &
                                                    /(4.0_r8*dp_gll(gi,gj,ilyr)   &
                                                      *elem(ie)%spheremp(gi,gj) ) &
                                                    +elem(ie)%state%q(:,:,ilyr,m)
                else
                  elem(ie)%derived%FQ(:,:,ilyr,m) = q_tmp(icol,ilyr,m,ie)
                end if
              end do
            end if ! fv_nphys == 1
            !-------------------------------------------------------------------
            !-------------------------------------------------------------------
            ! for pg2 we need to copy the FV physics state to quadrants of GLL grid
            if (fv_nphys == 2) then
              do jj = 1,2
                do ii = 1,2
                  if (i==1) gi = ii
                  if (j==1) gj = jj
                  if (i==2) gi = ii+2
                  if (j==2) gj = jj+2
                  elem(ie)%derived%FT(gi,gj,  ilyr) =  T_tmp(icol,  ilyr,ie)
                  elem(ie)%derived%FM(gi,gj,1,ilyr) = uv_tmp(icol,1,ilyr,ie)
                  elem(ie)%derived%FM(gi,gj,2,ilyr) = uv_tmp(icol,2,ilyr,ie)
                  do m = 1,pcnst
                    if ( ftype==2 .or. ftype==4 ) then
                      ! subtract initial phys state and add previous dyn state
                      elem(ie)%derived%FQ(gi,gj,ilyr,m) = ( q_tmp(icol,ilyr,m,ie)       &
                                                           -qo_phys(icol,ilyr,m) )      &
                                                          *dp_fvm_sum(icol,ilyr)        &
                                                          /(4.0_r8*dp_gll(gi,gj,ilyr)   &
                                                            *elem(ie)%spheremp(gi,gj) ) &
                                                          +elem(ie)%state%q(gi,gj,ilyr,m)
                    else
                      elem(ie)%derived%FQ(gi,gj,ilyr,m) = q_tmp(icol,ilyr,m,ie)
                    end if
                  end do
                end do ! ii
              end do ! jj
            end if ! fv_nphys == 2
            !-------------------------------------------------------------------
            !-------------------------------------------------------------------
          end do ! ilyr
        end do ! i
      end do ! j
    end do ! ie
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

  end subroutine fv_phys_to_dyn
  !=================================================================================================
  !=================================================================================================
  subroutine fv_phys_to_dyn_topo(elem,phis_tmp)
    ! Purpose: topo is initially defined on phys grid, 
    !          so this routine copys it to the dynamics grid
    implicit none
    !---------------------------------------------------------------------------
    ! interface arguments
    type(element_t), intent(inout) :: elem(:)        ! dynamics element structure
    real(r8),        intent(inout) :: phis_tmp(:,:)  ! temp array to hold PHIS field from file
    ! local variables
    integer(i4) :: ie, i, j, icol  ! loop iterators
    integer(i4) :: ii, jj, gi, gj  ! GLL loop iterator and indices for pg2
    !---------------------------------------------------------------------------
    ! Copy topography on the physics grid over to the dynamics grid (GLL)
    !---------------------------------------------------------------------------
    do ie = 1,nelemd
      icol = 0
      do j = 1,fv_nphys
        do i = 1,fv_nphys 
          icol = icol + 1
          !-------------------------------------------------------------------
          ! Store topo data in fv_physgrid to avoid mapping back and forth
          !-------------------------------------------------------------------
          fv_physgrid(ie)%topo(i,j) = phis_tmp(icol,ie)
          !-------------------------------------------------------------------
          !-------------------------------------------------------------------
          ! pg1 case 
          if (fv_nphys == 1) then
            elem(ie)%state%phis(:,:) = phis_tmp(icol,ie)
          end if ! fv_nphys == 1
          !-------------------------------------------------------------------
          !-------------------------------------------------------------------
          ! for pg2 we need to copy the FV state to quadrants of GLL grid
          if (fv_nphys == 2) then
            do jj = 1,2
              do ii = 1,2
                if (i==1) gi = ii
                if (j==1) gj = jj
                if (i==2) gi = ii+2
                if (j==2) gj = jj+2
                elem(ie)%state%phis(gi,gj) = phis_tmp(icol,ie)
              end do ! ii
            end do ! jj
          end if ! fv_nphys == 2
          !-------------------------------------------------------------------
          !-------------------------------------------------------------------
        end do ! i
      end do ! j
      ! Weight topo field for boundary exchange in read_inidat()
      elem(ie)%state%phis(:,:) = elem(ie)%state%phis(:,:) &
                                *elem(ie)%spheremp(:,:)   &
                                *elem(ie)%rspheremp(:,:)
    end do ! ie
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
  end subroutine fv_phys_to_dyn_topo
  !=================================================================================================
  !=================================================================================================
  subroutine dyn_to_fv_phys(elem,ps_tmp,zs_tmp,T_tmp,uv_tmp,om_tmp,Q_tmp)
    ! Purpose: average dynamics state over subcells and assign to physics state
    use derivative_mod,     only: subcell_integration
    use dyn_comp,           only: TimeLevel, hvcoord
    use element_ops,        only: get_temperature
    use time_manager,       only: is_first_step
    implicit none
    !---------------------------------------------------------------------------
    ! interface arguments
    type(element_t), intent(inout) :: elem(:)          ! dynamics element structure
    real(r8),        intent(inout) :: ps_tmp(:,:)      ! temp array to hold ps
    real(r8),        intent(inout) :: zs_tmp(:,:)      ! temp array to hold phis  
    real(r8),        intent(inout) :: T_tmp (:,:,:)    ! temp array to hold T
    real(r8),        intent(inout) :: uv_tmp(:,:,:,:)  ! temp array to hold u and v
    real(r8),        intent(inout) :: om_tmp(:,:,:)    ! temp array to hold omega
    real(r8),        intent(inout) :: Q_tmp (:,:,:,:)  ! temp to hold advected constituents
    ! local variables
    integer(i4) :: ie, m, icol, ilyr             ! loop iterators
    integer     :: tl_f                          ! time level
    integer     :: ncol
    real(r8), dimension(np,np,nlev)        :: temperature   ! Temperature from dynamics
    real(r8), dimension(np,np)             :: tmp_area      ! area for weighting
    real(r8), dimension(fv_nphys,fv_nphys) :: inv_area      ! inverse area for weighting
    real(r8), dimension(np,np)             :: dp_gll        ! pressure thickness on GLL grid
    real(r8), dimension(np,np)             :: dp_gll_in     ! pressure thickness on GLL grid
    real(r8), dimension(fv_nphys,fv_nphys) :: inv_dp_fvm    ! inverted pressure thickness on FV grid
    real(r8), dimension(fv_nphys,fv_nphys) :: inv_dp_fvm_in ! inverted pressure thickness on FV grid
    real(r8), dimension(npsq)              :: T_tmp_in      ! temp array to hold previous dyn state T 
    real(r8), dimension(npsq,pcnst)        :: Q_tmp_in      ! temp array to hold previous dyn state q 
    real(r8), dimension(npsq,2)            :: uv_tmp_in     ! temp array to hold previous dyn state uv
    !---------------------------------------------------------------------------
    ! Integrate dynamics field with appropriate weighting 
    ! to get average state in each physics cell
    !---------------------------------------------------------------------------
    tl_f = TimeLevel%n0
    ncol = fv_nphys*fv_nphys
    tmp_area(:,:) = 1.0_r8

    do ie = 1,nelemd

      inv_area(:,:) = 1.0_r8/subcell_integration(tmp_area,np,fv_nphys,elem(ie)%metdet(:,:))

      ps_tmp(:,ie) = RESHAPE( subcell_integration(                  &
                     elem(ie)%state%ps_v(:,:,tl_f),                 &
                     np, fv_nphys, elem(ie)%metdet(:,:) )           &
                     *inv_area , (/ncol/) )

      zs_tmp(:,ie) = RESHAPE( fv_physgrid(ie)%topo(:,:), (/ncol/) )

      call get_temperature(elem(ie),temperature,hvcoord,tl_f)

      do ilyr = 1,pver

        dp_gll(:,:) = elem(ie)%state%dp3d(:,:,ilyr,tl_f)
        inv_dp_fvm = 1.0 / subcell_integration(dp_gll,np,fv_nphys,elem(ie)%metdet(:,:))
        
        T_tmp(:ncol,ilyr,ie)      = RESHAPE( subcell_integration(             &
                                    temperature(:,:,ilyr)*dp_gll,             &
                                    np, fv_nphys, elem(ie)%metdet(:,:) )      &
                                    *inv_dp_fvm, (/ncol/) )

        om_tmp(:ncol,ilyr,ie)     = RESHAPE( subcell_integration(             &
                                    elem(ie)%derived%omega_p(:,:,ilyr),       &
                                    np, fv_nphys, elem(ie)%metdet(:,:) )      &
                                    *inv_area , (/ncol/) )
        do m = 1,2
          uv_tmp(:ncol,m,ilyr,ie) = RESHAPE( subcell_integration(             &
                                    elem(ie)%state%V(:,:,m,ilyr,tl_f),        &
                                    np, fv_nphys, elem(ie)%metdet(:,:) )      &
                                    *inv_area , (/ncol/) )
        end do
        do m = 1,pcnst
          Q_tmp(:ncol,ilyr,m,ie)  = RESHAPE( subcell_integration(             &
                                    elem(ie)%state%Q(:,:,ilyr,m)*dp_gll,      &
                                    np, fv_nphys, elem(ie)%metdet(:,:) )      &
                                    *inv_dp_fvm, (/ncol/) )
        end do

        !-----------------------------------------------------------------------
        ! Map previous dynamics state to physgrid for CRM forcing
        !-----------------------------------------------------------------------
        if (.not.is_first_step()) then

          dp_gll_in(:,:) = elem(ie)%state%dp_in(:,:,ilyr)
          inv_dp_fvm_in = 1.0 / subcell_integration(dp_gll_in,np,fv_nphys,elem(ie)%metdet(:,:))

          T_tmp_in(:ncol)         = RESHAPE( subcell_integration(             &
                                    elem(ie)%state%T_in(:,:,ilyr)*dp_gll_in,  &
                                    np, fv_nphys, elem(ie)%metdet(:,:) )      &
                                    *inv_dp_fvm_in, (/ncol/) )

          do m = 1,pcnst
            Q_tmp_in(:ncol,m)     = RESHAPE( subcell_integration(             &
                                    elem(ie)%state%Q_in(:,:,ilyr,m)*dp_gll,   &
                                    np, fv_nphys, elem(ie)%metdet(:,:) )      &
                                    *inv_dp_fvm_in, (/ncol/) )
          end do

          do m = 1,2
            uv_tmp_in(:ncol,m)    = RESHAPE( subcell_integration(             &
                                    elem(ie)%state%V_in(:,:,m,ilyr),          &
                                    np, fv_nphys, elem(ie)%metdet(:,:) )      &
                                    *inv_area , (/ncol/) )
          end do

          ! Calculate tendency from mapped states
          T_tmp (:ncol,ilyr,ie)   =  T_tmp(:ncol,ilyr,ie)   -  T_tmp_in(:ncol)
          Q_tmp (:ncol,ilyr,:,ie) =  Q_tmp(:ncol,ilyr,:,ie) -  Q_tmp_in(:ncol,:)
          uv_tmp(:ncol,:,ilyr,ie) = uv_tmp(:ncol,:,ilyr,ie) - uv_tmp_in(:ncol,:)

        end if ! not is_first_step
        !-----------------------------------------------------------------------
        !-----------------------------------------------------------------------

      end do ! ilyr
    end do ! ie
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
  end subroutine dyn_to_fv_phys
  !=================================================================================================
  !=================================================================================================
end module fv_physics_coupling_mod