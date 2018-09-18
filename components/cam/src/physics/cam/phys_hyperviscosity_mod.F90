module phys_hyperviscosity_mod

#if defined( DIFFUSE_PHYS_TEND )

   use kinds,           only: real_kind
   use hybrid_mod,      only: hybrid_t
   use element_mod,     only: element_t
   use derivative_mod,  only: derivative_t
   use time_mod,        only: timelevel_t
   use hybvcoord_mod,   only: hvcoord_t
   use domain_mod,      only: domain1d_t
   use edgetype_mod,    only: EdgeBuffer_t
   use bndry_mod,       only: bndry_exchangeV

   implicit none

   public :: phys_hyperviscosity
   
   contains

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------

subroutine phys_hyperviscosity(ptend)
   !----------------------------------------------
   ! Main interface for hyperviscosity calculation on state after being modifed by physics.
   ! This mimics what happens in dynamics hypervis calculation, but is not identical.
   ! Original motivation was to effectively smooth the output tendencies from the CRM.
   ! Author: Walter Hannah (LLNL)
   !----------------------------------------------
   use parallel_mod,    only: par
   use time_mod,        only: tstep
   use hybrid_mod,      only: hybrid_create
   use dyn_comp,        only: hvcoord, TimeLevel
   use thread_mod,      only: hthreads, omp_get_thread_num
   use prim_driver_base,only: prim_init1
   use phys_grid,       only: get_ncols_p, get_gcol_all_p
   use dimensions_mod,  only: np, npsq, nelemd, nlev
   use ppgrid,          only: begchunk, endchunk, pcols, pver, pverp
   use dof_mod,         only: putUniquePoints
   use physics_types,   only: physics_ptend 
   !!! Interface arguments
   type(physics_ptend), intent(inout) :: ptend(begchunk:endchunk)
   !!! Local variables
   integer :: nt
   integer :: ithr
   integer :: icol, ncols, lchnk
   integer :: nets, nete
   integer :: ie, ioff, k
   integer :: pgcols(pcols), idmb1(1), idmb2(1), idmb3(1)
   type (element_t),  pointer :: phys_elem(:)        ! dummy element structure to hold tendencies from physics
   type (domain1d_t), pointer :: dom_mt(:) => null()
   type (hybrid_t)            :: hybrid
   type (derivative_t)        :: deriv

   real (kind=real_kind), dimension(npsq,pver,nelemd) :: T_tmp       ! temporary array to physics tend
   
   !----------------------------------------------
   ! setup element and TimeLevel structure
   !----------------------------------------------   
   call f(phys_elem,par,dom_mt,TimeLevel)

   nets = dom_mt(ithr)%start
   nete = dom_mt(ithr)%end

   nt = TimeLevel%n0

   !----------------------------------------------
   ! initialize hybrid and derivative structure
   !----------------------------------------------
   ithr = omp_get_thread_num()
   hybrid = hybrid_create(par,ithr,hthreads)

   call derivinit(deriv)

   !----------------------------------------------
   ! copy physics tendencies to element structure 
   ! based on cam/src/dynamics/se/dp_coupling.F90
   !----------------------------------------------
   do lchnk = begchunk,endchunk
      ncols = get_ncols_p(lchnk)
      call get_gcol_all_p(lchnk,pcols,pgcols)
      do icol = 1,ncols
         call get_gcol_block_d(pgcols(icol),1,idmb1,idmb2,idmb3)
         ie   = idmb3(1)
         ioff = idmb2(1)
         do k = 1,pver
            T_tmp(ioff,k,ie) = ptend(lchnk)%s(icol,k)
         end do ! k
      end do ! icol
   end do ! lchnk

   ! if (par%dynproc) then
      do ie = 1,nelemd
         ncols = phys_elem(ie)%idxP%NumUniquePts
         ! call putUniquePoints( elem(ie)%idxP, nlev, T_tmp(1:ncols,:,ie), elem(ie)%derived%fT(:,:,:) )
         call putUniquePoints( phys_elem(ie)%idxP, nlev, T_tmp(1:ncols,:,ie), phys_elem(ie)%state%T(:,:,:,nt) )
      end do
   ! end if ! par%dynproc
   
   !----------------------------------------------
   ! diffuse temperature tendency
   !----------------------------------------------
   call phys_hyperviscosity_T( phys_elem, hvcoord, hybrid, deriv, nt, nets, nete, dt_in=tstep )

   !----------------------------------------------
   ! deallocate element structure
   !----------------------------------------------
   deallocate(phys_elem)

   !----------------------------------------------
   !----------------------------------------------
end subroutine phys_hyperviscosity

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------

subroutine phys_hyperviscosity_T( elem, hvcoord, hybrid, deriv, nt, nets, nete, dt_in )
   !----------------------------------------------
   ! hyperviscsoity operator for foward-in-time scheme
   ! take one timestep of: T(:,:,:,np) = T(:,:,:,np) +  dt2*nu_s*laplacian**order ( T )
   ! For correct scaling, dt2 should be the same 'dt2' used in the leapfrog scheme
   ! based on :
   !     advance_hypervis_dp()      from homme/src/preqx/share/prim_advance_mod.F90
   !     advance_hypervis_scalar()  from homme/src/share/prim_advection_base.F90
   ! Author: Walter Hannah (LLNL)
   !----------------------------------------------
   use dimensions_mod,  only: np, nlev
   use edge_mod,        only: edgevpack_nlyr
   use bndry_mod,       only: bndry_exchangev
   use control_mod,     only: hypervis_order, nu, nu_q, nu_s, nu_p, hypervis_subcycle
   !!! Interface arguments
   type (element_t)     , intent(inout) :: elem(:)
   type (hvcoord_t)     , intent(in   ) :: hvcoord
   type (hybrid_t)      , intent(in   ) :: hybrid
   type (derivative_t)  , intent(in   ) :: deriv
   integer              , intent(in   ) :: nt
   integer              , intent(in   ) :: nets
   integer              , intent(in   ) :: nete
   real (kind=real_kind), intent(in   ) :: dt_in
   !!! Local variables
   real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: tensor 
   real (kind=real_kind), dimension(np,np,nlev)           :: dp
   type (EdgeBuffer_t)   :: edge_buffer
   real (kind=real_kind) :: dt
   integer :: k, i, j, ie, ic, q

   !----------------------------------------------
   ! Exit conditions
   !----------------------------------------------
   if ( nu_q == 0 ) return
   if ( hypervis_order /= 2 ) return
   if ( nu_s==0 .and. nu==0 .and. nu_p==0 ) return

   !----------------------------------------------
   ! hyper viscosity calculation
   !----------------------------------------------
   ! nu_p = 0 :
   !   scale T dissipaton by dp  (conserve IE, dissipate T^2)
   ! nu_p > 0 :
   !   dont scale: T equation IE dissipation matches (to truncation error)
   !               IE dissipation from continuity equation (1 deg: to about 0.1 W/m^2)

   dt = dt_in / hypervis_subcycle

   do ic = 1, hypervis_subcycle

      !------------------------------------
      ! compute laplacian
      !------------------------------------
      !!! Note: tensor = input and output
      ! call biharmonic_wk_scalar( elem, tensor, deriv, edge_g, hybrid, nets, nete )

      call phys_biharmonic_wk_scalar(hybrid, deriv, nets, nete, elem, edge_buffer, tensor)

      !------------------------------------
      ! Apply hyperviscosity
      !------------------------------------
      do ie = nets, nete
         do k = 1, nlev
            do j = 1, np
               do i = 1, np
                  ! elem(ie)%state%T(i,j,k,nt) = elem(ie)%state%T(i,j,k,nt) * elem(ie)%spheremp(i,j) - dt * nu_s * tensor(i,j)
                  elem(ie)%state%T(i,j,k,nt) = elem(ie)%state%T(i,j,k,nt) - dt * nu_s * tensor(i,j,k,ie) * elem(ie)%rspheremp(i,j)
               enddo ! i = 1, np
            enddo ! j = 1, np
         enddo ! k = 1, nlev

         !!! prepare edge buffer
         call edgeVpack_nlyr(edge_buffer, elem(ie)%desc, elem(ie)%state%T(:,:,:,nt), nlev, 0, nlev )
      enddo ! ie = nets, nete

      !------------------------------------
      ! update boundaries
      !------------------------------------
      call bndry_exchangeV(hybrid, edge_buffer )

      !------------------------------------
      !------------------------------------
      do ie = nets, nete
         call edgeVunpack_nlyr(edge_buffer, elem(ie)%desc, elem(ie)%state%T(:,:,:,nt), nlev, 0, nlev)
         
         !!! apply inverse mass matrix
         ! do k = 1, nlev
         !    elem(ie)%state%T(:,:,k,nt) = elem(ie)%rspheremp(:,:) * elem(ie)%state%T(:,:,k,q,nt)
         ! enddo ! k
         
      enddo ! ie = nets, nete
      !------------------------------------
      !------------------------------------

   enddo ! ic = 1, hypervis_subcycle
   !----------------------------------------------
   !----------------------------------------------
  
end subroutine phys_hyperviscosity_T

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------

! subroutine physics_hypervis_q( elem, hvcoord, hybrid, deriv, nt_qdp, nets, nete, dt2 )
!   ! hyperviscsoity operator for foward-in-time scheme
!   ! take one timestep of:
!   !         Q(:,:,:,np) = Q(:,:,:,np) +  dt2*nu*laplacian**order ( Q )
!   ! For correct scaling, dt2 should be the same 'dt2' used in the leapfrog scheme
!   use kinds          , only : real_kind
!   use dimensions_mod , only : np, nlev
!   use hybrid_mod     , only : hybrid_t
!   use element_mod    , only : element_t
!   use derivative_mod , only : derivative_t
!   use bndry_mod      , only : bndry_exchangev
!   use viscosity_mod,   only: biharmonic_wk_scalar
!   implicit none
!   type (element_t)     , intent(inout), target :: elem(:)
!   type (hvcoord_t)     , intent(in   )         :: hvcoord
!   type (hybrid_t)      , intent(in   )         :: hybrid
!   type (derivative_t)  , intent(in   )         :: deriv
!   integer              , intent(in   )         :: nt_qdp
!   integer              , intent(in   )         :: nets
!   integer              , intent(in   )         :: nete
!   real (kind=real_kind), intent(in   )         :: dt2

!   ! local
!   real (kind=real_kind), dimension(np,np,nlev,qsize,nets:nete) :: Qtens 
!   real (kind=real_kind), dimension(np,np,nlev                ) :: dp
!   real (kind=real_kind) :: dt
!   integer :: k, i, j, ie, ic, q

!   if ( nu_q           == 0 ) return
!   if ( hypervis_order /= 2 ) return

!   !----------------------------------------------
!   ! hyper viscosity calculation
!   !----------------------------------------------

!   ! Note: Qtens = Q/dp   (apply hyperviscsoity to (dp0*Qdp/dp), not Qdp)
!   ! various options:
!   !   1)  biharmonic( Qdp )
!   !   2)  dp0 * biharmonic( Qdp/dp )
!   !   3)  dpave * biharmonic(Q/dp)      (where dpave is avg MF from nu_p contribution from dynamics)
!   ! For trace mass / mass consistenciy, we use #2 when nu_p=0 and #3 when nu_p>0, 

!   dt = dt2 / hypervis_subcycle_q

!   do ic = 1, hypervis_subcycle_q

!     !------------------------------------
!     ! Calculate Qtens
!     !------------------------------------
!     do ie = nets, nete
!       if (nu_p>0) then
!         !!! use option #3
!         do q = 1, qsize
!           do k = 1, nlev
!             dp(:,:,k) = elem(ie)%derived%dp(:,:,k) - dt2*elem(ie)%derived%divdp_proj(:,:,k)
!             Qtens(:,:,k,q,ie) = elem(ie)%derived%dpdiss_ave(:,:,k)*&
!                                 elem(ie)%state%Qdp(:,:,k,q,nt_qdp) / dp(:,:,k)
!           enddo ! k = 1, nlev
!         enddo ! q = 1, qsize

!       else
!         !!! use option #2
!         do q = 1, qsize
!           do k = 1, nlev
!             dp(:,:,k) = elem(ie)%derived%dp(:,:,k) - dt2 * elem(ie)%derived%divdp_proj(:,:,k)
!             Qtens(:,:,k,q,ie) = hvcoord%dp0(k) * elem(ie)%state%Qdp(:,:,k,q,nt_qdp) / dp(:,:,k)
!           enddo ! k = 1, nlev
!         enddo ! q = 1, qsize

!       endif ! nu_p > 0
!     enddo ! ie = nets, nete

!     !------------------------------------
!     ! compute laplacian
!     !------------------------------------
!     !!! Note: Qtens = input and output
!     call biharmonic_wk_scalar( elem, Qtens, deriv, edge_g, hybrid, nets, nete )

!     !------------------------------------
!     ! Apply hyperviscosity
!     !------------------------------------
!     do ie = nets, nete
!       do q = 1, qsize
!         do k = 1, nlev
!           do j = 1, np
!             do i = 1, np
!               ! advection Qdp.  For mass advection consistency:
!               ! DIFF( Qdp) ~   dp0 DIFF (Q)  =  dp0 DIFF ( Qdp/dp )
!               elem(ie)%state%Qdp(i,j,k,q,nt_qdp) = elem(ie)%state%Qdp(i,j,k,q,nt_qdp) * elem(ie)%spheremp(i,j) - dt * nu_q * Qtens(i,j,k,q,ie)
!             enddo ! i = 1, np
!           enddo ! j = 1, np
!         enddo ! k = 1, nlev
        
!         if (limiter_option .ne. 0 ) then
!            !!! smooth the negativities introduced by diffusion:
!            call limiter2d_zero( elem(ie)%state%Qdp(:,:,:,q,nt_qdp) )
!         endif

!       enddo ! q = 1, qsize
!       !!! prepare edge buffer
!       call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%Qdp(:,:,:,:,nt_qdp), qsize*nlev, 0, qsize*nlev )
!     enddo ! ie = nets, nete

!     !------------------------------------
!     ! update boundaries
!     !------------------------------------
!     call bndry_exchangeV( hybrid, edge_g )

!     !------------------------------------
!     !------------------------------------
!     do ie = nets, nete
!       call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%Qdp(:,:,:,:,nt_qdp), qsize*nlev, 0, qsize*nlev)
!       do q = 1, qsize
!         !!! apply inverse mass matrix
!         do k = 1, nlev
!           elem(ie)%state%Qdp(:,:,k,q,nt_qdp) = elem(ie)%rspheremp(:,:) * elem(ie)%state%Qdp(:,:,k,q,nt_qdp)
!         enddo ! k
!       enddo ! q
!     enddo ! ie = nets, nete
!     !------------------------------------
!     !------------------------------------
    
!   enddo ! ic = 1, hypervis_subcycle_q
!   !----------------------------------------------
!   !----------------------------------------------
  
! end subroutine physics_hypervis_q

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------

subroutine phys_biharmonic_wk_scalar(hybrid, deriv, nets, nete, elem, edge_buffer, qtens)
   !----------------------------------------------
   ! compute weak biharmonic operator (for physics)
   !    input:  qtens = Q
   !    output: qtens = weak biharmonic of Q
   ! Adapted from homme/src/share/viscosity_base.F90
   ! Author: Walter Hannah (LLNL)
   !----------------------------------------------
   use dimensions_mod,      only: np, nlev
   use control_mod,         only: hypervis_scaling
   use derivative_mod_base, only: laplace_sphere_wk
   !!! Interface arguments
   type (hybrid_t)      , intent(in   ) :: hybrid
   type (derivative_t)  , intent(in   ) :: deriv
   integer              , intent(in   ) :: nets
   integer              , intent(in   ) :: nete
   type (element_t)     , intent(inout) :: elem(:)
   type (EdgeBuffer_t)  , intent(inout) :: edge_buffer
   real (kind=real_kind), intent(inout),  dimension(np,np,nlev,nets:nete) :: qtens

   !!! Local Variables
   integer :: i,j,k,ie
   real (kind=real_kind), dimension(np,np) :: lap_p
   logical var_coef1

   !!! if tensor hyperviscosity with tensor V is used, 
   !!! then biharmonic operator is (\grad\cdot V\grad) (\grad \cdot \grad) 
   !!! so tensor is only used on second call to phys_laplace_sphere_wk
   var_coef1 = .true.
   if(hypervis_scaling > 0) var_coef1 = .false.

   !----------------------------------------------
   ! Calculate laplacian tensor
   !----------------------------------------------
   do ie=nets,nete
      do k=1,nlev
         lap_p(:,:) = qtens(:,:,k,ie)
         ! qtens(:,:,k,ie) = phys_laplace_sphere_wk(lap_p, deriv, elem(ie), var_coef=var_coef1)
         qtens(:,:,k,ie) = laplace_sphere_wk(lap_p, deriv, elem(ie), var_coef=var_coef1)
      enddo
      !!! prepare edge buffer
      call edgeVpack_nlyr(edge_buffer, elem(ie)%desc, qtens(:,:,:,ie),nlev, 0, nlev)
   enddo

   !----------------------------------------------
   ! update boundaries
   !----------------------------------------------
   call bndry_exchangeV(hybrid,edge_buffer)
   
   !----------------------------------------------
   ! apply inverse mass matrix, then apply laplace again
   !----------------------------------------------
   do ie=nets,nete
      call edgeVunpack_nlyr(edge_buffer, elem(ie)%desc, qtens(:,:,:,ie), nlev, 0, nlev)
      do k=1,nlev 
         lap_p(:,:) = elem(ie)%rspheremp(:,:)*qtens(:,:,k,ie)
         ! qtens(:,:,k,ie) = phys_laplace_sphere_wk(lap_p, deriv, elem(ie), var_coef=.true.)
         qtens(:,:,k,ie) = laplace_sphere_wk(lap_p, deriv, elem(ie), var_coef=.true.)
      enddo
   enddo

  !----------------------------------------------
  !----------------------------------------------
end subroutine phys_biharmonic_wk_scalar


!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------

! function phys_laplace_sphere_wk(s, deriv, elem, var_coef) result(laplace)
!    !----------------------------------------------
!    ! Calculate laplacian on sphere
!    !   input:  s = scalar
!    !   ouput:  -< grad(PHI), grad(s) >   = weak divergence of grad(s)
!    ! note: for this form of the operator, grad(s) does not need to be made C0 
!    ! Adapted from homme/src/share/derivative_mod_base.F90
!    ! Author: Walter Hannah (LLNL)
!    !----------------------------------------------
!    use derivative_mod_base, only: divergence_sphere_wk
!    !!! Interface arguments
!    real (kind=real_kind), intent(in)       :: s(np,np) 
!    type (derivative_t),   intent(in)       :: deriv
!    type (element_t),      intent(in)       :: elem
!    logical,               intent(in)       :: var_coef
!    real (kind=real_kind), dimension(np,np) :: laplace
!    !!! Local variables
!    integer i,j
!    real(kind=real_kind), dimension(np,np,2) :: grad_sph_1
!    real(kind=real_kind), dimension(np,np,2) :: grad_sph_2
!    !----------------------------------------------
!    !----------------------------------------------
!    grad_sph_1 = gradient_sphere(s,deriv,elem%Dinv)

!    if (var_coef) then
!       if (hypervis_power/=0) then
!          !!! scalar viscosity with variable coefficient
!          grad_sph_1(:,:,1) = grad_sph_1(:,:,1) * elem%variable_hyperviscosity(:,:)
!          grad_sph_1(:,:,2) = grad_sph_1(:,:,2) * elem%variable_hyperviscosity(:,:)
!       else if (hypervis_scaling /=0 ) then
!          !!! tensor hv, (3)
!          grad_sph_2 = grad_sph_1
!          do j = 1,np
!             do i = 1,np
!                grad_sph_1(i,j,1) = grad_sph_2(i,j,1) * elem%tensorVisc(i,j,1,1) + &
!                                    grad_sph_2(i,j,2) * elem%tensorVisc(i,j,1,2)
!                grad_sph_1(i,j,2) = grad_sph_2(i,j,1) * elem%tensorVisc(i,j,2,1) + &
!                                    grad_sph_2(i,j,2) * elem%tensorVisc(i,j,2,2)
!             end do
!          end do
!       endif ! hypervis_power/=0
!    endif ! var_coef
!    !----------------------------------------------
!    !----------------------------------------------
!    !!! divergnece_sphere and divergence_sphere_wk are identical *after* bndry_exchange
!    !!! if input is C_0.  Here input is not C_0, so we should use divergence_sphere_wk().
!    laplace = divergence_sphere_wk(grad_sph_1,deriv,elem)
!    !----------------------------------------------
!    !----------------------------------------------
! end function phys_laplace_sphere_wk

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------

! subroutine biharmonic_wk_dp3d(elem,dptens,ptens,vtens,deriv,edge3,hybrid,nt,nets,nete)
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ! compute weak biharmonic operator
!    !    input:  h,v (stored in elem()%, in lat-lon coordinates
!    !    output: ptens,vtens  overwritten with weak biharmonic of h,v (output in lat-lon coordinates)
!    !
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    type (hybrid_t)      , intent(in) :: hybrid
!    type (element_t)     , intent(inout), target :: elem(:)
!    integer              , intent(in)  :: nt,nets,nete
!    real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)  :: vtens
!    real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: ptens,dptens
!    type (EdgeBuffer_t)  , intent(inout) :: edge3
!    type (derivative_t)  , intent(in) :: deriv

!    ! local
!    integer :: i,j,k,kptr,ie
!    real (kind=real_kind), dimension(:,:), pointer :: rspheremv
!    real (kind=real_kind), dimension(np,np) :: tmp
!    real (kind=real_kind), dimension(np,np) :: tmp2
!    real (kind=real_kind), dimension(np,np,2) :: v
!    real (kind=real_kind) :: nu_ratio1, nu_ratio2
!    logical var_coef1

!    !if tensor hyperviscosity with tensor V is used, then biharmonic operator is (\grad\cdot V\grad) (\grad \cdot \grad) 
!    !so tensor is only used on second call to laplace_sphere_wk
!    var_coef1 = .true.
!    if(hypervis_scaling > 0)    var_coef1 = .false.

!    ! note: there is a scaling bug in the treatment of nu_div
!    ! nu_ratio is applied twice, once in each laplace operator
!    ! so in reality:   nu_div_actual = (nu_div/nu)**2 nu
!    ! We should fix this, but it requires adjusting all CAM defaults
!    nu_ratio1=1
!    nu_ratio2=1
!    if (nu_div/=nu) then
!       if(hypervis_scaling /= 0) then
!          ! we have a problem with the tensor in that we cant seperate
!          ! div and curl components.  So we do, with tensor V:
!          ! nu * (del V del ) * ( nu_ratio * grad(div) - curl(curl))
!          nu_ratio1=(nu_div/nu)**2   ! preserve buggy scaling
!          nu_ratio2=1
!       else
!          nu_ratio1=nu_div/nu
!          nu_ratio2=nu_div/nu
!       endif
!    endif


!    do ie=nets,nete

!       do k=1,nlev
!          tmp=elem(ie)%state%T(:,:,k,nt) 
!          ptens(:,:,k,ie)=laplace_sphere_wk(tmp,deriv,elem(ie),var_coef=var_coef1)
!       enddo
!       kptr=0
!       call edgeVpack_nlyr(edge3, elem(ie)%desc, ptens(1,1,1,ie),nlev,kptr,4*nlev)

!    enddo
   
!    call bndry_exchangeV(hybrid,edge3)
   
!    do ie=nets,nete
!       rspheremv     => elem(ie)%rspheremp(:,:)
      
!       kptr=0
!       call edgeVunpack_nlyr(edge3, elem(ie)%desc, ptens(1,1,1,ie), nlev, kptr, 4*nlev)
      
      
!       ! apply inverse mass matrix, then apply laplace again
!       do k=1,nlev
!          tmp(:,:)=rspheremv(:,:)*ptens(:,:,k,ie)
!          ptens(:,:,k,ie)=laplace_sphere_wk(tmp,deriv,elem(ie),var_coef=.true.)
!       enddo
!    enddo
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! end subroutine

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------

#endif /* DIFFUSE_PHYS_TEND */

end module phys_hyperviscosity_mod