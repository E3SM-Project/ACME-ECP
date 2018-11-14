module phys_hyperviscosity_mod

#if defined( DIFFUSE_PHYS_TEND )

   use shr_kind_mod,    only: r8 => shr_kind_r8
   use kinds,           only: real_kind
   use hybrid_mod,      only: hybrid_t
   use element_mod,     only: element_t
   use derivative_mod,  only: derivative_t
   use time_mod,        only: timelevel_t
   use hybvcoord_mod,   only: hvcoord_t
   use domain_mod,      only: domain1d_t
   use edgetype_mod,    only: EdgeBuffer_t
   use dyn_comp,        only: dyn_import_t
   use spmd_dyn,        only: local_dp_map, block_buf_nrecs, chunk_buf_nrecs
   use spmd_utils,      only: mpicom, iam 
   use bndry_mod,       only: bndry_exchangeV
   use edge_mod_base,   only: edgevpack_nlyr, edgevunpack_nlyr, initEdgeBuffer
   use dimensions_mod,  only: nelem, nelemd
   use edge_mod,        only: edge_g
   use perf_mod,        only: t_startf, t_stopf, t_barrierf

   use ieee_arithmetic

   implicit none

   public :: phys_hyperviscosity_init
   public :: phys_hyperviscosity

   ! type (element_t),  public, dimension(nelem), target :: phys_elem
   type (element_t),  public, pointer :: phys_elem(:) => null()        ! dummy element structure to hold tendencies from physics
   type (domain1d_t), pointer :: dom_mt(:)    => null()

   ! type (EdgeBuffer_t)   :: edge_buffer
   
   contains

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
subroutine phys_hyperviscosity_init()
   use parallel_mod,       only: par
   use dyn_comp,           only: hvcoord, TimeLevel
   use prim_driver_base,   only: prim_init1, prim_init2
   use hybrid_mod,         only: hybrid_create
   use thread_mod,         only: hthreads, omp_get_thread_num
   use element_mod,        only: allocate_element_desc, setup_element_pointers

   integer :: nets, nete   ! start/end of element index
   integer :: ithr
   type (hybrid_t) :: hybrid

   integer, dimension(nelemd) :: gid

   !----------------------------------------------
   ! setup element and TimeLevel structure
   !----------------------------------------------

   ! if (par%dynproc) then

      !!! Took this from prim_init1() in prim_driver_base.F90
      ! allocate( phys_elem(nelemd) )
      ! call setup_element_pointers( phys_elem )
      ! call allocate_element_desc( phys_elem )

      call prim_init1(phys_elem,par,dom_mt,TimeLevel)
      ithr = omp_get_thread_num()
      nets = dom_mt(ithr)%start
      nete = dom_mt(ithr)%end 
      hybrid = hybrid_create(par,ithr,hthreads)
      call prim_init2(phys_elem, hybrid, nets, nete, TimeLevel, hvcoord)

      ! do ie = 1,nelemd
      !    global_id(ie) = phys_elem(ie)%GlobalID
      ! end do
      
      !!! Took this from prim_init1() in prim_driver_base.F90
      ! allocate( phys_elem(nelemd) )
      ! call setup_element_pointers( phys_elem )
      ! call allocate_element_desc( phys_elem )  

      ! do ie = 1,nelemd
      !    phys_elem(ie)%GlobalID = global_id(ie)
      ! end do    

      ! call prim_init2(phys_elem, hybrid, nets, nete, TimeLevel, hvcoord)

   ! end if ! par%dynproc

   !----------------------------------------------
   !----------------------------------------------

end subroutine phys_hyperviscosity_init

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
subroutine phys_hyperviscosity(ptend)
   !----------------------------------------------
   ! Main interface for hyperviscosity calculation on state after being modifed by physics.
   ! This mimics what happens in dynamics hypervis calculation, but is not identical.
   ! Original motivation was to effectively smooth the output tendencies from the CRM.
   ! Author: Walter Hannah (LLNL)
   !----------------------------------------------
   use parallel_mod,       only: par
   use time_mod,           only: tstep
   use hybrid_mod,         only: hybrid_create
   use dyn_comp,           only: hvcoord, TimeLevel
   use thread_mod,         only: hthreads, omp_get_thread_num
   use prim_driver_base,   only: prim_init1, prim_init2
   use derivative_mod_base,only: derivinit
   use dyn_grid,           only: get_gcol_block_d
   use phys_grid,          only: get_ncols_p, get_gcol_all_p,                          &
                                 transpose_block_to_chunk, transpose_chunk_to_block,   &
                                 block_to_chunk_send_pters, block_to_chunk_recv_pters, &
                                 chunk_to_block_send_pters, chunk_to_block_recv_pters
   use dimensions_mod,     only: np, npsq, nelemd, nlev
   use ppgrid,             only: begchunk, endchunk, pcols, pver, pverp
   use dof_mod,            only: PutUniquePoints, UniquePoints
   use physics_types,      only: physics_ptend 
   !!! Interface arguments
   type(physics_ptend), intent(inout) :: ptend(begchunk:endchunk)
   ! type(dyn_import_t),  intent(inout) :: dyn_dum
   !!! Local variables
   ! type (element_t), dimension(nelem), target :: phys_elem
   integer :: nt
   integer :: ithr
   integer :: icol, ncols, lchnk
   integer :: nets, nete            ! start/end of element index
   integer :: ie, ioff, k, i
   integer :: pgcols(pcols), idmb1(1), idmb2(1), idmb3(1)
   type (hybrid_t)            :: hybrid
   type (derivative_t)        :: deriv
   
   integer :: tsize                                                  ! number of variables for transpose buffer
   integer, dimension(pcols,0:pver) :: cpter                         ! offsets into chunk buffer for packing data
   integer, dimension(npsq ,0:pver) :: bpter                         ! offsets into block buffer for packing data
   ! integer, dimension(pcols,pver+1) :: cpter                         ! offsets into chunk buffer for packing data
   ! integer, dimension(npsq ,pver+1) :: bpter                         ! offsets into block buffer for packing data
   real (kind=real_kind), allocatable, dimension(:) :: bbuffer       ! transpose buffer
   real (kind=real_kind), allocatable, dimension(:) :: cbuffer       ! transpose buffer

   real (kind=real_kind), dimension(npsq,pver,nelemd) :: T_tmp       ! temporary array to physics temperature tend
   real (kind=real_kind), dimension(npsq,pver,nelemd) :: q_tmp       ! temporary array to physics sp humidity tend

   !----------------------------------------------
   ! setup element and TimeLevel structure
   !----------------------------------------------

   ithr = omp_get_thread_num()
   nets = dom_mt(ithr)%start
   nete = dom_mt(ithr)%end 

   hybrid = hybrid_create(par,ithr,hthreads)

   ! nt = TimeLevel%n0
   nt = TimeLevel%np1

   !----------------------------------------------
   ! initialize hybrid and derivative structure
   !----------------------------------------------
   hybrid = hybrid_create(par,ithr,hthreads)

   call derivinit(deriv)

   !----------------------------------------------
   ! copy physics columns to element structure
   ! based on cam/src/dynamics/se/dp_coupling.F90
   !----------------------------------------------
   ! if (local_dp_map) then

      ! do lchnk = begchunk,endchunk
      !    ncols = get_ncols_p(lchnk)
      !    call get_gcol_all_p(lchnk,pcols,pgcols)
      !    do icol = 1,ncols
      !       call get_gcol_block_d(pgcols(icol),1,idmb1,idmb2,idmb3)
      !       ie   = idmb3(1)
      !       ioff = idmb2(1)
      !       do k = 1,pver
      !          T_tmp(ioff,k,ie) = ptend(lchnk)%s(icol,k)
      !          q_tmp(ioff,k,ie) = ptend(lchnk)%q(icol,k,1)
      !       end do ! k
      !    end do ! icol
      ! end do ! lchnk

   ! else ! .not. local_dp_map

      tsize = 2
      allocate( bbuffer(tsize*block_buf_nrecs) )
      allocate( cbuffer(tsize*chunk_buf_nrecs) )
      do lchnk = begchunk,endchunk
         ncols = get_ncols_p(lchnk)
         call chunk_to_block_send_pters(lchnk,pcols,pver+1,tsize,cpter)
         do i = 1,ncols
            cbuffer(cpter(i,0):cpter(i,0)+tsize-1) = 0.0_r8
         end do
         do icol = 1,ncols
            do k = 1,pver
               cbuffer (cpter(icol,k)  ) = ptend(lchnk)%s(icol,k)
               cbuffer (cpter(icol,k)+1) = ptend(lchnk)%q(icol,k,1)
            end do
         end do
      end do
      call t_barrierf('sync_chk_to_blk', mpicom)
      call transpose_chunk_to_block(tsize, cbuffer, bbuffer)
      if (par%dynproc) then
         do ie = 1,nelemd
            call chunk_to_block_recv_pters(phys_elem(ie)%GlobalID,npsq,pver+1,tsize,bpter)
            do icol = 1,phys_elem(ie)%idxP%NumUniquePts
               do k = 1,pver
                  T_tmp(icol,k,ie) = bbuffer( bpter(icol,k)   )
                  q_tmp(icol,k,ie) = bbuffer( bpter(icol,k)+1 )
               end do
            end do
         end do
      end if
      deallocate( bbuffer )
      deallocate( cbuffer )

   ! end if ! local_dp_map

   ! if (par%dynproc) then
      do ie = 1,nelemd
         ncols = phys_elem(ie)%idxP%NumUniquePts
         ! ncols = get_ncols_p(lchnk)
         call PutUniquePoints( phys_elem(ie)%idxP, nlev, T_tmp(1:ncols,:,ie), phys_elem(ie)%state%T(:,:,:,nt) )
         call PutUniquePoints( phys_elem(ie)%idxP, nlev, q_tmp(1:ncols,:,ie), phys_elem(ie)%state%q(:,:,:,nt) )
      end do ! ie
   ! end if ! par%dynproc
   
   !----------------------------------------------
   ! diffuse temperature tendency
   !----------------------------------------------

   call phys_hyperviscosity_Tq( phys_elem, hvcoord, hybrid, deriv, nt, nets, nete, dt_in=tstep )

   !----------------------------------------------
   ! copy smoothed tendencies back to physics columns
   !----------------------------------------------
   do ie = 1,nelemd
      call UniquePoints( phys_elem(ie)%idxP, nlev, phys_elem(ie)%state%T(:,:,:,nt), T_tmp(1:ncols,:,ie) )
      call UniquePoints( phys_elem(ie)%idxP, nlev, phys_elem(ie)%state%q(:,:,:,nt), q_tmp(1:ncols,:,ie) )
   end do ! ie
   
   ! if (local_dp_map) then

      ! do lchnk = begchunk,endchunk
      !    ncols = get_ncols_p(lchnk)
      !    call get_gcol_all_p(lchnk,pcols,pgcols)
      !    do icol = 1,ncols
      !       call get_gcol_block_d(pgcols(icol),1,idmb1,idmb2,idmb3)
      !       ie   = idmb3(1)
      !       ioff = idmb2(1)
      !       do k = 1,pver
      !          ptend(lchnk)%s(icol,k)   = T_tmp(ioff,k,ie) 
      !          ptend(lchnk)%q(icol,k,1) = q_tmp(ioff,k,ie)
      !       end do ! k
      !    end do ! icol
      ! end do ! lchnk

   ! else ! .not. local_dp_map

      allocate(bbuffer(tsize*block_buf_nrecs))
      allocate(cbuffer(tsize*chunk_buf_nrecs))
      if (par%dynproc) then
         do ie = 1,nelemd
            call block_to_chunk_send_pters(phys_elem(ie)%GlobalID,npsq,pver+1,tsize,bpter)
            do icol = 1,phys_elem(ie)%idxP%NumUniquePts
               bbuffer(bpter(icol,0):bpter(icol,0)+tsize-1) = 0.0_r8
               do k = 1,pver
                  bbuffer(bpter(icol,k)  ) = T_tmp(icol,k,ie)
                  bbuffer(bpter(icol,k)+1) = q_tmp(icol,k,ie)
               end do ! k
            end do ! icol
         end do ! ie
      else
         bbuffer(:) = 0._r8
      end if
      call t_barrierf ('sync_blk_to_chk', mpicom)
      call transpose_block_to_chunk(tsize, bbuffer, cbuffer)
      do lchnk = begchunk,endchunk
         ncols = get_ncols_p(lchnk)
         call block_to_chunk_recv_pters(lchnk,pcols,pver+1,tsize,cpter)
         do icol = 1,ncols
            do k = 1,pver
               ptend(lchnk)%s (icol,k)   = cbuffer(cpter(icol,k)  )
               ptend(lchnk)%q (icol,k,1) = cbuffer(cpter(icol,k)+1)
            end do ! k
         end do ! icol
      end do ! lchnk
      deallocate( bbuffer )
      deallocate( cbuffer )

   ! end if ! local_dp_map
   !----------------------------------------------
   !----------------------------------------------
end subroutine phys_hyperviscosity

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
subroutine phys_hyperviscosity_Tq( elem, hvcoord, hybrid, deriv, nt, nets, nete, dt_in )
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
   use control_mod,     only: hypervis_order, nu, nu_q, nu_s, nu_p, hypervis_subcycle
   use derivative_mod_base, only: laplace_sphere_wk
   !!! Interface arguments
   type (element_t)     , intent(inout) :: elem(:)
   type (hvcoord_t)     , intent(in   ) :: hvcoord
   type (hybrid_t)      , intent(in   ) :: hybrid
   type (derivative_t)  , intent(in   ) :: deriv
   integer              , intent(in   ) :: nt
   integer              , intent(in   ) :: nets    ! start of element index
   integer              , intent(in   ) :: nete    ! end of element index
   real (kind=real_kind), intent(in   ) :: dt_in
   !!! Local variables
   real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: tensor_T 
   real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: tensor_q 
   real (kind=real_kind), dimension(np,np,nlev)           :: dp
   real (kind=real_kind) :: dt
   integer :: k, i, j, ie, ic, q

   integer                :: factor_subcycle
   real (kind=real_kind)  :: factor_nu 

   !----------------------------------------------
   !----------------------------------------------
   factor_subcycle = 1
   factor_nu       = 1

#if defined( PHYS_HYPERVIS_FACTOR_5X5 )
   factor_subcycle = 5
   factor_nu       = 5.
#endif

   !----------------------------------------------
   ! Exit conditions
   !----------------------------------------------
   ! if ( nu_q == 0 ) return
   ! if ( hypervis_order /= 2 ) return
   ! if ( nu_s==0 .and. nu==0 .and. nu_p==0 ) return

   !----------------------------------------------
   ! hyper viscosity calculation
   !----------------------------------------------
   !!! nu_p = 0 :
   !!!   scale T dissipaton by dp  (conserve IE, dissipate T^2)
   !!! nu_p > 0 :
   !!!   dont scale: T equation IE dissipation matches (to truncation error)
   !!!               IE dissipation from continuity equation (1 deg: to about 0.1 W/m^2)


   dt = dt_in / ( hypervis_subcycle * factor_subcycle )

   do ic = 1, ( hypervis_subcycle * factor_subcycle )

      !------------------------------------
      ! Populate the tensor variables
      !------------------------------------
      do ie = nets, nete
         ! do k = 1, nlev
         !    tensor_T(:,:,k,ie) = elem(ie)%state%T(:,:,k,nt)
         !    tensor_q(:,:,k,ie) = elem(ie)%state%q(:,:,k,nt)
         ! end do ! k
         tensor_T(:,:,:,ie) = elem(ie)%state%T(:,:,:,nt)
         tensor_q(:,:,:,ie) = elem(ie)%state%q(:,:,:,nt)
      end do ! ie

      !------------------------------------
      ! compute laplacian
      !------------------------------------
      call phys_biharmonic_wk_scalar(hybrid, deriv, nets, nete, elem, tensor_T)
      call phys_biharmonic_wk_scalar(hybrid, deriv, nets, nete, elem, tensor_q)

      !------------------------------------
      ! Apply hyperviscosity
      !------------------------------------
      do ie = nets, nete

         !!! Apply hyperviscosity and mass matrix weighting
         do k = 1, nlev
            !!! regular viscosity tensor - this results in too much diffusion and is unstable
            ! tensor_T(:,:,k,ie) = laplace_sphere_wk( elem(ie)%state%T(:,:,k,nt), deriv, elem(ie), var_coef=.false. )
            ! tensor_q(:,:,k,ie) = laplace_sphere_wk( elem(ie)%state%q(:,:,k,nt), deriv, elem(ie), var_coef=.false. )
            do j = 1, np
               do i = 1, np
                  elem(ie)%state%T(i,j,k,nt) = elem(ie)%state%T(i,j,k,nt) * elem(ie)%spheremp(i,j) - dt*nu_s*factor_nu*tensor_T(i,j,k,ie)
                  elem(ie)%state%q(i,j,k,nt) = elem(ie)%state%q(i,j,k,nt) * elem(ie)%spheremp(i,j) - dt*nu_s*factor_nu*tensor_q(i,j,k,ie)
               enddo ! i = 1, np
            enddo ! j = 1, np
         enddo ! k = 1, nlev

         !!! prepare edge buffer
         call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%T(:,:,1:nlev,nt), nlev, nlev*0, nlev*2 )
         call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%q(:,:,1:nlev,nt), nlev, nlev*1, nlev*2 )

      enddo ! ie = nets, nete

      !------------------------------------
      ! update boundaries
      !------------------------------------
      call bndry_exchangeV(hybrid, edge_g )

      !------------------------------------
      ! Unpack and apply inverse weighting
      !------------------------------------
      do ie = nets, nete

         !!! unpack edge buffer
         call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%T(:,:,1:nlev,nt), nlev, nlev*0, nlev*2)
         call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%q(:,:,1:nlev,nt), nlev, nlev*1, nlev*2)
         
         !!! apply inverse mass matrix
         do k = 1, nlev
            elem(ie)%state%T(:,:,k,nt) = elem(ie)%rspheremp(:,:) * elem(ie)%state%T(:,:,k,nt)
            elem(ie)%state%q(:,:,k,nt) = elem(ie)%rspheremp(:,:) * elem(ie)%state%q(:,:,k,nt)
         enddo ! k
         
      enddo ! ie = nets, nete
      !------------------------------------
      !------------------------------------

   enddo ! ic = 1, hypervis_subcycle
   !----------------------------------------------
   !----------------------------------------------
  
end subroutine phys_hyperviscosity_Tq

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
subroutine phys_biharmonic_wk_scalar(hybrid, deriv, nets, nete, elem, qtens)
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
   integer              , intent(in   ) :: nets    ! start of element index
   integer              , intent(in   ) :: nete    ! end of element index
   type (element_t)     , intent(inout) :: elem(:)
   ! type (EdgeBuffer_t)  , intent(inout) :: edge_buffer_in
   real (kind=real_kind), intent(inout),  dimension(np,np,nlev,nets:nete) :: qtens

   !!! Local Variables
   integer :: i,j,k,ie
   real (kind=real_kind), dimension(np,np) :: lap_tmp
   logical var_coef1

   !!! if tensor hyperviscosity with tensor V is used, 
   !!! then biharmonic operator is ( \grad \cdot V \grad ) ( \grad \cdot \grad ) 
   !!! so tensor is only used on second call to phys_laplace_sphere_wk
   var_coef1 = .true.
   if(hypervis_scaling > 0) var_coef1 = .false.

   !----------------------------------------------
   ! Calculate Laplacian and pack edge buffer
   !----------------------------------------------
   do ie = nets,nete
      do k = 1,nlev
         lap_tmp(:,:) = qtens(:,:,k,ie)
         qtens(:,:,k,ie) = laplace_sphere_wk(lap_tmp, deriv, elem(ie), var_coef=var_coef1)
      end do ! k
      !!! prepare edge buffer
      call edgeVpack_nlyr(edge_g, elem(ie)%desc, qtens(:,:,:,ie),nlev, 0, nlev)
   end do ! ie

   !----------------------------------------------
   ! Update boundaries
   !----------------------------------------------
   call bndry_exchangeV(hybrid,edge_g)
   
   !----------------------------------------------
   ! Unpack edge buffer and reapply Laplacian
   !----------------------------------------------
   do ie=nets,nete
      call edgeVunpack_nlyr(edge_g, elem(ie)%desc, qtens(:,:,:,ie), nlev, 0, nlev)
      do k = 1,nlev 
         !!! apply inverse mass matrix
         lap_tmp(:,:) = elem(ie)%rspheremp(:,:)*qtens(:,:,k,ie)
         !!! Apply Laplacian again
         qtens(:,:,k,ie) = laplace_sphere_wk(lap_tmp, deriv, elem(ie), var_coef=.true.)
      end do ! k
   end do ! ie

  !----------------------------------------------
  !----------------------------------------------
end subroutine phys_biharmonic_wk_scalar

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------

#endif /* DIFFUSE_PHYS_TEND */

end module phys_hyperviscosity_mod