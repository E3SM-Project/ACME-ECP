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
subroutine phys_hyperviscosity(ptend,elem)
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
   use constituents,       only: cnst_get_ind
   !!! Interface arguments
   type(physics_ptend), intent(inout) :: ptend(begchunk:endchunk)
   type(element_t),     intent(inout) :: elem(:)
   !!! Local variables
   integer :: nt
   integer :: ithr
   integer :: icol, ncols, lchnk
   integer :: nets, nete            ! start/end of element index
   integer :: ie, ioff, k, i
   integer :: ixcldliq, ixcldice
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

   real (kind=real_kind), dimension(npsq,pver,nelemd) :: T_tend_tmp       ! temporary array to physics temperature tend
   real (kind=real_kind), dimension(npsq,pver,nelemd) :: qv_tend_tmp       ! temporary array to physics sp humidity tend
   real (kind=real_kind), dimension(npsq,pver,nelemd) :: ql_tend_tmp       ! temporary array to physics sp humidity tend
   real (kind=real_kind), dimension(npsq,pver,nelemd) :: qi_tend_tmp       ! temporary array to physics sp humidity tend

   !----------------------------------------------
   !----------------------------------------------
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)
   !----------------------------------------------
   ! setup element and TimeLevel structure
   !----------------------------------------------
   ithr = omp_get_thread_num()
   nets = dom_mt(ithr)%start
   nete = dom_mt(ithr)%end 

   hybrid = hybrid_create(par,ithr,hthreads)

   nt = TimeLevel%n0
   ! nt = TimeLevel%np1

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
      !          T_tend_tmp(ioff,k,ie) = ptend(lchnk)%s(icol,k)
      !          q_tend_tmp(ioff,k,ie) = ptend(lchnk)%q(icol,k,1)
      !       end do ! k
      !    end do ! icol
      ! end do ! lchnk

   ! else ! .not. local_dp_map

      tsize = 4

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
               cbuffer (cpter(icol,k)+2) = ptend(lchnk)%q(icol,k,ixcldliq)
               cbuffer (cpter(icol,k)+3) = ptend(lchnk)%q(icol,k,ixcldice)
            end do
         end do
      end do
      call t_barrierf('sync_chk_to_blk', mpicom)
      call transpose_chunk_to_block(tsize, cbuffer, bbuffer)
      if (par%dynproc) then
         do ie = 1,nelemd
            call chunk_to_block_recv_pters(elem(ie)%GlobalID,npsq,pver+1,tsize,bpter)
            do icol = 1,elem(ie)%idxP%NumUniquePts
            ! call chunk_to_block_recv_pters(phys_elem(ie)%GlobalID,npsq,pver+1,tsize,bpter)
            ! do icol = 1,phys_elem(ie)%idxP%NumUniquePts
               do k = 1,pver
                  T_tend_tmp(icol,k,ie)  = bbuffer( bpter(icol,k)   )
                  qv_tend_tmp(icol,k,ie) = bbuffer( bpter(icol,k)+1 )
                  ql_tend_tmp(icol,k,ie) = bbuffer( bpter(icol,k)+2 )
                  qi_tend_tmp(icol,k,ie) = bbuffer( bpter(icol,k)+3 )
               end do
            end do
         end do
      end if
      deallocate( bbuffer )
      deallocate( cbuffer )

   ! end if ! local_dp_map

   ! if (par%dynproc) then
      do ie = 1,nelemd
         ! ncols = get_ncols_p(lchnk)

         ncols = elem(ie)%idxP%NumUniquePts
         call PutUniquePoints( elem(ie)%idxP, nlev, T_tend_tmp(1:ncols,:,ie),  elem(ie)%T_tend(:,:,:) )
         call PutUniquePoints( elem(ie)%idxP, nlev, qv_tend_tmp(1:ncols,:,ie), elem(ie)%qv_tend(:,:,:) )
         call PutUniquePoints( elem(ie)%idxP, nlev, ql_tend_tmp(1:ncols,:,ie), elem(ie)%ql_tend(:,:,:) )
         call PutUniquePoints( elem(ie)%idxP, nlev, qi_tend_tmp(1:ncols,:,ie), elem(ie)%qi_tend(:,:,:) )

         ! ncols = phys_elem(ie)%idxP%NumUniquePts
         ! call PutUniquePoints( phys_elem(ie)%idxP, nlev, T_tend_tmp(1:ncols,:,ie),  phys_elem(ie)%T_tend(:,:,:) )
         ! call PutUniquePoints( phys_elem(ie)%idxP, nlev, qv_tend_tmp(1:ncols,:,ie), phys_elem(ie)%qv_tend(:,:,:) )
         ! call PutUniquePoints( phys_elem(ie)%idxP, nlev, ql_tend_tmp(1:ncols,:,ie), phys_elem(ie)%ql_tend(:,:,:) )
         ! call PutUniquePoints( phys_elem(ie)%idxP, nlev, qi_tend_tmp(1:ncols,:,ie), phys_elem(ie)%qi_tend(:,:,:) )
         
      end do ! ie
   ! end if ! par%dynproc
   
   !----------------------------------------------
   ! diffuse temperature tendency
   !----------------------------------------------
   call phys_hyperviscosity_Tq( elem, hvcoord, hybrid, deriv, nt, nets, nete, dt_in=tstep )
   ! call phys_hyperviscosity_Tq( phys_elem, hvcoord, hybrid, deriv, nt, nets, nete, dt_in=tstep )

   !----------------------------------------------
   ! copy smoothed tendencies back to physics columns
   !----------------------------------------------
   do ie = 1,nelemd
      call UniquePoints( elem(ie)%idxP, nlev, elem(ie)%T_tend(:,:,:),  T_tend_tmp(1:ncols,:,ie) )
      call UniquePoints( elem(ie)%idxP, nlev, elem(ie)%qv_tend(:,:,:), qv_tend_tmp(1:ncols,:,ie) )
      call UniquePoints( elem(ie)%idxP, nlev, elem(ie)%ql_tend(:,:,:), ql_tend_tmp(1:ncols,:,ie) )
      call UniquePoints( elem(ie)%idxP, nlev, elem(ie)%qi_tend(:,:,:), qi_tend_tmp(1:ncols,:,ie) )

      ! call UniquePoints( phys_elem(ie)%idxP, nlev, phys_elem(ie)%T_tend(:,:,:),  T_tend_tmp(1:ncols,:,ie) )
      ! call UniquePoints( phys_elem(ie)%idxP, nlev, phys_elem(ie)%qv_tend(:,:,:), qv_tend_tmp(1:ncols,:,ie) )
      ! call UniquePoints( phys_elem(ie)%idxP, nlev, phys_elem(ie)%ql_tend(:,:,:), ql_tend_tmp(1:ncols,:,ie) )
      ! call UniquePoints( phys_elem(ie)%idxP, nlev, phys_elem(ie)%qi_tend(:,:,:), qi_tend_tmp(1:ncols,:,ie) )
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
      !          ptend(lchnk)%s(icol,k)   = T_tend_tmp(ioff,k,ie) 
      !          ptend(lchnk)%q(icol,k,1) = q_tend_tmp(ioff,k,ie)
      !       end do ! k
      !    end do ! icol
      ! end do ! lchnk

   ! else ! .not. local_dp_map

      allocate(bbuffer(tsize*block_buf_nrecs))
      allocate(cbuffer(tsize*chunk_buf_nrecs))
      if (par%dynproc) then
         do ie = 1,nelemd
            call block_to_chunk_send_pters(elem(ie)%GlobalID,npsq,pver+1,tsize,bpter)
            do icol = 1,elem(ie)%idxP%NumUniquePts
            ! call block_to_chunk_send_pters(phys_elem(ie)%GlobalID,npsq,pver+1,tsize,bpter)
            ! do icol = 1,phys_elem(ie)%idxP%NumUniquePts
               bbuffer(bpter(icol,0):bpter(icol,0)+tsize-1) = 0.0_r8
               do k = 1,pver
                  bbuffer(bpter(icol,k)  ) = T_tend_tmp(icol,k,ie)
                  bbuffer(bpter(icol,k)+1) = qv_tend_tmp(icol,k,ie)
                  bbuffer(bpter(icol,k)+2) = ql_tend_tmp(icol,k,ie)
                  bbuffer(bpter(icol,k)+3) = qi_tend_tmp(icol,k,ie)
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
               ptend(lchnk)%s (icol,k)          = cbuffer(cpter(icol,k)  )
               ptend(lchnk)%q (icol,k,1)        = cbuffer(cpter(icol,k)+1)
               ptend(lchnk)%q (icol,k,ixcldliq) = cbuffer(cpter(icol,k)+2)
               ptend(lchnk)%q (icol,k,ixcldice) = cbuffer(cpter(icol,k)+3)
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
   real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: tensor_qv 
   real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: tensor_ql 
   real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: tensor_qi 
   ! real (kind=real_kind), dimension(np,np,nlev)           :: dp
   real (kind=real_kind) :: dp_tmp
   real (kind=real_kind) :: dt
   integer :: k, i, j, ie, ic, q
   
   integer                :: factor_subcycle
   real (kind=real_kind)  :: factor_nu    



   !----------------------------------------------
   !----------------------------------------------
   factor_subcycle = 1
   factor_nu       = 1.

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

#ifndef SP_DUMMY_HYPERVIS
      
      !------------------------------------
      ! apply pressure weighting
      !------------------------------------
      ! do ie = nets, nete
      !    do k = 1, nlev
      !       do j = 1, np
      !          do i = 1, np
      !             dp_tmp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
      !                      ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,nt)
      !             elem(ie)%qv_tend(i,j,k) = elem(ie)%qv_tend(i,j,k) * dp_tmp
      !             elem(ie)%ql_tend(i,j,k) = elem(ie)%ql_tend(i,j,k) * dp_tmp
      !             elem(ie)%qi_tend(i,j,k) = elem(ie)%qi_tend(i,j,k) * dp_tmp
      !          enddo ! i = 1, np
      !       enddo ! j = 1, np
      !    enddo ! k = 1, nlev
      ! enddo ! ie = nets, nete

      !------------------------------------
      ! Populate the tensor variables
      !------------------------------------
      do ie = nets, nete
         tensor_T(:,:,:,ie)  = elem(ie)%T_tend(:,:,:)
         tensor_qv(:,:,:,ie) = elem(ie)%qv_tend(:,:,:)
         tensor_ql(:,:,:,ie) = elem(ie)%ql_tend(:,:,:)
         tensor_qi(:,:,:,ie) = elem(ie)%qi_tend(:,:,:)
      end do ! ie

      !------------------------------------
      ! compute laplacian
      !------------------------------------
      call phys_biharmonic_wk_scalar(hybrid, deriv, nets, nete, elem, tensor_T)
      call phys_biharmonic_wk_scalar(hybrid, deriv, nets, nete, elem, tensor_qv)
      call phys_biharmonic_wk_scalar(hybrid, deriv, nets, nete, elem, tensor_ql)
      call phys_biharmonic_wk_scalar(hybrid, deriv, nets, nete, elem, tensor_qi)

      !!! regular viscosity tensor - this results in too much diffusion and is unstable
      ! tensor_T(:,:,k,ie) = laplace_sphere_wk( elem(ie)%T_tend(:,:,k),  deriv, elem(ie), var_coef=.false. )
      ! tensor_q(:,:,k,ie) = laplace_sphere_wk( elem(ie)%qv_tend(:,:,k), deriv, elem(ie), var_coef=.false. )

      !------------------------------------
      ! Apply hyperviscosity
      !------------------------------------
      do ie = nets, nete
         !!! Apply hyperviscosity and mass matrix weighting
         do k = 1, nlev
            do j = 1, np
               do i = 1, np
                  ! dp_tmp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                  !          ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,nt)
                  elem(ie)%T_tend(i,j,k)  = elem(ie)%T_tend(i,j,k)  * elem(ie)%spheremp(i,j) - dt*nu_s*factor_nu*tensor_T(i,j,k,ie)
                  elem(ie)%qv_tend(i,j,k) = elem(ie)%qv_tend(i,j,k) * elem(ie)%spheremp(i,j) - dt*nu_s*factor_nu*tensor_qv(i,j,k,ie)
                  elem(ie)%ql_tend(i,j,k) = elem(ie)%ql_tend(i,j,k) * elem(ie)%spheremp(i,j) - dt*nu_s*factor_nu*tensor_ql(i,j,k,ie)
                  elem(ie)%qi_tend(i,j,k) = elem(ie)%qi_tend(i,j,k) * elem(ie)%spheremp(i,j) - dt*nu_s*factor_nu*tensor_qi(i,j,k,ie)
               enddo ! i = 1, np
            enddo ! j = 1, np
         enddo ! k = 1, nlev
         !!! prepare edge buffer
         call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%T_tend(:,:,1:nlev),  nlev, nlev*0, nlev*4 )
         call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%qv_tend(:,:,1:nlev), nlev, nlev*1, nlev*4 )
         call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%ql_tend(:,:,1:nlev), nlev, nlev*2, nlev*4 )
         call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%qi_tend(:,:,1:nlev), nlev, nlev*3, nlev*4 )
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
         call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%T_tend(:,:,1:nlev),  nlev, nlev*0, nlev*4)
         call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%qv_tend(:,:,1:nlev), nlev, nlev*1, nlev*4)
         call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%ql_tend(:,:,1:nlev), nlev, nlev*2, nlev*4)
         call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%qi_tend(:,:,1:nlev), nlev, nlev*3, nlev*4)
         
         !!! apply inverse mass matrix
         do k = 1, nlev
            elem(ie)%T_tend(:,:,k)  = elem(ie)%rspheremp(:,:) * elem(ie)%T_tend(:,:,k)
            elem(ie)%qv_tend(:,:,k) = elem(ie)%rspheremp(:,:) * elem(ie)%qv_tend(:,:,k)
            elem(ie)%ql_tend(:,:,k) = elem(ie)%rspheremp(:,:) * elem(ie)%ql_tend(:,:,k)
            elem(ie)%qi_tend(:,:,k) = elem(ie)%rspheremp(:,:) * elem(ie)%qi_tend(:,:,k)
         enddo ! k

         ! !!! apply inverse mass matrix and undo pressure weighting
         ! do k = 1, nlev
         !    do j = 1, np
         !       do i = 1, np
         !          dp_tmp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
         !                   ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,nt)
         !          elem(ie)%T_tend(i,j,k)  = elem(ie)%T_tend(i,j,k)  * elem(ie)%rspheremp(i,j)
         !          elem(ie)%qv_tend(i,j,k) = elem(ie)%qv_tend(i,j,k) * elem(ie)%rspheremp(i,j) / dp_tmp
         !          elem(ie)%ql_tend(i,j,k) = elem(ie)%ql_tend(i,j,k) * elem(ie)%rspheremp(i,j) / dp_tmp
         !          elem(ie)%qi_tend(i,j,k) = elem(ie)%qi_tend(i,j,k) * elem(ie)%rspheremp(i,j) / dp_tmp
         !       enddo ! i = 1, np
         !    enddo ! j = 1, np
         ! enddo ! k = 1, nlev
         
      enddo ! ie = nets, nete
      !------------------------------------
      !------------------------------------

#else
      !------------------------------------
      ! DUMMY TESTING VERSION - NO HYPERVIS
      !------------------------------------
      do ie = nets, nete
         !!! apply mass matrix
         do k = 1, nlev
            ! elem(ie)%qv_tend(:,:,k) = elem(ie)%qv_tend(:,:,k) * elem(ie)%spheremp(:,:) 
            
            ! elem(ie)%qv_tend(:,:,k) = elem(ie)%qv_tend(:,:,k) * elem(ie)%spheremp(:,:) * elem(ie)%derived%dpdiss_ave(:,:,k) / hvcoord%dp0(k)

            do j = 1, np
               do i = 1, np
                  dp_tmp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                           ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,nt)
                  elem(ie)%qv_tend(i,j,k) = elem(ie)%qv_tend(i,j,k) * elem(ie)%spheremp(i,j) * dp_tmp
               end do
            end do
            
            ! elem(ie)%qv_tend(:,:,k) = elem(ie)%qv_tend(:,:,k) * elem(ie)%spheremp(:,:) * elem(ie)%state%dp3d(:,:,k,nt)
            ! elem(ie)%qv_tend(:,:,k) = elem(ie)%qv_tend(:,:,k) * elem(ie)%spheremp(:,:) * elem(ie)%derived%dp(:,:,k)

         enddo ! k
         !!! prepare edge buffer
         call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%qv_tend(:,:,1:nlev), nlev, nlev*0, nlev )

      enddo ! ie = nets, nete
      
      !!! boundary exchange
      call bndry_exchangeV(hybrid,edge_g)
      
      do ie = nets, nete
         !!! unpack edge buffer
         call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%qv_tend(:,:,1:nlev), nlev, nlev*0, nlev )
         !!! apply inverse mass matrix
         do k = 1, nlev
            ! elem(ie)%qv_tend(:,:,k) = elem(ie)%qv_tend(:,:,k) * elem(ie)%rspheremp(:,:)

            do j = 1, np
               do i = 1, np
                  dp_tmp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                           ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,nt)
                  elem(ie)%qv_tend(i,j,k) = elem(ie)%qv_tend(i,j,k) * elem(ie)%rspheremp(i,j) / dp_tmp
               end do
            end do

            ! elem(ie)%qv_tend(:,:,k) = elem(ie)%qv_tend(:,:,k) * elem(ie)%rspheremp(:,:) / elem(ie)%state%dp3d(:,:,k,nt)
            ! elem(ie)%qv_tend(:,:,k) = elem(ie)%qv_tend(:,:,k) * elem(ie)%rspheremp(:,:) / elem(ie)%derived%dp(:,:,k)

         enddo ! k
      enddo ! ie = nets, nete
      !------------------------------------
      ! DUMMY TESTING VERSION - NO HYPERVIS
      !------------------------------------
#endif

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