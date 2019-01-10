module gll_grid_mod
!----------------------------------------------------------------------- 
! 
! Purpose: Definition of GLL grid for output of dynamics fields 
!          used when physics is on alternate (non-GLL) grid
!
! Method: Variables are private; interface routines used to extract
!         information for use in user code.
! 
! Entry points:
!      gll_grid_init       initialize chunk'ed data structure
!
! Author: Walter Hannah - adapted from phys_grid.F90
! 
!-----------------------------------------------------------------------
#if ( defined SPMD )
   use spmd_dyn,         only: block_buf_nrecs, chunk_buf_nrecs, local_dp_map
   use mpishorthand
#endif
   use shr_kind_mod,     only: r8 => shr_kind_r8, r4 => shr_kind_r4
   use physconst,        only: pi
   use ppgrid,           only: pcols, pver, begchunk, endchunk
   use spmd_utils,       only: iam, masterproc, npes, proc_smp_map, nsmps
   use m_MergeSorts,     only: IndexSet, IndexSort
   use cam_abortutils,   only: endrun
   use perf_mod
   use cam_logfile,      only: iulog
   use scamMod,          only: single_column, scmlat, scmlon
   use shr_const_mod,    only: SHR_CONST_PI
   use dyn_grid,         only: dyn_grid_find_gcols
   use dycore,           only: dycore_is

   use phys_grid,        only: phys_decomp

   implicit none
   save

#if ( ! defined SPMD )
   integer, private :: block_buf_nrecs
   integer, private :: chunk_buf_nrecs
   logical, private :: local_dp_map=.true. 
#endif

   !!! The identifier for the physics grid
   ! integer, parameter, public :: gll_decomp = phys_decomp
   integer, parameter, public :: gll_decomp = 102

   !!! dynamics field grid information
   integer, private :: hdim1_d, hdim2_d   ! dimensions of rectangular horz grid data structure, If 1D, then hdim2_d == 1

   !!! physics field data structures
   integer, private :: ngcols           ! global column count in gll grid (all)
   integer, private :: ngcols_p         ! global column count in gll grid 
                                       ! (without holes)
   integer, dimension(:), allocatable, private :: dyn_to_latlon_gcol_map   ! map from unsorted (dynamics) to lat/lon sorted grid indices
   ! integer, dimension(:), allocatable, private :: latlon_to_dyn_gcol_map   ! map from lat/lon sorted grid to unsorted (dynamics) indices (not used here)
   ! integer, dimension(:), allocatable, private :: lonlat_to_dyn_gcol_map   ! map from lon/lat sorted grid to unsorted (dynamics) indices (not used here)

   integer, private :: clat_gll_tot ! number of unique latitudes
   integer, private :: clon_gll_tot ! number of unique longitudes

   integer, dimension(:), allocatable, private :: clat_gll_cnt ! number of repeats for each latitude index in latlon ordering 
   integer, dimension(:), allocatable, private :: clat_gll_idx ! for first occurence of lat corresponding to given lat index

   integer, dimension(:), allocatable, private :: clon_gll_cnt ! number of repeats for each longitude index in lonlat ordering
   integer, dimension(:), allocatable, private :: clon_gll_idx ! for first  occurrence of lon corresponding to  given lat index

   real(r8), dimension(:), allocatable :: clat_gll           ! unique latitudes (radians, increasing)
   real(r8), dimension(:), allocatable :: clon_gll           ! unique longitudes (radians, increasing)

   integer, dimension(:), allocatable, private :: lat_gll      ! index into list of unique column latitudes
   integer, dimension(:), allocatable, private :: lon_gll      ! index into list of unique column longitudes

   !!! chunk data structures
   type gll_chunk
     integer  :: ncols                 ! number of vertical columns
     integer  :: gcol(pcols)           ! global physics column indices
     integer  :: lon(pcols)            ! global longitude indices
     integer  :: lat(pcols)            ! global latitude indices
     integer  :: owner                 ! id of process where chunk assigned
     integer  :: lcid                  ! local chunk index
   end type gll_chunk

   integer, private :: nchunks                  ! global chunk count

   type(gll_chunk), dimension(:), allocatable, private :: chunks     ! global computational grid
   integer,         dimension(:), allocatable, private :: npchunks   ! number of chunks assigned to each process

   type gll_lchunk
     integer  :: ncols                 ! number of vertical columns
     integer  :: cid                   ! global chunk index
     integer  :: gcol(pcols)           ! global physics column indices
     real(r8) :: area(pcols)           ! column surface area (from dynamics)
     real(r8) :: wght(pcols)           ! column integration weight (from dynamics)
   end type gll_lchunk

   integer, private :: nlchunks        ! local chunk count
   type (gll_lchunk), dimension(:), allocatable, private :: gll_lchunks  
                                       ! local chunks

   !!! column mapping data structures
   type gll_column_map
     integer  :: chunk                 ! global chunk index
     integer  :: ccol                  ! column ordering in chunk
   end type gll_column_map

   integer, private :: nlcols           ! local column count
   type (gll_column_map), dimension(:), allocatable, private :: pgcols  ! ordered list of columns (for use in gather/scatter)
                                                                        ! NOTE: consistent with local ordering

   !!! column remap data structures
   integer, dimension(:), allocatable, private :: gs_col_num      ! number of columns scattered to each process in field_to_chunk scatter
   integer, dimension(:), allocatable, private :: gs_col_offset   ! offset of columns (-1) in pgcols scattered to each process in field_to_chunk scatter
   integer, dimension(:), allocatable, private :: btofc_blk_num   ! number of grid points scattered to each process in block_to_chunk alltoallv, 
                                                                  ! and gathered from each process in chunk_to_block alltoallv
   integer, dimension(:), allocatable, private :: btofc_chk_num   ! number of grid points gathered from each process in block_to_chunk alltoallv, 
                                                                  ! and scattered to each process in chunk_to_block alltoallv
   type btofc_pters
     integer :: ncols                  ! number of columns in block
     integer :: nlvls                  ! number of levels in columns
     integer, dimension(:,:), pointer :: pter 
   end type btofc_pters
   type (btofc_pters), dimension(:), allocatable, private :: btofc_blk_offset
                                       ! offset in btoc send array (-1) where 
                                       ! (blockid, bcid, k) column should be packed in
                                       ! block_to_chunk alltoallv, AND
                                       ! offset in ctob receive array (-1) from which
                                       ! (blockid, bcid, k) column should be unpacked in
                                       ! chunk_to_block alltoallv

   type (btofc_pters), dimension(:), allocatable, private :: btofc_chk_offset
                                       ! offset in btoc receive array (-1) from which
                                       ! (lcid, i, k) data should be unpacked in
                                       ! block_to_chunk alltoallv, AND
                                       ! offset in ctob send array (-1) where
                                       ! (lcid, i, k) data should be packed in
                                       ! chunk_to_block alltoallv

   !!! miscellaneous phys_grid data
   integer, private :: dp_coup_steps   ! number of swaps in transpose algorithm
   integer, dimension(:), private, allocatable :: dp_coup_proc
                                       ! swap partner in each step of 
                                       !  transpose algorithm
   logical :: gll_grid_set = .false.   ! flag indicates physics grid has been set
   integer, private :: max_nproc_smpx  ! maximum number of processes assigned to a
                                       !  single virtual SMP used to define physics 
                                       !  load balancing
   integer, private :: nproc_busy_d    ! number of processes active during the dynamics
                                       !  (assigned a dynamics block)

   ! Physics grid decomposition options:  
   ! -1: each chunk is a dynamics block
   !  0: chunk definitions and assignments do not require interprocess comm.
   !  1: chunk definitions and assignments do not require internode comm.
   !  2: chunk definitions and assignments may require communication between all processes
   !  3: chunk definitions and assignments only require communication with one other process
   !  4: concatenated blocks, no load balancing, no interprocess communication
   integer, private, parameter :: min_lbal_opt = -1
   integer, private, parameter :: max_lbal_opt = 5
   integer, private, parameter :: def_lbal_opt = 2               ! default 
   ! integer, private :: lbal_opt = def_lbal_opt
   integer, private :: lbal_opt = -1

   ! Physics grid load balancing options:  
   !  0: assign columns to chunks as single columns, wrap mapped across chunks
   !  1: use (day/night; north/south) twin algorithm to determine load-balanced pairs of 
   !       columns and assign columns to chunks in pairs, wrap mapped
   integer, private, parameter :: min_twin_alg = 0
   integer, private, parameter :: max_twin_alg = 1
   integer, private, parameter :: def_twin_alg_lonlat = 1         ! default
   integer, private, parameter :: def_twin_alg_unstructured = 0
   integer, private :: twin_alg = def_twin_alg_lonlat

   !!! target number of chunks per thread
   integer, private, parameter :: min_chunks_per_thread = 1
   integer, private, parameter :: def_chunks_per_thread = &
                                    min_chunks_per_thread         ! default
   integer, private :: chunks_per_thread = def_chunks_per_thread

   ! Dynamics/physics transpose method for nonlocal load-balance:
   ! -1: use "0" if max_nproc_smpx and nproc_busy_d are both > npes/2; otherwise use "1"
   !  0: use mpi_alltoallv
   !  1: use point-to-point MPI-1 two-sided implementation
   !  2: use point-to-point MPI-2 one-sided implementation if supported, 
   !       otherwise use MPI-1 implementation
   !  3: use Co-Array Fortran implementation if supported, 
   !       otherwise use MPI-1 implementation
   !  11-13: use mod_comm, choosing any of several methods internal to mod_comm.
   !      The method within mod_comm (denoted mod_method) has possible values 0,1,2 and
   !      is set according to mod_method = phys_alltoall - modmin_alltoall, where
   !      modmin_alltoall is 11.
   integer, private, parameter :: min_alltoall = -1
   integer, private, parameter :: max_alltoall = 3
# if defined(MODCM_DP_TRANSPOSE)
   integer, private, parameter :: modmin_alltoall = 11
   integer, private, parameter :: modmax_alltoall = 13
# endif
   integer, private, parameter :: def_alltoall = 1                ! default
   integer, private :: phys_alltoall = def_alltoall

contains

!==================================================================================================
!==================================================================================================
integer function get_ncols_gll(lcid)
   !!! Purpose: Return number of columns in chunk given the local chunk id.
   integer, intent(in)  :: lcid      ! local chunk id
   get_ncols_gll = gll_lchunks(lcid)%ncols
   return
end function get_ncols_gll
!==================================================================================================
!==================================================================================================
subroutine get_gcol_all_gll(lcid, latdim, gcols)
   ! Purpose: Return all global column indices for chunk
   integer, intent(in)  :: lcid           ! local chunk id
   integer, intent(in)  :: latdim         ! declared size of output array
   integer, intent(out) :: gcols(pcols)   ! array of global latitude indices
   integer :: i                           ! loop index
   gcols = -1
   do i = 1,gll_lchunks(lcid)%ncols
      gcols(i) = gll_lchunks(lcid)%gcol(i)
   end do
   return
end subroutine get_gcol_all_gll
!==================================================================================================
!==================================================================================================
real(r8) function get_lat_gll(lcid, col)
   !!! Purpose: Return latitude (in radians) for chunk column
   ! use ppgrid
   integer, intent(in)  :: lcid          ! local chunk id
   integer, intent(in)  :: col           ! column index
   integer :: cid                        ! global chunk id
   cid = gll_lchunks(lcid)%cid
   get_lat_gll = chunks(cid)%lat(col)
   return
end function get_lat_gll
!==================================================================================================
!==================================================================================================
real(r8) function get_lon_gll(lcid, col)
   !!! Purpose: Return longitude (in radians) for chunk column
   ! use ppgrid
   integer, intent(in)  :: lcid          ! local chunk id
   integer, intent(in)  :: col           ! column index
   integer :: cid                        ! global chunk id
   integer :: lat                        ! latitude index
   integer :: gcol                       ! global column id in latlon ordering
   cid = gll_lchunks(lcid)%cid
   lat = chunks(cid)%lat(col)
   gcol = dyn_to_latlon_gcol_map(chunks(cid)%gcol(col))
   get_lon_gll = (gcol - clat_gll_idx(lat)) + 1
   return
end function get_lon_gll
!==================================================================================================
!==================================================================================================
real(r8) function get_rlat_gll(lcid, col)
   !!! Purpose: Return latitude (in radians) for chunk column
   ! use ppgrid
   integer, intent(in)  :: lcid          ! local chunk id
   integer, intent(in)  :: col           ! column index
   integer :: cid                        ! global chunk id
   cid = gll_lchunks(lcid)%cid
   get_rlat_gll = clat_gll( chunks(cid)%lat(col) )
   return
end function get_rlat_gll
!==================================================================================================
!==================================================================================================
real(r8) function get_rlon_gll(lcid, col)
   !!! Purpose: Return longitude (in radians) for chunk column
   ! use ppgrid
   integer, intent(in)  :: lcid          ! local chunk id
   integer, intent(in)  :: col           ! column index
   integer :: cid                        ! global chunk id
   cid = gll_lchunks(lcid)%cid
   get_rlon_gll = clon_gll( chunks(cid)%lon(col) )
   return
end function get_rlon_gll
!==================================================================================================
!==================================================================================================
subroutine get_rlat_gll_all(lcid, rlatdim, rlat)
   !!! Purpose: Return all latitudes (in radians) for chunk
   ! use ppgrid
   integer,  intent(in)  :: lcid             ! local chunk id
   integer,  intent(in)  :: rlatdim          ! declared size of output array
   real(r8), intent(out) :: rlat(rlatdim)    ! array of longitudes
   integer :: i                              ! loop index
   integer :: cid                            ! global chunk id
   cid = gll_lchunks(lcid)%cid
   do i = 1,chunks(cid)%ncols
      rlat(i) = clat_gll( chunks(cid)%lat(i) )
   end do
   return
end subroutine get_rlat_gll_all
!==================================================================================================
!==================================================================================================
subroutine get_rlon_gll_all(lcid, rlondim, rlon)
   !!! Purpose: Return all longitudes (in radians) for chunk
   ! use ppgrid
   integer,  intent(in)  :: lcid             ! local chunk id
   integer,  intent(in)  :: rlondim          ! declared size of output array
   real(r8), intent(out) :: rlon(rlondim)    ! array of longitudes
   integer :: i                              ! loop index
   integer :: cid                            ! global chunk id
   cid = gll_lchunks(lcid)%cid
   do i = 1,chunks(cid)%ncols
      rlon(i) = clon_gll( chunks(cid)%lon(i) )
   end do
   return
end subroutine get_rlon_gll_all
!==================================================================================================
!==================================================================================================
!!! This duplicates init_geo_unique() from physics_types.F90 with a minor change
subroutine gll_state_init_geo_unique( ncol, phys_state )
   use physics_types, only: physics_state
   integer,             intent(in)    :: ncol
   type(physics_state), intent(inout) :: phys_state
   logical :: match
   integer :: i, j, ulatcnt, uloncnt

   phys_state%ulat = -999._r8
   phys_state%ulon = -999._r8
   phys_state%latmapback = 0
   phys_state%lonmapback = 0
   ulatcnt = 0
   uloncnt = 0

   do i = 1,ncol
      
      match = .false.
      do j = 1,ulatcnt
         if ( phys_state%lat(i) .eq. phys_state%ulat(j) ) then
         match=.true.
         phys_state%latmapback(i) = j
         end if
      end do
      if (.not.match) then
         ulatcnt=ulatcnt+1
         phys_state%ulat(ulatcnt) = phys_state%lat(i)
         phys_state%latmapback(i) = ulatcnt
      end if

      match = .false.
      do j = 1,uloncnt
         if ( phys_state%lon(i) .eq. phys_state%ulon(j) ) then
            match=.true.
            phys_state%lonmapback(i) = j
         end if
      end do
      if (.not.match) then
         uloncnt=uloncnt+1
         phys_state%ulon(uloncnt) = phys_state%lon(i)
         phys_state%lonmapback(i) = uloncnt
      end if

   end do ! i

   phys_state%uloncnt=uloncnt
   phys_state%ulatcnt=ulatcnt

   call get_gcol_all_gll(phys_state%lchnk,pcols,phys_state%cid)

end subroutine gll_state_init_geo_unique
!==================================================================================================
!==================================================================================================
subroutine gll_grid_init( )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Physics mapping initialization routine:  
    ! 
    ! Method: 
    ! 
    ! Author: John Drake and Patrick Worley
    ! 
    !-----------------------------------------------------------------------
    use pmgrid, only: plev
    use dycore, only: dycore_is
    use dyn_grid, only: get_block_bounds_d, &
         get_block_gcol_d, get_block_gcol_cnt_d, &
         get_block_levels_d, get_block_lvl_cnt_d, &
         get_block_owner_d, &
         get_gcol_block_d, get_gcol_block_cnt_d, &
#if defined( PHYS_GRID_1x1_TEST )
         get_horiz_grid_e, get_horiz_grid_dim_e, &
#endif
         get_horiz_grid_dim_d, get_horiz_grid_d, physgrid_copy_attributes_d
    use spmd_utils, only: pair, ceil2
    use cam_grid_support, only: cam_grid_register, iMap, max_hcoordname_len
    use cam_grid_support, only: horiz_coord_t, horiz_coord_create
    use cam_grid_support, only: cam_grid_attribute_copy
#if defined( PHYS_GRID_1x1_TEST )
    use cam_grid_support, only: cam_grid_attribute_register
    use dimensions_mod,   only: nelem, nelemd, ne
#endif

    !
    !------------------------------Arguments--------------------------------
    !
    !
    !---------------------------Local workspace-----------------------------
    !
    integer :: i, j, jb, k, p             ! loop indices
    integer :: pre_i                      ! earlier index in loop iteration
    integer :: clat_gll_dex, clon_gll_dex     ! indices into unique lat. and lon. arrays
    integer :: maxblksiz                  ! maximum number of columns in a dynamics block
    integer :: beg_dex, end_dex           ! index range
    integer :: cid, lcid                  ! global and local chunk ids
    integer :: max_ncols                  ! upper bound on number of columns in a block
    integer :: ncols                      ! number of columns in current chunk
    integer :: curgcol, curgcol_d         ! current global column index
    integer :: firstblock, lastblock      ! global block indices
    integer :: blksiz                     ! current block size
    integer :: glbcnt, curcnt             ! running grid point counts
    integer :: curp                       ! current process id
    integer :: block_cnt                  ! number of blocks containing data
    ! for a given vertical column
    integer :: numlvl                     ! number of vertical levels in block 
    ! column
    integer :: levels(plev+1)             ! vertical level indices
    integer :: owner_d                    ! process owning given block column
    integer :: owner_p                    ! process owning given chunk column
    integer :: blockids(plev+1)           ! block indices
    integer :: bcids(plev+1)              ! block column indices
    real(r8), parameter :: deg2rad = SHR_CONST_PI/180.0


    ! column surface area (from dynamics)
    real(r8), dimension(:), allocatable :: area_d

    ! column integration weight (from dynamics)
    real(r8), dimension(:), allocatable :: wght_d

    ! chunk global ordering
    integer, dimension(:), allocatable :: pchunkid

    ! permutation array used in physics column sorting;
    ! reused later as work space in (lbal_opt == -1) logic
    integer, dimension(:), allocatable :: cdex

    ! latitudes and longitudes and column area for dynamics columns
    real(r8), dimension(:), allocatable :: clat_d
    real(r8), dimension(:), allocatable :: clon_d
    real(r8), dimension(:), allocatable :: lat_d
    real(r8), dimension(:), allocatable :: lon_d
    real(r8) :: clat_gll_tmp
    real(r8) :: clon_gll_tmp

    ! Maps and values for physics grid
    real(r8),       pointer             :: lonvals(:)
    real(r8),       pointer             :: latvals(:)
    real(r8),               allocatable :: latdeg_p(:)
    real(r8),               allocatable :: londeg_p(:)
    integer(iMap),  pointer             :: grid_map(:,:)
    integer(iMap),  pointer             :: coord_map(:)
    type(horiz_coord_t), pointer        :: lat_coord
    type(horiz_coord_t), pointer        :: lon_coord
    integer                             :: gcols(pcols)
    character(len=max_hcoordname_len), pointer :: copy_attributes(:)
    character(len=max_hcoordname_len)   :: copy_gridname
    logical                             :: unstructured
    real(r8)                            :: lonmin, latmin

    real(r8), dimension(:), allocatable :: local_pe_area


    nullify(lonvals)
    nullify(latvals)
    nullify(grid_map)
    nullify(coord_map)
    nullify(lat_coord)
    nullify(lon_coord)

    ! call t_adj_detailf(-2)
    ! call t_startf("phys_grid_init")

    !-----------------------------------------------------------------------
    !
    ! Initialize physics grid, using dynamics grid
    ! a) column coordinates
    if (single_column .and. dycore_is ('SE')) lbal_opt = -1 !+PAB make this default option for SCM
! #if defined( PHYS_GRID_1x1_TEST )
!     call get_horiz_grid_dim_e(hdim1_d,hdim2_d)
! #else
    call get_horiz_grid_dim_d(hdim1_d,hdim2_d)
! #endif
    ! if (single_column .and. dycore_is('SE')) then
    !   ngcols = 1
    ! else
      ngcols = hdim1_d*hdim2_d
    ! endif
#if defined( PHYS_GRID_1x1_TEST )
    lbal_opt = -1     ! need local mapping of dyn and phys columns (i.e. local_dp_map = .true. )
#endif


    allocate( clat_d(1:ngcols) )
    allocate( clon_d(1:ngcols) )
    allocate( lat_d(1:ngcols) )
    allocate( lon_d(1:ngcols) )
    allocate( cdex(1:ngcols) )

    clat_d = 100000.0_r8
    clon_d = 100000.0_r8
    ! if (single_column .and. dycore_is('SE')) then
    !   lat_d = scmlat
    !   lon_d = scmlon
    !   clat_d = scmlat * deg2rad
    !   clon_d = scmlon * deg2rad
    ! else
! #if defined( PHYS_GRID_1x1_TEST )
!       call get_horiz_grid_e(ngcols, lat_rad_out=clat_d, lon_rad_out=clon_d, lat_deg_out=lat_d, lon_deg_out=lon_d)
! #else
      call get_horiz_grid_d(ngcols, clat_d_out=clat_d, clon_d_out=clon_d, lat_d_out=lat_d, lon_d_out=lon_d)
! #endif /* PHYS_GRID_1x1_TEST */
    ! endif
    latmin = MINVAL(ABS(lat_d))
    lonmin = MINVAL(ABS(lon_d))
!!XXgoldyXX: To do: replace collection above with local physics points

    ! count number of "real" column indices
    ngcols_p = 0
    do i=1,ngcols
       if (clon_d(i) < 100000.0_r8) then
          ngcols_p = ngcols_p + 1
       endif
    enddo

    ! sort over longitude and identify unique longitude coordinates
    call IndexSet(ngcols,cdex)
    call IndexSort(ngcols,cdex,clon_d,descend=.false.)
    clon_gll_tmp = clon_d(cdex(1))
    clon_gll_tot = 1

    do i=2,ngcols_p
       if (clon_d(cdex(i)) > clon_gll_tmp) then
          clon_gll_tot = clon_gll_tot + 1
          clon_gll_tmp = clon_d(cdex(i))
       endif
    enddo

    allocate( clon_gll(1:clon_gll_tot) )
    allocate( clon_gll_cnt(1:clon_gll_tot) )
    allocate( londeg_p(1:clon_gll_tot) )

    pre_i = 1
    clon_gll_tot = 1
    clon_gll(1) = clon_d(cdex(1))
    londeg_p(1) = lon_d(cdex(1))
    do i=2,ngcols_p
       if (clon_d(cdex(i)) > clon_gll(clon_gll_tot)) then
          clon_gll_cnt(clon_gll_tot) = i-pre_i
          pre_i = i
          clon_gll_tot = clon_gll_tot + 1
          clon_gll(clon_gll_tot) = clon_d(cdex(i))
          londeg_p(clon_gll_tot) = lon_d(cdex(i))
       endif
    enddo
    clon_gll_cnt(clon_gll_tot) = (ngcols_p+1)-pre_i

    ! sort over latitude and identify unique latitude coordinates
    call IndexSet(ngcols,cdex)
    call IndexSort(ngcols,cdex,clat_d,descend=.false.)
    clat_gll_tmp = clat_d(cdex(1))
    clat_gll_tot = 1
    do i=2,ngcols_p
       if (clat_d(cdex(i)) > clat_gll_tmp) then
          clat_gll_tot = clat_gll_tot + 1
          clat_gll_tmp = clat_d(cdex(i))
       endif
    enddo

    allocate( clat_gll(1:clat_gll_tot) )
    allocate( clat_gll_cnt(1:clat_gll_tot) )
    allocate( clat_gll_idx(1:clat_gll_tot) )
    allocate( latdeg_p(1:clat_gll_tot) )

    pre_i = 1
    clat_gll_tot = 1
    clat_gll(1) = clat_d(cdex(1))
    latdeg_p(1) = lat_d(cdex(1))
    do i=2,ngcols_p
       if (clat_d(cdex(i)) > clat_gll(clat_gll_tot)) then
          clat_gll_cnt(clat_gll_tot) = i-pre_i
          pre_i = i
          clat_gll_tot = clat_gll_tot + 1
          clat_gll(clat_gll_tot) = clat_d(cdex(i))
          latdeg_p(clat_gll_tot) = lat_d(cdex(i))
       endif
    enddo
    clat_gll_cnt(clat_gll_tot) = (ngcols_p+1)-pre_i

    clat_gll_idx(1) = 1
    do j=2,clat_gll_tot
       clat_gll_idx(j) = clat_gll_idx(j-1) + clat_gll_cnt(j-1)
    enddo

    deallocate(lat_d)
    deallocate(lon_d)

    ! sort by longitude within latitudes
    end_dex = 0
    do j=1,clat_gll_tot
       beg_dex = end_dex + 1
       end_dex = end_dex + clat_gll_cnt(j)
       call IndexSort(cdex(beg_dex:end_dex),clon_d,descend=.false.)
    enddo

    !!! Early clean-up, to minimize memory high water mark (not executing find_partner or find_twin)
    if (((twin_alg .ne. 1) .and. (lbal_opt .ne. 3)) .or. (lbal_opt .eq. -1)) deallocate( clat_gll_cnt)

    !!! save "longitude within latitude" column ordering and determine mapping 
    !!! from unsorted global column index to  unique latitude/longitude indices
    allocate( lat_gll(1:ngcols) )
    allocate( lon_gll(1:ngcols) )
    allocate( dyn_to_latlon_gcol_map(1:ngcols) )
    ! if (lbal_opt .ne. -1) allocate( latlon_to_dyn_gcol_map(1:ngcols_p) )

    clat_gll_dex = 1
    lat_gll = -1
    dyn_to_latlon_gcol_map = -1
    do i=1,ngcols_p
       ! if (lbal_opt .ne. -1) latlon_to_dyn_gcol_map(i) = cdex(i)
       dyn_to_latlon_gcol_map(cdex(i)) = i

       do while ((clat_gll(clat_gll_dex) < clat_d(cdex(i))) .and. &
                 (clat_gll_dex < clat_gll_tot))
          clat_gll_dex = clat_gll_dex + 1
       enddo
       lat_gll(cdex(i)) = clat_gll_dex
    enddo

    !!! sort by latitude within longitudes
    call IndexSet(ngcols,cdex)
    call IndexSort(ngcols,cdex,clon_d,descend=.false.)
    end_dex = 0
    do i=1,clon_gll_tot
       beg_dex = end_dex + 1
       end_dex = end_dex + clon_gll_cnt(i)
       call IndexSort(cdex(beg_dex:end_dex),clat_d,descend=.false.)
    enddo

    !!! Early clean-up, to minimize memory high water mark (not executing find_twin)
    if ((twin_alg .ne. 1) .or. (lbal_opt .eq. -1)) deallocate( clon_gll_cnt )

    !!! save "latitude within longitude" column ordering (only need in find_twin)
    ! if ((twin_alg .eq. 1) .and. (lbal_opt .ne. -1)) allocate( lonlat_to_dyn_gcol_map(1:ngcols_p) )

    clon_gll_dex = 1
    lon_gll = -1
    do i=1,ngcols_p
       ! if ((twin_alg .eq. 1) .and. (lbal_opt .ne. -1)) lonlat_to_dyn_gcol_map(i) = cdex(i)
       do while ((clon_gll(clon_gll_dex) < clon_d(cdex(i))) .and. &
                 (clon_gll_dex < clon_gll_tot))
          clon_gll_dex = clon_gll_dex + 1
       enddo
       lon_gll(cdex(i)) = clon_gll_dex
    enddo

    !!! Clean-up
    deallocate( clat_d )
    deallocate( clon_d )
    deallocate( cdex )

    !!!
    !!! Determine block index bounds
    !!!
    call get_block_bounds_d(firstblock,lastblock)

    !!! Allocate storage to save number of chunks and columns assigned to each
    !!! process during chunk creation and assignment
    allocate( npchunks(0:npes-1) )
    allocate( gs_col_num(0:npes-1) )
    npchunks(:) = 0
    gs_col_num(:) = 0

    !!!
    !!! Option -1: each dynamics block is a single chunk
    !!!          
    if (lbal_opt == -1) then

       !!!
       !!! Check that pcols >= maxblksiz
       !!!
       maxblksiz = 0
       do jb=firstblock,lastblock
          ! if (single_column .and. dycore_is('SE')) then
          !   maxblksiz = 1
          ! else
            maxblksiz = max(maxblksiz,get_block_gcol_cnt_d(jb))
          ! endif
       enddo
       if (pcols < maxblksiz) then
           write(iulog,*) 'pcols = ',pcols, ' maxblksiz=',maxblksiz
          call endrun ('GLL_GRID_INIT error: phys_loadbalance -1 specified but PCOLS < MAXBLKSIZ')
       endif

       !!!
       !!! Determine total number of chunks
       !!!
       ! if (single_column .and. dycore_is('SE')) then
       !   nchunks = 1
       ! else
          nchunks = (lastblock-firstblock+1)
       ! endif

       !
       ! Set max virtual SMP node size
       !
       max_nproc_smpx = 1

       !
       ! Allocate and initialize chunks data structure
       !
       allocate( cdex(1:maxblksiz) )
       allocate( chunks(1:nchunks) )

       do cid=1,nchunks
          ! get number of global column indices in block
          ! if (single_column .and. dycore_is('SE')) then
          !   max_ncols = 1
          ! else
            max_ncols = get_block_gcol_cnt_d(cid+firstblock-1)
          ! endif
! #if defined( PHYS_GRID_1x1_TEST )
!           max_ncols = 1
! #endif
          ! fill cdex array with global indices from current block
          call get_block_gcol_d(cid+firstblock-1,max_ncols,cdex)

          ncols = 0
          do i=1,max_ncols
             ! check whether global index is for a column that dynamics
             ! intends to pass to the physics
             curgcol_d = cdex(i)
             if (dyn_to_latlon_gcol_map(curgcol_d) .ne. -1) then
                ! yes - then save the information
                ncols = ncols + 1
                chunks(cid)%gcol(ncols) = curgcol_d
                chunks(cid)%lat(ncols) = lat_gll(curgcol_d)
                chunks(cid)%lon(ncols) = lon_gll(curgcol_d)      
             endif
          enddo
          chunks(cid)%ncols = ncols
       enddo

       if (masterproc) write(iulog,*) 'gll_grid_init: ncols = ',ncols

       ! Clean-up
       deallocate( cdex )
       deallocate( lat_gll )
       deallocate( lon_gll )

       !
       ! Specify parallel decomposition 
       !
       do cid=1,nchunks
#if (defined SPMD)
          p = get_block_owner_d(cid+firstblock-1)
#else
          p = 0
#endif
          chunks(cid)%owner = p
          npchunks(p)       = npchunks(p) + 1
          gs_col_num(p)     = gs_col_num(p) + chunks(cid)%ncols
       enddo

       !
       ! Set flag indicating columns in physics and dynamics 
       ! decompositions reside on the same processes
       !
       local_dp_map = .true. 
       !
    ! else ! (lbal_opt == -1)
    !    !
    !    ! Option == 0: split local blocks into chunks,
    !    !               while attempting to create load-balanced chunks.
    !    !               Does not work with vertically decomposed blocks.
    !    !               (default)
    !    ! Option == 1: split SMP-local blocks into chunks,
    !    !               while attempting to create load-balanced chunks.
    !    !               Does not work with vertically decomposed blocks.
    !    ! Option == 2: load balance chunks with respect to diurnal and
    !    !               seaonsal cycles and wth respect to latitude, 
    !    !               and assign chunks to processes
    !    !               in a way that attempts to minimize communication costs
    !    ! Option == 3: divide processes into pairs and split 
    !    !               blocks assigned to these pairs into 
    !    !               chunks, attempting to create load-balanced chunks.
    !    !               The process pairs are chosen to maximize load balancing
    !    !               opportunities.
    !    !               Does not work with vertically decomposed blocks.
    !    ! Option == 4: concatenate local blocks, then
    !    !               divide into chunks.
    !    !               Does not work with vertically decomposed blocks.
    !    ! Option == 5: split indiviudal blocks into chunks,
    !    !               assigning columns using block ordering
    !    !
    !    !
    !    ! Allocate and initialize chunks data structure, then
    !    ! assign chunks to processes.
    !    !

    !    if  (twin_alg .eq. 1) then
    !       ! precompute clon_gll_idx: index in lonlat ordering for first 
    !       ! occurrence of longitude corresponding to given latitude index,
    !       ! used in twin option in create_chunks; used in create_chunks
    !       allocate( clon_gll_idx(1:clon_gll_tot) )
    !       clon_gll_idx(1) = 1
    !       do i=2,clon_gll_tot
    !          clon_gll_idx(i) = clon_gll_idx(i-1) + clon_gll_cnt(i-1)
    !       enddo
    !    endif

    !    call create_chunks(lbal_opt, chunks_per_thread)

    !    ! Early clean-up, to minimize memory high water mark
    !    deallocate( lat_gll )
    !    deallocate( lon_gll )
    !    !deallocate( latlon_to_dyn_gcol_map ) !do not deallocate as it is being used in RRTMG radiation.F90
    !    if  (twin_alg .eq. 1) deallocate( lonlat_to_dyn_gcol_map )
    !    if  (twin_alg .eq. 1) deallocate( clon_gll_cnt )
    !    if  (twin_alg .eq. 1) deallocate( clon_gll_idx )
    !    if ((twin_alg .eq. 1) .or. (lbal_opt .eq. 3)) deallocate( clat_gll_cnt )

    !    !
    !    ! Determine whether dynamics and physics decompositions
    !    ! are colocated, not requiring any interprocess communication
    !    ! in the coupling.
    !    local_dp_map = .true.   
    !    do cid=1,nchunks
    !       do i=1,chunks(cid)%ncols
    !          curgcol_d = chunks(cid)%gcol(i)
    !          block_cnt = get_gcol_block_cnt_d(curgcol_d)
    !          call get_gcol_block_d(curgcol_d,block_cnt,blockids,bcids)
    !          do jb=1,block_cnt
    !             owner_d = get_block_owner_d(blockids(jb)) 
    !             if (owner_d .ne. chunks(cid)%owner) then
    !                local_dp_map = .false.   
    !             endif
    !          enddo
    !       enddo
    !    enddo
    end if ! lbal_opt == -1

    !!!
    !!! Allocate and initialize data structures for gather/scatter
    !!!  

    allocate( pgcols(1:ngcols_p) )
    allocate( gs_col_offset(0:npes) )
    allocate( pchunkid(0:npes) )

    ! Initialize pchunkid and gs_col_offset by summing 
    ! number of chunks and columns per process, respectively
    pchunkid(0) = 0
    gs_col_offset(0) = 0
    do p=1,npes-1
       pchunkid(p)      = pchunkid(p-1)      + npchunks(p-1)
       gs_col_offset(p) = gs_col_offset(p-1) + gs_col_num(p-1)
    enddo
    
    ! Determine local ordering via "process id" bin sort
    do cid=1,nchunks
       p = chunks(cid)%owner
       pchunkid(p) = pchunkid(p) + 1

       chunks(cid)%lcid = pchunkid(p) + lastblock

       curgcol = gs_col_offset(p)
       do i=1,chunks(cid)%ncols
          curgcol = curgcol + 1
          pgcols(curgcol)%chunk = cid
          pgcols(curgcol)%ccol = i
       enddo
       gs_col_offset(p) = curgcol
    enddo

    ! Reinitialize pchunkid and gs_col_offset (for real)
    pchunkid(0) = 1
    gs_col_offset(0) = 1
    do p=1,npes-1
       pchunkid(p)      = pchunkid(p-1)      + npchunks(p-1)
       gs_col_offset(p) = gs_col_offset(p-1) + gs_col_num(p-1)
    enddo
    pchunkid(npes)      = pchunkid(npes-1)      + npchunks(npes-1)
    gs_col_offset(npes) = gs_col_offset(npes-1) + gs_col_num(npes-1)

    ! Save local information
    ! (Local chunk index range chosen so that it does not overlap 
    !  {begblock,...,endblock})
    ! 
    nlcols   = gs_col_num(iam)
    nlchunks = npchunks(iam)
    begchunk = pchunkid(iam)   + lastblock
    endchunk = pchunkid(iam+1) + lastblock - 1
    !
    allocate( gll_lchunks(begchunk:endchunk) )
    do cid=1,nchunks
       if (chunks(cid)%owner == iam) then
          lcid = chunks(cid)%lcid
          gll_lchunks(lcid)%ncols = chunks(cid)%ncols
          gll_lchunks(lcid)%cid   = cid
          do i=1,chunks(cid)%ncols
             gll_lchunks(lcid)%gcol(i) = chunks(cid)%gcol(i)
          enddo
       endif
    enddo

    deallocate( pchunkid )
    deallocate( npchunks )

    !!!
    !!!-----------------------------------------------------------------------
    !!!
    !!! Initialize physics grid, using dynamics grid
    !!! b) column area and integration weight

    allocate( area_d(1:ngcols) )
    allocate( wght_d(1:ngcols) )
    area_d = 0.0_r8
    wght_d = 0.0_r8

    ! if (single_column .and. dycore_is('SE')) then
    !   area_d = 4.0_r8*pi
    !   wght_d = 4.0_r8*pi
    ! else
! #if defined( PHYS_GRID_1x1_TEST )
!       call get_horiz_grid_e(ngcols, area_out=area_d, wght_out=wght_d)
! #else
      call get_horiz_grid_d(ngcols, area_d_out=area_d, wght_d_out=wght_d)
! #endif /* PHYS_GRID_1x1_TEST */
    ! endif

    if ( abs(sum(area_d) - 4.0_r8*pi) > 1.e-10_r8 ) then
       write(iulog,*) ' ERROR: sum of areas on globe does not equal 4*pi'
       write(iulog,*) ' sum of areas = ', sum(area_d), sum(area_d)-4.0_r8*pi
       call endrun('phys_grid')
    end if

    if ( abs(sum(wght_d) - 4.0_r8*pi) > 1.e-10_r8 ) then
       write(iulog,*) ' ERROR: sum of integration weights on globe does not equal 4*pi'
       write(iulog,*) ' sum of weights = ', sum(wght_d), sum(wght_d)-4.0_r8*pi
       call endrun('phys_grid')
    end if

    do lcid=begchunk,endchunk
       do i=1,gll_lchunks(lcid)%ncols
          gll_lchunks(lcid)%area(i) = area_d(gll_lchunks(lcid)%gcol(i))
          gll_lchunks(lcid)%wght(i) = wght_d(gll_lchunks(lcid)%gcol(i))
       enddo
    enddo
    deallocate( area_d )
    deallocate( wght_d )

    ! if (.not. local_dp_map) then
    !    !
    !    ! allocate and initialize data structures for transposes
    !    !  
    !    allocate( btofc_blk_num(0:npes-1) )
    !    btofc_blk_num = 0
    !    allocate( btofc_blk_offset(firstblock:lastblock) )
    !    do jb = firstblock,lastblock
    !       nullify( btofc_blk_offset(jb)%pter )
    !    enddo
    !    !
    !    glbcnt = 0
    !    curcnt = 0
    !    curp = 0
    !    do curgcol=1,ngcols_p
    !       cid = pgcols(curgcol)%chunk
    !       i   = pgcols(curgcol)%ccol
    !       owner_p   = chunks(cid)%owner
    !       do while (curp < owner_p)
    !          btofc_blk_num(curp) = curcnt
    !          curcnt = 0
    !          curp = curp + 1
    !       enddo
    !       curgcol_d = chunks(cid)%gcol(i)
    !       block_cnt = get_gcol_block_cnt_d(curgcol_d)
    !       call get_gcol_block_d(curgcol_d,block_cnt,blockids,bcids)
    !       do jb = 1,block_cnt
    !          owner_d = get_block_owner_d(blockids(jb))
    !          if (iam == owner_d) then
    !             if (.not. associated(btofc_blk_offset(blockids(jb))%pter)) then
    !                blksiz = get_block_gcol_cnt_d(blockids(jb))
    !                numlvl = get_block_lvl_cnt_d(blockids(jb),bcids(jb))
    !                btofc_blk_offset(blockids(jb))%ncols = blksiz
    !                btofc_blk_offset(blockids(jb))%nlvls = numlvl
    !                allocate( btofc_blk_offset(blockids(jb))%pter(blksiz,numlvl) )
    !             endif
    !             do k=1,btofc_blk_offset(blockids(jb))%nlvls
    !                btofc_blk_offset(blockids(jb))%pter(bcids(jb),k) = glbcnt
    !                curcnt = curcnt + 1
    !                glbcnt = glbcnt + 1
    !             enddo
    !          endif
    !       enddo
    !    enddo
       
       ! btofc_blk_num(curp) = curcnt
       ! block_buf_nrecs = glbcnt
       ! !  
       ! allocate( btofc_chk_num(0:npes-1) )
       ! btofc_chk_num = 0
       ! allocate( btofc_chk_offset(begchunk:endchunk) )
       ! do lcid=begchunk,endchunk
       !    ncols = lchunks(lcid)%ncols
       !    btofc_chk_offset(lcid)%ncols = ncols
       !    btofc_chk_offset(lcid)%nlvls = pver+1
       !    allocate( btofc_chk_offset(lcid)%pter(ncols,pver+1) )
       ! enddo
       ! !
       ! curcnt = 0
       ! glbcnt = 0
       ! do p=0,npes-1
       !    do curgcol=gs_col_offset(iam),gs_col_offset(iam+1)-1
       !       cid  = pgcols(curgcol)%chunk
       !       owner_p  = chunks(cid)%owner
       !       if (iam == owner_p) then
       !          i    = pgcols(curgcol)%ccol
       !          lcid = chunks(cid)%lcid
       !          curgcol_d = chunks(cid)%gcol(i)
       !          block_cnt = get_gcol_block_cnt_d(curgcol_d)
       !          call get_gcol_block_d(curgcol_d,block_cnt,blockids,bcids)
       !          do jb = 1,block_cnt
       !             owner_d = get_block_owner_d(blockids(jb))
       !             if (p == owner_d) then
       !                numlvl = get_block_lvl_cnt_d(blockids(jb),bcids(jb))
       !                call get_block_levels_d(blockids(jb),bcids(jb),numlvl,levels)
       !                do k=1,numlvl
       !                   btofc_chk_offset(lcid)%pter(i,levels(k)+1) = glbcnt
       !                   curcnt = curcnt + 1
       !                   glbcnt = glbcnt + 1
       !                enddo
       !             endif
       !          enddo
       !       endif
       !    enddo
       !    btofc_chk_num(p) = curcnt
       !    curcnt = 0
       ! enddo
       ! chunk_buf_nrecs = glbcnt

       !!!
       !!! Precompute swap partners and number of steps in point-to-point
       !!! implementations of alltoall algorithm.
       !!! First, determine number of swaps.
       !!!

       ! dp_coup_steps = 0
       ! do i=1,ceil2(npes)-1
       !    p = pair(npes,i,iam)
       !    if (p >= 0) then
       !       if ((btofc_blk_num(p) > 0 .or. btofc_chk_num(p) > 0)) then
       !          dp_coup_steps = dp_coup_steps + 1
       !       end if
       !    end if
       ! end do

       !!!
       !!! Second, determine swap partners.
       !!!

       ! allocate( dp_coup_proc(dp_coup_steps) )
       ! dp_coup_steps = 0
       ! do i=1,ceil2(npes)-1
       !    p = pair(npes,i,iam)
       !    if (p >= 0) then
       !       if ((btofc_blk_num(p) > 0 .or. btofc_chk_num(p) > 0)) then
       !          dp_coup_steps = dp_coup_steps + 1
       !          dp_coup_proc(dp_coup_steps) = p
       !       end if
       !    end if
       ! end do
       ! !
    ! end if ! .not. local_dp_map

    !!! Final clean-up
    deallocate( gs_col_offset )
    !!! (if eliminate get_lon_xxx, can also deallocate clat_gll_idx, and grid_latlon?))

    ! Add physics-package grid to set of CAM grids
    ! physgrid always uses 'lat' and 'lon' as coordinate names; 
    ! If dynamics grid is different, it will use different coordinate names

    ! First, create a map for the physics grid
    ! It's structure will depend on whether or not the physics grid is unstructured
    unstructured = dycore_is('UNSTRUCTURED')
    if (unstructured) then
      allocate(grid_map(3, pcols * (endchunk - begchunk + 1)))
    else
      allocate(grid_map(4, pcols * (endchunk - begchunk + 1)))
    end if
    grid_map = 0
    allocate(latvals(size(grid_map, 2)))
    allocate(lonvals(size(grid_map, 2)))
    p = 0
    do lcid = begchunk, endchunk
      ncols = gll_lchunks(lcid)%ncols
      call get_gcol_all_gll(lcid, pcols, gcols)
      ! collect latvals and lonvals
      cid = gll_lchunks(lcid)%cid
      do i = 1, chunks(cid)%ncols
        latvals(p + i) = latdeg_p(chunks(cid)%lat(i))
        lonvals(p + i) = londeg_p(chunks(cid)%lon(i))
      end do
      if (pcols > ncols) then
        ! Need to set these to detect unused columns
        latvals(p+ncols+1:p+pcols) = 1000.0_r8
        lonvals(p+ncols+1:p+pcols) = 1000.0_r8
      end if

      ! Set grid values for this chunk
      do i = 1, pcols
        p = p + 1
        grid_map(1, p) = i
        grid_map(2, p) = lcid
        if ((i <= ncols) .and. (gcols(i) > 0)) then
          if (unstructured) then
            grid_map(3, p) = gcols(i)
          else
            grid_map(3, p) = get_lon_gll(lcid, i)
            grid_map(4, p) = get_lat_gll(lcid, i)
          end if
        else
          if (i <= ncols) then
            call endrun("gll_grid_init: unmapped column")
          end if
        end if
      end do ! i 
    end do ! lcid

    ! Note that if the dycore is using the same points as the physics grid,
    !      it will have already set up 'lat' and 'lon' axes for the physics grid
    !      However, these will be in the dynamics decomposition

    if (unstructured) then
      coord_map => grid_map(3,:)
      lon_coord => horiz_coord_create('lon', 'ncol', ngcols_p, 'longitude',   &
           'degrees_east', 1, size(lonvals), lonvals, map=coord_map)
      lat_coord => horiz_coord_create('lat', 'ncol', ngcols_p, 'latitude',    &
           'degrees_north', 1, size(latvals), latvals, map=coord_map)
    ! else
    !   ! Create a lon coord map which only writes from one of each unique lon
    !   allocate(coord_map(size(grid_map, 2)))
    !   where(latvals == latmin)
    !     coord_map(:) = grid_map(3, :)
    !   elsewhere
    !     coord_map(:) = 0_iMap
    !   end where
    !   lon_coord => horiz_coord_create('lon', 'lon', hdim1_d, 'longitude',     &
    !        'degrees_east', 1, size(lonvals), lonvals, map=coord_map)
    !   nullify(coord_map)
    !   ! Create a lat coord map which only writes from one of each unique lat
    !   allocate(coord_map(size(grid_map, 2)))
    !   where(lonvals == lonmin)
    !     coord_map(:) = grid_map(4, :)
    !   elsewhere
    !     coord_map(:) = 0_iMap
    !   end where
    !   lat_coord => horiz_coord_create('lat', 'lat', hdim2_d, 'latitude',      &
    !        'degrees_north', 1, size(latvals), latvals, map=coord_map)
    end if
    ! nullify(coord_map)
    ! call cam_grid_register('physgrid', phys_decomp, lat_coord, lon_coord,     &
    !      grid_map, unstruct=unstructured, block_indexed=.true.)
! #if defined( PHYS_GRID_1x1_TEST )
    !----------------------------------------------------------------------------
    ! Normally grid attributes are registered in dyn_grid and copied to physgrid
    ! here, but for the 1x1 grid we need to directly regsiter attributes here.
    !----------------------------------------------------------------------------
    allocate( local_pe_area( pcols*(endchunk-begchunk+1) ) )
    p = 0
    do lcid = begchunk, endchunk
      do i = 1,gll_lchunks(lcid)%ncols
        p = p+1
        local_pe_area(p) = gll_lchunks(lcid)%area(i)
      end do ! i
    end do ! lcid
    ! call cam_grid_attribute_register('physgrid','area','physics grid areas','ncol', local_pe_area, map=coord_map)
    ! call cam_grid_attribute_register('physgrid','ne','',ne)
    ! call cam_grid_attribute_register('physgrid','pg','',1)
    !----------------------------------------------------------------------------
    ! In order to get output on the dynamics grid (GLL), we need to 
    ! register a 3rd grid consisting of the unique GLL nodes.
    !----------------------------------------------------------------------------
    call cam_grid_register('gll_grid', gll_decomp, lat_coord, lon_coord,     &
         grid_map, unstruct=unstructured, block_indexed=.true.)
    
    !!! Copy required attributes from the dynamics array
    ! nullify(copy_attributes)
    ! call physgrid_copy_attributes_d(copy_gridname, copy_attributes)
    ! do i = 1, size(copy_attributes)
    !   call cam_grid_attribute_copy(copy_gridname, 'gll_grid', copy_attributes(i))
    ! end do

    !!! register unique attributes
    call cam_grid_attribute_register('gll_grid','gll_area','physics grid areas','ncol', local_pe_area, map=coord_map)
    call cam_grid_attribute_register('gll_grid','gll_ne','',ne)
    call cam_grid_attribute_register('gll_grid','gll_pg','',1)
    nullify(coord_map)

    !----------------------------------------------------------------------------
    !----------------------------------------------------------------------------
! #else /* PHYS_GRID_1x1_TEST */
    !!! Copy required attributes from the dynamics array
    ! nullify(copy_attributes)
    ! call physgrid_copy_attributes_d(copy_gridname, copy_attributes)
    ! do i = 1, size(copy_attributes)
    !   call cam_grid_attribute_copy(copy_gridname, 'physgrid', copy_attributes(i))
    ! end do
! #endif /* PHYS_GRID_1x1_TEST */
    ! Cleanup pointers (they belong to the grid now)
    nullify(grid_map)
    deallocate(latvals)
    nullify(latvals)
    deallocate(lonvals)
    nullify(lonvals)
    ! Cleanup, we are responsible for copy attributes
    if (associated(copy_attributes)) then
      deallocate(copy_attributes)
      nullify(copy_attributes)
    end if

    gll_grid_set = .true.   ! Set flag indicating physics grid is now set
    
    ! if (masterproc) then
    !    write(iulog,*) 'GLL_GRID_INIT:  Using PCOLS=',pcols,     &
    !         '  phys_loadbalance=',lbal_opt,            &
    !         '  phys_twin_algorithm=',twin_alg,         &
    !         '  phys_alltoall=',phys_alltoall,          &
    !         '  chunks_per_thread=',chunks_per_thread
    ! endif
    

    ! call t_stopf("gll_grid_init")
    ! call t_adj_detailf(+2)

    return
end subroutine gll_grid_init
!==================================================================================================
!==================================================================================================
end module gll_grid_mod
