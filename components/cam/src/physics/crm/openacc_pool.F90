
module openacc_pool
  use iso_c_binding, only: c_ptr, c_f_pointer, c_loc
  implicit none
  private

  !All size / location variables are integer(8) simply to ensure they never overflow
  !Pool is in integer(8) to keep double precision alignment for vectorization

  !No more than a million allocated variables (seems reasonable to me)
  integer(8), parameter :: stack_tot = 1000000
  !Pool hands out data in 8-byte chunks to maintain double alignment for vectorization
  integer(8), pointer :: pool(:)
  !Size of the pool (specified by user during init)
  integer(8) :: pool_size
  !Current head of the stack-based pool (in terms of index, not byte)
  integer(8) :: pool_head
  !Lengths of the entries in the stack
  integer(8) :: stack_lengths(stack_tot)
  !Current head of the stack entries
  integer(8) :: stack_head
  !Dummy variable for sizeof() to know how many bytes are in one index of the pool
  integer(8) :: dummy

  !Interface to map a pointer to various Fortran types / shapes
  interface pool_push
     module procedure pool_push_log_1d
     module procedure pool_push_log_2d
     module procedure pool_push_log_3d
     module procedure pool_push_log_4d
     module procedure pool_push_log_5d
     module procedure pool_push_log_6d
     module procedure pool_push_log_7d

     module procedure pool_push_int_1d
     module procedure pool_push_int_2d
     module procedure pool_push_int_3d
     module procedure pool_push_int_4d
     module procedure pool_push_int_5d
     module procedure pool_push_int_6d
     module procedure pool_push_int_7d

     module procedure pool_push_int8_1d
     module procedure pool_push_int8_2d
     module procedure pool_push_int8_3d
     module procedure pool_push_int8_4d
     module procedure pool_push_int8_5d
     module procedure pool_push_int8_6d
     module procedure pool_push_int8_7d

     module procedure pool_push_real_1d
     module procedure pool_push_real_2d
     module procedure pool_push_real_3d
     module procedure pool_push_real_4d
     module procedure pool_push_real_5d
     module procedure pool_push_real_6d
     module procedure pool_push_real_7d

     module procedure pool_push_real8_1d
     module procedure pool_push_real8_2d
     module procedure pool_push_real8_3d
     module procedure pool_push_real8_4d
     module procedure pool_push_real8_5d
     module procedure pool_push_real8_6d
     module procedure pool_push_real8_7d

     module procedure pool_push_real8_nonstd_1d
     module procedure pool_push_real8_nonstd_2d
     module procedure pool_push_real8_nonstd_3d
     module procedure pool_push_real8_nonstd_4d
     module procedure pool_push_real8_nonstd_5d
     module procedure pool_push_real8_nonstd_6d
     module procedure pool_push_real8_nonstd_7d

     module procedure pool_push_real_nonstd_1d
     module procedure pool_push_real_nonstd_2d
     module procedure pool_push_real_nonstd_3d
     module procedure pool_push_real_nonstd_4d
     module procedure pool_push_real_nonstd_5d
     module procedure pool_push_real_nonstd_6d
     module procedure pool_push_real_nonstd_7d

     module procedure pool_push_int8_nonstd_1d
     module procedure pool_push_int8_nonstd_2d
     module procedure pool_push_int8_nonstd_3d
     module procedure pool_push_int8_nonstd_4d
     module procedure pool_push_int8_nonstd_5d
     module procedure pool_push_int8_nonstd_6d
     module procedure pool_push_int8_nonstd_7d

     module procedure pool_push_int_nonstd_1d
     module procedure pool_push_int_nonstd_2d
     module procedure pool_push_int_nonstd_3d
     module procedure pool_push_int_nonstd_4d
     module procedure pool_push_int_nonstd_5d
     module procedure pool_push_int_nonstd_6d
     module procedure pool_push_int_nonstd_7d

     module procedure pool_push_log_nonstd_1d
     module procedure pool_push_log_nonstd_2d
     module procedure pool_push_log_nonstd_3d
     module procedure pool_push_log_nonstd_4d
     module procedure pool_push_log_nonstd_5d
     module procedure pool_push_log_nonstd_6d
     module procedure pool_push_log_nonstd_7d
  end interface

  !Specify the size of the pool, allocate it on host and device (using OpenACC)
  !call pool_init(pool_size)
  !integer(8), intent(in) pool_size_in !size of the pool in bytes
  public :: pool_init

  !Push an element of a certain size onto the pool stack, and assign the pointer's location
  !Two ways to call this routine. The first assumes indices begin with 1. The second allows
  !arbitrary beginnings to the incides (such as -1:5, for instance)
  !
  !call pool_push(ptr,dims)
  ![real,real(8),integer,integer(8),logical], pointer, intent(out) :: ptr(:[,:,:,:,:,:,:])
  !integer, intent(in) :: dims(:) !Shape of the array
  !
  !call pool_push(ptr,dims1,dims2)
  ![real,real(8),integer,integer(8),logical], pointer, intent(out) :: ptr(:[,:,:,:,:,:,:])
  !integer, intent(in) :: dims1(:) !starting arry indices for each dimension
  !integer, intent(in) :: dims2(:) !ending array indices for each dimension
  public :: pool_push

  !Pop one element from the stack
  !call pool_pop()
  public :: pool_pop

  !Pop multiple elements from the stack
  !call pool_pop(num)
  !integer, intent(in) :: number of elements to pop from the stack
  public :: pool_pop_multiple

  !Deallocate the pool on host and device
  !No parameters needed
  public :: pool_finalize


contains


  !Specify the size of the pool, allocate it on host and device (using OpenACC)
  !INPUT: pool_size_in (size of the pool in bytes)
  subroutine pool_init(pool_size_in)
    implicit none
    integer(8), intent(in) :: pool_size_in
    !Calculate number of integers taken up by [size] bytes
    pool_size = ceiling( pool_size_in / real( sizeof(dummy) ) )
    allocate(pool(pool_size))
    !write(*,*) 'Initializing pool of size: ',pool_size*sizeof(dummy),' bytes'
    !$acc enter data pcreate(pool)
    pool_head = 1
    stack_head = 1
  end subroutine pool_init


  !Push an element of a certain size onto the pool stack, and assign the pointer's location
  !INPUT : size (size of the element to push to the stack in bytes)
  !OUTPUT: ptr (Assigned pointer to the chunk of size bytes in the stack)
  subroutine pool_push_flat(ptr,size)
    implicit none
    integer(8), pointer, intent(  out) :: ptr(:)
    integer(8)         , intent(in   ) :: size
    integer(8) :: locsize

    !Calculate number of integers taken up by [size] bytes
    locsize = ceiling( size / real( sizeof(dummy) ) )

    !Complain and exit if the pool is too small
    if (pool_head + locsize - 1 > pool_size) then
      write(*,*) 'ERROR: Pool size is too small. Exiting...'
      stop
    endif

    !Complain and exit if the total number of stack entries allowed is too small
    if (stack_head + 1 > stack_tot) then
      write(*,*) 'ERROR: Total stack entries is too small. Exiting...'
      stop
    endif

    ptr => pool( pool_head:(pool_head+locsize-1) )
    pool_head = pool_head + locsize

    stack_lengths(stack_head) = locsize
    stack_head = stack_head + 1

    !write(*,*) 'Pushing size: ', size
  end subroutine pool_push_flat


  !Pop one element from the stack
  subroutine pool_pop()
    implicit none
    if (stack_head-1 < 1) then
      write(*,*) 'ERROR: Popping an element that does not exist.'
      return
    endif
    !write(*,*) 'Popping size: ',stack_lengths(stack_head-1)*sizeof(dummy)
    pool_head = pool_head - stack_lengths(stack_head-1)
    stack_head = stack_head - 1
  end subroutine pool_pop


  !Pop multiple elements from the stack
  !INPUT: num (number of elements to pop)
  subroutine pool_pop_multiple(num)
    implicit none
    integer, intent(in) :: num
    integer :: i
    do i = 1 , num
      call pool_pop()
    enddo
  end subroutine pool_pop_multiple


  !Deallocate the pool on host and device
  subroutine pool_finalize()
    implicit none
    !$acc exit data delete(pool)
    deallocate(pool)
  end subroutine pool_finalize


  subroutine pool_push_log_1d(ptr,dims)
    implicit none
    logical, pointer, intent(  out) :: ptr(:)
    integer         , intent(in   ) :: dims(:)
    logical :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_log_1d
  subroutine pool_push_log_2d(ptr,dims)
    implicit none
    logical, pointer, intent(  out) :: ptr(:,:)
    integer         , intent(in   ) :: dims(:)
    logical :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_log_2d
  subroutine pool_push_log_3d(ptr,dims)
    implicit none
    logical, pointer, intent(  out) :: ptr(:,:,:)
    integer         , intent(in   ) :: dims(:)
    logical :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_log_3d
  subroutine pool_push_log_4d(ptr,dims)
    implicit none
    logical, pointer, intent(  out) :: ptr(:,:,:,:)
    integer         , intent(in   ) :: dims(:)
    logical :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_log_4d
  subroutine pool_push_log_5d(ptr,dims)
    implicit none
    logical, pointer, intent(  out) :: ptr(:,:,:,:,:)
    integer         , intent(in   ) :: dims(:)
    logical :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_log_5d
  subroutine pool_push_log_6d(ptr,dims)
    implicit none
    logical, pointer, intent(  out) :: ptr(:,:,:,:,:,:)
    integer         , intent(in   ) :: dims(:)
    logical :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_log_6d
  subroutine pool_push_log_7d(ptr,dims)
    implicit none
    logical, pointer, intent(  out) :: ptr(:,:,:,:,:,:,:)
    integer         , intent(in   ) :: dims(:)
    logical :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_log_7d


  subroutine pool_push_int8_1d(ptr,dims)
    implicit none
    integer(8), pointer, intent(  out) :: ptr(:)
    integer         , intent(in   ) :: dims(:)
    integer(8) :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_int8_1d
  subroutine pool_push_int8_2d(ptr,dims)
    implicit none
    integer(8), pointer, intent(  out) :: ptr(:,:)
    integer         , intent(in   ) :: dims(:)
    integer(8) :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_int8_2d
  subroutine pool_push_int8_3d(ptr,dims)
    implicit none
    integer(8), pointer, intent(  out) :: ptr(:,:,:)
    integer         , intent(in   ) :: dims(:)
    integer(8) :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_int8_3d
  subroutine pool_push_int8_4d(ptr,dims)
    implicit none
    integer(8), pointer, intent(  out) :: ptr(:,:,:,:)
    integer         , intent(in   ) :: dims(:)
    integer(8) :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_int8_4d
  subroutine pool_push_int8_5d(ptr,dims)
    implicit none
    integer(8), pointer, intent(  out) :: ptr(:,:,:,:,:)
    integer         , intent(in   ) :: dims(:)
    integer(8) :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_int8_5d
  subroutine pool_push_int8_6d(ptr,dims)
    implicit none
    integer(8), pointer, intent(  out) :: ptr(:,:,:,:,:,:)
    integer         , intent(in   ) :: dims(:)
    integer(8) :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_int8_6d
  subroutine pool_push_int8_7d(ptr,dims)
    implicit none
    integer(8), pointer, intent(  out) :: ptr(:,:,:,:,:,:,:)
    integer         , intent(in   ) :: dims(:)
    integer(8) :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_int8_7d


  subroutine pool_push_int_1d(ptr,dims)
    implicit none
    integer, pointer, intent(  out) :: ptr(:)
    integer         , intent(in   ) :: dims(:)
    integer :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_int_1d
  subroutine pool_push_int_2d(ptr,dims)
    implicit none
    integer, pointer, intent(  out) :: ptr(:,:)
    integer         , intent(in   ) :: dims(:)
    integer :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_int_2d
  subroutine pool_push_int_3d(ptr,dims)
    implicit none
    integer, pointer, intent(  out) :: ptr(:,:,:)
    integer         , intent(in   ) :: dims(:)
    integer :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_int_3d
  subroutine pool_push_int_4d(ptr,dims)
    implicit none
    integer, pointer, intent(  out) :: ptr(:,:,:,:)
    integer         , intent(in   ) :: dims(:)
    integer :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_int_4d
  subroutine pool_push_int_5d(ptr,dims)
    implicit none
    integer, pointer, intent(  out) :: ptr(:,:,:,:,:)
    integer         , intent(in   ) :: dims(:)
    integer :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_int_5d
  subroutine pool_push_int_6d(ptr,dims)
    implicit none
    integer, pointer, intent(  out) :: ptr(:,:,:,:,:,:)
    integer         , intent(in   ) :: dims(:)
    integer :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_int_6d
  subroutine pool_push_int_7d(ptr,dims)
    implicit none
    integer, pointer, intent(  out) :: ptr(:,:,:,:,:,:,:)
    integer         , intent(in   ) :: dims(:)
    integer :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_int_7d


  subroutine pool_push_real_1d(ptr,dims)
    implicit none
    real   , pointer, intent(  out) :: ptr(:)
    integer         , intent(in   ) :: dims(:)
    real    :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_real_1d
  subroutine pool_push_real_2d(ptr,dims)
    implicit none
    real   , pointer, intent(  out) :: ptr(:,:)
    integer         , intent(in   ) :: dims(:)
    real    :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_real_2d
  subroutine pool_push_real_3d(ptr,dims)
    implicit none
    real   , pointer, intent(  out) :: ptr(:,:,:)
    integer         , intent(in   ) :: dims(:)
    real    :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_real_3d
  subroutine pool_push_real_4d(ptr,dims)
    implicit none
    real   , pointer, intent(  out) :: ptr(:,:,:,:)
    integer         , intent(in   ) :: dims(:)
    real    :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_real_4d
  subroutine pool_push_real_5d(ptr,dims)
    implicit none
    real   , pointer, intent(  out) :: ptr(:,:,:,:,:)
    integer         , intent(in   ) :: dims(:)
    real    :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_real_5d
  subroutine pool_push_real_6d(ptr,dims)
    implicit none
    real   , pointer, intent(  out) :: ptr(:,:,:,:,:,:)
    integer         , intent(in   ) :: dims(:)
    real    :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_real_6d
  subroutine pool_push_real_7d(ptr,dims)
    implicit none
    real   , pointer, intent(  out) :: ptr(:,:,:,:,:,:,:)
    integer         , intent(in   ) :: dims(:)
    real    :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_real_7d


  subroutine pool_push_real8_1d(ptr,dims)
    implicit none
    real(8), pointer, intent(  out) :: ptr(:)
    integer         , intent(in   ) :: dims(:)
    real(8) :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_real8_1d
  subroutine pool_push_real8_2d(ptr,dims)
    implicit none
    real(8), pointer, intent(  out) :: ptr(:,:)
    integer         , intent(in   ) :: dims(:)
    real(8) :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_real8_2d
  subroutine pool_push_real8_3d(ptr,dims)
    implicit none
    real(8), pointer, intent(  out) :: ptr(:,:,:)
    integer         , intent(in   ) :: dims(:)
    real(8) :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_real8_3d
  subroutine pool_push_real8_4d(ptr,dims)
    implicit none
    real(8), pointer, intent(  out) :: ptr(:,:,:,:)
    integer         , intent(in   ) :: dims(:)
    real(8) :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_real8_4d
  subroutine pool_push_real8_5d(ptr,dims)
    implicit none
    real(8), pointer, intent(  out) :: ptr(:,:,:,:,:)
    integer         , intent(in   ) :: dims(:)
    real(8) :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_real8_5d
  subroutine pool_push_real8_6d(ptr,dims)
    implicit none
    real(8), pointer, intent(  out) :: ptr(:,:,:,:,:,:)
    integer         , intent(in   ) :: dims(:)
    real(8) :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_real8_6d
  subroutine pool_push_real8_7d(ptr,dims)
    implicit none
    real(8), pointer, intent(  out) :: ptr(:,:,:,:,:,:,:)
    integer         , intent(in   ) :: dims(:)
    real(8) :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims)
  end subroutine pool_push_real8_7d







  subroutine pool_push_real8_nonstd_1d(ptr,dims1,dims2)
    implicit none
    real(8), pointer, intent(  out) :: ptr(:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    real(8) :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):dims2(1)) => ptr(:)
  end subroutine pool_push_real8_nonstd_1d
  subroutine pool_push_real8_nonstd_2d(ptr,dims1,dims2)
    implicit none
    real(8), pointer, intent(  out) :: ptr(:,:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    real(8) :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):,dims1(2):) => ptr(:,:)
  end subroutine pool_push_real8_nonstd_2d
  subroutine pool_push_real8_nonstd_3d(ptr,dims1,dims2)
    implicit none
    real(8), pointer, intent(  out) :: ptr(:,:,:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    real(8) :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):,dims1(2):,dims1(3):) => ptr(:,:,:)
  end subroutine pool_push_real8_nonstd_3d
  subroutine pool_push_real8_nonstd_4d(ptr,dims1,dims2)
    implicit none
    real(8), pointer, intent(  out) :: ptr(:,:,:,:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    real(8) :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):,dims1(2):,dims1(3):,dims1(4):) => ptr(:,:,:,:)
  end subroutine pool_push_real8_nonstd_4d
  subroutine pool_push_real8_nonstd_5d(ptr,dims1,dims2)
    implicit none
    real(8), pointer, intent(  out) :: ptr(:,:,:,:,:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    real(8) :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):,dims1(2):,dims1(3):,dims1(4):,dims1(5):) => ptr(:,:,:,:,:)
  end subroutine pool_push_real8_nonstd_5d
  subroutine pool_push_real8_nonstd_6d(ptr,dims1,dims2)
    implicit none
    real(8), pointer, intent(  out) :: ptr(:,:,:,:,:,:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    real(8) :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):,dims1(2):,dims1(3):,dims1(4):,dims1(5):,dims1(6):) => ptr(:,:,:,:,:,:)
  end subroutine pool_push_real8_nonstd_6d
  subroutine pool_push_real8_nonstd_7d(ptr,dims1,dims2)
    implicit none
    real(8), pointer, intent(  out) :: ptr(:,:,:,:,:,:,:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    real(8) :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):,dims1(2):,dims1(3):,dims1(4):,dims1(5):,dims1(6):,dims1(7):) => ptr(:,:,:,:,:,:,:)
  end subroutine pool_push_real8_nonstd_7d


  subroutine pool_push_real_nonstd_1d(ptr,dims1,dims2)
    implicit none
    real   , pointer, intent(  out) :: ptr(:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    real    :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):dims2(1)) => ptr(:)
  end subroutine pool_push_real_nonstd_1d
  subroutine pool_push_real_nonstd_2d(ptr,dims1,dims2)
    implicit none
    real   , pointer, intent(  out) :: ptr(:,:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    real    :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):,dims1(2):) => ptr(:,:)
  end subroutine pool_push_real_nonstd_2d
  subroutine pool_push_real_nonstd_3d(ptr,dims1,dims2)
    implicit none
    real   , pointer, intent(  out) :: ptr(:,:,:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    real    :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):,dims1(2):,dims1(3):) => ptr(:,:,:)
  end subroutine pool_push_real_nonstd_3d
  subroutine pool_push_real_nonstd_4d(ptr,dims1,dims2)
    implicit none
    real   , pointer, intent(  out) :: ptr(:,:,:,:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    real    :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):,dims1(2):,dims1(3):,dims1(4):) => ptr(:,:,:,:)
  end subroutine pool_push_real_nonstd_4d
  subroutine pool_push_real_nonstd_5d(ptr,dims1,dims2)
    implicit none
    real   , pointer, intent(  out) :: ptr(:,:,:,:,:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    real    :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):,dims1(2):,dims1(3):,dims1(4):,dims1(5):) => ptr(:,:,:,:,:)
  end subroutine pool_push_real_nonstd_5d
  subroutine pool_push_real_nonstd_6d(ptr,dims1,dims2)
    implicit none
    real   , pointer, intent(  out) :: ptr(:,:,:,:,:,:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    real    :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):,dims1(2):,dims1(3):,dims1(4):,dims1(5):,dims1(6):) => ptr(:,:,:,:,:,:)
  end subroutine pool_push_real_nonstd_6d
  subroutine pool_push_real_nonstd_7d(ptr,dims1,dims2)
    implicit none
    real   , pointer, intent(  out) :: ptr(:,:,:,:,:,:,:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    real    :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):,dims1(2):,dims1(3):,dims1(4):,dims1(5):,dims1(6):,dims1(7):) => ptr(:,:,:,:,:,:,:)
  end subroutine pool_push_real_nonstd_7d


  subroutine pool_push_int8_nonstd_1d(ptr,dims1,dims2)
    implicit none
    integer(8), pointer, intent(  out) :: ptr(:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    integer(8) :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):dims2(1)) => ptr(:)
  end subroutine pool_push_int8_nonstd_1d
  subroutine pool_push_int8_nonstd_2d(ptr,dims1,dims2)
    implicit none
    integer(8), pointer, intent(  out) :: ptr(:,:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    integer(8) :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):,dims1(2):) => ptr(:,:)
  end subroutine pool_push_int8_nonstd_2d
  subroutine pool_push_int8_nonstd_3d(ptr,dims1,dims2)
    implicit none
    integer(8), pointer, intent(  out) :: ptr(:,:,:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    integer(8) :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):,dims1(2):,dims1(3):) => ptr(:,:,:)
  end subroutine pool_push_int8_nonstd_3d
  subroutine pool_push_int8_nonstd_4d(ptr,dims1,dims2)
    implicit none
    integer(8), pointer, intent(  out) :: ptr(:,:,:,:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    integer(8) :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):,dims1(2):,dims1(3):,dims1(4):) => ptr(:,:,:,:)
  end subroutine pool_push_int8_nonstd_4d
  subroutine pool_push_int8_nonstd_5d(ptr,dims1,dims2)
    implicit none
    integer(8), pointer, intent(  out) :: ptr(:,:,:,:,:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    integer(8) :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):,dims1(2):,dims1(3):,dims1(4):,dims1(5):) => ptr(:,:,:,:,:)
  end subroutine pool_push_int8_nonstd_5d
  subroutine pool_push_int8_nonstd_6d(ptr,dims1,dims2)
    implicit none
    integer(8), pointer, intent(  out) :: ptr(:,:,:,:,:,:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    integer(8) :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):,dims1(2):,dims1(3):,dims1(4):,dims1(5):,dims1(6):) => ptr(:,:,:,:,:,:)
  end subroutine pool_push_int8_nonstd_6d
  subroutine pool_push_int8_nonstd_7d(ptr,dims1,dims2)
    implicit none
    integer(8), pointer, intent(  out) :: ptr(:,:,:,:,:,:,:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    integer(8) :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):,dims1(2):,dims1(3):,dims1(4):,dims1(5):,dims1(6):,dims1(7):) => ptr(:,:,:,:,:,:,:)
  end subroutine pool_push_int8_nonstd_7d


  subroutine pool_push_int_nonstd_1d(ptr,dims1,dims2)
    implicit none
    integer, pointer, intent(  out) :: ptr(:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    integer :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):dims2(1)) => ptr(:)
  end subroutine pool_push_int_nonstd_1d
  subroutine pool_push_int_nonstd_2d(ptr,dims1,dims2)
    implicit none
    integer, pointer, intent(  out) :: ptr(:,:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    integer :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):,dims1(2):) => ptr(:,:)
  end subroutine pool_push_int_nonstd_2d
  subroutine pool_push_int_nonstd_3d(ptr,dims1,dims2)
    implicit none
    integer, pointer, intent(  out) :: ptr(:,:,:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    integer :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):,dims1(2):,dims1(3):) => ptr(:,:,:)
  end subroutine pool_push_int_nonstd_3d
  subroutine pool_push_int_nonstd_4d(ptr,dims1,dims2)
    implicit none
    integer, pointer, intent(  out) :: ptr(:,:,:,:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    integer :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):,dims1(2):,dims1(3):,dims1(4):) => ptr(:,:,:,:)
  end subroutine pool_push_int_nonstd_4d
  subroutine pool_push_int_nonstd_5d(ptr,dims1,dims2)
    implicit none
    integer, pointer, intent(  out) :: ptr(:,:,:,:,:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    integer :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):,dims1(2):,dims1(3):,dims1(4):,dims1(5):) => ptr(:,:,:,:,:)
  end subroutine pool_push_int_nonstd_5d
  subroutine pool_push_int_nonstd_6d(ptr,dims1,dims2)
    implicit none
    integer, pointer, intent(  out) :: ptr(:,:,:,:,:,:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    integer :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):,dims1(2):,dims1(3):,dims1(4):,dims1(5):,dims1(6):) => ptr(:,:,:,:,:,:)
  end subroutine pool_push_int_nonstd_6d
  subroutine pool_push_int_nonstd_7d(ptr,dims1,dims2)
    implicit none
    integer, pointer, intent(  out) :: ptr(:,:,:,:,:,:,:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    integer :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):,dims1(2):,dims1(3):,dims1(4):,dims1(5):,dims1(6):,dims1(7):) => ptr(:,:,:,:,:,:,:)
  end subroutine pool_push_int_nonstd_7d


  subroutine pool_push_log_nonstd_1d(ptr,dims1,dims2)
    implicit none
    logical, pointer, intent(  out) :: ptr(:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    logical :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):dims2(1)) => ptr(:)
  end subroutine pool_push_log_nonstd_1d
  subroutine pool_push_log_nonstd_2d(ptr,dims1,dims2)
    implicit none
    logical, pointer, intent(  out) :: ptr(:,:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    logical :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):,dims1(2):) => ptr(:,:)
  end subroutine pool_push_log_nonstd_2d
  subroutine pool_push_log_nonstd_3d(ptr,dims1,dims2)
    implicit none
    logical, pointer, intent(  out) :: ptr(:,:,:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    logical :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):,dims1(2):,dims1(3):) => ptr(:,:,:)
  end subroutine pool_push_log_nonstd_3d
  subroutine pool_push_log_nonstd_4d(ptr,dims1,dims2)
    implicit none
    logical, pointer, intent(  out) :: ptr(:,:,:,:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    logical :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):,dims1(2):,dims1(3):,dims1(4):) => ptr(:,:,:,:)
  end subroutine pool_push_log_nonstd_4d
  subroutine pool_push_log_nonstd_5d(ptr,dims1,dims2)
    implicit none
    logical, pointer, intent(  out) :: ptr(:,:,:,:,:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    logical :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):,dims1(2):,dims1(3):,dims1(4):,dims1(5):) => ptr(:,:,:,:,:)
  end subroutine pool_push_log_nonstd_5d
  subroutine pool_push_log_nonstd_6d(ptr,dims1,dims2)
    implicit none
    logical, pointer, intent(  out) :: ptr(:,:,:,:,:,:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    logical :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):,dims1(2):,dims1(3):,dims1(4):,dims1(5):,dims1(6):) => ptr(:,:,:,:,:,:)
  end subroutine pool_push_log_nonstd_6d
  subroutine pool_push_log_nonstd_7d(ptr,dims1,dims2)
    implicit none
    logical, pointer, intent(  out) :: ptr(:,:,:,:,:,:,:)
    integer         , intent(in   ) :: dims1(:)
    integer         , intent(in   ) :: dims2(:)
    logical :: dummy
    integer(8), pointer :: iptr(:)
    type(c_ptr) :: cptr
    integer(8) :: nbytes
    nbytes = product(dims2-dims1+1)*sizeof(dummy)
    call pool_push_flat(iptr,nbytes)
    cptr = c_loc(iptr)
    call c_f_pointer(cptr,ptr,shape=dims2-dims1+1)
    ptr(dims1(1):,dims1(2):,dims1(3):,dims1(4):,dims1(5):,dims1(6):,dims1(7):) => ptr(:,:,:,:,:,:,:)
  end subroutine pool_push_log_nonstd_7d


end module openacc_pool
