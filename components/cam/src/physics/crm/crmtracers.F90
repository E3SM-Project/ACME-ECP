module crmtracers

  ! This module serves as a template for adding tracer transport in the model. The tracers can be
  ! chemical tracers, or bin microphysics drop/ice categories, etc.
  ! The number of tracers is set by the parameter ntracers which is set in domain.f90.
  ! Also, the logical flag dotracers should be set to .true. in namelist (default is .false.).
  ! The model will transport the tracers around automatically (advection and SGS diffusion).
  ! The user must supply the initialization in the subroutine tracers_init() in this module.
  ! By default, the surface flux of all tracers is zero. Nonzero values can be set in tracers_flux().
  ! The local sinks/sources of tracers should be supplied in tracers_physics().



  use grid
  use params, only: crm_rknd
  use utils,  only: lenstr
  implicit none

  real(crm_rknd), allocatable :: tracer  (:,:,:,:,:) !REDIM  
  real(crm_rknd), allocatable :: fluxbtr (:,:,:,:)   !REDIM ! surface flux of tracers
  real(crm_rknd), allocatable :: fluxttr (:,:,:,:)   !REDIM ! top boundary flux of tracers
  real(crm_rknd), allocatable :: trwle   (:,:,:)     !REDIM ! resolved vertical flux
  real(crm_rknd), allocatable :: trwsb   (:,:,:)     !REDIM ! SGS vertical flux
  real(crm_rknd), allocatable :: tradv   (:,:,:)     !REDIM ! tendency due to vertical advection
  real(crm_rknd), allocatable :: trdiff  (:,:,:)     !REDIM ! tendency due to vertical diffusion
  real(crm_rknd), allocatable :: trphys  (:,:,:)     !REDIM ! tendency due to physics
  character *4 tracername(0:ntracers)
  character *10 tracerunits(0:ntracers)

CONTAINS

  subroutine allocate_tracers(ncrms)
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) :: zero
    allocate( tracer  (ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 0:ntracers) )
    allocate( fluxbtr (ncrms,nx, ny, 0:ntracers) ) ! surface flux of tracers
    allocate( fluxttr (ncrms,nx, ny, 0:ntracers) ) ! top boundary flux of tracers
    allocate( trwle   (ncrms,nz,0:ntracers) )  ! resolved vertical flux
    allocate( trwsb   (ncrms,nz,0:ntracers) )  ! SGS vertical flux
    allocate( tradv   (ncrms,nz,0:ntracers) )  ! tendency due to vertical advection
    allocate( trdiff  (ncrms,nz,0:ntracers) )  ! tendency due to vertical diffusion
    allocate( trphys  (ncrms,nz,0:ntracers) )  ! tendency due to physics

    zero = 0

    tracer  = zero 
    fluxbtr = zero 
    fluxttr = zero 
    trwle   = zero 
    trwsb   = zero 
    tradv   = zero 
    trdiff  = zero 
    trphys  = zero 
  end subroutine allocate_tracers

  subroutine deallocate_tracers
    implicit none
    deallocate( tracer  )
    deallocate( fluxbtr ) ! surface flux of tracers
    deallocate( fluxttr ) ! top boundary flux of tracers
    deallocate( trwle   ) ! resolved vertical flux
    deallocate( trwsb   ) ! SGS vertical flux
    deallocate( tradv   ) ! tendency due to vertical advection
    deallocate( trdiff  ) ! tendency due to vertical diffusion
    deallocate( trphys  ) ! tendency due to physics
  end subroutine deallocate_tracers

  subroutine tracers_init(ncrms,icrm)
    implicit none
    integer, intent(in) :: ncrms,icrm
    integer k,ntr
    character *2 ntrchar
    ! integer, external :: lenstr

    tracer (icrm,:,:,:,:) = 0.
    fluxbtr(icrm,:,:,:) = 0.
    fluxttr(icrm,:,:,:) = 0.

    ! Add your initialization code here. Default is to set to 0 in setdata.f90.

    if(nrestart.eq.0) then

      ! here ....

    end if

    ! Specify te tracers' default names:

    ! Default names are TRACER01, TRACER02, etc:

    do ntr = 1,ntracers
      write(ntrchar,'(i2)') ntr
      do k=1,3-lenstr(ntrchar)-1
        ntrchar(k:k)='0'
      end do
      tracername(ntr) = 'TR'//ntrchar(1:2)
      tracerunits(ntr) = '[TR]'
    end do

  end subroutine tracers_init



  subroutine tracers_flux()

    ! Set surface and top fluxes of tracers. Default is 0 set in setdata.f90

  end subroutine tracers_flux



  subroutine tracers_physics(ncrms,icrm)
    implicit none
    integer, intent(in) :: ncrms,icrm

    ! add here a call to a subroutine that does something to tracers besides advection and diffusion.
    ! The transport is done automatically.

    trphys(icrm,:,:) = 0. ! Default tendency due to physics. You code should compute this to output statistics.

  end subroutine tracers_physics



  subroutine tracers_hbuf_init(namelist,deflist,unitlist,status,average_type,count,trcount)

    ! Initialize the list of tracers statistics variables written in statistics.f90

    character(*) namelist(*), deflist(*), unitlist(*)
    integer status(*),average_type(*),count,trcount
    integer ntr


    do ntr=1,ntracers

      count = count + 1
      trcount = trcount + 1
      namelist(count) = trim(tracername(ntr))
      deflist(count) = trim(tracername(ntr))
      unitlist(count) = trim(tracerunits(ntr))
      status(count) = 1
      average_type(count) = 0
      count = count + 1
      trcount = trcount + 1
      namelist(count) = trim(tracername(ntr))//'FLX'
      deflist(count) = 'Total flux of '//trim(tracername(ntr))
      unitlist(count) = trim(tracerunits(ntr))//' kg/m2/s'
      status(count) = 1
      average_type(count) = 0
      count = count + 1
      trcount = trcount + 1
      namelist(count) = trim(tracername(ntr))//'FLXS'
      deflist(count) = 'SGS flux of '//trim(tracername(ntr))
      unitlist(count) = trim(tracerunits(ntr))//' kg/m2/s'
      status(count) = 1
      average_type(count) = 0
      count = count + 1
      trcount = trcount + 1
      namelist(count) = trim(tracername(ntr))//'ADV'
      deflist(count) = 'Tendency of '//trim(tracername(ntr)//'due to vertical advection')
      unitlist(count) = trim(tracerunits(ntr))//'/day'
      status(count) = 1
      average_type(count) = 0
      count = count + 1
      trcount = trcount + 1
      namelist(count) = trim(tracername(ntr))//'DIFF'
      deflist(count) = 'Tendency of '//trim(tracername(ntr)//'due to vertical SGS transport')
      unitlist(count) = trim(tracername(ntr))//'/day'
      status(count) = 1
      average_type(count) = 0
      count = count + 1
      trcount = trcount + 1
      namelist(count) = trim(tracername(ntr))//'PHYS'
      deflist(count) = 'Tendency of '//trim(tracername(ntr)//'due to physics')
      unitlist(count) = trim(tracername(ntr))//'/day'
      status(count) = 1
      average_type(count) = 0
    end do

  end subroutine tracers_hbuf_init

end module crmtracers
