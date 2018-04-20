module cam_fluxes
   use shr_kind_mod, only: r8=>shr_kind_r8
   use assertions, only: assert, assert_valid, assert_range
   use cam_abortutils, only: endrun
   implicit none

   ! Class to store fluxes on the CAM vertical grid and handle functions specific
   ! to CAM fluxes, like sending values to the physics buffer, exporting surface
   ! fluxes to the surface exchange fields, and sending outputs to the history
   ! buffer.
   type cam_flux_type

      ! RRTMGP may operate on a subset of vertical levels in the case that the CAM
      ! vertical grid contains levels with very small or very large values of
      ! pressure that are outside of those assumed reasonable by RRTMGP (this might
      ! happen in the case of a model that extends into the stratosphere maybe?).
      ! To handle these cases, we allow RRTMGP to operate on its own vertical grid.
      ! The variables ktop and kbot are the indices to the top-most and bottom-most
      ! CAM level used in the RRTMGP grid. That is, the RRTMGP grid will consist of
      ! all CAM levels between ktop and kbot, and we will have a mapping between
      ! the CAM and RRTMGP grid where cam_level(ktop) = rrtmgp_level(1),
      ! cam_level(kbot) = rrtmgp_level(nlevels_rrtmgp).
      integer :: ktop  ! Index of top-most level fluxes were calculated at
      integer :: kbot  ! Index of bottom-most level fluxes were calculated at

      ! Fluxes mapped back to CAM vertical grid. Note that those with a level
      ! dimension are defined at the layer *interfaces*, not at the model level
      ! midpoints.
      real(r8), allocatable :: flux_dn(:,:)
      real(r8), allocatable :: flux_up(:,:)
      real(r8), allocatable :: flux_net(:,:)
      real(r8), allocatable :: bnd_flux_dn(:,:,:)
      real(r8), allocatable :: bnd_flux_up(:,:,:)
      real(r8), allocatable :: bnd_flux_net(:,:,:)
      real(r8), allocatable :: flux_net_bot(:)
      real(r8), allocatable :: flux_net_top(:)

      ! Heating rate, to be calculated from the fluxes at model level interfaces.
      real(r8), allocatable :: heating_rate(:,:)

      ! Variable to hold the "type" of flux a particlar instance of this object
      ! contains, being either "shortwave" or "longwave". This will be used to
      ! decide how to calculate net fluxes (because net is defined as down minus up
      ! for shortwave but up minus down for longwave) and how to map fluxes in this
      ! object to outputs.
      character(len=32) :: flux_type

   contains

      ! Type-bound procedures for this class.
      procedure :: initialize=>cam_fluxes_initialize
      procedure :: finalize=>cam_fluxes_finalize
      procedure :: copy_from_rrtmgp=>cam_fluxes_copy_from_rrtmgp
      procedure :: calculate_heating_rate=>cam_fluxes_calculate_heating_rate
      procedure :: to_pbuf=>cam_fluxes_to_pbuf
      procedure :: export_surface_fluxes=>cam_fluxes_export_surface_fluxes

   end type cam_flux_type

contains
   !-------------------------------------------------------------------------------
   ! Type-bound procedures for cam_fluxes_type class
   !-------------------------------------------------------------------------------
   subroutine cam_fluxes_export_surface_fluxes(this, cam_out)
      use camsrfexch, only: cam_out_t

      class(cam_flux_type), intent(in) :: this
      type(cam_out_t), intent(inout) :: cam_out
      integer :: i

      ! Export surface fluxes
      ! The RRTMGP output isn't broken down into direct and diffuse.  For first cut
      ! put entire flux into direct and leave the diffuse set to zero.  To break the
      ! fluxes into the UV/vis and near-IR bands use the same scheme as for the albedos
      ! which is hardcoded for 14 spectral bands.
      !
      ! sols(pcols)      Direct solar rad on surface (< 0.7)
      ! soll(pcols)      Direct solar rad on surface (>= 0.7)
      !
      ! Near-IR bands (1-9 and 14), 820-16000 cm-1, 0.625-12.195 microns
      !
      ! Put half of band 9 in each of the UV/visible and near-IR values,
      ! since this band straddles 0.7 microns:
      !
      ! UV/visible bands 10-13, 16000-50000 cm-1, 0.200-0.625 micron
      if (trim(this%flux_type) == 'shortwave') then
         cam_out%soll = 0
         cam_out%sols = 0
         cam_out%solld = 0
         cam_out%solsd = 0
         do i = 1,size(this%bnd_flux_dn, 1)
            cam_out%soll(i) &
               = sum(this%bnd_flux_dn(i,this%kbot,1:8)) &
               + 0.5_r8 * this%bnd_flux_dn(i,this%kbot,9) &
               + this%bnd_flux_dn(i,this%kbot,14)

            cam_out%sols(i) &
               = 0.5_r8 * this%bnd_flux_dn(i,this%kbot,9) &
               + sum(this%bnd_flux_dn(i,this%kbot,10:13))

            cam_out%netsw(i) = this%flux_net_bot(i)
         end do
      else if (trim(this%flux_type) == 'longwave') then
         do i = 1,size(this%flux_dn, 1)
            cam_out%flwds(i) = this%flux_dn(i,this%kbot)
         end do
      else
         call endrun('flux_type ' // this%flux_type // ' not known.')
      end if

   end subroutine cam_fluxes_export_surface_fluxes


   subroutine cam_fluxes_copy_from_rrtmgp(this, flux_input, cam_indices)
      use mo_fluxes_byband, only:  ty_fluxes_byband
      class(cam_flux_type), intent(inout) :: this
      type(ty_fluxes_byband), intent(in) :: flux_input
      integer, intent(in), optional :: cam_indices(:)

      integer :: k_cam, k_rad, i
      character(len=32) :: sub_name = 'cam_fluxes_copy_from_rrtmgp'

      ! DEBUG checks on inputs 
      call assert_valid(flux_input%flux_dn, sub_name // ': flux_input%flux_dn')
      call assert_valid(flux_input%flux_up, sub_name // ': flux_input%flux_up')
      call assert_valid(flux_input%flux_net, sub_name // ': flux_input%flux_net')
      call assert_valid(flux_input%bnd_flux_dn, sub_name // ': flux_input%bnd_flux_dn')
      call assert_valid(flux_input%bnd_flux_up, sub_name // ': flux_input%bnd_flux_up')
      call assert_valid(flux_input%bnd_flux_net, sub_name // ': flux_input%bnd_flux_net')

      ! Map fluxes on radiation grid to cam grid
      if (present(cam_indices)) then

         ! Sanity check on cam_indices
         if (size(cam_indices) /= size(flux_input%flux_dn, 1)) then
            call endrun(sub_name // 'indices do not conform.')
         end if

         do i = 1,size(cam_indices)
            k_rad = 1
            do k_cam = this%ktop,this%kbot
               ! Broadband fluxes
               this%flux_dn(cam_indices(i),k_cam) = flux_input%flux_dn(i,k_rad)
               this%flux_up(cam_indices(i),k_cam) = flux_input%flux_up(i,k_rad)
               this%flux_net(cam_indices(i),k_cam) = flux_input%flux_net(i,k_rad)

               ! Band-by-band fluxes
               this%bnd_flux_dn(cam_indices(i),k_cam,:) = flux_input%bnd_flux_dn(i,k_rad,:)
               this%bnd_flux_up(cam_indices(i),k_cam,:) = flux_input%bnd_flux_up(i,k_rad,:)
               this%bnd_flux_net(cam_indices(i),k_cam,:) = flux_input%bnd_flux_net(i,k_rad,:)

               ! Increment rad level index
               k_rad = k_rad + 1
            end do
         end do
      else
         k_rad = 1
         do k_cam = this%ktop,this%kbot
            ! Broadband fluxes
            this%flux_dn(:,k_cam) = flux_input%flux_dn(:,k_rad)
            this%flux_up(:,k_cam) = flux_input%flux_up(:,k_rad)
            this%flux_net(:,k_cam) = flux_input%flux_net(:,k_rad)

            ! Band-by-band fluxes
            this%bnd_flux_dn(:,k_cam,:) = flux_input%bnd_flux_dn(:,k_rad,:)
            this%bnd_flux_up(:,k_cam,:) = flux_input%bnd_flux_up(:,k_rad,:)
            this%bnd_flux_net(:,k_cam,:) = flux_input%bnd_flux_net(:,k_rad,:)

            ! Increment rad level index
            k_rad = k_rad + 1
         end do
      end if

      ! Compute net fluxes at surface and top of model
      if (trim(this%flux_type) == 'longwave') then
         ! Net fluxes are always down minus up in RRTMGP, but CAM expects upward to
         ! be positive for longwave, so we invert the net flux here
         this%flux_net(:,:) = - this%flux_net(:,:)
         this%flux_net_bot(:) = this%flux_up(:,this%kbot) - this%flux_dn(:,this%kbot)
         this%flux_net_top(:) = this%flux_up(:,this%ktop) - this%flux_dn(:,this%ktop)
      else if (trim(this%flux_type) == 'shortwave') then
         this%flux_net_bot(:) = this%flux_dn(:,this%kbot) - this%flux_up(:,this%kbot)
         this%flux_net_top(:) = this%flux_dn(:,this%ktop) - this%flux_up(:,this%ktop)
      else
         call endrun('flux_type ' // this%flux_type // ' not known')
      end if

      ! Check outputs to make sure we set things right
      call assert_valid(this%flux_up, trim(sub_name) // ': flux_up')
      call assert_valid(this%flux_dn, trim(sub_name) // ': flux_dn')
      call assert_valid(this%flux_net, trim(sub_name) // ': flux_net')
      call assert_valid(this%bnd_flux_up, trim(sub_name) // ': bnd_flux_up')
      call assert_valid(this%bnd_flux_dn, trim(sub_name) // ': bnd_flux_dn')
      call assert_valid(this%bnd_flux_net, trim(sub_name) // ': bnd_flux_net')
      call assert_valid(this%flux_net_bot, trim(sub_name) // ': flux_net_bot')
      call assert_valid(this%flux_net_top, trim(sub_name) // ': flux_net_top')

   end subroutine cam_fluxes_copy_from_rrtmgp


   subroutine cam_fluxes_calculate_heating_rate(this, pint)

      use physconst, only: gravit

      ! Inputs
      class(cam_flux_type), intent(inout) :: this 
      real(r8), intent(in) :: pint(:,:)

      ! Loop indices
      integer :: i, k

      ! Everyone needs a name
      character(len=32) :: sub_name = 'calculate_heating_rate'

      ! Check inputs
      call assert_valid(this%flux_net, 'calculate_heating_rate: flux_net')

      ! Loop over levels and calculate heating rates; note that the fluxes *should*
      ! be defined at interfaces, so the loop ktop,kbot and grabbing the current
      ! and next value of k should be safe. ktop should be the top interface, and
      ! kbot + 1 should be the bottom interface.
      !
      ! NOTE: to get heating rate in K/day, normally we would use:
      !
      !     H = dF / dp * g * (sec/day) * (1e-5) / (cpair)
      !
      ! Here we just use
      !
      !     H = dF / dp * g
      !
      ! Why? Something to do with convenience with applying the fluxes to the
      ! heating tendency?
      if (trim(this%flux_type) == 'longwave') then
         do i = 1,size(this%flux_net,1)
            do k = 1,size(this%flux_net,2)-1
               this%heating_rate(i,k) = ( &
                  this%flux_net(i,k+1) - this%flux_net(i,k) &
               ) * gravit / (pint(i,k+1) - pint(i,k)) !pdel(i,k)
            end do
         end do
      else if (trim(this%flux_type) == 'shortwave') then
         do i = 1,size(this%flux_net,1)
            do k = 1,size(this%flux_net,2)-1
               this%heating_rate(i,k) = ( &
                  this%flux_net(i,k) - this%flux_net(i,k+1) &
               ) * gravit / (pint(i,k+1) - pint(i,k)) !pdel(i,k)
            end do
         end do
      else
         call endrun(sub_name // ': flux_type ' // trim(this%flux_type) // ' not defined.')
      end if

   end subroutine cam_fluxes_calculate_heating_rate
         

   subroutine cam_fluxes_to_pbuf(this, pbuf)
      use physics_buffer, only: physics_buffer_desc, pbuf_get_field, &
                                pbuf_get_index

      class(cam_flux_type), intent(in) :: this
      type(physics_buffer_desc), pointer :: pbuf(:)

      ! Pointers to pbuf fields
      real(r8), pointer :: fsds(:)
      real(r8), pointer :: fsns(:)
      real(r8), pointer :: fsnt(:)
      real(r8), pointer :: flns(:)
      real(r8), pointer :: flnt(:)

      integer :: ncol, nlevels
      character(len=32) :: sub_name = 'cam_fluxes_to_pbuf'

      ! Check inputs
      call assert_valid(this%flux_dn, trim(sub_name) // ': this%flux_dn')
      call assert_valid(this%flux_net_bot, trim(sub_name) // ': this%flux_net_bot')
      call assert_valid(this%flux_net_top, trim(sub_name) // ': this%flux_net_top')

      call assert_valid(this%flux_up, trim(sub_name) // ': flux_up')
      call assert_valid(this%flux_dn, trim(sub_name) // ': flux_dn')
      call assert_valid(this%flux_net, trim(sub_name) // ': flux_net')
      call assert_valid(this%bnd_flux_up, trim(sub_name) // ': bnd_flux_up')
      call assert_valid(this%bnd_flux_dn, trim(sub_name) // ': bnd_flux_dn')
      call assert_valid(this%bnd_flux_net, trim(sub_name) // ': bnd_flux_net')
      call assert_valid(this%flux_net_bot, trim(sub_name) // ': flux_net_bot')
      call assert_valid(this%flux_net_top, trim(sub_name) // ': flux_net_top')

      ncol = size(this%flux_up, 1)
      nlevels = size(this%flux_up, 2)
      if (trim(this%flux_type) == 'shortwave') then
         ! Associate pointers
         call pbuf_get_field(pbuf, pbuf_get_index('FSDS'), fsds)
         call pbuf_get_field(pbuf, pbuf_get_index('FSNS'), fsns)
         call pbuf_get_field(pbuf, pbuf_get_index('FSNT'), fsnt)

         ! Copy data
         if (size(this%flux_dn, 2) >= this%kbot) then
            fsds(:ncol) = this%flux_dn(:ncol,this%kbot)
         else
            call endrun(trim(sub_name) // ': number of levels < kbot')
         end if
         fsns(:ncol) = this%flux_net_bot(:ncol)
         fsnt(:ncol) = this%flux_net_top(:ncol)

         ! Check that we set them right
         call assert_valid(fsds(:ncol), trim(sub_name) // ': fsds')
         call assert_valid(fsns(:ncol), trim(sub_name) // ': fsns')
         call assert_valid(fsnt(:ncol), trim(sub_name) // ': fsnt')
      else if (trim(this%flux_type) == 'longwave') then
         ! Associate pointers
         call pbuf_get_field(pbuf, pbuf_get_index('FLNS'), flns)
         call pbuf_get_field(pbuf, pbuf_get_index('FLNT'), flnt)

         ! Copy data
         flns(:ncol) = this%flux_net_bot(:ncol)
         flnt(:ncol) = this%flux_net_top(:ncol)

         ! Check that we set them right
         call assert_valid(flns(:ncol), trim(sub_name) // ': flns')
         call assert_valid(flnt(:ncol), trim(sub_name) // ': flnt')
      else
         call endrun(trim(sub_name) // ': flux_type ' // this%flux_type // ' not known.')
      end if

   end subroutine cam_fluxes_to_pbuf


   subroutine cam_fluxes_initialize(this, ncol, nlevels, nbands, ktop, kbot, flux_type)
      class(cam_flux_type), intent(inout) :: this
      integer, intent(in) :: ncol, nlevels, nbands, ktop, kbot
      character(len=*), intent(in) :: flux_type

      this%ktop = ktop
      this%kbot = kbot
      this%flux_type = flux_type

      ! allocate CAM fluxes
      allocate(this%flux_up(ncol,nlevels), &
               this%flux_dn(ncol,nlevels), &
               this%flux_net(ncol,nlevels), &
               this%bnd_flux_up(ncol,nlevels,nbands), &
               this%bnd_flux_dn(ncol,nlevels,nbands), &
               this%bnd_flux_net(ncol,nlevels,nbands), &
               this%flux_net_bot(ncol), &
               this%flux_net_top(ncol))

      ! allocate CAM heating rate
      allocate(this%heating_rate(ncol,nlevels-1))

      ! initialize to zero
      this%flux_up = 0.0
      this%flux_dn = 0.0
      this%flux_net = 0.0
      this%bnd_flux_up = 0.0
      this%bnd_flux_dn = 0.0
      this%bnd_flux_net = 0.0
      this%flux_net_bot = 0.0
      this%flux_net_top = 0.0
      this%heating_rate = 0.0

   end subroutine cam_fluxes_initialize

   subroutine cam_fluxes_finalize(this)
      class(cam_flux_type), intent(inout) :: this
      deallocate(this%flux_up, this%flux_dn, this%flux_net, &
                 this%bnd_flux_up, this%bnd_flux_dn, this%bnd_flux_net, &
                 this%flux_net_bot, this%flux_net_top, &
                 this%heating_rate)
   end subroutine cam_fluxes_finalize

   subroutine cam_fluxes_aggregate_average(this, that)
      class(cam_flux_type), intent(inout) :: this
      class(cam_flux_type), intent(in) :: that

      ! Average fluxes from two flux type instances
      this%flux_up = (this%flux_up + that%flux_up) * 0.5
      this%flux_dn = (this%flux_dn + that%flux_dn) * 0.5
      this%flux_net = (this%flux_net + that%flux_net) * 0.5
      this%bnd_flux_up = (this%bnd_flux_up + that%bnd_flux_up) * 0.5
      this%bnd_flux_dn = (this%bnd_flux_dn + that%bnd_flux_dn) * 0.5
      this%bnd_flux_net = (this%bnd_flux_net + that%bnd_flux_net) * 0.5
      this%flux_net_bot = (this%flux_net_bot + that%flux_net_bot) * 0.5
      this%flux_net_top = (this%flux_net_top + that%flux_net_top) * 0.5
      this%heating_rate = (this%heating_rate + that%heating_rate) * 0.5
   end subroutine cam_fluxes_aggregate_average

end module cam_fluxes
