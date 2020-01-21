module cam_optics

   use shr_kind_mod, only: r8=>shr_kind_r8
   use assertions, only: assert, assert_valid, assert_range
   use radiation_state, only: nlev_rad, ktop, kbot
   use radiation_utils, only: handle_error
   use radconstants, only: nswbands, nswgpts, &
                           nlwbands, nlwgpts
   use perf_mod, only: t_startf, t_stopf

   implicit none
   private

   public ::                 &
      free_optics_sw,        &
      free_optics_lw,        &
      set_cloud_optics_sw,   &
      set_cloud_optics_lw,   &
      set_aerosol_optics_sw, &
      set_aerosol_optics_lw

   ! Mapping from old RRTMG sw bands to new band ordering in RRTMGP
   integer, dimension(14) :: map_rrtmg_to_rrtmgp_swbands = (/ &
      14, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 &
   /)
   integer, dimension(14) :: map_rrtmgp_to_rrtmg_swbands = (/ &
      2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 1 &
   /)

contains

   subroutine get_cloud_optics_sw(cld, cldfsnow, iclwp, iciwp, icswp, &
                                  rel, rei, dei, des, lambdac, mu, &
                                  tau_out, ssa_out, asm_out)

      use physics_types, only: physics_state
      use physics_buffer, only: physics_buffer_desc, pbuf_get_field, &
                                pbuf_get_index, pbuf_old_tim_idx
      use cloud_rad_props, only: mitchell_ice_optics_sw, &
                                 gammadist_liquid_optics_sw
      use ebert_curry, only: ec_ice_optics_sw
      use slingo, only: slingo_liq_optics_sw
      use phys_control, only: phys_getopts
      use cam_abortutils, only: endrun
      use rad_constituents, only: icecldoptics, liqcldoptics

      ! Pointers to fields on the physics buffer
      real(r8), intent(in), dimension(:,:) :: &
         iciwp, iclwp, icswp, rei, rel, dei, des, cld, cldfsnow, lambdac, mu

      ! Outputs are shortwave cloud optical properties *by band*. Dimensions should
      ! be ncol,pver,nswbands. Generally, this should be able to handle cases were
      ! ncol might be something like nday, and pver could be arbitrary so long as
      ! corresponding pbuf/state fields were defined for all indices of pver. This
      ! isn't the case right now I don't think, as cloud_rad_props makes explicit
      ! assumptions about array sizes.
      real(r8), intent(inout), dimension(:,:,:) :: tau_out, ssa_out, asm_out

      ! Temporary variables to hold cloud optical properties before combining into
      ! output arrays. Same shape as output arrays, so get shapes from output.
      real(r8), dimension(nswbands,size(cld, 1),size(cld, 2)) :: &
            liq_tau, liq_tau_ssa, liq_tau_ssa_g, liq_tau_ssa_f, &
            ice_tau, ice_tau_ssa, ice_tau_ssa_g, ice_tau_ssa_f, &
            snow_tau, snow_tau_ssa, snow_tau_ssa_g, snow_tau_ssa_f, &
            cloud_tau, cloud_tau_ssa, cloud_tau_ssa_g, cloud_tau_ssa_f, &
            combined_tau, combined_tau_ssa, combined_tau_ssa_g, combined_tau_ssa_f, &
            tau_tmp, tau_ssa_tmp, tau_ssa_g_tmp, tau_ssa_f_tmp

      integer :: ncol, nlev, ibnd, icol, ilev


      ! Initialize
      ice_tau = 0
      ice_tau_ssa = 0
      ice_tau_ssa_g = 0
      ice_tau_ssa_f = 0
      liq_tau = 0
      liq_tau_ssa = 0
      liq_tau_ssa_g = 0
      liq_tau_ssa_f = 0
      snow_tau = 0
      snow_tau_ssa = 0
      snow_tau_ssa_g = 0
      snow_tau_ssa_f = 0
      combined_tau = 0
      combined_tau_ssa = 0
      combined_tau_ssa_g = 0
      combined_tau_ssa_f = 0

      ncol = size(tau_out, 1)
      nlev = size(tau_out, 2)

      ! Get ice cloud optics
      if (trim(icecldoptics) == 'mitchell') then
         call mitchell_ice_optics_sw(ncol, iciwp, dei, &
                                tau_tmp, tau_ssa_tmp, &
                                tau_ssa_g_tmp, tau_ssa_f_tmp)

         ! Mitchell optics hard-coded for RRTMG band-ordering
         do ibnd = 1,nswbands
            do ilev = 1,nlev
               do icol = 1,ncol
                  ice_tau      (ibnd,icol,ilev) = tau_tmp      (map_rrtmg_to_rrtmgp_swbands(ibnd),icol,ilev)
                  ice_tau_ssa  (ibnd,icol,ilev) = tau_ssa_tmp  (map_rrtmg_to_rrtmgp_swbands(ibnd),icol,ilev)
                  ice_tau_ssa_g(ibnd,icol,ilev) = tau_ssa_g_tmp(map_rrtmg_to_rrtmgp_swbands(ibnd),icol,ilev)
                  ice_tau_ssa_f(ibnd,icol,ilev) = tau_ssa_f_tmp(map_rrtmg_to_rrtmgp_swbands(ibnd),icol,ilev)
               end do
            end do
         end do
      else if (trim(icecldoptics) == 'ebertcurry') then
         call ec_ice_optics_sw( &
            ncol, iciwp, rei, cld, &
            ice_tau, ice_tau_ssa, ice_tau_ssa_g, ice_tau_ssa_f &
         )
            
      else
         call endrun('Ice optics scheme ' // trim(icecldoptics) // ' not recognized.')
      end if
      
      ! Get liquid cloud optics
      if (trim(liqcldoptics) == 'gammadist') then
         call gammadist_liquid_optics_sw( &
            ncol, iclwp, lambdac, mu, &
            tau_tmp, tau_ssa_tmp, &
            tau_ssa_g_tmp, tau_ssa_f_tmp &
         )

         ! Gammadist optics hard-coded for RRTMG band-ordering
         do ibnd = 1,nswbands
            do ilev = 1,nlev
               do icol = 1,ncol
                  liq_tau      (ibnd,icol,ilev) = tau_tmp      (map_rrtmg_to_rrtmgp_swbands(ibnd),icol,ilev)
                  liq_tau_ssa  (ibnd,icol,ilev) = tau_ssa_tmp  (map_rrtmg_to_rrtmgp_swbands(ibnd),icol,ilev)
                  liq_tau_ssa_g(ibnd,icol,ilev) = tau_ssa_g_tmp(map_rrtmg_to_rrtmgp_swbands(ibnd),icol,ilev)
                  liq_tau_ssa_f(ibnd,icol,ilev) = tau_ssa_f_tmp(map_rrtmg_to_rrtmgp_swbands(ibnd),icol,ilev)
               end do
            end do
         end do
      else if (trim(liqcldoptics) == 'slingo') then
         call slingo_liq_optics_sw( &
            ncol, iclwp, rel, cld, &
            liq_tau, liq_tau_ssa, liq_tau_ssa_g, liq_tau_ssa_f &
         )
      else
         call endrun('Liquid optics scheme ' // trim(liqcldoptics) // ' not recognized.')
      end if 

      ! Combine all cloud optics from CAM routines
      cloud_tau = ice_tau + liq_tau
      cloud_tau_ssa = ice_tau_ssa + liq_tau_ssa
      cloud_tau_ssa_g = ice_tau_ssa_g + liq_tau_ssa_g

      ! Get snow cloud optics?
      if (do_snow_optics()) then
         ! Doing snow optics; call procedure to get these from CAM state and pbuf
         call mitchell_ice_optics_sw( &
            ncol, icswp, des, &
            tau_tmp, tau_ssa_tmp, &
            tau_ssa_g_tmp, tau_ssa_f_tmp &
         )
         ! Mitchell optics hard-coded for RRTMG band-ordering
         do ibnd = 1,nswbands
            do ilev = 1,nlev
               do icol = 1,ncol
                  snow_tau      (ibnd,icol,ilev) = tau_tmp      (map_rrtmg_to_rrtmgp_swbands(ibnd),icol,ilev)
                  snow_tau_ssa  (ibnd,icol,ilev) = tau_ssa_tmp  (map_rrtmg_to_rrtmgp_swbands(ibnd),icol,ilev)
                  snow_tau_ssa_g(ibnd,icol,ilev) = tau_ssa_g_tmp(map_rrtmg_to_rrtmgp_swbands(ibnd),icol,ilev)
                  snow_tau_ssa_f(ibnd,icol,ilev) = tau_ssa_f_tmp(map_rrtmg_to_rrtmgp_swbands(ibnd),icol,ilev)
               end do
            end do
         end do

         call combine_properties( &
            cld, cloud_tau, cldfsnow, snow_tau, combined_tau &
         )
         call combine_properties( &
            cld, cloud_tau_ssa, cldfsnow, snow_tau_ssa, combined_tau_ssa &
         )
         call combine_properties( &
            cld, cloud_tau_ssa_g, cldfsnow, snow_tau_ssa_g, combined_tau_ssa_g &
         )
      else
         combined_tau = cloud_tau
         combined_tau_ssa = cloud_tau_ssa
         combined_tau_ssa_g = cloud_tau_ssa_g
      end if
     
      ! Copy to output arrays, converting to optical depth, single scattering
      ! albedo, and assymmetry parameter from the products that the CAM routines
      ! return. Make sure we do not try to divide by zero.
      tau_out = 0
      ssa_out = 0
      asm_out = 0
      do ibnd = 1,nswbands
         do ilev = 1,nlev
            do icol = 1,ncol
               tau_out(icol,ilev,ibnd) = combined_tau(ibnd,icol,ilev)
               if (combined_tau(ibnd,icol,ilev) > 0) then
                  ssa_out(icol,ilev,ibnd) = combined_tau_ssa(ibnd,icol,ilev) / combined_tau(ibnd,icol,ilev)
               end if
               if (combined_tau_ssa(ibnd,icol,ilev) > 0) then
                  asm_out(icol,ilev,ibnd) = combined_tau_ssa_g(ibnd,icol,ilev) / combined_tau_ssa(ibnd,icol,ilev)
               end if
            end do
         end do
      end do
                 
   end subroutine get_cloud_optics_sw

   !----------------------------------------------------------------------------

   ! Provide a procedure to combine cloud optical properties by weighting
   ! contributions by fraction present. I.e., for combining cloud and snow
   ! optical properties, we weight the cloud and snow properties by the fraction
   ! of cloud and snow present.
   subroutine combine_properties(fraction1, property1, &
                                 fraction2, property2, &
                                 combined_property)
      
      ! Input fractions of each type of constituent
      real(r8), intent(in) :: fraction1(:,:), fraction2(:,:)

      ! Individual optical properties for each constituent
      real(r8), intent(in) :: property1(:,:,:), property2(:,:,:)

      ! Combined optical property
      real(r8), intent(out) :: combined_property(:,:,:)

      ! Combined fraction (max of property 1 and 2)
      real(r8) :: combined_fraction

      ! Loop variables
      integer :: nbnd, ncol, nlev, ibnd, icol, ilev

      ! Dimension sizes
      nbnd = size(property1, 1)
      ncol = min(size(property1, 2), size(combined_property, 2))
      nlev = size(property1, 3)

      ! Combine optical properties by weighting by amount of cloud and snow
      !$acc parallel loop collapse(3)
      do ilev = 1,nlev
         do icol = 1,ncol
            do ibnd = 1,nbnd
               combined_fraction = max(fraction1(icol,ilev), fraction2(icol,ilev))
               if (combined_fraction > 0) then
                  combined_property(ibnd,icol,ilev) = ( &
                     fraction1(icol,ilev) * property1(ibnd,icol,ilev) &
                   + fraction2(icol,ilev) * property2(ibnd,icol,ilev) &
                  ) / combined_fraction
               else
                  combined_property(ibnd,icol,ilev) = 0
               end if
            end do
         end do
      end do
   end subroutine combine_properties

   !----------------------------------------------------------------------------

   subroutine set_cloud_optics_sw(map_gpt_to_bnd, pmid, cld, cldfsnow, &
                                  iclwp, iciwp, icswp, rel, rei, dei, des, &
                                  lambdac, mu, &
                                  tau_gpt, ssa_gpt, asm_gpt)
      
      use ppgrid, only: pcols, pver, pverp
      use physics_types, only: physics_state
      use physics_buffer, only: physics_buffer_desc, &
                                pbuf_get_field, &
                                pbuf_get_index, pbuf_old_tim_idx
      use mo_optical_props, only: ty_optical_props_2str
      use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp
      use mcica_subcol_gen, only: mcica_subcol_mask

      integer, intent(in), dimension(:) :: map_gpt_to_bnd
      real(r8), intent(in), dimension(:,:) :: &
         pmid, iciwp, iclwp, icswp, rei, rel, dei, des, cld, cldfsnow, lambdac, mu
      real(r8), intent(inout), dimension(:,:,:) :: tau_gpt, ssa_gpt, asm_gpt

      ! For MCICA sampling routine
      integer, parameter :: changeseed = 1

      ! Dimension sizes
      integer :: ncol, nlev

      ! Loop variables
      integer :: icol, ilev, igpt, ibnd, ilev_cam, ilev_rad

      real(r8), allocatable, dimension(:,:) :: combined_cloud_fraction
      real(r8), allocatable, dimension(:,:,:) :: tau_bnd, ssa_bnd, asm_bnd
      logical, allocatable, dimension(:,:,:) :: iscloudy

      ! Set a name for this subroutine to write to error messages
      character(len=32) :: subname = 'set_cloud_optics_sw'

      ncol = size(pmid, 1)
      nlev = size(pmid, 2)

      call t_startf('optics_allocate_sw')
      allocate( &
         combined_cloud_fraction(ncol,nlev), &
         iscloudy(ncol,nlev,nswgpts), &
         tau_bnd(ncol,nlev,nswbands), &
         ssa_bnd(ncol,nlev,nswbands), &
         asm_bnd(ncol,nlev,nswbands) &
      )
      call t_stopf('optics_allocate_sw')

      ! Get optics by band
      call t_startf('optics_get_cloud_optics_sw')
      call get_cloud_optics_sw( &
         cld, cldfsnow, &
         iclwp, iciwp, icswp, &
         rel, rei, dei, des, lambdac, mu, &
         tau_bnd, ssa_bnd, asm_bnd &
      )
      call t_stopf('optics_get_cloud_optics_sw')

      ! Get combined cloud fraction (cloud + snow) for MCICA sampling
      call t_startf('optics_combine_cloud_fraction_sw')
      if (do_snow_optics()) then
         do ilev = 1,nlev
            do icol = 1,ncol
               combined_cloud_fraction(icol,ilev) = max(cld(icol,ilev), cldfsnow(icol,ilev))
            end do
         end do
      else
         do ilev = 1,nlev
            do icol = 1,ncol
               combined_cloud_fraction(icol,ilev) = cld(icol,ilev)
            end do
         end do
      end if
      call t_stopf('optics_combine_cloud_fraction_sw')

      ! Do MCICA sampling of optics here. This will map bands to gpoints,
      ! while doing stochastic sampling of cloud state
      call t_startf('optics_mcica_subcol_mask_sw')
      call mcica_subcol_mask( &
         nswgpts, ncol, nlev, changeseed, &
         pmid, combined_cloud_fraction, iscloudy &
      )
      call t_stopf('optics_mcica_subcol_mask_sw')

      ! Generate subcolumns for clouds
      call t_startf('optics_sample_optics_sw')
      do igpt = 1,size(tau_gpt, 3)
         do ilev = 1,size(tau_gpt, 2)
            do icol = 1,size(tau_gpt, 1)
               tau_gpt(icol,ilev,igpt) = 0
               ssa_gpt(icol,ilev,igpt) = 1
               asm_gpt(icol,ilev,igpt) = 0
            end do
         end do
      end do
      call sample_optics_sw( &
         map_gpt_to_bnd, iscloudy, &
         tau_bnd, ssa_bnd, asm_bnd, &
         tau_gpt, ssa_gpt, asm_gpt &
      )
      call t_stopf('optics_sample_optics_sw')

      call t_startf('optics_deallocate_sw')
      deallocate(combined_cloud_fraction, iscloudy, tau_bnd, ssa_bnd, asm_bnd)
      call t_stopf('optics_deallocate_sw')

   end subroutine set_cloud_optics_sw

   !----------------------------------------------------------------------------

   subroutine sample_optics_sw(gpt2bnd, iscloudy, &
                               tau_bnd, ssa_bnd, asm_bnd, &
                               tau_gpt, ssa_gpt, asm_gpt)
      integer, intent(in), dimension(:) :: gpt2bnd
      logical, intent(in), dimension(:,:,:) :: iscloudy
      real(r8), intent(in), dimension(:,:,:) :: tau_bnd, ssa_bnd, asm_bnd
      real(r8), intent(inout), dimension(:,:,:) :: tau_gpt, ssa_gpt, asm_gpt
      integer :: nlev1, nlev2, ncol, ngpt, igpt, ibnd, icol, ilev1, ilev2

      ! Loop over g-points and map bands to g-points; each subcolumn
      ! corresponds to a single g-point. This is how this code implements the
      ! McICA assumptions: simultaneously sampling over cloud state and
      ! g-point.
      ncol  = size(tau_gpt, 1)
      nlev1 = size(tau_bnd, 2)
      nlev2 = size(tau_gpt, 2)
      ngpt  = size(tau_gpt, 3)
      !$acc parallel loop collapse(3) private(ilev2)
      do igpt = 1,ngpt
         do ilev1 = 1,nlev1
            do icol = 1,ncol
               ! Outputs may have larger vertical dimension
               ilev2 = ilev1 + (nlev2 - nlev1)
               if (iscloudy(icol,ilev1,igpt)) then
                  tau_gpt(icol,ilev2,igpt) = tau_bnd(icol,ilev1,gpt2bnd(igpt))
                  ssa_gpt(icol,ilev2,igpt) = ssa_bnd(icol,ilev1,gpt2bnd(igpt))
                  asm_gpt(icol,ilev2,igpt) = asm_bnd(icol,ilev1,gpt2bnd(igpt))
               else
                  tau_gpt(icol,ilev2,igpt) = 0._r8
                  ssa_gpt(icol,ilev2,igpt) = 1._r8
                  asm_gpt(icol,ilev2,igpt) = 0._r8
               end if  ! iscloudy
            end do  ! icol
         end do  ! ilev
      end do  ! igpt
   end subroutine sample_optics_sw

   !----------------------------------------------------------------------------

   subroutine set_cloud_optics_lw(map_gpt_to_bnd, pmid, cld, cldfsnow, &
                                  iclwp, iciwp, icswp, rei, dei, des, &
                                  lambdac, mu, tau_out)
      
      use physics_types, only: physics_state
      use physics_buffer, only: physics_buffer_desc, &
                                pbuf_get_field, pbuf_get_index, pbuf_old_tim_idx
      use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp
      use mcica_subcol_gen, only: mcica_subcol_mask
      use cloud_rad_props, only: gammadist_liquid_optics_lw, &
                                 mitchell_ice_optics_lw
      use ebert_curry, only: ec_ice_optics_lw
      use slingo, only: slingo_liq_optics_lw
      use radconstants, only: nlwbands
      use cam_abortutils, only: endrun
      use rad_constituents, only: icecldoptics, liqcldoptics

      integer, intent(in), dimension(:) :: map_gpt_to_bnd
      real(r8), intent(in), dimension(:,:) :: &
         pmid, rei, dei, des, iclwp, iciwp, icswp, cld, cldfsnow, lambdac, mu

      real(r8), intent(inout), dimension(:,:,:) :: tau_out

      real(r8), allocatable :: combined_cloud_fraction(:,:)

      ! For MCICA sampling routine
      integer, parameter :: changeseed = 1

      ! Temporary arrays to hold mcica-sampled cloud optics
      logical, allocatable :: iscloudy(:,:,:)

      real(r8), allocatable, dimension(:,:,:) :: tau_bnd

      ! Temporary variables to hold absorption optical depth
      real(r8), allocatable, dimension(:,:,:) :: &
            ice_tau, liq_tau, snow_tau, cloud_tau, combined_tau

      ! Loop variables
      integer :: ncol, nlev, icol, ilev, igpt, ibnd

      ncol = size(pmid, 1)
      nlev = size(pmid, 2)

      allocate(combined_cloud_fraction(ncol,nlev))
      allocate(iscloudy(ncol,nlev,nlwgpts))
      allocate(tau_bnd(ncol,nlev,nlwbands))
      allocate( &
         ice_tau(nlwbands,ncol,nlev), liq_tau(nlwbands,ncol,nlev), &
         snow_tau(nlwbands,ncol,nlev), cloud_tau(nlwbands,ncol,nlev), &
         combined_tau(nlwbands,ncol,nlev) &
      )

      ! initialize
      !$acc parallel loop collapse(3)
      do ilev = 1,nlev
         do icol = 1,ncol
            do ibnd = 1,nlwbands
               ice_tau(ibnd,icol,ilev) = 0.0
               liq_tau(ibnd,icol,ilev) = 0.0
               snow_tau(ibnd,icol,ilev) = 0.0
               cloud_tau(ibnd,icol,ilev) = 0.0
               combined_tau(ibnd,icol,ilev) = 0.0
            end do
         end do
      end do

      ! Get ice optics
      if (trim(icecldoptics) == 'mitchell') then
         call mitchell_ice_optics_lw(ncol, iciwp, dei, ice_tau)
      else if (trim(icecldoptics) == 'ebertcurry') then
         call ec_ice_optics_lw(ncol, iciwp, iclwp, rei, cld, ice_tau)
      else
         call endrun('Ice optics scheme ' // trim(icecldoptics) // ' not recognized.')
      end if

      ! Get liquid optics
      if (trim(liqcldoptics) == 'gammadist') then
         call gammadist_liquid_optics_lw(ncol, iclwp, lambdac, mu, liq_tau)
      else if (trim(liqcldoptics) == 'slingo') then
         call slingo_liq_optics_lw(ncol, iclwp, iciwp, rei, cld, liq_tau)
      else
         call endrun('Liquid optics scheme ' // trim(liqcldoptics) // ' not recognized.')
      end if

      ! Combine liquid and ice
      !$acc parallel loop collapse(3)
      do ilev = 1,nlev
         do icol = 1,ncol
            do ibnd = 1,nlwbands
               cloud_tau(ibnd,icol,ilev) = liq_tau(ibnd,icol,ilev) + ice_tau(ibnd,icol,ilev)
            end do
         end do
      end do

      ! Get snow optics?
      if (do_snow_optics()) then
         call mitchell_ice_optics_lw(ncol, icswp, des, snow_tau)
         call combine_properties(cld, cloud_tau, cldfsnow, snow_tau, combined_tau)
      else
         !$acc parallel loop collapse(3)
         do ilev = 1,nlev
            do icol = 1,ncol
               do ibnd = 1,nlwbands
                  combined_tau(ibnd,icol,ilev) = cloud_tau(ibnd,icol,ilev)
               end do
            end do
         end do
      end if

      ! Re-order array dimensions for outputs
      !$acc parallel loop collapse(3)
      do ibnd = 1,nlwbands
         do ilev = 1,nlev
            do icol = 1,ncol
               tau_bnd(icol,ilev,ibnd) = combined_tau(ibnd,icol,ilev)
            end do
         end do
      end do

      ! Combine cloud and snow fractions for MCICA sampling
      if (do_snow_optics()) then
         !$acc parallel loop collapse(2)
         do ilev = 1,nlev
            do icol = 1,ncol
               combined_cloud_fraction(icol,ilev) = max(cld(icol,ilev), cldfsnow(icol,ilev))
            end do
         end do  
      else
         !$acc parallel loop collapse(2)
         do ilev = 1,nlev
            do icol = 1,ncol
               combined_cloud_fraction(icol,ilev) = cld(icol,ilev)
            end do
         end do  
      end if

      ! Do MCICA sampling of optics here. This will map bands to gpoints,
      ! while doing stochastic sampling of cloud state
      !
      ! First, just get the stochastic subcolumn cloudy mask...
      call t_startf('mcica_subcol_mask_lw')
      call mcica_subcol_mask( &
         nlwgpts, ncol, nlev, changeseed, &
         pmid, combined_cloud_fraction, iscloudy &
      )
      call t_stopf('mcica_subcol_mask_lw')

      ! ... and now map optics to g-points, selecting a single subcolumn for each
      ! g-point. This implementation generates homogeneous clouds, but it would be
      ! straightforward to extend this to handle horizontally heterogeneous clouds
      ! as well.
      call sample_optics_lw(map_gpt_to_bnd, iscloudy, tau_bnd, tau_out)
     
      deallocate( &
         iscloudy, tau_bnd, combined_cloud_fraction, &
         ice_tau, liq_tau, snow_tau, cloud_tau, combined_tau &
      )

   end subroutine set_cloud_optics_lw

   !----------------------------------------------------------------------------

   ! Map bands to g-points; each subcolumn corresponds to a single g-point.
   ! This is how this code implements the McICA assumptions: simultaneously
   ! sampling over cloud state and g-point.
   subroutine sample_optics_lw(gpt2bnd, iscloudy, tau_bnd, tau_gpt)
      integer, intent(in), dimension(:) :: gpt2bnd
      logical, intent(in), dimension(:,:,:) :: iscloudy
      real(r8), intent(in), dimension(:,:,:) :: tau_bnd
      real(r8), intent(inout), dimension(:,:,:) :: tau_gpt
      integer :: nlev1, nlev2, ncol, ngpt, igpt, ibnd, icol, ilev1, ilev2

      ncol  = size(tau_gpt, 1)
      nlev1 = size(tau_bnd, 2)
      nlev2 = size(tau_gpt, 2)
      ngpt  = size(tau_gpt, 3)
      !$acc parallel loop collapse(3)
      do igpt = 1,ngpt
         do ilev2 = 1,nlev2
            do icol = 1,ncol
               tau_gpt(icol,ilev2,igpt) = 0
            end do  ! icol
         end do  ! ilev
      end do  ! igpt
      !$acc parallel loop collapse(3) private(ilev2)
      do igpt = 1,ngpt
         do ilev1 = 1,nlev1
            do icol = 1,ncol
               ! Outputs may have larger vertical dimension
               ilev2 = ilev1 + (nlev2 - nlev1)
               if (iscloudy(icol,ilev1,igpt)) then
                  tau_gpt(icol,ilev2,igpt) = tau_bnd(icol,ilev1,gpt2bnd(igpt))
               else
                  tau_gpt(icol,ilev2,igpt) = 0._r8
               end if  ! iscloudy
            end do  ! icol
         end do  ! ilev
      end do  ! igpt
   end subroutine sample_optics_lw

   !----------------------------------------------------------------------------

   subroutine set_aerosol_optics_lw(icall, state, pbuf, is_cmip6_volc, optics_out)
     
      use ppgrid, only: pcols, pver
      use physics_types, only: physics_state
      use physics_buffer, only: physics_buffer_desc, pbuf_get_index, &
                                pbuf_get_field, pbuf_old_tim_idx
      use aer_rad_props, only: aer_rad_props_lw
      use radconstants, only: nlwbands
      use mo_optical_props, only: ty_optical_props_1scl

      integer, intent(in) :: icall
      type(physics_state), intent(in) :: state
      type(physics_buffer_desc), pointer :: pbuf(:)
      logical, intent(in) :: is_cmip6_volc
      type(ty_optical_props_1scl), intent(inout) :: optics_out
      real(r8) :: tau(pcols,pver,nlwbands)
      integer :: ncol

      ! Subroutine name for error messages
      character(len=*), parameter :: subroutine_name = 'set_aerosol_optics_lw'

      ! Get aerosol absorption optical depth from CAM routine
      tau = 0.0
      call aer_rad_props_lw(is_cmip6_volc, icall, state, pbuf, tau)

      ! Populate the RRTMGP optical properties object with CAM optical depth
      ncol = size(optics_out%tau, 1)
      optics_out%tau(:,:,:) = 0.0
      optics_out%tau(1:ncol,ktop:kbot,1:nlwbands) = tau(1:ncol,1:pver,1:nlwbands)

   end subroutine set_aerosol_optics_lw

   !----------------------------------------------------------------------------

   subroutine set_aerosol_optics_sw(icall, state, pbuf, &
                                    night_indices, &
                                    is_cmip6_volc, &
                                    optics_out)

      use ppgrid, only: pcols, pver
      use physics_types, only: physics_state
      use physics_buffer, only: physics_buffer_desc
      use aer_rad_props, only: aer_rad_props_sw
      use mo_optical_props, only: ty_optical_props_2str

      integer, intent(in) :: icall
      type(physics_state), intent(in) :: state
      type(physics_buffer_desc), pointer :: pbuf(:)
      integer, intent(in) :: night_indices(:)
      logical, intent(in) :: is_cmip6_volc
      type(ty_optical_props_2str), intent(inout) :: optics_out

      real(r8) :: tau_out(pcols,pver,nswbands)
      real(r8) :: ssa_out(pcols,pver,nswbands)
      real(r8) :: asm_out(pcols,pver,nswbands)

      ! NOTE: aer_rad_props expects 0:pver indexing on these! It appears this is to
      ! account for the extra layer added above model top, but it is not entirely
      ! clear. This is not done for the longwave, and it is not really documented
      ! anywhere that I can find. Regardless, optical properties for the zero index
      ! are set to zero in aer_rad_props_sw as far as I can tell.
      !
      ! NOTE: dimension ordering is different than for cloud optics!
      real(r8), dimension(pcols,0:pver,nswbands) :: tau, tau_w, tau_w_g, tau_w_f

      ! Loop indices
      integer :: ncol, icol, ilev, ibnd, ilev_rad

      ncol = state%ncol

      ! Get aerosol absorption optical depth from CAM routine
      tau = 0.0
      tau_w = 0.0
      tau_w_g = 0.0
      tau_w_f = 0.0
      call aer_rad_props_sw(icall, state, pbuf, &
                            count(night_indices > 0), night_indices, is_cmip6_volc, &
                            tau, tau_w, tau_w_g, tau_w_f)

      ! Convert from products to optical properties
      tau_out = 0
      ssa_out = 0
      asm_out = 0
      do ibnd = 1,nswbands
         do ilev = 1,pver
            do icol = 1,ncol 
               tau_out(icol,ilev,ibnd) = tau(icol,ilev,ibnd)
               if (tau(icol,ilev,ibnd) > 0) then
                  ssa_out(icol,ilev,ibnd) = tau_w(icol,ilev,ibnd) / tau(icol,ilev,ibnd)
               end if
               if (tau_w(icol,ilev,ibnd) > 0) then
                  asm_out(icol,ilev,ibnd) = tau_w_g(icol,ilev,ibnd) / tau_w(icol,ilev,ibnd)
               end if
            end do
         end do
      end do

      ! Reset outputs (also handles case where radiation grid contains an extra
      ! layer above CAM grid)
      optics_out%tau = 0
      optics_out%ssa = 1
      optics_out%g = 0

      ! Set values
      ! We need to fix band ordering because the old input files assume RRTMG band
      ! ordering, but this has changed in RRTMGP.
      ! TODO: fix the input files themselves!
      do ibnd = 1,nswbands
         do icol = 1,ncol
            do ilev = 1,pver
               ilev_rad = ilev + (nlev_rad - pver)
               optics_out%tau(icol,ilev_rad,ibnd) = tau_out(icol,ilev,map_rrtmg_to_rrtmgp_swbands(ibnd))
               optics_out%ssa(icol,ilev_rad,ibnd) = ssa_out(icol,ilev,map_rrtmg_to_rrtmgp_swbands(ibnd))
               optics_out%g  (icol,ilev_rad,ibnd) = asm_out(icol,ilev,map_rrtmg_to_rrtmgp_swbands(ibnd))
            end do
         end do
      end do

      ! Check values
      call handle_error(optics_out%validate())

   end subroutine set_aerosol_optics_sw

   !----------------------------------------------------------------------------

   subroutine free_optics_sw(optics)
      use mo_optical_props, only: ty_optical_props_2str
      type(ty_optical_props_2str), intent(inout) :: optics
      if (allocated(optics%tau)) deallocate(optics%tau)
      if (allocated(optics%ssa)) deallocate(optics%ssa)
      if (allocated(optics%g)) deallocate(optics%g)
      call optics%finalize()
   end subroutine free_optics_sw

   !----------------------------------------------------------------------------

   subroutine free_optics_lw(optics)
      use mo_optical_props, only: ty_optical_props_1scl
      type(ty_optical_props_1scl), intent(inout) :: optics
      if (allocated(optics%tau)) deallocate(optics%tau)
      call optics%finalize()
   end subroutine free_optics_lw

   !----------------------------------------------------------------------------

   subroutine output_cloud_optics_sw(state, tau, ssa, asm)

      use ppgrid, only: pver
      use physics_types, only: physics_state
      use cam_history, only: outfld
      use radconstants, only: idx_sw_diag

      type(physics_state), intent(in) :: state
      real(r8), intent(in), dimension(:,:,:) :: tau, ssa, asm
      character(len=*), parameter :: subname = 'output_cloud_optics_sw'

      ! Send outputs to history buffer
      call outfld('CLOUD_TAU_SW', &
                  tau(1:state%ncol,1:pver,1:nswbands), &
                  state%ncol, state%lchnk)
      call outfld('CLOUD_SSA_SW', &
                  ssa(1:state%ncol,1:pver,1:nswbands), &
                  state%ncol, state%lchnk)
      call outfld('CLOUD_G_SW', &
                  asm(1:state%ncol,1:pver,1:nswbands), &
                  state%ncol, state%lchnk)
      call outfld('TOT_ICLD_VISTAU', &
                  tau(1:state%ncol,1:pver,idx_sw_diag), &
                  state%ncol, state%lchnk)

   end subroutine output_cloud_optics_sw

   !----------------------------------------------------------------------------

   subroutine output_cloud_optics_lw(state, tau)

      use ppgrid, only: pver
      use physics_types, only: physics_state
      use cam_history, only: outfld

      type(physics_state), intent(in) :: state
      real(r8), intent(in), dimension(:,:,:) :: tau

      ! Output
      call outfld('CLOUD_TAU_LW', &
                  tau(1:state%ncol,1:pver,1:nlwbands), &
                  state%ncol, state%lchnk)

   end subroutine output_cloud_optics_lw

   !----------------------------------------------------------------------------

   ! Should we do snow optics? Check for existence of "cldfsnow" variable
   logical function do_snow_optics()
      use physics_buffer, only: physics_buffer_desc, pbuf_get_index
      use phys_control, only: phys_getopts
      use cam_abortutils, only: endrun
      real(r8), pointer :: pbuf(:)
      integer :: err, idx
      logical :: use_SPCAM
      character(len=16) :: SPCAM_microp_scheme

      idx = pbuf_get_index('CLDFSNOW', errcode=err)
      if (idx > 0) then
         do_snow_optics = .true.
      else
         do_snow_optics = .false.
      end if

      ! Reset to false if using SPCAM with 1-mom scheme
      call phys_getopts(use_SPCAM_out           = use_SPCAM          )
      call phys_getopts(SPCAM_microp_scheme_out = SPCAM_microp_scheme)
      if (use_SPCAM .and. (trim(SPCAM_microp_scheme) == 'sam1mom')) then
         do_snow_optics = .false.
      end if

      return
   end function do_snow_optics

end module cam_optics
