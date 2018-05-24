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
! Description:  Unit test for gas optical depth calculation.
!

program test_gas_optics
  use mo_rte_kind,   only: wp
  use mo_test_files_io, only: read_atmos, read_lw_bc, &
                              write_spectral_disc, write_optical_props, write_direction, &
                              write_lw_Planck_sources, write_sw_solar_sources
  use mo_gas_concentrations, &
                        only: ty_gas_concs
  use mo_optical_props, only: ty_optical_props_arry, &
                              ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
  use mo_gas_optics, only: ty_gas_optics_specification
  use mo_load_coefficients,        only: load_and_init

  implicit none

  real(wp), dimension(:,:),   allocatable :: p_lay, t_lay, p_lev, t_lev
  real(wp), dimension(:),     allocatable :: t_sfc
  real(wp), dimension(:,:),   allocatable :: emis_sfc
  real(wp), dimension(:,:),   allocatable :: col_dry

  class (ty_optical_props_arry), &
                              allocatable :: optical_props, optical_props_subset
  type(ty_gas_optics_specification)       :: kdist
  type(ty_gas_concs)                      :: gas_concs, gas_concs_subset

  ! Optional fields for internal source gas optics
  !
  real(wp), dimension(:,:  ), allocatable :: sfc_src
  real(wp), dimension(:,:,:), allocatable :: lay_src, lev_src_inc, lev_src_dec
  !
  ! Optional fields for external source gas optics
  !
  real(wp), dimension(:,:  ), allocatable :: toa_src ! TOA incident radiation  (ncol, ngpt)

  integer, parameter :: blockSize = 4, nmom = 3
  ! dimensions
  integer :: ncol, nlay, ngpt, nbnd
  integer :: b, items_to_write, nBlocks, colS, colE
  logical :: top_at_1

  character(len=64), parameter :: fileName = 'rrtmgp-inputs-outputs.nc'
  ! ==========================================================
  ! get started
  print *, 'test start ...'
  print *, ' Using blocks of size ', blockSize

  ! load profile
  call read_atmos(fileName, p_lay, t_lay, p_lev, t_lev, gas_concs, col_dry)

  ncol = size(p_lay, 1)
  nlay = size(p_lay, 2)

  call load_and_init(kdist,'coefficients.nc', gas_concs)
  ngpt = kdist%get_ngpt()
  nbnd = kdist%get_nband()

  allocate(sfc_src(ncol, ngpt))
  allocate(lay_src(ncol, nlay, ngpt))
  allocate(lev_src_inc(ncol, nlay, ngpt))
  allocate(lev_src_dec(ncol, nlay, ngpt))
  allocate(toa_src(ncol, ngpt))

  if(kdist%is_internal_source_present()) then
    print *, " Compute longwave gas optical depths"
    allocate(ty_optical_props_1scl::optical_props)
    allocate(ty_optical_props_1scl::optical_props_subset)
    call read_lw_bc(fileName, t_sfc, emis_sfc)
  else
    print *, " Compute shortwave gas optical depths"
    allocate(ty_optical_props_2str::optical_props)
    allocate(ty_optical_props_2str::optical_props_subset)
  end if

  select type (optical_props)
    type is (ty_optical_props_1scl) ! two-stream calculation
        call stop_on_err(optical_props%init_1scl(ncol, nlay, ngpt))
    type is (ty_optical_props_2str) ! two-stream calculation
        call stop_on_err(optical_props%init_2str(ncol, nlay, ngpt))
    type is (ty_optical_props_nstr)
        call stop_on_err(optical_props%init_nstr(nmom, ncol, nlay, ngpt))
  end select

  nBlocks = ncol/blockSize ! Integer division
  print *, "Doing ", nBlocks, "blocks"
  do b = 1, nBlocks
    colS = (b-1) * blockSize + 1
    colE =  b    * blockSize
    call stop_on_err(    gas_concs%get_subset(colS, colE-colS+1, gas_concs_subset))
    call stop_on_err(optical_props%get_subset(colS, colE-colS+1, optical_props_subset))
    !
    ! The if statement for the presence of col_dry makes this look long, but I expect most
    !   users will decide ahead of time whether they will compute col_dry themselves or
    !   let us do so
    !
    if(kdist%is_internal_source_present()) then
      if(allocated(col_dry)) then
        call stop_on_err(kdist%gas_optics(p_lay(colS:colE,:), p_lev(colS:colE,:), &
                                          t_lay(colS:colE,:), t_sfc(colS:colE  ), &
                                          gas_concs_subset,                       &
                                          optical_props_subset,                   &
                                          lay_src    (colS:colE,:,:), &
                                          lev_src_inc(colS:colE,:,:), &
                                          lev_src_dec(colS:colE,:,:), &
                                          sfc_src    (colS:colE,  :), &
                                          tlev=t_lev(colS:colE,:  ), &
                                          col_dry=col_dry(colS:colE, :)) )
      else
        call stop_on_err(kdist%gas_optics(p_lay(colS:colE,:), p_lev(colS:colE,:), &
                                          t_lay(colS:colE,:), t_sfc(colS:colE  ), &
                                          gas_concs_subset,                       &
                                          optical_props_subset,                   &
                                          lay_src    (colS:colE,:,:), &
                                          lev_src_inc(colS:colE,:,:), &
                                          lev_src_dec(colS:colE,:,:), &
                                          sfc_src    (colS:colE,  :), &
                                          tlev=t_lev(colS:colE,:  ))  )
      end if
    else
      if(allocated(col_dry)) then
        call stop_on_err(kdist%gas_optics(p_lay(colS:colE,:), p_lev(colS:colE,:), &
                                          t_lay(colS:colE,:),                     &
                                          gas_concs_subset,                       &
                                          optical_props_subset,                   &
                                          toa_src        (colS:colE,:), &
                                          col_dry=col_dry(colS:colE, :)) )
      else
        call stop_on_err(kdist%gas_optics(p_lay(colS:colE,:), p_lev(colS:colE,:), &
                                          t_lay(colS:colE,:),                     &
                                          gas_concs_subset,                       &
                                          optical_props_subset,                   &
                                          toa_src    (colS:colE,:)) )
      end if
    end if
    call stop_on_err(assign_subset(optical_props_subset, colS, colE, optical_props))
  end do

  if(mod(ncol, blockSize) /= 0) then
    colS = ncol/blockSize * blockSize + 1  ! Integer arithmetic
    colE = ncol
    print *, "Doing ", colE-colS+1, "extra columns"
    call stop_on_err(    gas_concs%get_subset(colS, colE-colS+1, gas_concs_subset))
    call stop_on_err(optical_props%get_subset(colS, colE-colS+1, optical_props_subset))
     if(kdist%is_internal_source_present()) then
      if(allocated(col_dry)) then
        call stop_on_err(kdist%gas_optics(p_lay(colS:colE,:), p_lev(colS:colE,:), &
                                          t_lay(colS:colE,:), t_sfc(colS:colE  ), &
                                          gas_concs_subset,                       &
                                          optical_props_subset,                   &
                                          lay_src    (colS:colE,:,:), &
                                          lev_src_inc(colS:colE,:,:), &
                                          lev_src_dec(colS:colE,:,:), &
                                          sfc_src    (colS:colE,:  ), &
                                          tlev=t_lev(colS:colE,:  ), &
                                          col_dry=col_dry(colS:colE, :)) )
      else
        call stop_on_err(kdist%gas_optics(p_lay(colS:colE,:), p_lev(colS:colE,:), &
                                          t_lay(colS:colE,:), t_sfc(colS:colE  ), &
                                          gas_concs_subset,                       &
                                          optical_props_subset,                   &
                                          lay_src    (colS:colE,:,:), &
                                          lev_src_inc(colS:colE,:,:), &
                                          lev_src_dec(colS:colE,:,:), &
                                          sfc_src    (colS:colE,:  ), &
                                          tlev=t_lev(colS:colE,:  ))  )
      end if
    else
      if(allocated(col_dry)) then
        call stop_on_err(kdist%gas_optics(p_lay(colS:colE,:), p_lev(colS:colE,:), &
                                          t_lay(colS:colE,:),                     &
                                          gas_concs_subset,                       &
                                          optical_props_subset,                   &
                                          toa_src        (colS:colE,:), &
                                          col_dry=col_dry(colS:colE, :)) )
      else
        call stop_on_err(kdist%gas_optics(p_lay(colS:colE,:), p_lev(colS:colE,:), &
                                          t_lay(colS:colE,:),                     &
                                          gas_concs_subset,                       &
                                          optical_props_subset,                   &
                                          toa_src    (colS:colE,:)) )
      end if
    end if
    call stop_on_err(assign_subset(optical_props_subset, colS, colE, optical_props))
  end if

  !
  ! Write fields out
  !
  call write_spectral_disc(fileName, kdist)
  call write_optical_props(fileName, optical_props)
  top_at_1 = p_lay(1, 1) < p_lay(1, nlay)
  call write_direction(fileName, top_at_1)
  if(kdist%is_internal_source_present()) then
    call write_lw_Planck_sources(fileName, lay_src, lev_src_inc, lev_src_dec, sfc_src)
  else
    call write_sw_solar_sources(fileName, toa_src)
  end if


  ! all done
  print *, 'test end.'

contains
! -----------------------------------------------------------------------------------
  subroutine stop_on_err(error_msg)
    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: error_msg

    if(error_msg /= "") then
      write (error_unit,*) trim(error_msg)
      write (error_unit,*) "test_gas_optics stopping"
      stop
    end if

  end subroutine stop_on_err
! -----------------------------------------------------------------------------------
!
! Assign optical properties from a set of columns to a specified position in a larger set
!   Could be bound to ty_optical_props but this would require more careful error checking and
!   it doesn't seem likely users will need this
!
  function assign_subset(subset, colS, colE, full) result(error_msg)
    class(ty_optical_props_arry), intent(in   ) :: subset
    integer,                      intent(in   ) :: colS, colE
    class(ty_optical_props_arry), intent(inout) :: full
    character(len=128)                          :: error_msg

    real(wp), dimension(colE-colS+1, size(subset%tau,2), size(subset%tau,3)) :: tau, g, ssa
    real(wp), dimension(:,:,:,:), allocatable :: p
    integer :: nmom
    ! ------------------------------------------
    error_msg = ""
    if(colS > size(full%tau,1) .or.  colE > size(full%tau,1) .or. &
       colS < 1 .or.  colE <1) then
      error_msg = "  Subset, colS, colE not consistent with full optical properties arrays???"
      return
    end if

    nmom = 0
    select type (subset)
      type is (ty_optical_props_nstr) ! n-stream calculation
        nmom = size(subset%p,1)
        allocate(p(nmom, colE-colS+1, nlay, ngpt))
    end select

    full%tau(colS:colE,:,:) = subset%tau
    ! For whatever reason the Intel compiler, at least, can't tell that full and subset
    !   have the same type, so we copy to intermediate storage.
    select type (subset)
      type is (ty_optical_props_2str) ! two-stream calculation
        ssa = subset%ssa
        g   = subset%g
      type is (ty_optical_props_nstr) ! n-stream calculation
        ssa = subset%ssa
        p   = subset%p
    end select
    select type (full)
      type is (ty_optical_props_2str) ! two-stream calculation
        full%ssa(colS:colE,:,:) = ssa
        full%g  (colS:colE,:,:) = g
      type is (ty_optical_props_nstr) ! two-stream calculation
        full%ssa(colS:colE,:,:) = ssa
        full%p(1:nmom, &
                 colS:colE,:,:) = p
    end select

  end function assign_subset
! -----------------------------------------------------------------------------------
end program test_gas_optics
