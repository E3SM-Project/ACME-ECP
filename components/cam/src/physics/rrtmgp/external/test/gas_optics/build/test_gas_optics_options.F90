! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2018,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!
! Description:  Unit test for gas optical depth calculation.
!

program test_gas_optics
  use mo_rte_kind,      only: wp
  use mo_test_files_io, only: read_atmos, write_atmos
  use mo_gas_concentrations, &
                        only: ty_gas_concs
  use mo_gas_optics,    only: get_col_dry

  implicit none

  real(wp), dimension(:,:),   allocatable :: p_lay, t_lay, p_lev, t_lev
  real(wp), dimension(:),     allocatable :: t_sfc
  real(wp), dimension(:,:),   allocatable :: emis_sfc
  real(wp), dimension(:,:),   allocatable :: col_dry, vmrh2o

  type(ty_gas_concs)                      :: gas_concs, gas_concs_subset

  ! dimensions
  integer :: icol, ilay, ncol, nlay

  character(len=64), parameter :: fileName = 'rrtmgp-inputs-outputs.nc'
  ! ==========================================================
  ! get started
  print *, 'test start ...'

  ! load profile
  call read_atmos(fileName, p_lay, t_lay, p_lev, t_lev, gas_concs, col_dry)

  ncol = size(p_lay, 1)
  nlay = size(p_lay, 2)
  allocate(vmrh2o(ncol, nlay))
  call stop_on_err(gas_concs%get_vmr('h2o', vmrh2o))
  col_dry =  get_col_dry(vmrh2o, p_lev, t_lay)

  !
  ! Code from ~line 285 of mo_gas_optics.F90
  !   Interpolation of level temperatures from layer Ts, pressures
  !
  do icol = 1, ncol
    t_lev(icol,1) = t_lay(icol,1) &
                      + (p_lev(icol,1)-p_lay(icol,1))*(t_lay(icol,2)-t_lay(icol,1))  &
         &                                           / (p_lay(icol,2)-p_lay(icol,1))
  end do
  do ilay = 2, nlay
    do icol = 1, ncol
      t_lev(icol,ilay) = (p_lay(icol,ilay-1)*t_lay(icol,ilay-1)*(p_lev(icol,ilay  )-p_lay(icol,ilay)) &
                           +  p_lay(icol,ilay  )*t_lay(icol,ilay  )*(p_lay(icol,ilay-1)-p_lev(icol,ilay))) /  &
                             (p_lev(icol,ilay)*(p_lay(icol,ilay-1) - p_lay(icol,ilay)))
    end do
  end do
  do icol = 1, ncol
    t_lev(icol,nlay+1) = t_lay(icol,nlay)                                                             &
                           + (p_lev(icol,nlay+1)-p_lay(icol,nlay))*(t_lay(icol,nlay)-t_lay(icol,nlay-1))  &
                                                                 / (p_lay(icol,nlay)-p_lay(icol,nlay-1))
  end do


  call write_atmos(fileName, t_lev, col_dry)
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
end program test_gas_optics
