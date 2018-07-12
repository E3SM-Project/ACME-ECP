program test_cld_samp
! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2015-2016,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!
! Description: Unit test for cloud sampling/overlap
!  Tests against analytic results
!
  use mo_rte_kind,      only: wp
  use mo_rng,              only: ty_rng
  use mo_rng_MT19937,      only: ty_rng_mt
  use mo_cloud_sampling,   only: set_overlap, sample_clouds, MAX_OVERLAP, RAN_OVERLAP, MAX_RAN_OVERLAP
  implicit none 
  
  integer, parameter :: ncol = 6, nlay = 5
  integer, parameter :: nsamps = 1000
  
  real(wp), dimension(ncol,nlay) :: cld_frac_in, cld_frac_samp, cld_frac_err
  logical,  dimension(nsamps,nlay,ncol) :: cld_mask
  real(wp), dimension(ncol) :: proj_cld_frac_in, proj_cld_frac_samp

  
  integer, parameter               :: seed_length = 10
  integer, dimension(seed_length)  :: seed
  type(ty_rng_mt), dimension(ncol) :: rngs
  
  character(len=128)               :: error_msg

  
  logical :: top_at_1
  integer :: icol, ilay 
  ! --------------------------------------------------------
  cld_frac_in = RESHAPE([1.00_wp, 0.00_wp, 0.00_wp, 0.00_wp, 0.00_wp, 0.00_wp, &
                         0.00_wp, 1.00_wp, 0.00_wp, 0.00_wp, 0.00_wp, 0.00_wp, &
                         0.00_wp, 0.00_wp, 0.50_wp, 0.50_wp, 0.00_wp, 0.00_wp, &
                         0.00_wp, 0.00_wp, 0.00_wp, 0.50_wp, 0.25_wp, 0.00_wp, &
                         0.00_wp, 0.00_wp, 0.00_wp, 0.00_wp, 0.25_wp, 0.00_wp], & 
                         shape = [ncol, nlay])
  
  do icol = 1, ncol
    seed(:) = icol
    call rngs(icol)%init_rng(seed) 
  end do 

  ! -----------------
  print *, "Testing max overlap"
  error_msg = set_overlap(MAX_OVERLAP)
  if(len_trim(error_msg) > 0) then 
    print *, error_msg
    stop
  end if 
  ! Raisanen et al., 2004, https://dx.doi.org/10.1256/qj.03.99, Eq 4. 
  proj_cld_frac_in(1:ncol) = maxval(cld_frac_in(1:ncol,1:nlay),dim=2)

  print *, "  top_at_1 = false"             
  top_at_1 = .false. 
  call compute_cld_frac_error
  call report_local_cloud_frac_error
  call report_proj_cloud_frac_error
  
  print *, "  top_at_1 = true"             
  top_at_1 = .true. 
  call compute_cld_frac_error
  call report_local_cloud_frac_error
  call report_proj_cloud_frac_error
  ! -----------------
  print *
  print *, "Testing ran overlap"
  error_msg = set_overlap(RAN_OVERLAP)
  if(len_trim(error_msg) > 0) then 
    print *, error_msg
    stop
  end if 
  proj_cld_frac_in(1:ncol) = cld_frac_in(:, 1) 
  ! Raisanen et al., 2004, https://dx.doi.org/10.1256/qj.03.99, Eq 5 
  do ilay = 1, nlay
    proj_cld_frac_in(1:ncol) = proj_cld_frac_in(1:ncol) + cld_frac_in(:, ilay) - & 
                               proj_cld_frac_in(1:ncol) * cld_frac_in(:, ilay)
  end do 
  print *, "  top_at_1 = false"             
  top_at_1 = .false. 
  call compute_cld_frac_error
  call report_local_cloud_frac_error
  call report_proj_cloud_frac_error
  
  print *, "  top_at_1 = true"             
  top_at_1 = .true. 
  call compute_cld_frac_error
  call report_local_cloud_frac_error
  call report_proj_cloud_frac_error
  ! -----------------
  print *
  print *, "Testing max-ran overlap"
  error_msg = set_overlap(MAX_RAN_OVERLAP)
  if(len_trim(error_msg) > 0) then 
    print *, error_msg
    stop
  end if 
  ! Raisanen et al., 2004, https://dx.doi.org/10.1256/qj.03.99, Eq 5 
  proj_cld_frac_in(1:ncol) = 1._wp
  do ilay = 2, nlay
    proj_cld_frac_in(1:ncol) = proj_cld_frac_in(1:ncol) *                                 &  
                               (1._wp - max(cld_frac_in(:,ilay-1),cld_frac_in(:,ilay))) / & 
                               (1._wp -     cld_frac_in(:,ilay-1)) 
  end do 
  proj_cld_frac_in(1:ncol) = 1._wp - (1._wp - cld_frac_in(:,1)) * proj_cld_frac_in(1:ncol)
  print *, "  top_at_1 = false"             
  top_at_1 = .false. 
  call compute_cld_frac_error
  call report_local_cloud_frac_error
  call report_proj_cloud_frac_error
  
  print *, "  top_at_1 = true"             
  top_at_1 = .true. 
  call compute_cld_frac_error
  call report_local_cloud_frac_error
  call report_proj_cloud_frac_error
          
contains 
  ! -----------------------------------------------------------------------------
  subroutine compute_cld_frac_error 
    error_msg = sample_clouds(ncol,nlay,nsamps,rngs,cld_frac_in,top_at_1,cld_mask)
    if(len_trim(error_msg) > 0) then 
      print *, "Ack! Sampling failed." 
      print *, error_msg
      stop
    end if 
    cld_frac_samp(1:ncol,1:nlay) = transpose(real(count(cld_mask, dim=1),wp)/real(nsamps,wp))
    cld_frac_err = cld_frac_in - cld_frac_samp
  end subroutine compute_cld_frac_error 
  ! -----------------------------------------------------------------------------
  subroutine report_local_cloud_frac_error 
    print *, "    Local cloud fraction errors:"  
    print *, "    range = ", minval(cld_frac_err), maxval(cld_frac_err)
    print *, "    RMS   = ", sqrt(sum(cld_frac_err*cld_frac_err)/real(ncol*nlay))
  end subroutine report_local_cloud_frac_error 
  ! -----------------------------------------------------------------------------
  subroutine report_proj_cloud_frac_error 
    
    proj_cld_frac_samp =  real(count(any(cld_mask, dim=2),dim=1))/real(nsamps) 
    print *, "    Projected cloud fraction errors:"  
    write(*, '("    Ref :", 16(f5.3, 2x))') proj_cld_frac_in  (1:ncol)
    write(*, '("    TST :", 16(f5.3, 2x))') proj_cld_frac_samp(1:ncol)
  end subroutine report_proj_cloud_frac_error 
  ! -----------------------------------------------------------------------------
end program test_cld_samp
