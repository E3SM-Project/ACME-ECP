subroutine qneg3 (subnam  ,idx     ,ncol    ,ncold   ,lver    ,lconst_beg  , &
                  lconst_end       ,qmin    ,q       ,lfix)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Check moisture and tracers for minimum value, reset any below
! minimum value to minimum value and return information to allow
! warning message to be printed. The global average is NOT preserved.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: J. Rosinski
! 
! Modifications: 
!
!  2016-08   Kai Zhang (kai.zhang@pnnl.gov) 
!            Added an option to print out the information about negative 
!            values without clipping the tracer concentration. This can
!            be used together with subroutine massborrow.  
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use cam_logfile,  only: iulog
   use phys_control, only: use_mass_borrower, print_fixer_message 

   implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   character*(*), intent(in) :: subnam ! name of calling routine

   integer, intent(in) :: idx          ! chunk/latitude index
   integer, intent(in) :: ncol         ! number of atmospheric columns
   integer, intent(in) :: ncold        ! declared number of atmospheric columns
   integer, intent(in) :: lver         ! number of vertical levels in column
   integer, intent(in) :: lconst_beg   ! beginning constituent
   integer, intent(in) :: lconst_end   ! ending    constituent
   logical, intent(in) :: lfix         ! if true, fix negative tracers

   real(r8), intent(in) :: qmin(lconst_beg:lconst_end)      ! Global minimum constituent concentration

!
! Input/Output arguments
!
   real(r8), intent(inout) :: q(ncold,lver,lconst_beg:lconst_end) ! moisture/tracer field
!
!---------------------------Local workspace-----------------------------
!
   integer indx(ncol,lver)  ! array of indices of points < qmin
   integer nval(lver)       ! number of points < qmin for 1 level
   integer nvals            ! number of values found < qmin
   integer nn
   integer iwtmp
   integer i,ii,k           ! longitude, level indices
   integer m                ! constituent index
   integer iw,kw            ! i,k indices of worst violator

   logical found            ! true => at least 1 minimum violator found

   real(r8) worst           ! biggest violator
!
!-----------------------------------------------------------------------
!

   do m=lconst_beg,lconst_end
      nvals = 0
      found = .false.
      worst = 1.e35_r8
      iw = -1
!
! Test all field values for being less than minimum value. Set q = qmin
! for all such points. Trace offenders and identify worst one.
!
!DIR$ preferstream
      do k=1,lver
         nval(k) = 0
!DIR$ prefervector
         nn = 0
         do i=1,ncol
            if (q(i,k,m) < qmin(m)) then
               nn = nn + 1
               indx(nn,k) = i
            end if
         end do
         nval(k) = nn
      end do

      do k=1,lver
         if (nval(k) > 0) then
            found = .true.
            nvals = nvals + nval(k)
            iwtmp = -1
!cdir nodep,altcode=loopcnt
            do ii=1,nval(k)
               i = indx(ii,k)
               if (q(i,k,m) < worst) then
                  worst = q(i,k,m)
                  iwtmp = ii
               end if
            end do
            if (iwtmp /= -1 ) kw = k
            if (iwtmp /= -1 ) iw = indx(iwtmp,k)
            if(lfix) then 
!cdir nodep,altcode=loopcnt
               do ii=1,nval(k)
                  i = indx(ii,k)
                  q(i,k,m) = qmin(m)
               end do
            end if
         end if
      end do

      if(print_fixer_message .and. found .and. abs(worst)>max(qmin(m),1.e-8_r8) ) then 
         if (lfix) then 
            write(iulog,9000)subnam,m,idx,nvals,qmin(m),worst,iw,kw
         else
            write(iulog,8000)subnam,m,idx,nvals,worst,iw,kw
         end if
      end if

   end do
!
   return
8000 format(' QNEG3 from ',a,':m=',i3,' lat/lchnk=',i7, &
            ' Min. mixing ratio violated at ',i4,' points. ', &
            ' Worst =',e8.1,' at i,k=',i4,i3)
9000 format(' QNEG3 from ',a,':m=',i3,' lat/lchnk=',i7, &
            ' Min. mixing ratio violated at ',i4,' points.  Reset to ', &
            1p,e8.1,' Worst =',e8.1,' at i,k=',i4,i3)
end subroutine qneg3



